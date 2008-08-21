/**
 * MLP encoder
 * Copyright (c) 2008 Ramiro Polla <ramiro@lisha.ufsc.br>
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "avcodec.h"
#include "bitstream.h"
#include "libavutil/crc.h"
#include "mlp.h"
#include "dsputil.h"
#include "lpc.h"

#define MAJOR_HEADER_INTERVAL 16

#define MLP_MIN_LPC_ORDER      0
#define MLP_MAX_LPC_ORDER      8
#define MLP_MAX_LPC_SHIFT     16

typedef struct {
    uint8_t         min_channel;         ///< The index of the first channel coded in this substream.
    uint8_t         max_channel;         ///< The index of the last channel coded in this substream.
    uint8_t         max_matrix_channel;  ///< The number of channels input into the rematrix stage.

    uint8_t         noise_shift;         ///< The left shift applied to random noise in 0x31ea substreams.
    uint32_t        noisegen_seed;       ///< The current seed value for the pseudorandom noise generator(s).

    int             data_check_present;  ///< Set if the substream contains extra info to check the size of VLC blocks.

    int32_t         lossless_check_data; ///< XOR of all output samples
} RestartHeader;

typedef struct {
    uint16_t        blocksize;              ///< number of PCM samples in current audio block
    uint8_t         quant_step_size[MAX_CHANNELS];  ///< left shift to apply to Huffman-decoded residuals

    uint8_t         num_primitive_matrices; ///< number of matrices to apply

    int32_t         matrix_coeff[MAX_MATRICES][MAX_CHANNELS+2];
    uint8_t         frac_bits[MAX_CHANNELS];

    int8_t          output_shift[MAX_CHANNELS]; ///< Left shift to apply to decoded PCM values to get final 24-bit output.

    uint8_t         param_presence_flags;       ///< Bitmask of which parameter sets are conveyed in a decoding parameter block.
#define PARAM_PRESENCE_FLAGS    (1 << 8)

#define PARAMS_DEFAULT      (0xFF)
#define PARAM_BLOCKSIZE     (1 << 7)
#define PARAM_MATRIX        (1 << 6)
#define PARAM_OUTSHIFT      (1 << 5)
#define PARAM_QUANTSTEP     (1 << 4)
#define PARAM_FIR           (1 << 3)
#define PARAM_IIR           (1 << 2)
#define PARAM_HUFFOFFSET    (1 << 1)

} DecodingParams;

#define HUFF_OFFSET_MIN    -16384
#define HUFF_OFFSET_MAX     16383

typedef struct {
    AVCodecContext *avctx;

    int             num_substreams; ///< Number of substreams contained within this stream.

    int             num_channels;   /**< Number of channels in major_frame_buffer.
                                     *   Normal channels + noise channels. */

    int             sample_fmt;     ///< sample format encoded for MLP
    int             mlp_sample_rate;    ///< sample rate encoded for MLP

    int32_t        *sample_buffer;  ///< Pointer to current access unit samples.
    int32_t        *major_frame_buffer; ///< Buffer with all data for one entire major frame interval.
    int32_t        *last_frame;     ///< Pointer to last frame with data to encode.

    unsigned int    major_frame_size;   ///< Number of samples in current major frame being encoded.
    unsigned int    next_major_frame_size;  ///< Counter of number of samples for next major frame.

    int32_t        *lossless_check_data;    ///< Array with lossless_check_data for each access unit.

    unsigned int   *frame_size;         ///< Array with number of samples/channel in each access unit.
    unsigned int    frame_index;        ///< Index of current frame being encoded.

    unsigned int    one_sample_buffer_size; ///< Number of samples*channel for one access unit.

    unsigned int    major_header_interval;  ///< Interval of access units in between two major frames.

    uint16_t        timestamp;      ///< Timestamp of current access unit.

    uint8_t         mlp_channels;   ///< channel arrangement for MLP streams

    uint8_t         mlp_channels2;  ///< 1 bit for each channel
    uint8_t         mlp_channels3;  /**< TODO unknown channel-related field
                                     *   These values are correct for mono and stereo. */

    ChannelParams   channel_params[MAX_CHANNELS];

    DecodingParams  decoding_params[MAX_SUBSTREAMS];
    RestartHeader   restart_header[MAX_SUBSTREAMS];

    DSPContext      dsp;
} MLPEncodeContext;

#define SYNC_MAJOR      0xf8726f

#define SYNC_MLP        0xbb
#define SYNC_TRUEHD     0xba

enum InputBitDepth {
    BITS_16,
    BITS_20,
    BITS_24,
};

/** Returns the coded sample_rate for MLP. */
static int mlp_sample_rate(int sample_rate)
{
    switch (sample_rate) {
    case 44100 << 0: return 0x8 + 0;
    case 44100 << 1: return 0x8 + 1;
    case 44100 << 2: return 0x8 + 2;
    case 48000 << 0: return 0x0 + 0;
    case 48000 << 1: return 0x0 + 1;
    case 48000 << 2: return 0x0 + 2;
    default:
        return -1;
    }
}

/** Writes a major sync header to the bitstream. */
static void write_major_sync(MLPEncodeContext *ctx, uint8_t *buf, int buf_size)
{
    PutBitContext pb;

    init_put_bits(&pb, buf, buf_size);

    put_bits(&pb, 24, SYNC_MAJOR       );
    put_bits(&pb,  8, SYNC_MLP         );
    put_bits(&pb,  4, ctx->sample_fmt  );
    put_bits(&pb,  4, ctx->sample_fmt  );
    put_bits(&pb,  4, ctx->mlp_sample_rate );
    put_bits(&pb,  4, ctx->mlp_sample_rate );
    put_bits(&pb, 11, 0                ); /* This value is 0 in all tested
                                           * MLP samples. */
    put_bits(&pb,  5, ctx->mlp_channels);

    /* These values seem to be constant for all tested MLP samples. */
    put_bits(&pb, 16, 0xb752);
    put_bits(&pb, 16, 0x4000);
    put_bits(&pb, 16, 0x0000);

    put_bits(&pb,  1, 1); /* This value is 1 in all tested MLP samples.
                           * I suppose it would be 0 only when no filters
                           * or codebooks are used. */
    put_bits(&pb, 15, 0); /* TODO peak_bitrate: Most MLP samples tested encode
                           * a value that evaluates peak_bitrate to 9600000 or
                           * a little bit less. */
    put_bits(&pb,  4, 1); /* TODO Support more num_substreams. */

    put_bits(&pb, 20, 0x1054c); /* TODO These values have something to do with
                                 * the sample rate. The ones used here come
                                 * from samples that are stereo and have
                                 * 44100Hz. */
    put_bits(&pb,  8, ctx->mlp_channels2);
    put_bits(&pb, 32, 0x00008080); /* These values seem */
    put_bits(&pb,  8, 0x00      ); /* to be constants.  */
    put_bits(&pb,  8, ctx->mlp_channels3); /* TODO Finish understanding this field. */

    flush_put_bits(&pb);

    AV_WL16(buf+26, ff_mlp_checksum16(buf, 26));
}

/** Writes a restart header to the bitstream. Damaged streams can start being
 *  decoded losslessly again after such a header and the subsequent decoding
 *  params header.
 */
static void write_restart_header(MLPEncodeContext *ctx,
                                 PutBitContext *pb, int substr)
{
    RestartHeader *rh = &ctx->restart_header[substr];
    int32_t lossless_check = xor_32_to_8(rh->lossless_check_data);
    unsigned int start_count = put_bits_count(pb);
    PutBitContext tmpb;
    uint8_t checksum;
    unsigned int ch;

    put_bits(pb, 14, 0x31ea                ); /* TODO 0x31eb */
    put_bits(pb, 16, 0                     ); /* TODO I don't know what this is. Ask Ian. */
    put_bits(pb,  4, rh->min_channel       );
    put_bits(pb,  4, rh->max_channel       );
    put_bits(pb,  4, rh->max_matrix_channel);
    put_bits(pb,  4, rh->noise_shift       );
    put_bits(pb, 23, rh->noisegen_seed     );
    put_bits(pb, 19, 0                     ); /* TODO What the hell is this? */
    put_bits(pb,  1, rh->data_check_present);
    put_bits(pb,  8, lossless_check        );
    put_bits(pb, 16, 0                     ); /* This is zero =) */

    for (ch = 0; ch <= rh->max_matrix_channel; ch++)
        put_bits(pb, 6, ch);

    /* Data must be flushed for the checksum to be correct. */
    tmpb = *pb;
    flush_put_bits(&tmpb);

    checksum = ff_mlp_restart_checksum(pb->buf, put_bits_count(pb) - start_count);

    put_bits(pb,  8, checksum);
}

/** Encodes the third type of channel information for the sync headers.
 *  TODO This field is not yet fully understood. These values are just copied
 *       from some samples out in the wild.
 */
static uint8_t get_channels3_code(int channels)
{
    switch (channels) {
    case 1: return 0x1f;
    case 2: return 0x1b;
    case 6: return 0x00;
    default:
        return 0x1b;
    }
}

/** Clears a DecodingParams struct the way it should be after a restart header. */
static void clear_decoding_params(DecodingParams decoding_params[MAX_SUBSTREAMS])
{
    unsigned int substr;

    for (substr = 0; substr < MAX_SUBSTREAMS; substr++) {
        DecodingParams *dp = &decoding_params[substr];

        dp->param_presence_flags   = 0xff;
        dp->num_primitive_matrices = 0;
        dp->blocksize              = 8;

        memset(dp->frac_bits      , 0, sizeof(dp->frac_bits      ));
        memset(dp->matrix_coeff   , 0, sizeof(dp->matrix_coeff   ));
        memset(dp->output_shift   , 0, sizeof(dp->output_shift   ));
        memset(dp->quant_step_size, 0, sizeof(dp->quant_step_size));
    }
}

/** Clears a ChannelParams struct the way it should be after a restart header. */
static void clear_channel_params(ChannelParams channel_params[MAX_CHANNELS])
{
    unsigned int channel;

    for (channel = 0; channel < MAX_CHANNELS; channel++) {
        ChannelParams *cp = &channel_params[channel];

        memset(&cp->filter_params, 0, sizeof(cp->filter_params));

        /* Default audio coding is 24-bit raw PCM. */
        cp->huff_offset      =  0;
        cp->codebook         =  0;
        cp->huff_lsbs        = 24;
    }
}

/** Sets default vales in our encoder for a DecodingParams struct. */
static void default_decoding_params(MLPEncodeContext *ctx,
     DecodingParams decoding_params[MAX_SUBSTREAMS])
{
    unsigned int substr;

    clear_decoding_params(decoding_params);

    for (substr = 0; substr < MAX_SUBSTREAMS; substr++) {
        DecodingParams *dp = &decoding_params[substr];
        uint8_t param_presence_flags = 0;

        dp->num_primitive_matrices = ctx->avctx->channels - 1;

        param_presence_flags |= PARAM_BLOCKSIZE;
        param_presence_flags |= PARAM_MATRIX;
        param_presence_flags |= PARAM_OUTSHIFT;
        param_presence_flags |= PARAM_QUANTSTEP;
        param_presence_flags |= PARAM_FIR;
/*      param_presence_flags |= PARAM_IIR; */
        param_presence_flags |= PARAM_HUFFOFFSET;

        dp->param_presence_flags = param_presence_flags;
    }
}

static av_cold int mlp_encode_init(AVCodecContext *avctx)
{
    MLPEncodeContext *ctx = avctx->priv_data;
    unsigned int major_frame_buffer_size;
    unsigned int lossless_check_data_size;
    unsigned int frame_size_size;
    unsigned int substr;

    ctx->avctx = avctx;

    ctx->mlp_sample_rate = mlp_sample_rate(avctx->sample_rate);
    if (ctx->mlp_sample_rate < 0) {
        av_log(avctx, AV_LOG_ERROR, "Unsupported sample rate %d. Supported "
                            "sample rates are 44100, 88200, 176400, 48000, "
                            "96000, and 192000.\n", avctx->sample_rate);
        return -1;
    }

    /* TODO support more channels. */
    if (avctx->channels > 2) {
        av_log(avctx, AV_LOG_ERROR,
               "Only mono and stereo are supported at the moment.\n");
        return -1;
    }

    switch (avctx->sample_fmt) {
    case SAMPLE_FMT_S16: ctx->sample_fmt = BITS_16; break;
    /* TODO 20 bits: */
    case SAMPLE_FMT_S24: ctx->sample_fmt = BITS_24; break;
    default:
        av_log(avctx, AV_LOG_ERROR, "Sample format not supported. "
               "Only 16- and 24-bit samples are supported.\n");
        return -1;
    }

    avctx->frame_size               = 40 << (ctx->mlp_sample_rate & 0x7);
    avctx->coded_frame              = avcodec_alloc_frame();
    avctx->coded_frame->key_frame   = 1;

    ctx->num_channels = avctx->channels + 2; /* +2 noise channels */
    ctx->one_sample_buffer_size = avctx->frame_size
                                * ctx->num_channels;
    /* TODO Let user pass major header interval as parameter. */
    ctx->major_header_interval = MAJOR_HEADER_INTERVAL;

    /* TODO Let user pass parameters for LPC filter. */

    major_frame_buffer_size = ctx->one_sample_buffer_size
                           * ctx->major_header_interval * sizeof(int32_t);

    ctx->major_frame_buffer = av_malloc(major_frame_buffer_size);
    if (!ctx->major_frame_buffer) {
        av_log(avctx, AV_LOG_ERROR,
               "Not enough memory for buffering samples.\n");
        return -1;
    }

    ff_mlp_init_crc();
    ff_mlp_init_crc2D(NULL);

    /* TODO mlp_channels is more complex, but for now
     * we only accept mono and stereo. */
    ctx->mlp_channels   = avctx->channels - 1;
    ctx->mlp_channels2  = (1 << avctx->channels) - 1;
    ctx->mlp_channels3  = get_channels3_code(avctx->channels);
    ctx->num_substreams = 1;

    frame_size_size = sizeof(unsigned int)
                    * ctx->major_header_interval;

    ctx->frame_size = av_malloc(frame_size_size);
    if (!ctx->frame_size)
        return -1;

    lossless_check_data_size = sizeof(int32_t) * ctx->num_substreams
                             * ctx->major_header_interval;

    ctx->lossless_check_data = av_malloc(lossless_check_data_size);
    if (!ctx->lossless_check_data)
        return -1;

    for (substr = 0; substr < ctx->num_substreams; substr++) {
        RestartHeader  *rh = &ctx->restart_header [substr];

        /* TODO see if noisegen_seed is really worth it. */
        rh->noisegen_seed      = 0;

        rh->min_channel        = 0;
        rh->max_channel        = avctx->channels - 1;
        rh->max_matrix_channel = 1;
    }

    default_decoding_params(ctx, ctx->decoding_params);
    clear_channel_params(ctx->channel_params);

    dsputil_init(&ctx->dsp, avctx);

    return 0;
}

/** Calculates the smallest number of bits it takes to encode a given signed
 *  value in two's complement.
 */
static int inline number_sbits(int number)
{
    if (number < 0)
        number++;

    return av_log2(FFABS(number)) + 1 + !!number;
}

/** Determines the smallest number of bits needed to encode the filter
 *  coefficients, and if it's possible to right-shift their values without
 *  losing any precision.
 */
static void code_filter_coeffs(MLPEncodeContext *ctx,
                               FilterParams *fp,
                               int *pcoeff_shift, int *pcoeff_bits)
{
    int min = INT_MAX, max = INT_MIN;
    int bits, shift;
    int or = 0;
    int order;

    for (order = 0; order < fp->order; order++) {
        int coeff = fp->coeff[order];

        if (coeff < min)
            min = coeff;
        if (coeff > max)
            max = coeff;

        or |= coeff;
    }

    bits = FFMAX(number_sbits(min), number_sbits(max));

    for (shift = 0; shift < 7 && !(or & (1<<shift)); shift++);

    *pcoeff_bits  = bits;
    *pcoeff_shift = shift;
}

/** Writes filter parameters for one filter to the bitstream. */
static void write_filter_params(MLPEncodeContext *ctx, PutBitContext *pb,
                                unsigned int channel, unsigned int filter)
{
    FilterParams *fp = &ctx->channel_params[channel].filter_params[filter];

    put_bits(pb, 4, fp->order);

    if (fp->order > 0) {
        int coeff_shift;
        int coeff_bits;
        int i;

        code_filter_coeffs(ctx, fp, &coeff_shift, &coeff_bits);

        put_bits(pb, 4, fp->shift  );
        put_bits(pb, 5, coeff_bits );
        put_bits(pb, 3, coeff_shift);

        for (i = 0; i < fp->order; i++) {
            put_sbits(pb, coeff_bits, fp->coeff[i] >> coeff_shift);
        }

        /* TODO state data for IIR filter. */
        put_bits(pb, 1, 0);
    }
}

/** Writes matrix params for all primitive matrices to the bitstream. */
static void write_matrix_params(MLPEncodeContext *ctx, PutBitContext *pb,
                                unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    unsigned int mat;

    put_bits(pb, 4, dp->num_primitive_matrices);

    for (mat = 1; mat <= dp->num_primitive_matrices; mat++) {
        unsigned int channel;

        put_bits(pb, 4, mat               ); /* matrix_out_ch */
        put_bits(pb, 4, dp->frac_bits[mat]);
        put_bits(pb, 1, 0                 ); /* lsb_bypass */

        for (channel = 0; channel < ctx->num_channels; channel++) {
            int32_t coeff = dp->matrix_coeff[mat][channel];

            if (coeff) {
                put_bits(pb, 1, 1);

                coeff >>= 14 - dp->frac_bits[mat];

                put_sbits(pb, dp->frac_bits[mat] + 2, coeff);
            } else {
                put_bits(pb, 1, 0);
            }
        }
    }
}

/** Writes decoding parameters to the bitstream. These change very often,
 *  usually at almost every frame.
 */
static void write_decoding_params(MLPEncodeContext *ctx, PutBitContext *pb,
                                  unsigned int substr, int params_changed)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    RestartHeader  *rh = &ctx->restart_header [substr];
    unsigned int ch;

    if (dp->param_presence_flags != PARAMS_DEFAULT &&
        params_changed & PARAM_PRESENCE_FLAGS) {
        put_bits(pb, 1, 1);
        put_bits(pb, 8, dp->param_presence_flags);
    } else {
        put_bits(pb, 1, 0);
    }

    if (dp->param_presence_flags & PARAM_BLOCKSIZE) {
        if (params_changed       & PARAM_BLOCKSIZE) {
            put_bits(pb, 1, 1);
            put_bits(pb, 9, dp->blocksize);
        } else {
            put_bits(pb, 1, 0);
        }
    }

    if (dp->param_presence_flags & PARAM_MATRIX) {
        if (params_changed       & PARAM_MATRIX) {
            put_bits(pb, 1, 1);
            write_matrix_params(ctx, pb, substr);
        } else {
            put_bits(pb, 1, 0);
        }
    }

    if (dp->param_presence_flags & PARAM_OUTSHIFT) {
        if (params_changed       & PARAM_OUTSHIFT) {
            put_bits(pb, 1, 1);
            for (ch = 0; ch <= rh->max_matrix_channel; ch++)
                put_sbits(pb, 4, dp->output_shift[ch]);
        } else {
            put_bits(pb, 1, 0);
        }
    }

    if (dp->param_presence_flags & PARAM_QUANTSTEP) {
        if (params_changed       & PARAM_QUANTSTEP) {
            put_bits(pb, 1, 1);
            for (ch = 0; ch <= rh->max_channel; ch++)
                put_bits(pb, 4, dp->quant_step_size[ch]);
        } else {
            put_bits(pb, 1, 0);
        }
    }

    for (ch = rh->min_channel; ch <= rh->max_channel; ch++) {
        ChannelParams *cp = &ctx->channel_params[ch];

        if (dp->param_presence_flags & 0xF) {
            put_bits(pb, 1, 1);

            if (dp->param_presence_flags & PARAM_FIR) {
                if (params_changed       & PARAM_FIR) {
                    put_bits(pb, 1, 1);
                    write_filter_params(ctx, pb, ch, FIR);
                } else {
                    put_bits(pb, 1, 0);
                }
            }

            if (dp->param_presence_flags & PARAM_IIR) {
                if (params_changed       & PARAM_IIR) {
                    put_bits(pb, 1, 1);
                    write_filter_params(ctx, pb, ch, IIR);
                } else {
                    put_bits(pb, 1, 0);
                }
            }

            if (dp->param_presence_flags & PARAM_HUFFOFFSET) {
                if (params_changed       & PARAM_HUFFOFFSET) {
                    put_bits(pb,  1, 1);
                    put_sbits(pb, 15, cp->huff_offset);
                } else {
                    put_bits(pb, 1, 0);
                }
            }

            put_bits(pb, 2, cp->codebook );
            put_bits(pb, 5, cp->huff_lsbs);
        } else {
            put_bits(pb, 1, 0);
        }
    }
}

/** Inputs data from the samples passed by lavc into the context, shifts them
 *  appropriately depending on the bit-depth, and calculates the
 *  lossless_check_data that will be written to the restart header.
 */
static void input_data_internal(MLPEncodeContext *ctx, const uint8_t *samples,
                                int is24)
{
    int32_t *lossless_check_data = ctx->lossless_check_data;
    const int32_t *samples_32 = (const int32_t *) samples;
    const int16_t *samples_16 = (const int16_t *) samples;
    unsigned int substr;

    lossless_check_data += ctx->frame_index * ctx->num_substreams;

    for (substr = 0; substr < ctx->num_substreams; substr++) {
        RestartHeader  *rh = &ctx->restart_header [substr];
        int32_t *sample_buffer = ctx->sample_buffer;
        int32_t temp_lossless_check_data = 0;
        unsigned int channel;
        int i;

        for (i = 0; i < ctx->frame_size[ctx->frame_index]; i++) {
            for (channel = 0; channel <= rh->max_channel; channel++) {
                int32_t sample;

                if (is24) sample = *samples_32++ >> 8;
                else      sample = *samples_16++ << 8;

                temp_lossless_check_data ^= (sample & 0x00ffffff) << channel;
                *sample_buffer++ = sample;
            }
            sample_buffer += 2; /* noise channels */
        }

        *lossless_check_data++ = temp_lossless_check_data;
    }
}

/** Wrapper function for inputting data in two different bit-depths. */
static void input_data(MLPEncodeContext *ctx, void *samples)
{
    if (ctx->avctx->sample_fmt == SAMPLE_FMT_S24)
        input_data_internal(ctx, samples, 1);
    else
        input_data_internal(ctx, samples, 0);
}

/** Counts the number of trailing zeroes in a value */
static int number_trailing_zeroes(int32_t sample)
{
    int bits;

    for (bits = 0; bits < 24 && !(sample & (1<<bits)); bits++);

    /* All samples are 0. TODO Return previous quant_step_size to avoid
     * writing a new header. */
    if (bits == 24)
        return 0;

    return bits;
}

/** Determines how many bits are zero at the end of all samples so they can be
 *  shifted out for the huffman coder.
 */
static void determine_quant_step_size(MLPEncodeContext *ctx, unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    RestartHeader  *rh = &ctx->restart_header [substr];
    int32_t *sample_buffer = ctx->sample_buffer;
    int32_t sample_mask[MAX_CHANNELS];
    unsigned int channel;
    int i;

    memset(sample_mask, 0x00, sizeof(sample_mask));

    for (i = 0; i < ctx->major_frame_size; i++) {
        for (channel = 0; channel <= rh->max_channel; channel++)
            sample_mask[channel] |= *sample_buffer++;

        sample_buffer += 2; /* noise channels */
    }

    for (channel = 0; channel <= rh->max_channel; channel++)
        dp->quant_step_size[channel] = number_trailing_zeroes(sample_mask[channel]);
}

/** Determines the best filter parameters for the given data and writes the
 *  necessary information to the context.
 *  TODO Add IIR filter predictor!
 */
static void set_filter_params(MLPEncodeContext *ctx,
                              unsigned int channel, unsigned int filter,
                              int clear_filter)
{
    FilterParams *fp = &ctx->channel_params[channel].filter_params[filter];

    if (clear_filter) {
        fp->order    =  0;
        return;
    }

    if (filter == FIR) {
        int32_t *sample_buffer = ctx->sample_buffer + channel;
        int32_t coefs[MAX_LPC_ORDER][MAX_LPC_ORDER];
        /* TODO Should ctx->major_frame_buffer be reorded in channels, so that
         *      this buffer becomes unnecessary (but then every access to the
         *      same offset in all channels will span over several Kbs), or
         *      should this be a new buffer allocated in the context, or...?
         *      If it stays on the stack, there will be a limit to
         *      major_header_interval.
         */
        int32_t samples[ctx->major_frame_size];
        int32_t *lpc_samples = samples;
        int shift[MLP_MAX_LPC_ORDER];
        unsigned int i;
        int order;

        for (i = 0; i < ctx->major_frame_size; i++) {
            *lpc_samples++ = *sample_buffer;
            sample_buffer += ctx->num_channels;
        }

        order = ff_lpc_calc_coefs(&ctx->dsp, samples, ctx->major_frame_size,
                                  MLP_MIN_LPC_ORDER, MLP_MAX_LPC_ORDER, 7,
                                  coefs, shift, 1,
                                  ORDER_METHOD_EST, MLP_MAX_LPC_SHIFT, 0);

        fp->order = order;
        fp->shift = shift[order-1];

        for (i = 0; i < order; i++)
            fp->coeff[i] = coefs[order-1][i];
    } else { /* IIR */
        fp->order    =  0;
        fp->shift    =  0;
    }
}

#define INT24_MAX ((1 << 23) - 1)
#define INT24_MIN (~INT24_MAX)

#define MSB_MASK(bits)  (-1u << bits)

/** Applies the filter to the current samples, and saves the residual back
 *  into the samples buffer. If the filter is too bad and overflows the
 *  maximum amount of bits allowed (24), the samples buffer is left as is and
 *  the function returns -1.
 */
static int apply_filter(MLPEncodeContext *ctx, unsigned int substr,
                        unsigned int channel)
{
    FilterParams *fp[NUM_FILTERS] = { &ctx->channel_params[channel].filter_params[FIR],
                                      &ctx->channel_params[channel].filter_params[IIR], };
    int32_t filter_state_buffer[NUM_FILTERS][ctx->major_frame_size];
    int32_t mask = MSB_MASK(ctx->decoding_params[substr].quant_step_size[channel]);
    int32_t *sample_buffer = ctx->sample_buffer + channel;
    unsigned int major_frame_size = ctx->major_frame_size;
    unsigned int filter_shift = fp[FIR]->shift;
    int filter;
    int i;

    for (i = 0; i < 8; i++) {
        filter_state_buffer[FIR][i] = *sample_buffer;
        filter_state_buffer[IIR][i] = *sample_buffer;

        sample_buffer += ctx->num_channels;
    }

    for (i = 8; i < major_frame_size; i++) {
        int32_t sample = *sample_buffer;
        unsigned int order;
        int64_t accum = 0;
        int32_t residual;

        for (filter = 0; filter < NUM_FILTERS; filter++)
            for (order = 0; order < fp[filter]->order; order++)
                accum += (int64_t)filter_state_buffer[filter][i - 1 - order] *
                         fp[filter]->coeff[order];

        accum  >>= filter_shift;
        residual = sample - (accum & mask);

        if (residual < INT24_MIN || residual > INT24_MAX)
            return -1;

        filter_state_buffer[FIR][i] = sample;
        filter_state_buffer[IIR][i] = residual;

        sample_buffer += ctx->num_channels;
    }

    sample_buffer = ctx->sample_buffer + channel;
    for (i = 0; i < major_frame_size; i++) {
        *sample_buffer = filter_state_buffer[IIR][i];

        sample_buffer += ctx->num_channels;
    }

    return 0;
}

/** Generates two noise channels worth of data. */
static void generate_2_noise_channels(MLPEncodeContext *ctx, unsigned int substr)
{
    int32_t *sample_buffer = ctx->sample_buffer + ctx->num_channels - 2;
    RestartHeader *rh = &ctx->restart_header[substr];
    unsigned int i;
    uint32_t seed = rh->noisegen_seed;

    for (i = 0; i < ctx->major_frame_size; i++) {
        uint16_t seed_shr7 = seed >> 7;
        *sample_buffer++ = ((int8_t)(seed >> 15)) << rh->noise_shift;
        *sample_buffer++ = ((int8_t) seed_shr7)   << rh->noise_shift;

        seed = (seed << 16) ^ seed_shr7 ^ (seed_shr7 << 5);

        sample_buffer += ctx->num_channels - 2;
    }

    rh->noisegen_seed = seed & ((1 << 24)-1);
}

/** Determines how many fractional bits are needed to encode matrix
 *  coefficients. Also shifts the coefficients to fit within 2.14 bits.
 */
static int code_matrix_coeffs(MLPEncodeContext *ctx,
                               unsigned int substr, unsigned int mat)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    int32_t min = INT32_MAX, max = INT32_MIN;
    int32_t coeff_mask = 0;
    unsigned int channel;
    unsigned int shift;
    unsigned int bits;

    /* No decorrelation for mono. */
    if (ctx->num_channels - 2 == 1)
        return 0;

    for (channel = 0; channel < ctx->num_channels; channel++) {
        int32_t coeff = dp->matrix_coeff[mat][channel];

        if (coeff < min)
            min = coeff;
        if (coeff > max)
            max = coeff;

        coeff_mask |= coeff;
    }

    shift = FFMAX(0, FFMAX(number_sbits(min), number_sbits(max)) - 16);

    if (shift) {
        for (channel = 0; channel < ctx->num_channels; channel++)
            dp->matrix_coeff[mat][channel] >>= shift;

        coeff_mask >>= shift;
    }

    for (bits = 0; bits < 14 && !(coeff_mask & (1<<bits)); bits++);

    dp->frac_bits   [mat] = 14 - bits;

    for (channel = 0; channel < ctx->num_channels; channel++)
        dp->output_shift[channel] = shift;

    return ctx->num_channels - 3;
}

/** Determines best coefficients to use for the lossless matrix. */
static void lossless_matrix_coeffs(MLPEncodeContext *ctx, unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];

    generate_2_noise_channels(ctx, substr);

    /* TODO actual decorrelation. */

    dp->matrix_coeff[1][0] =  1 << 14;
    dp->matrix_coeff[1][1] = -1 << 14;
    dp->matrix_coeff[1][2] =  0 << 14;
    dp->matrix_coeff[1][3] =  0 << 14;

    dp->num_primitive_matrices = code_matrix_coeffs(ctx, substr, 1);
}

/** Applies output_shift to all channels when it is needed because of shifted
 *  matrix coefficients.
 */
static void output_shift_channels(MLPEncodeContext *ctx, unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    int32_t *sample_buffer = ctx->sample_buffer;
    unsigned int i;

    for (i = 0; i < ctx->major_frame_size; i++) {
        unsigned int channel;

        for (channel = 0; channel < ctx->num_channels - 2; channel++) {
            *sample_buffer++ >>= dp->output_shift[channel];
        }

        sample_buffer += 2;
    }
}

/** Rematrixes all channels using chosen coefficients. */
static void rematrix_channels(MLPEncodeContext *ctx, unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    int32_t *sample_buffer = ctx->sample_buffer;
    unsigned int mat, i, maxchan;

    maxchan = ctx->num_channels;

    for (mat = 1; mat <= dp->num_primitive_matrices; mat++) {
        unsigned int msb_mask_bits = (ctx->avctx->sample_fmt == SAMPLE_FMT_S16 ? 8 : 0) - dp->output_shift[mat];
        int32_t mask = MSB_MASK(msb_mask_bits);

        sample_buffer = ctx->sample_buffer;
        for (i = 0; i < ctx->major_frame_size; i++) {
            unsigned int src_ch;
            int64_t accum = 0;

            for (src_ch = 0; src_ch < maxchan; src_ch++) {
                int32_t sample = *(sample_buffer + src_ch);
                accum += (int64_t) sample * dp->matrix_coeff[mat][src_ch];
            }
            sample_buffer[mat] = (accum >> 14) & mask;

            sample_buffer += ctx->num_channels;
        }
    }
}

/** Min and max values that can be encoded with each codebook. The values for
 *  the third codebook take into account the fact that the sign shift for this
 *  codebook is outside the coded value, so it has one more bit of precision.
 *  It should actually be -7 -> 7, shifted down by 0.5.
 */
static int codebook_extremes[3][2] = {
    {-9, 8}, {-8, 7}, {-15, 14},
};

typedef struct BestOffset {
    int16_t offset;
    int bitcount;
    int lsb_bits;
} BestOffset;

/** Determines the least amount of bits needed to encode the samples using no
 *  codebooks.
 */
static void no_codebook_bits(MLPEncodeContext *ctx, unsigned int substr,
                             unsigned int channel,
                             int32_t min, int32_t max,
                             BestOffset *bo)
{
    ChannelParams  *cp = &ctx->channel_params[channel];
    DecodingParams *dp = &ctx->decoding_params[substr];
    int16_t offset;
    int32_t unsign;
    uint32_t diff;
    int lsb_bits;

    /* Set offset inside huffoffset's boundaries by adjusting extremes
     * so that more bits are used, thus shifting the offset. */
    if (min < HUFF_OFFSET_MIN)
        max = FFMAX(max, 2 * HUFF_OFFSET_MIN - min + 1);
    if (max > HUFF_OFFSET_MAX)
        min = FFMIN(min, 2 * HUFF_OFFSET_MAX - max - 1);

    /* Determine offset and minimum number of bits. */
    diff = max - min;

    lsb_bits = number_sbits(diff) - 1;

    unsign = 1 << (lsb_bits - 1);

    /* If all samples are the same (lsb_bits == 0), offset must be
     * adjusted because of sign_shift. */
    offset = min + diff / 2 + !!lsb_bits;

    /* Check if we can use the same offset as last access_unit to save
     * on writing a new header. */
    if (lsb_bits + dp->quant_step_size[channel] == cp->huff_lsbs) {
        int16_t cur_offset = cp->huff_offset;
        int32_t cur_max    = cur_offset + unsign - 1;
        int32_t cur_min    = cur_offset - unsign;

        if (min > cur_min && max < cur_max)
            offset = cur_offset;
    }

    bo->offset   = offset;
    bo->lsb_bits = lsb_bits;
    bo->bitcount = lsb_bits * dp->blocksize;
}

/** Determines the least amount of bits needed to encode the samples using a
 *  given codebook and a given offset.
 */
static inline void codebook_bits_offset(MLPEncodeContext *ctx, unsigned int substr,
                                 unsigned int channel, int codebook,
                                 int32_t min, int32_t max, int16_t offset,
                                 BestOffset *bo, int *pnext, int up)
{
    int32_t codebook_min = codebook_extremes[codebook][0];
    int32_t codebook_max = codebook_extremes[codebook][1];
    int32_t *sample_buffer = ctx->sample_buffer + channel;
    DecodingParams *dp = &ctx->decoding_params[substr];
    int codebook_offset  = 7 + (2 - codebook);
    int32_t unsign_offset = offset;
    int lsb_bits = 0, bitcount = 0;
    int next = INT_MAX;
    int unsign, mask;
    int i;

    min -= offset;
    max -= offset;

    while (min < codebook_min || max > codebook_max) {
        lsb_bits++;
        min >>= 1;
        max >>= 1;
    }

    unsign = 1 << lsb_bits;
    mask   = unsign - 1;

    if (codebook == 2) {
        unsign_offset -= unsign;
        lsb_bits++;
    }

    for (i = 0; i < dp->blocksize; i++) {
        int32_t sample = *sample_buffer >> dp->quant_step_size[channel];
        int temp_next;

        sample -= unsign_offset;

        if (up)
            temp_next = unsign - (sample & mask);
        else
            temp_next = (sample & mask) + 1;

        if (temp_next < next)
            next = temp_next;

        sample >>= lsb_bits;

        bitcount += ff_mlp_huffman_tables[codebook][sample + codebook_offset][1];

        sample_buffer += ctx->num_channels;
    }

    bo->offset   = offset;
    bo->lsb_bits = lsb_bits;
    bo->bitcount = lsb_bits * dp->blocksize + bitcount;

    *pnext       = next;
}

/** Determines the least amount of bits needed to encode the samples using a
 *  given codebook. Searches for the best offset to minimize the bits.
 */
static inline void codebook_bits(MLPEncodeContext *ctx, unsigned int substr,
                          unsigned int channel, int codebook,
                          int average, int32_t min, int32_t max,
                          BestOffset *bo, int direction)
{
    int offset = av_clip(average, HUFF_OFFSET_MIN, HUFF_OFFSET_MAX);
    int previous_count = INT_MAX;
    int offset_min, offset_max;
    int is_greater = 0;
    int next;

    offset_min = FFMAX(min, HUFF_OFFSET_MIN);
    offset_max = FFMIN(max, HUFF_OFFSET_MAX);

    for (;;) {
        BestOffset temp_bo;

        codebook_bits_offset(ctx, substr, channel, codebook,
                             min, max, offset,
                             &temp_bo, &next, direction);

        if (temp_bo.bitcount < previous_count) {
            if (temp_bo.bitcount < bo->bitcount)
                *bo = temp_bo;

            is_greater = 0;
        } else if (++is_greater >= 3)
            break;

        previous_count = temp_bo.bitcount;

        if (direction) {
            offset += next;
            if (offset > offset_max)
                break;
        } else {
            offset -= next;
            if (offset < offset_min)
                break;
        }
    }
}

/** Determines the least amount of bits needed to encode the samples using
 *  any or no codebook.
 */
static void determine_bits(MLPEncodeContext *ctx, unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    RestartHeader  *rh = &ctx->restart_header [substr];
    unsigned int channel;

    for (channel = 0; channel <= rh->max_channel; channel++) {
        int32_t *sample_buffer = ctx->sample_buffer + channel;
        ChannelParams *cp = &ctx->channel_params[channel];
        int32_t min = INT32_MAX, max = INT32_MIN;
        int best_codebook = 0;
        int average = 0;
        BestOffset bo;
        int i;

        /* Determine extremes and average. */
        for (i = 0; i < dp->blocksize; i++) {
            int32_t sample = *sample_buffer >> dp->quant_step_size[channel];
            if (sample < min)
                min = sample;
            if (sample > max)
                max = sample;
            average += sample;
            sample_buffer += ctx->num_channels;
        }
        average /= dp->blocksize;

        no_codebook_bits(ctx, substr, channel, min, max, &bo);

        for (i = 1; i < 4; i++) {
            BestOffset temp_bo = { 0, INT_MAX, 0, };

            codebook_bits(ctx, substr, channel, i - 1, average,
                          min, max, &temp_bo, 0);
            codebook_bits(ctx, substr, channel, i - 1, average,
                          min, max, &temp_bo, 1);

            if (temp_bo.bitcount < bo.bitcount) {
                bo = temp_bo;
                best_codebook = i;
            }
        }

        /* Update context. */
        cp->huff_offset = bo.offset;
        cp->huff_lsbs   = bo.lsb_bits + dp->quant_step_size[channel];
        cp->codebook    = best_codebook;
    }
}

/** Writes the residuals to the bitstream. That is, the VLC codes from the
 *  codebooks (if any is used), and then the residual.
 */
static void write_block_data(MLPEncodeContext *ctx, PutBitContext *pb,
                             unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    RestartHeader  *rh = &ctx->restart_header [substr];
    int32_t *sample_buffer = ctx->sample_buffer;
    int32_t sign_huff_offset[MAX_CHANNELS];
    int codebook            [MAX_CHANNELS];
    int lsb_bits            [MAX_CHANNELS];
    unsigned int i, ch;

    for (ch = rh->min_channel; ch <= rh->max_channel; ch++) {
        ChannelParams *cp = &ctx->channel_params[ch];
        int sign_shift;

        lsb_bits        [ch] = cp->huff_lsbs - dp->quant_step_size[ch];
        codebook        [ch] = cp->codebook  - 1;
        sign_huff_offset[ch] = cp->huff_offset;

        sign_shift = lsb_bits[ch] - 1;

        if (codebook[ch] >= 0) {
            sign_huff_offset[ch] -= 7 << lsb_bits[ch];
            sign_shift += 2 - codebook[ch];
        }

        /* Unsign if needed. */
        if (sign_shift >= 0)
            sign_huff_offset[ch] -= 1 << sign_shift;
    }

    for (i = 0; i < dp->blocksize; i++) {
        for (ch = rh->min_channel; ch <= rh->max_channel; ch++) {
            int32_t sample = *sample_buffer++ >> dp->quant_step_size[ch];

            sample -= sign_huff_offset[ch];

            if (codebook[ch] >= 0) {
                int vlc = sample >> lsb_bits[ch];
                put_bits(pb, ff_mlp_huffman_tables[codebook[ch]][vlc][1],
                             ff_mlp_huffman_tables[codebook[ch]][vlc][0]);
            }

            put_sbits(pb, lsb_bits[ch], sample);
        }
        sample_buffer += 2; /* noise channels */
    }
}

/** Compares two FilterParams structures and returns 1 if anything has
 *  changed. Returns 0 if they are both equal.
 */
static int compare_filter_params(FilterParams *prev, FilterParams *fp)
{
    int i;

    if (prev->order != fp->order)
        return 1;

    if (!prev->order)
        return 0;

    if (prev->shift != fp->shift)
        return 1;

    for (i = 0; i < fp->order; i++)
        if (prev->coeff[i] != fp->coeff[i])
            return 1;

    return 0;
}

/** Compare two primitive matrices and returns 1 if anything has changed.
 *  Returns 0 if they are both equal.
 */
static int compare_primitive_matrices(DecodingParams *prev, DecodingParams *dp)
{
    unsigned int channel, mat;

    if (prev->num_primitive_matrices != dp->num_primitive_matrices)
        return 1;

    if (!prev->num_primitive_matrices)
        return 0;

    for (channel = 0; channel < MAX_CHANNELS; channel++)
        if (prev->frac_bits[channel] != dp->frac_bits[channel])
            return 1;

    for (mat = 0; mat < MAX_MATRICES; mat++)
        for (channel = 0; channel < MAX_CHANNELS + 2; channel++)
            if (prev->matrix_coeff[mat][channel] != dp->matrix_coeff[mat][channel])
                return 1;

    return 0;
}

/** Compares two DecodingParams and ChannelParams structures to decide if a
 *  new decoding params header has to be written.
 */
static int decoding_params_diff(MLPEncodeContext *ctx, DecodingParams *prev,
                                ChannelParams channel_params[MAX_CHANNELS],
                                unsigned int substr)
{
    DecodingParams *dp = &ctx->decoding_params[substr];
    RestartHeader  *rh = &ctx->restart_header [substr];
    unsigned int ch;
    int retval = 0;

    if (prev->param_presence_flags != dp->param_presence_flags)
        retval |= PARAM_PRESENCE_FLAGS;

    if (prev->blocksize != dp->blocksize)
        retval |= PARAM_BLOCKSIZE;

    if (compare_primitive_matrices(prev, dp))
        retval |= PARAM_MATRIX;

    for (ch = 0; ch <= rh->max_matrix_channel; ch++)
        if (prev->output_shift[ch] != dp->output_shift[ch]) {
            retval |= PARAM_OUTSHIFT;
            break;
        }

    for (ch = 0; ch <= rh->max_channel; ch++)
        if (prev->quant_step_size[ch] != dp->quant_step_size[ch]) {
            retval |= PARAM_QUANTSTEP;
            break;
        }

    for (ch = rh->min_channel; ch <= rh->max_channel; ch++) {
        ChannelParams *prev_cp = &channel_params[ch];
        ChannelParams *cp = &ctx->channel_params[ch];

        if (!(retval & PARAM_FIR) &&
            compare_filter_params(&prev_cp->filter_params[FIR],
                                  &     cp->filter_params[FIR]))
            retval |= PARAM_FIR;

        if (!(retval & PARAM_IIR) &&
            compare_filter_params(&prev_cp->filter_params[IIR],
                                  &     cp->filter_params[IIR]))
            retval |= PARAM_IIR;

        if (prev_cp->huff_offset != cp->huff_offset)
            retval |= PARAM_HUFFOFFSET;

        if (prev_cp->codebook    != cp->codebook  ||
            prev_cp->huff_lsbs   != cp->huff_lsbs  )
            retval |= 0x1;
    }

    return retval;
}

/** Writes the access unit and substream headers to the bitstream. */
static void write_frame_headers(MLPEncodeContext *ctx, uint8_t *frame_header,
                                uint8_t *substream_headers, unsigned int length,
                                uint16_t substream_data_len[MAX_SUBSTREAMS])
{
    uint16_t access_unit_header = 0;
    uint16_t parity_nibble = 0;
    unsigned int substr;

    parity_nibble  = ctx->timestamp;
    parity_nibble ^= length;

    for (substr = 0; substr < ctx->num_substreams; substr++) {
        uint16_t substr_hdr = 0;

        substr_hdr |= (0 << 15); /* extraword */
        substr_hdr |= (0 << 14); /* ??? */
        substr_hdr |= (1 << 13); /* checkdata */
        substr_hdr |= (0 << 12); /* ??? */
        substr_hdr |= (substream_data_len[substr] / 2) & 0x0FFF;

        AV_WB16(substream_headers, substr_hdr);

        parity_nibble ^= *substream_headers++;
        parity_nibble ^= *substream_headers++;
    }

    parity_nibble ^= parity_nibble >> 8;
    parity_nibble ^= parity_nibble >> 4;
    parity_nibble &= 0xF;

    access_unit_header |= (parity_nibble ^ 0xF) << 12;
    access_unit_header |= length & 0xFFF;

    AV_WB16(frame_header  , access_unit_header);
    AV_WB16(frame_header+2, ctx->timestamp    );
}

/** Tries to determine a good prediction filter, and applies it to the samples
 *  buffer if the filter is good enough. Sets the filter data to be cleared if
 *  no good filter was found.
 */
static void determine_filters(MLPEncodeContext *ctx, unsigned int substr)
{
    RestartHeader *rh = &ctx->restart_header[substr];
    int channel, filter;

    for (channel = rh->min_channel; channel <= rh->max_channel; channel++) {
        for (filter = 0; filter < NUM_FILTERS; filter++)
            set_filter_params(ctx, channel, filter, 0);
        if (apply_filter(ctx, substr, channel) < 0) {
            /* Filter is horribly wrong.
             * Clear filter params and update state. */
            set_filter_params(ctx, channel, FIR, 1);
            set_filter_params(ctx, channel, IIR, 1);
            apply_filter(ctx, substr, channel);
        }
    }
}

/** Writes the substreams data to the bitstream. */
static uint8_t *write_substrs(MLPEncodeContext *ctx, uint8_t *buf, int buf_size,
                             int restart_frame,
                             DecodingParams decoding_params[MAX_SUBSTREAMS],
                             uint16_t substream_data_len[MAX_SUBSTREAMS],
                             ChannelParams channel_params[MAX_CHANNELS])
{
    int32_t *lossless_check_data = ctx->lossless_check_data;
    unsigned int substr;
    int end = 0;

    lossless_check_data += ctx->frame_index * ctx->num_substreams;

    for (substr = 0; substr < ctx->num_substreams; substr++) {
        unsigned int subblock, num_subblocks = restart_frame;
        DecodingParams *dp = &ctx->decoding_params[substr];
        RestartHeader  *rh = &ctx->restart_header [substr];
        ChannelParams backup_cp[MAX_CHANNELS];
        int32_t *backup_sample_buffer;
        uint8_t parity, checksum;
        PutBitContext pb, tmpb;
        int params_changed;
        int last_block = 0;

        init_put_bits(&pb, buf, buf_size);

        for (subblock = 0; subblock <= num_subblocks; subblock++) {

            if (num_subblocks) {
                if (!subblock) {
                    dp->blocksize = 8;

                    backup_sample_buffer = ctx->sample_buffer;

                    memcpy(backup_cp, ctx->channel_params, sizeof(backup_cp));
                    memcpy(ctx->channel_params, channel_params, sizeof(ctx->channel_params));
                } else {
                    ctx->sample_buffer += ctx->num_channels * dp->blocksize;
                    dp->blocksize = ctx->frame_size[ctx->frame_index] - dp->blocksize;

                    memcpy(channel_params, ctx->channel_params, sizeof(ctx->channel_params));
                    memcpy(ctx->channel_params, backup_cp, sizeof(ctx->channel_params));

                    restart_frame = 0;
                }
            }

            if (ctx->frame_size[ctx->frame_index] < dp->blocksize) {
                dp->blocksize = ctx->frame_size[ctx->frame_index];
                last_block = 1;
            }

            determine_bits(ctx, substr);

            params_changed = decoding_params_diff(ctx, &decoding_params[substr],
                                                channel_params,
                                                substr);

            if (restart_frame || params_changed) {
                put_bits(&pb, 1, 1);

                if (restart_frame) {
                    put_bits(&pb, 1, 1);

                    write_restart_header(ctx, &pb, substr);
                    rh->lossless_check_data = 0;
                } else {
                    put_bits(&pb, 1, 0);
                }

                write_decoding_params(ctx, &pb, substr, params_changed);
            } else {
                put_bits(&pb, 1, 0);
            }

            if (!restart_frame)
                rh->lossless_check_data ^= *lossless_check_data++;

            write_block_data(ctx, &pb, substr);

            put_bits(&pb, 1, !restart_frame);
        }

        if (num_subblocks)
            ctx->sample_buffer = backup_sample_buffer;

        put_bits(&pb, (-put_bits_count(&pb)) & 15, 0);

        if (last_block) {
            /* TODO find a sample and implement shorten_by. */
            put_bits(&pb, 32, END_OF_STREAM);
        }

        /* Data must be flushed for the checksum and parity to be correct. */
        tmpb = pb;
        flush_put_bits(&tmpb);

        parity   = ff_mlp_calculate_parity(buf, put_bits_count(&pb) >> 3) ^ 0xa9;
        checksum = ff_mlp_checksum8       (buf, put_bits_count(&pb) >> 3);

        put_bits(&pb, 8, parity  );
        put_bits(&pb, 8, checksum);

        flush_put_bits(&pb);

        end += put_bits_count(&pb) >> 3;
        substream_data_len[substr] = end;

        buf += put_bits_count(&pb) >> 3;
    }

    return buf;
}

static int mlp_encode_frame(AVCodecContext *avctx, uint8_t *buf, int buf_size,
                            void *data)
{
    DecodingParams decoding_params[MAX_SUBSTREAMS];
    uint16_t substream_data_len[MAX_SUBSTREAMS];
    ChannelParams channel_params[MAX_CHANNELS];
    MLPEncodeContext *ctx = avctx->priv_data;
    uint8_t *buf2, *buf1, *buf0 = buf;
    int total_length = 0;
    unsigned int substr;
    int restart_frame;

    ctx->frame_index = avctx->frame_number % ctx->major_header_interval;

    ctx->sample_buffer = ctx->major_frame_buffer
                       + ctx->frame_index * ctx->one_sample_buffer_size;

    if (ctx->last_frame == ctx->sample_buffer) {
        return 0;
    }

    if (avctx->frame_number < ctx->major_header_interval) {
        if (data) {
            goto input_and_return;
        } else {
            /* There are less frames than the requested major header interval.
             * Update the context to reflect this.
             */
            ctx->major_header_interval = avctx->frame_number;
            ctx->frame_index = 0;

            ctx->sample_buffer = ctx->major_frame_buffer;
        }
    }

    if (ctx->frame_size[ctx->frame_index] > MAX_BLOCKSIZE) {
        av_log(avctx, AV_LOG_ERROR, "Invalid frame size (%d > %d)\n",
               ctx->frame_size[ctx->frame_index], MAX_BLOCKSIZE);
        return -1;
    }

    if (buf_size < 4)
        return -1;

    /* Frame header will be written at the end. */
    buf      += 4;
    buf_size -= 4;

    restart_frame = !(avctx->frame_number & (ctx->major_header_interval - 1));

    if (restart_frame) {
        if (buf_size < 28)
            return -1;
        write_major_sync(ctx, buf, buf_size);
        buf      += 28;
        buf_size -= 28;

        clear_decoding_params(decoding_params);
        clear_channel_params (channel_params );
        clear_channel_params(ctx->channel_params);

        ctx->major_frame_size = ctx->next_major_frame_size;
        ctx->next_major_frame_size = 0;

        for (substr = 0; substr < ctx->num_substreams; substr++) {
            lossless_matrix_coeffs   (ctx, substr);
            output_shift_channels    (ctx, substr);
            rematrix_channels        (ctx, substr);
            determine_quant_step_size(ctx, substr);
            determine_filters        (ctx, substr);
        }
    } else {
        /* TODO Should these be a (DecodingParams *) in the context instead of
         * memcpy'ing things around?
         */
        memcpy(decoding_params, ctx->decoding_params, sizeof(decoding_params));
        memcpy(channel_params, ctx->channel_params, sizeof(channel_params));
    }

    buf1 = buf;

    /* Substream headers will be written at the end. */
    for (substr = 0; substr < ctx->num_substreams; substr++) {
        DecodingParams *dp = &ctx->decoding_params[substr];

        dp->blocksize = ctx->frame_size[ctx->frame_index];

        buf      += 2;
        buf_size -= 2;
    }

    buf2 = buf;

    total_length = buf - buf0;

    buf = write_substrs(ctx, buf, buf_size, restart_frame, decoding_params,
                        substream_data_len, channel_params);

    total_length += buf - buf2;

    write_frame_headers(ctx, buf0, buf1, total_length / 2, substream_data_len);

    ctx->timestamp += ctx->frame_size[ctx->frame_index];

input_and_return:

    if (data) {
        ctx->frame_size[ctx->frame_index] = avctx->frame_size;
        ctx->next_major_frame_size += avctx->frame_size;
        input_data(ctx, data);
    } else if (!ctx->last_frame) {
        ctx->last_frame = ctx->sample_buffer;
    }

    return total_length;
}

static av_cold int mlp_encode_close(AVCodecContext *avctx)
{
    MLPEncodeContext *ctx = avctx->priv_data;

    av_freep(&ctx->lossless_check_data);
    av_freep(&ctx->major_frame_buffer);
    av_freep(&avctx->coded_frame);
    av_freep(&ctx->frame_size);

    return 0;
}

AVCodec mlp_encoder = {
    "mlp",
    CODEC_TYPE_AUDIO,
    CODEC_ID_MLP,
    sizeof(MLPEncodeContext),
    mlp_encode_init,
    mlp_encode_frame,
    mlp_encode_close,
    .capabilities = CODEC_CAP_SMALL_LAST_FRAME | CODEC_CAP_DELAY,
    .sample_fmts = (enum SampleFormat[]){SAMPLE_FMT_S16,SAMPLE_FMT_S24,SAMPLE_FMT_NONE},
    .long_name = NULL_IF_CONFIG_SMALL("Meridian Lossless Packing"),
};
