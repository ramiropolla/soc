/*
 * MPEG-4 ALS decoder
 * Copyright (c) 2009 Thilo Borgmann <thilo.borgmann _at_ googlemail.com>
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

/**
 * @file
 * MPEG-4 ALS decoder
 * @author Thilo Borgmann <thilo.borgmann _at_ googlemail.com>
 */


//#define DEBUG


#include "avcodec.h"
#include "get_bits.h"
#include "unary.h"
#include "mpeg4audio.h"
#include "bytestream.h"

#include "als_data.h"

enum RA_Flag {
    RA_FLAG_NONE,
    RA_FLAG_FRAMES,
    RA_FLAG_HEADER
};


typedef struct {
    int resolution;           ///< 000 = 8-bit; 001 = 16-bit; 010 = 24-bit; 011 = 32-bit
    int floating;             ///< 1 = IEEE 32-bit floating-point, 0 = integer
    int frame_length;         ///< frame length for each frame (last frame may differ)
    int ra_distance;          ///< distance between RA frames (in frames, 0...255)
    enum RA_Flag ra_flag;     ///< indicates where the size of ra units is stored
    int adapt_order;          ///< adaptive order: 1 = on, 0 = off
    int coef_table;           ///< table index of Rice code parameters
    int long_term_prediction; ///< long term prediction (LTP): 1 = on, 0 = off
    int max_order;            ///< maximum prediction order (0..1023)
    int block_switching;      ///< number of block switching levels
    int bgmc;                 ///< "Block Gilbert-Moore Code": 1 = on, 0 = off (Rice coding only)
    int sb_part;              ///< sub-block partition
    int joint_stereo;         ///< joint stereo: 1 = on, 0 = off
    int mc_coding;            ///< extended inter-channel coding (multi channel coding): 1 = on, 0 = off
    int chan_config;          ///< indicates that a chan_config_info field is present
    int chan_sort;            ///< channel rearrangement: 1 = on, 0 = off
    int rlslms;               ///< use "Recursive Least Square-Least Mean Square" predictor: 1 = on, 0 = off
    int chan_config_info;     ///< mapping of channels to loudspeaker locations. Unused until setting channel configuration is implemented.
    int *chan_pos;            ///< original channel positions
    uint32_t header_size;     ///< header size of original audio file in bytes, provided for debugging
    uint32_t trailer_size;    ///< trailer size of original audio file in bytes, provided for debugging
} ALSSpecificConfig;


typedef struct {
    int stop_flag;
    int master_channel;
    int time_diff_flag;
    int time_diff_sign;
    int time_diff_index;
    int weighting[6];
} ALSChannelData;


typedef struct {
    AVCodecContext *avctx;
    ALSSpecificConfig sconf;
    GetBitContext gb;
    unsigned int num_frames;        ///< number of frames to decode, 0 if unknown
    unsigned int cur_frame_length;  ///< length of the current frame to decode
    unsigned int last_frame_length; ///< length of the last frame to decode, 0 if unknown
    unsigned int frame_id;          ///< the frame ID / number of the current frame
    unsigned int js_switch;         ///< if true, joint-stereo decoding is enforced
    unsigned int num_blocks;        ///< number of blocks used in the current frame
    int ltp_lag_length;             ///< number of bits used for ltp lag value
    int *use_ltp;                   ///< contains use_ltp flags for all channels
    int *ltp_lag;                   ///< contains ltp lag values for all channels
    int **ltp_gain;                 ///< gain values for ltp 5-tap filter for a channel
    int *ltp_gain_buffer;           ///< contains all gain values for ltp 5-tap filter
    int32_t **quant_cof;            ///< quantized parcor coefficients for a channel
    int32_t *quant_cof_buffer;      ///< contains all quantized parcor coefficients
    int32_t **lpc_cof;              ///< coefficients of the direct form prediction filter for a channel
    int32_t *lpc_cof_buffer;        ///< contains all coefficients of the direct form prediction filter
    ALSChannelData **chan_data;     ///< channel data for multi-channel correlation
    ALSChannelData *chan_data_buffer; ///< contains channel data for all channels
    int32_t *prev_raw_samples;      ///< contains unshifted raw samples from the previous block
    int32_t **raw_samples;          ///< decoded raw samples for each channel
    int32_t *raw_buffer;            ///< contains all decoded raw samples including carryover samples
} ALSDecContext;


typedef struct {
    unsigned int block_length;      ///< number of samples within the block
    unsigned int ra_block;          ///< if true, this is a random access block
    int          const_block;       ///< if true, this is a constant value block
    int32_t      const_val;         ///< the sample value of a constant block
    int          js_blocks;         ///< true if this block contains a difference signal
    unsigned int shift_lsbs;        ///< shift of values for this block
    unsigned int opt_order;         ///< prediction order of this block
    int          store_prev_samples;///< if true, carryover samples have to be stored
    int          *use_ltp;          ///< if true, long-term prediction is used
    int          *ltp_lag;          ///< lag value for long-term prediction
    int          *ltp_gain;         ///< gain values for ltp 5-tap filter
    int32_t      *quant_cof;        ///< quantized parcor coefficients
    int32_t      *lpc_cof;          ///< coefficients of the direct form prediction
    int32_t      *raw_samples;      ///< decoded raw samples / residuals for this block
    int32_t      *prev_raw_samples; ///< contains unshifted raw samples from the previous block
    int32_t      *raw_other;        ///< decoded raw samples of the other channel of a channel pair
} ALSBlockData;


static av_cold void dprint_specific_config(ALSDecContext *ctx)
{
#ifdef DEBUG
    AVCodecContext *avctx    = ctx->avctx;
    ALSSpecificConfig *sconf = &ctx->sconf;

    dprintf(avctx, "resolution = %i\n",           sconf->resolution);
    dprintf(avctx, "floating = %i\n",             sconf->floating);
    dprintf(avctx, "frame_length = %i\n",         sconf->frame_length);
    dprintf(avctx, "ra_distance = %i\n",          sconf->ra_distance);
    dprintf(avctx, "ra_flag = %i\n",              sconf->ra_flag);
    dprintf(avctx, "adapt_order = %i\n",          sconf->adapt_order);
    dprintf(avctx, "coef_table = %i\n",           sconf->coef_table);
    dprintf(avctx, "long_term_prediction = %i\n", sconf->long_term_prediction);
    dprintf(avctx, "max_order = %i\n",            sconf->max_order);
    dprintf(avctx, "block_switching = %i\n",      sconf->block_switching);
    dprintf(avctx, "bgmc = %i\n",                 sconf->bgmc);
    dprintf(avctx, "sb_part = %i\n",              sconf->sb_part);
    dprintf(avctx, "joint_stereo = %i\n",         sconf->joint_stereo);
    dprintf(avctx, "mc_coding = %i\n",            sconf->mc_coding);
    dprintf(avctx, "chan_config = %i\n",          sconf->chan_config);
    dprintf(avctx, "chan_sort = %i\n",            sconf->chan_sort);
    dprintf(avctx, "RLSLMS = %i\n",               sconf->rlslms);
    dprintf(avctx, "chan_config_info = %i\n",     sconf->chan_config_info);
    dprintf(avctx, "header_size = %i\n",          sconf->header_size);
    dprintf(avctx, "trailer_size = %i\n",         sconf->trailer_size);

    dprintf(avctx, " num_frames = %i\n",          ctx->num_frames);
    dprintf(avctx, " last_frame_length = %i\n",   ctx->last_frame_length);
#endif
}


/** Reads an ALSSpecificConfig from a buffer into the output struct.
 */
static av_cold int read_specific_config(ALSDecContext *ctx)
{
    GetBitContext gb;
    uint64_t ht_size;
    int i, config_offset, crc_enabled;
    MPEG4AudioConfig m4ac;
    ALSSpecificConfig *sconf = &ctx->sconf;
    AVCodecContext *avctx    = ctx->avctx;
    uint32_t samples, als_id;

    init_get_bits(&gb, avctx->extradata, avctx->extradata_size * 8);

    config_offset = ff_mpeg4audio_get_config(&m4ac, avctx->extradata,
                                             avctx->extradata_size);

    if (config_offset < 0)
        return -1;

    skip_bits_long(&gb, config_offset);

    if (get_bits_left(&gb) < (30 << 3))
        return -1;

    // read the fixed items
    als_id                      = get_bits_long(&gb, 32);
    avctx->sample_rate          = m4ac.sample_rate;
    skip_bits_long(&gb, 32); // sample rate already known
    samples                     = get_bits_long(&gb, 32);
    avctx->channels             = m4ac.channels;
    skip_bits(&gb, 16);      // number of channels already knwon
    skip_bits(&gb, 3);       // skip file_type
    sconf->resolution           = get_bits(&gb, 3);
    sconf->floating             = get_bits1(&gb);
    skip_bits1(&gb);         // skip msb_first
    sconf->frame_length         = get_bits(&gb, 16) + 1;
    sconf->ra_distance          = get_bits(&gb, 8);
    sconf->ra_flag              = get_bits(&gb, 2);
    sconf->adapt_order          = get_bits1(&gb);
    sconf->coef_table           = get_bits(&gb, 2);
    sconf->long_term_prediction = get_bits1(&gb);
    sconf->max_order            = get_bits(&gb, 10);
    sconf->block_switching      = get_bits(&gb, 2);
    sconf->bgmc                 = get_bits1(&gb);
    sconf->sb_part              = get_bits1(&gb);
    sconf->joint_stereo         = get_bits1(&gb);
    sconf->mc_coding            = get_bits1(&gb);
    sconf->chan_config          = get_bits1(&gb);
    sconf->chan_sort            = get_bits1(&gb);
    crc_enabled                 = get_bits1(&gb);
    sconf->rlslms               = get_bits1(&gb);
    skip_bits(&gb, 5);       // skip 5 reserved bits
    skip_bits1(&gb);         // skip aux_data_enabled


    // check for ALSSpecificConfig struct
    if (als_id != MKBETAG('A','L','S','\0'))
        return -1;

    ctx->cur_frame_length = sconf->frame_length;

    // calculate total number of frames to decode if possible
    if (samples != 0xFFFFFFFF) {
        ctx->num_frames        = (samples - 1) / sconf->frame_length + 1;
        ctx->last_frame_length = (samples - 1) % sconf->frame_length + 1;
    } else {
        ctx->num_frames        = 0;
        ctx->last_frame_length = 0;
    }

    // read channel config
    if (sconf->chan_config)
        sconf->chan_config_info = get_bits(&gb, 16);
    // TODO: use this to set avctx->channel_layout


    // read channel sorting
    if (sconf->chan_sort && avctx->channels > 1) {
        int chan_pos_bits = av_ceil_log2(avctx->channels);
        int bits_needed  = avctx->channels * chan_pos_bits + 7;
        if (get_bits_left(&gb) < bits_needed)
            return -1;

        if (!(sconf->chan_pos = av_malloc(avctx->channels * sizeof(*sconf->chan_pos))))
            return AVERROR(ENOMEM);

        for (i = 0; i < avctx->channels; i++)
            sconf->chan_pos[i] = get_bits(&gb, chan_pos_bits);

        align_get_bits(&gb);
        // TODO: use this to actually do channel sorting
    } else {
        sconf->chan_sort = 0;
    }


    // read fixed header and trailer sizes,
    // if size = 0xFFFFFFFF then there is no data field!
    if (get_bits_left(&gb) < 64)
        return -1;

    sconf->header_size  = get_bits_long(&gb, 32);
    sconf->trailer_size = get_bits_long(&gb, 32);
    if (sconf->header_size  == 0xFFFFFFFF)
        sconf->header_size  = 0;
    if (sconf->trailer_size == 0xFFFFFFFF)
        sconf->trailer_size = 0;

    ht_size = ((int64_t)(sconf->header_size) + (int64_t)(sconf->trailer_size)) << 3;


    // skip the header and trailer data
    if (get_bits_left(&gb) < ht_size)
        return -1;

    if (ht_size > INT32_MAX)
        return -1;

    skip_bits_long(&gb, ht_size);


    // skip the crc data
    if (crc_enabled) {
        if (get_bits_left(&gb) < 32)
            return -1;

        skip_bits_long(&gb, 32);
    }


    // no need to read the rest of ALSSpecificConfig (ra_unit_size & aux data)

    dprint_specific_config(ctx);

    return 0;
}


/** Checks the ALSSpecificConfig for unsupported features.
 */
static int check_specific_config(ALSDecContext *ctx)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    int error = 0;

    // report unsupported feature and set error value
    #define MISSING_ERR(cond, str, errval)              \
    {                                                   \
        if (cond) {                                     \
            av_log_missing_feature(ctx->avctx, str, 0); \
            error = errval;                             \
        }                                               \
    }

    MISSING_ERR(sconf->floating,             "Floating point decoding",     -1);
    MISSING_ERR(sconf->bgmc,                 "BGMC entropy decoding",       -1);
    MISSING_ERR(sconf->rlslms,               "Adaptive RLS-LMS prediction", -1);
    MISSING_ERR(sconf->chan_sort,            "Channel sorting",              0);

    return error;
}


/** Parses the bs_info field to extract the block partitioning used in
 *  block switching mode, refer to ISO/IEC 14496-3, section 11.6.2.
 */
static void parse_bs_info(const uint32_t bs_info, unsigned int n,
                          unsigned int div, unsigned int **div_blocks,
                          unsigned int *num_blocks)
{
    if (n < 31 && ((bs_info << n) & 0x40000000)) {
        // if the level is valid and the investigated bit n is set
        // then recursively check both children at bits (2n+1) and (2n+2)
        n   *= 2;
        div += 1;
        parse_bs_info(bs_info, n + 1, div, div_blocks, num_blocks);
        parse_bs_info(bs_info, n + 2, div, div_blocks, num_blocks);
    } else {
        // else the bit is not set or the last level has been reached
        // (bit implicitly not set)
        **div_blocks = div;
        (*div_blocks)++;
        (*num_blocks)++;
    }
}


/** Reads and decodes a Rice codeword.
 */
static int32_t decode_rice(GetBitContext *gb, unsigned int k)
{
    int max = gb->size_in_bits - get_bits_count(gb) - k;
    int q   = get_unary(gb, 0, max);
    int r   = k ? get_bits1(gb) : !(q & 1);

    if (k > 1) {
        q <<= (k - 1);
        q  += get_bits_long(gb, k - 1);
    } else if (!k) {
        q >>= 1;
    }
    return r ? q : ~q;
}


/** Converts PARCOR coefficient k to direct filter coefficient.
 */
static void parcor_to_lpc(unsigned int k, const int32_t *par, int32_t *cof)
{
    int i, j;

    for (i = 0, j = k - 1; i < j; i++, j--) {
        int tmp1 = ((MUL64(par[k], cof[j]) + (1 << 19)) >> 20);
        cof[j]  += ((MUL64(par[k], cof[i]) + (1 << 19)) >> 20);
        cof[i]  += tmp1;
    }
    if (i == j)
        cof[i] += ((MUL64(par[k], cof[j]) + (1 << 19)) >> 20);

    cof[k] = par[k];
}


/** Reformat block sizes from log2 format to direct form. Also assure that the
 *  block sizes of the last frame correspond to the actual number of samples.
 */
static void reconstruct_block_sizes(ALSDecContext *ctx, uint32_t *div_blocks)
{
    unsigned int b;

    // The last frame may have an overdetermined block structure given in
    // the bitstream. In that case the defined block structure would need
    // more samples than available to be consistent.
    // The block structure is actually used but the block sizes are adapted
    // to fit the actual number of available samples.
    // Example: 5 samples, 2nd level block sizes: 2 2 2 2.
    // This results in the actual block sizes:    2 2 1 0.
    // This is not specified in 14496-3 but actually done by the reference
    // codec RM22 revision 2.
    // This appears to happen in case of an odd number of samples in the last
    // frame which is actually not allowed by the block length switching part
    // of 14496-3.
    // The ALS conformance files feature an odd number of samples in the last
    // frame.

    for (b = 0; b < ctx->num_blocks; b++)
        div_blocks[b] = ctx->sconf.frame_length >> div_blocks[b];

    if (ctx->cur_frame_length == ctx->last_frame_length) {
        unsigned int remaining = ctx->cur_frame_length;

        for (b = 0; b < ctx->num_blocks; b++) {
            if (remaining < div_blocks[b]) {
                div_blocks[b] = remaining;
                ctx->num_blocks = b + 1;
                break;
            }

            remaining -= div_blocks[b];
        }
    }
}


/** Reads block switching field if necessary and sets actual block sizes.
 */
static void get_block_sizes(ALSDecContext *ctx, unsigned int *div_blocks,
                            uint32_t *bs_info)
{
    ALSSpecificConfig *sconf     = &ctx->sconf;
    GetBitContext *gb            = &ctx->gb;
    unsigned int *ptr_div_blocks = div_blocks;

    if (sconf->block_switching) {
        unsigned int bs_info_len = 1 << (sconf->block_switching + 2);
        *bs_info = get_bits_long(gb, bs_info_len);
        *bs_info <<= (32 - bs_info_len);
    }

    ctx->num_blocks = 0;
    parse_bs_info(*bs_info, 0, 0, &ptr_div_blocks, &ctx->num_blocks);
    reconstruct_block_sizes(ctx, div_blocks);
}


/** Reads the block data for a constant block
 */
static void read_const_block_data(ALSDecContext *ctx, ALSBlockData *bd)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    AVCodecContext *avctx    = ctx->avctx;
    GetBitContext *gb        = &ctx->gb;

    bd->const_val    = 0;
    bd->const_block  = get_bits1(gb);    // 1 = constant value, 0 = zero block (silence)
    bd->js_blocks    = get_bits1(gb);

    // skip 5 reserved bits
    skip_bits(gb, 5);

    if (bd->const_block) {
        unsigned int const_val_bits = sconf->floating ? 24 : avctx->bits_per_raw_sample;
        bd->const_val = get_sbits_long(gb, const_val_bits);
    }

    // ensure constant block decoding by reusing this field
    bd->const_block = 1;
}


/** Decodes the block data for a constant block
 */
static void decode_const_block_data(ALSDecContext *ctx, ALSBlockData *bd)
{
    int smp;

    // write raw samples into buffer
    for (smp = 0; smp < bd->block_length; smp++)
        bd->raw_samples[smp] = bd->const_val;
}


/** Reads the block data for a non-constant block
 */
static int read_var_block_data(ALSDecContext *ctx, ALSBlockData *bd)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    AVCodecContext *avctx    = ctx->avctx;
    GetBitContext *gb        = &ctx->gb;
    unsigned int k;
    unsigned int s[8];
    unsigned int sub_blocks, log2_sub_blocks, sb_length;
    unsigned int start      = 0;
    int          sb;

    // ensure variable block decoding by reusing this field
    bd->const_block = 0;

    bd->opt_order   = 1;
    bd->js_blocks   = get_bits1(gb);

    // determine the number of subblocks for entropy decoding
    if (!sconf->bgmc && !sconf->sb_part) {
        log2_sub_blocks = 0;
    } else {
        if (sconf->bgmc && sconf->sb_part)
            log2_sub_blocks = get_bits(gb, 2);
        else
            log2_sub_blocks = 2 * get_bits1(gb);
    }

    sub_blocks = 1 << log2_sub_blocks;

    // do not continue in case of a damaged stream since
    // block_length must be evenly divisible by sub_blocks
    if (bd->block_length & (sub_blocks - 1)) {
        av_log(avctx, AV_LOG_WARNING,
               "Block length is not evenly divisible by the number of subblocks.\n");
        return -1;
    }

    sb_length = bd->block_length >> log2_sub_blocks;


    if (sconf->bgmc) {
        // TODO: BGMC mode
    } else {
        s[0] = get_bits(gb, 4 + (sconf->resolution > 1));
        for (k = 1; k < sub_blocks; k++)
            s[k] = s[k - 1] + decode_rice(gb, 0);
    }

    if (get_bits1(gb))
        bd->shift_lsbs = get_bits(gb, 4) + 1;

    bd->store_prev_samples = (bd->js_blocks && bd->raw_other) ||
                                      bd->shift_lsbs;

    if (!sconf->rlslms) {
        if (sconf->adapt_order) {
            int opt_order_length  = av_ceil_log2(av_clip((bd->block_length >> 3) - 1,
                                                 2, sconf->max_order + 1));
            bd->opt_order = get_bits(gb, opt_order_length);
        } else {
            bd->opt_order = sconf->max_order;
        }

        if (bd->opt_order) {
            int add_base;

            if (sconf->coef_table == 3) {
                add_base = 0x7F;

                // read coefficient 0
                bd->quant_cof[0]  = parcor_scaled_values[get_bits(gb, 7)];
                bd->quant_cof[0] *= 32;

                // read coefficient 1
                bd->quant_cof[1]  = -parcor_scaled_values[get_bits(gb, 7)];
                bd->quant_cof[0] *= 32;

                // read coefficients 2 to opt_order
                for (k = 2; k < bd->opt_order; k++)
                    bd->quant_cof[k] = get_bits(gb, 7);
            } else {
                int k_max;
                add_base = 1;

                // read coefficient 0 to 19
                k_max = FFMIN(bd->opt_order, 20);
                for (k = 0; k < k_max; k++) {
                    int rice_param = parcor_rice_table[sconf->coef_table][k][1];
                    int offset     = parcor_rice_table[sconf->coef_table][k][0];
                    bd->quant_cof[k] = decode_rice(gb, rice_param) + offset;
                }

                // read coefficients 20 to 126
                k_max = FFMIN(bd->opt_order, 127);
                for (; k < k_max; k++)
                    bd->quant_cof[k] = decode_rice(gb, 2) + (k & 1);

                // read coefficients 127 to opt_order
                for (; k < bd->opt_order; k++)
                    bd->quant_cof[k] = decode_rice(gb, 1);

                bd->quant_cof[0]  =  parcor_scaled_values[bd->quant_cof[0] + 64];
                bd->quant_cof[0] *= 32;
                bd->quant_cof[1]  = -parcor_scaled_values[bd->quant_cof[1] + 64];
                bd->quant_cof[1] *= 32;
            }

        for (k = 2; k < bd->opt_order; k++)
            bd->quant_cof[k] = (bd->quant_cof[k] << 14) + (add_base << 13);
        }
    }

    // read LTP gain and lag values
    if (sconf->long_term_prediction) {
        *bd->use_ltp = get_bits1(gb);

        if (*bd->use_ltp) {
            bd->ltp_gain[0]   = decode_rice(gb, 1) << 3;
            bd->ltp_gain[1]   = decode_rice(gb, 2) << 3;

            bd->ltp_gain[2]   = get_unary(gb, 0, 4);
            bd->ltp_gain[2] <<= 2;
            bd->ltp_gain[2]  += get_bits(gb, 2);
            bd->ltp_gain[2]   = ltp_gain_values[bd->ltp_gain[2]];

            bd->ltp_gain[3]   = decode_rice(gb, 2) << 3;
            bd->ltp_gain[4]   = decode_rice(gb, 1) << 3;

            *bd->ltp_lag      = get_bits(gb, ctx->ltp_lag_length);
            *bd->ltp_lag     += FFMAX(4, bd->opt_order + 1);
        }
    }

    // read first value and residuals in case of a random access block
    if (bd->ra_block) {
        if (bd->opt_order)
            bd->raw_samples[0] = decode_rice(gb, avctx->bits_per_raw_sample - 4);
        if (bd->opt_order > 1)
            bd->raw_samples[1] = decode_rice(gb, s[0] + 3);
        if (bd->opt_order > 2)
            bd->raw_samples[2] = decode_rice(gb, s[0] + 1);

        start = FFMIN(bd->opt_order, 3);
    }

    // read all residuals
    if (sconf->bgmc) {
        // TODO: BGMC mode
    } else {
        int32_t *current_res = bd->raw_samples + start;

        for (sb = 0; sb < sub_blocks; sb++, start = 0)
            for (; start < sb_length; start++)
                *current_res++ = decode_rice(gb, s[sb]);
    }

    if (!sconf->mc_coding || ctx->js_switch)
        align_get_bits(gb);

    return 0;
}


/** Decodes the block data for a non-constant block
 */
static int decode_var_block_data(ALSDecContext *ctx, ALSBlockData *bd)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    unsigned int block_length = bd->block_length;
    unsigned int smp = 0;
    unsigned int k;
    int sb;
    int64_t y;

    // reverse long-term prediction
    if (*bd->use_ltp) {
        int ltp_smp;

        for (ltp_smp = FFMAX(*bd->ltp_lag - 2, 0); ltp_smp < block_length; ltp_smp++) {
            int center = ltp_smp - *bd->ltp_lag;
            int begin  = FFMAX(0, center - 2);
            int end    = center + 3;
            int tab    = 5 - (end - begin);
            int base;

            y = 1 << 6;

            for (base = begin; base < end; base++, tab++)
                y += MUL64(bd->ltp_gain[tab], bd->raw_samples[base]);

            bd->raw_samples[ltp_smp] += y >> 7;
        }
    }

    // reconstruct all samples from residuals
    if (bd->ra_block) {
        for (smp = 0; smp < bd->opt_order; smp++) {
            y = 1 << 19;

            for (sb = 0; sb < smp; sb++)
                y += MUL64(bd->lpc_cof[sb], bd->raw_samples[smp - (sb + 1)]);

            bd->raw_samples[smp] -= y >> 20;
            parcor_to_lpc(smp, bd->quant_cof, bd->lpc_cof);
        }
    } else {
        for (k = 0; k < bd->opt_order; k++)
            parcor_to_lpc(k, bd->quant_cof, bd->lpc_cof);

        // store previous samples in case that they have to be altered
        if (bd->store_prev_samples)
            memcpy(bd->prev_raw_samples, bd->raw_samples - sconf->max_order,
                   sizeof(*bd->prev_raw_samples) * sconf->max_order);

        // reconstruct difference signal for prediction (joint-stereo)
        if (bd->js_blocks && bd->raw_other) {
            int32_t *left, *right;

            if (bd->raw_other > bd->raw_samples) {  // D = R - L
                left  = bd->raw_samples;
                right = bd->raw_other;
            } else {                                // D = R - L
                left  = bd->raw_other;
                right = bd->raw_samples;
            }

            for (sb = -1; sb >= -sconf->max_order; sb--)
                bd->raw_samples[sb] = right[sb] - left[sb];
        }

        // reconstruct shifted signal
        if (bd->shift_lsbs)
            for (sb = -1; sb >= -sconf->max_order; sb--)
                bd->raw_samples[sb] >>= bd->shift_lsbs;
    }

    // reconstruct raw samples
    for (; smp < bd->block_length; smp++) {
        y = 1 << 19;

        for (sb = 0; sb < bd->opt_order; sb++)
            y += MUL64(bd->lpc_cof[sb], bd->raw_samples[smp - (sb + 1)]);

        bd->raw_samples[smp] -= y >> 20;
    }

    // restore previous samples in case that they have been altered
    if (bd->store_prev_samples)
        memcpy(bd->raw_samples - sconf->max_order, bd->prev_raw_samples,
               sizeof(*bd->raw_samples) * sconf->max_order);

    return 0;
}


/** Reads the block data.
 */
static int read_block(ALSDecContext *ctx, ALSBlockData *bd)
{
    GetBitContext *gb        = &ctx->gb;

    // read block type flag and read the samples accordingly
    if (get_bits1(gb))
        read_var_block_data(ctx, bd);
    else
        read_const_block_data(ctx, bd);

    return 0;
}


/** Decodes the block data.
 */
static int decode_block(ALSDecContext *ctx, ALSBlockData *bd)
{
    unsigned int smp;

    // read block type flag and read the samples accordingly
    if (bd->const_block)
        decode_const_block_data(ctx, bd);
    else if (decode_var_block_data(ctx, bd))
        return -1;

    // TODO: read RLSLMS extension data

    if (bd->shift_lsbs)
        for (smp = 0; smp < bd->block_length; smp++)
            bd->raw_samples[smp] <<= bd->shift_lsbs;

    return 0;
}


/** Computes the number of samples left to decode for the current frame and
 *  sets these samples to zero.
 */
static void zero_remaining(unsigned int b, unsigned int b_max,
                           const unsigned int *div_blocks, int32_t *buf)
{
    unsigned int count = 0;

    while (b < b_max)
        count += div_blocks[b];

    memset(buf, 0, sizeof(*buf) * count);
}


/** Decodes blocks independently.
 */
static int decode_blocks_ind(ALSDecContext *ctx, unsigned int ra_frame,
                             unsigned int c, const unsigned int *div_blocks,
                             unsigned int *js_blocks)
{
    unsigned int b;
    ALSBlockData bd;

    memset(&bd, 0, sizeof(ALSBlockData));

    bd.ra_block         = ra_frame;
    bd.use_ltp          = ctx->use_ltp;
    bd.ltp_lag          = ctx->ltp_lag;
    bd.ltp_gain         = ctx->ltp_gain[0];
    bd.quant_cof        = ctx->quant_cof[0];
    bd.lpc_cof          = ctx->lpc_cof[0];
    bd.prev_raw_samples = ctx->prev_raw_samples;
    bd.raw_samples      = ctx->raw_samples[c];


    for (b = 0; b < ctx->num_blocks; b++) {
        bd.shift_lsbs       = 0;
        bd.block_length     = div_blocks[b];

        if (read_block(ctx, &bd) ||
            decode_block(ctx, &bd)) {
            zero_remaining(b, ctx->num_blocks, div_blocks, bd.raw_samples);
            return -1;
            }

        bd.raw_samples += div_blocks[b];
        bd.ra_block     = 0;
    }

    *js_blocks      = bd.js_blocks;

    return 0;
}


/** Decodes blocks dependently.
 */
static int decode_blocks(ALSDecContext *ctx, unsigned int ra_frame,
                         unsigned int c, const unsigned int *div_blocks,
                         unsigned int *js_blocks)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    unsigned int offset = 0;
    unsigned int b;
    ALSBlockData bd[2];

    memset(bd, 0, 2 * sizeof(ALSBlockData));

    bd[0].ra_block         = ra_frame;
    bd[0].use_ltp          = ctx->use_ltp;
    bd[0].ltp_lag          = ctx->ltp_lag;
    bd[0].ltp_gain         = ctx->ltp_gain[0];
    bd[0].quant_cof        = ctx->quant_cof[0];
    bd[0].lpc_cof          = ctx->lpc_cof[0];
    bd[0].prev_raw_samples = ctx->prev_raw_samples;
    bd[0].js_blocks        = *js_blocks;

    bd[1].ra_block         = ra_frame;
    bd[1].use_ltp          = ctx->use_ltp;
    bd[1].ltp_lag          = ctx->ltp_lag;
    bd[1].ltp_gain         = ctx->ltp_gain[0];
    bd[1].quant_cof        = ctx->quant_cof[0];
    bd[1].lpc_cof          = ctx->lpc_cof[0];
    bd[1].prev_raw_samples = ctx->prev_raw_samples;
    bd[1].js_blocks        = *(js_blocks + 1);

    // decode all blocks
    for (b = 0; b < ctx->num_blocks; b++) {
        unsigned int s;

        bd[0].shift_lsbs   = 0;
        bd[1].shift_lsbs   = 0;

        bd[0].block_length = div_blocks[b];
        bd[1].block_length = div_blocks[b];

        bd[0].raw_samples  = ctx->raw_samples[c    ] + offset;
        bd[1].raw_samples  = ctx->raw_samples[c + 1] + offset;

        bd[0].raw_other    = bd[1].raw_samples;
        bd[1].raw_other    = bd[0].raw_samples;

        if (read_block(ctx, &bd[0]) || decode_block(ctx, &bd[0]) ||
            read_block(ctx, &bd[1]) || decode_block(ctx, &bd[1])) {
            // damaged block, write zero for the rest of the frame
            zero_remaining(b, ctx->num_blocks, div_blocks, bd[0].raw_samples);
            zero_remaining(b, ctx->num_blocks, div_blocks, bd[1].raw_samples);
            return -1;
        }

        // reconstruct joint-stereo blocks
        if (bd[0].js_blocks) {
            if (bd[1].js_blocks)
                av_log(ctx->avctx, AV_LOG_WARNING, "Invalid channel pair!\n");

            for (s = 0; s < div_blocks[b]; s++)
                bd[0].raw_samples[s] = bd[1].raw_samples[s] - bd[0].raw_samples[s];
        } else if (bd[1].js_blocks) {
            for (s = 0; s < div_blocks[b]; s++)
                bd[1].raw_samples[s] = bd[1].raw_samples[s] + bd[0].raw_samples[s];
        }

        offset  += div_blocks[b];
        bd[0].ra_block = 0;
        bd[1].ra_block = 0;
    }

    // store carryover raw samples,
    // the others channel raw samples are stored by the calling function.
    memmove(ctx->raw_samples[c] - sconf->max_order,
            ctx->raw_samples[c] - sconf->max_order + sconf->frame_length,
            sizeof(*ctx->raw_samples[c]) * sconf->max_order);

    return 0;
}


/** Reads the channel data
  */
static void read_channel_data(ALSDecContext *ctx, ALSChannelData *cd, int c)
{
    GetBitContext *gb       = &ctx->gb;
    ALSChannelData *current = cd;

    while (!(current->stop_flag = get_bits1(gb))) {
        current->master_channel = get_bits_long(gb, av_ceil_log2(ctx->avctx->channels));

        if (current->master_channel != c) {
            current->time_diff_flag = get_bits1(gb);
            current->weighting[0]   = mcc_weightings[av_clip(decode_rice(gb, 1) + 16, 0, 32)];
            current->weighting[1]   = mcc_weightings[av_clip(decode_rice(gb, 2) + 14, 0, 32)];
            current->weighting[2]   = mcc_weightings[av_clip(decode_rice(gb, 1) + 16, 0, 32)];

            if (current->time_diff_flag) {
                current->weighting[3] = mcc_weightings[av_clip(decode_rice(gb, 1) + 16, 0, 32)];
                current->weighting[4] = mcc_weightings[av_clip(decode_rice(gb, 1) + 16, 0, 32)];
                current->weighting[5] = mcc_weightings[av_clip(decode_rice(gb, 1) + 16, 0, 32)];

                current->time_diff_sign  = get_bits1(gb);
                current->time_diff_index = get_bits(gb, ctx->ltp_lag_length - 3) + 3;
                if (current->time_diff_sign)
                    current->time_diff_index = -current->time_diff_index;
            }
        }

        current++;
    }

    align_get_bits(gb);
}


static void revert_channel_correlation(ALSDecContext *ctx, ALSBlockData *bd,
                                       ALSChannelData **cd, int *reverted,
                                       unsigned int offset, int c)
{
    ALSChannelData *ch = cd[c];
    unsigned int   dep = 0;

    if (reverted[c])
        return;

    while (!ch[dep].stop_flag) {
        if (!reverted[ch[dep].master_channel])
            revert_channel_correlation(ctx, bd, cd, reverted, offset,
                                       ch[dep].master_channel);

        dep++;
    }

    bd->use_ltp     = ctx->use_ltp + c;
    bd->ltp_lag     = ctx->ltp_lag + c;
    bd->ltp_gain    = ctx->ltp_gain[c];
    bd->lpc_cof     = ctx->lpc_cof[c];
    bd->quant_cof   = ctx->quant_cof[c];
    bd->raw_samples = ctx->raw_samples[c] + offset;

    dep = 0;
    while (!ch[dep].stop_flag) {
        unsigned int smp;
        unsigned int begin = 1;
        unsigned int end   = bd->block_length - 1;
        int64_t y;
        int32_t *master = ctx->raw_samples[ch[dep].master_channel] + offset;

        if (ch[dep].time_diff_flag) {
            int t = ch[dep].time_diff_index;

            if (ch[dep].time_diff_sign)
                begin -= t;
            else
                end   -= t;

            for (smp = begin; smp < end; smp++) {
                y  = (1 << 6) +
                     MUL64(ch[dep].weighting[0], master[smp - 1    ]) +
                     MUL64(ch[dep].weighting[1], master[smp        ]) +
                     MUL64(ch[dep].weighting[2], master[smp + 1    ]) +
                     MUL64(ch[dep].weighting[3], master[smp - 1 + t]) +
                     MUL64(ch[dep].weighting[4], master[smp     + t]) +
                     MUL64(ch[dep].weighting[5], master[smp + 1 + t]);

                bd->raw_samples[smp] += y >> 7;
            }
        } else {
            for (smp = begin; smp < end; smp++) {
                y  = (1 << 6) +
                     MUL64(ch[dep].weighting[0], master[smp - 1]) +
                     MUL64(ch[dep].weighting[1], master[smp    ]) +
                     MUL64(ch[dep].weighting[2], master[smp + 1]);

                bd->raw_samples[smp] += y >> 7;
            }
        }

        dep++;
    }

    reverted[c] = 1;
}


/** Reads the frame data.
 */
static int read_frame_data(ALSDecContext *ctx, unsigned int ra_frame)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    AVCodecContext *avctx    = ctx->avctx;
    GetBitContext *gb = &ctx->gb;
    unsigned int div_blocks[32];                ///< block sizes.
    unsigned int b, c;
    unsigned int js_blocks[2];

    uint32_t bs_info = 0;

    // skip the size of the ra unit if present in the frame
    if (sconf->ra_flag == RA_FLAG_FRAMES && ra_frame)
        skip_bits_long(gb, 32);

    if (sconf->mc_coding && sconf->joint_stereo) {
        ctx->js_switch = get_bits1(gb);
        align_get_bits(gb);
    }

    if (!sconf->mc_coding || ctx->js_switch) {
        int independent_bs = !sconf->joint_stereo;

        for (c = 0; c < avctx->channels; c++) {
            js_blocks[0] = 0;
            js_blocks[1] = 0;

            get_block_sizes(ctx, div_blocks, &bs_info);

            // if joint_stereo and block_switching is set, independent decoding
            // is signaled via the first bit of bs_info
            if (sconf->joint_stereo && sconf->block_switching)
                if (bs_info >> 31)
                    independent_bs = 2;

            // if this is the last channel, it has to be decoded independently
            if (c == avctx->channels - 1)
                independent_bs = 1;

            if (independent_bs) {
                if (decode_blocks_ind(ctx, ra_frame, c, div_blocks, js_blocks))
                    return -1;

                if (independent_bs)
                    independent_bs--;
            } else {
                if (decode_blocks(ctx, ra_frame, c, div_blocks, js_blocks))
                    return -1;

                c++;
            }

        // store carryover raw samples
        memmove(ctx->raw_samples[c] - sconf->max_order,
                ctx->raw_samples[c] - sconf->max_order + sconf->frame_length,
                sizeof(*ctx->raw_samples[c]) * sconf->max_order);
        }
    } else { // multi-channel coding
        ALSBlockData   bd;
        int            reverted_channels[avctx->channels];
        unsigned int   offset = 0;

        memset(&bd,                0, sizeof(ALSBlockData));
        memset(&reverted_channels, 0, sizeof(reverted_channels) * avctx->channels);

        bd.ra_block         = ra_frame;
        bd.prev_raw_samples = ctx->prev_raw_samples;

        get_block_sizes(ctx, div_blocks, &bs_info);

        for (b = 0; b < ctx->num_blocks; b++) {
            bd.shift_lsbs   = 0;
            bd.block_length = div_blocks[b];

            for (c = 0; c < avctx->channels; c++) {
                bd.use_ltp     = ctx->use_ltp + c;
                bd.ltp_lag     = ctx->ltp_lag + c;
                bd.ltp_gain    = ctx->ltp_gain[c];
                bd.lpc_cof     = ctx->lpc_cof[c];
                bd.quant_cof   = ctx->quant_cof[c];
                bd.raw_samples = ctx->raw_samples[c] + offset;
                bd.raw_other   = NULL;

                read_block(ctx, &bd);
                read_channel_data(ctx, ctx->chan_data[c], c);
            }

            for (c = 0; c < avctx->channels; c++)
                revert_channel_correlation(ctx, &bd, ctx->chan_data,
                                           reverted_channels, offset, c);

            for (c = 0; c < avctx->channels; c++) {
                bd.use_ltp     = ctx->use_ltp + c;
                bd.ltp_lag     = ctx->ltp_lag + c;
                bd.ltp_gain    = ctx->ltp_gain[c];
                bd.lpc_cof     = ctx->lpc_cof[c];
                bd.quant_cof   = ctx->quant_cof[c];
                bd.raw_samples = ctx->raw_samples[c] + offset;
                decode_block(ctx, &bd);
            }

            memset(&reverted_channels, 0, avctx->channels * sizeof(reverted_channels));
            offset      += div_blocks[b];
            bd.ra_block  = 0;
        }

        // store carryover raw samples
        for (c = 0; c < avctx->channels; c++)
            memmove(ctx->raw_samples[c] - sconf->max_order,
                    ctx->raw_samples[c] - sconf->max_order + sconf->frame_length,
                    sizeof(*ctx->raw_samples[c]) * sconf->max_order);
    }

    // TODO: read_diff_float_data

    return 0;
}


/** Decodes an ALS frame.
 */
static int decode_frame(AVCodecContext *avctx,
                        void *data, int *data_size,
                        AVPacket *avpkt)
{
    ALSDecContext *ctx       = avctx->priv_data;
    ALSSpecificConfig *sconf = &ctx->sconf;
    const uint8_t *buffer    = avpkt->data;
    int buffer_size          = avpkt->size;
    int invalid_frame, size;
    unsigned int c, sample, ra_frame, bytes_read, shift;

    init_get_bits(&ctx->gb, buffer, buffer_size * 8);

    // In the case that the distance between random access frames is set to zero
    // (sconf->ra_distance == 0) no frame is treated as a random access frame.
    // For the first frame, if prediction is used, all samples used from the
    // previous frame are assumed to be zero.
    ra_frame = sconf->ra_distance && !(ctx->frame_id % sconf->ra_distance);

    // the last frame to decode might have a different length
    if (ctx->num_frames && ctx->num_frames - 1 == ctx->frame_id)
        ctx->cur_frame_length = ctx->last_frame_length;

    // decode the frame data
    if ((invalid_frame = read_frame_data(ctx, ra_frame) < 0))
        av_log(ctx->avctx, AV_LOG_WARNING,
               "Reading frame data failed. Skipping RA unit.\n");

    ctx->frame_id++;

    // check for size of decoded data
    size = ctx->cur_frame_length * avctx->channels *
           (av_get_bits_per_sample_format(avctx->sample_fmt) >> 3);

    if (size > *data_size) {
        av_log(avctx, AV_LOG_ERROR, "Decoded data exceeds buffer size.\n");
        return -1;
    }

    *data_size = size;

    // transform decoded frame into output format
    #define INTERLEAVE_OUTPUT(bps)                                 \
    {                                                              \
        int##bps##_t *dest = (int##bps##_t*) data;                 \
        shift = bps - ctx->avctx->bits_per_raw_sample;             \
        for (sample = 0; sample < ctx->cur_frame_length; sample++) \
            for (c = 0; c < avctx->channels; c++)                  \
                *dest++ = ctx->raw_samples[c][sample] << shift;    \
    }

    if (ctx->avctx->bits_per_raw_sample <= 16) {
        INTERLEAVE_OUTPUT(16)
    } else {
        INTERLEAVE_OUTPUT(32)
    }

    bytes_read = invalid_frame ? buffer_size :
                                 (get_bits_count(&ctx->gb) + 7) >> 3;

    return bytes_read;
}


/** Uninitializes the ALS decoder.
 */
static av_cold int decode_end(AVCodecContext *avctx)
{
    ALSDecContext *ctx = avctx->priv_data;

    av_freep(&ctx->sconf.chan_pos);

    av_freep(&ctx->use_ltp);
    av_freep(&ctx->ltp_lag);
    av_freep(&ctx->ltp_gain);
    av_freep(&ctx->ltp_gain_buffer);
    av_freep(&ctx->quant_cof);
    av_freep(&ctx->lpc_cof);
    av_freep(&ctx->quant_cof_buffer);
    av_freep(&ctx->lpc_cof_buffer);
    av_freep(&ctx->prev_raw_samples);
    av_freep(&ctx->raw_samples);
    av_freep(&ctx->raw_buffer);

    if (ctx->sconf.mc_coding) {
        av_freep(&ctx->chan_data);
        av_freep(&ctx->chan_data_buffer);
    }

    return 0;
}


/** Initializes the ALS decoder.
 */
static av_cold int decode_init(AVCodecContext *avctx)
{
    unsigned int c;
    unsigned int channel_size;
    int num_buffers;
    ALSDecContext *ctx = avctx->priv_data;
    ALSSpecificConfig *sconf = &ctx->sconf;
    ctx->avctx = avctx;

    if (!avctx->extradata) {
        av_log(avctx, AV_LOG_ERROR, "Missing required ALS extradata.\n");
        return -1;
    }

    if (read_specific_config(ctx)) {
        av_log(avctx, AV_LOG_ERROR, "Reading ALSSpecificConfig failed.\n");
        decode_end(avctx);
        return -1;
    }

    if (check_specific_config(ctx)) {
        decode_end(avctx);
        return -1;
    }

    if (sconf->floating) {
        avctx->sample_fmt          = SAMPLE_FMT_FLT;
        avctx->bits_per_raw_sample = 32;
    } else {
        avctx->sample_fmt          = sconf->resolution > 1
                                     ? SAMPLE_FMT_S32 : SAMPLE_FMT_S16;
        avctx->bits_per_raw_sample = (sconf->resolution + 1) * 8;
    }

    // set lag value for long-term prediction
    if      (avctx->sample_rate >= 192000)
        ctx->ltp_lag_length = 10;
    else if (avctx->sample_rate >=  96000)
        ctx->ltp_lag_length = 9;
    else
        ctx->ltp_lag_length = 8;

    // allocate quantized parcor coefficient buffer
    num_buffers = sconf->mc_coding ? avctx->channels : 1;

    ctx->quant_cof        = av_malloc(sizeof(*ctx->quant_cof) * num_buffers);
    ctx->lpc_cof          = av_malloc(sizeof(*ctx->lpc_cof)   * num_buffers);
    ctx->quant_cof_buffer = av_malloc(sizeof(*ctx->quant_cof_buffer) *
                                      num_buffers * sconf->max_order);
    ctx->lpc_cof_buffer   = av_malloc(sizeof(*ctx->lpc_cof_buffer) *
                                      num_buffers * sconf->max_order);

    if (!ctx->quant_cof        || !ctx->lpc_cof       ||
        !ctx->quant_cof_buffer || !ctx->lpc_cof_buffer) {
        av_log(avctx, AV_LOG_ERROR, "Allocating buffer memory failed.\n");
        return AVERROR(ENOMEM);
    }

    // assign quantized parcor coefficient buffers
    for (c = 0; c < num_buffers; c++) {
        ctx->quant_cof[c] = ctx->quant_cof_buffer + c * sconf->max_order;
        ctx->lpc_cof[c]   = ctx->lpc_cof_buffer   + c * sconf->max_order;
    }

    // allocate and assign lag and gain data buffer for ltp mode
    ctx->use_ltp         = av_mallocz(sizeof(*ctx->use_ltp)  * num_buffers);
    ctx->ltp_lag         = av_malloc (sizeof(*ctx->ltp_lag)  * num_buffers);
    ctx->ltp_gain        = av_malloc (sizeof(*ctx->ltp_gain) * num_buffers);
    ctx->ltp_gain_buffer = av_malloc (sizeof(*ctx->ltp_gain_buffer) *
                                      num_buffers * 5);

    if (!ctx->use_ltp  || !ctx->ltp_lag ||
        !ctx->ltp_gain || !ctx->ltp_gain_buffer) {
        av_log(avctx, AV_LOG_ERROR, "Allocating buffer memory failed.\n");
        decode_end(avctx);
        return AVERROR(ENOMEM);
    }

    for (c = 0; c < num_buffers; c++)
        ctx->ltp_gain[c] = ctx->ltp_gain_buffer + c * 5;

    // allocate and assign channel data buffer for mcc mode
    if (sconf->mc_coding) {
        ctx->chan_data_buffer = av_malloc(sizeof(*ctx->chan_data_buffer) *
                                          num_buffers);
        ctx->chan_data        = av_malloc(sizeof(ALSChannelData) *
                                          num_buffers);

        if (!ctx->chan_data_buffer || !ctx->chan_data) {
            av_log(avctx, AV_LOG_ERROR, "Allocating buffer memory failed.\n");
            decode_end(avctx);
            return AVERROR(ENOMEM);
        }

        for (c = 0; c < num_buffers; c++)
            ctx->chan_data[c] = ctx->chan_data_buffer + c;
    } else {
        ctx->chan_data        = NULL;
        ctx->chan_data_buffer = NULL;
    }

    avctx->frame_size = sconf->frame_length;
    channel_size      = sconf->frame_length + sconf->max_order;

    ctx->prev_raw_samples = av_malloc (sizeof(*ctx->prev_raw_samples) * sconf->max_order);
    ctx->raw_buffer       = av_mallocz(sizeof(*ctx->     raw_buffer)  * avctx->channels * channel_size);
    ctx->raw_samples      = av_malloc (sizeof(*ctx->     raw_samples) * avctx->channels);

    // allocate previous raw sample buffer
    if (!ctx->prev_raw_samples || !ctx->raw_buffer|| !ctx->raw_samples) {
        av_log(avctx, AV_LOG_ERROR, "Allocating buffer memory failed.\n");
        decode_end(avctx);
        return AVERROR(ENOMEM);
    }

    // assign raw samples buffers
    ctx->raw_samples[0] = ctx->raw_buffer + sconf->max_order;
    for (c = 1; c < avctx->channels; c++)
        ctx->raw_samples[c] = ctx->raw_samples[c - 1] + channel_size;

    return 0;
}


/** Flushes (resets) the frame ID after seeking.
 */
static av_cold void flush(AVCodecContext *avctx)
{
    ALSDecContext *ctx = avctx->priv_data;

    ctx->frame_id = 0;
}


AVCodec als_decoder = {
    "als",
    CODEC_TYPE_AUDIO,
    CODEC_ID_MP4ALS,
    sizeof(ALSDecContext),
    decode_init,
    NULL,
    decode_end,
    decode_frame,
    .flush = flush,
    .capabilities = CODEC_CAP_SUBFRAMES,
    .long_name = NULL_IF_CONFIG_SMALL("MPEG-4 Audio Lossless Coding (ALS)"),
};

