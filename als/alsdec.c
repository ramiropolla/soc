/*
 * ALS decoder
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
 * @file libavcodec/alsdec.c
 * MPEG-4 ALS decoder
 * @author Thilo Borgmann <thilo.borgmann _at_ googlemail.com>
 */


//#define DEBUG


#include "avcodec.h"
#include "get_bits.h"


typedef struct {
    uint32_t als_id;                 ///< ALS identifier
    uint32_t samp_freq;              ///< Sampling frequency in Hz
    uint32_t samples;                ///< Number of samples (per channel). = 0xFFFFFFFF if unknown.
    int      channels;               ///< Number of channels
    int      file_type;              ///< Not used. Provided for debugging.
    int      resolution;             ///< 000 = 8-bit; 001 = 16-bit; 010 = 24-bit; 011 = 32-bit
    int      floating;               ///< 1 = IEEE 32-bit floating-point, 0 = integer
    int      msb_first;              ///< Original byte order of the input audio data
    int      frame_length;           ///< Frame Length
    int      random_access;          ///< Distance between RA frames (in frames, 0...255)
    int      ra_flag;                ///< Indicates where the size of ra units is stored.
    int      adapt_order;            ///< Adaptive order: 1 = on, 0 = off
    int      coef_table;             ///< Table index of rice code parameters
    int      long_term_prediction;   ///< Long term prediction (LTP): 1 = on, 0 = off
    int      max_order;              ///< Maximum prediction order (0..1023)
    int      block_switching;        ///< Number of block switching levels
    int      bgmc_mode;              ///< BGMC Mode: 1 = on, 0 = off (Rice coding only)
    int      sb_part;                ///< Sub-block partition
    int      joint_stereo;           ///< Joint Stereo: 1 = on, 0 = off
    int      mc_coding;              ///< Extended inter-channel coding: 1 = on, 0 = off
    int      chan_config;            ///< Indicates that a chan_config_info field is present
    int      chan_sort;              ///< Channel rearrangement: 1 = on, 0 = off
    int      crc_enabled;            ///< Indicates that the crc field is present
    int      RLSLMS;                 ///< Use RLS-LMS predictor: 1 = on, 0 = off
    int      aux_data_enabled;       ///< Indicates that auxiliary data is present
    int      chan_config_info;       ///< Mapping of channels to loudspeaker locations
    int      *chan_pos;              ///< Original channel positions
    uint32_t header_size;            ///< Header size of original audio file in bytes. Provided for debugging.
    uint32_t trailer_size;           ///< Trailer size of original audio file in bytes. Provided for debugging.
    uint32_t crc;                    ///< 32-bit CCITT-32 CRC checksum
} ALSSpecificConfig;


typedef struct {
    AVCodecContext    *avctx;
    ALSSpecificConfig sconf;
    GetBitContext     gb;                ///< A bit reader context.
    unsigned int      num_frames;        ///< Number of frames to decode. 0 if unknown.
    unsigned int      last_frame_length; ///< Length of the last frame to decode. 0 if unknown.
    unsigned int      frame_id;          ///< The frame id / number of the current frame.
    unsigned int      js_switch;         ///< If true, joint-stereo decoding is enforced.
} ALSDecContext;


#ifdef DEBUG
static av_cold void dprint_specific_config(ALSDecContext *ctx)
{
    AVCodecContext *avctx    = ctx->avctx;
    ALSSpecificConfig *sconf = &ctx->sconf;

    dprintf(avctx, "als_id = %x\n",               sconf->als_id);
    dprintf(avctx, "samp_freq = %i\n",            sconf->samp_freq);
    dprintf(avctx, "samples = %i\n",              sconf->samples);
    dprintf(avctx, "channels = %i\n",             sconf->channels);
    dprintf(avctx, "file_type = %i\n",            sconf->file_type);
    dprintf(avctx, "resolution = %i\n",           sconf->resolution);
    dprintf(avctx, "floating = %i\n",             sconf->floating);
    dprintf(avctx, "msb_first = %i\n",            sconf->msb_first);
    dprintf(avctx, "frame_length = %i\n",         sconf->frame_length);
    dprintf(avctx, "ra_flag = %i\n",              sconf->ra_flag);
    dprintf(avctx, "adapt_order = %i\n",          sconf->adapt_order);
    dprintf(avctx, "coef_table = %i\n",           sconf->coef_table);
    dprintf(avctx, "long_term_prediction = %i\n", sconf->long_term_prediction);
    dprintf(avctx, "max_order = %i\n",            sconf->max_order);
    dprintf(avctx, "block_switching = %i\n",      sconf->block_switching);
    dprintf(avctx, "bgmc_mode = %i\n",            sconf->bgmc_mode);
    dprintf(avctx, "sb_part = %i\n",              sconf->sb_part);
    dprintf(avctx, "joint_stereo = %i\n",         sconf->joint_stereo);
    dprintf(avctx, "mc_coding = %i\n",            sconf->mc_coding);
    dprintf(avctx, "chan_config = %i\n",          sconf->chan_config);
    dprintf(avctx, "chan_sort = %i\n",            sconf->chan_sort);
    dprintf(avctx, "crc_enabled = %i\n",          sconf->crc_enabled);
    dprintf(avctx, "RLSLMS = %i\n",               sconf->RLSLMS);
    dprintf(avctx, "aux_data_enabled = %i\n",     sconf->aux_data_enabled);
    dprintf(avctx, "chan_config_info = %i\n",     sconf->chan_config_info);
    dprintf(avctx, "header_size = %i\n",          sconf->header_size);
    dprintf(avctx, "trailer_size = %i\n",         sconf->trailer_size);
    dprintf(avctx, "crc_enabled = %i\n",          sconf->crc_enabled);

    dprintf(avctx, " num_frames = %i\n",          ctx->num_frames);
    dprintf(avctx, " last_frame_length = %i\n",   ctx->last_frame_length);
}
#endif


/** Reads an ALSSpecificConfig from a buffer into the output struct.
 */
static av_cold int read_specific_config(ALSDecContext *ctx,
                                        const uint8_t *buffer, int buffer_size,
                                        ALSSpecificConfig *sconf)
{
    GetBitContext gb;
    uint64_t ht_size;
    int i;

    if (buffer_size < 22)
        return -1;

    init_get_bits(&gb, buffer, buffer_size * 8);

    // read the fixed items
    sconf->als_id               = get_bits_long(&gb, 32);
    sconf->samp_freq            = get_bits_long(&gb, 32);
    sconf->samples              = get_bits_long(&gb, 32);
    sconf->channels             = get_bits(&gb, 16) + 1;
    sconf->file_type            = get_bits(&gb, 3);
    sconf->resolution           = get_bits(&gb, 3);
    sconf->floating             = get_bits1(&gb);
    sconf->msb_first            = get_bits1(&gb);
    sconf->frame_length         = get_bits(&gb, 16) + 1;
    sconf->random_access        = get_bits(&gb, 8);
    sconf->ra_flag              = get_bits(&gb, 2);
    sconf->adapt_order          = get_bits1(&gb);
    sconf->coef_table           = get_bits(&gb, 2);
    sconf->long_term_prediction = get_bits1(&gb);
    sconf->max_order            = get_bits(&gb, 10);
    sconf->block_switching      = get_bits(&gb, 2);
    sconf->bgmc_mode            = get_bits1(&gb);
    sconf->sb_part              = get_bits1(&gb);
    sconf->joint_stereo         = get_bits1(&gb);
    sconf->mc_coding            = get_bits1(&gb);
    sconf->chan_config          = get_bits1(&gb);
    sconf->chan_sort            = get_bits1(&gb);
    sconf->crc_enabled          = get_bits1(&gb);
    sconf->RLSLMS               = get_bits1(&gb);
    skip_bits(&gb, 5);                                      // skip 5 reserved bits
    sconf->aux_data_enabled     = get_bits1(&gb);
    buffer_size -= 22;


    // check for ALSSpecificConfig struct
    if (sconf->als_id != MKBETAG('A','L','S','\0'))
        return -1;


    // calculate total number of frames to decode if possible
    if (sconf->samples != 0xFFFFFFFF) {
        ctx->num_frames        = ((sconf->samples - 1) / sconf->frame_length) + 1;
        ctx->last_frame_length = sconf->samples % ctx->sconf.frame_length;
        if (!ctx->last_frame_length) {
            ctx->last_frame_length = sconf->frame_length;
        }
    } else {
        ctx->num_frames        = 0;
        ctx->last_frame_length = 0;
    }


    // read channel config
    if (sconf->chan_config) {
        if (buffer_size < 2)
            return -1;

        sconf->chan_config_info = get_bits(&gb, 16);
        buffer_size -= 2;
        // TODO: use this to set avctx->channel_layout
    }


    // read channel sorting
    if (sconf->chan_sort && sconf->channels > 1) {
        int chan_pos_bits = av_log2((sconf->channels - 1) << 1);
        int bytes_needed  = (sconf->channels * chan_pos_bits + 7) / 8;
        if (buffer_size < bytes_needed)
            return -1;

        if(!(sconf->chan_pos = av_malloc(sconf->channels * sizeof(int))))
            return -1;

        for (i = 0; i < sconf->channels; i++) {
            sconf->chan_pos[i] = get_bits(&gb, chan_pos_bits);
        }

        align_get_bits(&gb);
        buffer_size -= bytes_needed;
    } else {
        sconf->chan_sort = 0;
    }


    // read fixed header and trailer sizes, if size = 0xFFFFFFFF then there is no data field!
    if (buffer_size < 8)
        return -1;

    sconf->header_size  = get_bits_long(&gb, 32);
    sconf->trailer_size = get_bits_long(&gb, 32);
    if (sconf->header_size  == 0xFFFFFFFF)
        sconf->header_size  = 0;
    if (sconf->trailer_size == 0xFFFFFFFF)
        sconf->trailer_size = 0;

    ht_size = sconf->header_size + sconf->trailer_size;

    buffer_size -= 8;


    // skip the header and trailer data
    if (buffer_size < ht_size)
        return -1;

    ht_size *= 8;

    while (ht_size > 0) {
        int len = FFMIN(ht_size, INT32_MAX);
        skip_bits_long(&gb, len);
        ht_size -= len;
    }

    buffer_size -= ht_size;


    // read the crc data
    if (sconf->crc_enabled) {
        if (buffer_size < 4)
            return -1;

        sconf->crc = get_bits_long(&gb, 32);
    }


    // no need to read the rest of ALSSpecificConfig (ra_unit_size & aux data)
#ifdef DEBUG
    dprint_specific_config(ctx);
#endif
    return 0;
}


/** Parses the bs_info item to extract the block partitioning.
 */
static void parse_bs_info(uint32_t bs_info, unsigned int n, unsigned int div,
                          unsigned int **div_blocks, unsigned int *num_blocks)
{
    if (n < 32 && ((bs_info >> (30 - n)) & 1)) {
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


/** Reads the frame data.
 */
static int read_frame_data(ALSDecContext *ctx)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    GetBitContext *gb = &ctx->gb;
    int c, b;
    uint32_t bs_info = 0;
    unsigned int num_blocks;
    unsigned int div_blocks[32];
    unsigned int *ptr_div_blocks = &div_blocks[0];

    // skip ra_unit_size if present
    if (sconf->ra_flag == 1 && !(ctx->frame_id % sconf->random_access))
        skip_bits_long(gb, 32);

    if (sconf->mc_coding && sconf->joint_stereo) {
        ctx->js_switch = get_bits1(gb);
        align_get_bits(gb);
    }

    if (!sconf->mc_coding || ctx->js_switch) {
        int independent_bs = !sconf->joint_stereo;

        for (c = 0; c < sconf->channels; c++) {
            if (sconf->block_switching) {
                unsigned int bs_info_len = 1 << (sconf->block_switching + 2);
                bs_info = get_bits_long(gb, bs_info_len);
                bs_info <<= (32 - bs_info_len);
            }

            num_blocks = 0;
            parse_bs_info(bs_info, 0, 0, &ptr_div_blocks, &num_blocks);
#ifdef DEBUG
            dprintf(ctx->avctx, "bs_info = %x, block sizes:", bs_info);
            for (b = 0; b < num_blocks; b++)
                dprintf(ctx->avctx, " %i", div_blocks[b]);
            dprintf(ctx->avctx, "\n");
#endif
            // if this is the last channel, it has to be decoded independently
            if (c == sconf->channels - 1)
                independent_bs = 1;

            // if joint_stereo and block_switching is set, independent decoding
            // is signaled via the first bit of bs_info
            if(sconf->joint_stereo && sconf->block_switching)
                independent_bs = bs_info >> 31;

            if (independent_bs) {
                for (b = 0; b < num_blocks; b++) {
                    dprintf(ctx->avctx, "reading block A\n");
                    read_block_data(ctx);
                }
            } else {
                for (b = 0; b < num_blocks; b++) {
                    dprintf(ctx->avctx, "reading block B\n");
                    read_block_data(ctx);
                    dprintf(ctx->avctx, "reading block C\n");
                    read_block_data(ctx);
                }
                c++;
            }
        }
    } else {
        if (sconf->block_switching) {
            unsigned int bs_info_len = 1 << (sconf->block_switching + 2);
            bs_info = get_bits_long(gb, bs_info_len);
            bs_info <<= (32 - bs_info_len);
        }

        num_blocks = 0;
        parse_bs_info(bs_info, 0, 0, &ptr_div_blocks, &num_blocks);
#ifdef DEBUG
            dprintf(ctx->avctx, "bs_info = %x, block sizes:", bs_info);
            for (b = 0; b < num_blocks; b++)
                dprintf(ctx->avctx, " %i", div_blocks[b]);
            dprintf(ctx->avctx, "\n");
#endif

        for (b = 0; b < num_blocks; b++) {
            dprintf(ctx->avctx, "reading block D\n");
            read_block_data(ctx);
            // TODO: read_channel_data
        }
    }

    if (sconf->floating) {
        unsigned int num_bytes_diff_float = get_bits_long(gb, 32);
        // TODO: read_diff_float_data
    }

    return 0;
}


/** Decodes an ALS frame.
 */
static int decode_frame(AVCodecContext *avctx,
                        void *data, int *data_size,
                        AVPacket *avpkt)
{
    ALSDecContext *ctx       = avctx->priv_data;
    const uint8_t *buffer    = avpkt->data;
    int buffer_size          = avpkt->size;

    init_get_bits(&ctx->gb, buffer, buffer_size * 8);

    if(read_frame_data(ctx)) {
        av_log(ctx->avctx, AV_LOG_ERROR, "Reading frame data failed.\n");
        return -1;
    }

    return 0;
}


/** Initializes the ALS decoder.
 */
static av_cold int decode_init(AVCodecContext *avctx)
{
    ALSDecContext *ctx = avctx->priv_data;
    ctx->avctx = avctx;

    if (!avctx->extradata) {
        av_log(avctx, AV_LOG_ERROR, "Missing required ALS extradata.\n");
        return -1;
    }

    if (read_specific_config(ctx, avctx->extradata + 6,
                             avctx->extradata_size - 6, &ctx->sconf)) {
        av_log(avctx, AV_LOG_ERROR, "Reading ALSSpecificConfig failed.\n");
        return -1;
    }

    avctx->sample_rate = ctx->sconf.samp_freq;
    avctx->channels    = ctx->sconf.channels;

    if (ctx->sconf.floating) {
        avctx->sample_fmt          = SAMPLE_FMT_FLT;
        avctx->bits_per_raw_sample = 32;
    } else {
        avctx->sample_fmt          = ctx->sconf.resolution > 1 ? SAMPLE_FMT_S32 : SAMPLE_FMT_S16;
        avctx->bits_per_raw_sample = (ctx->sconf.resolution + 1) * 8;
    }

    return 0;
}


/** Uninitializes the ALS decoder.
 */
static av_cold int decode_end(AVCodecContext *avctx)
{
    ALSDecContext *ctx = avctx->priv_data;

    av_freep(&ctx->sconf.chan_pos);

    return 0;
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
    .long_name = NULL_IF_CONFIG_SMALL("MPEG-4 Audio Lossless Coding"),
};

