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
#include "unary.h"

#include "als_data.h"

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
    int      coef_table;             ///< Table index of Rice code parameters
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
    unsigned int      num_blocks;        ///< Number of blocks used in the current frame.
    int64_t           *residuals;        ///< Decoded residuals of the current block.
    int64_t           **raw_samples;     ///< Decoded raw samples for each channel.
    int64_t           *raw_buffer;       ///< Contains all decoded raw samples including carryover samples.
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


/** Computes ceil(log2(x)) using av_log2.
 */
static inline int ceil_log2(int x) {
    if (x <= 0)
        return 0;
    return av_log2((x - 1) << 1);
}


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
        int chan_pos_bits = ceil_log2(sconf->channels);
        int bytes_needed  = (sconf->channels * chan_pos_bits + 7) >> 3;
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


/** Reads and decodes a Rice codeword.
 */
static int64_t decode_rice(GetBitContext *gb, unsigned int k)
{
    int64_t value = 0;
    int64_t q = 0;
    int     max = gb->size_in_bits - get_bits_count(gb) - k;

    if (!k) {
        q = get_unary(gb, 0, max);
        return (q & 1) ? -((q + 1) >> 1) : ((q + 1) >> 1);
    } else if (k == 1) {
        q = get_unary(gb, 0, max);
        return get_bits1(gb) ? q : -(q + 1);
    } else {
        unsigned int r, sub_sign;

        q         = get_unary(gb, 0, max);
        sub_sign  = get_bits1(gb);
        r         = get_bits_long(gb, k - 1);

        value = (q << (k - 1)) + r;

        return sub_sign ? value : -(value + 1);
    }
}


/** Converts PARCOR coefficient k to direct filter coefficient.
 */
static void parcor_to_lpc(unsigned int k, int64_t *par, int64_t *cof)
{
    int i;
    int64_t tmp1, tmp2;

    for (i = 0; i < ((k+1) >> 1); i++) {
        tmp1 = cof[    i    ] + ((par[k] * cof[k - i - 1] + (1 << 19)) >> 20);
        tmp2 = cof[k - i - 1] + ((par[k] * cof[    i    ] + (1 << 19)) >> 20);
        cof[k - i - 1] = tmp2;
        cof[    i    ] = tmp1;
    }

    cof[k] = par[k];
}


/** Converts all PARCOR coefficients to direct filter coefficients.
 */
static void all_parcor_to_lpc(unsigned int num, int64_t *par, int64_t *cof)
{
    int k;

    for (k = 0; k < num; k++)
        parcor_to_lpc(k, par, cof);
}


/** Reads the block data.
 */
static int read_block_data(ALSDecContext *ctx, unsigned int ra_block,
                            int64_t *raw_samples, unsigned int block_length,
                            uint32_t *js_blocks, int64_t *raw_other)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    AVCodecContext *avctx    = ctx->avctx;
    GetBitContext *gb        = &ctx->gb;
    unsigned int block_type;
    unsigned int k;

    *js_blocks <<= 1;

    // read block type
    block_type = get_bits1(gb);

    if (block_type == 0) {
        unsigned int const_block;
        int32_t      const_val = 0;

        const_block  = get_bits1(gb);    // 1 = constant value, 0 = zero block (silence)
        *js_blocks  |= get_bits1(gb);

        // skip 5 reserved bits
        skip_bits(gb, 5);

        if (const_block) {
            unsigned int const_val_bits;

            if (sconf->resolution == 2 || sconf->floating) {
                const_val_bits = 24;
            } else {
                const_val_bits = avctx->bits_per_raw_sample;
            }

            const_val = get_bits_long(gb, const_val_bits);
        }

        // write raw samples into buffer
        for (k = 0; k < block_length; k++) {
            raw_samples[k] = const_val;
        }
    } else {
        unsigned int s[8];
        unsigned int sub_blocks, sb_length, shift_lsbs;
        unsigned int opt_order = 1;
        int64_t      quant_cof[sconf->max_order];
        int64_t      lpc_cof[sconf->max_order];
        unsigned int start = 0;
        int64_t      *res = ctx->residuals;
        unsigned int sb, smp;
        int64_t      y;

        *js_blocks |= get_bits1(gb);

        // Determine the number of sub blocks for entropy decoding
        if (!sconf->bgmc_mode && !sconf->sb_part)
            sub_blocks = 1;
        else if (sconf->bgmc_mode && sconf->sb_part)
            sub_blocks = 1 << get_bits(gb, 2);
        else
            sub_blocks = get_bits1(gb) ? 4 : 1;

        // Do not continue in case of a damaged stream
        if (block_length % sub_blocks)
            return -1;

        sb_length = block_length / sub_blocks;


        if (!sconf->bgmc_mode) {
            s[0] = get_bits(gb, (sconf->resolution > 1) ? 5 : 4);
            for (k = 1; k < sub_blocks; k++)
                s[k] = s[k - 1] + decode_rice(gb, 0);
        } else {
            // TODO: BGMC mode
        }

        shift_lsbs = get_bits1(gb);

        if (shift_lsbs) {
            // TODO: LSBS shifts
        }


        if (!sconf->RLSLMS) {
            int64_t quant_index;

            if (sconf->adapt_order) {
                int opt_order_length =
                        FFMIN(
                        ceil_log2(sconf->max_order+1),
                        FFMAX(ceil_log2((block_length >> 3) - 1), 1)
                        );
                opt_order = get_bits(gb, opt_order_length);
            } else {
                opt_order = sconf->max_order;
            }

            if (opt_order) {
                if (sconf->coef_table == 3) {
                    // read coefficient 0
                    quant_index = get_bits(gb, 7) - 64;
                    quant_cof[0] = parcor_scaled_values[quant_index + 64];

                    // read coefficient 1
                    quant_index = get_bits(gb, 7) - 64;
                    quant_cof[1] = -parcor_scaled_values[quant_index + 64];

                    // read coefficients 2 - opt_order
                    for (k = 2; k < opt_order; k++) {
                        quant_index = get_bits(gb, 7) - 64;
                        quant_cof[k] = (quant_index << 14) + (1 << 13);
                    }
                } else {
                    int offset, rice_param, k_max;

                    // read coefficient 0
                    offset       = parcor_rice_table[sconf->coef_table][0][0];
                    rice_param   = parcor_rice_table[sconf->coef_table][0][1];
                    quant_index  = decode_rice(gb, rice_param) + offset;
                    quant_cof[0] = parcor_scaled_values[quant_index + 64];

                    // read coefficient 1
                    offset       = parcor_rice_table[sconf->coef_table][1][0];
                    rice_param   = parcor_rice_table[sconf->coef_table][1][1];
                    quant_index  = decode_rice(gb, rice_param) + offset;
                    quant_cof[1] = -parcor_scaled_values[quant_index + 64];

                    // read coefficients 2 - 19
                    k_max = FFMIN(20, opt_order);
                    for (k = 2; k < k_max; k++) {
                        offset       = parcor_rice_table[sconf->coef_table][k][0];
                        rice_param   = parcor_rice_table[sconf->coef_table][k][1];
                        quant_index  = decode_rice(gb, rice_param) + offset;
                        quant_cof[k] = (quant_index << 14) + (1 << 13);
                    }

                    // read coefficients 20 - 126
                    k_max = FFMIN(127, opt_order);
                    for (k = 20; k < k_max; k++) {
                        offset       = k & 1;
                        rice_param   = 2;
                        quant_index  = decode_rice(gb, rice_param) + offset;
                        quant_cof[k] = (quant_index << 14) + (1 << 13);
                    }

                    // read coefficients 127 - opt_order
                    for (k = 127; k < opt_order; k++) {
                        offset       = 0;
                        rice_param   = 1;
                        quant_index  = decode_rice(gb, rice_param) + offset;
                        quant_cof[k] = (quant_index << 14) + (1 << 13);
                    }
                }
            }
        }

        if (sconf->long_term_prediction) {
            // TODO: LTP mode
        }

        start = 0;

        // read first value and residuals in case of a random access block
        if (ra_block) {
            if (opt_order) {
                raw_samples[0] = decode_rice(gb, avctx->bits_per_raw_sample - 4);
                res[0] = raw_samples[0];
            }
            if (opt_order > 1)
                res[1] = decode_rice(gb, s[0] + 3);
            if (opt_order > 2)
                res[2] = decode_rice(gb, s[0] + 1);

            start = FFMIN(opt_order, 3);
        } else {
            // TODO: check if this has to be a function after features are implemented.
            all_parcor_to_lpc(opt_order, quant_cof, lpc_cof);
        }

        // read all residuals
        // TODO: decode directly into ctx->raw_samples[] instead of storing the residuals
        if (sconf->bgmc_mode) {
            // TODO: BGMC mode
        } else {
            int64_t *current_res = res;

            for (sb = 0; sb < sub_blocks; sb++) {
                for (k = start; k < sb_length; k++) {
                    current_res[k] = decode_rice(gb, s[sb]);
                }
                current_res += sb_length;
                start = 0;
            }
         }

        // reconstruct all samples from residuals
        if (ra_block) {
            unsigned int progressive = FFMIN(block_length, opt_order);

            for (smp = 0; smp < progressive; smp++) {
                y = 1 << 19;

                for (sb = 0; sb < smp; sb++)
                    y += lpc_cof[sb] * raw_samples[smp - (sb + 1)];

                raw_samples[smp] = res[smp] - (y >> 20);
                parcor_to_lpc(smp, quant_cof, lpc_cof);
            }

            for (; smp < block_length; smp++) {
                y = 1 << 19;

                for (sb = 0; sb < progressive; sb++)
                    y += lpc_cof[sb] * raw_samples[smp - (sb + 1)];

                raw_samples[smp] = res[smp] - (y >> 20);
            }
        } else {
            // reconstruct difference signal for prediction (joint-stereo)
            if (*js_blocks & 1 && raw_other) {
                int i;
                if (raw_other > raw_samples) {          // L = R - D
                    for (i = -1; i >= -sconf->max_order; i--)
                        raw_samples[i] = raw_other[i] - raw_samples[i];
                } else {                                // R = D + L
                    for (i = -1; i >= -sconf->max_order; i--)
                        raw_samples[i] = raw_samples[i] - raw_other[i];
                }
            }

            // reconstruct raw samples
            for (smp = 0; smp < block_length; smp++) {
                y = 1 << 19;

                for (sb = 0; sb < opt_order; sb++)
                    y += lpc_cof[sb] * raw_samples[smp - (sb + 1)];

                raw_samples[smp] = res[smp] - (y >> 20);
            }
        }
    }

    if (sconf->RLSLMS) {
        // TODO: read RLSLMS extension data
    }

    if (!sconf->mc_coding || ctx->js_switch) {
        align_get_bits(gb);
    }

    return 0;
}


/** Reads the frame data.
 */
static int read_frame_data(ALSDecContext *ctx, unsigned int ra_frame)
{
    ALSSpecificConfig *sconf = &ctx->sconf;
    GetBitContext *gb = &ctx->gb;
    unsigned int div_blocks[32];                ///< Block sizes.
    unsigned int c, b, ra_block, block_length;
    int64_t *raw_samples;
    uint32_t js_blocks[2];

    uint32_t bs_info = 0;
    unsigned int *ptr_div_blocks = &div_blocks[0];

    // skip ra_unit_size if present
    if (sconf->ra_flag == 1 && ra_frame)
        skip_bits_long(gb, 32);

    if (sconf->mc_coding && sconf->joint_stereo) {
        ctx->js_switch = get_bits1(gb);
        align_get_bits(gb);
    }

    if (!sconf->mc_coding || ctx->js_switch) {
        int independent_bs = !sconf->joint_stereo;

        for (c = 0; c < sconf->channels; c++) {
            js_blocks[0] = 0;
            js_blocks[1] = 0;

            if (sconf->block_switching) {
                unsigned int bs_info_len = 1 << (sconf->block_switching + 2);
                bs_info = get_bits_long(gb, bs_info_len);
                bs_info <<= (32 - bs_info_len);
            }

            ctx->num_blocks = 0;
            parse_bs_info(bs_info, 0, 0, &ptr_div_blocks, &ctx->num_blocks);

            // if this is the last channel, it has to be decoded independently
            if (c == sconf->channels - 1)
                independent_bs = 1;

            // if joint_stereo and block_switching is set, independent decoding
            // is signaled via the first bit of bs_info
            if(sconf->joint_stereo && sconf->block_switching)
                independent_bs = bs_info >> 31;

            if (independent_bs) {
                raw_samples = ctx->raw_samples[c];

                for (b = 0; b < ctx->num_blocks; b++) {
                    ra_block = !b && ra_frame;
                    block_length = sconf->frame_length >> div_blocks[b];
                    read_block_data(ctx, ra_block, raw_samples, block_length,
                                    &js_blocks[0], NULL);
                    raw_samples += block_length;
                }

                memmove((ctx->raw_samples[c]) - sconf->max_order,
                        (ctx->raw_samples[c]) - sconf->max_order + sconf->frame_length,
                        sizeof(int64_t) * sconf->max_order);
            } else {
                unsigned int offset = 0;

                for (b = 0; b < ctx->num_blocks; b++) {
                    ra_block = !b && ra_frame;
                    block_length = sconf->frame_length >> div_blocks[b];

                    raw_samples = ctx->raw_samples[c] + offset;
                    read_block_data(ctx, ra_block, raw_samples, block_length,
                                    &js_blocks[0], ctx->raw_samples[c + 1] + offset);

                    raw_samples = ctx->raw_samples[c + 1] + offset;
                    read_block_data(ctx, ra_block, raw_samples, block_length,
                                    &js_blocks[1], ctx->raw_samples[c] + offset);

                    offset += block_length;
                }

                if (js_blocks[0] || js_blocks[1]) {
                    int64_t *raw_samples_L, *raw_samples_R;
                    unsigned int s;
                    b = ctx->num_blocks - 1;

                    raw_samples_L = ctx->raw_samples[c    ] + sconf->frame_length;
                    raw_samples_R = ctx->raw_samples[c + 1] + sconf->frame_length;

                    while (js_blocks[0] || js_blocks[1]) {
                        unsigned int diff_l, diff_r;
                        block_length   = sconf->frame_length >> div_blocks[b];
                        raw_samples_L -= block_length;
                        raw_samples_R -= block_length;

                        diff_l = js_blocks[0] & 1;
                        diff_r = js_blocks[1] & 1;

                        if (diff_l) {                     // L = R - D
                            if (diff_r)
                                av_log(ctx->avctx, AV_LOG_WARNING, "Invalid channel pair!");

                            for (s = 0; s < block_length; s++)
                                raw_samples_L[s] = raw_samples_R[s] - raw_samples_L[s];
                        } else if (diff_r) {                // R = D + L
                            for (s = 0; s < block_length; s++)
                                raw_samples_R[s] = raw_samples_R[s] + raw_samples_L[s];
                        }

                        b--;
                        js_blocks[0] >>= 1;
                        js_blocks[1] >>= 1;
                    }
                }

                // store carryover raw samples
                memmove((ctx->raw_samples[c]) - sconf->max_order,
                        (ctx->raw_samples[c]) - sconf->max_order + sconf->frame_length,
                        sizeof(int64_t) * sconf->max_order);

                memmove((ctx->raw_samples[c + 1]) - sconf->max_order,
                        (ctx->raw_samples[c + 1]) - sconf->max_order + sconf->frame_length,
                        sizeof(int64_t) * sconf->max_order);

                c++;
            }
        }
    } else {
        if (sconf->block_switching) {
            unsigned int bs_info_len = 1 << (sconf->block_switching + 2);
            bs_info = get_bits_long(gb, bs_info_len);
            bs_info <<= (32 - bs_info_len);
        }

        ctx->num_blocks = 0;
        parse_bs_info(bs_info, 0, 0, &ptr_div_blocks, &ctx->num_blocks);

        // TODO: multi channel coding might use a temporary buffer instead as
        //       the actual channel is not known when read_block-data is called
        raw_samples = ctx->raw_samples[0];

        for (b = 0; b < ctx->num_blocks; b++) {
            ra_block = !b && ra_frame;
            block_length = sconf->frame_length >> div_blocks[b];
            read_block_data(ctx, ra_block, raw_samples, block_length,
                            &js_blocks[0], NULL);
            raw_samples += block_length;
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
    ALSSpecificConfig *sconf = &ctx->sconf;
    const uint8_t *buffer    = avpkt->data;
    int buffer_size          = avpkt->size;
    int16_t *dest            = (int16_t*)data;
    unsigned int c, sample, ra_frame, bytes_read;

    init_get_bits(&ctx->gb, buffer, buffer_size * 8);
    ra_frame = !ctx->frame_id ||
               (sconf->random_access && !(ctx->frame_id % sconf->random_access));

    // the last frame to decode might have a different length
    if (ctx->num_frames && ctx->num_frames - 1 == ctx->frame_id) {
        sconf->frame_length = ctx->last_frame_length;
    }

    // decode the frame data
    if (read_frame_data(ctx, ra_frame)) {
        av_log(ctx->avctx, AV_LOG_ERROR, "Reading frame data failed.\n");
        return -1;
    }

    // increment the frame counter
    ctx->frame_id++;

    // transform decoded frame into output format
    // TODO: Support other resolutions than 16 bit
    for (sample = 0; sample < sconf->frame_length; sample++) {
        for (c = 0; c < sconf->channels; c++) {
            *(dest++) = (int16_t) (ctx->raw_samples[c][sample]);
        }
    }

    *data_size = sconf->frame_length * sconf->channels
                                     * (avctx->sample_fmt == SAMPLE_FMT_S16 ?
                                        2 : 4);

    bytes_read = (get_bits_count(&ctx->gb) + 7) >> 3;

    return bytes_read;
}


/** Uninitializes the ALS decoder.
 */
static av_cold int decode_end(AVCodecContext *avctx)
{
    ALSDecContext *ctx       = avctx->priv_data;

    av_freep(&ctx->sconf.chan_pos);
    av_freep(&ctx->residuals);

    av_freep(&ctx->raw_samples);
    av_freep(&ctx->raw_buffer);

    return 0;
}


/** Initializes the ALS decoder.
 */
static av_cold int decode_init(AVCodecContext *avctx)
{
    unsigned int c;
    unsigned int channel_size;
    ALSDecContext *ctx = avctx->priv_data;
    ctx->avctx = avctx;

    if (!avctx->extradata) {
        av_log(avctx, AV_LOG_ERROR, "Missing required ALS extradata.\n");
        return -1;
    }

    if (read_specific_config(ctx, avctx->extradata + 6,
                             avctx->extradata_size - 6, &ctx->sconf)) {
        av_log(avctx, AV_LOG_ERROR, "Reading ALSSpecificConfig failed.\n");
        decode_end(avctx);
        return -1;
    }

    avctx->sample_rate = ctx->sconf.samp_freq;
    avctx->channels    = ctx->sconf.channels;

    if (ctx->sconf.floating) {
        avctx->sample_fmt          = SAMPLE_FMT_FLT;
        avctx->bits_per_raw_sample = 32;
    } else {
        avctx->sample_fmt          = ctx->sconf.resolution > 1
                                     ? SAMPLE_FMT_S32 : SAMPLE_FMT_S16;
        avctx->bits_per_raw_sample = (ctx->sconf.resolution + 1) * 8;
    }

    avctx->frame_size = ctx->sconf.frame_length;
    channel_size      = ctx->sconf.frame_length + ctx->sconf.max_order;

    // allocate residual buffer
    if (!(ctx->residuals = av_malloc(sizeof(int64_t) * ctx->sconf.frame_length))) {
        av_log(avctx, AV_LOG_ERROR, "Allocating buffer memory failed.\n");
        decode_end(avctx);
        return AVERROR_NOMEM;
    }

    // allocate raw and carried sample buffer
    if (!(ctx->raw_buffer = av_malloc(sizeof(int64_t) *
                                      avctx->channels * channel_size))) {
        av_log(avctx, AV_LOG_ERROR, "Allocating buffer memory failed.\n");
        decode_end(avctx);
        return AVERROR_NOMEM;
    }

    // allocate raw sample array buffer
    if (!(ctx->raw_samples = av_malloc(sizeof(int64_t*) * avctx->channels))) {
        av_log(avctx, AV_LOG_ERROR, "Allocating buffer array failed.\n");
        decode_end(avctx);
        return AVERROR_NOMEM;
    }

    // allocate raw and carried samples buffers
    for (c = 0; c < avctx->channels; c++) {
        ctx->raw_samples[c] = ctx->raw_buffer + ctx->sconf.max_order +
                              c * channel_size;
    }

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

