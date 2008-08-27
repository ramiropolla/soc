/*
 * E-AC-3 decoder
 * Copyright (c) 2007 Bartlomiej Wolowiec <bartek.wolowiec@gmail.com>
 * Copyright (c) 2008 Justin Ruggles
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
#include "ac3.h"
#include "ac3_parser.h"
#include "ac3dec.h"
#include "ac3dec_data.h"

/** gain adaptive quantization mode */
typedef enum {
    EAC3_GAQ_NO =0,
    EAC3_GAQ_12,
    EAC3_GAQ_14,
    EAC3_GAQ_124
} EAC3GaqMode;

#define EAC3_SR_CODE_REDUCED  3

/** lrint(M_SQRT2*cos(2*M_PI/12)*(1<<23)) */
#define COEFF_0 10273905LL

/** lrint(M_SQRT2*cos(0*M_PI/12)*(1<<23)) = lrint(M_SQRT2*(1<<23)) */
#define COEFF_1 11863283LL

/** lrint(M_SQRT2*cos(5*M_PI/12)*(1<<23)) */
#define COEFF_2  3070444LL

/**
 * Calculate 6-point IDCT of the pre-mantissas.
 * All calculations are 24-bit fixed-point.
 */
static void idct6(int pre_mant[6])
{
    int tmp;
    int even0, even1, even2, odd0, odd1, odd2;

    odd1 = pre_mant[1] - pre_mant[3] - pre_mant[5];

    even2 = ( pre_mant[2]                * COEFF_0) >> 23;
    tmp   = ( pre_mant[4]                * COEFF_1) >> 23;
    odd0  = ((pre_mant[1] + pre_mant[5]) * COEFF_2) >> 23;

    even0 = pre_mant[0] + (tmp >> 1);
    even1 = pre_mant[0] - tmp;

    tmp = even0;
    even0 = tmp + even2;
    even2 = tmp - even2;

    tmp = odd0;
    odd0 = tmp + pre_mant[1] + pre_mant[3];
    odd2 = tmp + pre_mant[5] - pre_mant[3];

    pre_mant[0] = even0 + odd0;
    pre_mant[1] = even1 + odd1;
    pre_mant[2] = even2 + odd2;
    pre_mant[3] = even2 - odd2;
    pre_mant[4] = even1 - odd1;
    pre_mant[5] = even0 - odd0;
}

void ff_eac3_decode_transform_coeffs_aht_ch(AC3DecodeContext *s, int ch)
{
    int bin, blk, gs;
    int end_bap, gaq_mode;
    GetBitContext *gbc = &s->gbc;
    int gaq_gain[AC3_MAX_COEFS];

    gaq_mode = get_bits(gbc, 2);
    end_bap = (gaq_mode < 2) ? 12 : 17;

    /* if GAQ gain is used, decode gain codes for bins with hebap between
       8 and end_bap */
    gs = 0;
    if (gaq_mode == EAC3_GAQ_12 || gaq_mode == EAC3_GAQ_14) {
        /* read 1-bit GAQ gain codes */
        for (bin = s->start_freq[ch]; bin < s->end_freq[ch]; bin++) {
            if (s->bap[ch][bin] > 7 && s->bap[ch][bin] < end_bap)
                gaq_gain[gs++] = get_bits1(gbc) << (gaq_mode-1);
        }
    } else if (gaq_mode == EAC3_GAQ_124) {
        /* read 1.67-bit GAQ gain codes (3 codes in 5 bits) */
        int gc = 2;
        for (bin = s->start_freq[ch]; bin < s->end_freq[ch]; bin++) {
            if (s->bap[ch][bin] > 7 && s->bap[ch][bin] < end_bap) {
                if (gc++ == 2) {
                    int group_gain = get_bits(gbc, 5);
                    if (group_gain > 26) {
                        av_log(s->avctx, AV_LOG_WARNING, "GAQ gain value out-of-range\n");
                        group_gain = 26;
                    }
                    gaq_gain[gs++] = ff_ac3_ungroup_3_in_5_bits_tab[group_gain][0];
                    gaq_gain[gs++] = ff_ac3_ungroup_3_in_5_bits_tab[group_gain][1];
                    gaq_gain[gs++] = ff_ac3_ungroup_3_in_5_bits_tab[group_gain][2];
                    gc = 0;
                }
            }
        }
    }

    gs=0;
    for (bin = s->start_freq[ch]; bin < s->end_freq[ch]; bin++) {
        int hebap = s->bap[ch][bin];
        int bits = ff_eac3_bits_vs_hebap[hebap];
        if (!hebap) {
            /* zero-mantissa dithering */
            for (blk = 0; blk < 6; blk++) {
                s->pre_mantissa[ch][bin][blk] = (av_lfg_get(&s->dith_state) & 0x7FFFFF) - 0x400000;
            }
        } else if (hebap < 8) {
            /* Vector Quantization */
            int v = get_bits(gbc, bits);
            for (blk = 0; blk < 6; blk++) {
                s->pre_mantissa[ch][bin][blk] = ff_eac3_vq_hebap[hebap][v][blk] << 8;
            }
        } else {
            /* Gain Adaptive Quantization */
            int gbits, log_gain;
            if (gaq_mode != EAC3_GAQ_NO && hebap < end_bap) {
                log_gain = gaq_gain[gs++];
            } else {
                log_gain = 0;
            }
            gbits = bits - log_gain;

            for (blk = 0; blk < 6; blk++) {
                int mant = get_sbits(gbc, gbits);
                if (mant == -(1 << (gbits-1))) {
                    /* large mantissa */
                    int b;
                    mant = get_sbits(gbc, bits-2+log_gain) << (26-log_gain-bits);
                    /* remap mantissa value to correct for asymmetric quantization */
                    if (mant >= 0)
                        b = 32768 >> (log_gain+8);
                    else
                        b = ff_eac3_gaq_remap_2_4_b[hebap-8][log_gain-1];
                    mant += (ff_eac3_gaq_remap_2_4_a[hebap-8][log_gain-1] * (mant>>8) + b) >> 7;
                } else {
                    /* small mantissa, no GAQ, or Gk=1 */
                    mant <<= 24 - bits;
                    if (!log_gain) {
                        /* remap mantissa value for no GAQ or Gk=1 */
                        mant += (ff_eac3_gaq_remap_1[hebap-8] * (mant>>8)) >> 7;
                    }
                }
                s->pre_mantissa[ch][bin][blk] = mant;
            }
        }
        idct6(s->pre_mantissa[ch][bin]);
    }
}

int ff_eac3_parse_header(AC3DecodeContext *s)
{
    int i, blk, ch;
    int ac3_exponent_strategy, parse_aht_info, parse_spx_atten_data;
    int parse_transient_proc_info;
    int num_cpl_blocks;
    GetBitContext *gbc = &s->gbc;

    /* An E-AC-3 stream can have multiple independent streams which the
       application can select from. each independent stream can also contain
       dependent streams which are used to add or replace channels. */
    if (s->frame_type == EAC3_FRAME_TYPE_DEPENDENT) {
        av_log_missing_feature(s->avctx, "Dependent substream decoding", 1);
        return AC3_PARSE_ERROR_FRAME_TYPE;
    } else if (s->frame_type == EAC3_FRAME_TYPE_RESERVED) {
        av_log(s->avctx, AV_LOG_ERROR, "Reserved frame type\n");
        return AC3_PARSE_ERROR_FRAME_TYPE;
    }

    /* The substream id indicates which substream this frame belongs to. each
       independent stream has its own substream id, and the dependent streams
       associated to an independent stream have matching substream id's. */
    if (s->substreamid) {
        /* only decode substream with id=0. skip any additional substreams. */
        av_log_missing_feature(s->avctx, "Additional substreams", 1);
        return AC3_PARSE_ERROR_FRAME_TYPE;
    }

    if (s->bit_alloc_params.sr_code == EAC3_SR_CODE_REDUCED) {
        /* The E-AC-3 specification does not tell how to handle reduced sample
           rates in bit allocation.  The best assumption would be that it is
           handled like AC-3 DolbyNet, but we cannot be sure until we have a
           sample which utilizes this feature. */
        av_log_missing_feature(s->avctx, "Reduced sampling rates", 1);
        return -1;
    }
    skip_bits(gbc, 5); // skip bitstream id

    /* volume control params */
    for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
        skip_bits(gbc, 5); // skip dialog normalization
        if (get_bits1(gbc)) {
            skip_bits(gbc, 8); // skip compression gain word
        }
    }

    /* dependent stream channel map */
    if (s->frame_type == EAC3_FRAME_TYPE_DEPENDENT) {
        if (get_bits1(gbc)) {
            skip_bits(gbc, 16); // skip custom channel map
        }
    }

    /* mixing metadata */
    if (get_bits1(gbc)) {
        /* center and surround mix levels */
        if (s->channel_mode > AC3_CHMODE_STEREO) {
            skip_bits(gbc, 2);  // skip preferred stereo downmix mode
            if (s->channel_mode & 1) {
                /* if three front channels exist */
                skip_bits(gbc, 3); //skip Lt/Rt center mix level
                s->center_mix_level = get_bits(gbc, 3);
            }
            if (s->channel_mode & 4) {
                /* if a surround channel exists */
                skip_bits(gbc, 3); //skip Lt/Rt surround mix level
                s->surround_mix_level = get_bits(gbc, 3);
            }
        }

        /* lfe mix level */
        if (s->lfe_on && get_bits1(gbc)) {
            // TODO: use LFE mix level
            skip_bits(gbc, 5); // skip LFE mix level code
        }

        /* info for mixing with other streams and substreams */
        if (s->frame_type == EAC3_FRAME_TYPE_INDEPENDENT) {
            for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
                // TODO: apply program scale factor
                if (get_bits1(gbc)) {
                    skip_bits(gbc, 6);  // skip program scale factor
                }
            }
            if (get_bits1(gbc)) {
                skip_bits(gbc, 6);  // skip external program scale factor
            }
            /* skip mixing parameter data */
            switch(get_bits(gbc, 2)) {
                case 1: skip_bits(gbc, 5);  break;
                case 2: skip_bits(gbc, 12); break;
                case 3: {
                    int mix_data_size = (get_bits(gbc, 5) + 2) << 3;
                    skip_bits_long(gbc, mix_data_size);
                    break;
                }
            }
            /* skip pan information for mono or dual mono source */
            if (s->channel_mode < AC3_CHMODE_STEREO) {
                for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
                    if (get_bits1(gbc)) {
                        /* note: this is not in the ATSC A/52B specification
                           reference: ETSI TS 102 366 V1.1.1
                                      section: E.1.3.1.25 */
                        skip_bits(gbc, 8);  // skip pan mean direction index
                        skip_bits(gbc, 6);  // skip reserved paninfo bits
                    }
                }
            }
            /* skip mixing configuration information */
            if (get_bits1(gbc)) {
                if (s->num_blocks == 1) {
                    skip_bits(gbc, 5);
                } else {
                    for (blk = 0; blk < s->num_blocks; blk++) {
                        if (get_bits1(gbc)) {
                            skip_bits(gbc, 5);
                        }
                    }
                }
            }
        }
    }

    /* informational metadata */
    if (get_bits1(gbc)) {
        skip_bits(gbc, 3); // skip bit stream mode
        skip_bits(gbc, 2); // skip copyright bit and original bitstream bit
        if (s->channel_mode == AC3_CHMODE_STEREO) {
            skip_bits(gbc, 4); // skip Dolby surround and headphone mode
        }
        if (s->channel_mode >= AC3_CHMODE_2F2R) {
            skip_bits(gbc, 2); // skip Dolby surround EX mode
        }
        for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
            if (get_bits1(gbc)) {
                skip_bits(gbc, 8); // skip mix level, room type, and A/D converter type
            }
        }
        if (s->bit_alloc_params.sr_code != EAC3_SR_CODE_REDUCED) {
            skip_bits1(gbc); // skip source sample rate code
        }
    }

    /* converter synchronization flag
       If frames are less than six blocks, this bit should be turned on
       once every 6 blocks to indicate the start of a frame set.
       reference: RFC 4598, Section 2.1.3  Frame Sets */
    if (s->frame_type == EAC3_FRAME_TYPE_INDEPENDENT && s->num_blocks != 6) {
        skip_bits1(gbc); // skip converter synchronization flag
    }

    /* original frame size code if this stream was converted from AC-3 */
    if (s->frame_type == EAC3_FRAME_TYPE_AC3_CONVERT &&
            (s->num_blocks == 6 || get_bits1(gbc))) {
        skip_bits(gbc, 6); // skip frame size code
    }

    /* additional bitstream info */
    if (get_bits1(gbc)) {
        int addbsil = get_bits(gbc, 6);
        for (i = 0; i < addbsil + 1; i++) {
            skip_bits(gbc, 8); // skip additional bit stream info
        }
    }

    /* audio frame syntax flags, strategy data, and per-frame data */

    if (s->num_blocks == 6) {
        ac3_exponent_strategy = get_bits1(gbc);
        parse_aht_info = get_bits1(gbc);
    } else {
        /* less than 6 blocks, so use AC-3-style exponent strategy syntax, and
           do not use AHT */
        ac3_exponent_strategy = 1;
        parse_aht_info = 0;
    }

    s->snr_offset_strategy    = get_bits(gbc, 2);
    parse_transient_proc_info = get_bits1(gbc);

    s->block_switch_syntax = get_bits1(gbc);
    if (!s->block_switch_syntax)
        memset(s->block_switch, 0, sizeof(s->block_switch));

    s->dither_flag_syntax = get_bits1(gbc);
    if (!s->dither_flag_syntax) {
        s->dither_all = 1;
        for (ch = 1; ch <= s->fbw_channels; ch++)
            s->dither_flag[ch] = 1;
    }
    s->dither_flag[CPL_CH] = s->dither_flag[s->lfe_ch] = 0;

    s->bit_allocation_syntax = get_bits1(gbc);
    if (!s->bit_allocation_syntax) {
        /* set default bit allocation parameters */
        s->bit_alloc_params.slow_decay = ff_ac3_slow_decay_tab[2];
        s->bit_alloc_params.fast_decay = ff_ac3_fast_decay_tab[1];
        s->bit_alloc_params.slow_gain  = ff_ac3_slow_gain_tab [1];
        s->bit_alloc_params.db_per_bit = ff_ac3_db_per_bit_tab[2];
        s->bit_alloc_params.floor      = ff_ac3_floor_tab     [7];
    }

    s->fast_gain_syntax  = get_bits1(gbc);
    s->dba_syntax        = get_bits1(gbc);
    s->skip_syntax       = get_bits1(gbc);
    parse_spx_atten_data = get_bits1(gbc);

    /* coupling strategy occurance and coupling use per block */
    num_cpl_blocks = 0;
    if (s->channel_mode > 1) {
        for (blk = 0; blk < s->num_blocks; blk++) {
            s->cpl_strategy_exists[blk] = (!blk || get_bits1(gbc));
            if (s->cpl_strategy_exists[blk]) {
                s->cpl_in_use[blk] = get_bits1(gbc);
            } else {
                s->cpl_in_use[blk] = s->cpl_in_use[blk-1];
            }
            num_cpl_blocks += s->cpl_in_use[blk];
        }
    } else {
        memset(s->cpl_in_use, 0, sizeof(s->cpl_in_use));
    }

    /* exponent strategy data */
    if (ac3_exponent_strategy) {
        /* AC-3-style exponent strategy syntax */
        for (blk = 0; blk < s->num_blocks; blk++) {
            for (ch = !s->cpl_in_use[blk]; ch <= s->fbw_channels; ch++) {
                s->exp_strategy[blk][ch] = get_bits(gbc, 2);
            }
        }
    } else {
        /* LUT-based exponent strategy syntax */
        int frmchexpstr;
        for (ch = !((s->channel_mode > 1) && num_cpl_blocks); ch <= s->fbw_channels; ch++) {
            frmchexpstr = get_bits(gbc, 5);
            for (blk = 0; blk < 6; blk++) {
                s->exp_strategy[blk][ch] = ff_eac3_frm_expstr[frmchexpstr][blk];
            }
        }
    }
    /* LFE exponent strategy */
    if (s->lfe_on) {
        for (blk = 0; blk < s->num_blocks; blk++) {
            s->exp_strategy[blk][s->lfe_ch] = get_bits1(gbc);
        }
    }
    /* original exponent strategies if this stream was converted from AC-3 */
    if (s->frame_type == EAC3_FRAME_TYPE_INDEPENDENT &&
            (s->num_blocks == 6 || get_bits1(gbc))) {
        for (ch = 1; ch <= s->fbw_channels; ch++) {
            skip_bits(gbc, 5); // skip converter channel exponent strategy
        }
    }

    /* determine which channels use AHT */
    if (parse_aht_info) {
        /* AHT is only available when there are 6 blocks in the frame.
           The coupling channel can only use AHT when coupling is in use for
           all blocks.
           reference: Section E3.3.2 Bit Stream Helper Variables */
        s->channel_uses_aht[CPL_CH]=0;
        for (ch = (num_cpl_blocks != 6); ch <= s->channels; ch++) {
            int nchregs = 0;
            for (blk = 0; blk < 6; blk++) {
                if (ch)
                    nchregs += (s->exp_strategy[blk][ch] != EXP_REUSE);
                else
                    nchregs += s->cpl_strategy_exists[blk] ||
                               (s->exp_strategy[blk][CPL_CH] != EXP_REUSE);
            }
            s->channel_uses_aht[ch] = (nchregs == 1) && get_bits1(gbc);
        }
    } else {
        memset(s->channel_uses_aht, 0, sizeof(s->channel_uses_aht));
    }

    /* per-frame SNR offset */
    if (!s->snr_offset_strategy) {
        int csnroffst = (get_bits(gbc, 6) - 15) << 4;
        int snroffst = (csnroffst + get_bits(gbc, 4)) << 2;
        for (ch = 0; ch <= s->channels; ch++)
            s->snr_offset[ch] = snroffst;
    }

    /* transient pre-noise processing data */
    if (parse_transient_proc_info) {
        for (ch = 1; ch <= s->fbw_channels; ch++) {
            if (get_bits1(gbc)) { // channel in transient processing
                skip_bits(gbc, 10); // skip transient processing location
                skip_bits(gbc, 8);  // skip transient processing length
            }
        }
    }

    /* block start information */
    if (s->num_blocks > 1 && get_bits1(gbc)) {
        /* reference: Section E2.3.2.27
           nblkstrtbits = (numblks - 1) * (4 + ceiling(log2(words_per_frame)))
           The spec does not say what this data is or what it's used for.
           It is likely the offset of each block within the frame. */
        int block_start_bits = (s->num_blocks-1) * (4 + av_log2(s->frame_size-2));
        skip_bits(gbc, block_start_bits);
    }

    /* syntax state initialization */
    for (ch = 1; ch <= s->fbw_channels; ch++) {
        s->first_cpl_coords[ch] = 1;
    }
    s->first_cpl_leak = 1;

    return 0;
}
