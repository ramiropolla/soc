/*
 * AAC Spectral Band Replication decoding functions
 * Copyright (c) 2008-2009 Robert Swain ( rob opendot cl )
 * Copyright (c) 2009 Alex Converse <alex.converse@gmail.com>
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
 * AAC Spectral Band Replication decoding functions
 * @author Robert Swain ( rob opendot cl )
 */

#include "aac.h"
#include "aacsbr.h"
#include "aacsbrdata.h"

#include <stdint.h>

static VLC vlc_sbr[10];
static const int8_t vlc_sbr_lav[10] =
    { 60, 60, 24, 24, 31, 31, 12, 12, 31, 12 };

av_cold void ff_aac_sbr_init(void)
{
    static const struct {
        const void *sbr_codes, *sbr_bits;
        const unsigned int table_size, elem_size;
    } sbr_tmp[] = {
        { t_huffman_env_1_5dB_codes,       t_huffman_env_1_5dB_bits,
            sizeof(t_huffman_env_1_5dB_codes),       sizeof(t_huffman_env_1_5dB_codes[0]) },
        { f_huffman_env_1_5dB_codes,       f_huffman_env_1_5dB_bits,
            sizeof(f_huffman_env_1_5dB_codes),       sizeof(f_huffman_env_1_5dB_codes[0]) },
        { t_huffman_env_bal_1_5dB_codes,   t_huffman_env_bal_1_5dB_bits,
            sizeof(t_huffman_env_bal_1_5dB_codes),   sizeof(t_huffman_env_bal_1_5dB_codes[0]) },
        { f_huffman_env_bal_1_5dB_codes,   f_huffman_env_bal_1_5dB_bits,
            sizeof(f_huffman_env_bal_1_5dB_codes),   sizeof(f_huffman_env_bal_1_5dB_codes[0]) },
        { t_huffman_env_3_0dB_codes,       t_huffman_env_3_0dB_bits,
            sizeof(t_huffman_env_3_0dB_codes),       sizeof(t_huffman_env_3_0dB_codes[0]) },
        { f_huffman_env_3_0dB_codes,       f_huffman_env_3_0dB_bits,
            sizeof(f_huffman_env_3_0dB_codes),       sizeof(f_huffman_env_3_0dB_codes[0]) },
        { t_huffman_env_bal_3_0dB_codes,   t_huffman_env_bal_3_0dB_bits,
            sizeof(t_huffman_env_bal_3_0dB_codes),   sizeof(t_huffman_env_bal_3_0dB_codes[0]) },
        { f_huffman_env_bal_3_0dB_codes,   f_huffman_env_bal_3_0dB_bits,
            sizeof(f_huffman_env_bal_3_0dB_codes),   sizeof(f_huffman_env_bal_3_0dB_codes[0]) },
        { t_huffman_noise_3_0dB_codes,     t_huffman_noise_3_0dB_bits,
            sizeof(t_huffman_noise_3_0dB_codes),     sizeof(t_huffman_noise_3_0dB_codes[0]) },
        { t_huffman_noise_bal_3_0dB_codes, t_huffman_noise_bal_3_0dB_bits,
            sizeof(t_huffman_noise_bal_3_0dB_codes), sizeof(t_huffman_noise_bal_3_0dB_codes[0]) },
    };

    // SBR VLC table initialization
    SBR_INIT_VLC_STATIC(0, 1098);
    SBR_INIT_VLC_STATIC(1, 1092);
    SBR_INIT_VLC_STATIC(2, 768);
    SBR_INIT_VLC_STATIC(3, 1026);
    SBR_INIT_VLC_STATIC(4, 1058);
    SBR_INIT_VLC_STATIC(5, 1052);
    SBR_INIT_VLC_STATIC(6, 544);
    SBR_INIT_VLC_STATIC(7, 544);
    SBR_INIT_VLC_STATIC(8, 592);
    SBR_INIT_VLC_STATIC(9, 512);
}

av_cold void ff_aac_sbr_ctx_init(SpectralBandReplication *sbr)
{
    sbr->k[4] = sbr->k[3] = 32; //Typo in spec, kx' inits to 32
}

static unsigned int sbr_header(SpectralBandReplication *sbr, GetBitContext *gb)
{
    unsigned int cnt = get_bits_count(gb);
    uint8_t bs_header_extra_1;
    uint8_t bs_header_extra_2;

    sbr->start = 1;

    // Save last spectrum parameters variables to compare to new ones
    memcpy(&sbr->spectrum_params[0], &sbr->spectrum_params[1], sizeof(SpectrumParameters));

    sbr->bs_amp_res_header                 = get_bits1(gb);
    sbr->spectrum_params[1].bs_start_freq  = get_bits(gb, 4);
    sbr->spectrum_params[1].bs_stop_freq   = get_bits(gb, 4);
    sbr->spectrum_params[1].bs_xover_band  = get_bits(gb, 3);
                                             skip_bits(gb, 2); // bs_reserved

    bs_header_extra_1 = get_bits1(gb);
    bs_header_extra_2 = get_bits1(gb);

    if (bs_header_extra_1) {
        sbr->spectrum_params[1].bs_freq_scale  = get_bits(gb, 2);
        sbr->spectrum_params[1].bs_alter_scale = get_bits1(gb);
        sbr->spectrum_params[1].bs_noise_bands = get_bits(gb, 2);
    } else {
        sbr->spectrum_params[1].bs_freq_scale  = 2;
        sbr->spectrum_params[1].bs_alter_scale = 1;
        sbr->spectrum_params[1].bs_noise_bands = 2;
    }

    // Check if spectrum parameters changed
    if (memcmp(&sbr->spectrum_params[0], &sbr->spectrum_params[1],
               sizeof(SpectrumParameters)))
        sbr->reset = 1;

    if (bs_header_extra_2) {
        sbr->bs_limiter_bands  = get_bits(gb, 2);
        sbr->bs_limiter_gains  = get_bits(gb, 2);
        sbr->bs_interpol_freq  = get_bits1(gb);
        sbr->bs_smoothing_mode = get_bits1(gb);
    } else {
        sbr->bs_limiter_bands  = 2;
        sbr->bs_limiter_gains  = 2;
        sbr->bs_interpol_freq  = 1;
        sbr->bs_smoothing_mode = 1;
    }

    return get_bits_count(gb) - cnt;
}

static int array_min_int16(int16_t *array, int nel)
{
    int i, min = array[0];
    for (i = 1; i < nel; i++)
        if (array[i] < min)
            min = array[i];
    return min;
}

static int qsort_comparison_function_int(const void *a, const void *b)
{
    return *(const int *)a - *(const int *)b;
}

static int qsort_comparison_function_int16(const void *a, const void *b)
{
    return *(const int16_t *)a - *(const int16_t *)b;
}

// Master Frequency Band Table (14496-3 sp04 p194)
static int sbr_make_f_master(AACContext *ac, SpectralBandReplication *sbr,
                             SpectrumParameters *spectrum)
{
    unsigned int temp;
    unsigned int start_min, stop_min;
    int k;
    const uint8_t *sbr_offset_ptr;
    int stop_dk[13];

    if (sbr->sample_rate < 32000) {
        temp = 3000;
    } else if (sbr->sample_rate < 64000) {
        temp = 4000;
    } else
        temp = 5000;

    start_min = (unsigned int)lroundf((temp << 7) / (float)sbr->sample_rate);
    stop_min  = (unsigned int)lroundf((temp << 8) / (float)sbr->sample_rate);

    if (sbr->sample_rate == 16000) {
        sbr_offset_ptr = sbr_offset[0];
    } else if (sbr->sample_rate == 22050) {
        sbr_offset_ptr = sbr_offset[1];
    } else if (sbr->sample_rate == 24000) {
        sbr_offset_ptr = sbr_offset[2];
    } else if (sbr->sample_rate == 32000) {
        sbr_offset_ptr = sbr_offset[3];
    } else if ((sbr->sample_rate >= 44100) &&
               (sbr->sample_rate <= 64000)) {
        sbr_offset_ptr = sbr_offset[4];
    } else if (sbr->sample_rate > 64000) {
        sbr_offset_ptr = sbr_offset[5];
    } else {
        av_log(ac->avccontext, AV_LOG_ERROR, "Unsupported sample rate for SBR: %d\n", sbr->sample_rate);
        return -1;
    }

    sbr->k[0] = start_min + sbr_offset_ptr[spectrum->bs_start_freq];

    if (spectrum->bs_stop_freq < 14) {
        sbr->k[2] = stop_min;
        for (k = 0; k < 13; k++) {
            stop_dk[k] = lroundf(stop_min * powf(64.0f / (float)stop_min, (k + 1) / 13.0f))  -
                         lroundf(stop_min * powf(64.0f / (float)stop_min,  k      / 13.0f));
        }
        qsort(stop_dk, 13, sizeof(stop_dk[0]), qsort_comparison_function_int);
        for (k = 0; k < spectrum->bs_stop_freq; k++)
            sbr->k[2] += stop_dk[k];
    } else if (spectrum->bs_stop_freq == 14) {
        sbr->k[2] = 2*sbr->k[0];
    } else if (spectrum->bs_stop_freq == 15) {
        sbr->k[2] = 3*sbr->k[0];
    } else {
        av_log(ac->avccontext, AV_LOG_ERROR, "Invalid bs_stop_freq: %d\n", spectrum->bs_stop_freq);
        return -1;
    }
    sbr->k[2] = FFMIN(64, sbr->k[2]);

    if (!spectrum->bs_freq_scale) {
        unsigned int dk;
        int k2diff;

        if (!spectrum->bs_alter_scale) {
            dk = 1;
            sbr->n_master = ((unsigned int)((sbr->k[2] - sbr->k[0]) / (float)(dk << 1))) << 1;
        } else {
            dk = 2;
            sbr->n_master =         lroundf((sbr->k[2] - sbr->k[0]) / (float)(dk << 1))  << 1;
        }

        for (k = 1; k <= sbr->n_master; k++)
            sbr->f_master[k] = dk;

        k2diff = sbr->k[2] - sbr->k[0] - sbr->n_master * dk;
        if (k2diff) {
            int incr;
            if (k2diff < 0) {
                incr = 1;
                k    = 1;
            } else {
                incr = -1;
                k    = sbr->n_master;
            }

            while (k2diff) {
                sbr->f_master[k] -= incr;
                k                += incr;
                k2diff           += incr;
            }
        }

        sbr->f_master[0] = sbr->k[0];
        for (k = 1; k <= sbr->n_master; k++)
            sbr->f_master[k] += sbr->f_master[k - 1];

        // Requirements (14496-3 sp04 p205)
        if (sbr->n_master <= 0) {
            av_log(ac->avccontext, AV_LOG_ERROR, "Invalid n_master: %d\n", sbr->n_master);
            return -1;
        }
    } else {
        unsigned int bands = 14 - (spectrum->bs_freq_scale << 1); // bs_freq_scale  = {1,2,3}
        float warp = spectrum->bs_alter_scale ? 1.3 : 1.0;        // bs_alter_scale = {0,1}
        unsigned int two_regions, num_bands_0;
        int vdk0_max, vdk1_min;
        int16_t vk0[49];

        if (sbr->k[2] / (float)sbr->k[0] > 2.2449) {
            two_regions = 1;
            sbr->k[1] = sbr->k[0] << 1;
        } else {
            two_regions = 0;
            sbr->k[1] = sbr->k[2];
        }

        num_bands_0 = lroundf(bands * logf(sbr->k[1] / (float)sbr->k[0]) / (2.0f * logf(2.0f))) << 1;

        if (num_bands_0 <= 0) { // Requirements (14496-3 sp04 p205)
            av_log(ac->avccontext, AV_LOG_ERROR, "Invalid num_bands_0: %d\n", num_bands_0);
            return -1;
        }

        vk0[0] = 0;

        for (k = 0; k < num_bands_0; k++) {
            vk0[k + 1] = lroundf(sbr->k[0] * powf(sbr->k[1] / (float)sbr->k[0], (k + 1) / (float)num_bands_0)) -
                         lroundf(sbr->k[0] * powf(sbr->k[1] / (float)sbr->k[0],  k      / (float)num_bands_0));
        }

        qsort(vk0 + 1, num_bands_0, sizeof(vk0[1]), qsort_comparison_function_int16);
        vdk0_max = vk0[num_bands_0];

        vk0[0] = sbr->k[0];
        for (k = 1; k <= num_bands_0; k++) {
            if (vk0[k] <= 0) { // Requirements (14496-3 sp04 p205)
                av_log(ac->avccontext, AV_LOG_ERROR, "Invalid vDk0[%d]: %d\n", k, vk0[k]);
                return -1;
            }
            vk0[k] += vk0[k-1];
        }

        if (two_regions) {
            int16_t vk1[49];
            unsigned int num_bands_1 = lroundf(bands * logf(sbr->k[2] / (float)sbr->k[1]) /
                                              (2.0f * logf(2.0f) * warp)) << 1;

            for (k = 0; k < num_bands_1; k++) {
                vk1[k + 1] = lroundf(sbr->k[1] * powf(sbr->k[2] / (float)sbr->k[1], (k + 1) / (float)num_bands_1)) -
                             lroundf(sbr->k[1] * powf(sbr->k[2] / (float)sbr->k[1],  k      / (float)num_bands_1));
            }

            vdk1_min = array_min_int16(vk1 + 1, num_bands_1);

            if (vdk1_min < vdk0_max) {
                int change;
                qsort(vk1 + 1, num_bands_1, sizeof(vk1[1]), qsort_comparison_function_int16);
                change = FFMIN(vdk0_max - vk1[1], (vk1[num_bands_1] - vk1[1]) >> 1);
                vk1[1]           += change;
                vk1[num_bands_1] -= change;
            }

            qsort(vk1 + 1, num_bands_1, sizeof(vk1[1]), qsort_comparison_function_int16);

            vk1[0] = sbr->k[1];
            for (k = 1; k <= num_bands_1; k++) {
                if (vk1[k] <= 0) { // Requirements (14496-3 sp04 p205)
                    av_log(ac->avccontext, AV_LOG_ERROR, "Invalid vDk1[%d]: %d\n", k, vk1[k]);
                    return -1;
                }
                vk1[k] += vk1[k-1];
            }

            sbr->n_master = num_bands_0 + num_bands_1;
            memcpy(&sbr->f_master[0],               vk0,    (num_bands_0 + 1) * sizeof(sbr->f_master[0]));
            memcpy(&sbr->f_master[num_bands_0 + 1], vk1 + 1, num_bands_1      * sizeof(sbr->f_master[0]));

        } else {
            sbr->n_master = num_bands_0;
            memcpy(sbr->f_master, vk0, (num_bands_0 + 1) * sizeof(sbr->f_master[0]));
        }
    }
    // Requirements (14496-3 sp04 p205)
    if (sbr->spectrum_params[1].bs_xover_band >= sbr->n_master) {
        av_log(ac->avccontext, AV_LOG_ERROR,
               "Invalid bitstream, crossover band index beyond array bounds: %d\n",
               sbr->spectrum_params[1].bs_xover_band);
        return -1;
    }
    // temp == max number of QMF subbands
    if (sbr->sample_rate <= 32000) {
        temp = 48;
    } else if (sbr->sample_rate == 44100) {
        temp = 35;
    } else if (sbr->sample_rate >= 48000)
        temp = 32;

    if (sbr->k[2] - sbr->k[0] > temp) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Invalid bitstream, too many QMF subbands: %d\n", sbr->k[2] - sbr->k[0]);
        return -1;
    }

    return 0;
}

// High Frequency Generation - Patch Construction (14496-3 sp04 p216 fig. 4.46)
static int sbr_hf_calc_npatches(AACContext *ac, SpectralBandReplication *sbr)
{
    int i, k, sb = 0;
    int msb = sbr->k[0];
    int usb = sbr->k[3];
    int goal_sb = lroundf((1 << 11) * 1000 / (float)sbr->sample_rate);

    sbr->num_patches = 0;

    if (goal_sb < sbr->k[3] + sbr->m) {
        for (k = 0; sbr->f_master[k] < goal_sb; k++);
    } else
        k = sbr->n_master;

    do {
        int odd = 0;
        for (i = k; i == k || sb > (sbr->k[0] - 1 + msb - odd); i--) {
            sb = sbr->f_master[i];
            odd = (sb - 2 + sbr->k[0]) & 1;
        }

        sbr->patch_num_subbands[sbr->num_patches]  = FFMAX(sb - usb, 0);
        sbr->patch_start_subband[sbr->num_patches] = sbr->k[0] - odd - sbr->patch_num_subbands[sbr->num_patches];

        if (sbr->patch_num_subbands[sbr->num_patches] > 0) {
            usb = sb;
            msb = sb;
            sbr->num_patches++;
        } else
            msb = sbr->k[3];

        if (sbr->f_master[k] - sb < 3)
            k = sbr->n_master;
    } while (sb != sbr->k[3] + sbr->m);

    if ((sbr->patch_num_subbands[sbr->num_patches-1] < 3) && (sbr->num_patches > 1))
        sbr->num_patches--;

    if (sbr->num_patches > 5) { // Requirements (14496-3 sp04 p205)
        av_log(ac->avccontext, AV_LOG_ERROR, "Too many patches: %d\n", sbr->num_patches);
        return -1;
    }

    return 0;
}

static inline void remove_table_element(void *table, uint8_t *last_el, int el_size,
                                        int el)
{
    memmove((uint8_t *)table + el_size*el, (uint8_t *)table + el_size*(el + 1), (*last_el - el)*el_size);
    (*last_el)--;
}

static inline int in_table(void *table, int last_el, int el_size, void *needle)
{
    int i;
    uint8_t *table_ptr = table; // avoids a warning with void * ptr arith
    for (i = 0; i <= last_el; i++, table_ptr += el_size)
        if (!memcmp(table_ptr, needle, el_size))
            return 1;
    return 0;
}

// Derived Frequency Band Tables (14496-3 sp04 p197)
static int sbr_make_f_derived(AACContext *ac, SpectralBandReplication *sbr)
{
    int k, temp;

    sbr->n[1] = sbr->n_master - sbr->spectrum_params[1].bs_xover_band;
    sbr->n[0] = (sbr->n[1] + 1) >> 1;

    memcpy(sbr->f_tablehigh, &sbr->f_master[sbr->spectrum_params[1].bs_xover_band],
           (sbr->n[1] + 1) * sizeof(sbr->f_master[0]));
    sbr->mold = sbr->m;
    sbr->m    = sbr->f_tablehigh[sbr->n[1]] - sbr->f_tablehigh[0];
    sbr->k[4] = sbr->k[3];
    sbr->k[3] = sbr->f_tablehigh[0];

    // Requirements (14496-3 sp04 p205)
    if (sbr->k[3] + sbr->m > 64) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Stop frequency border too high: %d\n", sbr->k[3] + sbr->m);
        return -1;
    }
    if (sbr->k[3] > 32) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Start frequency border too high: %d\n", sbr->k[3]);
        return -1;
    }

    sbr->f_tablelow[0] = sbr->f_tablehigh[0];
    temp = sbr->n[1] & 1;
    for (k = 1; k <= sbr->n[0]; k++)
        sbr->f_tablelow[k] = sbr->f_tablehigh[(k << 1) - temp];

    sbr->n_q = FFMAX(1, lroundf(sbr->spectrum_params[1].bs_noise_bands * logf(sbr->k[2] / (float)sbr->k[3]) / logf(2.0f))); // 0 <= bs_noise_bands <= 3
    if (sbr->n_q > 5) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Too many noise floor scale factors: %d\n", sbr->n_q);
        return -1;
    }

    sbr->f_tablenoise[0] = sbr->f_tablelow[0];
    temp = 0;
    for (k = 1; k <= sbr->n_q; k++) {
        temp += (sbr->n[0] - temp) / (sbr->n_q + 1 - k);
        sbr->f_tablenoise[k] = sbr->f_tablelow[temp];
    }

    sbr_hf_calc_npatches(ac, sbr);

    // Limiter Frequency Band Table (14496-3 sp04 p198)
    if (sbr->bs_limiter_bands > 0) {
        const float lim_bands_per_octave[3] = {1.2, 2, 3};
        int16_t patch_borders[sbr->num_patches + 1];

        patch_borders[0] = sbr->k[3];
        for (k=1; k <= sbr->num_patches; k++)
            patch_borders[k] = patch_borders[k-1] + sbr->patch_num_subbands[k-1];

        memcpy( sbr->f_tablelim,                  sbr->f_tablelow, (sbr->n[0]        + 1) * sizeof(sbr->f_tablelow[0]));
        memcpy(&sbr->f_tablelim[sbr->n[0] + 1], &patch_borders[1], (sbr->num_patches - 1) * sizeof(patch_borders[0]));

        qsort(sbr->f_tablelim, sbr->num_patches + sbr->n[0], sizeof(sbr->f_tablelim[0]), qsort_comparison_function_int16);

        k = 1;
        sbr->n_lim = sbr->n[0] + sbr->num_patches - 1;
        while (k <= sbr->n_lim) {
            // if ( nOctaves * limBands >= 0.49) ...
            if (log2(sbr->f_tablelim[k] / (float)sbr->f_tablelim[k-1]) * lim_bands_per_octave[sbr->bs_limiter_bands - 1] >= 0.49) {
                k++;
                continue;
            }
            if (sbr->f_tablelim[k] == sbr->f_tablelim[k-1] ||
                !in_table(patch_borders, sbr->num_patches, sizeof(patch_borders[0]), &sbr->f_tablelim[k]))
                remove_table_element(sbr->f_tablelim, &sbr->n_lim, sizeof(sbr->f_tablelim[0]), k);
            else if (!in_table(patch_borders, sbr->num_patches, sizeof(patch_borders[0]), &sbr->f_tablelim[k-1]))
                remove_table_element(sbr->f_tablelim, &sbr->n_lim, sizeof(sbr->f_tablelim[0]), k-1);
            else
                k++;
        };
    } else {
        sbr->f_tablelim[0] = sbr->f_tablelow[0];
        sbr->f_tablelim[1] = sbr->f_tablelow[sbr->n[0]];
        sbr->n_lim = 1;
    }

    sbr->f_indexnoise[0] = 0;
    sbr->f_indexnoise[1] = 0;

    return 0;
}

static int8_t ceil_log2[] = {
    0, 0, 1, 2, 2, 3, 3,
};

static int sbr_grid(AACContext *ac, SpectralBandReplication *sbr,
                    GetBitContext *gb, SBRData *ch_data)
{
    int i;

    ch_data->bs_freq_res[0] = ch_data->bs_freq_res[ch_data->bs_num_env[1]];
    ch_data->bs_num_env[0] = ch_data->bs_num_env[1];
    ch_data->bs_amp_res = sbr->bs_amp_res_header;

    switch (ch_data->bs_frame_class = get_bits(gb, 2)) {
    case FIXFIX:
        ch_data->bs_num_env[1] = 1 << get_bits(gb, 2);
        if (ch_data->bs_num_env[1] == 1)
            ch_data->bs_amp_res = 0;

        ch_data->bs_freq_res[1] = get_bits1(gb);
        for (i = 1; i < ch_data->bs_num_env[1]; i++)
            ch_data->bs_freq_res[i + 1] = ch_data->bs_freq_res[1];
        break;
    case FIXVAR:
        ch_data->bs_var_bord[1] = get_bits(gb, 2);
        ch_data->bs_num_rel[1]  = get_bits(gb, 2);
        ch_data->bs_num_env[1] = ch_data->bs_num_rel[1] + 1;

        for (i = 0; i < ch_data->bs_num_rel[1]; i++)
            ch_data->bs_rel_bord[1][i] = (get_bits(gb, 2) << 1) + 2;

        ch_data->bs_pointer = get_bits(gb, ceil_log2[ch_data->bs_num_env[1] + 1]);

        for (i = 0; i < ch_data->bs_num_env[1]; i++)
            ch_data->bs_freq_res[ch_data->bs_num_env[1] - i] = get_bits1(gb);
        break;
    case VARFIX:
        ch_data->bs_var_bord[0] = get_bits(gb, 2);
        ch_data->bs_num_rel[0]  = get_bits(gb, 2);
        ch_data->bs_num_env[1]  = ch_data->bs_num_rel[0] + 1;

        for (i = 0; i < ch_data->bs_num_rel[0]; i++)
            ch_data->bs_rel_bord[0][i] = (get_bits(gb, 2) << 1) + 2;

        ch_data->bs_pointer = get_bits(gb, ceil_log2[ch_data->bs_num_env[1] + 1]);

        for (i = 0; i < ch_data->bs_num_env[1]; i++)
            ch_data->bs_freq_res[i + 1] = get_bits1(gb);
        break;
    case VARVAR:
        ch_data->bs_var_bord[0] = get_bits(gb, 2);
        ch_data->bs_var_bord[1] = get_bits(gb, 2);
        ch_data->bs_num_rel[0]  = get_bits(gb, 2);
        ch_data->bs_num_rel[1]  = get_bits(gb, 2);
        ch_data->bs_num_env[1]  = ch_data->bs_num_rel[0] + ch_data->bs_num_rel[1] + 1;

        for (i = 0; i < ch_data->bs_num_rel[0]; i++)
            ch_data->bs_rel_bord[0][i] = (get_bits(gb, 2) << 1) + 2;
        for (i = 0; i < ch_data->bs_num_rel[1]; i++)
            ch_data->bs_rel_bord[1][i] = (get_bits(gb, 2) << 1) + 2;

        ch_data->bs_pointer = get_bits(gb, ceil_log2[ch_data->bs_num_env[1] + 1]);

        for (i = 0; i < ch_data->bs_num_env[1]; i++)
            ch_data->bs_freq_res[i + 1] = get_bits1(gb);
        break;
    default:
        break;
    }

    if (ch_data->bs_frame_class == FIXFIX && ch_data->bs_num_env[1] > 4) {
        av_log(ac->avccontext, AV_LOG_ERROR,
               "Invalid bitstream, too many SBR envelopes in FIXFIX type SBR frame: %d\n", ch_data->bs_num_env[1]);
        return -1;
    }
    if (ch_data->bs_frame_class == VARVAR && ch_data->bs_num_env[1] > 5) {
        av_log(ac->avccontext, AV_LOG_ERROR,
               "Invalid bitstream, too many SBR envelopes in VARVAR type SBR frame: %d\n", ch_data->bs_num_env[1]);
        return -1;
    }

    ch_data->bs_num_noise = ch_data->bs_num_env[1] > 1 ? 2 : 1;

    return 0;
}

static void sbr_dtdf(SpectralBandReplication *sbr, GetBitContext *gb,
                     SBRData *ch_data)
{
    int i;

    for (i = 0; i < ch_data->bs_num_env[1]; i++)
        ch_data->bs_df_env[i] = get_bits1(gb);
    for (i = 0; i < ch_data->bs_num_noise; i++)
        ch_data->bs_df_noise[i] = get_bits1(gb);
}

static void sbr_invf(SpectralBandReplication *sbr, GetBitContext *gb,
                     SBRData *ch_data)
{
    int i;

    memcpy(ch_data->bs_invf_mode[1], ch_data->bs_invf_mode[0], 5 * sizeof(uint8_t));
    for (i = 0; i < sbr->n_q; i++)
        ch_data->bs_invf_mode[0][i] = get_bits(gb, 2);
}

static void sbr_envelope(SpectralBandReplication *sbr, GetBitContext *gb,
                         SBRData *ch_data, int ch)
{
    int bits, max_depth;
    int i, j;
    VLC_TYPE (*t_huff)[2], (*f_huff)[2];
    int t_lav, f_lav;

    if (sbr->bs_coupling && ch) {
        max_depth = 2;
        if (ch_data->bs_amp_res) {
            bits = 5;
            t_huff = vlc_sbr[T_HUFFMAN_ENV_BAL_3_0DB].table;
            t_lav  = vlc_sbr_lav[T_HUFFMAN_ENV_BAL_3_0DB];
            f_huff = vlc_sbr[F_HUFFMAN_ENV_BAL_3_0DB].table;
            f_lav  = vlc_sbr_lav[F_HUFFMAN_ENV_BAL_3_0DB];
        } else {
            bits = 6;
            t_huff = vlc_sbr[T_HUFFMAN_ENV_BAL_1_5DB].table;
            t_lav  = vlc_sbr_lav[T_HUFFMAN_ENV_BAL_1_5DB];
            f_huff = vlc_sbr[F_HUFFMAN_ENV_BAL_1_5DB].table;
            f_lav  = vlc_sbr_lav[F_HUFFMAN_ENV_BAL_1_5DB];
        }
    } else {
        max_depth = 3;
        if (ch_data->bs_amp_res) {
            bits = 6;
            t_huff = vlc_sbr[T_HUFFMAN_ENV_3_0DB].table;
            t_lav  = vlc_sbr_lav[T_HUFFMAN_ENV_3_0DB];
            f_huff = vlc_sbr[F_HUFFMAN_ENV_3_0DB].table;
            f_lav  = vlc_sbr_lav[F_HUFFMAN_ENV_3_0DB];
        } else {
            bits = 7;
            t_huff = vlc_sbr[T_HUFFMAN_ENV_1_5DB].table;
            t_lav  = vlc_sbr_lav[T_HUFFMAN_ENV_1_5DB];
            f_huff = vlc_sbr[F_HUFFMAN_ENV_1_5DB].table;
            f_lav  = vlc_sbr_lav[F_HUFFMAN_ENV_1_5DB];
        }
    }

    for (i = 0; i < ch_data->bs_num_env[1]; i++) {
        if (!ch_data->bs_df_env[i]) {
            ch_data->bs_data_env[i][0] = get_bits(gb, bits); // bs_env_start_value_balance
            for (j = 1; j < sbr->n[ch_data->bs_freq_res[i + 1]]; j++)
                ch_data->bs_data_env[i][j] = get_vlc2(gb, f_huff, 9, max_depth) - f_lav;
        } else {
            for (j = 0; j < sbr->n[ch_data->bs_freq_res[i + 1]]; j++)
                ch_data->bs_data_env[i][j] = get_vlc2(gb, t_huff, 9, max_depth) - t_lav;
        }
    }
}

static void sbr_noise(SpectralBandReplication *sbr, GetBitContext *gb,
                      SBRData *ch_data, int ch)
{
    int max_depth;
    int i, j;
    VLC_TYPE (*t_huff)[2], (*f_huff)[2];
    int t_lav, f_lav;

    if (sbr->bs_coupling && ch) {
        max_depth = 1;
        t_huff = vlc_sbr[T_HUFFMAN_NOISE_BAL_3_0DB].table;
        t_lav  = vlc_sbr_lav[T_HUFFMAN_NOISE_BAL_3_0DB];
        f_huff = vlc_sbr[F_HUFFMAN_ENV_BAL_3_0DB].table;
        f_lav  = vlc_sbr_lav[F_HUFFMAN_ENV_BAL_3_0DB];
    } else {
        max_depth = 2;
        t_huff = vlc_sbr[T_HUFFMAN_NOISE_3_0DB].table;
        t_lav  = vlc_sbr_lav[T_HUFFMAN_NOISE_3_0DB];
        f_huff = vlc_sbr[F_HUFFMAN_ENV_3_0DB].table;
        f_lav  = vlc_sbr_lav[F_HUFFMAN_ENV_3_0DB];
    }

    for (i = 0; i < ch_data->bs_num_noise; i++) {
        if (!ch_data->bs_df_noise[i]) {
            ch_data->bs_data_noise[i][0] = get_bits(gb, 5); // bs_noise_start_value_balance or bs_noise_start_value_level
            for (j = 1; j < sbr->n_q; j++)
                ch_data->bs_data_noise[i][j] = get_vlc2(gb, f_huff, 9, max_depth + 1) - f_lav;
        } else {
            for (j = 0; j < sbr->n_q; j++)
                ch_data->bs_data_noise[i][j] = get_vlc2(gb, t_huff, 9, max_depth) - t_lav;
        }
    }
}

static void sbr_sinusoidal_coding(SpectralBandReplication *sbr,
                                  GetBitContext *gb, SBRData *ch_data)
{
    int i;
    for (i = 0; i < sbr->n[1]; i++)
        ch_data->bs_add_harmonic[i] = get_bits1(gb);
}

static void sbr_extension(SpectralBandReplication *sbr, GetBitContext *gb,
                          int bs_extension_id, int *num_bits_left)
{
/* FIXME - implement ps_data for parametric stereo parsing
    switch (bs_extension_id) {
    case EXTENSION_ID_PS:
        num_bits_left -= ps_data(sbr, gb);
        break;
    default:
*/
        skip_bits(gb, *num_bits_left); // bs_fill_bits
        *num_bits_left = 0;
/*
        break;
    }
*/
}

static void sbr_single_channel_element(AACContext *ac,
                                       SpectralBandReplication *sbr,
                                       GetBitContext *gb)
{
    if (get_bits1(gb)) // bs_data_extra
        skip_bits(gb, 4); // bs_reserved

    sbr_grid(ac, sbr, gb, &sbr->data[0]);
    sbr_dtdf(sbr, gb, &sbr->data[0]);
    sbr_invf(sbr, gb, &sbr->data[0]);
    sbr_envelope(sbr, gb, &sbr->data[0], 0);
    sbr_noise(sbr, gb, &sbr->data[0], 0);

    if ((sbr->data[0].bs_add_harmonic_flag = get_bits1(gb)))
        sbr_sinusoidal_coding(sbr, gb, &sbr->data[0]);
}

static void sbr_channel_pair_element(AACContext *ac,
                                     SpectralBandReplication *sbr,
                                     GetBitContext *gb)
{
    if (get_bits1(gb))    // bs_data_extra
        skip_bits(gb, 8); // bs_reserved

    if ((sbr->bs_coupling = get_bits1(gb))) {
        sbr_grid(ac, sbr, gb, &sbr->data[0]);
        memcpy(&sbr->data[1], &sbr->data[0], sizeof(sbr->data[1]));
        sbr_dtdf(sbr, gb, &sbr->data[0]);
        sbr_dtdf(sbr, gb, &sbr->data[1]);
        sbr_invf(sbr, gb, &sbr->data[0]);
        memcpy(sbr->data[1].bs_invf_mode, sbr->data[0].bs_invf_mode, sizeof(sbr->data[1].bs_invf_mode));
        sbr_envelope(sbr, gb, &sbr->data[0], 0);
        sbr_noise(sbr, gb, &sbr->data[0], 0);
        sbr_envelope(sbr, gb, &sbr->data[1], 1);
        sbr_noise(sbr, gb, &sbr->data[1], 1);
    } else {
        sbr_grid(ac, sbr, gb, &sbr->data[0]);
        sbr_grid(ac, sbr, gb, &sbr->data[1]);
        sbr_dtdf(sbr, gb, &sbr->data[0]);
        sbr_dtdf(sbr, gb, &sbr->data[1]);
        sbr_invf(sbr, gb, &sbr->data[0]);
        sbr_invf(sbr, gb, &sbr->data[1]);
        sbr_envelope(sbr, gb, &sbr->data[0], 0);
        sbr_envelope(sbr, gb, &sbr->data[1], 1);
        sbr_noise(sbr, gb, &sbr->data[0], 0);
        sbr_noise(sbr, gb, &sbr->data[1], 1);
    }

    if ((sbr->data[0].bs_add_harmonic_flag = get_bits1(gb)))
        sbr_sinusoidal_coding(sbr, gb, &sbr->data[0]);
    if ((sbr->data[1].bs_add_harmonic_flag = get_bits1(gb)))
        sbr_sinusoidal_coding(sbr, gb, &sbr->data[1]);
}

static unsigned int sbr_data(AACContext *ac, SpectralBandReplication *sbr,
                             GetBitContext *gb, int id_aac)
{
    unsigned int cnt = get_bits_count(gb);

    if (id_aac == TYPE_SCE || id_aac == TYPE_CCE) {
        sbr_single_channel_element(ac, sbr, gb);
    } else if (id_aac == TYPE_CPE) {
        sbr_channel_pair_element(ac, sbr, gb);
    } else {
        av_log(ac->avccontext, AV_LOG_ERROR,
            "Invalid bitstream - cannot apply SBR to element type %d\n", id_aac);
        return -1;
    }
    if (get_bits1(gb)) { // bs_extended_data
        int num_bits_left = get_bits(gb, 4); // bs_extension_size
        if (num_bits_left == 15)
            num_bits_left += get_bits(gb, 8); // bs_esc_count

        num_bits_left <<= 3;
        while (num_bits_left > 7) {
            num_bits_left -= 2;
            sbr_extension(sbr, gb, get_bits(gb, 2), &num_bits_left); // bs_extension_id
        }
    }

    return get_bits_count(gb) - cnt;
}

static void sbr_reset(AACContext *ac, SpectralBandReplication *sbr)
{
    int err;
    err = sbr_make_f_master(ac, sbr, &sbr->spectrum_params[1]);
    if (err >= 0)
        err = sbr_make_f_derived(ac, sbr);
    if (err < 0) {
        av_log(ac->avccontext, AV_LOG_ERROR, "SBR reset failed. Switching SBR to pure upsampling mode.\n");
        sbr->start = 0;
    }
}

/**
 * Decode Spectral Band Replication extension data; reference: table 4.55.
 *
 * @param   crc flag indicating the presence of CRC checksum
 * @param   cnt length of TYPE_FIL syntactic element in bytes
 *
 * @return  Returns number of bytes consumed from the TYPE_FIL element.
 */
int ff_decode_sbr_extension(AACContext *ac, SpectralBandReplication *sbr,
                            GetBitContext *gb_host, int crc, int cnt, int id_aac)
{
    unsigned int num_sbr_bits = 0, num_align_bits;
    unsigned bytes_read;
    GetBitContext gbc = *gb_host;
    GetBitContext *gb = &gbc;
    skip_bits_long(gb_host, cnt*8 - 4);

    sbr->reset = 0;

    if (!sbr->sample_rate)
        sbr->sample_rate = 2 * ac->m4ac.sample_rate; //TODO use the nominal sample rate for arbitrary sample rate support

    if (crc) {
        skip_bits(gb, 10); // bs_sbr_crc_bits; FIXME - implement CRC check
        num_sbr_bits += 10;
    }

    num_sbr_bits++;
    if (get_bits1(gb)) // bs_header_flag
        num_sbr_bits += sbr_header(sbr, gb);

    if (sbr->reset)
        sbr_reset(ac, sbr);

    if (sbr->start)
    num_sbr_bits  += sbr_data(ac, sbr, gb, id_aac);
    num_align_bits = ((cnt << 3) - 4 - num_sbr_bits) & 7;

    bytes_read = ((num_sbr_bits + num_align_bits + 4) / 8);
    if (bytes_read > cnt) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Expected to read %d SBR bytes actually read %d.\n", cnt, bytes_read);
    }
    return cnt;
}

// Time/frequency Grid (14496-3 sp04 p200)
static int sbr_time_freq_grid(AACContext *ac, SpectralBandReplication *sbr,
                              SBRData *ch_data, int ch)
{
    int abs_bord_lead  =  ch_data->bs_frame_class >= 2 ? ch_data->bs_var_bord[0] : 0;
    // frameLengthFlag ? 15 : 16; 960 sample length frames unsupported; this value is numTimeSlots
    int abs_bord_trail = (ch_data->bs_frame_class & 1 ? ch_data->bs_var_bord[1] : 0) + 16;
    int n_rel_lead, n_rel_trail;
    int i;

    if (ch_data->bs_frame_class == FIXFIX) {
        n_rel_lead = ch_data->bs_num_env[1] - 1;
    } else if (ch_data->bs_frame_class == FIXVAR) {
        n_rel_lead = 0;
    } else if (ch_data->bs_frame_class < 4) { // VARFIX or VARVAR
        n_rel_lead = ch_data->bs_num_rel[0];
    } else {
        av_log(ac->avccontext, AV_LOG_ERROR, "Invalid bs_frame_class for SBR: %d\n", ch_data->bs_frame_class);
        return -1;
    }

    n_rel_trail = ch_data->bs_frame_class & 1 ? ch_data->bs_num_rel[1] : 0;

    sbr->t_env_num_env_old[ch] = sbr->t_env[ch][ch_data->bs_num_env[0]]; //FIXME move me into a setup next frame area
    sbr->t_env[ch][0]                      = abs_bord_lead;
    sbr->t_env[ch][ch_data->bs_num_env[1]] = abs_bord_trail;

    if (ch_data->bs_frame_class == FIXFIX) {
        unsigned int temp = (unsigned int)lroundf(abs_bord_trail / (float)ch_data->bs_num_env[1]);
        for (i = 0; i < n_rel_lead; i++)
            sbr->t_env[ch][i + 1] = sbr->t_env[ch][i] + temp;
    } else if (ch_data->bs_frame_class > 1) { // VARFIX or VARVAR
        for (i = 0; i < n_rel_lead; i++)
            sbr->t_env[ch][i + 1] = sbr->t_env[ch][i] + ch_data->bs_rel_bord[0][i];
    } else { // FIXVAR
        for (i = 0; i < n_rel_lead; i++)
            sbr->t_env[ch][i + 1] = abs_bord_lead;
    }

    if (ch_data->bs_frame_class & 1) { // FIXVAR or VARVAR
        for (i = ch_data->bs_num_env[1] - 1; i > n_rel_lead; i--)
            sbr->t_env[ch][i] = sbr->t_env[ch][i + 1] - ch_data->bs_rel_bord[1][ch_data->bs_num_env[1] - 1 - i];
    } else { // FIXFIX or VARFIX
        for (i = n_rel_lead; i < ch_data->bs_num_env[1]; i++)
            sbr->t_env[ch][i + 1] = abs_bord_trail;
    }

    sbr->t_q[ch][0] = sbr->t_env[ch][0];
    if (ch_data->bs_num_noise > 1) { // typo in spec bases this on bs_num_env...
        unsigned int idx;
        if (ch_data->bs_frame_class == FIXFIX) {
            idx = ch_data->bs_num_env[1] >> 1;
        } else if (ch_data->bs_frame_class & 1) { // FIXVAR or VARVAR
            idx = ch_data->bs_num_env[1] - FFMAX(ch_data->bs_pointer - 1, 1);
        } else { // VARFIX
            if (!ch_data->bs_pointer)
                idx = 1;
            else if (ch_data->bs_pointer == 1)
                idx = ch_data->bs_num_env[1] - 1;
            else // bs_pointer > 1
                idx = ch_data->bs_pointer - 1;
        }
        sbr->t_q[ch][1] = sbr->t_env[ch][idx];
        sbr->t_q[ch][2] = sbr->t_env[ch][ch_data->bs_num_env[1]];
    } else
        sbr->t_q[ch][1] = sbr->t_env[ch][ch_data->bs_num_env[1]];

    return 0;
}

// SBR Envelope and Noise Floor Decoding (14496-3 sp04 p201)
static void sbr_env_noise_floors(SpectralBandReplication *sbr, SBRData *ch_data,
                                 int ch)
{
    int delta = (ch == 1 && sbr->bs_coupling == 1) ? 2 : 1;
    int i, k, l;
    const int temp = sbr->n[1] & 1;
    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        if (ch_data->bs_df_env[l]) {
            // bs_freq_res[0] == bs_freq_res[bs_num_env[1]] from prev frame
            if (ch_data->bs_freq_res[l + 1] == ch_data->bs_freq_res[l]) {
                for (k = 0; k < sbr->n[ch_data->bs_freq_res[l + 1]]; k++)
                    sbr->env_facs[ch][l + 1][k] = sbr->env_facs[ch][l][k] + delta * ch_data->bs_data_env[l][k];
            } else if (ch_data->bs_freq_res[l + 1]) {
                for (k = 0; k < sbr->n[ch_data->bs_freq_res[l + 1]]; k++) {
                    i = (k + temp) >> 1; // find i such that f_tablelow[i] <= f_tablehigh[k] < f_tablelow[i + 1]
                    sbr->env_facs[ch][l + 1][k] = sbr->env_facs[ch][l][i] + delta * ch_data->bs_data_env[l][k];
                }
            } else {
                for (k = 0; k < sbr->n[ch_data->bs_freq_res[l + 1]]; k++) {
                    i = k ? 2*k - temp : 0; // find i such that f_tablehigh[i] == f_tablelow[k]
                    sbr->env_facs[ch][l + 1][k] = sbr->env_facs[ch][l][i] + delta * ch_data->bs_data_env[l][k];
                }
            }
        } else {
            sbr->env_facs[ch][l + 1][0] = delta * ch_data->bs_data_env[l][0];
            for (k = 1; k < sbr->n[ch_data->bs_freq_res[l + 1]]; k++)
                sbr->env_facs[ch][l + 1][k] = sbr->env_facs[ch][l + 1][k - 1] + delta * ch_data->bs_data_env[l][k];
        }
    }

    for (l = 0; l < ch_data->bs_num_noise; l++) {
        if (ch_data->bs_df_noise[l])
            for (k = 0; k < sbr->n_q; k++)
                sbr->noise_facs[ch][l + 1][k] = sbr->noise_facs[ch][l][k] + delta * ch_data->bs_data_noise[l][k];
        else {
            sbr->noise_facs[ch][l + 1][0] = delta * ch_data->bs_data_noise[l][0];
            for (k = 1; k < sbr->n_q; k++)
                sbr->noise_facs[ch][l + 1][k] = sbr->noise_facs[ch][l + 1][k - 1] + delta * ch_data->bs_data_noise[l][k];
        }
    }

    //assign 0th elements of (env|noise)_facs from last elements
    memcpy(  sbr->env_facs[ch][0],   sbr->env_facs[ch][ch_data->bs_num_env[1]], sizeof(  sbr->env_facs[ch][0]));
    memcpy(sbr->noise_facs[ch][0], sbr->noise_facs[ch][ch_data->bs_num_noise ], sizeof(sbr->noise_facs[ch][0]));
}

// Dequantisation and stereo decoding (14496-3 sp04 p203)
static void sbr_dequant(SpectralBandReplication *sbr, int id_aac, int ch)
{
    int k, l;
    float alpha = sbr->data[ch].bs_amp_res ? 1.0f : 0.5f;

    if (id_aac == TYPE_CPE && sbr->bs_coupling) {
        float pan_offset = sbr->data[ch].bs_amp_res ? 12.0f : 24.0f;
        for (l = 1; l <= sbr->data[ch].bs_num_env[1]; l++) {
            for (k = 0; k < sbr->n[sbr->data[ch].bs_freq_res[l]]; k++) {
                float temp1 = powf(2.0f, sbr->env_facs[0][l][k] * alpha + 7.0f);
                float temp2 = (pan_offset - sbr->env_facs[1][l][k]) * alpha;
                sbr->env_facs[0][l][k] = temp1 / (1.0f + powf(2.0f,  temp2));
                sbr->env_facs[1][l][k] = temp1 / (1.0f + powf(2.0f, -temp2));
            }
        }
        for (l = 1; l <= sbr->data[ch].bs_num_noise; l++) {
            for (k = 0; k < sbr->n_q; k++) {
                float temp1 = powf(2.0f, NOISE_FLOOR_OFFSET - sbr->noise_facs[0][l][k] + 1);
                float temp2 = 12 - sbr->noise_facs[1][l][k];
                sbr->noise_facs[0][l][k] = temp1 / (1.0f + powf(2.0f,  temp2));
                sbr->noise_facs[1][l][k] = temp1 / (1.0f + powf(2.0f, -temp2));
            }
        }
    } else { // SCE or one non-coupled CPE
        for (l = 1; l <= sbr->data[ch].bs_num_env[1]; l++)
            for (k = 0; k < sbr->n[sbr->data[ch].bs_freq_res[l]]; k++)
                sbr->env_facs[ch][l][k] = powf(2.0f, alpha * sbr->env_facs[ch][l][k] + 6.0f);
        for (l = 1; l <= sbr->data[ch].bs_num_noise; l++)
            for (k = 0; k < sbr->n_q; k++)
                sbr->noise_facs[ch][l][k] = powf(2.0f, NOISE_FLOOR_OFFSET - sbr->noise_facs[ch][l][k]);
    }
}

  /**
 * Analysis QMF Bank (14496-3 sp04 p206)
 *
 * @param   x       pointer to the beginning of the first sample window
 * @param   W       array of complex-valued samples split into subbands
 */
static void sbr_qmf_analysis(const float *in, float *x, float W[2][32][32][2])
{
    int i, k, l, n;
    memcpy(W[1], W[0], sizeof(W[0]));
    memcpy(x    , x+1024, (320-32)*sizeof(x[0]));
    memcpy(x+288, in    ,     1024*sizeof(x[0]));
    x += 319;
    for (l = 0; l < 32; l++) { // 32 = numTimeSlots*RATE = 16*2 as 960 sample frames are not supported
        float z[320], u[64];
        for (i = 0; i < 320; i++)
            z[i] = x[-i] * sbr_qmf_window[i << 1];
        for (i = 0; i < 64; i++)
            u[i] = z[i] + z[i + 64] + z[i + 128] + z[i + 192] + z[i + 256];
        for (k = 0; k < 32; k++) {
            float temp1 = u[0] * 2.0f;
            float temp2 = -(k + 0.5f) * M_PI / 128.0f;
            W[0][k][l][0] = temp1 * cosf(temp2);
            W[0][k][l][1] = temp1 * sinf(temp2);
            for (n = 1; n < 64; n++) {
                temp1 = u[n] * 2.0f;
                temp2 = (n - 0.25f) * (k + 0.5f) * M_PI / 32.0f;
                W[0][k][l][0] += temp1 * cosf(temp2);
                W[0][k][l][1] += temp1 * sinf(temp2);
            }
        }
        x += 32;
    }
}

// Synthesis QMF Bank (14496-3 sp04 p206)
// Downsampled Synthesis QMF Bank (14496-3 sp04 p206)
static void sbr_qmf_synthesis(float *out, float X[64][40][2],
                              float *v, const unsigned int div)
{
    int k, l, n;
    float w[640];
    for (l = 0; l < 32; l++) {
        memmove(&v[128 / div], v, (1280 - 128) / div * sizeof(float));
        for (n = 0; n < 128 / div; n++) {
            v[n] = X[0][l][0] * cosf((2.0f * n - 255.0f / div) * M_PI / (256.0f / div)) -
                   X[0][l][1] * sinf((2.0f * n - 255.0f / div) * M_PI / (256.0f / div));
            for (k = 1; k < 64 / div; k++) {
                v[n] += X[k][l][0] * cosf((k + 0.5f) * (2.0f * n - 255.0f / div) * M_PI / (128.0f / div)) -
                        X[k][l][1] * sinf((k + 0.5f) * (2.0f * n - 255.0f / div) * M_PI / (128.0f / div));
            }
            v[n] /= 64.0f / div;
        }
        for (n = 0; n <= 4; n++) {
            int temp1 = 128 / div * n, temp2 = temp1 << 1;
            for (k = 0; k < 64 / div; k++) {
                w[temp1 + k]            = v[temp2 + k]             * sbr_qmf_window[temp1 + k];
                w[temp1 + k + 64 / div] = v[temp2 + k + 192 / div] * sbr_qmf_window[temp1 + k + 64 / div];
            }
        }
        for (k = 0; k < 64 / div; k++) {
            out[k] = w[k]             + w[64  / div + k] + w[128 / div + k] + w[192 / div + k] + w[256 / div + k]
                   + w[320 / div + k] + w[384 / div + k] + w[448 / div + k] + w[512 / div + k] + w[576 / div + k];
        }
        out += 64 / div;
    }
}

// High Frequency Generation (14496-3 sp04 p214+)

// Inverse Filtering (14496-3 sp04 p214)
static void sbr_hf_inverse_filter(float (*alpha0)[2], float (*alpha1)[2],
                                  float x_low[32][40][2], int k0)
{
    int i, j, k, n;
    for (k = 0; k < k0; k++) {
        float phi[3][2][2], dk;

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 2; j++) {
                unsigned int idxtmp1 = ENVELOPE_ADJUSTMENT_OFFSET - i;
                unsigned int idxtmp2 = ENVELOPE_ADJUSTMENT_OFFSET - (j + 1);

                phi[i][j][0] = 0.0f;
                phi[i][j][1] = 0.0f;

                if (i <= j + 1)
                for (n = 0; n < 16 * 2 + 6; n++) {
                    unsigned int idx1 = n + idxtmp1;
                    unsigned int idx2 = n + idxtmp2;
                    phi[i][j][0] += x_low[k][idx1][0] * x_low[k][idx2][0] +
                                    x_low[k][idx1][1] * x_low[k][idx2][1];
                    if (i != j + 1) { // imaginary part
                        phi[i][j][1] += x_low[k][idx1][1] * x_low[k][idx2][0] -
                                        x_low[k][idx1][0] * x_low[k][idx2][1];
                    }
                }
            }
        }

        dk = phi[2][1][0] * phi[1][0][0] -
               (phi[1][1][0] * phi[1][1][0] + phi[1][1][1] * phi[1][1][1]) / 1.000001f;

        if (!dk) {
            alpha1[k][0] = 0;
            alpha1[k][1] = 0;
        } else {
            float temp_real, temp_im;
            temp_real = phi[0][0][0] * phi[1][1][0] -
                        phi[0][0][1] * phi[1][1][1] -
                        phi[0][1][0] * phi[1][0][0];
            temp_im   = phi[0][0][0] * phi[1][1][1] +
                        phi[0][0][1] * phi[1][1][0] -
                        phi[0][1][1] * phi[1][0][0];

            alpha1[k][0] = temp_real / dk;
            alpha1[k][1] = temp_im   / dk;
        }

        if (!phi[1][0][0]) {
            alpha0[k][0] = 0;
            alpha0[k][1] = 0;
        } else {
            float temp_real, temp_im;
            temp_real = phi[0][0][0] + alpha1[k][0] * phi[1][1][0] +
                                       alpha1[k][1] * phi[1][1][1];
            temp_im   = phi[0][0][1] + alpha1[k][1] * phi[1][1][0] -
                                       alpha1[k][0] * phi[1][1][1];

            alpha0[k][0] = -temp_real / phi[1][0][0];
            alpha0[k][0] = -temp_im   / phi[1][0][0];
        }

        if (alpha1[k][0] * alpha1[k][0] + alpha1[k][1] * alpha1[k][1] >= 16.0f ||
           alpha0[k][0] * alpha0[k][0] + alpha0[k][1] * alpha0[k][1] >= 16.0f) {
            alpha1[k][0] = 0;
            alpha1[k][1] = 0;
            alpha0[k][0] = 0;
            alpha0[k][1] = 0;
        }
    }
}

// Chirp Factors (14496-3 sp04 p214)
static void sbr_chirp(SpectralBandReplication *sbr, SBRData *ch_data)
{
    int i;
    float new_bw;
    float temp_bw;

    for (i = 0; i < sbr->n_q; i++) {
        switch (ch_data->bs_invf_mode[0][i]) {
        case 0:
            if (ch_data->bs_invf_mode[1][i] == 1) {
                new_bw = 0.6f;
            } else
                new_bw = 0.0f;
            break;
        case 1:
            if (!ch_data->bs_invf_mode[1][i]) {
                new_bw = 0.6f;
            } else
                new_bw = 0.75f;
            break;
        case 2:
            new_bw = 0.9f;
            break;
        case 3:
            new_bw = 0.98f;
            break;
        default:
            break;
        }

        if (new_bw < sbr->bw_array[1][i]) {
            temp_bw = 0.75f    * new_bw + 0.25f    * sbr->bw_array[1][i];
        } else
            temp_bw = 0.90625f * new_bw + 0.09375f * sbr->bw_array[1][i];
        sbr->bw_array[0][i] = temp_bw < 0.015625f ? 0.0f : temp_bw;
    }

    // update previous bw_array values
    memcpy(sbr->bw_array[1], sbr->bw_array[0], 5 * sizeof(float));
}

static inline int find_freq_subband(uint16_t *table, int nel, int needle)
{
    int i;
    for (i = 0; i < nel; i++) {
        if (needle >= table[i] && needle < table[i + 1])
            return i;
    }
    return -1;
}

// Generate the subband filtered lowband
static int sbr_lf_gen(AACContext *ac, SpectralBandReplication *sbr,
                      float x_low[32][40][2], float W[2][32][32][2]) {
    int k, l;
    const int t_HFGen = 8;
    const int l_f = 32;
    memset(x_low, 0, 32*sizeof(*x_low));
    for (k = 0; k < sbr->k[3]; k++) {
        for (l = t_HFGen; l < l_f + t_HFGen; l++) {
            x_low[k][l][0] = W[0][k][l - t_HFGen][0];
            x_low[k][l][1] = W[0][k][l - t_HFGen][1];
        }
    }
    for (k = 0; k < sbr->k[4]; k++) {
        for (l = 0; l < t_HFGen; l++) {
            x_low[k][l][0] = W[1][k][l + l_f - t_HFGen][0];
            x_low[k][l][1] = W[1][k][l + l_f - t_HFGen][1];
        }
    }
    return 0;
}

// High Frequency Generator (14496-3 sp04 p215)
static int sbr_hf_gen(AACContext *ac, SpectralBandReplication *sbr,
                      float x_high[64][40][2], float x_low[32][40][2], float (*alpha0)[2],
                      float (*alpha1)[2], float bw_array[2][5], uint8_t *t_env,
                      int bs_num_env)
{
    int i, x, l;
    int ktmp = sbr->k[3];
    for (i = 0; i < sbr->num_patches; i++) {
        if (i >= 1)
            ktmp += sbr->patch_num_subbands[i-1];
        for (x = 0; x < sbr->patch_num_subbands[i]; x++) {
            const int k = ktmp + x;
            const int g = find_freq_subband(sbr->f_tablenoise, sbr->n_q + 1, k);
            const int p = sbr->patch_start_subband[i] + x;

            if (g < 0) {
                av_log(ac->avccontext, AV_LOG_ERROR, "ERROR : no subband found for frequency %d\n", k);
                return -1;
            }

            for (l = t_env[0] << 1; l < t_env[bs_num_env] << 1; l++) {
                const int idx = l + ENVELOPE_ADJUSTMENT_OFFSET;
                x_high[k][idx][0] =
                    (x_low[p][idx - 2][0] * alpha1[p][0] -
                     x_low[p][idx - 2][1] * alpha1[p][1]) * bw_array[0][g] * bw_array[0][g] +
                    (x_low[p][idx - 1][0] * alpha0[p][0] -
                     x_low[p][idx - 1][1] * alpha0[p][1]) * bw_array[0][g] +
                     x_low[p][idx][0];
                x_high[k][idx][1] =
                    (x_low[p][idx - 2][1] * alpha1[p][0] +
                     x_low[p][idx - 2][0] * alpha1[p][1]) * bw_array[0][g] * bw_array[0][g] +
                    (x_low[p][idx - 1][1] * alpha0[p][0] +
                     x_low[p][idx - 1][0] * alpha0[p][1]) * bw_array[0][g] +
                     x_low[p][idx][1];
            }
        }
    }

    return 0;
}

// Generate the subband filtered lowband
static int sbr_x_gen(SpectralBandReplication *sbr,
                      float x[64][40][2], float x_low[32][40][2], float Y[2][64][40][2], int ch) {
    int k, l;
    const int t_HFAdj = ENVELOPE_ADJUSTMENT_OFFSET;
    const int l_f = 32;
    const int l_Temp = FFMAX(2*sbr->t_env_num_env_old[ch] - l_f, 0); //FIXME hack to make l_Temp initialize to zero
    memset(x, 0, 64*sizeof(*x));
    for (k = 0; k < sbr->k[4]; k++) {
        for (l = 0; l < l_Temp; l++) {
            x[k][l][0] = x_low[k][l + t_HFAdj][0];
            x[k][l][1] = x_low[k][l + t_HFAdj][1];
        }
    }
    for (; k < sbr->k[4] + sbr->mold; k++) {
        for (l = 0; l < l_Temp; l++) {
            x[k][l][0] = Y[1][k][l + t_HFAdj + l_f][0];
            x[k][l][1] = Y[1][k][l + t_HFAdj + l_f][1];
        }
    }

    for (k = 0; k < sbr->k[3]; k++) {
        for (l = l_Temp; l < l_f; l++) {
            x[k][l][0] = x_low[k][l + t_HFAdj][0];
            x[k][l][1] = x_low[k][l + t_HFAdj][1];
        }
    }
    for (; k < sbr->k[3] + sbr->m; k++) {
        for (l = l_Temp; l < l_f; l++) {
            x[k][l][0] = Y[0][k][l + t_HFAdj][0];
            x[k][l][1] = Y[0][k][l + t_HFAdj][1];
        }
    }
    return 0;
}

// High Frequency Adjustment (14496-3 sp04 p217)

// Mapping (14496-3 sp04 p217)
static void sbr_mapping(AACContext *ac, SpectralBandReplication *sbr,
                        SBRData *ch_data, int ch, int l_a[2])
{
    int i, l, m;

    // The following is used for
    l_a[0] = l_a[1]; // update previous frame's l_a value
    l_a[1] = -1;
    if ((ch_data->bs_frame_class & 1) && ch_data->bs_pointer) { // FIXVAR or VARVAR and bs_pointer != 0
        l_a[1] = ch_data->bs_num_env[1] + 1 - ch_data->bs_pointer;
    } else if ((ch_data->bs_frame_class == 2) && (ch_data->bs_pointer > 1)) // VARFIX and bs_pointer > 1
        l_a[1] = ch_data->bs_pointer - 1;

    memset(sbr->s_indexmapped[1], 0, 7*sizeof(sbr->s_indexmapped[1]));
    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        const unsigned int ilim = sbr->n[ch_data->bs_freq_res[l + 1]];
        uint16_t *table = ch_data->bs_freq_res[l + 1] ? sbr->f_tablehigh : sbr->f_tablelow;
        int k;

        for (i = 0; i < ilim; i++)
            for (m = table[i]; m < table[i + 1]; m++)
                sbr->e_origmapped[l][m - sbr->k[3]] = sbr->env_facs[ch][l+1][i];

        // ch_data->bs_num_noise > 1 => 2 noise floors
        k = (ch_data->bs_num_noise > 1) && (sbr->t_env[ch][l] >= sbr->t_q[ch][1]);
        for (i = 0; i < sbr->n_q; i++)
            for (m = sbr->f_tablenoise[i]; m < sbr->f_tablenoise[i + 1]; m++)
                sbr->q_mapped[l][m - sbr->k[3]] = sbr->noise_facs[ch][k+1][i];

        for (i = 0; i < sbr->n[1]; i++) {
            if (ch_data->bs_add_harmonic_flag) {
                const unsigned int m_midpoint =
                    (sbr->f_tablehigh[i] + sbr->f_tablehigh[i + 1]) >> 1;

                sbr->s_indexmapped[l + 1][m_midpoint - sbr->k[3]] = ch_data->bs_add_harmonic[i] *
                    (l >= l_a[1] || (sbr->s_indexmapped[0][m_midpoint - sbr->k[3]] == 1));
            }
        }

        for (i = 0; i < ilim; i++) {
            int additional_sinusoid_present = 0;
            for (m = table[i]; m < table[i + 1]; m++) {
                if (sbr->s_indexmapped[l + 1][m - sbr->k[3]]) {
                    additional_sinusoid_present = 1;
                    break;
                }
            }
            memset(&sbr->s_mapped[l][table[i] - sbr->k[3]], additional_sinusoid_present,
                   (table[i + 1] - table[i]) * sizeof(sbr->s_mapped[l][0]));
        }
    }

    memcpy(sbr->s_indexmapped[0], sbr->s_indexmapped[ch_data->bs_num_env[1]], sizeof(sbr->s_indexmapped[0]));
}

// Estimation of current envelope (14496-3 sp04 p218)
static void sbr_env_estimate(float (*e_curr)[48], float x_high[64][40][2],
                             SpectralBandReplication *sbr, SBRData *ch_data,
                             int ch)
{
    int i, l, m;

    if (sbr->bs_interpol_freq) {
        for (l = 0; l < ch_data->bs_num_env[1]; l++) {
            const int env_size = (sbr->t_env[ch][l + 1] - sbr->t_env[ch][l]) << 1;
            int ilb = sbr->t_env[ch][l]     * 2 + ENVELOPE_ADJUSTMENT_OFFSET;
            int iub = sbr->t_env[ch][l + 1] * 2 + ENVELOPE_ADJUSTMENT_OFFSET;

            for (m = 0; m < sbr->m; m++) {
                float sum = 0.0f;

                for (i = ilb; i < iub; i++) {
                    sum += x_high[m + sbr->k[3]][i][0] * x_high[m + sbr->k[3]][i][0] +
                           x_high[m + sbr->k[3]][i][1] * x_high[m + sbr->k[3]][i][1];
                }
                e_curr[l][m] = sum / env_size;
            }
        }
    } else {
        int k, p;

        for (l = 0; l < ch_data->bs_num_env[1]; l++) {
            const int env_size = (sbr->t_env[ch][l + 1] - sbr->t_env[ch][l]) << 1;
            int ilb = sbr->t_env[ch][l]     * 2 + ENVELOPE_ADJUSTMENT_OFFSET;
            int iub = sbr->t_env[ch][l + 1] * 2 + ENVELOPE_ADJUSTMENT_OFFSET;
            const uint16_t *table = ch_data->bs_freq_res[l + 1] ? sbr->f_tablehigh : sbr->f_tablelow;

            for (p = 0; p < sbr->n[ch_data->bs_freq_res[l + 1]]; p++) {
                float sum = 0.0f;
                const int den = env_size * (table[p + 1] - table[p]);

                for (k = table[p]; k < table[p + 1]; k++) {
                    for (i = ilb; i < iub; i++) {
                        sum += x_high[k][i][0] * x_high[k][i][0] +
                               x_high[k][i][1] * x_high[k][i][1];
                    }
                }
                sum /= den;
                for (k = table[p]; k < table[p + 1]; k++) {
                    e_curr[l][k - sbr->k[3]] = sum;
                }
            }
        }
    }
}

// Calculation of levels of additional HF signal components (14496-3 sp04 p219)
static void sbr_hf_additional_levels(SpectralBandReplication *sbr,
                                     SBRData *ch_data)
{
    int l, m;

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (m = 0; m < sbr->m; m++) {
            const float temp = sbr->e_origmapped[l][m] / (1.0f + sbr->q_mapped[l][m]);
            sbr->q_m[l][m] = sqrtf(temp * sbr->q_mapped[l][m]);
            sbr->s_m[l][m] = sqrtf(temp * sbr->s_indexmapped[l + 1][m]);
        }
    }
}

// Calculation of gain (14496-3 sp04 p219)
static void sbr_gain_calc(AACContext * ac, SpectralBandReplication *sbr,
                          SBRData *ch_data, int l_a[2])
{
    int i, k, l, m;
    float gain_boost_temp[7][48];
    float gain_max_temp[7][48];
    // max gain limits : -3dB, 0dB, 3dB, inf dB (limiter off)
    const float limgain[4] = { 0.70795, 1.0, 1.41254, 10000000000 };

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        int delta = !((l == l_a[1]) || (l == -(l_a[0] != ch_data->bs_num_env[0])));
        for (m = 0; m < sbr->m; m++) {
            if (!sbr->s_mapped[l][m]) {
                sbr->gain[l][m] = sqrtf(sbr->e_origmapped[l][m] /
                                        ((1.0f + sbr->e_curr[l][m]) *
                                         (1.0f + sbr->q_mapped[l][m] * delta)));
            } else {
                sbr->gain[l][m] = sqrtf(sbr->e_origmapped[l][m] * sbr->q_mapped[l][m] /
                                        ((1.0f + sbr->e_curr[l][m]) *
                                         (1.0f + sbr->q_mapped[l][m])));
            }
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (k = 0; k < sbr->n_lim; k++) {
            float sum[2] = { 0.0f, 0.0f };
            for (i = sbr->f_tablelim[k] - sbr->k[3]; i < sbr->f_tablelim[k + 1] - sbr->k[3]; i++) {
                sum[0] += sbr->e_origmapped[l][i];
                sum[1] += sbr->e_curr[l][i];
            }
            gain_max_temp[l][k] = FFMIN(100000,
                                        limgain[sbr->bs_limiter_gains] * sqrtf((EPS0 + sum[0]) / (EPS0 + sum[1])));
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (m = 0; m < sbr->m; m++) {
            if ((k = find_freq_subband(sbr->f_tablelim, sbr->n_lim, m + sbr->k[3])) < 0) {
                av_log(ac->avccontext, AV_LOG_ERROR,
                       "No subband found for frequency %d\n", m + sbr->k[3]);
            }
            sbr->gain_max[l][m] = gain_max_temp[l][k];
            sbr->q_m_lim[l][m]  = FFMIN(sbr->q_m[l][m],  sbr->q_m[l][m] * sbr->gain_max[l][m] / sbr->gain[l][m]);
            sbr->gain_lim[l][m] = FFMIN(sbr->gain[l][m], sbr->gain_max[l][m]);
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        int delta = !((l == l_a[1]) || (l == -(l_a[0] != ch_data->bs_num_env[0])));
        for (k = 0; k < sbr->n_lim; k++) {
            float sum[2] = { 0.0f, 0.0f };
            for (i = sbr->f_tablelim[k] - sbr->k[3]; i < sbr->f_tablelim[k + 1] - sbr->k[3]; i++) {
                sum[0] += sbr->e_origmapped[l][i];
                sum[1] += sbr->e_curr[l][i] * sbr->gain_lim[l][i] * sbr->gain_lim[l][i]
                          + sbr->s_m[l][i] * sbr->s_m[l][i]
                          + (delta && !sbr->s_m[l][i]) * sbr->q_m_lim[l][i] * sbr->q_m_lim[l][i];
            }
            gain_boost_temp[l][k] = FFMIN(1.584893192, sqrtf((EPS0 + sum[0]) / (EPS0 + sum[1])));
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (m = 0; m < sbr->m; m++) {
            if ((k = find_freq_subband(sbr->f_tablelim, sbr->n_lim, m + sbr->k[3])) < 0) {
                av_log(ac->avccontext, AV_LOG_ERROR,
                       "No subband found for frequency %d\n", m + sbr->k[3]);
            }
            sbr->gain_boost[l][m] = gain_boost_temp[l][k];
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (m = 0; m < sbr->m; m++) {
            sbr->gain_limboost[l][m] = sbr->gain_lim[l][m] * sbr->gain_boost[l][m];
            sbr->q_m_limboost[l][m]  = sbr->q_m_lim[l][m]  * sbr->gain_boost[l][m];
            sbr->s_m_boost[l][m]     = sbr->s_m[l][m]      * sbr->gain_boost[l][m];
        }
    }
}

// Assembling HF Signals (14496-3 sp04 p220)
static void sbr_hf_assemble(float y[2][64][40][2], float x_high[64][40][2],
                            SpectralBandReplication *sbr, SBRData *ch_data,
                            int ch, int l_a[2])
{
    int i, j, l, m;
    const int h_SL = sbr->bs_smoothing_mode ? 0 : 4;
    const float h_smooth[5] = {
        0.33333333333333,
        0.30150283239582,
        0.21816949906249,
        0.11516383427084,
        0.03183050093751,
    };
    const int8_t phi[2][4] = {
        {  1,  0, -1,  0}, // real
        {  0,  1,  0, -1}, // imaginary
    };
    float g_filt[42][48], q_filt[42][48], w_temp[42][48][2];
    float (*g_temp)[48] = sbr->g_temp, (*q_temp)[48] = sbr->q_temp;
    memcpy(y[1], y[0], sizeof(y[0]));

    if (sbr->reset) {
        for (i = 0; i < h_SL; i++) {
            memcpy(g_temp[i], sbr->gain_limboost[0], sbr->m * sizeof(sbr->gain_limboost[0][0]));
            memcpy(q_temp[i], sbr->q_m_limboost[0],  sbr->m * sizeof(sbr->q_m_limboost[0][0]));
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (i = sbr->t_env[ch][l] << 1; i < sbr->t_env[ch][l + 1] << 1; i++) {
            memcpy(g_temp[h_SL + i], sbr->gain_limboost[l], sbr->m * sizeof(sbr->gain_limboost[l][0]));
            memcpy(q_temp[h_SL + i], sbr->q_m_limboost[l],  sbr->m * sizeof(sbr->q_m_limboost[l][0]));
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        if (h_SL && l != l_a[0] && l != l_a[1]) {
            for (i = sbr->t_env[ch][l] << 1; i < sbr->t_env[ch][l + 1] << 1; i++) {
                for (m = 0; m < sbr->m; m++) {
                    const int idx1 = i + h_SL;
                    g_filt[i][m] = 0.0f;
                    for (j = 0; j <= h_SL; j++)
                        g_filt[i][m] += g_temp[idx1 - j][m] * h_smooth[j];
                }
            }
        } else {
            for (i = sbr->t_env[ch][l] << 1; i < sbr->t_env[ch][l + 1] << 1; i++)
                memcpy(g_filt[i], g_temp[i + h_SL], sbr->m * sizeof(g_temp[i + h_SL][0]));
        }
    }

    for (i = sbr->t_env[ch][0] << 1; i < sbr->t_env[ch][ch_data->bs_num_env[1]] << 1; i++) {
        const int idx2 = i + ENVELOPE_ADJUSTMENT_OFFSET;
        for (m = 0; m < sbr->m; m++) {
            const int idx1 = m + sbr->k[3];
            w_temp[i][m][0] = x_high[idx1][idx2][0] * g_filt[i][m];
            w_temp[i][m][1] = x_high[idx1][idx2][1] * g_filt[i][m];
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        if (l != l_a[0] && l != l_a[1]) {
            for (i = sbr->t_env[ch][l] << 1; i < sbr->t_env[ch][l + 1] << 1; i++) {
                for (m = 0; m < sbr->m; m++) {
                    if (sbr->s_m_boost[l][m])
                        q_filt[i][m] = 0.0f;
                    else if (h_SL) {
                        const int idx1 = i + h_SL;
                        q_filt[i][m] = 0.0f;
                        for (j = 0; j <= h_SL; j++)
                            q_filt[i][m] += q_temp[idx1 - j][m] * h_smooth[j];
                    } else
                        q_filt[i][m] = q_temp[i][m];
                }
            }
        } else {
            for (i = sbr->t_env[ch][l] << 1; i < sbr->t_env[ch][l + 1] << 1; i++)
                memset(q_filt[i], 0, sbr->m * sizeof(q_filt[i][0]));
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (i = sbr->t_env[ch][l] << 1; i < sbr->t_env[ch][l + 1] << 1; i++) {
            for (m = 0; m < sbr->m; m++) {
                sbr->f_indexnoise[ch] = (sbr->f_indexnoise[ch] + 1) & 0x1ff;
                w_temp[i][m][0] += q_filt[i][m] * sbr_noise_table[sbr->f_indexnoise[ch]][0];
                w_temp[i][m][1] += q_filt[i][m] * sbr_noise_table[sbr->f_indexnoise[ch]][1];
            }
        }
    }

    for (l = 0; l < ch_data->bs_num_env[1]; l++) {
        for (i = sbr->t_env[ch][l] << 1; i < sbr->t_env[ch][l + 1] << 1; i++) {
            sbr->f_indexsine[ch] = (sbr->f_indexsine[ch] + 1) & 3;
            for (m = 0; m < sbr->m; m++) {
                y[0][m + sbr->k[3]][i + ENVELOPE_ADJUSTMENT_OFFSET][0] =
                    w_temp[i][m][0] + sbr->s_m_boost[l][m] * phi[0][sbr->f_indexsine[ch]];
                y[0][m + sbr->k[3]][i + ENVELOPE_ADJUSTMENT_OFFSET][1] =
                    w_temp[i][m][1] + sbr->s_m_boost[l][m] * phi[1][sbr->f_indexsine[ch]] * (1 - 2*((m + sbr->k[3]) & 1));
            }
        }
    }

   memcpy(g_temp[0], g_temp[2*sbr->t_env[ch][ch_data->bs_num_env[1]] - 4], 4*sizeof(g_temp[0]));
   memcpy(q_temp[0], q_temp[2*sbr->t_env[ch][ch_data->bs_num_env[1]] - 4], 4*sizeof(q_temp[0]));
}

void ff_sbr_apply(AACContext *ac, SpectralBandReplication *sbr, int id_aac, int ch, float* in, float* out)
{
    int* l_a = sbr->data[ch].l_a;

    if (sbr->start) {
        sbr_time_freq_grid(ac, sbr, &sbr->data[ch], ch);
        sbr_env_noise_floors(sbr, &sbr->data[ch], ch);
        sbr_dequant(sbr, id_aac, ch);
    }

    /* decode channel */
    sbr_qmf_analysis(in, sbr->data[ch].analysis_filterbank_samples, sbr->W);
    sbr_lf_gen(ac, sbr, sbr->x_low, sbr->W);
    if (sbr->start) {
        sbr_hf_inverse_filter(sbr->alpha0, sbr->alpha1, sbr->x_low, sbr->k[0]);
        sbr_chirp(sbr, &sbr->data[ch]);
        sbr_hf_gen(ac, sbr, sbr->x_high, sbr->x_low, sbr->alpha0, sbr->alpha1,
                   sbr->bw_array, sbr->t_env[ch], sbr->data[ch].bs_num_env[1]);

    // hf_adj
        sbr_mapping(ac, sbr, &sbr->data[ch], ch, l_a);
        sbr_env_estimate(sbr->e_curr, sbr->x_high, sbr, &sbr->data[ch], ch);
        sbr_hf_additional_levels(sbr, &sbr->data[ch]);
        sbr_gain_calc(ac, sbr, &sbr->data[ch], l_a);
        sbr_hf_assemble(sbr->y, sbr->x_high, sbr, &sbr->data[ch], ch, l_a);
    }

    /* synthesis */
    sbr_x_gen(sbr, sbr->X, sbr->x_low, sbr->y, ch);
    sbr_qmf_synthesis(out, sbr->X, sbr->data[ch].synthesis_filterbank_samples, 1);
}
