/*
 * AAC Spectral Band Replication definitions and structures
 * Copyright (c) 2008-2009 Robert Swain ( rob opendot cl )
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
 * AAC Spectral Band Replication definitions and structures
 * @author Robert Swain ( rob opendot cl )
 */

#ifndef AVCODEC_AACSBR_H
#define AVCODEC_AACSBR_H

#include <stdint.h>

#define SBR_INIT_VLC_STATIC(num, size) \
    INIT_VLC_STATIC(&vlc_sbr[num], 9, sbr_tmp[num].table_size / sbr_tmp[num].elem_size,     \
                    sbr_tmp[num].sbr_bits ,                      1,                      1, \
                    sbr_tmp[num].sbr_codes, sbr_tmp[num].elem_size, sbr_tmp[num].elem_size, \
                    size);

#define NOISE_FLOOR_OFFSET 6.0f
#define ENVELOPE_ADJUSTMENT_OFFSET 2
// constant to avoid division by zero, e.g. 96 dB below maximum signal input
#define EPS0 0.000000000001

/**
 * SBR VLC tables
 */
enum {
    T_HUFFMAN_ENV_1_5DB,
    F_HUFFMAN_ENV_1_5DB,
    T_HUFFMAN_ENV_BAL_1_5DB,
    F_HUFFMAN_ENV_BAL_1_5DB,
    T_HUFFMAN_ENV_3_0DB,
    F_HUFFMAN_ENV_3_0DB,
    T_HUFFMAN_ENV_BAL_3_0DB,
    F_HUFFMAN_ENV_BAL_3_0DB,
    T_HUFFMAN_NOISE_3_0DB,
    T_HUFFMAN_NOISE_BAL_3_0DB,
};

/**
 * bs_frame_class - frame class of current SBR frame (14496-3 sp04 p98)
 */
enum {
    FIXFIX,
    FIXVAR,
    VARFIX,
    VARVAR,
};

/**
 * Spectral Band Replication header - spectrum parameters that invoke a reset if they differ from the previous header.
 */
typedef struct {
    uint8_t bs_start_freq;
    uint8_t bs_stop_freq;
    uint8_t bs_xover_band;

    // if (bs_header_extra_1)
    uint8_t bs_freq_scale;
    uint8_t bs_alter_scale;
    uint8_t bs_noise_bands;
} SpectrumParameters;

/**
 * Spectral Band Replication data - main bitstream data variables
 */
typedef struct {
    uint8_t            bs_frame_class;
    uint8_t            bs_add_harmonic_flag;
    uint8_t            bs_num_env[2];
    uint8_t            bs_freq_res[7];
    uint8_t            bs_var_bord[2];
    uint8_t            bs_num_rel[2];
    uint8_t            bs_rel_bord[2][3];
    uint8_t            bs_pointer;
    uint8_t            bs_num_noise;
    uint8_t            bs_df_env[5];
    uint8_t            bs_df_noise[2];
    uint8_t            bs_invf_mode[2][5];
    int32_t            bs_data_env[7][32];
    int32_t            bs_data_noise[2][5];
    uint8_t            bs_add_harmonic[32];
    uint8_t            bs_amp_res;
    float              synthesis_filterbank_samples[1280];
    float              analysis_filterbank_samples [1312];
    int                l_a[2];
} SBRData;

/**
 * Spectral Band Replication
 */
typedef struct {
    uint8_t            start;
    int32_t            sample_rate;
    // SBR header bitstream variables
    uint8_t            bs_amp_res_header;
    SpectrumParameters spectrum_params[2];
    // if (bs_header_extra_2)
    uint8_t            bs_limiter_bands;
    uint8_t            bs_limiter_gains;
    uint8_t            bs_interpol_freq;
    uint8_t            bs_smoothing_mode;
    // SBR data bitstream variables
    uint8_t            bs_coupling;
    SBRData            data[2];
    // SBR tool variables
    uint8_t            reset;
    uint8_t            k[5]; // k0, k1, k2, kx, and kx' respectively
    uint8_t            m;
    uint8_t            mold;
    uint8_t            n_master;
    uint8_t            n[2]; // n_low and n_high respectively
    uint8_t            n_q;
    uint8_t            n_lim;
    uint16_t           f_master[49];
    uint16_t           f_tablelow[25];
    uint16_t           f_tablehigh[49];
    uint16_t           f_tablenoise[6];
    uint16_t           f_tablelim[29];
    uint8_t            num_patches;
    uint8_t            patch_num_subbands[5];
    uint8_t            patch_start_subband[5];
    uint8_t            t_env[2][8];
    uint8_t            t_env_num_env_old[2];
    uint8_t            t_q[2][3];
    float              env_facs[2][6][48];
    float              noise_facs[2][3][5];
    float              W[2][32][32][2];
    float              x_low[32][40][2];
    float              x_high[64][40][2];
    float              y[2][64][40][2];
    float              X[64][40][2];
    float              alpha0[64][2];
    float              alpha1[64][2];
    float              bw_array[2][5];
    float              e_origmapped[7][48];
    float              q_mapped[7][48];
    uint8_t            s_indexmapped[8][48];
    uint8_t            s_mapped[7][48];
    float              e_curr[7][48];
    float              q_m[7][48];
    float              s_m[7][48];
    float              gain[7][48];
    float              gain_max[7][48];
    float              gain_lim[7][48];
    float              gain_boost[7][48];
    float              gain_limboost[7][48];
    float              q_m_lim[7][48];
    float              q_m_limboost[7][48];
    float              s_m_boost[7][48];
    float              g_temp[42][48];
    float              q_temp[42][48];
    uint16_t           f_indexnoise[2];
    uint8_t            f_indexsine[2];
} SpectralBandReplication;

#endif /* AVCODEC_AACSBR_H */
