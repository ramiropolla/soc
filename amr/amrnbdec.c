/*
 * AMR narrowband decoder
 * Copyright (c) 2006-2007 Robert Swain
 * Copyright (c) 2009 Colin McQuillan
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
 * @file libavcodec/amrnbdec.c
 * AMR narrowband decoder
 */


#include <string.h>
#include <math.h>

#include "avcodec.h"
#include "get_bits.h"
#include "libavutil/common.h"
#include "celp_math.h"
#include "celp_filters.h"
#include "acelp_filters.h"
#include "acelp_vectors.h"
#include "lsp.h"

#include "amrnbdata.h"

typedef struct AMRContext {
    AMRNBFrame                        frame; ///< decoded AMR parameters (lsf coefficients, codebook indexes, etc)
    uint8_t             bad_frame_indicator; ///< bad frame ? 1 : 0
    enum Mode                cur_frame_mode; ///< current frame mode

    float       prev_lsf_r[LP_FILTER_ORDER]; ///< residual LSF vector from previous subframe
    float           lsp[4][LP_FILTER_ORDER]; ///< lsp vectors from current frame
    float    prev_lsp_sub4[LP_FILTER_ORDER]; ///< lsp vector for the 4th subframe of the previous frame

    float         lsf_q[4][LP_FILTER_ORDER]; ///< Interpolated LSF vector for fixed gain smoothing
    float          lsf_avg[LP_FILTER_ORDER]; ///< vector of averaged lsf vector

    float           lpc[4][LP_FILTER_ORDER]; ///< lpc coefficient vectors for 4 subframes

    uint8_t                   pitch_lag_int; ///< integer part of pitch lag from current subframe

    float excitation_buf[PITCH_LAG_MAX + LP_FILTER_ORDER + 1 + AMR_SUBFRAME_SIZE]; ///< excitation buffer
    float                       *excitation; ///< pointer to the current excitation vector in excitation_buf

    float   pitch_vector[AMR_SUBFRAME_SIZE]; ///< adaptive code book (pitch) vector

    float               prediction_error[4]; ///< quantified prediction errors {20log10(^gamma_gc)} for previous four subframes
    float                     pitch_gain[5]; ///< quantified pitch gains for the current and previous four subframes
    float                     fixed_gain[5]; ///< quantified fixed gains for the current and previous four subframes

    float                              beta; ///< previous pitch_gain, bounded by [0.0,SHARP_MAX]
    uint8_t                      diff_count; ///< the number of subframes for which diff has been above 0.65
    uint8_t                      hang_count; ///< the number of subframes since a hangover period started

    float            prev_sparse_fixed_gain; ///< previous fixed gain; used by anti-sparseness processing to determine "onset"
    uint8_t         prev_ir_filter_strength; ///< previous impulse response filter strength; 0 - strong, 1 - medium, 2 - none
    uint8_t                 ir_filter_onset; ///< flag for impulse response filter strength

    float                postfilter_mem[10]; ///< previous intermediate values in the formant filter
    float                          tilt_mem; ///< previous input to tilt compensation filter
    float                    postfilter_agc; ///< previous factor used for adaptive gain control
    float                  high_pass_mem[2]; ///< previous intermediate values in the high-pass filter

    float samples_in[LP_FILTER_ORDER + AMR_SUBFRAME_SIZE]; ///< floating point samples

} AMRContext;

static void reset_state(AMRContext *p)
{
    int i;

    for (i = 0; i < LP_FILTER_ORDER; i++) {
        p->prev_lsp_sub4[i] = lsp_sub4_init[i] * 1000 / (float)(1 << 15);
        p->lsf_avg[i]       =
        p->lsf_q[3][i]      = lsp_avg_init[i]         / (float)(1 << 15);
    }

    for (i = 0; i < 4; i++)
        p->prediction_error[i] = MIN_ENERGY;
}

static av_cold int amrnb_decode_init(AVCodecContext *avctx)
{
    AMRContext *p = avctx->priv_data;

    avctx->sample_fmt = SAMPLE_FMT_FLT;

    // p->excitation always points to the same position in p->excitation_buf
    p->excitation = &p->excitation_buf[PITCH_LAG_MAX + LP_FILTER_ORDER + 1];

    reset_state(p);

    return 0;
}


/**
 * Unpack an RFC4867 speech frame into the AMR frame mode and parameters.
 *
 * The order of speech bits is specified by 3GPP TS 26.101.
 *
 * @param p the context
 * @param buf               pointer to the input buffer
 * @param buf_size          size of the input buffer
 *
 * @return the frame mode
 */
static enum Mode unpack_bitstream(AMRContext *p, const uint8_t *buf,
                                  int buf_size)
{
    GetBitContext gb;
    enum Mode mode;

    init_get_bits(&gb, buf, buf_size * 8);

    // Decode the first octet.
    skip_bits(&gb, 1);                        // padding bit
    mode = get_bits(&gb, 4);                  // frame type
    p->bad_frame_indicator = !get_bits1(&gb); // quality bit
    skip_bits(&gb, 2);                        // two padding bits

    if (mode <= MODE_DTX) {
        uint16_t *data = (uint16_t *)&p->frame;
        const AMROrder *order = amr_unpacking_bitmaps_per_mode[mode];
        int i;

        memset(&p->frame, 0, sizeof(AMRNBFrame));
        for (i = 0; i < mode_bits[mode]; i++)
            data[order[i].index] += get_bits1(&gb) << order[i].bit;
    }

    return mode;
}


/// @defgroup amr_lpc_decoding AMR pitch LPC coefficient decoding functions
/// @{

/**
 * Convert an lsf vector into an lsp vector.
 *
 * @param lsf               input lsf vector
 * @param lsp               output lsp vector
 */
static void lsf2lsp(float *lsf, float *lsp)
{
    int i;

    for (i = 0; i < LP_FILTER_ORDER; i++)
        lsp[i] = cos(lsf[i] * FREQ_LSP_FAC); // FREQ_LSP_FAC = 2*M_PI / 8000.0
}

/**
 * Interpolate the LSF vector (used for fixed gain smoothing).
 * The interpolation is done over all four subframes even in MODE_122.
 *
 * @param[in,out] lsf_q     LSFs in [0,1] for each subframe
 * @param[in]     lsf_new   New LSFs in Hertz for subframe 4
 */
static void interpolate_lsf(float lsf_q[4][LP_FILTER_ORDER], float *lsf_new)
{
    int i;

    for (i = 0; i < 4; i++)
        ff_weighted_vector_sumf(lsf_q[i], lsf_q[3], lsf_new,
                                0.25 * (3 - i), 0.25 * (i + 1) * FREQ_LSF,
                                LP_FILTER_ORDER);
}

/**
 * Adjust the quantized LSFs so they are increasing and not too close.
 *
 * This step isn't mentioned in the spec but is in the reference C decoder.
 * Omitting this step creates audible distortion on the sinusoidal sweep
 * test vectors in 3GPP TS 26.074.
 *
 * @param[in,out] lsf    LSFs in Hertz
 */
static void adjust_lsf(float *lsf)
{
    int i;
    float prev = 0.0;
    for (i = 0; i < LP_FILTER_ORDER; i++)
        prev = lsf[i] = FFMAX(lsf[i], prev + MIN_LSF_SPACING);
}

/**
 * Decode a set of 5 split-matrix quantized lsf indexes into an lsp vector.
 *
 * @param p the context
 * @param lsp output LSP vector
 * @param lsf_no_r LSF vector without the residual vector added
 * @param lsf_quantizer pointers to LSF dictionary tables
 * @param quantizer_offset offset in tables
 * @param sign for the 3 dictionary table
 * @param update store data for computing the next frame's LSFs
 */
static void lsf2lsp_for_mode122(AMRContext *p, float lsp[LP_FILTER_ORDER],
                                const float lsf_no_r[LP_FILTER_ORDER],
                                const float *lsf_quantizer[5],
                                const int quantizer_offset,
                                const int sign, const int update)
{
    float lsf[LP_FILTER_ORDER]; // used for both the residual and total LSFs
    int i;

    for (i = 0; i < LP_FILTER_ORDER >> 1; i++)
        memcpy(&lsf[i << 1], &lsf_quantizer[i][quantizer_offset],
               2 * sizeof(float));

    if (sign) {
        lsf[4] *= -1;
        lsf[5] *= -1;
    }

    if (update)
        memcpy(p->prev_lsf_r, lsf, LP_FILTER_ORDER * sizeof(float));

    for (i = 0; i < LP_FILTER_ORDER; i++)
        lsf[i] += lsf_no_r[i];

    adjust_lsf(lsf);

    if (update)
        interpolate_lsf(p->lsf_q, lsf);

    lsf2lsp(lsf, lsp);
}

/**
 * Decode a set of 5 split-matrix quantized lsf indexes into 2 lsp vectors.
 *
 * @param p                 pointer to the AMRContext
 */
static void lsf2lsp_5(AMRContext *p)
{
    const uint16_t *lsf_param = p->frame.lsf;
    float lsf_no_r[LP_FILTER_ORDER]; // LSFs without the residual vector
    const float *lsf_quantizer[5];
    int i;

    lsf_quantizer[0] = lsf_5_1[lsf_param[0]];
    lsf_quantizer[1] = lsf_5_2[lsf_param[1]];
    lsf_quantizer[2] = lsf_5_3[lsf_param[2] >> 1];
    lsf_quantizer[3] = lsf_5_4[lsf_param[3]];
    lsf_quantizer[4] = lsf_5_5[lsf_param[4]];

    for (i = 0; i < LP_FILTER_ORDER; i++)
        lsf_no_r[i] = p->prev_lsf_r[i] * PRED_FAC_MODE_122 + lsf_5_mean[i];

    lsf2lsp_for_mode122(p, p->lsp[1], lsf_no_r, lsf_quantizer, 0, lsf_param[2] & 1, 0);
    lsf2lsp_for_mode122(p, p->lsp[3], lsf_no_r, lsf_quantizer, 2, lsf_param[2] & 1, 1);

    // interpolate LSP vectors at subframes 1 and 3
    ff_weighted_vector_sumf(p->lsp[0], p->prev_lsp_sub4, p->lsp[1], 0.5, 0.5, LP_FILTER_ORDER);
    ff_weighted_vector_sumf(p->lsp[2], p->lsp[1]       , p->lsp[3], 0.5, 0.5, LP_FILTER_ORDER);
}

/**
 * Decode a set of 3 split-matrix quantized lsf indexes into an lsp vector.
 *
 * @param p                 pointer to the AMRContext
 */
static void lsf2lsp_3(AMRContext *p)
{
    const uint16_t *lsf_param = p->frame.lsf;
    float lsf_r[LP_FILTER_ORDER]; // residual LSF vector
    float lsf_q[LP_FILTER_ORDER]; // quantified LSF vector
    const float *lsf_quantizer;
    int i;

    lsf_quantizer = (p->cur_frame_mode == MODE_795 ? lsf_3_1_MODE_795 : lsf_3_1)[lsf_param[0]];
    memcpy(lsf_r, lsf_quantizer, 3 * sizeof(*lsf_r));

    lsf_quantizer = lsf_3_2[lsf_param[1] << (p->cur_frame_mode <= MODE_515)];
    memcpy(lsf_r + 3, lsf_quantizer, 3 * sizeof(*lsf_r));

    lsf_quantizer = (p->cur_frame_mode <= MODE_515 ? lsf_3_3_MODE_515 : lsf_3_3)[lsf_param[2]];
    memcpy(lsf_r + 6, lsf_quantizer, 4 * sizeof(*lsf_r));

    // calculate mean-removed LSF vector and add mean
    for (i = 0; i < LP_FILTER_ORDER; i++)
        lsf_q[i] = lsf_r[i] + p->prev_lsf_r[i] * pred_fac[i] + lsf_3_mean[i];

    adjust_lsf(lsf_q);

    // store data for computing the next frame's LSFs
    interpolate_lsf(p->lsf_q, lsf_q);
    memcpy(p->prev_lsf_r, lsf_r, LP_FILTER_ORDER * sizeof(*lsf_r));

    lsf2lsp(lsf_q, p->lsp[3]);

    // interpolate LSP vectors at subframes 1, 2 and 3
    for (i = 0; i < 3; i++)
        ff_weighted_vector_sumf(p->lsp[i], p->prev_lsp_sub4, p->lsp[3],
                                0.25 * (3 - i), 0.25 * (i + 1),
                                LP_FILTER_ORDER);
}


/**
 * Convert an lsp vector to lpc coefficients.
 *
 * @param lsp                 input lsp vector
 * @param lpc_coeffs          output lpc coefficients
 */
static void lsp2lpc(float *lsp, float *lpc_coeffs)
{
    double lsp_double[LP_FILTER_ORDER];
    int i;

    for (i = 0; i < LP_FILTER_ORDER; i++)
        lsp_double[i] = lsp[i];

    ff_acelp_lspd2lpc(lsp_double, lpc_coeffs);
}

/// @}


/// @defgroup amr_pitch_vector_decoding AMR pitch vector decoding functions
/// @{

/**
 * Decode the adaptive codebook index to the integer and fractional parts
 * of the pitch lag for one subframe at 1/6 resolution for MODE_122,
 * 1/3 for other modes.
 *
 * The choice of pitch lag is described in 3GPP TS 26.090 section 5.6.1.
 *
 * @param lag_int             integer part of pitch lag of the current subframe
 * @param lag_frac            fractional part of pitch lag of the current subframe
 * @param pitch_index         parsed adaptive codebook (pitch) index
 * @param prev_lag_int        integer part of pitch lag for the previous subframe
 * @param subframe            current subframe number
 * @param mode                mode of the current frame
 */
static void decode_pitch_lag(int *lag_int, int *lag_frac, int pitch_index,
                             const int prev_lag_int, const int subframe,
                             const enum Mode mode)
{
    /* Note n * 10923 >> 15 is floor(x/3) for 0 <= n <= 32767 */
    if (subframe == 0 ||
        (subframe == 2 && mode != MODE_475 && mode != MODE_515)) {
        if (mode == MODE_122) {
            if (pitch_index < 463) {
                *lag_int  = (pitch_index + 5) / 6 + 17;
                *lag_frac = pitch_index - *lag_int * 6 + 105;
            } else {
                *lag_int  = pitch_index - 368;
                *lag_frac = 0;
            }
        } else if (pitch_index < 197) {
            *lag_int  = ((pitch_index + 2) * 10923 >> 15) + 19;
            *lag_frac = pitch_index - *lag_int * 3 + 58;
        } else {
            *lag_int  = pitch_index - 112;
            *lag_frac = 0;
        }
    } else {
        if (mode == MODE_122) {
            *lag_int  = (pitch_index + 5) / 6 - 1;
            *lag_frac = pitch_index - *lag_int * 6 - 3;
            *lag_int += av_clip(prev_lag_int - 5, PITCH_LAG_MIN_MODE_122,
                                PITCH_LAG_MAX - 9);
        } else if (mode <= MODE_67) {
            int search_range_min = av_clip(prev_lag_int - 5, PITCH_LAG_MIN,
                                           PITCH_LAG_MAX - 9);

            // decoding with 4-bit resolution
            if (pitch_index < 4) {
                // integer only precision for [search_range_min, search_range_min+3]
                *lag_int  = pitch_index + search_range_min;
                *lag_frac = 0;
            } else if (pitch_index < 12) {
                // 1/3 fractional precision for [search_range_min+3 1/3, search_range_min+5 2/3]
                *lag_int  = (pitch_index + 1) * 10923 >> 15;
                *lag_frac = pitch_index - *lag_int * 3;
                *lag_int += search_range_min + 2;
            } else {
                // integer only precision for [search_range_min+6, search_range_min+9]
                *lag_int  = pitch_index + search_range_min - 6;
                *lag_frac = 0;
            }
        } else {
            // decoding with 5 or 6 bit resolution, 1/3 fractional precision
            *lag_int  = ((pitch_index + 2) * 10923 >> 15) - 1;
            *lag_frac = pitch_index - *lag_int * 3 - 2;
            if (mode == MODE_795) {
                *lag_int += av_clip(prev_lag_int - 10, PITCH_LAG_MIN,
                                    PITCH_LAG_MAX - 19);
            } else
                *lag_int += av_clip(prev_lag_int - 5, PITCH_LAG_MIN,
                                    PITCH_LAG_MAX - 9);
        }
    }
}

/**
 * Calculate the pitch vector by interpolating the past excitation at the pitch
 * lag using a b60 hamming windowed sinc function.
 *
 * @param pitch_vector buffer that must hold for the previous state of the filter in
 *                     pitch_vector[-PITCH_LAG_MAX-LP_FILTER_ORDER-1, -1]
 * @param lag_int             integer part of pitch lag
 * @param lag_frac            fractional part of pitch lag
 * @param mode                current frame mode
 */
static void interp_pitch_vector(float *pitch_vector, int lag_int,
                                int lag_frac, enum Mode mode)
{
    int n, i;
    const float *b60_idx1, *b60_idx2;
    float *exc_idx;

    lag_frac *= -1;
    if (mode != MODE_122) {
        lag_frac <<= 1;
    }

    if (lag_frac < 0) {
        lag_frac += 6;
        lag_int++;
    }

    b60_idx1 = &b60[    lag_frac];
    b60_idx2 = &b60[6 - lag_frac];
    exc_idx  = &pitch_vector[-lag_int];

    for (n = 0; n < AMR_SUBFRAME_SIZE; n++) {
        pitch_vector[n] = 0.0;
        for (i = 0; i < 10; i++)
            pitch_vector[n] += b60_idx1[6 * i] * exc_idx[-i];
        exc_idx++;
        for (i = 0; i < 10; i++)
            pitch_vector[n] += b60_idx2[6 * i] * exc_idx[ i];
    }
}

static void decode_pitch_vector(AMRContext *p,
                                const AMRNBSubframe *amr_subframe,
                                const int subframe)
{
    int pitch_lag_int, pitch_lag_frac;

    decode_pitch_lag(&pitch_lag_int, &pitch_lag_frac, amr_subframe->p_lag,
                     p->pitch_lag_int, subframe, p->cur_frame_mode);

    interp_pitch_vector(p->excitation, pitch_lag_int, pitch_lag_frac,
                        p->cur_frame_mode);

    p->pitch_lag_int = pitch_lag_int; // store previous lag in a uint8_t
    memcpy(p->pitch_vector, p->excitation, AMR_SUBFRAME_SIZE * sizeof(float));
}

/// @}


/// @defgroup amr_algebraic_code_book AMR algebraic code book (fixed) vector decoding functions
/// @{

/**
 * Reconstruct the algebraic codebook vector.
 *
 * @param pulse_position       vector of pulse positions
 * @param sign                 signs of the pulses
 * @param nr_pulses            number of pulses
 * @param fixed_vector         algebraic codebook vector
 */
static void reconstruct_fixed_vector(int *pulse_position, int sign,
                                     int nr_pulses, float *fixed_vector)
{
    int i;

    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE * sizeof(float));

    for (i = 0; i < nr_pulses; i++)
        fixed_vector[pulse_position[i]] = ((sign >> i) & 1) ? 1.0 : -1.0;
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and
 * construct the algebraic codebook vector for MODE_102.
 *
 * @param fixed_index          positions of the eight pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 */
static void decode_8_pulses_31bits(const int16_t *fixed_index,
                                   float *fixed_vector)
{
    int pulse_position[8];
    int i, temp;

    // coded using 7+3 bits with the 3 LSBs being, individually, the LSB of 1 of
    // the 3 pulses and the upper 7 bits being coded in base 5
    temp = fixed_index[4] >> 3;
    pulse_position[0] = (( temp       % 5) << 1) + ( fixed_index[4]       & 1);
    pulse_position[4] = (((temp /  5) % 5) << 1) + ((fixed_index[4] >> 1) & 1);
    pulse_position[1] = (((temp / 25) % 5) << 1) + ((fixed_index[4] >> 2) & 1);

    // coded using 7+3 bits with the 3 LSBs being, individually, the LSB of 1 of
    // the 3 pulses and the upper 7 bits being coded in base 5
    temp = fixed_index[5] >> 3;
    pulse_position[2] = (( temp       % 5) << 1) + ( fixed_index[5]       & 1);
    pulse_position[6] = (((temp /  5) % 5) << 1) + ((fixed_index[5] >> 1) & 1);
    pulse_position[5] = (((temp / 25) % 5) << 1) + ((fixed_index[5] >> 2) & 1);

    // coded using 5+2 bits with the 2 LSBs being, individually, the LSB of 1 of
    // the 2 pulses and the upper 5 bits being coded in base 5
    temp = ((fixed_index[6] >> 2) * 25 + 12) >> 5;
    pulse_position[3] = temp % 5;
    pulse_position[7] = temp / 5;
    if (pulse_position[7] & 1)
        pulse_position[3] = 4 - pulse_position[3];
    pulse_position[3] = (pulse_position[3] << 1) + ( fixed_index[6]       & 1);
    pulse_position[7] = (pulse_position[7] << 1) + ((fixed_index[6] >> 1) & 1);

    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE * sizeof(float));

    for (i = 0; i < TRACKS_MODE_102; i++) {
        const int pos1   = (pulse_position[i]     << 2) + i;
        const int pos2   = (pulse_position[i + 4] << 2) + i;
        const float sign = fixed_index[i] ? -1.0 : 1.0;
        fixed_vector[pos1]  = sign;
        fixed_vector[pos2] += pos2 < pos1 ? -sign : sign;
    }
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and
 * construct the algebraic codebook vector for MODE_122.
 *
 * @note: The positions and signs are explicitly coded in MODE_122.
 *
 * @param fixed_index          positions of the ten pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 */
static void decode_10_pulses_35bits(const int16_t *fixed_index,
                                    float *fixed_vector)
{
    int i;

    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE * sizeof(float));

    for (i = 0; i < TRACKS; i++) {
        const int pos1   = gray_decode[fixed_index[i    ] & 7] * TRACKS + i;
        const int pos2   = gray_decode[fixed_index[i + 5] & 7] * TRACKS + i;
        const float sign = (fixed_index[i] & 8) ? -1.0 : 1.0;
        fixed_vector[pos1]  = sign;
        fixed_vector[pos2] += pos2 < pos1 ? -sign : sign;
    }
}

/**
 * Decode the algebraic codebook index to pulse positions and signs,
 * then construct the algebraic codebook vector.
 *
 *                           nb of pulses | bits encoding pulses
 * For MODE_475 or MODE_515,            2 | 1-3, 4-6, 7
 *                  MODE_59,            2 | 1,   2-4, 5-6, 7-9
 *                  MODE_67,            3 | 1-3, 4,   5-7, 8,  9-11
 *      MODE_74 or MODE_795,            4 | 1-3, 4-6, 7-9, 10, 11-13
 *
 * @param fixed_vector pointer to the algebraic codebook vector
 * @param pulses       algebraic codebook indexes
 * @param mode         mode of the current frame
 * @param subframe     current subframe number
 */
static void decode_fixed_vector(float *fixed_vector, const uint16_t *pulses,
                                const enum Mode mode, const int subframe)
{
    assert(MODE_475 <= mode && mode <= MODE_122);

    if (mode == MODE_122) {
        decode_10_pulses_35bits(pulses, fixed_vector);
    } else if (mode == MODE_102) {
        decode_8_pulses_31bits(pulses, fixed_vector);
    } else {
        int pulse_position[4], pulse_subset;
        const int fixed_index = pulses[0];

        if (mode <= MODE_515) {
            pulse_subset      = ((fixed_index >> 3) & 8)     + (subframe << 1);
            pulse_position[0] = ( fixed_index       & 7) * 5 + track_position[pulse_subset];
            pulse_position[1] = ((fixed_index >> 3) & 7) * 5 + track_position[pulse_subset + 1];
        } else if (mode == MODE_59) {
            pulse_subset      = ((fixed_index & 1) << 1) + 1;
            pulse_position[0] = ((fixed_index >> 1) & 7) * 5 + pulse_subset;
            pulse_subset      = (fixed_index  >> 4) & 3;
            pulse_position[1] = ((fixed_index >> 6) & 7) * 5 + pulse_subset + (pulse_subset == 3 ? 1 : 0);
        } else if (mode == MODE_67) {
            pulse_position[0] = (fixed_index        & 7) * 5;
            pulse_subset      = (fixed_index  >> 2) & 2;
            pulse_position[1] = ((fixed_index >> 4) & 7) * 5 + pulse_subset + 1;
            pulse_subset      = (fixed_index  >> 6) & 2;
            pulse_position[2] = ((fixed_index >> 8) & 7) * 5 + pulse_subset + 2;
        } else { // mode <= MODE_795
            pulse_position[0] = gray_decode[ fixed_index        & 7] * 5;
            pulse_position[1] = gray_decode[(fixed_index >> 3)  & 7] * 5 + 1;
            pulse_position[2] = gray_decode[(fixed_index >> 6)  & 7] * 5 + 2;
            pulse_subset      = (fixed_index >> 9) & 1;
            pulse_position[3] = gray_decode[(fixed_index >> 10) & 7] * 5 + pulse_subset + 3;
        }
        reconstruct_fixed_vector(pulse_position, pulses[1],
                                 pulses_nb_per_mode[mode], fixed_vector);
    }
}

/**
 * Apply pitch lag to the fixed vector (section 6.1.2)
 *
 * @param p the context
 * @param subframe unpacked amr subframe
 * @param mode mode of the current frame
 * @param fixed_vector vector to be modified
 */
static void pitch_sharpening(AMRContext *p, int subframe, enum Mode mode,
                             float *fixed_vector)
{
    int i;

    // The spec suggests the current pitch gain is always used, but in other
    // modes the pitch and codebook gains are joinly quantized (sec 5.8.2)
    // so the codebook gain cannot depend on the quantized pitch gain.
    if (mode == MODE_122)
        p->beta = FFMIN(p->pitch_gain[4], 1.0);

    // conduct pitch sharpening as appropriate (section 6.1.2)
    if (p->pitch_lag_int < AMR_SUBFRAME_SIZE)
        for (i = p->pitch_lag_int; i < AMR_SUBFRAME_SIZE; i++)
            fixed_vector[i] += p->beta * fixed_vector[i - p->pitch_lag_int];

    // Save pitch sharpening factor for the next subframe
    // MODE_475 only updates on the 2nd and 4th subframes - this follows from
    // the fact that the gains for two subframes are jointly quantized.
    if (mode != MODE_475 || subframe & 1)
        p->beta = av_clipf(p->pitch_gain[4], 0.0, SHARP_MAX);
}

/// @}


/// @defgroup amr_gain_decoding AMR gain decoding functions
/// @{

/**
 * fixed gain smoothing
 * Note that where the spec specifies the "spectrum in the q domain"
 * in section 6.1.4, in fact frequencies should be used.
 *
 * @param p the context
 * @param lsf LSFs for the current subframe, in the range [0,1]
 * @param lsf_avg averaged LSFs
 * @param mode mode of the current frame
 *
 * @return fixed gain smoothed
 */
static float fixed_gain_smooth(AMRContext *p , const float *lsf,
                               const float *lsf_avg, const enum Mode mode)
{
    float diff = 0.0;
    int i;

    for (i = 0; i < LP_FILTER_ORDER; i++)
        diff += fabs(lsf_avg[i] - lsf[i]) / lsf_avg[i];

    // If diff is large for ten subframes, disable smoothing for a 40-subframe
    // hangover period.
    p->diff_count = diff > 0.65 ? p->diff_count + 1 : 0;

    if (p->diff_count > 10)
        p->hang_count = 0;

    if (p->hang_count < 40) {
        p->hang_count++;
    } else if (mode < MODE_74 || mode == MODE_102) {
        const float smoothing_factor = av_clipf(4.0 * diff - 1.6, 0.0, 1.0);
        const float fixed_gain_mean = (p->fixed_gain[0] + p->fixed_gain[1] +
                                       p->fixed_gain[2] + p->fixed_gain[3] +
                                       p->fixed_gain[4]) * 0.2;
        return smoothing_factor * p->fixed_gain[4] +
               (1.0 - smoothing_factor) * fixed_gain_mean;
    }
    return p->fixed_gain[4];
}

/**
 * Decode pitch gain and fixed gain factor (part of section 6.1.3).
 *
 * @param p the context
 * @param amr_subframe unpacked amr subframe
 * @param mode mode of the current frame
 * @param subframe current subframe number
 * @param fixed_gain_factor decoded gain correction factor
 */
static void decode_gains(AMRContext *p, const AMRNBSubframe *amr_subframe,
                         const enum Mode mode, const int subframe,
                         float *fixed_gain_factor)
{
    if (mode == MODE_122 || mode == MODE_795) {
        p->pitch_gain[4]   = qua_gain_pit [amr_subframe->p_gain];
        *fixed_gain_factor = qua_gain_code[amr_subframe->fixed_gain];
    } else {
        const float *gains =
            mode >= MODE_67  ? gains_high[amr_subframe->p_gain] :
            mode >= MODE_515 ? gains_low [amr_subframe->p_gain] :
                // gain index is only coded in subframes 0,2 for MODE_475
                gains_MODE_475[(p->frame.subframe[subframe & 2].p_gain << 1) +
                               (subframe & 1)];

        p->pitch_gain[4]   = gains[0];
        *fixed_gain_factor = gains[1];
    }
}

/**
 * Calculate fixed gain (part of section 6.1.3)
 *
 * @param p the context
 * @param mode mode of the current frame
 * @param fixed_gain_factor gain correction factor
 * @param fixed_energy decoded algebraic codebook vector energy
 */
static void set_fixed_gain(AMRContext *p, const enum Mode mode,
                           float fixed_gain_factor, float fixed_energy)
{
    // ^g_c = ^gamma_gc * g_c' (equation 69)
    p->fixed_gain[4] = fixed_gain_factor *
        // Eqn 67: gc' = 10^0.05 (predicted dB + mean dB - dB of fixed vector)
        exp2f(log2f(10.0) * 0.05
                        * (ff_dot_productf(energy_pred_fac,
                                           p->prediction_error,
                                           4) + // predicted fixed energy
                           energy_mean[mode])) /
        // 10^(0.05 * -10log(average x^2)) = 1/sqrt((average x^2))
        sqrtf(fixed_energy / AMR_SUBFRAME_SIZE);

    // update quantified prediction error energy history
    memmove(&p->prediction_error[0], &p->prediction_error[1],
            3 * sizeof(p->prediction_error[0]));
    p->prediction_error[3] = 20.0 * log10f(fixed_gain_factor);
}

/// @}


/// @defgroup amr_pre_processing AMR pre-processing functions
/// @{

/**
 * Reduce fixed vector sparseness by smoothing with one of three IR filters.
 *
 * This implements 3GPP TS 26.090 section 6.1(5).
 *
 * In the patent description "Method and device for coding speech in
 * analysis-by-synthesis speech coders" by Ari P. Heikkinen, this method
 * is called "adaptive phase dispersion". This name is also used in the
 * reference source, but not in the spec.
 *
 * @param p the context
 * @param fixed_vector algebraic codebook vector
 * @param fixed_gain smoothed gain
 * @param spare_vector space for modified vector if necessary
 */
static float *anti_sparseness(AMRContext *p, float *fixed_vector,
                              float fixed_gain, float *spare_vector)
{
    int ir_filter_strength;

    if (p->pitch_gain[4] < 0.6) {
        ir_filter_strength = 0;      // strong filtering
    } else if (p->pitch_gain[4] < 0.9) {
        ir_filter_strength = 1;      // medium filtering
    } else
        ir_filter_strength = 2;      // no filtering

    // detect 'onset'
    if (fixed_gain > 2.0 * p->prev_sparse_fixed_gain) {
        p->ir_filter_onset = 2;
    } else if (p->ir_filter_onset)
        p->ir_filter_onset--;

    if (!p->ir_filter_onset) {
        int i, count = 0;

        for (i = 0; i < 5; i++)
            if (p->pitch_gain[i] < 0.6)
                count++;
        if (count > 2)
            ir_filter_strength = 0;

        if (ir_filter_strength > p->prev_ir_filter_strength + 1)
            ir_filter_strength--;
    } else if (ir_filter_strength < 2)
        ir_filter_strength++;

    // Disable filtering for very low level of fixed_gain.
    // Note this step is not specified in the technical description but is in
    // the reference source in the function Ph_disp.
    if (fixed_gain < 5.0)
        ir_filter_strength = 2;

    if (p->cur_frame_mode != MODE_74 && p->cur_frame_mode < MODE_102
         && ir_filter_strength < 2) {
        const float **filters = p->cur_frame_mode == MODE_795 ?
            ir_filters_lookup_MODE_795 : ir_filters_lookup;

        ff_celp_convolve_circf(spare_vector, fixed_vector,
                               filters[ir_filter_strength], AMR_SUBFRAME_SIZE);
        fixed_vector = spare_vector;
    }

    // update ir filter strength history
    p->prev_ir_filter_strength = ir_filter_strength;
    p->prev_sparse_fixed_gain  = fixed_gain;

    return fixed_vector;
}

/// @}


/// @defgroup amr_synthesis AMR synthesis functions
/// @{

/**
 * Conduct 10th order linear predictive coding synthesis.
 *
 * @param p             pointer to the AMRContext
 * @param lpc           pointer to the LPC coefficients
 * @param fixed_gain    fixed codebook gain for synthesis
 * @param fixed_vector  algebraic codebook vector
 * @param samples       pointer to the output speech samples
 * @param overflow      16-bit overflow flag
 */
static int synthesis(AMRContext *p, float *lpc,
                     float fixed_gain, float *fixed_vector, float *samples,
                     uint8_t overflow)
{
    int i, overflow_temp = 0;
    float excitation[AMR_SUBFRAME_SIZE];

    // if an overflow has been detected, the pitch vector is scaled down by a
    // factor of 4
    if (overflow)
        for (i = 0; i < AMR_SUBFRAME_SIZE; i++)
            p->pitch_vector[i] *= 0.25;

    ff_weighted_vector_sumf(excitation, p->pitch_vector, fixed_vector,
                            p->pitch_gain[4], fixed_gain, AMR_SUBFRAME_SIZE);

    // emphasize pitch vector contribution
    if (p->pitch_gain[4] > 0.5 && !overflow) {
        float energy = ff_energyf(excitation, AMR_SUBFRAME_SIZE);
        float pitch_factor = (p->cur_frame_mode == MODE_122 ? 0.25 : 0.5)
            * FFMIN(p->pitch_gain[4],
                    p->cur_frame_mode == MODE_122 ? 1.0 : SHARP_MAX)
            * p->pitch_gain[4];

        for (i = 0; i < AMR_SUBFRAME_SIZE; i++)
            excitation[i] += pitch_factor * p->pitch_vector[i];

        ff_set_energyf(excitation, excitation, energy, AMR_SUBFRAME_SIZE);
    }

    ff_celp_lp_synthesis_filterf(samples, lpc, excitation, AMR_SUBFRAME_SIZE,
                                 LP_FILTER_ORDER);

    // detect overflow
    for (i = 0; i < AMR_SUBFRAME_SIZE; i++)
        if (fabsf(samples[i]) > AMR_SAMPLE_BOUND) {
            overflow_temp = 1;
            samples[i] = av_clipf(samples[i], -AMR_SAMPLE_BOUND,
                                               AMR_SAMPLE_BOUND);
        }

    return overflow_temp;
}

/// @}


/// @defgroup amr_update AMR update functions
/// @{

/**
 * Update buffers and history at the end of decoding a subframe.
 *
 * @param p             pointer to the AMRContext
 */
static void update_state(AMRContext *p)
{
    memcpy(p->prev_lsp_sub4, p->lsp[3], LP_FILTER_ORDER * sizeof(float));

    memmove(&p->excitation_buf[0], &p->excitation_buf[AMR_SUBFRAME_SIZE],
            (PITCH_LAG_MAX + LP_FILTER_ORDER + 1) * sizeof(float));

    memmove(&p->pitch_gain[0], &p->pitch_gain[1], 4 * sizeof(float));
    memmove(&p->fixed_gain[0], &p->fixed_gain[1], 4 * sizeof(float));

    memmove(&p->samples_in[0], &p->samples_in[AMR_SUBFRAME_SIZE],
            LP_FILTER_ORDER * sizeof(float));
}

/// @}


/// @defgroup amr_postproc AMR Post processing functions
/// @{

/**
 * Get the tilt factor of a formant filter from its transfer function
 *
 * @param lpc_n LP_FILTER_ORDER coefficients of the numerator
 * @param lpc_d LP_FILTER_ORDER coefficients of the denominator
 */
static float tilt_factor(float *lpc_n, float *lpc_d)
{
    float rh0, rh1; // autocorrelation at lag 0 and 1

    // LP_FILTER_ORDER prior zeros are needed for ff_celp_lp_synthesis_filterf
    float impulse_buffer[LP_FILTER_ORDER + AMR_TILT_RESPONSE] = { 0 };
    float *hf = impulse_buffer + LP_FILTER_ORDER; // start of impulse response

    hf[0] = 1.0;
    memcpy(hf + 1, lpc_n, sizeof(float) * LP_FILTER_ORDER);
    ff_celp_lp_synthesis_filterf(hf, lpc_d, hf, AMR_TILT_RESPONSE,
                                 LP_FILTER_ORDER);

    rh0 = ff_dot_productf(hf, hf,     AMR_TILT_RESPONSE);
    rh1 = ff_dot_productf(hf, hf + 1, AMR_TILT_RESPONSE - 1);

    // The spec only specifies this check for 12.2 and 10.2 kbit/s
    // modes. But in the ref source the tilt is always non-negative.
    return rh1 >= 0.0 ? rh1 / rh0 * AMR_TILT_GAMMA_T : 0.0;
}

/**
 * Apply tilt compensation filter, 1 - tilt * z^-1
 *
 * @param mem Pointer to one float to keep the filter's state
 * @param tilt Tilt factor
 * @param samples AMR_SUBFRAME_SIZE array where the filter is applied
 */
static void tilt_compensation(float *mem, float tilt, float *samples)
{
    float new_tilt_mem = samples[AMR_SUBFRAME_SIZE - 1];
    int i;

    for (i = AMR_SUBFRAME_SIZE - 1; i > 0; i--)
         samples[i] -= tilt * samples[i - 1];

    samples[0] -= tilt * *mem;
    *mem = new_tilt_mem;
}

/**
 * Perform adaptive post-filtering to enhance the quality of the speech.
 * See section 6.2.1.
 *
 * @param p             pointer to the AMRContext
 * @param lpc           interpolated LP coefficients for this subframe
 * @param buf_out       output of the filter
 */
static void postfilter(AMRContext *p, float *lpc, float *buf_out)
{
    int i;
    float *samples          = p->samples_in + LP_FILTER_ORDER; // Start of input

    float gain_scale_factor = 1.0;
    float speech_gain       = ff_energyf(samples, AMR_SUBFRAME_SIZE);
    float postfilter_gain;

    float pole_out[AMR_SUBFRAME_SIZE + LP_FILTER_ORDER];  // Output of pole filter
    const float *gamma_n, *gamma_d;                       // Formant filter factor table
    float lpc_n[LP_FILTER_ORDER], lpc_d[LP_FILTER_ORDER]; // Transfer function coefficients

    if (p->cur_frame_mode == MODE_122 || p->cur_frame_mode == MODE_102) {
        gamma_n = formant_high_n;
        gamma_d = formant_high_d;
    } else {
        gamma_n = formant_low_n;
        gamma_d = formant_low_d;
    }

    for (i = 0; i < LP_FILTER_ORDER; i++) {
         lpc_n[i] = lpc[i] * gamma_n[i];
         lpc_d[i] = lpc[i] * gamma_d[i];
    }

    memcpy(pole_out, p->postfilter_mem, sizeof(float) * LP_FILTER_ORDER);
    ff_celp_lp_synthesis_filterf(pole_out + LP_FILTER_ORDER, lpc_d, samples,
                                 AMR_SUBFRAME_SIZE, LP_FILTER_ORDER);
    memcpy(p->postfilter_mem, pole_out + AMR_SUBFRAME_SIZE,
           sizeof(float) * LP_FILTER_ORDER);

    ff_celp_lp_zero_synthesis_filterf(buf_out, lpc_n,
                                      pole_out + LP_FILTER_ORDER,
                                      AMR_SUBFRAME_SIZE, LP_FILTER_ORDER);

    tilt_compensation(&p->tilt_mem, tilt_factor(lpc_n, lpc_d), buf_out);

    // Adaptive gain control
    postfilter_gain = ff_energyf(buf_out, AMR_SUBFRAME_SIZE);
    if (postfilter_gain)
        gain_scale_factor = sqrt(speech_gain / postfilter_gain);

    for (i = 0; i < AMR_SUBFRAME_SIZE; i++) {
        p->postfilter_agc = AMR_AGC_ALPHA * p->postfilter_agc +
                            (1.0 - AMR_AGC_ALPHA) * gain_scale_factor;
        buf_out[i] *= p->postfilter_agc;
    }
}

/// @}

static int amrnb_decode_frame(AVCodecContext *avctx, void *data, int *data_size,
                              AVPacket *avpkt)
{

    AMRContext *p = avctx->priv_data;        // pointer to private data
    const uint8_t *buf = avpkt->data;
    int buf_size       = avpkt->size;
    float *buf_out = data;                   // pointer to the output data buffer
    int i, subframe;
    float fixed_gain_factor;
    float fixed_vector[AMR_SUBFRAME_SIZE];   // algebraic code book (fixed) vector
    float spare_vector[AMR_SUBFRAME_SIZE];   // extra stack space to hold result from anti-sparseness processing
    float synth_fixed_gain;                  // the fixed gain that synthesis should use
    float *synth_fixed_vector;               // pointer to the fixed vector that synthesis should use

    p->cur_frame_mode = unpack_bitstream(p, buf, buf_size);
    if (p->cur_frame_mode == MODE_DTX) {
        av_log_missing_feature(avctx, "dtx mode", 1);
        return -1;
    }

    if (p->cur_frame_mode == MODE_122) {
        lsf2lsp_5(p);
    } else
        lsf2lsp_3(p);

    for (i = 0; i < 4; i++)
        lsp2lpc(p->lsp[i], p->lpc[i]);

    for (subframe = 0; subframe < 4; subframe++) {
        const AMRNBSubframe *amr_subframe = &p->frame.subframe[subframe];

        decode_pitch_vector(p, amr_subframe, subframe);

        decode_fixed_vector(fixed_vector, amr_subframe->pulses,
                            p->cur_frame_mode, subframe);

        // The fixed gain (section 6.1.3) depends on the fixed vector
        // (section 6.1.2), but the fixed vector calculation uses
        // pitch sharpening based on the on the pitch gain (section 6.1.3).
        // So the correct order is: pitch gain, pitch sharpening, fixed gain.
        decode_gains(p, amr_subframe, p->cur_frame_mode, subframe,
                     &fixed_gain_factor);

        pitch_sharpening(p, subframe, p->cur_frame_mode, fixed_vector);

        set_fixed_gain(p, p->cur_frame_mode, fixed_gain_factor,
                       ff_energyf(fixed_vector, AMR_SUBFRAME_SIZE));

        // The excitation feedback is calculated without any processing such
        // as fixed gain smoothing. This isn't mentioned in the specification.
        ff_weighted_vector_sumf(p->excitation, p->excitation, fixed_vector,
                                p->pitch_gain[4], p->fixed_gain[4],
                                AMR_SUBFRAME_SIZE);

        // In the ref decoder, excitation is stored with no fractional bits.
        // This step prevents buzz in silent periods. The ref encoder can
        // emit long sequences with pitch factor greater than one. This
        // creates unwanted feedback if the excitation vector is nonzero.
        // (e.g. test sequence T19_795.COD in 3GPP TS 26.074)
        for (i = 0; i < AMR_SUBFRAME_SIZE; i++)
            p->excitation[i] = truncf(p->excitation[i]);

        // Smooth fixed gain.
        // The specification is ambiguous, but in the reference source, the
        // smoothed value is NOT fed back into later fixed gain smoothing.
        synth_fixed_gain = fixed_gain_smooth(p, p->lsf_q[subframe],
                                              p->lsf_avg, p->cur_frame_mode);

        synth_fixed_vector = anti_sparseness(p, fixed_vector, synth_fixed_gain,
                                             spare_vector);

        if (synthesis(p, p->lpc[subframe], synth_fixed_gain,
                      synth_fixed_vector, &p->samples_in[LP_FILTER_ORDER], 0))
            // overflow detected -> rerun synthesis scaling pitch vector down
            // by a factor of 4, skipping pitch vector contribution emphasis
            // and adaptive gain control
            synthesis(p, p->lpc[subframe], synth_fixed_gain,
                      synth_fixed_vector, &p->samples_in[LP_FILTER_ORDER], 1);

        postfilter(p, p->lpc[subframe], buf_out + subframe * AMR_SUBFRAME_SIZE);

        // update buffers and history
        update_state(p);
    }

    ff_acelp_high_pass_filterf(buf_out, p->high_pass_mem, AMR_BLOCK_SIZE);

    for (i = 0; i < AMR_BLOCK_SIZE; i++)
        buf_out[i] = av_clipf(buf_out[i] * AMR_SAMPLE_SCALE,
                              -1.0, 32767.0/32768.0);

    /* Update averaged lsf vector (used for fixed gain smoothing).
     *
     * Note that lsf_avg should not incorporate the current frame's LSFs
     * for fixed_gain_smooth.
     * The specification has an incorrect formula: the reference decoder uses
     * qbar(n-1) rather than qbar(n) in section 6.1(4) equation 71. */
    ff_weighted_vector_sumf(p->lsf_avg, p->lsf_avg, p->lsf_q[3],
                            0.84, 0.16, LP_FILTER_ORDER);

    /* report how many samples we got */
    *data_size = AMR_BLOCK_SIZE * sizeof(float);

    /* return the amount of bytes consumed if everything was OK */
    return (mode_bits[p->cur_frame_mode] + 15) >> 3; // +7 for rounding and +8 for TOC
}


AVCodec amrnb_decoder = {
    .name           = "amrnb",
    .type           = CODEC_TYPE_AUDIO,
    .id             = CODEC_ID_AMR_NB,
    .priv_data_size = sizeof(AMRContext),
    .init           = amrnb_decode_init,
    .decode         = amrnb_decode_frame,
    .long_name      = NULL_IF_CONFIG_SMALL("Adaptive Multi-Rate NarrowBand"),
    .sample_fmts    = (enum SampleFormat[]){SAMPLE_FMT_FLT,SAMPLE_FMT_NONE},
};
