/*
 * AMR narrowband decoder (floating point)
 * Copyright (c) 2006-2007 Robert Swain
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
 * @file amrnbfloatdec.c
 * AMR narrowband decoder (floating point)
 */


#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "avcodec.h"
#include "bitstream.h"
#include "libavutil/common.h"
#include "libavcodec/internal.h"
#include "amrnbfloatdata.h"

void ff_celp_lspf2lpc(const double *lspf, float *lpc);

typedef struct AMRContext {

    GetBitContext                        gb;

    AMRNBFrame                        frame; ///< decoded AMR parameters (lsf coefficients, codebook indexes, etc)
    int                 bad_frame_indicator; ///< bad frame ? 1 : 0
    int                      cur_frame_mode; ///< current frame mode
    int                      cur_frame_type; ///< current frame type

    float       prev_lsf_r[LP_FILTER_ORDER]; ///< residual LSF vector from previous subframe
    float           lsp[4][LP_FILTER_ORDER]; ///< lsp vectors from current frame
    float    prev_lsp_sub4[LP_FILTER_ORDER]; ///< lsp vector for the 4th subframe of the previous frame

    float          lsp_avg[LP_FILTER_ORDER]; ///< vector of averaged lsp coefficients

    float           lpc[4][LP_FILTER_ORDER]; ///< lpc coefficient vectors for 4 subframes

    int                    search_range_min; ///< minimum pitch lag search range
    int                    search_range_max; ///< maximum pitch lag search range
    int                       pitch_lag_int; ///< integer part of pitch lag from current subframe
    int                      pitch_lag_frac; ///< fractional part of pitch lag from current subframe
    int                  prev_pitch_lag_int; ///< integer part of pitch lag from previous subframe

    float excitation_buf[PITCH_LAG_MAX + LP_FILTER_ORDER + 1 + AMR_SUBFRAME_SIZE]; ///< excitation buffer
    float                       *excitation; ///< pointer to the current excitation vector in excitation_buf

    float   pitch_vector[AMR_SUBFRAME_SIZE]; ///< adaptive code book (pitch) vector
    float   fixed_vector[AMR_SUBFRAME_SIZE]; ///< algebraic code book (fixed) vector

    float               prediction_error[4]; ///< quantified prediction errors {20log10(^γ_gc)} for previous four subframes
    float                     pitch_gain[5]; ///< quantified pitch gains for the current and previous four subframes
    float                 fixed_gain_factor; ///< fixed gain correction factor {^γ_gc} for the current frame
    float                     fixed_gain[5]; ///< quantified fixed gains for the current and previous four subframes

    float                              beta; ///< beta = pitch_gain, bounded by [0.0,1.0] for 12.2 kbps or [0.0,0.8] for other modes
    int                          diff_count; ///< the number of subframes for which diff has been above 0.65

    uint8_t           ir_filter_strength[2]; ///< impulse response filter strength; 0 - strong, 1 - medium, 2 - none
    const float                  *ir_filter; ///< pointer to impulse response filter data

    float samples_in[LP_FILTER_ORDER + AMR_SUBFRAME_SIZE]; ///< floating point samples

} AMRContext;

static void weighted_vector_sumf(float *out, const float *in_a,
                                 const float *in_b, float weight_coeff_a,
                                 float weight_coeff_b, int length)
{
    int i;

    for(i=0; i<length; i++)
        out[i] = weight_coeff_a * in_a[i]
               + weight_coeff_b * in_b[i];
}

static void reset_state(AMRContext *p)
{
    int i;

    for(i=0; i<LP_FILTER_ORDER; i++) {
      p->prev_lsp_sub4[i] = lsp_sub4_init[i] * 1000 / (float)(1 << 15);
      p->lsp_avg[i]       = lsp_avg_init[i]         / (float)(1 << 15);
    }

    for(i=0; i<4; i++) {
        p->prediction_error[i] = MIN_ENERGY;
    }
}

static av_cold int amrnb_decode_init(AVCodecContext *avctx)
{
    AMRContext *p = avctx->priv_data;

    avctx->sample_fmt = SAMPLE_FMT_FLT;

    // p->excitation always points to the same position in p->excitation_buf
    p->excitation = &p->excitation_buf[PITCH_LAG_MAX + LP_FILTER_ORDER + 1];

    reset_state(p);

    /* return 0 for a successful init, -1 for failure */
    return 0;
}


/**
 * Decode the bitstream into the AMR parameters and discover the frame mode.
 *
 * @param buf               pointer to the input buffer
 * @param buf_size          size of the input buffer
 * @param speech_mode       pointer to the speech mode
 *
 * @return the frame mode
 */

enum Mode decode_bitstream(AMRContext *p, const uint8_t *buf, int buf_size,
                           enum Mode *speech_mode)
{
    enum Mode mode;
    int i;
    const AMROrder *order;

    // initialize get_bits
    init_get_bits(&p->gb, buf, buf_size*8);
    skip_bits(&p->gb, 1);
    // set the mode
    mode = get_bits(&p->gb ,4);
    // set the bad frame indicator based on the quality bit
    p->bad_frame_indicator = !get_bits1(&p->gb);
    skip_bits(&p->gb, 2);

    switch(mode) {
        case MODE_DTX:
            order = order_MODE_DTX;
            p->cur_frame_type = RX_SID_FIRST; // get SID type bit
        break;
        case NO_DATA:
            p->cur_frame_type = RX_NO_DATA;
        break;
        case MODE_475:
            order = order_MODE_475;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_515:
            order = order_MODE_515;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_59:
            order = order_MODE_59;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_67:
            order = order_MODE_67;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_74:
            order = order_MODE_74;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_795:
            order = order_MODE_795;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_102:
            order = order_MODE_102;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_122:
            order = order_MODE_122;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        default:
            p->cur_frame_type = RX_SPEECH_BAD;
        break;
    }

    // reorder the bitstream to match the bit allocation in the specification
    if((p->cur_frame_type != RX_NO_DATA) && (p->cur_frame_type != RX_SPEECH_BAD)) {
        uint16_t *data = (uint16_t *)&p->frame;

        memset(&p->frame, 0, sizeof(AMRNBFrame));
        for(i=0; i<mode_bits[mode]; i++) {
            data[order[i].array_element] += get_bits1(&p->gb) * (1<< order[i].bit_mask);
        }
    }

    if(mode == MODE_DTX) {
        skip_bits(&p->gb, 4); // skip to the next byte
        if(get_bits1(&p->gb)) // use the update if there is one
            p->cur_frame_type = RX_SID_UPDATE;
        *speech_mode = get_bits(&p->gb, 3); // speech mode indicator
    }

    return mode;
}


/// @defgroup amr_lpc_decoding pitch LPC coefficient decoding functions
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

    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsp[i] = cos(lsf[i]*FREQ_LSP_FAC); // FREQ_LSP_FAC = 2*M_PI/8000.0
    }
}

/**
 * Decode a set of 5 split-matrix quantized lsf indexes into an lsp vector.
 *
 * @param p the context
 * @param lsp output LSP vector
 * @param prev_lsf previous LSF vector
 * @param lsf_quantizer pointers to LSF dictionary tables
 * @param quantizer_offset offset in tables
 * @param sign for the 3 dictionary table
 * @param update_prev_lsf_r update the prev_lsf_r in the context if true
 */
static void lsf2lsp_for_mode122(AMRContext *p, float lsp[LP_FILTER_ORDER],
                                const float prev_lsf[LP_FILTER_ORDER],
                                const float *lsf_quantizer[5], const int quantizer_offset,
                                const int sign, const int update_prev_lsf_r)
{
    float lsf[LP_FILTER_ORDER];
    int i;

    for(i=0; i<LP_FILTER_ORDER>>1; i++)
        memcpy(&lsf[2*i], &lsf_quantizer[i][quantizer_offset], 2*sizeof(float));

    if(sign) {
        lsf[4] *= -1;
        lsf[5] *= -1;
    }

    if(update_prev_lsf_r)
        memcpy(p->prev_lsf_r, lsf, LP_FILTER_ORDER*sizeof(float));

    for(i=0; i<LP_FILTER_ORDER; i++)
        lsf[i] += prev_lsf[i];

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
    float prev_lsf[LP_FILTER_ORDER]; // previous quantized LSF vectors
    const float *lsf_quantizer[5];
    int i;

    lsf_quantizer[0] = lsf_5_1[lsf_param[0]];
    lsf_quantizer[1] = lsf_5_2[lsf_param[1]];
    lsf_quantizer[2] = lsf_5_3[lsf_param[2] >> 1];
    lsf_quantizer[3] = lsf_5_4[lsf_param[3]];
    lsf_quantizer[4] = lsf_5_5[lsf_param[4]];

    for(i=0; i<LP_FILTER_ORDER;i++)
        prev_lsf[i] = p->prev_lsf_r[i]*PRED_FAC_MODE_122 + lsf_5_mean[i];

    lsf2lsp_for_mode122(p, p->lsp[1], prev_lsf, lsf_quantizer, 0, lsf_param[2] & 1, 0);
    lsf2lsp_for_mode122(p, p->lsp[3], prev_lsf, lsf_quantizer, 2, lsf_param[2] & 1, 1);

    // interpolate LSP vectors at subframes 1 and 3
    weighted_vector_sumf(p->lsp[0], p->prev_lsp_sub4, p->lsp[1], 0.5, 0.5, LP_FILTER_ORDER);
    weighted_vector_sumf(p->lsp[2], p->lsp[1]       , p->lsp[3], 0.5, 0.5, LP_FILTER_ORDER);
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
    memcpy(lsf_r, lsf_quantizer, 3*sizeof(float));

    lsf_quantizer = lsf_3_2[lsf_param[1] << (p->cur_frame_mode <= MODE_515)];
    memcpy(lsf_r + 3, lsf_quantizer, 3*sizeof(float));

    lsf_quantizer = (p->cur_frame_mode <= MODE_515 ? lsf_3_3_MODE_515 : lsf_3_3)[lsf_param[2]];
    memcpy(lsf_r + 6, lsf_quantizer, 4*sizeof(float));

    // calculate mean-removed LSF vector and add mean
    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsf_q[i] = lsf_r[i] + p->prev_lsf_r[i]*pred_fac[i] + lsf_3_mean[i];
    }
    // update residual LSF vector from previous subframe
    memcpy(p->prev_lsf_r, lsf_r, LP_FILTER_ORDER*sizeof(float));

    // convert LSF vector to LSP vector
    lsf2lsp(lsf_q, p->lsp[3]);

    // interpolate LSP vectors at subframes 1, 2 and 3
    for(i=0; i<3; i++)
        weighted_vector_sumf(p->lsp[i], p->prev_lsp_sub4, p->lsp[3], 0.25*(3-i), 0.25*(i+1), LP_FILTER_ORDER);
}


/**
 * Convert an lsp vector to lpc coefficients.
 *
 * @param lsp                 input lsp vector
 * @param lpc                 output lpc coefficients
 */

static void lsp2lpc(float *lsp, float *lpc_coeffs)
{
    double lsp_double[LP_FILTER_ORDER];
    int i;

    for(i=0; i<LP_FILTER_ORDER; i++)
        lsp_double[i] = lsp[i];

    ff_celp_lspf2lpc(lsp_double, lpc_coeffs);
}

/// @}


/// @defgroup amr_pitch_vector_decoding pitch vector decoding functions
/// @{

/**
 * Decode the adaptive codebook index to the integer and fractional parts
 * of the pitch lag for one subframe at 1/3 resolution.
 *
 * @param p                   pointer to the AMRContext
 * @param pitch_index         parsed adaptive codebook (pitch) index
 * @param subframe            current subframe
 */

static void decode_pitch_lag_3(AMRContext *p, int pitch_index, int subframe)
{
    // subframe 1 or 3
    if(!(subframe & 1)) {
        if(pitch_index < 197) {
            // 10923>>15 is approximately 1/3
            p->pitch_lag_int = ( ((pitch_index + 2)*10923)>>15 ) + 19;
            p->pitch_lag_frac = pitch_index - p->pitch_lag_int*3 + 58;
        }else {
            p->pitch_lag_int = pitch_index - 112;
            p->pitch_lag_frac = 0;
        }
    // subframe 2 or 4
    }else {
        if( (p->cur_frame_mode == MODE_475) || (p->cur_frame_mode == MODE_515) ||
            (p->cur_frame_mode == MODE_59)  || (p->cur_frame_mode == MODE_67) ) {
            // decoding with 4-bit resolution
            int t1_temp = FFMAX(FFMIN(p->prev_pitch_lag_int, p->search_range_max-4), p->search_range_min+5);

            if(pitch_index < 4) {
                // integer only precision for [t1_temp-5, t1_temp-2]
                p->pitch_lag_int = pitch_index + (t1_temp - 5);
                p->pitch_lag_frac = 0;
            }else if(pitch_index < 12) {
                // 1/3 fractional precision for [t1_temp-1 2/3, t1_temp+2/3]
                p->pitch_lag_int = ( ((pitch_index - 5)*10923)>>15 ) - 1;
                p->pitch_lag_frac = pitch_index - p->pitch_lag_int*3 - 9;
                p->pitch_lag_int += t1_temp;
            }else {
                // integer only precision for [t1_temp+1, t1_temp+4]
                p->pitch_lag_int = pitch_index + t1_temp - 11;
                p->pitch_lag_frac = 0;
            }
        }else {
            // decoding with 5 or 6 bit resolution, 1/3 fractional precision
            // 10923>>15 is approximately 1/3
            int temp = ( ((pitch_index + 2)*10923)>>15 ) - 1;
            p->pitch_lag_int = temp + p->search_range_min;
            p->pitch_lag_frac = pitch_index - temp*3 - 2;
        }
    }
}

/**
 * Decode the adaptive codebook index to the integer and fractional parts
 * of the pitch lag for one subframe at 1/6 resolution.
 *
 * @param p                   pointer to the AMRContext
 * @param pitch_index         parsed adaptive codebook (pitch) index
 * @param subframe            current subframe
 */

static void decode_pitch_lag_6(AMRContext *p, int pitch_index, int subframe)
{
    // subframe 1 or 3
    if(!(subframe & 1)) {
        if(pitch_index < 463){
            p->pitch_lag_int = (pitch_index + 5)/6 + 17;
            p->pitch_lag_frac = pitch_index - p->pitch_lag_int*6 + 105;
        }else {
            p->pitch_lag_int = pitch_index - 368;
            p->pitch_lag_frac = 0;
        }
    // subframe 2 or 4
    }else {
        int temp;
        // calculate the pitch lag
        temp = (pitch_index + 5)/6 - 1;
        p->pitch_lag_int = temp + p->search_range_min;
        p->pitch_lag_frac = pitch_index - temp*6 - 3;
    }
}

/**
 * Calculate the pitch vector by interpolating the past excitation at the pitch
 * pitch lag using a b60 hamming windowed sinc function.
 *
 * @param prev_excitation     pointer to the element after the previous excitations
 * @param lag_int             integer part of pitch lag
 * @param lag_frac            fractional part of pitch lag
 * @param mode                current frame mode
 * @param pitch_vector        pointer to the pitch vector
 */

static void interp_pitch_vector(float *prev_excitation, int lag_int,
                                int lag_frac, enum Mode mode,
                                float *pitch_vector)
{
    int n, i;
    const float *b60_idx1, *b60_idx2;
    float *exc_idx;

    lag_frac *= -1;
    if(mode != MODE_122) {
        lag_frac <<= 1;
    }

    if(lag_frac < 0) {
        lag_frac += 6;
        lag_int++;
    }

    b60_idx1 = &b60[    lag_frac];
    b60_idx2 = &b60[6 - lag_frac];
    exc_idx = &prev_excitation[-lag_int];

    for(n=0; n<AMR_SUBFRAME_SIZE; n++) {
        pitch_vector[n] = 0.0;
        for(i=0; i<10; i++) {
            pitch_vector[n] += b60_idx1[6*i] * exc_idx[-i];
        }
        exc_idx++;
        for(i=0; i<10; i++) {
            pitch_vector[n] += b60_idx2[6*i] * exc_idx[ i];
        }
    }
}

/// @}


/// @defgroup amr_algebraic_code_book algebraic code book (fixed) vector decoding functions
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

    // reset the code
    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE*sizeof(float));

    for(i=0; i<nr_pulses; i++)
        fixed_vector[pulse_position[i]] = ((sign >> i) & 1) ? 1.0 : -1.0;
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and
 * construct the algebraic codebook vector for MODE_102.
 *
 * @param fixed_index          positions of the eight pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 */

static void decode_8_pulses_31bits(const int16_t *fixed_index, float *fixed_vector)
{
    int pulse_position[8];
    int i, pos1, pos2, sign, temp;

    // decode pulse positions
    // coded using 7+3 bits with the 3 LSBs being, individually, the LSB of 1 of
    // the 3 pulses and the upper 7 bits being coded in base 5
    temp = fixed_index[4] >> 3;
    pulse_position[0] = (( temp    %5)<<1) + ( fixed_index[4]    &1);
    pulse_position[4] = (((temp /5)%5)<<1) + ((fixed_index[4]>>1)&1);
    pulse_position[1] = (((temp/25)%5)<<1) + ((fixed_index[4]>>2)&1);

    // coded using 7+3 bits with the 3 LSBs being, individually, the LSB of 1 of
    // the 3 pulses and the upper 7 bits being coded in base 5
    temp = fixed_index[5] >> 3;
    pulse_position[2] = (( temp    %5)<<1) + ( fixed_index[5]    &1);
    pulse_position[6] = (((temp /5)%5)<<1) + ((fixed_index[5]>>1)&1);
    pulse_position[5] = (((temp/25)%5)<<1) + ((fixed_index[5]>>2)&1);

    // coded using 5+2 bits with the 2 LSBs being, individually, the LSB of 1 of
    // the 2 pulses and the upper 5 bits being coded in base 5
    temp = ((fixed_index[6] >> 2)*25)>>5;
    pulse_position[3] = temp%5;
    pulse_position[7] = temp/5;
    if(pulse_position[7]&1)
        pulse_position[3] = 4 - pulse_position[3];
    pulse_position[3] = (pulse_position[3]<<1) + ( fixed_index[6]    &1);
    pulse_position[7] = (pulse_position[7]<<1) + ((fixed_index[6]>>1)&1);

    // reset the code
    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE*sizeof(float));

    // reconstruct the fixed code
    for(i=0; i<TRACKS_MODE_102; i++) {
        pos1 = (pulse_position[i]   << 2) + i; // ith pulse position
        pos2 = (pulse_position[i+4] << 2) + i; // i+4th pulse position
        sign = fixed_index[i] ? -1.0 : 1.0; // sign of ith pulse
        // assign the ith pulse (+/-1) to its appropriate position
        fixed_vector[pos1] = sign;
        // sign of i+4th pulse is relative to sign of ith pulse
        if(pos2 < pos1) sign = -sign;
        // assign the i+4th pulse (+/-1) to its appropriate position
        fixed_vector[pos2] += sign;
    }
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and
 * construct the algebraic codebook vector for MODE_122.
 *
 * @param fixed_index          positions of the ten pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 */

static void decode_10_pulses_35bits(const int16_t *fixed_index, float *fixed_vector)
{
    int i, pos1, pos2, sign;

    // reset the code
    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE*sizeof(float));

    // the positions and signs are explicitly coded in MODE_122

    // reconstruct the fixed code
    for(i=0; i<TRACKS; i++) {
        pos1 = gray_decode[fixed_index[i  ] & 7]*TRACKS + i; // ith pulse position
        sign = (fixed_index[i] & 8) ? -1.0 : 1.0; // sign of ith pulse
        pos2 = gray_decode[fixed_index[i+5] & 7]*TRACKS + i; // i+5th pulse position
        // assign the ith pulse (+/-1) to its appropriate position
        fixed_vector[pos1] = sign;
        // sign of i+5th pulse is relative to sign of ith pulse
        if(pos2 < pos1) sign = -sign;
        // assign the i+5th pulse (+/-1) to its appropriate position
        fixed_vector[pos2] += sign;
    }
}

/**
 * Decode the algebraic codebook index to pulse positions and signs,
 * then construct the algebraic codebook vector.
 *
 *                           nb of pulses | sign index | bits encoding pulses
 * For MODE_475 or MODE_515,            2 |      bit 7 | 1-3, 4-7
 *                  MODE_59,            2 |      bit 1 | 2-4, 5-6, 7-9
 *                  MODE_67,            3 |            | 1-3, 4,   5-7, 8,  9-11
 *      MODE_74 or MODE_795,            4 |            | 1-3, 4-6, 7-9, 10, 11-13
 *
 * @param fixed_vector pointer to the algebraic codebook vector
 * @param pulses       algebraic codebook indexes
 * @param mode         mode of the current frame
 * @param subframe     current subframe
 */
static void decode_fixed_vector(float *fixed_vector, const uint16_t *pulses,
                                const enum Mode mode, const int subframe)
{
    assert(MODE_475 <= mode && mode <= MODE_122);

    if (mode == MODE_122) {
        decode_10_pulses_35bits(pulses, fixed_vector);
    }else if(mode == MODE_102) {
        decode_8_pulses_31bits(pulses, fixed_vector);
    }else {
        int pulse_position[4], pulse_subset, pulse_nb;
        const int fixed_index = pulses[0], sign = pulses[1];

        if(mode <= MODE_515) {
            pulse_nb = 2;
            pulse_subset = (fixed_index & 0x40)>>6;
            pulse_position[0] = ( fixed_index       & 7)*5 + track_position[ (pulse_subset<<3) + (subframe<<1) ];
            pulse_position[1] = ((fixed_index >> 3) & 7)*5 + track_position[ (pulse_subset<<3) + (subframe<<1) + 1 ];
        }else if(mode == MODE_59) {
            pulse_nb = 2;
            pulse_subset = fixed_index & 1;
            pulse_position[0] = ((fixed_index >> 1) & 7)*5 + (pulse_subset<<1) + 1;
            pulse_subset = (fixed_index >> 4) & 3;
            pulse_position[1] = ((fixed_index >> 6) & 7)*5 + pulse_subset + (pulse_subset == 3 ? 1 : 0);
        }else if(mode == MODE_67) {
            pulse_nb = 3;
            pulse_position[0] = ( fixed_index       & 7)*5;
            pulse_subset = (fixed_index >> 3) & 1;
            pulse_position[1] = ((fixed_index >> 4) & 7)*5 + (pulse_subset<<1) + 1;
            pulse_subset = (fixed_index >> 7) & 1;
            pulse_position[2] = ((fixed_index >> 8) & 7)*5 + (pulse_subset<<1) + 2;
        }else { // mode <= MODE_795
            pulse_nb = 4;
            pulse_position[0] = gray_decode[ fixed_index        & 7]*5;
            pulse_position[1] = gray_decode[(fixed_index >> 3)  & 7]*5 + 1;
            pulse_position[2] = gray_decode[(fixed_index >> 6)  & 7]*5 + 2;
            pulse_subset = (fixed_index >> 9) & 1;
            pulse_position[3] = gray_decode[(fixed_index >> 10) & 7]*5 + pulse_subset + 3;
        }
        reconstruct_fixed_vector(pulse_position, sign, pulse_nb, fixed_vector);
    }
}

/// @}


/// @defgroup amr_gain_decoding gain decoding functions
/// @{

/**
 * Predict the fixed gain.
 *
 * @param fixed_vector         pointer to the algebraic codebook vector
 * @param prev_pred_error      pointer to the quantified prediction errors from the previous four subframes
 *
 * @return the predicted fixed gain
 */

static float fixed_gain_prediction(float *fixed_vector, float *prev_pred_error,
                                   enum Mode mode)
{
    int i;
    float energy_pred = 0.0, energy_fixed_mean = 0.0;

    // Calculate the predicted energy
    for(i=0; i<4; i++) {
        energy_pred += energy_pred_fac[i]*prev_pred_error[3-i];
    }

    // Calculate the mean fixed vector energy
    for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
        energy_fixed_mean += fixed_vector[i]*fixed_vector[i];
    }
    energy_fixed_mean = 10.0*log10f(energy_fixed_mean/(float)AMR_SUBFRAME_SIZE);

    // predicted fixed gain =
    // 10^(0.05 * (predicted energy + desired mean energy - mean fixed vector energy))
    return powf(10.0, 0.05*(energy_pred + energy_mean[mode] - energy_fixed_mean));
}

/// @}


/// @defgroup amr_pre_processing pre-processing functions
/// @{

/**
 * Comparison function for use with qsort.
 *
 * @param a             first value for comparison
 * @param b             second value for comparison
 * @return a-b : the result of the comparison
 */

int qsort_compare(const void *a, const void *b)
{
    float diff = *(const float *)a - *(const float *)b;
    if(diff > 0.0f)
        return 1;
    if(diff < 0.0f)
        return -1;
    return 0;
}

/**
 * Find the median of some float values.
 *
 * @param values        pointer to the values of which to find the median
 * @param n             number of values
 * @return the median value
 */

static float medianf(float *values, int n)
{
    float temp[9]; // largest n used for median calculation is 9

    memcpy(temp, values, n * sizeof(float));

    qsort(temp, n, sizeof(float), qsort_compare);

    if(n&1) {
        return                     temp[ n>>1 ];
    }else {
        return (temp[ (n>>1)-1 ] + temp[ n>>1 ])/2.0;
    }
}

/**
 * Circularly convolve the fixed vector with a phase dispersion impulse response
 * filter.
 *
 * @param fixed_vector  pointer to the fixed vector
 * @param ir_filter     pointer to the impulse response filter
 */

static void convolve_circ(float *fixed_vector, const float *ir_filter)
{
    int i, j, k;
    int npulses = 0, pulse_positions[AMR_SUBFRAME_SIZE];
    float fixed_vector_temp[AMR_SUBFRAME_SIZE];

    memcpy(fixed_vector_temp, fixed_vector, AMR_SUBFRAME_SIZE*sizeof(float));
    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE*sizeof(float));

    // Find non-zero pulses (most are zero)
    for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
        if(fixed_vector_temp[i]) {
            pulse_positions[npulses] = i;
            npulses++;
        }
    }

    for(i=0; i<npulses; i++) {
        k = 0;
        for(j=pulse_positions[i]; j<AMR_SUBFRAME_SIZE; j++) {
            fixed_vector[j] += fixed_vector_temp[pulse_positions[i]]*ir_filter[k++];
        }
        for(j=0; j<pulse_positions[i]; j++) {
            fixed_vector[j] += fixed_vector_temp[pulse_positions[i]]*ir_filter[k++];
        }
    }
}

/// @}


/// @defgroup amr_synthesis synthesis functions
/// @{

/**
 * Conduct 10th order linear predictive coding synthesis.
 *
 * @param p             pointer to the AMRContext
 * @param excitation    pointer to the excitation vector
 * @param lpc           pointer to the LPC coefficients
 * @param samples       pointer to the output speech samples
 * @param overflow      16-bit overflow flag
 */

static int synthesis(AMRContext *p, float *excitation, float *lpc,
                     float *samples, uint8_t overflow)
{
    int i, j, overflow_temp = 0;

    // if an overflow has been detected, the pitch vector is scaled down by a
    // factor of 4
    if(overflow) {
        for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
            p->pitch_vector[i] /= 4.0;
        }
    }

    // construct the excitation vector
    for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
        excitation[i] = p->pitch_gain[4]*p->pitch_vector[i] + p->fixed_gain[4]*p->fixed_vector[i];
    }

    // if an overflow has been detected, pitch vector contribution emphasis and
    // adaptive gain control are skipped
    if(p->pitch_gain[4] > 0.5 && !overflow) {
        float excitation_temp[AMR_SUBFRAME_SIZE];
        float pitch_factor = (p->cur_frame_mode == MODE_122 ? 0.25 : 0.5)*p->beta*p->pitch_gain[4];
        float eta, temp1 = 0.0, temp2 = 0.0;

        for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
            // emphasize pitch vector contribution
            excitation_temp[i] = excitation[i] + pitch_factor*p->pitch_vector[i];
            // find gain scale
            temp1 +=      excitation[i]*excitation[i];
            temp2 += excitation_temp[i]*excitation_temp[i];
        }

        // adaptive gain control by gain scaling
        eta = sqrt(temp1/temp2);
        for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
            excitation[i] = eta*excitation_temp[i];
        }
    }

    for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
        samples[i] = excitation[i];
        for(j=0; j<LP_FILTER_ORDER; j++) {
            samples[i] -= lpc[j]*samples[i-j-1];
        }
        // detect overflow
        if(fabsf(samples[i])>1.0) {
            overflow_temp = 1;
            samples[i] = av_clipf(samples[i], -1.0, 1.0);
        }
    }

    return overflow_temp;
}

/// @}


/// @defgroup amr_update update functions
/// @{

/**
 * Update buffers and history at the end of decoding a subframe.
 *
 * @param p             pointer to the AMRContext
 */

static void update_state(AMRContext *p)
{
    // update the previous frame's fourth subframe LSP vector
    memcpy(p->prev_lsp_sub4, p->lsp[3], LP_FILTER_ORDER * sizeof(float));

    // update the excitation buffer moving the current values into the buffer
    // pushing out those no longer needed
    memmove(&p->excitation_buf[0], &p->excitation_buf[AMR_SUBFRAME_SIZE],
        (PITCH_LAG_MAX + LP_FILTER_ORDER + 1)*sizeof(float));

    // update quantified prediction error energy history
    p->prediction_error[0] = p->prediction_error[1];
    p->prediction_error[1] = p->prediction_error[2];
    p->prediction_error[2] = p->prediction_error[3];
    p->prediction_error[3] = 20.0*log10f(p->fixed_gain_factor);

    // update pitch lag history
    p->prev_pitch_lag_int = p->pitch_lag_int;

    // update gain history
    memmove(&p->pitch_gain[0], &p->pitch_gain[1], 4*sizeof(float));
    memmove(&p->fixed_gain[0], &p->fixed_gain[1], 4*sizeof(float));

    // update ir filter strength history
    p->ir_filter_strength[0] = p->ir_filter_strength[1];

    // update speech sample history
    memmove(&p->samples_in[0], &p->samples_in[AMR_SUBFRAME_SIZE],
        LP_FILTER_ORDER*sizeof(float));
}

/// @}



static int amrnb_decode_frame(AVCodecContext *avctx, void *data, int *data_size,
                              const uint8_t *buf, int buf_size)
{

    AMRContext *p = avctx->priv_data;        // pointer to private data
    float *buf_out = data;                   // pointer to the output data buffer
    int i, subframe;                         // counters
    int gains_index_MODE_475 = 0;            // MODE_475 gains index coded every other subframe
    enum Mode speech_mode = MODE_475;        // ???

    // decode the bitstream to AMR parameters
    p->cur_frame_mode = decode_bitstream(p, buf, buf_size, &speech_mode);
    if(p->cur_frame_mode == MODE_DTX) {
        ff_log_missing_feature(avctx, "dtx mode", 1);
        return -1;
    }
/*** LPC coefficient decoding ***/

    if(p->cur_frame_mode == MODE_122) {
        // decode split-matrix quantized lsf vector indexes to lsp vectors
        lsf2lsp_5(p);
    }else {
        // decode split-matrix quantized lsf vector indexes to an lsp vector
        lsf2lsp_3(p);
    }

    // convert LSP vectors to LPC coefficient vectors
    for(i=0; i<4; i++) {
        lsp2lpc(p->lsp[i], p->lpc[i]);
    }

    // update averaged lsp vector (used for fixed gain smoothing)
    weighted_vector_sumf(p->lsp_avg, p->lsp_avg, p->prev_lsp_sub4, 0.84, 0.16, LP_FILTER_ORDER);

/*** end of LPC coefficient decoding ***/

    for(subframe = 0; subframe < 4; subframe++) {
        const AMRNBSubframe *amr_subframe = &p->frame.subframe[subframe];
/*** adaptive code book (pitch) vector decoding ***/

        // find the search range
        p->search_range_min = FFMAX(p->prev_pitch_lag_int - 5, p->cur_frame_mode == MODE_122 ? PITCH_LAG_MIN_MODE_122 : PITCH_LAG_MIN);
        p->search_range_max = p->search_range_min + 9;
        if(p->search_range_max > PITCH_LAG_MAX) {
            p->search_range_max = PITCH_LAG_MAX;
            p->search_range_min = p->search_range_max - 9;
        }

        // decode integer and fractional parts of pitch lag from parsed pitch
        // index
        if(p->cur_frame_mode == MODE_122) {
            decode_pitch_lag_6(p, amr_subframe->p_lag, subframe);
        }else {
            decode_pitch_lag_3(p, amr_subframe->p_lag, subframe);
        }

        // interpolate the past excitation at the pitch lag to obtain the pitch
        // vector
        interp_pitch_vector(p->excitation, p->pitch_lag_int, p->pitch_lag_frac, p->cur_frame_mode, p->pitch_vector);

/*** end of adaptive code book (pitch) vector decoding ***/

        decode_fixed_vector(p->fixed_vector, amr_subframe->pulses, p->cur_frame_mode, subframe);

/*** gain decoding ***/

        // calculate the predicted fixed gain g_c'
        p->fixed_gain[4] = fixed_gain_prediction(p->fixed_vector, p->prediction_error, p->cur_frame_mode);

        // decode pitch gain and fixed gain correction factor
        if(p->cur_frame_mode == MODE_122 || p->cur_frame_mode == MODE_795) {
            p->pitch_gain[4] =     qua_gain_pit[amr_subframe->p_gain];
            p->fixed_gain_factor = qua_gain_code[amr_subframe->fixed_gain];
        }else if(p->cur_frame_mode == MODE_67 || p->cur_frame_mode == MODE_74 ||
                 p->cur_frame_mode == MODE_102) {
            p->pitch_gain[4] =     gains_high[amr_subframe->p_gain][0];
            p->fixed_gain_factor = gains_high[amr_subframe->p_gain][1];
        }else if(p->cur_frame_mode == MODE_515 || p->cur_frame_mode == MODE_59) {
            p->pitch_gain[4] =     gains_low[amr_subframe->p_gain][0];
            p->fixed_gain_factor = gains_low[amr_subframe->p_gain][1];
        }else {
            // gain index is only coded in subframes 0,2
            if(!(subframe&1)) {
                gains_index_MODE_475 = amr_subframe->p_gain<<1;
            }
            p->pitch_gain[4] =     gains_MODE_475[gains_index_MODE_475 + (subframe&1)][0];
            p->fixed_gain_factor = gains_MODE_475[gains_index_MODE_475 + (subframe&1)][1];
        }

        // ^g_c = g_c' * ^gamma_gc
        p->fixed_gain[4] *= p->fixed_gain_factor;

/*** end of gain decoding ***/

/*** pre-processing ***/

        p->beta = av_clipf(p->pitch_gain[4], 0.0, p->cur_frame_mode == MODE_122 ? 1.0 : 0.8);

        // conduct pitch sharpening as appropriate
        if(p->pitch_lag_int < AMR_SUBFRAME_SIZE) {
            for(i=p->pitch_lag_int; i<AMR_SUBFRAME_SIZE; i++) {
                p->fixed_vector[i] += p->beta*p->fixed_vector[i-p->pitch_lag_int];
            }
        }

        // smooth fixed gain
        if(p->cur_frame_mode < MODE_74 || p->cur_frame_mode == MODE_102) {
            float diff = 0.0;
            float smoothing_factor = 0.0;

            for(i=0; i<LP_FILTER_ORDER; i++) {
                // calculate diff
                diff += fabs(p->lsp_avg[i]-p->lsp[subframe][i])/p->lsp_avg[i];
            }

            // if diff has been >0.65 for 10 frames (40 subframes) no smoothing is applied
            if(diff > 0.65) {
                p->diff_count++;
            }else {
                p->diff_count = 0;
            }

            if(p->diff_count < 40) {
                float fixed_gain_mean;
                // calculate the fixed gain smoothing factor (k_m)
                smoothing_factor = FFMIN(0.25, FFMAX(0.0, diff - 0.4))/0.25;
                // calculate the mean fixed gain for the current subframe
                fixed_gain_mean = (p->fixed_gain[0] + p->fixed_gain[1] + p->fixed_gain[2] + p->fixed_gain[3] + p->fixed_gain[4])/5.0;
                // calculate the smoothed fixed gain
                p->fixed_gain[4] = smoothing_factor*p->fixed_gain[4] + (1.0 - smoothing_factor)*fixed_gain_mean;
            }
        }

        // anti-sparseness processing
        if(p->pitch_gain[4] < 0.6) {
            // strong filtering
            p->ir_filter_strength[1] = 0;
        }else if(p->pitch_gain[4] < 0.9) {
            // medium filtering
            p->ir_filter_strength[1] = 1;
        }else {
            // no filtering
            p->ir_filter_strength[1] = 2;
        }

        // detect 'onset'
        if(p->fixed_gain[4] > 2.0*p->fixed_gain[3]) {
            p->ir_filter_strength[1] = FFMIN(p->ir_filter_strength[1] + 1, 2);
        }else if(p->ir_filter_strength[1] == 0 && medianf(p->pitch_gain, 5) >= 0.6 &&
                    p->ir_filter_strength[1] > p->ir_filter_strength[0] + 1) {
            p->ir_filter_strength[1] = p->ir_filter_strength[0] + 1;
        }

        if(p->cur_frame_mode != MODE_74 && p->cur_frame_mode != MODE_102 &&
                p->cur_frame_mode != MODE_122 && p->ir_filter_strength[1] < 2) {
            // assign the correct impulse response
            if(p->ir_filter_strength[1] == 1) {
                p->ir_filter = ir_filter_medium;
            }else {
                if(p->cur_frame_mode != MODE_795) {
                    p->ir_filter = ir_filter_strong;
                }else {
                    p->ir_filter = ir_filter_strong_MODE_795;
                }
            }

            // circularly convolve the fixed vector with the impulse response
            convolve_circ(p->fixed_vector, p->ir_filter);
        }

/*** end of pre-processing ***/

/*** synthesis ***/

        if(synthesis(p, p->excitation, p->lpc[subframe], &p->samples_in[LP_FILTER_ORDER], 0)) {
            // overflow detected -> rerun synthesis scaling pitch vector down by
            // a factor of 4, skipping pitch vector contribution emphasis and
            // adaptive gain control
            synthesis(p, p->excitation, p->lpc[subframe], &p->samples_in[LP_FILTER_ORDER], 1);
        }

/*** end of synthesis ***/

        // update buffers and history
        update_state(p);

        memcpy(&buf_out[subframe*AMR_SUBFRAME_SIZE], &p->samples_in[LP_FILTER_ORDER],
               AMR_SUBFRAME_SIZE*sizeof(float));
    }

    /* report how many samples we got */
    *data_size = AMR_BLOCK_SIZE * sizeof(float);

    /* return the amount of bytes consumed if everything was OK */
    return (mode_bits[p->cur_frame_mode] + 15)>>3; // +7 for rounding and +8 for TOC
}


AVCodec amrnb_decoder = {
    .name = "amrnb",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_AMR_NB,
    .priv_data_size = sizeof(AMRContext),
    .init = amrnb_decode_init,
    .decode = amrnb_decode_frame,
    .long_name = NULL_IF_CONFIG_SMALL("Adaptive Multi-Rate NarrowBand"),
};

