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
#include "common.h"
#include "amrnbfloatdata.h"

typedef struct AMRContext {

    GetBitContext                        gb;
    int16_t                  *sample_buffer;

    int16_t                       *amr_prms; ///< pointer to the decoded amr parameters (lsf coefficients, codebook indices, etc)
    int                 bad_frame_indicator; ///< bad frame ? 1 : 0
    int                      cur_frame_mode; ///< current frame mode
    int                      cur_frame_type; ///< current frame type

    float       prev_lsf_r[LP_FILTER_ORDER];
    float           lsp[4][LP_FILTER_ORDER]; ///< lsp vectors from current frame
    float    prev_lsp_sub4[LP_FILTER_ORDER]; ///< lsp vector for the 4th subframe of the previous frame

    float          lsp_avg[LP_FILTER_ORDER]; ///< vector of averaged lsp coefficients

    float           lpc[4][LP_FILTER_ORDER]; ///< vectors of lpc coefficients for 4 subframes

    int                       pitch_lag_int; ///< integer part of pitch lag from current subframe
    int                      pitch_lag_frac; ///< fractional part of pitch lag from current subframe
    int                  prev_pitch_lag_int; ///< integer part of pitch lag from previous subframe

    float excitation_buf[PITCH_LAG_MAX + LP_FILTER_ORDER + 1 + AMR_SUBFRAME_SIZE]; ///< excitation buffer
    float                       *excitation; ///< pointer to the current excitation vector in excitation_buf

    float                      pitch_vector[AMR_SUBFRAME_SIZE]; ///< adaptive code book (pitch) vector

    float                      fixed_vector[AMR_SUBFRAME_SIZE]; ///< algebraic code book (fixed) vector

    float               prediction_error[4]; ///< quantified prediction errors {20log10(^γ_gc)} for previous four subframes
    float                     pitch_gain[5]; ///< quantified pitch gains for the current and previous four subframes
    float                 fixed_gain_factor;
    float                     fixed_gain[5]; ///< quantified fixed gains for the current and previous four subframes

    float                              beta; ///< beta = pitch_gain, bounded by [0.0,1.0] for 12.2 kbps or [0.0,0.8] for other modes
    int                          diff_count; ///< the number of subframes for which diff has been above 0.65

    uint8_t           ir_filter_strength[2]; ///< impulse response filter strength; 0 - strong, 1 - medium, 2 - none
    float                        *ir_filter; ///< pointer to impulse response filter data

} AMRContext;


static int amrnb_decode_init(AVCodecContext *avctx) {
    AMRContext *p = avctx->priv_data;

    // allocate and zero the 16-bit mono sample buffer
    p->sample_buffer = av_mallocz(sizeof(int16_t)*AMR_BLOCK_SIZE);
    // allocate and zero the amr parameters
    p->amr_prms = av_mallocz(sizeof(int16_t)*PRMS_MODE_122);

    /* Check if the allocation was successful */
    if(p->sample_buffer == NULL)
        return -1;
    // Check amr_prms allocation
    if(p->amr_prms == NULL)
        return -1;

    // p->excitation always points to the same position in p->excitation_buf
    p->excitation = &p->excitation_buf[PITCH_LAG_MAX + LP_FILTER_ORDER + 1];

    /* return 0 for a successful init, -1 for failure */
    return 0;
}


/**
 * Decode the bitstream into the AMR parameters and discover the frame mode
 *
 * @param buf               pointer to the input buffer
 * @param buf_size          size of the input buffer
 * @param speech_mode       pointer to the speech mode
 *
 * @return Returns the frame mode
 */

enum Mode decode_bitstream(AMRContext *p, uint8_t *buf, int buf_size, enum Mode *speech_mode) {
    enum Mode mode;
    int i;
    const AMROrder *order;

    // initialise get_bits
    init_get_bits(&p->gb, buf, buf_size*8);
    skip_bits1(&p->gb);
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
        for(i=0; i<mode_bits[mode]; i++) {
            p->amr_prms[ order[i].array_element ] += get_bits1(&p->gb) * (1<< order[i].bit_mask);
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


/*** LPC coefficient decoding functions ***/

/**
 * Convert an lsf vector into an lsp vector
 *
 * @param lsf               input lsf vector
 * @param lsp               output lsp vector
 *
 * @return void
 */

static void lsf2lsp(float *lsf, float *lsp) {
    int i;

    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsp[i] = cos(lsf[i]*FREQ_LSP_FAC); // FREQ_LSP_FAC = 2*M_PI/8000.0
    }
}

/**
 * Decode a set of 5 split-matrix quantised lsf indices into an lsp vector
 *
 * @param p                 pointer to the AMRContext
 *
 * @return void
 */

static void lsf2lsp_5(AMRContext *p) {
    float lsf_r[2][LP_FILTER_ORDER]; // residual LSF vectors
    float lsf_q[2][LP_FILTER_ORDER]; // quantified LSF vectors
    float sign;
    int i, idx;

    // decode split-matrix quantized residual LSF vectors

    idx = p->amr_prms[0];
    lsf_r[0][0] = lsf_5_1[ idx ][0];
    lsf_r[0][1] = lsf_5_1[ idx ][1];
    lsf_r[1][0] = lsf_5_1[ idx ][2];
    lsf_r[1][1] = lsf_5_1[ idx ][3];

    idx = p->amr_prms[1];
    lsf_r[0][2] = lsf_5_2[ idx ][0];
    lsf_r[0][3] = lsf_5_2[ idx ][1];
    lsf_r[1][2] = lsf_5_2[ idx ][2];
    lsf_r[1][3] = lsf_5_2[ idx ][3];

    // lsb of p->amr_prms[2] is the sign bit
    sign = (p->amr_prms[2] & 1) ? -1.0 : 1.0;
    idx = p->amr_prms[2]>>1;
    lsf_r[0][4] = lsf_5_3[ idx ][0]*sign;
    lsf_r[0][5] = lsf_5_3[ idx ][1]*sign;
    lsf_r[1][4] = lsf_5_3[ idx ][2]*sign;
    lsf_r[1][5] = lsf_5_3[ idx ][3]*sign;

    idx = p->amr_prms[3];
    lsf_r[0][6] = lsf_5_4[ idx ][0];
    lsf_r[0][7] = lsf_5_4[ idx ][1];
    lsf_r[1][6] = lsf_5_4[ idx ][2];
    lsf_r[1][7] = lsf_5_4[ idx ][3];

    idx = p->amr_prms[4];
    lsf_r[0][8] = lsf_5_5[ idx ][0];
    lsf_r[0][9] = lsf_5_5[ idx ][1];
    lsf_r[1][8] = lsf_5_5[ idx ][2];
    lsf_r[1][9] = lsf_5_5[ idx ][3];

    // calculate mean-removed LSF vectors and add mean
    for(i=0; i<LP_FILTER_ORDER; i++) {
        float temp = p->prev_lsf_r[i]*PRED_FAC_MODE_122 + lsf_5_mean[i];
        lsf_q[0][i] = lsf_r[0][i] + temp;
        lsf_q[1][i] = lsf_r[1][i] + temp;
    }
    // update residual LSD vector from previous subframe
    memcpy(p->prev_lsf_r, lsf_r[1], LP_FILTER_ORDER*sizeof(float));

    // convert LSF vectors to LSP vectors
    lsf2lsp(lsf_q[0], p->lsp[1]);
    lsf2lsp(lsf_q[1], p->lsp[3]);
}

/**
 * Decode a set of 3 split-matrix quantised lsf indices into an lsp vector
 *
 * @param p                 pointer to the AMRContext
 *
 * @return void
 */

static void lsf2lsp_3(AMRContext *p) {
    float lsf_r[LP_FILTER_ORDER]; // residual LSF vector
    float lsf_q[LP_FILTER_ORDER]; // quantified LSF vector
    const float (*lsf_3_1_tmp)[3], (*lsf_3_3_tmp)[4]; // temp ptrs for switching tables depending on mode
    float sign;
    int idx, i;

    // assign lsf tables according to mode
    if((p->cur_frame_mode == MODE_475) || (p->cur_frame_mode == MODE_515)) {
        lsf_3_1_tmp = lsf_3_1;
        lsf_3_3_tmp = lsf_3_MODE_515;
    }else if(p->cur_frame_mode == MODE_795) {
        lsf_3_1_tmp = lsf_3_MODE_795;
        lsf_3_3_tmp = lsf_3_3;
    }else {
        lsf_3_1_tmp = lsf_3_1;
        lsf_3_3_tmp = lsf_3_3;
    }

    // decode split-matrix quantized residual LSF vector

    idx = p->amr_prms[0];
    lsf_r[0] = lsf_3_1_tmp[ idx ][0];
    lsf_r[1] = lsf_3_1_tmp[ idx ][1];
    lsf_r[2] = lsf_3_1_tmp[ idx ][2];

    idx = p->amr_prms[1];
    // MODE_475, MODE_515 only use every other entry as their indices are stored
    // using 1 less bit (8-bits vs 9-bits)
    if((p->cur_frame_mode == MODE_475) || (p->cur_frame_mode == MODE_515)) {
        idx <<= 1;
    }
    lsf_r[3] = lsf_3_2[ idx ][0];
    lsf_r[4] = lsf_3_2[ idx ][1];
    lsf_r[5] = lsf_3_2[ idx ][2];

    lsf_r[6] = lsf_3_3_tmp[ idx ][0];
    lsf_r[7] = lsf_3_3_tmp[ idx ][1];
    lsf_r[8] = lsf_3_3_tmp[ idx ][2];
    lsf_r[9] = lsf_3_3_tmp[ idx ][3];

    // calculate mean-removed LSF vector and add mean
    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsf_q[i] = lsf_r[i] + p->prev_lsf_r[i]*pred_fac[i] + lsf_3_mean[i];
    }
    // update residual LSF vector from previous subframe
    memcpy(p->prev_lsf_r, lsf_r, LP_FILTER_ORDER*sizeof(float));

    // convert LSF vector to LSP vector
    lsf2lsp(lsf_q, p->lsp[3]);
}

/**
 * Interpolate lsp vectors for subframes 1 and 3
 *
 * @param p                 pointer to the AMRContext
 *
 * @return void
 */

static void interp_lsp_13(AMRContext *p) {
    int i;

    for(i=0; i<LP_FILTER_ORDER; i++) {
        p->lsp[0][i] = 0.5*(p->prev_lsp_sub4[i] + p->lsp[1][i]);
        p->lsp[2][i] = 0.5*(       p->lsp[1][i] + p->lsp[3][i]);
    }
}

/**
 * Interpolate lsp vectors for subframes 1, 2 and 3
 *
 * @param p                 pointer to the AMRContext
 *
 * @return void
 */

static void interp_lsp_123(AMRContext *p) {
    int i;

    for(i=0; i<LP_FILTER_ORDER; i++) {
        p->lsp[0][i] = 0.75*p->prev_lsp_sub4[i] + 0.25*p->lsp[3][i];
        p->lsp[1][i] = 0.5*(p->prev_lsp_sub4[i] +      p->lsp[3][i]);
        p->lsp[2][i] = 0.25*p->prev_lsp_sub4[i] + 0.75*p->lsp[3][i];
    }
}

/**
 * Find the polynomial F1(z) or F2(z) from the lsp vectors
 *
 * @param lsp               input lsp vector
 * @param f                 pointer to the polynomial F1(z) or F2(z)
 *
 * @return void
 */

static void lsp2poly(float *lsp, float *f) {
    int i, j;

    f[0] = 0.0;
    f[1] = 1.0;

    for(i=2; i<7; i++) {
        int idx_lsp = 2*i-4;
        f[i] = 2.0*(-lsp[idx_lsp]*f[i-1] + f[i-2]);
        for(j=i-1; j>0; j--) {
            f[j] += -2.0*lsp[idx_lsp]*f[j-1] + f[j-2];
        }
    }
}

/**
 * Convert an lsp vector to lpc coefficients
 *
 * @param lsp                 input lsp vector
 * @param lpc                 output lpc coefficients
 *
 * @return void
 */

static void lsp2lpc(float *lsp, float *lpc_coeffs) {
    float f1[6], f2[6];
    int i;

    // find F1(z) and F2(z) from the lsps
    lsp2poly(&lsp[0], f1);
    lsp2poly(&lsp[1], f2);

    // multiply F1(z) by 1+z^{-1} and F2(z) by 1-z^{-1} to obtain F1'(z) and F2'(z)
    for(i=5; i>0; i--) {
        f1[i] += f1[i-1];
        f2[i] -= f2[i-1];
    }

    // A(z) = ( F1'(z) + F2'(z) )/2
    // note f1 and f2 are actually f1' and f2'
    for(i=0; i<5; i++) {
        lpc_coeffs[i]   = 0.5*(f1[i+1] + f2[i+1]); // lpc 0..4 uses indices to f, 1..5
        lpc_coeffs[i+5] = 0.5*(f1[5-i] - f2[5-i]); // lpc 5..9 uses indices to f, 5..1
    }
}

/*** end of LPC coefficient decoding functions ***/


/*** pitch vector decoding functions ***/

/**
 * Decode the adaptive codebook index to the integer and fractional parts of the
 * pitch lag for one subframe at 1/3 resolution
 *
 * @param p                   pointer to the AMRContext
 * @param pitch_index         parsed adaptive codebook (pitch) index
 *
 * @return void
 */

static void decode_pitch_lag_3(AMRContext *p, int pitch_index) {
    // subframe 1 or 3
    if(p->cur_subframe & 1) {
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
            // decoding with 4 bit resolution
            int t1_temp = clip(p->prev_pitch_lag_int, p->search_range_max-4, p->search_range_min+5);

            if(pitch_index < 4) {
                // integer only precision for [t1_temp-5, t1_temp-2]
                p->pitch_lag_int = pitch_index + (t1_temp - 5);
                p->pitch_lag_frac = 0;
            }else if(pitch_index < 12) {
                // 1/3 fractional precision for [t1_temp-1 2/3, t1_temp+2/3]
                p->pitch_lag_int = ( ((pitch_index - 5)*10923)>>15 ) + t1_temp - 1;
                p->pitch_lag_frac = pitch_index - p->pitch_lag_int*3 - 9;
            }else {
                // integer only precision for [t1_temp+1, t1_temp+4]
                pitch_lag_int = pitch_index + t1_temp - 11;
                pitch_lag_frac = 0;
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
 * Decode the adaptive codebook index to the integer and fractional parts of the
 * pitch lag for one subframe at 1/6 resolution
 *
 * @param p                   pointer to the AMRContext
 * @param pitch_index         parsed adaptive codebook (pitch) index
 *
 * @return void
 */

static void decode_pitch_lag_6(AMRContext *p, int pitch_index) {
    // subframe 1 or 3
    if(p->cur_subframe & 1) {
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
        // find the search range
        p->search_range_min = FFMAX(p->pitch_lag_int - 5, PITCH_LAG_MIN_MODE_122);
        p->search_range_max = p->search_range_min + 9;
        if(p->search_range_max > PITCH_LAG_MAX) {
            p->search_range_max = PITCH_LAG_MAX;
            p->search_range_min = p->search_range_max - 9;
        }
        // calculate the pitch lag
        temp = (pitch_index + 5)/6 - 1;
        p->pitch_lag_int = temp + p->search_range_min;
        p->pitch_lag_frac = pitch_index - temp*6 - 3;
    }
}

/**
 * Calculate the pitch vector by interpolating the past excitation at the pitch
 * pitch lag using a b60 hamming windowed sinc function
 *
 * @param prev_excitation     pointer to the element after the previous excitations
 * @param lag_int             integer part of pitch lag
 * @param lag_frac            fractional part of pitch lag
 * @param mode                current frame mode
 * @param pitch_vector        pointer to the pitch vector
 *
 * @return void
 */

static void interp_pitch_vector(float *prev_excitation, int lag_int, int lag_frac, enum Mode mode, float *pitch_vector) {
    int n, i;
    float *b60_idx1, *b60_idx2, *exc_idx;

    lag_frac *= -1;
    if(mode != MODE_122) {
        lag_frac <<= 1;
    }

    if(lag_frac < 0) {
        lag_frac += 6;
        lag_int--;
    }

    b60_idx1 = &b60[    lag_frac];
    b60_idx2 = &b60[6 - lag_frac];
    exc_idx = &prev_excitation[-lag_int];

    for(n=0; n<AMR_SUBFRAME_SIZE; n++) {
        for(i=0; i<10; i++) {
            pitch_vector[n] += b60_idx1[6*i] * exc_idx[-i];
        }
        exc_idx++;
        for(i=0; i<10; i++) {
            pitch_vector[n] += b60_idx2[6*i] * exc_idx[ i];
        }
        exc_idx++;
    }
}

/*** end of pitch vector decoding functions ***/


/*** algebraic code book (fixed) vector decoding functions ***/

/**
 * Reconstruct the algebraic codebook vector
 *
 * @param pulse_position       vector of pulse positions
 * @param sign                 signs of the pulses
 * @param nr_pulses            number of pulses
 * @param fixed_vector         algebraic codebook vector
 *
 * @return void
 */

static void reconstruct_fixed_vector(float *pulse_position, int sign, int nr_pulses, float *fixed_vector) {
    int i;

    // reset the code
    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE*sizeof(float));

    for(i=0; i<nr_pulses; i++)
        fixed_vector[pulse_position[i]] = ((sign >> i) & 1) ? 1.0 : -1.0;
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and construct
 * the algebraic codebook vector for MODE_475 and MODE_515
 *
 * @param fixed_index          positions of the two pulses
 * @param sign                 signs of the two pulses
 * @param subframe             current subframe
 * @param fixed_vector         pointer to the algebraic codebook vector
 *
 * @return void
 */

static void decode_2_pulses_9bits(int fixed_index, int sign, int subframe, float *fixed_vector) {
    int pulse_position[2];
    int pulse_subset;

    // pulse subset is the msb (bit 7) of the 7-bits used to code the 2 pulses
    pulse_subset = (fixed_index & 0x40)>>6;
    // first pulse position is coded in bits 1-3
    pulse_position[0] = ( fixed_index       & 7)*5 + track_position[ (pulse_subset<<3) + (subframe<<1) ];
    // second pulse position is coded in bits 4-6
    pulse_position[1] = ((fixed_index >> 3) & 7)*5 + track_position[ (pulse_subset<<3) + (subframe<<1) + 1 ];

    // reconstruct the fixed code
    reconstruct_fixed_vector(pulse_position, sign, 2, fixed_vector);
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and construct
 * the algebraic codebook vector for MODE_59
 *
 * @param fixed_index          positions of the two pulses
 * @param sign                 signs of the two pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 *
 * @return void
 */

static void decode_2_pulses_11bits(int fixed_index, int sign, float *fixed_vector) {
    int pulse_position[2];
    int pulse_subset;

    // pulse subset for the first pulse is the lsb (bit 1) of the 9-bits used
    // to code the 2 pulses
    pulse_subset = fixed_index & 1;
    // first pulse position is coded in bits 2-4
    pulse_position[0] = ((fixed_index >> 1) & 7)*5 + pulse_subset<<1 + 1;
    // pulse subset for the second pulse is coded in bits 5-6
    pulse_subset = (fixed_index >> 4) & 3;
    // second pulse position is coded in bits 7-9
    pulse_position[1] = ((fixed_index >> 6) & 7)*5 + pulse_subset + (pulse_subset == 3 ? 1 : 0);

    // reconstruct the fixed code
    reconstruct_fixed_vector(pulse_position, sign, 2, fixed_vector);
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and construct
 * the algebraic codebook vector for MODE_67
 *
 * @param fixed_index          positions of the three pulses
 * @param sign                 signs of the three pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 *
 * @return void
 */

static void decode_3_pulses_14bits(int fixed_index, int sign, float *fixed_vector) {
    int pulse_position[3];
    int pulse_subset;

    // first pulse position is coded in bits 1-3
    pulse_position[0] = ( fixed_index       & 7)*5;
    // pulse subset for the second pulse is coded in bit 4
    pulse_subset = (fixed_index >> 3) & 1;
    // second pulse position is coded in bits 5-7
    pulse_position[1] = ((fixed_index >> 4) & 7)*5 + pulse_subset<<1 + 1;
    // pulse subset for the second pulse is coded in bit 8
    pulse_subset = (fixed_index >> 7) & 1;
    // third pulse position is coded in bits 9-11
    pulse_position[2] = ((fixed_index >> 8) & 7)*5 + pulse_subset<<1 + 2;

    // reconstruct the fixed code
    reconstruct_fixed_vector(pulse_position, sign, 3, fixed_vector);
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and construct
 * the algebraic codebook vector for MODE_74 and MODE_795
 *
 * @param fixed_index          positions of the four pulses
 * @param sign                 signs of the four pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 *
 * @return void
 */

static void decode_4_pulses_17bits(int fixed_index, int sign, float *fixed_vector) {
    int pulse_position[4];
    int pulse_subset;

    // first pulse position is Gray coded in bits 1-3
    pulse_position[0] = gray_decode[ fixed_index        & 7]*5;
    // second pulse position is Gray coded in bits 4-6
    pulse_position[1] = gray_decode[(fixed_index >> 3)  & 7]*5 + 1;
    // third pulse position is Gray coded in bits 7-9
    pulse_position[2] = gray_decode[(fixed_index >> 6)  & 7]*5 + 2;
    // pulse subset for the fourth pulse is coded in bit 10
    pulse_subset = (fixed_index >> 9) & 1;
    // third pulse position is Gray coded in bits 11-13
    pulse_position[3] = gray_decode[(fixed_index >> 10) & 7]*5 + pulse_subset + 3;

    // reconstruct the fixed code
    reconstruct_fixed_vector(pulse_position, sign, 4, fixed_vector);
}

/**
 * Decode the algebraic codebook index to pulse positions and signs and construct
 * the algebraic codebook vector for MODE_102
 *
 * @param fixed_index          positions of the eight pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 *
 * @return void
 */

static void decode_8_pulses_31bits(int16_t *fixed_index, float *fixed_vector) {
    int pulse_position[8];
    int i, pos1, pos2, sign, temp;

    // decode pulse positions
    // coded using 7+3 bits with the 3 LSBs being, individually, the LSB of 1 of
    // the 3 pulses and the upper 7 bits being coded in base 5
    temp = fixed_index[4] >> 3;
    pulse_position[0] = (temp    %5)<<1 +  fixed_index[4]    &1;
    pulse_position[4] = (temp /5)%5)<<1 + (fixed_index[4]>>1)&1;
    pulse_position[1] = (temp/25)%5)<<1 + (fixed_index[4]>>2)&1;

    // coded using 7+3 bits with the 3 LSBs being, individually, the LSB of 1 of
    // the 3 pulses and the upper 7 bits being coded in base 5
    temp = fixed_index[5] >> 3;
    pulse_position[2] = (temp    %5)<<1 +  fixed_index[5]    &1;
    pulse_position[6] = (temp /5)%5)<<1 + (fixed_index[5]>>1)&1;
    pulse_position[5] = (temp/25)%5)<<1 + (fixed_index[5]>>2)&1;

    // coded using 5+2 bits with the 2 LSBs being, individually, the LSB of 1 of
    // the 2 pulses and the upper 5 bits being coded in base 5
    temp = ((fixed_index[6] >> 2)*25)>>5;
    pulse_position[3] = temp%5;
    pulse_position[7] = temp/5;
    if(pulse_position[7]&1)
        pulse_position[3] = 4 - pulse_position[3];
    pulse_position[3] = pulse_position[3]<<1 +  fixed_index[6]    &1;
    pulse_position[7] = pulse_position[7]<<1 + (fixed_index[6]>>1)&1;

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
 * Decode the algebraic codebook index to pulse positions and signs and construct
 * the algebraic codebook vector for MODE_122
 *
 * @param fixed_index          positions of the ten pulses
 * @param fixed_vector         pointer to the algebraic codebook vector
 *
 * @return void
 */

static void decode_10_pulses_35bits(int16_t *fixed_index, float *fixed_vector) {
    int i, pos1, pos2, sign;

    // reset the code
    memset(fixed_vector, 0, AMR_SUBFRAME_SIZE*sizeof(float));

    // the positions and signs are explicitly coded in MODE_122

    // reconstruct the fixed code
    for(i=0; i<5; i++) {
        pos1 = gray_decode[fixed_index[i  ] & 7]*5 + i; // ith pulse position
        sign = (fixed_index[i] & 8) ? -1.0 : 1.0; // sign of ith pulse
        pos2 = gray_decode[fixed_index[i+5] & 7]*5 + i; // i+5th pulse position
        // assign the ith pulse (+/-1) to its appropriate position
        fixed_vector[pos1] = sign;
        // sign of i+5th pulse is relative to sign of ith pulse
        if(pos2 < pos1) sign = -sign;
        // assign the i+5th pulse (+/-1) to its appropriate position
        fixed_vector[pos2] += sign;
    }
}

/*** end of algebraic code book (fixed) vector decoding functions ***/


/*** gain decoding functions ***/

/**
 * Predict the fixed gain
 *
 * @param fixed_vector         pointer to the algebraic codebook vector
 * @param prev_pred_error      pointer to the quantified prediction errors from the previous four subframes
 *
 * @return Returns the predicted fixed gain
 */

static float fixed_gain_prediction(float *fixed_vector, float *prev_pred_error) {
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

/*** end of gain decoding functions ***/


/*** pre-processing functions ***/

static inline float av_clipf(float a, float min, float max) {
    if(a>max)
        return max;
    else if(a<min)
        return min;
    else
        return a;
}

/**
 * Comparison function for use with qsort
 *
 * @param a             First value for comparison
 * @param b             Second value for comparison
 * @return a-b : the result of the comparison
 */

float qsort_compare(const float *a, const float *b) {
    return (float)(*a - *b);
}

/**
 * Find the median of some float values
 *
 * @param values        pointer to the values of which to find the median
 * @param n             number of values
 * @return Returns the median value
 */

static float medianf(float *values, int n) {
    float temp[9]; // largest n used for median calculation is 9

    memcpy(values, temp, n*sizeof(float));

    qsort(temp, n, sizeof(float), qsort_compare);

    if(n&1) {
        return                     temp[ n>>1 ];
    }else {
        return (temp[ (n>>1)-1 ] + temp[ n>>1 ])/2.0;
    }
}

/**
 * Circularly convolve the fixed vector with a phase dispersion impulse response
 * filter
 *
 * @param fixed_vector  pointer to the fixed vector
 * @param ir_filter     pointer to the impulse response filter
 */

static void convolve_circ(float *fixed_vector, float *ir_filter) {
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

/*** end of pre-processing functions ***/


static int amrnb_decode_frame(AVCodecContext *avctx,
        void *data, int *data_size, uint8_t *buf, int buf_size) {

    AMRContext *p = avctx->priv_data;        // pointer to private data
    int16_t *buf_out = data;                 // pointer to the output data buffer
    int i, subframe;                         // counters
    int index = 0;                           // index counter (different modes
                                             // advance through amr_prms differently)
    enum Mode speech_mode = MODE_475;        // ???

    // decode the bitstream to amr parameters
    p->cur_frame_mode = decode_bitstream(p, buf, buf_size, &speech_mode);

/*** LPC coefficient decoding ***/

    if(p->cur_frame_mode == MODE_122) {
        // decode split-matrix quantised lsf vector indices to lsp vectors
        lsf2lsp_5(p);
        // interpolate LSP vectors at subframes 1 and 3
        interp_lsp_13(p);
        // advance index into amr_prms
        index += 5;
    }else {
        // decode split-matrix quantised lsf vector indices to an lsp vector
        lsf2lsp_3(p);
        // interpolate LSP vectors at subframes 1, 2 and 3
        interp_lsp_123(p);
        // advance index into amr_prms
        index += 3;
    }

    // convert LSP vectors to LPC coefficient vectors
    for(i=0; i<4; i++) {
        lsp2lpc(p->lsp[i], p->lpc[i]);
    }

    // update averaged lsp vector (used for fixed gain smoothing)
    for(i=0; i<LP_FILTER_ORDER; i++) {
        // calculate averaged lsp vector
        p->lsp_avg[i] = 0.84*p->lsp_avg[i] + 0.16*p->prev_lsp_sub4[i];
    }

/*** end of LPC coefficient decoding ***/

    for(subframe = 0; subframe < 5; subframe++) {

/*** adaptive code book (pitch) vector decoding ***/

        // decode integer and fractional parts of pitch lag from parsed pitch
        // index
        if(p->cur_frame_mode == MODE_122) {
            decode_pitch_lag_6(p, p->amr_prms[index]);
        }else {
            decode_pitch_lag_3(p, p->amr_prms[index]);
        }

        // interpolate the past excitation at the pitch lag to obtain the pitch
        // vector
        interp_pitch_vector(p->excitation, p->pitch_lag_int, p->pitch_lag_frac, p->cur_frame_mode, p->pitch_vector);

/*** end of adaptive code book (pitch) vector decoding ***/

/*** algebraic code book (fixed) vector decoding ***/

        switch(p->cur_frame_mode) {
            case MODE_475:
            case MODE_515:
                decode_2_pulses_9bits(*index, *(index+1), subframe, &p->fixed_vector);
                index += 2;
            break;
            case MODE_59:
                decode_2_pulses_11bits(*index, *(index+1), &p->fixed_vector);
                index += 2;
            break;
            case MODE_67:
                decode_3_pulses_14bits(*index, *(index+1), &p->fixed_vector);
                index += 2;
            break;
            case MODE_74:
            case MODE_795:
                decode_4_pulses_17bits(*index, *(index+1), &p->fixed_vector);
                index += 2;
            break;
            case MODE_102:
                decode_8_pulses_31bits(index, &p->fixed_vector);
                index += 7;
            break;
            case MODE_122:
                // decode pitch gain
                p->pitch_gain[4] = qua_gain_pit[*index++];
                decode_10_pulses_35bits(index, &p->fixed_vector);
                index += 10;
            break;
            default:
            break;
        }

/*** end of algebraic code book (fixed) vector decoding ***/

/*** gain decoding ***/

        // calculate the predicted fixed gain g_c'
        p->fixed_gain[4] = fixed_gain_prediction(&p->fixed_vector, &p->prediction_error);

        // decode pitch gain and fixed gain correction factor
        if(p->cur_frame_mode == MODE_122) {
            p->fixed_gain_factor = qua_gain_code[*index++];
        }else if(p->cur_frame_mode == MODE_795) {
            p->pitch_gain[4] =     qua_gain_pit[*index++];
            p->fixed_gain_factor = qua_gain_code[*index++];
        }else if(p->cur_frame_mode == MODE_67 || p->cur_frame_mode == MODE_74 ||
                 p->cur_frame_mode == MODE_102) {
            p->pitch_gain[4] =     gains_high[index][0];
            p->fixed_gain_factor = gains_high[index][1];
            *index++;
        }else if(p->cur_frame_mode == MODE_515 || p->cur_frame_mode == MODE_59) {
            p->pitch_gain[4] =     gains_low[index][0];
            p->fixed_gain_factor = gains_low[index][1];
            *index++;
        }else {
            p->pitch_gain[4] =     gains_MODE_475[index + (subframe&1)<<1][0];
            p->fixed_gain_factor = gains_MODE_475[index + (subframe&1)<<1][1];
            *index++;
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
        if(p->pitch_gain < 0.6) {
            // strong filtering
            p->ir_filter_strength[1] = 0;
        }else if(p->pitch_gain < 0.9) {
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

    }

    /* Report how many samples we got */
    *data_size = buf_size;

    /* Return the amount of bytes consumed if everything was ok */
    return *data_size*sizeof(int16_t);
}


static int amrnb_decode_close(AVCodecContext *avctx) {
    AMRContext *p = avctx->priv_data;

    /* Free allocated memory */
    av_free(p->sample_buffer);
    av_free(p->amr_prms);

    /* Return 0 if everything is ok, -1 if not */
    return 0;
}


AVCodec amrnb_decoder =
{
    .name = "amrnb",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_AMR_NB,
    .priv_data_size = sizeof(AMRContext),
    .init = amrnb_decode_init,
    .close = amrnb_decode_close,
    .decode = amrnb_decode_frame,
};
