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

    float       prev_lsf_r[LP_FILTER_ORDER];
    float           lsp[4][LP_FILTER_ORDER]; ///< lsp vectors from current frame
    float    prev_lsp_sub4[LP_FILTER_ORDER]; ///< lsp vector for the 4th subframe of the previous frame

    float           lpc[4][LP_FILTER_ORDER]; ///< vectors of lpc coefficients for 4 subframes

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
    int idx;

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


static int amrnb_decode_frame(AVCodecContext *avctx,
        void *data, int *data_size, uint8_t *buf, int buf_size) {

    AMRContext *p = avctx->priv_data;        // pointer to private data
    int16_t *buf_out = data;                 // pointer to the output data buffer
    int i;                                   // counter
    enum Mode speech_mode = MODE_475;        // ???

    // decode the bitstream to amr parameters
    p->cur_frame_mode = decode_bitstream(p, buf, buf_size, &speech_mode);

/*** LPC coefficient decoding ***/

    if(p->cur_frame_mode == MODE_122) {
        // decode split-matrix quantised lsf vector indices to lsp vectors
        lsf2lsp_5(p);
        // interpolate LSP vectors at subframes 1 and 3
        interp_lsp_13(p);
    }else {
        // decode split-matrix quantised lsf vector indices to an lsp vector
        lsf2lsp_3(p);
        // interpolate LSP vectors at subframes 1, 2 and 3
        interp_lsp_123(p);
    }

    // convert LSP vectors to LPC coefficient vectors
    for(i=0; i<4; i++) {
        lsp2lpc(p->lsp[i], p->lpc[i]);
    }

/*** end of LPC coefficient decoding ***/

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
