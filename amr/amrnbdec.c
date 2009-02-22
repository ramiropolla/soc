/*
 * AMR narrowband decoder
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
 * @file amrnbdec.c
 * AMR narrowband decoder
 */


#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

#include "avcodec.h"
#include "bitstream.h"
#include "common.h"
#include "amrdata.h"

typedef struct AMRContext {

    GetBitContext                        gb;
    int16_t                  *sample_buffer;

    int16_t                       *amr_prms; ///< pointer to the decoded AMR parameters (lsf coefficients, codebook indexes, etc)
    int                 bad_frame_indicator; ///< bad frame ? 1 : 0

    int                   prev_frame_homing; ///< previous frame was a homing frame ? 1 : 0
    enum RXFrameType        prev_frame_type; ///< frame type of previous frame
    enum Mode               prev_frame_mode; ///< mode of previous frame

    int                        *prev_lsf_dq; ///< previous dequantized lsfs
    int                    *prev_residual_q; ///< previous quantized residual

    int             lsp1_q[LP_FILTER_ORDER]; ///< vector of quantized lsps
    int             lsp2_q[LP_FILTER_ORDER]; ///< vector of quantized lsps

    int                    cur_frame_homing; ///< current frame homing ? 1 : 0
    enum RXFrameType         cur_frame_type; ///< current frame type
    enum Mode                cur_frame_mode; ///< current frame mode

    int    lpc_coeffs[4][LP_FILTER_ORDER+1]; ///< lpc coefficients, A(z), for all four subframes
    int      prev_lsp_sub4[LP_FILTER_ORDER]; ///< vector of lsps from subframe 4 of the previous frame
    int       cur_lsp_sub2[LP_FILTER_ORDER]; ///< vector of lsps from subframe 2 of the current frame
    int       cur_lsp_sub4[LP_FILTER_ORDER]; ///< vector of lsps from subframe 4 of the current frame

    int                        cur_subframe; ///< subframe number (1-4)
    int                    search_range_min; ///< minimum of search range
    int                    search_range_max; ///< maximum of search range
    int                  prev_pitch_lag_int; ///< integer pitch lag from previous subframe (used in 2nd and 4th subframes)
    int                   cur_pitch_lag_int; ///< integer part of pitch lag from current subframe
    int                  cur_pitch_lag_frac; ///< fractional part of pitch lag from current subframe
    int                         *excitation; ///< excitation buffer

    AMRDecoderState                  *state; ///< current decoder state

} AMRContext;


static av_cold int amrnb_decode_init(AVCodecContext *avctx)
{
    AMRContext *p = avctx->priv_data;

    // variables needed for cos table generation
    int i, temp;
    double pi_step = M_PI/64.0;

    // generate the cos table
    for(i=0; i<65; i++) {
        // avoid (int16_t)32768 -> -32768
        if( (temp = lrint( 32768.0*cos(pi_step*(double)i) )) > 32767 ) {
            cos_table[i] = 32767;
        }else {
            cos_table[i] = (int16_t)temp;
        }
    }

    // allocate and zero the sample buffer
    p->sample_buffer = av_mallocz(sizeof(int16_t)*AMR_BLOCK_SIZE);
    // allocate and zero the AMR parameters
    p->amr_prms = av_mallocz(sizeof(int16_t)*PRMS_MODE_122);
    // allocate and zero the decoder state
    p->state = av_mallocz(sizeof(AMRDecoderState));

    /* Check if the allocation was successful */
    if(p->sample_buffer == NULL)
        return -1;
    // Check amr_prms allocation
    if(p->amr_prms == NULL)
        return -1;
    // Check state allocation
    if(p->state == NULL)
        return -1;

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

enum Mode decode_bitstream(AVCodecContext *avctx, uint8_t *buf, int buf_size,
                           enum Mode *speech_mode)
{
    AMRContext *p = avctx->priv_data;
    enum Mode mode;
    int i;
    const AMROrder *order;

    // initialize get_bits
    init_get_bits(&p->gb, buf, buf_size*8);
    skip_bits1(&p->gb);
    // set the mode
    mode = get_bits(&p->gb ,4);
    // set the bad frame indicator based on the quality bit
    p->bad_frame_indicator = !get_bits1(&p->gb);
    skip_bits(&p->gb, 2);

    switch(mode) {
        case MODE_DTX:
            order             = order_MODE_DTX;
            p->cur_frame_type = RX_SID_FIRST; // get SID type bit
        break;
        case NO_DATA:
            p->cur_frame_type = RX_NO_DATA;
        break;
        case MODE_475:
            order             = order_MODE_475;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_515:
            order             = order_MODE_515;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_59:
            order             = order_MODE_59;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_67:
            order             = order_MODE_67;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_74:
            order             = order_MODE_74;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_795:
            order             = order_MODE_795;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_102:
            order             = order_MODE_102;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_122:
            order             = order_MODE_122;
            p->cur_frame_type = RX_SPEECH_GOOD;
        break;
        default:
            p->cur_frame_type = RX_SPEECH_BAD;
        break;
    }

    // reorder the bitstream to match the bit allocation in the specification
    if((p->cur_frame_type != RX_NO_DATA) && (p->cur_frame_type != RX_SPEECH_BAD)) {
        for(i=0; i<mode_bits[mode]; i++) {
            p->amr_prms[ order[i].array_element ] += get_bits1(&p->gb) * (1 << order[i].bit_mask);
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


/**
 * DESCRIPTION FROM REFERENCE SOURCE:
 * Make sure that the LSFs are properly ordered and keep a certain minimum
 * distance between adjacent LSFs.
 *
 * @param lsf               a vector of lsfs (range: 0<=val<=0.5)
 * @param min_dist          minimum required separation of lsfs
 */

static void reorder_lsf(int *lsf, int min_dist)
{
    int i;
    int lsf_min = min_dist;

    for(i=0; i<LP_FILTER_ORDER; i++) {
        if(lsf[i] < lsf_min) {
            lsf[i] = lsf_min;
        }
        lsf_min = lsf[i] + min_dist;
    }
}


/**
 * Convert a vector of lsfs into the corresponding cosine domain lsps.
 *
 * @param lsf               a vector of lsfs
 * @param lsp               a vector of lsps
 */

static void lsf2lsp(int *lsf, int *lsp)
{
    int i;
    int index, offset;

    for(i=0; i<LP_FILTER_ORDER; i++) {
        index  = lsf[i] >> 8;     // bits 8 to 15 of lsf[i]
        offset = lsf[i] & 0x00ff; // bits 0 to  7 of lsf[i]
        lsp[i] = cos_table[index] + ( (( cos_table[index+1]-cos_table[index] )*offset)>>8 );
    }
}


/**
 * Convert a set of 3 lsf quantization indexes into lsp parameters.
 *
 * @param avctx             pointer to the AVCodecContext for AMR
 */

static void decode_lsf2lsp_3(AVCodecContext *avctx)
{
    AMRContext *p = avctx->priv_data;

    int lsf1_r[LP_FILTER_ORDER]; // vector of residual lsfs
    int lsf1_q[LP_FILTER_ORDER]; // vector of quantized lsfs
    const int16_t (*lsf_3_temp1)[3], (*lsf_3_temp3)[4]; // temp ptrs for switching tables depending on mode
    int index_temp; // temp lsf index
    int i; // counter

    /* if the current frame is bad, estimate the past quantized residual
     * based on the past lsf shifted slightly towards the mean */
    if(p->bad_frame_indicator) {
        for(i=0; i<LP_FILTER_ORDER; i++) {
            lsf1_q[i] = ( (p->prev_lsf_dq[i]*ALPHA)>>15 ) + ( (lsf_3_mean[i]*ONE_ALPHA)>>15 );
        }
        if(p->cur_frame_mode != MODE_DTX) {
            for(i=0; i<LP_FILTER_ORDER; i++) {
                p->prev_residual_q[i] = lsf1_q[i] - lsf_3_mean[i] - ( (p->prev_residual_q[i] * pred_fac[i])>>15 );
             }
        }else {
            for(i=0; i<LP_FILTER_ORDER; i++) {
                p->prev_residual_q[i] = lsf1_q[i] - lsf_3_mean[i] - p->prev_residual_q[i];
            }
        }
    }else {
        // assign lsf tables according to mode
        if((p->cur_frame_mode == MODE_475) || (p->cur_frame_mode == MODE_515)) {
            lsf_3_temp1 = lsf_3_1;
            lsf_3_temp3 = lsf_3_MODE_515;
        }else if(p->cur_frame_mode == MODE_795) {
            lsf_3_temp1 = lsf_3_MODE_795;
            lsf_3_temp3 = lsf_3_3;
        }else {
            lsf_3_temp1 = lsf_3_1;
            lsf_3_temp3 = lsf_3_3;
        }

        // decode lsf residuals from lsf tables using the bitstream indexes
        lsf1_r[0] = lsf_3_temp1[ p->amr_prms[0] ][0];
        lsf1_r[1] = lsf_3_temp1[ p->amr_prms[0] ][1];
        lsf1_r[2] = lsf_3_temp1[ p->amr_prms[0] ][2];

        index_temp = p->amr_prms[1];
        /* MODE_475, MODE_515 only use every second entry
         * NOTE: unsure what they mean and the spec doesn't mention this */
        if((p->cur_frame_mode == MODE_475) || (p->cur_frame_mode == MODE_515)) {
            index_temp = index_temp << 1;
        }
        lsf1_r[3] = lsf_3_2[index_temp][0];
        lsf1_r[4] = lsf_3_2[index_temp][1];
        lsf1_r[5] = lsf_3_2[index_temp][2];

        lsf1_r[6] = lsf_3_temp3[ p->amr_prms[2] ][0];
        lsf1_r[7] = lsf_3_temp3[ p->amr_prms[2] ][1];
        lsf1_r[8] = lsf_3_temp3[ p->amr_prms[2] ][2];
        lsf1_r[9] = lsf_3_temp3[ p->amr_prms[2] ][3];

        /* Compute quantized LSFs and update the past quantized residual */
        if(p->cur_frame_mode != MODE_DTX) {
            for(i=0; i<LP_FILTER_ORDER; i++) {
                lsf1_q[i] = lsf1_r[i] + lsf_3_mean[i] + ( (p->prev_residual_q[i]*pred_fac[i])>>15 );
            }
            memcpy( p->prev_residual_q, lsf1_r, LP_FILTER_ORDER <<2 );
        }else {
            for(i=0; i<LP_FILTER_ORDER; i++) {
                lsf1_q[i] = lsf1_r[i] + lsf_3_mean[i] + p->prev_residual_q[i];
            }
            memcpy(p->prev_residual_q, lsf1_r, LP_FILTER_ORDER<<2);
        }
    }

    /* verification that LSFs have minimum distance of LSF_GAP Hz */
    reorder_lsf(lsf1_q, LSF_GAP);
    memcpy(p->prev_lsf_dq, lsf1_q, LP_FILTER_ORDER<<2);

    /*  convert LSFs to the cosine domain */
    lsf2lsp(lsf1_q, p->lsp1_q);
}


/**
 * Convert a set of 5 lsf quantization indexes into lsp parameters.
 *
 * @param avctx             pointer to the AVCodecContext for AMR
 */

static void decode_lsf2lsp_5(AVCodecContext *avctx)
{
    AMRContext *p = avctx->priv_data;

    int lsf1_r[LP_FILTER_ORDER], lsf2_r[LP_FILTER_ORDER]; // vectors of residual lsfs
    int lsf1_q[LP_FILTER_ORDER], lsf2_q[LP_FILTER_ORDER]; // vectors of quantized lsfs
    int temp;
    int i, sign; // counter and sign of 3rd lsf table

    /* if the current frame is bad, estimate the past quantized residual
     * based on the past lsf shifted slightly towards the mean */
    if(p->bad_frame_indicator) {
        for(i=0; i<LP_FILTER_ORDER; i++) {
            lsf1_q[i] = ( (p->prev_lsf_dq[i]*ALPHA_122)>>15 ) + ( (lsf_5_mean[i]*ONE_ALPHA_122)>>15 );
        }
        memcpy(lsf2_q, lsf1_q, LP_FILTER_ORDER<<1);
        for(i=0; i<LP_FILTER_ORDER; i++) {
            p->prev_residual_q[i] = lsf2_q[i] - lsf_5_mean[i] - ( (p->prev_residual_q[i]*LSP_PRED_FAC_MODE_122)>>15 );
        }
    }else {
        /* decode prediction residuals from 5 received indexes */
        lsf1_r[0] = lsf_5_1[ p->amr_prms[0] ][0];
        lsf1_r[1] = lsf_5_1[ p->amr_prms[0] ][1];
        lsf2_r[0] = lsf_5_1[ p->amr_prms[0] ][2];
        lsf2_r[1] = lsf_5_1[ p->amr_prms[0] ][3];

        lsf1_r[2] = lsf_5_2[ p->amr_prms[1] ][0];
        lsf1_r[3] = lsf_5_2[ p->amr_prms[1] ][1];
        lsf2_r[2] = lsf_5_2[ p->amr_prms[1] ][2];
        lsf2_r[3] = lsf_5_2[ p->amr_prms[1] ][3];

        sign = (p->amr_prms[2] & 1) ? -1 : 1;
        // lsb of p->amr_prms[2] is the sign bit so p->amr_prms[2]>>1
        lsf1_r[4] = lsf_5_3[ p->amr_prms[2]>>1 ][0]*sign;
        lsf1_r[5] = lsf_5_3[ p->amr_prms[2]>>1 ][1]*sign;
        lsf2_r[4] = lsf_5_3[ p->amr_prms[2]>>1 ][2]*sign;
        lsf2_r[5] = lsf_5_3[ p->amr_prms[2]>>1 ][3]*sign;

        lsf1_r[6] = lsf_5_4[ p->amr_prms[3] ][0];
        lsf1_r[7] = lsf_5_4[ p->amr_prms[3] ][1];
        lsf2_r[6] = lsf_5_4[ p->amr_prms[3] ][2];
        lsf2_r[7] = lsf_5_4[ p->amr_prms[3] ][3];

        lsf1_r[8] = lsf_5_5[ p->amr_prms[4] ][0];
        lsf1_r[9] = lsf_5_5[ p->amr_prms[4] ][1];
        lsf2_r[8] = lsf_5_5[ p->amr_prms[4] ][2];
        lsf2_r[9] = lsf_5_5[ p->amr_prms[4] ][3];

        /* Compute quantized LSFs and update the past quantized residual */
        for(i=0; i<LP_FILTER_ORDER; i++) {
            temp = lsf_5_mean[i] + ( (p->prev_residual_q[i]*LSP_PRED_FAC_MODE_122)>>15 );
            lsf1_q[i] = lsf1_r[i] + temp;
            lsf2_q[i] = lsf2_r[i] + temp;
            p->prev_residual_q[i] = lsf2_r[i];
        }
    }

    /* verification that LSFs have minimum distance of LSF_GAP Hz */
    reorder_lsf(lsf1_q, LSF_GAP);
    reorder_lsf(lsf2_q, LSF_GAP);
    memcpy(p->prev_lsf_dq, lsf2_q, LP_FILTER_ORDER<<2);

    /*  convert LSFs to the cosine domain */
    lsf2lsp(lsf1_q, p->lsp1_q);
    lsf2lsp(lsf2_q, p->lsp2_q);
}


/**
 * Find the polynomial F1(z) or F2(z) from the lsps
 *
 * The ref source description of how to find the polynomials was as follows:
 *    Find the polynomial F1(z) or F2(z) from the LSPs.
 *
 *    F1(z) = product ( 1 - 2 * lsp[i] * z^-1 + z^-2 )
 *             i=0,2,4,6,8
 *    F2(z) = product ( 1 - 2 * lsp[i] * z^-1 + z^-2 )
 *             i=1,3,5,7,9
 *
 * @param lsp               vector of lsps
 * @param f                 pointer to the polynomial F1(z) or F2(z)
 */

static void lsp2poly(int *lsp, int *f)
{
    int i, j;

    f[0] = 1<<24;
    f[1] = -(lsp[0]<<10);
    for(i=2; i<6; i++) {
        f[i] = (f[i-2]<<1) - (( (f[i-1]>>16)*lsp[2*i-2] + ( ((f[i-1] & 0xFFFE)*lsp[2*i-2]) >>16) )<<2);
        for(j=i-1; j>1; j--) {
            f[j] += f[j-2] - (( (f[j-1]>>16)*lsp[2*i-2] + ( ((f[j-1] & 0xFFFE)*lsp[2*i-2]) >>16) )<<2);
        }
        f[1] -= lsp[2*i-2]<<10;
    }
}
/* FIXME - check rounding sensitivity of F1(z) and F2(z)
 * as the code below is faster but rounds differently
static void lsp2poly(int *lsp, int *f) {
    int i, j;
    f[0] = 1<<24;
    f[1] = -lsp[0]<<10;
    for(i=2; i<6; i++) {
        f[i] = f[i-2];
        for(j=i; j>1; j--)
            f[j] +=  f[j-2] - ((lsp[2*i-2] * (int64_t)f[j-1] + (1<<13))>>14);
        f[1] -= lsp[2*i-2]<<10;
    }
}
*/


/**
 * Convert a vector of lsps to lpc coefficients for a 10th order filter.
 *
 * @param lsp                 vector of lsps
 * @param lpc_coeffs          pointer to a subframe of lpc coefficients
 */

static void lsp2lpc(int *lsp, int *lpc_coeffs)
{
    int f1[6], f2[6];
    int temp, i;

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
    // lpc_coeffs[i] = 0.5*f1[i]    + 0.5*f2[i]       for i = 1,...,5
    //                 0.5*f1[11-i] - 0.5*f2[11-i]    for i = 6,...,10
    lpc_coeffs[0] = 0x1000;
    for(i=1; i<6; i++) {
        temp = f1[i] + f2[i];
        lpc_coeffs[i] = (int16_t)(temp>>13);  // ref source: emulate fixed point bug
        if(temp & 0x1000) {
            lpc_coeffs[i]++;
        }
        temp = f1[11-i] - f2[11-i];
        lpc_coeffs[11-i] = (int16_t)(temp>>13);  // ref source: emulate fixed point bug
        if(temp & 0x1000) {
            lpc_coeffs[11-i]++;
        }
    }
}


/**
 * Interpolate lsps for subframes 1 and 3 and convert them to lpc coefficients.
 *
 * @param avctx                 pointer to AVCodecContext
 * @param lpc_coeffs            array of lpc coefficients for the four subframes
 */

static void lpc_interp_13(AVCodecContext *avctx, int **lpc_coeffs)
{
    AMRContext *p = avctx->priv_data;
    int lsp[LP_FILTER_ORDER];
    int i;

    // subframe 1: interpolate lsp vector
    // q1(n) = 0.5q4(n-1) + 0.5q2(n)
    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsp[i] = (p->cur_lsp_sub2[i]>>1) + (p->prev_lsp_sub4[i]>>1);
    }
    // subframe 1: convert the lsps to lpc coefficients
    lsp2lpc(lsp, lpc_coeffs[0]);

    // subframe 2: convert the lsps to lpc coefficients
    lsp2lpc(p->cur_lsp_sub2, lpc_coeffs[1]);

    // subframe 3: interpolate lsp vector
    // q3(n) = 0.5q2(n) + 0.5q4(n)
    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsp[i] = (p->cur_lsp_sub2[i]>>1) + (p->cur_lsp_sub4[i]>>1);
    }
    // subframe 3: convert the lsps to lpc coefficients
    lsp2lpc(lsp, lpc_coeffs[2]);

    // subframe 4: convert the lsps to lpc coefficients
    lsp2lpc(lsp, lpc_coeffs[3]);
}


/**
 * Interpolate lsps for subframes 1, 2 and 3 and convert them to lpc coefficients.
 *
 * @param avctx                 pointer to AVCodecContext
 * @param lpc_coeffs            array of lpc coefficients for the four subframes
 */

static void lpc_interp_123(AVCodecContext *avctx, int **lpc_coeffs)
{
    AMRContext *p = avctx->priv_data;
    int lsp[LP_FILTER_ORDER];
    int i;

    // subframe 1: interpolate lsp vector
    // q1(n) = 0.75q4(n-1) + 0.25q4(n)
    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsp[i] = ((p->cur_lsp_sub4[i] - p->prev_lsp_sub4[i])>>2) + p->prev_lsp_sub4[i];
    }
    // subframe 1: convert the lsps to lpc coefficients
    lsp2lpc(lsp, lpc_coeffs[0]);

    // subframe 2: interpolate lsp vector
    // q2(n) = 0.5q4(n-1) + 0.5q4(n)
    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsp[i] = (p->prev_lsp_sub4[i]>>1) + (p->cur_lsp_sub4[i]>>1);
    }
    // subframe 2: convert the lsps to lpc coefficients
    lsp2lpc(lsp, lpc_coeffs[1]);

    // subframe 3: interpolate lsp vector
    // q3(n) = 0.25q4(n-1) + 0.75q4(n)
    for(i=0; i<LP_FILTER_ORDER; i++) {
        lsp[i] = (p->prev_lsp_sub4[i]>>2) + ( p->cur_lsp_sub4[i] - (p->cur_lsp_sub4[i]>>2) );
    }
    // subframe 3: convert the lsps to lpc coefficients
    lsp2lpc(lsp, lpc_coeffs[2]);

    // subframe 4: convert the lsps to lpc coefficients
    lsp2lpc(p->cur_lsp_sub4, lpc_coeffs[3]);
}


/**
 * Decode the adaptive codebook index to the integer and fractional
 * pitch lags for one subframe at 1/3 resolution.
 *
 * Description from reference source:
 *    The fractional lag in 1st and 3rd subframes is encoded with 8 bits
 *    while that in 2nd and 4th subframes is relatively encoded with 4, 5
 *    and 6 bits depending on the mode.
 *
 * @param avctx              pointer to AVCodecContext
 * @param pitch_index        received adaptive codebook (pitch) index
 * @param pitch_lag_int      pointer to integer part of the subframe pitch lag
 * @param pitch_lag_frac     pointer to fractional part of the subframe pitch lag
 */

static void decode_pitch_lag_3(AVCodecContext *avctx, int pitch_index,
                               int *pitch_lag_int, int *pitch_lag_frac)
{
    AMRContext *p = avctx->priv_data;
    int tmp_lag;

    // FIXME: find out where these constants come from and add appropriate
    // comments such that people reading the code can understand why these
    // particular values are used

    // subframe 1 or 3
    if(p->cur_subframe & 1) {
        if(pitch_index < 197) {
            *pitch_lag_int  = (( (pitch_index + 2)*10923 )>>15) + 19;
            *pitch_lag_frac = pitch_index - *pitch_lag_int*3 + 58;
            // For the sake of deciphering where this comes from,
            // the above code probably means:
            // *pitch_lag_int  = (pitch_index + 2)/3 + 19;
            // *pitch_lag_frac = (pitch_index + 2)%3 + 56;
        }else {
            *pitch_lag_int  = pitch_index - 112;
            *pitch_lag_frac = 0;
        }
    // subframe 2 or 4
    }else {
        if( (p->cur_frame_mode == MODE_475) || (p->cur_frame_mode == MODE_515) ||
            (p->cur_frame_mode == MODE_59)  || (p->cur_frame_mode == MODE_67) ) {

            // decoding with 4 bit resolution
            tmp_lag= clip(p->prev_pitch_lag_int, p->search_range_max-4, p->search_range_min+5);

            if(pitch_index < 4) {
                *pitch_lag_int  = pitch_index + tmp_lag - 5;
                *pitch_lag_frac = 0;
            }else {
                if(pitch_index < 12) {
                    *pitch_lag_int  = ( ((pitch_index - 5)*10923)>>15 ) + tmp_lag - 1;
                    *pitch_lag_frac = pitch_index - *pitch_lag_int*3 - 9;
                }else {
                    *pitch_lag_int  = pitch_index + tmp_lag - 11;
                    *pitch_lag_frac = 0;
                }
            }
        }else {
            // decoding with 5 or 6 bit resolution
            *pitch_lag_int  = ( ((pitch_index + 2)*10923 )>>15) - 1 + p->search_range_min;
            *pitch_lag_frac = pitch_index - (*pitch_lag_int - p->search_range_min)*3 - 2;
        }
    }
}


/**
 * Decode the adaptive codebook index to the integer and fractional
 * pitch lag for one subframe at 1/6 resolution.
 *
 * Description from reference source:
 *    The fractional lag in 1st and 3rd subframes is encoded with 9 bits
 *    while that in 2nd and 4th subframes is relatively encoded with 6 bits.
 *    Note that in relative encoding only 61 values are used. If the
 *    decoder receives 61, 62, or 63 as the relative pitch index, it means
 *    that a transmission error occurred. In this case, the pitch lag from
 *    the previous subframe (actually from the previous frame) is used.
 *
 * @param avctx              pointer to AVCodecContext
 * @param pitch_index        received adaptive codebook (pitch) index
 * @param pitch_lag_int      pointer to integer part of the subframe pitch lag
 * @param pitch_lag_frac     pointer to fractional part of the subframe pitch lag
 */

static void decode_pitch_lag_6(AVCodecContext *avctx, int pitch_index,
                               int *pitch_lag_int, int *pitch_lag_frac)
{
    AMRContext *p = avctx->priv_data;
    int temp;

    // FIXME - find out where these constants come from and add appropriate
    // comments such that people reading the code can understand why these
    // particular values are used.

    // subframe 1 or 3
    if(p->cur_subframe & 1) {
        if(pitch_index < 463){
            *pitch_lag_int  = (pitch_index + 5) / 6 + 17;
            *pitch_lag_frac = pitch_index - *pitch_lag_int*6 + 105;
        }else {
            *pitch_lag_int  = pitch_index - 368;
            *pitch_lag_frac = 0;
        }
    // subframe 2 or 4
    }else {
        // find the search range
        // FIXME - set pitch_lag_int = 0 before each frame
        p->search_range_min = FFMAX(*pitch_lag_int - 5, PITCH_LAG_MIN_MODE_122);
        p->search_range_max = p->search_range_min + 9;
        if(p->search_range_max > PITCH_LAG_MAX) {
            p->search_range_max = PITCH_LAG_MAX;
            p->search_range_min = p->search_range_max - 9;
        }
        // calculate the pitch lag
        temp = (pitch_index + 5)/6 - 1;
        *pitch_lag_int  = temp + p->search_range_min;
        *pitch_lag_frac = (pitch_index - 3) - temp*6;
    }
}


/**
 * Find the adaptive codebook (pitch) vector from interpolating the past
 * excitation at the decoded pitch lag with corresponding resolution
 * of 1/3 or 1/6.
 *
 * See the specification for details of the underlying mathematics and the used
 * FIR filter.
 *
 * @param avctx              pointer to AVCodecContext
 * @param excitation         pointer to the excitation buffer
 */

static void decode_pitch_vector(AVCodecContext *avctx, int *excitation)
{
    AMRContext *p = avctx->priv_data;
    int i, j, temp;
    int *excitation_temp;

    excitation_temp = &excitation[-p->cur_pitch_lag_int];
    p->cur_pitch_lag_frac *= -1;
    if(p->cur_frame_mode != MODE_122) {
        p->cur_pitch_lag_frac <<= 1;
    }
    if(p->cur_pitch_lag_frac < 0) {
        p->cur_pitch_lag_frac += 6;
        excitation_temp--;
    }
    for(i=0; i<AMR_SUBFRAME_SIZE; i++) {
        // reset temp
        temp = 0;
        for(j=0; j<10; j++) {
            temp += excitation_temp[i-j  ] * b60[ j   *6 + p->cur_pitch_lag_frac];
            temp += excitation_temp[i+j+1] * b60[(j+1)*6 - p->cur_pitch_lag_frac];
        }
        excitation[i] = (temp + 0x4000)>>15;
    }
}


/**
 * Reconstruct the algebraic codebook vector.
 *
 * @param fixed_code           algebraic codebook vector
 * @param pulse_position       vector of pulse positions
 * @param sign                 signs of the pulses
 * @param nr_pulses            number of pulses
 */

static void reconstruct_fixed_code(int *fixed_code, int *pulse_position,
                                   int sign, int nr_pulses)
{
    int i;

    // reset the code
    memset(fixed_code, 0, AMR_SUBFRAME_SIZE*sizeof(int));

    // assign the pulse values (+/-1) to their appropriate positions
    for(i=0; i<nr_pulses; i++)
        fixed_code[pulse_position[i]] = ((sign >> i) & 1) ? 8191 : -8192;
}


/**
 * Decode the algebraic codebook index to pulse position indexes for MODE_102.
 *
 * @param fixed_index          positions of the eight pulses (encoded)
 * @param position_index       position index of the eight pulses
 */

static void fixed2position(int16_t *fixed_index, int *position_index)
{
    int MSBs, LSBs, MSBs0_24, divMSB;

    // indexes from track 1 (7+3 bits)
    MSBs= FFMIN(fixed_index[0] >> 3, 124);
    // LSBs = fixed_index[0]%8;
    LSBs = fixed_index[0] & 0x7;
    // position_index[0] = ((MSBs-25*(MSBs/25))%5)*2 + (LSBs-4*(LSBs/4))%2;
    // position_index[4] = ((MSBs-25*(MSBs/25))/5)*2 + (LSBs-4*(LSBs/4))/2;
    // position_index[1] = (MSBs/25)*2 + LSBs/4;
    divMSB = MSBs/25;
    position_index[0] = (( (MSBs - 25*divMSB)%5 )<<1) +  (LSBs & 0x1);
    position_index[4] = (( (MSBs - 25*divMSB)/5 )<<1) + ((LSBs & 0x2)>>1);
    position_index[1] = (divMSB<<1) + (LSBs>>2);

    // indexes from track 2 (7+3 bits)
    MSBs= FFMIN(fixed_index[1] >> 3, 124);
    // LSBs = fixed_index[1]%8;
    LSBs = fixed_index[1] & 0x7;
    // position_index[2] = ((MSBs-25*(MSBs/25))%5)*2 + (LSBs-4*(LSBs/4))%2;
    // position_index[6] = ((MSBs-25*(MSBs/25))/5)*2 + (LSBs-4*(LSBs/4))/2;
    // position_index[5] = (MSBs/25)*2 + LSBs/4;
    divMSB = MSBs/25;
    position_index[2] = (( (MSBs - 25*divMSB)%5 )<<1) +  (LSBs & 0x1);
    position_index[6] = (( (MSBs - 25*divMSB)/5 )<<1) + ((LSBs & 0x2)>>1);
    position_index[5] = (divMSB<<1) + (LSBs>>2);

    // indexes from track 3 (5+2 bits)
    MSBs = fixed_index[2] >> 2;
    // LSBs = fixed_index[2]%4;
    LSBs = fixed_index[2] & 0x3;
    MSBs0_24 = (( MSBs*25 + 12 )>>5);
    // FIXME - check equivalence of the following with the ref source as their code
    // didn't appear to match with their description at first glance on which the
    // code below is based
    if((MSBs0_24/5)%2==1)
        position_index[3] = ( (4-(MSBs0_24%5))<<1 ) + ( LSBs & 0x1 );
    else
        position_index[3] = (    (MSBs0_24%5) <<1 ) + ( LSBs & 0x1 );
    position_index[7] = ((MSBs0_24 / 5)<<1) + ( LSBs >> 1 );
}


/**
 * Decode the algebraic codebook index to pulse positions and signs and
 * construct the algebraic codebook vector for MODE_475 and MODE_515.
 *
 * @param avctx                pointer to AVCodecContext
 * @param sign                 signs of the two pulses
 * @param fixed_index          positions of the two pulses
 * @param fixed_code           pointer to the algebraic codebook vector
 */

static void decode_2_pulses_9bits(AVCodecContext *avctx, int sign,
                                  int fixed_index, int *fixed_code)
{
    AMRContext *p = avctx->priv_data;
    int pulse_position[2];
    int pulse_subset;

    // find the subset of pulses used
    pulse_subset = (fixed_index & 0x40)>>6;
    // find the position of the first pulse
    pulse_position[0] = ( fixed_index       & 7)*5 + track_position[ (pulse_subset<<3) + (p->cur_subframe<<1) ];
    // find the position of the second pulse
    pulse_position[1] = ((fixed_index >> 3) & 7)*5 + track_position[ (pulse_subset<<3) + (p->cur_subframe<<1) + 1 ];

    // reconstruct the fixed code
    reconstruct_fixed_code(fixed_code, pulse_position, sign, 2);
}


/**
 * Decode the algebraic codebook index to pulse positions and signs
 * and construct the algebraic codebook vector for MODE_59.
 *
 * @param sign                 signs of the two pulses
 * @param fixed_index          positions of the two pulses
 * @param fixed_code           pointer to the algebraic codebook vector
 */

static void decode_2_pulses_11bits(int sign, int fixed_index, int *fixed_code)
{
    int pulse_position[2];
    int pulse_subset;

    // find the subset of pulses used for the first pulse
    pulse_subset = fixed_index & 1;
    // find the position of the first pulse
    pulse_position[0] =     ((fixed_index >> 1) & 7)*5 + pulse_subset*2 + 1;
    // find the subset of pulses used for the second pulse
    pulse_subset = (fixed_index >> 4) & 3;
    // find the position of the second pulse
    if(pulse_subset == 3) {
        pulse_position[1] = ((fixed_index >> 6) & 7)*5 + 4;
    }else {
        pulse_position[1] = ((fixed_index >> 6) & 7)*5 + pulse_subset;
    }

    // reconstruct the fixed code
    reconstruct_fixed_code(fixed_code, pulse_position, sign, 2);
}


/**
 * Decode the algebraic codebook index to pulse positions and signs
 * and construct the algebraic codebook vector for MODE_67.
 *
 * @param sign                 signs of the three pulses
 * @param fixed_index          positions of the three pulses
 * @param fixed_code           pointer to the algebraic codebook vector
 */

static void decode_3_pulses_14bits(int sign, int fixed_index, int *fixed_code)
{
    int pulse_position[3];
    int pulse_subset;

    // find the position of the first pulse
    pulse_position[0] = ( fixed_index       & 7)*5;
    // find the subset of pulses used for the second pulse
    pulse_subset = (fixed_index >> 3) & 1;
    // find the position of the second pulse
    pulse_position[1] = ((fixed_index >> 4) & 7)*5 + pulse_subset*2 + 1;
    // find the subset of pulses used for the third pulse
    pulse_subset = (fixed_index >> 7) & 1;
    // find the position of the third pulse
    pulse_position[2] = ((fixed_index >> 8) & 7)*5 + pulse_subset*2 + 2;

    // reconstruct the fixed code
    reconstruct_fixed_code(fixed_code, pulse_position, sign, 3);
}


/**
 * Decode the algebraic codebook index to pulse positions and signs and
 * construct the algebraic codebook vector for MODE_74 and MODE_795.
 *
 * @param sign                 signs of the four pulses
 * @param fixed_index          positions of the four pulses
 * @param fixed_code           pointer to the algebraic codebook vector
 */

static void decode_4_pulses_17bits(int sign, int fixed_index, int *fixed_code)
{
    int pulse_position[4];
    int pulse_subset;

    // find the position of the first pulse
    pulse_position[0] = dgray[ fixed_index        & 7]*5;
    // find the position of the second pulse
    pulse_position[1] = dgray[(fixed_index >> 3)  & 7]*5 + 1;
    // find the position of the third pulse
    pulse_position[2] = dgray[(fixed_index >> 6)  & 7]*5 + 2;
    // find the subset of pulses used for the fourth pulse
    pulse_subset = (fixed_index >> 9) & 1;
    // find the position of the fourth pulse
    pulse_position[3] = dgray[(fixed_index >> 10) & 7]*5 + pulse_subset + 3;

    // reconstruct the fixed code
    reconstruct_fixed_code(fixed_code, pulse_position, sign, 4);
}


/**
 * Decode the algebraic codebook index to pulse positions and signs
 * and construct the algebraic codebook vector for MODE_102.
 *
 * @param fixed_index          positions of the eight pulses
 * @param fixed_code           pointer to the algebraic codebook vector
 */

static void decode_8_pulses_31bits(int16_t *fixed_index, int *fixed_code)
{
    int position_index[8];
    int i, pos1, pos2, sign;

    // reset the code
    memset(fixed_code, 0, AMR_SUBFRAME_SIZE*sizeof(int));
    fixed2position(&fixed_index[TRACKS_MODE_102], position_index);

    for(i=0; i<TRACKS_MODE_102; i++) {
        // find the position of the ith pulse
        pos1 = (position_index[i]   << 2) + i;
        // find the position of the i+4th pulse
        pos2 = (position_index[i+4] << 2) + i;
        // find the sign of the ith pulse
        sign = !fixed_index[i] ? 8191 : -8191; // +/-1 : -8191 is used here in the ref source
        // assign the ith pulse (+/-1) to its appropriate position
        fixed_code[pos1] = sign;
        // find the sign of the i+4th pulse (relative to the sign of the ith pulse)
        if(pos2 < pos1) sign = -sign;
        // assign the i+4th pulse (+/-1) to its appropriate position
        fixed_code[pos2] += sign;
    }
}


/**
 * Decode the algebraic codebook index to pulse positions and signs
 * and construct the algebraic codebook vector for MODE_122.
 *
 * @param fixed_index          positions of the ten pulses
 * @param fixed_code           pointer to the algebraic codebook vector
 */

static void decode_10_pulses_35bits(int16_t *fixed_index, int *fixed_code)
{
    int i, pos1, pos2, sign;

    // reset the code
    memset(fixed_code, 0, AMR_SUBFRAME_SIZE*sizeof(int));

    for(i=0; i<5; i++) {
        // find the position of the ith pulse
        pos1 = dgray[fixed_index[i  ] & 7]*5 + i;
        // find the sign of the ith pulse
        sign = (fixed_index[i] & 8) ? -4096 : 4096; // +/-1 : 4096 is used here in the reference source
        // find the position of the i+5th pulse
        pos2 = dgray[fixed_index[i+5] & 7]*5 + i;
        // assign the ith pulse (+/-1) to its appropriate position
        fixed_code[pos1] = sign;
        // find the sign of the i+5th pulse (relative to the sign of the ith pulse)
        if(pos2 < pos1) sign = -sign;
        // assign the i+5th pulse (+/-1) to its appropriate position
        fixed_code[pos2] += sign;
    }
}


// general functions FIXME - useful enough to put into libavutil?

/**
 * Comparison function for use with qsort.
 *
 * @param a             first value for comparison
 * @param b             second value for comparison
 * @return a-b : the result of the comparison
 */

int qsort_compare(const int *a, const int *b)
{
    return (int)(*a - *b);
}

/**
 * Find the median of some values.
 *
 * @param values        pointer to the values of which to find the median
 * @param n             number of values
 * @return the median value
 */

static int median(int *values, int n)
{
    int temp[9]; // largest n used for median calculation is 9

    memcpy(values, temp, n*sizeof(int));

    qsort(temp, n, sizeof(int), qsort_compare);

    return(temp[ n>>1 ]);
}


// gain functions

/**
 * Calculate the pitch gain from previous values.
 *
 * @param state_ptr             pointer to the current state
 * @return the pitch gain
 */

static int find_pitch_gain(AMRDecoderState *state_ptr)
{
    int temp_median;

    // find the median of the previous five pitch gains
    temp_median = median(state_ptr->prev_pitch_gains, 5);

    // clip the median pitch gain to the previous pitch gain
    if(temp_median > state_ptr->prev_pitch_gain) {
        temp_median = state_ptr->prev_pitch_gain;
    }
    return ( (temp_median*pitch_gain_attenuation[ state_ptr->state ])>>15 );
}


/**
 * Decode the pitch gain using the received index.
 *
 * @param mode              current mode
 * @param index             quantization index
 * @return the pitch gain
 */

static int decode_pitch_gain(enum Mode mode, int index)
{
    int gain;

    if(mode == MODE_122) {
        // zero the two least significant bits
        // gain = ( pitch_gain_quant[index]>>2 )<<2;
        gain = pitch_gain_quant[index] & 0xFFFC;
    }else {
        gain = pitch_gain_quant[index];
    }
    return gain;
}


/**
 * Update the pitch gain and limit pitch_gain if the previous frame was bad.
 *
 * @param state_ptr             pointer to the current state
 * @param bad_frame_indicator   bad frame indicator
 * @param pitch_gain            pointer to the pitch gain
 */

static void pitch_gain_update(AMRDecoderState *state_ptr,
                              int bad_frame_indicator, int *pitch_gain)
{
    if(bad_frame_indicator == 0) {
        if(state_ptr->prev_frame_bad != 0) {
            // if the previous frame was bad, limit the current pitch gain to
            // the previous good pitch gain
            if(*pitch_gain > state_ptr->prev_good_pitch_gain) {
                *pitch_gain = state_ptr->prev_good_pitch_gain;
            }
        }
        // if the current frame is good, update the previous good pitch gain
        state_ptr->prev_good_pitch_gain = *pitch_gain;
    }
    state_ptr->prev_pitch_gain = *pitch_gain;

    // clip the previous pitch gain to 1.0
    if(state_ptr->prev_pitch_gain > 16384) {
        state_ptr->prev_pitch_gain = 16384;
    }

    // update the array of the previous five pitch gains
    state_ptr->prev_pitch_gains[0] = state_ptr->prev_pitch_gains[1];
    state_ptr->prev_pitch_gains[1] = state_ptr->prev_pitch_gains[2];
    state_ptr->prev_pitch_gains[2] = state_ptr->prev_pitch_gains[3];
    state_ptr->prev_pitch_gains[3] = state_ptr->prev_pitch_gains[4];
    state_ptr->prev_pitch_gains[4] = state_ptr->prev_pitch_gain;
}


/**
 * Reset the AMR frame parameters.
 *
 * @param avctx             pointer to the AVCodecContext for AMR
 */

void decode_reset(AVCodecContext *avctx)
{
    AMRContext *p = avctx->priv_data;

    p->prev_frame_homing = 1;
    p->prev_frame_type   = RX_SPEECH_GOOD;
    p->prev_frame_mode   = MODE_475;
    // FIXME reset AMRDecoderState too!
}


static int amrnb_decode_frame(AVCodecContext *avctx, void *data,
                              int *data_size, uint8_t *buf, int buf_size)
{
    AMRContext *p = avctx->priv_data;        // pointer to private data
    int16_t *outbuffer = data;               // pointer to the output data buffer
    int i;                                   // counter
    enum Mode speech_mode = MODE_475;        // ???
    const int16_t *homing_frame;             // pointer to the homing frame
    int homing_frame_size;                   // homing frame size

    // decode the bitstream to AMR parameters
    p->cur_frame_mode = decode_bitstream(avctx, buf, buf_size, &speech_mode);

    // guess the mode from the previous frame if no data or bad data
    if(p->cur_frame_type == RX_SPEECH_BAD) {
        if(p->prev_frame_type > RX_SPEECH_BAD) {
            p->cur_frame_type = RX_SID_BAD;
            p->cur_frame_mode = MODE_DTX;
        }else {
            p->cur_frame_mode = p->prev_frame_mode;
        }
    }else if(p->cur_frame_type == RX_NO_DATA) {
        p->cur_frame_mode = p->prev_frame_mode;
    }

    if(p->bad_frame_indicator) {
        if(p->cur_frame_mode < MODE_DTX) {
            p->cur_frame_type = RX_SPEECH_BAD;
        }else if(p->cur_frame_mode != NO_DATA) {
            p->cur_frame_type = RX_SID_BAD;
        }
    }

    if(p->prev_frame_homing) {
        switch(p->cur_frame_mode) {
            case MODE_475:
                homing_frame      = dhf_MODE_475;
                homing_frame_size = 7;
            break;
            case MODE_515:
                homing_frame      = dhf_MODE_515;
                homing_frame_size = 7;
            break;
            case MODE_59:
                homing_frame      = dhf_MODE_59;
                homing_frame_size = 7;
            break;
            case MODE_67:
                homing_frame      = dhf_MODE_67;
                homing_frame_size = 7;
            break;
            case MODE_74:
                homing_frame      = dhf_MODE_74;
                homing_frame_size = 7;
            break;
            case MODE_795:
                homing_frame      = dhf_MODE_795;
                homing_frame_size = 8;
            break;
            case MODE_102:
                homing_frame      = dhf_MODE_102;
                homing_frame_size = 12;
            break;
            case MODE_122:
                homing_frame      = dhf_MODE_122;
                homing_frame_size = 18;
            break;
            default:
                homing_frame      = NULL;
                homing_frame_size = 0;
            break;
        }
        for(i=0; i<homing_frame_size; i++) {
            // check if the frame is homing
            if( (p->cur_frame_homing = p->amr_prms[i] ^ homing_frame[i]) ) break;
        }
    }

    if(!p->cur_frame_homing && p->prev_frame_homing) {
        for(i=0; i<AMR_BLOCK_SIZE; i++) {
            p->sample_buffer[i] = EHF_MASK;
        }
    }else {
        // decode frame (ref funcn): Speech_Decode_Frame( p, p->cur_frame_mode, p->amr_prms, p->cur_frame_type, p->sample_buffer )
    }

    if(!p->prev_frame_homing) {
        switch(p->cur_frame_mode) {
            case MODE_475:
                homing_frame      = dhf_MODE_475;
                homing_frame_size = PRMS_MODE_475;
            break;
            case MODE_515:
                homing_frame      = dhf_MODE_515;
                homing_frame_size = PRMS_MODE_515;
            break;
            case MODE_59:
                homing_frame      = dhf_MODE_59;
                homing_frame_size = PRMS_MODE_59;
            break;
            case MODE_67:
                homing_frame      = dhf_MODE_67;
                homing_frame_size = PRMS_MODE_67;
            break;
            case MODE_74:
                homing_frame      = dhf_MODE_74;
                homing_frame_size = PRMS_MODE_74;
            break;
            case MODE_795:
                homing_frame      = dhf_MODE_795;
                homing_frame_size = PRMS_MODE_795;
            break;
            case MODE_102:
                homing_frame      = dhf_MODE_102;
                homing_frame_size = PRMS_MODE_102;
            break;
            case MODE_122:
                homing_frame      = dhf_MODE_122;
                homing_frame_size = PRMS_MODE_122;
            break;
            default:
                homing_frame      = NULL;
                homing_frame_size = 0;
            break;
        }
        for(i=0; i<homing_frame_size; i++) {
            // check if the frame is homing
            if( (p->cur_frame_homing = p->amr_prms[i] ^ homing_frame[i]) ) break;
        }
    }

    if(p->cur_frame_homing) {
        decode_reset(avctx);
    }
    p->prev_frame_homing = !p->cur_frame_homing;
    p->prev_frame_type   =  p->cur_frame_type;
    p->prev_frame_mode   =  p->cur_frame_mode;

    /* To make it easy the stream can only be 16 bits mono,
     * so let's convert it to that */
    for (i=0 ; i<buf_size; i++)
        outbuffer[i] = (int16_t)p->sample_buffer[i];

    /* Report how many samples we got */
    *data_size = buf_size;

    /* Return the amount of bytes consumed if everything was OK */
    return *data_size*sizeof(int16_t);
}


static av_cold int amrnb_decode_close(AVCodecContext *avctx)
{
    AMRContext *p = avctx->priv_data;

    /* Free allocated memory */
    av_free(p->sample_buffer);
    av_free(p->amr_prms);
    av_free(p->state);

    /* Return 0 if everything is OK, -1 if not */
    return 0;
}


AVCodec amrnb_decoder = {
    .name           = "amrnb",
    .type           = CODEC_TYPE_AUDIO,
    .id             = CODEC_ID_AMR_NB,
    .priv_data_size = sizeof(AMRContext),
    .init           = amrnb_decode_init,
    .close          = amrnb_decode_close,
    .decode         = amrnb_decode_frame,
    .long_name      = NULL_IF_CONFIG_SMALL("Adaptive Multi-Rate NarrowBand"),
};
