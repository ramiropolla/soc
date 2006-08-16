/*
 * AMR audio codec
 * Copyright (c) 2006 Robert Swain
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */


/**
 * @file amr.c
 * AMR codec.
 */


#include <math.h>
#include <stddef.h>
#include <stdio.h>

#include "avcodec.h"
#include "bitstream.h"
#include "amrdata.h"

// #define DEBUG_BITSTREAM

typedef struct AMRContext {

    GetBitContext                        gb;
    float                    *sample_buffer;

    int16_t                       *amr_prms; // pointer to the decoded amr parameters (lsf coefficients, codebook indices, etc)
    int                 bad_frame_indicator;

    int                   prev_frame_homing; // previous frame was a homing frame ? 1 : 0
    enum RXFrameType        prev_frame_type; // frame type of previous frame
    enum Mode               prev_frame_mode; // mode of previous frame

    int                        *prev_lsf_dq; // previous dequantised lsfs
    int                    *prev_residual_q; // previous quantised residual

    int             lsp1_q[LP_FILTER_ORDER]; // vector of quantised lsps
    int             lsp2_q[LP_FILTER_ORDER]; // vector of quantised lsps

    int                    cur_frame_homing; // current frame homing ? 1 : 0
    enum RXFrameType         cur_frame_type; // current frame type
    enum Mode                cur_frame_mode; // current frame mode

    int    lpc_coeffs[4][LP_FILTER_ORDER+1]; // lpc coefficients, A(z), for all four subframes
    int      prev_lsp_sub4[LP_FILTER_ORDER]; // vector of lsps from subframe 4 of the previous frame
    int       cur_lsp_sub2[LP_FILTER_ORDER]; // vector of lsps from subframe 2 of the current frame
    int       cur_lsp_sub4[LP_FILTER_ORDER]; // vector of lsps from subframe 4 of the current frame

//    struct AMRDecoderState {                 // struct to hold current decoder state
//    }; AMRDecoderState

} AMRContext;


static int amr_nb_decode_init(AVCodecContext *avctx) {

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
    p->sample_buffer = av_mallocz(sizeof(float)*1024);
    // allocate and zero the amr parameters
    p->amr_prms = av_mallocz(sizeof(int16_t)*PRMS_MODE_122);

    /* Check if the allocation was successful */
    if(p->sample_buffer == NULL)
        return -1;

    /* return 0 for a successful init, -1 for failure */
    return 0;
}


static int amr_nb_decode_frame(AVCodecContext *avctx,
        void *data, int *data_size, uint8_t *buf, int buf_size) {

    AMRContext *p = avctx->priv_data;        // pointer to private data
    int16_t *outbuffer = data;               // pointer to the output data buffer
    int i;                                   // counter
    enum Mode speech_mode = MODE_475;        // ???
    const int16_t *homing_frame;             // pointer to the homing frame
    int homing_frame_size;                   // homing frame size

#ifdef DEBUG_BITSTREAM
    init_get_bits(&p->gb, buf, buf_size*8);
    av_log(NULL, AV_LOG_ERROR, "\n\n\nBits from one frame (%d):\n", buf_size*8);
    for(i=0; i<buf_size*8; i++) {
        av_log(NULL, AV_LOG_ERROR, "%d",get_bits1(&p->gb));
    }
    av_log(NULL, AV_LOG_ERROR, "\n\n\n");
    return -1;
#endif // DEBUG_BITSTREAM

    // decode the bitstream to amr parameters
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
                homing_frame = dhf_MODE_475;
                homing_frame_size = 7;
            break;
            case MODE_515:
                homing_frame = dhf_MODE_515;
                homing_frame_size = 7;
            break;
            case MODE_59:
                homing_frame = dhf_MODE_59;
                homing_frame_size = 7;
            break;
            case MODE_67:
                homing_frame = dhf_MODE_67;
                homing_frame_size = 7;
            break;
            case MODE_74:
                homing_frame = dhf_MODE_74;
                homing_frame_size = 7;
            break;
            case MODE_795:
                homing_frame = dhf_MODE_795;
                homing_frame_size = 8;
            break;
            case MODE_102:
                homing_frame = dhf_MODE_102;
                homing_frame_size = 12;
            break;
            case MODE_122:
                homing_frame = dhf_MODE_122;
                homing_frame_size = 18;
            break;
            default:
                homing_frame = NULL;
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
                homing_frame = dhf_MODE_475;
                homing_frame_size = PRMS_MODE_475;
            break;
            case MODE_515:
                homing_frame = dhf_MODE_515;
                homing_frame_size = PRMS_MODE_515;
            break;
            case MODE_59:
                homing_frame = dhf_MODE_59;
                homing_frame_size = PRMS_MODE_59;
            break;
            case MODE_67:
                homing_frame = dhf_MODE_67;
                homing_frame_size = PRMS_MODE_67;
            break;
            case MODE_74:
                homing_frame = dhf_MODE_74;
                homing_frame_size = PRMS_MODE_74;
            break;
            case MODE_795:
                homing_frame = dhf_MODE_795;
                homing_frame_size = PRMS_MODE_795;
            break;
            case MODE_102:
                homing_frame = dhf_MODE_102;
                homing_frame_size = PRMS_MODE_102;
            break;
            case MODE_122:
                homing_frame = dhf_MODE_122;
                homing_frame_size = PRMS_MODE_122;
            break;
            default:
                homing_frame = NULL;
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
    p->prev_frame_type   = p->cur_frame_type;
    p->prev_frame_mode   = p->cur_frame_mode;

    /* To make it easy the stream can only be 16 bits mono, so let's convert it to that */
    for (i=0 ; i<buf_size; i++)
        outbuffer[i] = (int16_t)p->sample_buffer[i];

    /* Report how many samples we got */
    *data_size = buf_size;

    /* Return the amount of bytes consumed if everything was ok */
    return *data_size*sizeof(int16_t);
}


static int amr_nb_decode_close(AVCodecContext *avctx) {

    AMRContext *p = avctx->priv_data;

    /* Free allocated memory */
    av_free(p->sample_buffer);
    av_free(p->amr_prms);

    /* Return 0 if everything is ok, -1 if not */
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

enum Mode decode_bitstream(AVCodecContext *avctx, uint8_t *buf, int buf_size, enum Mode *speech_mode) {

    AMRContext *p = avctx->priv_data;
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


/**
 * Convert a set of 3 lsf quantisation indices into lsp parameters
 *
 * @param avctx             pointer to the AVCodecContext for AMR
 */

static void decode_lsf2lsp_3(AVCodecContext *avctx) {

    AMRContext *p = avctx->priv_data;

    int lsf1_r[LP_FILTER_ORDER]; // vector of residual lsfs
    int lsf1_q[LP_FILTER_ORDER]; // vector of quantised lsfs
    const int16_t (*lsf_3_temp1)[3], (*lsf_3_temp3)[4]; // temp ptrs for switching tables depending on mode
    int index_temp; // temp lsf index
    int i; // counter

    // if the current frame is bad estimate the past quantised residual based on the past lsf shifted slightly towards the mean
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

        // decode lsf residuals from lsf tables using the indices from the bitstream
        lsf1_r[0] = lsf_3_temp1[ p->amr_prms[0] ][0];
        lsf1_r[1] = lsf_3_temp1[ p->amr_prms[0] ][1];
        lsf1_r[2] = lsf_3_temp1[ p->amr_prms[0] ][2];

        index_temp = p->amr_prms[1];
        // MODE_475, MODE_515 only use every second entry - NOTE : not quite sure what they mean and the spec doesn't mention this
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

    /* verification that LSFs has minimum distance of LSF_GAP Hz */
    reorder_lsf(lsf1_q, LSF_GAP);
    memcpy(p->prev_lsf_dq, lsf1_q, LP_FILTER_ORDER<<2);

    /*  convert LSFs to the cosine domain */
    lsf2lsp(lsf1_q, p->lsp1_q);
}


/**
 * Convert a set of 5 lsf quantisation indices into lsp parameters
 *
 * @param avctx             pointer to the AVCodecContext for AMR
 */

static void decode_lsf2lsp_5(AVCodecContext *avctx) {

    AMRContext *p = avctx->priv_data;

    int lsf1_r[LP_FILTER_ORDER], lsf2_r[LP_FILTER_ORDER]; // vectors of residual lsfs
    int lsf1_q[LP_FILTER_ORDER], lsf2_q[LP_FILTER_ORDER]; // vectors of quantised lsfs
    int temp;
    int i, sign; // counter and sign of 3rd lsf table

    // if the current frame is bad estimate the past quantised residual based on the past lsf shifted slightly towards the mean
    if(p->bad_frame_indicator) {
        for(i=0; i<LP_FILTER_ORDER; i++) {
            lsf1_q[i] = ( (p->prev_lsf_dq[i]*ALPHA_122)>>15 ) + ( (lsf_5_mean[i]*ONE_ALPHA_122)>>15 );
        }
        memcpy(lsf2_q, lsf1_q, LP_FILTER_ORDER<<1);
        for(i=0; i<LP_FILTER_ORDER; i++) {
            p->prev_residual_q[i] = lsf2_q[i] - lsf_5_mean[i] - ( (p->prev_residual_q[i]*LSP_PRED_FAC_MODE_122)>>15 );
        }
    }else {
        /* decode prediction residuals from 5 received indices */
        lsf1_r[0] = lsf_5_1[ p->amr_prms[0] ][0];
        lsf1_r[1] = lsf_5_1[ p->amr_prms[0] ][1];
        lsf2_r[0] = lsf_5_1[ p->amr_prms[0] ][2];
        lsf2_r[1] = lsf_5_1[ p->amr_prms[0] ][3];

        lsf1_r[2] = lsf_5_2[ p->amr_prms[1] ][0];
        lsf1_r[3] = lsf_5_2[ p->amr_prms[1] ][1];
        lsf2_r[2] = lsf_5_2[ p->amr_prms[1] ][2];
        lsf2_r[3] = lsf_5_2[ p->amr_prms[1] ][3];

        sign = (p->amr_prms[2] & 1) ? -1 : 1;
        // I don't know why p->amr_prms[2]>>1 but that's how it is in the ref source
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
 * DESCRIPTION FROM REF SOURCE:
 * Make sure that the LSFs are properly ordered and to keep a certain minimum
 * distance between adjacent LSFs.
 *
 * @param lsf               a vector of lsfs (range: 0<=val<=0.5)
 * @param min_dist          minimum required separation of lsfs
 */

static void reorder_lsf(int *lsf, int min_dist) {
    int i;
    int lsf_min;

    for(i=0; i<LP_FILTER_ORDER; i++) {
        if(lsf[i] < lsf_min) {
            lsf[i] = lsf_min;
        }
        lsf_min = lsf[i] + min_dist;
    }
}


/**
 * Convert a vector of lsfs into the corresponding, cosine domain lsps
 *
 * @param lsf               a vector of lsfs
 * @param lsp               a vector of lsps
 */

static void lsf2lsp(int *lsf, int *lsp) {
    int i;
    int index, offset;

    for(i=0; i<LP_FILTER_ORDER; i++) {
        index = lsf[i] >> 8;      // bits 8 to 15 of lsf[i]
        offset = lsf[i] & 0x00ff; // bits 0 to  7 of lsf[i]
        lsp[i] = cos_table[index] + ( (( cos_table[index+1]-cos_table[index] )*offset)>>8 );
    }
}


/**
 * Find the polynomial F1(z) or F2(z) from the lsps
 *
 * The ref source description of how to find the polynomials was as follows:
 *    Find the polynomial F1(z) or F2(z) from the LSPs.
 *
 *    F1(z) = product ( 1 - 2 lsp[i] z^-1 + z^-2 )
 *             i=0,2,4,6,8
 *    F2(z) = product   ( 1 - 2 lsp[i] z^-1 + z^-2 )
 *             i=1,3,5,7,9
 *
 * @param lsp               vector of lsps
 * @param f                 pointer to the polynomial F1(z) or F2(z)
 *
 * @return void
 */

static void lsp2poly(int *lsp, int *f) {
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
/* FIXME - check rounding sensitivity of F1(z) and F2(z) as the code below is faster but rounds differently
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
 * Convert a vector of lsps to lpc coefficients for a 10th order filter
 *
 * @param lsp                 vector of lsps
 * @param lpc_coeffs          pointer to a subframe of lpc coefficients
 *
 * @return void
 */

static void lsp2lpc(int *lsp, int *lpc_coeffs) {
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
        if(( temp & 0x1000 )) {
            lpc_coeffs[i]++;
        }
        temp = f1[11-i] - f2[11-i];
        lpc_coeffs[11-i] = (int16_t)(temp>>13);  // ref source: emulate fixed point bug
        if(( temp & 0x1000 )) {
            lpc_coeffs[11-i]++;
        }
    }
}


/**
 * Interpolate lsps for subframes 1 and 3 and convert them to lpc coefficients
 *
 * @param avctx                 pointer to AVCodecContext
 * @param lpc_coeffs            array of lpc coefficients for the four subframes
 *
 * @return void
 */

static void lpc_interp_13(AVCodecContext *avctx, int **lpc_coeffs) {
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
 * Interpolate lsps for subframes 1, 2 and 3 and convert them to lpc coefficients
 *
 * @param avctx                 pointer to AVCodecContext
 * @param lpc_coeffs            array of lpc coefficients for the four subframes
 *
 * @return void
 */

static void lpc_interp_123(AVCodecContext *avctx, int **lpc_coeffs) {
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
 * Reset the AMR frame parameters
 *
 * @param avctx             pointer to the AVCodecContext for AMR
 */

void decode_reset(AVCodecContext *avctx) {
    AMRContext *p = avctx->priv_data;

    p->prev_frame_homing = 1;
    p->prev_frame_type = RX_SPEECH_GOOD;
    p->prev_frame_mode = MODE_475;
    // FIXME reset AMRDecoderState too!
}

AVCodec amr_nb_decoder =
{
    .name = "amr_nb",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_AMR_NB,
    .priv_data_size = sizeof(AMRContext),
    .init = amr_nb_decode_init,
    .close = amr_nb_decode_close,
    .decode = amr_nb_decode_frame,
};
