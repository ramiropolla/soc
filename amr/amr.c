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

    int                    cur_frame_homing; // current frame homing ? 1 : 0
    enum RXFrameType         cur_frame_type; // current frame type
    enum Mode                cur_frame_mode; // current frame mode

//    struct AMRDecoderState {                 // struct to hold current decoder state
//    }; AMRDecoderState

} AMRContext;


static int amr_nb_decode_init(AVCodecContext *avctx) {

    AMRContext *p = avctx->priv_data;

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
    int16_t q_bit;                           // FIXME rename q_bit when I know what it means
    enum Mode speech_mode = MODE_475;        // ???
    const int16_t *homing_frame;             // pointer to the homing frame
    int16_t homing_frame_size;               // homing frame size

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
    p->cur_frame_mode = decode_bitstream(avctx, buf, buf_size, &speech_mode, &q_bit);

    // set the bad frame indicator depending on q_bit
    if(!p->bad_frame_indicator) p->bad_frame_indicator = !q_bit;

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
            if(p->cur_frame_homing = p->amr_prms[i] ^ homing_frame[i]) break;
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
            if(p->cur_frame_homing = p->amr_prms[i] ^ homing_frame[i]) break;
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
        outbuffer[i] = (int16_t)p->sample_buffer[i]; // FIXME check possible output data type(s)

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
 * @param speech_mode       pointer to the speech mode
 * @param q_bit             pointer to q_bit which is used to decide if the frame
 *                          is good or bad
 *
 * @return Returns the frame mode
 */

enum Mode decode_bitstream(AVCodecContext *avctx, uint8_t *buf, int buf_size, enum Mode *speech_mode, int16_t *q_bit) {

    AMRContext *p = avctx->priv_data;
    enum Mode mode;
    int i;
    AMROrder *order;

    // initialise get_bits
    init_get_bits(&p->gb, buf, buf_size*8);
    skip_bits1(&p->gb);
    mode = get_bits(&p->gb ,4);
    *q_bit = get_bits1(&p->gb); // FIXME rename q_bit to something more meaningful when i understand what it is
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
