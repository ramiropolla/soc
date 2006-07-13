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
    float*                    sample_buffer;

    int                   prev_frame_homing; // previous frame was a homing frame ? 1 : 0
    enum RXFrameType        prev_frame_type; // frame type of previous frame
    enum Mode               prev_frame_mode; // mode of previous frame

    struct AMRDecoderState {                 // struct to hold current decoder state
    }; AMRDecoderState

} AMRContext;


static int amr_nb_decode_init(AVCodecContext *avctx) {

    AMRContext *p = avctx->priv_data;

    /* We also need to allocate a sample buffer */
    p->sample_buffer = av_mallocz(sizeof(float)*1024);  // here we used av_mallocz instead of av_malloc
                                                        // av_mallocz memsets the whole buffer to 0

    /* Check if the allocation was successful */
    if(p->sample_buffer == NULL)
        return -1;

    /* return 0 for a successful init, -1 for failure */
    return 0;
}


static int amr_nb_decode_frame(AVCodecContext *avctx,
        void *data, int *data_size, uint8_t *buf, int buf_size) {

    AMRContext *p = avctx->priv_data;
    int16_t *outbuffer = data; // FIXME check possible output data type(s)
    int i;
    int16_t amr_prms, q_bit;

    enum RXFrameType frame_type;
    enum Mode mode, speech_mode;

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
    mode = decode_bitstream(avctx, &amr_prms, buf, buf_size, &frame_type, &speech_mode, &q_bit);

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

    /* Free allocated memory buffer */
    av_free(p->sample_buffer);

    /* Return 0 if everything is ok, -1 if not */
    return 0;
}

/**
 * Decode the bitstream into the AMR parameters and discover the frame mode
 *
 * @param amr_prms          pointer to the AMR parameters
 * @param buf               pointer to the input buffer
 * @param frame_type        pointer to the frame type
 * @param speech_mode       pointer to the speech mode
 * @param q_bit             pointer to q_bit which is used to decide if the frame
 *                          is good or bad
 *
 * @return Returns the frame mode
 */

enum Mode decode_bitstream(AVCodecContext *avctx, int16_t *amr_prms, uint8_t *buf, int buf_size,
                           enum RXFrameType *frame_type, enum Mode *speech_mode, int16_t *q_bit) {

    AMRContext *p = avctx->priv_data;
    enum Mode mode;
    int i;
    int16_t *mask;

    // initialise get_bits
    init_get_bits(&p->gb, buf, buf_size*8);
    memset(amr_prms, 0, PRMS_MODE_122 << 1);
    skip_bits1(&p->gb);
    mode = get_bits(&p->gb ,4);
    *q_bit = get_bits1(&p->gb); // FIXME rename q_bit to something more meaningful when i understand what it is
    skip_bits(&p->gb, 2);

    switch(mode) {
        case MODE_DTX:
            mask = order_MODE_DTX;
            *frame_type = RX_SID_FIRST; // get SID type bit
        break;
        case NO_DATA:
            *frame_type = RX_NO_DATA;
        break;
        case MODE_475:
            mask = order_MODE_475;
            *frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_515:
            mask = order_MODE_515;
            *frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_59:
            mask = order_MODE_59;
            *frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_67:
            mask = order_MODE_67;
            *frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_74:
            mask = order_MODE_74;
            *frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_795:
            mask = order_MODE_795;
            *frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_102:
            mask = order_MODE_102;
            *frame_type = RX_SPEECH_GOOD;
        break;
        case MODE_122:
            mask = order_MODE_122;
            *frame_type = RX_SPEECH_GOOD;
        break;
        default:
            *frame_type = RX_SPEECH_BAD;
        break;
    }

    if((*frame_type != RX_NO_DATA) && (*frame_type != RX_SPEECH_BAD)) {
        for(i=1; i<mode_bits[mode]; i++) {
            amr_prms[*mask] += get_bits1(&p->gb) * mask[1];
            mask += 2;
        }
    }

    if(mode == MODE_DTX) {
        skip_bits(&p->gb, 4); // skip to the next byte
        if(get_bits1(&p->gb)) // use the update if there is one
            *frame_type = RX_SID_UPDATE;
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
