/*
 * QCELP Decoder
 * Copyright (c) 2007 Reynaldo H. Verdejo Pinochet
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
 * @file qcelpdec.c
 * QCELP decoder.
 */

/* First we include some default includes */
#include <math.h>
#include <stddef.h>
#include <stdio.h>

/* The following includes have the bitstream reader, various dsp functions and the various defaults */
#define ALT_BITSTREAM_READER
#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"

#include "qcelp.h"

#define DEBUG 1

typedef struct {
    GetBitContext gb;
    QCELPFrame    *frame;
    uint8_t       erasure_count;
    uint8_t       ifq_count;
} QCELPContext;

static int qcelp_decode_init(AVCodecContext *avctx);
static int qcelp_decode_frame(AVCodecContext *avctx, void *data,
           int *data_size, uint8_t *buf, int buf_size);
static int qcelp_decode_close(AVCodecContext *avctx);


static int qcelp_decode_init(AVCodecContext *avctx)
{
    QCELPContext *q = (QCELPContext *) avctx->priv_data;

    avctx->sample_rate = 8000;
    avctx->channels = 1;

    q->frame = av_mallocz(sizeof(QCELPFrame));


    if(q->frame == NULL)
        return -1;

    return 0;
}

static int qcelp_parse_pkt_full(uint8_t *buf, QCELPFrame *frame)
{
    return 0;
}

static int qcelp_decode_frame(AVCodecContext *avctx, void *data,
           int *data_size, uint8_t *buf, int buf_size)
{
    QCELPContext *q = avctx->priv_data;
    int16_t *outbuffer = data;
    int8_t samples;
    int8_t bitcount;
    int16_t first16 = 0; /* needed for rate 1/8 particularities */

    QCELPBitmap *order;

    order = NULL;

    init_get_bits(&q->gb, buf, buf_size*8);

    /*
     * FIXME this comment should actually make some sence ..
     *
     * Here we try to identify each frame's rate by its byte size,
     * then, after setting a few utility vars we point 'order'
     * to start at the location of the rate's reference _slice_
     * inside the big REFERECE_FRAME array. We then proceed with
     * the bit reordering that will leave a full raw frame's data
     * ordered in our 'universal frame'
     */

    switch(buf_size)
    {
        case 34:
            q->frame->type = RATE_FULL;
            q->frame->bits = qcelp_bits_per_type[RATE_FULL];
            order = QCELP_REFERENCE_FRAME + QCELP_FULLPKT_REFERENCE_POS;
            break;
        case 16:
            q->frame->type = RATE_HALF;
            q->frame->bits = qcelp_bits_per_type[RATE_HALF];
            order = QCELP_REFERENCE_FRAME + QCELP_HALFPKT_REFERENCE_POS;
            break;
        case 7:
            q->frame->type = RATE_QUARTER;
            q->frame->bits = qcelp_bits_per_type[RATE_QUARTER];
            order = QCELP_REFERENCE_FRAME + QCELP_4THRPKT_REFERENCE_POS;
            break;
        case 3:
            q->frame->type = RATE_OCTAVE;
            q->frame->bits = qcelp_bits_per_type[RATE_OCTAVE];
            order = QCELP_REFERENCE_FRAME + QCELP_8THRPKT_REFERENCE_POS;
            break;
        case 0: /* FIXME */
            q->frame->type = BLANK;
            q->frame->bits = 0;
            break;
        default:
            q->frame->type = RATE_UNKNOWN;
            q->frame->bits = 0;
            /*
            printf("UNKNOWN PACKET RATE\n");
            */
            break;
    }

    /*
     * reordering loop
     */

    bitcount=0;
    while(bitcount < q->frame->bits)
    {
        /*
         * order[bitcount]->index holds the placement of this
         * input stream bit in the universal frame.
         *
         * order[bitcount]->pos holds the bit pos inside this value
         * byte.
         *
         */

        q->frame->data[ order[bitcount].index ] |=
        get_bits1(&q->gb)>>(order[bitcount].bitpos);

        /*
         * viral sample! :D
         *
         * just needed for rate 1/8 packets
         *
         */

        if(bitcount<16)
        {
            first16 |= q->frame->data[ order[bitcount].index ]>>bitcount
        }

        bitcount++;
    }

    /* DONE REORDERING */

    /*
     * check for erasures/blanks on rates 1, 1/4 and 1/8
     *
     */

    if(q->frame->type != RATE_HALF)
    {
        if(!q->frame->data[QCELP_RSRVD_POS])
        {
            /*
             * flag aproach: set flag for ifq/blank/incorrect
             * decoding
             */
        }
    }

    /* particularities for rate 1/8 */
    if(q->frame->type == RATE_OCTAVE)
    {
        if(first16==0xFFFF)
        {
            /*
             * flag aproach: set flag for ifq/blank/incorrect
             * decoding
             */
        }
    }

    /*
     * check for badly received packets
     * for rate 1, 1/2 and 1/4
     */

    if(q->frame->type != RATE_OCTAVE)
    {
        /*
         * flag aproach: set flag for ifq/blank/incorrect
         * decoding
         */

    }



     /*
      * decode loop
      *
      */




    return 1;
}

AVCodec qcelp_decoder =
{
    .name = "qcelp",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_QCELP,
    .priv_data_size = sizeof(QCELPContext),
    .init = qcelp_decode_init,
    .close = qcelp_decode_close,
    .decode = qcelp_decode_frame,
};
