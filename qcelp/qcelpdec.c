/*
 * QCELP decoder
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
 * QCELP decoder
 */

#include <math.h>
#include <stddef.h>

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"

#include "qcelpdata.h"

#define DEBUG 1

typedef struct
{
    qcelp_packet_rate rate;
    uint8_t data[76];       /*!< data from a _parsed_ frame */
    uint8_t bits;
} QCELPFrame;

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

static int qcelp_decode_close(AVCodecContext *avctx)
{
    QCELPContext *q = avctx->priv_data;

    av_free(q->frame);

    return 0;
}

static int qcelp_decode_frame(AVCodecContext *avctx, void *data,
           int *data_size, uint8_t *buf, int buf_size)
{
    QCELPContext *q    = avctx->priv_data;
    const QCELPBitmap *order = NULL;
    int16_t  *outbuffer = data;
    int8_t   samples;
    int      n;
    uint16_t first16 = 0; /*!< needed for rate 1/8 peculiarities */
    int      is_ifq = 0;

    init_get_bits(&q->gb, buf, buf_size*8);

    /**
     * figure out frame's rate by its size, set up a few utility vars
     * and point 'order' to the rate's reference _slice_ inside the
     * big REFERENCE_FRAME array.
     */

    switch(buf_size)
    {
        case 34:
            q->frame->rate = RATE_FULL;
            q->frame->bits = qcelp_bits_per_rate[RATE_FULL];
            order = QCELP_REFERENCE_FRAME + QCELP_FULLPKT_REFERENCE_POS;
            break;
        case 16:
            q->frame->rate = RATE_HALF;
            q->frame->bits = qcelp_bits_per_rate[RATE_HALF];
            order = QCELP_REFERENCE_FRAME + QCELP_HALFPKT_REFERENCE_POS;
            break;
        case 7:
            q->frame->rate = RATE_QUARTER;
            q->frame->bits = qcelp_bits_per_rate[RATE_QUARTER];
            order = QCELP_REFERENCE_FRAME + QCELP_4THRPKT_REFERENCE_POS;
            break;
        case 3:
            q->frame->rate = RATE_OCTAVE;
            q->frame->bits = qcelp_bits_per_rate[RATE_OCTAVE];
            order = QCELP_REFERENCE_FRAME + QCELP_8THRPKT_REFERENCE_POS;
            break;
        case 0: /* FIXME */
            q->frame->rate = BLANK;
            q->frame->bits = 0;
            break;
        default:
            q->frame->rate = RATE_UNKNOWN;
            q->frame->bits = 0;
            /*
            printf("UNKNOWN PACKET RATE\n");
            */
            break;
    }

    /**
     * reordering loop
     */

    for(n=0; n < q->frame->bits; n++)
    {
        q->frame->data[ order[n].index ] |=
        get_bits1(&q->gb)<<order[n].bitpos;

        if(n<16)
        {
            first16 |= q->frame->data[ order[n].index ]>>n;
        }

    }

    /* DONE REORDERING */

    /**
     * check for erasures/blanks on rates 1, 1/4 and 1/8
     */

    if(q->frame->rate != RATE_HALF && !q->frame->data[QCELP_RSRVD_POS])
        is_ifq=1;

    if(q->frame->rate == RATE_OCTAVE && first16==0xFFFF)
        is_ifq=1;

    /* check for badly received packets */

    if(q->frame->rate != RATE_OCTAVE)
    {
        /* check for outbound LSP freqs and codebook gain params */
        if(q->frame->rate != RATE_QUARTER)
        {
               /* magic here */
        }

    }

     /*
      * decode loop
      *
      */


    return 1;
}

AVCodec qcelp_decoder =
{
    .name   = "qcelp",
    .type   = CODEC_TYPE_AUDIO,
    .id     = CODEC_ID_QCELP,
    .init   = qcelp_decode_init,
    .close  = qcelp_decode_close,
    .decode = qcelp_decode_frame,
    .priv_data_size = sizeof(QCELPContext),
};
