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

/**
 * Decodes the five |R2 LSPVi vectors to get the 10
 * quantized LSP frequencies from any packet rate but 1/8
 *
 * For details see TIA/EIA/IS-733 2.4.3.2.6.2-2
 */
void qcelp_decode_lspf(const QCELPFrame *frame, float *lspf)
{
    /* FIXME a loop is wanted here */
    /* WIP implement I_F_Q handling? */
    const uint8_t *lspv;
    int i;

    switch(frame->rate)
    {
        case RATE_FULL:
        case RATE_HALF:
        case RATE_QUARTER:
            lspv=frame->data+QCELP_LSPV0_POS;
            av_log(NULL, AV_LOG_DEBUG,"LSPV:%5d,%5d,%5d,%5d,%5d\n",
                   lspv[0],lspv[1],lspv[2],
                   lspv[3],lspv[4]);
            lspf[0]=        qcelp_lspvq1[lspv[0]].x;
            lspf[1]=lspf[0]+qcelp_lspvq1[lspv[0]].y;
            lspf[2]=lspf[1]+qcelp_lspvq2[lspv[1]].x;
            lspf[3]=lspf[2]+qcelp_lspvq2[lspv[1]].y;
            lspf[4]=lspf[3]+qcelp_lspvq3[lspv[2]].x;
            lspf[5]=lspf[4]+qcelp_lspvq3[lspv[2]].y;
            lspf[6]=lspf[5]+qcelp_lspvq4[lspv[3]].x;
            lspf[7]=lspf[6]+qcelp_lspvq4[lspv[3]].y;
            lspf[8]=lspf[7]+qcelp_lspvq5[lspv[4]].x;
            lspf[9]=lspf[8]+qcelp_lspvq5[lspv[4]].y;
            break;
        case RATE_OCTAVE:
            lspv=frame->data+QCELP_LSP0_POS;
            for(i=0; i<10; i++)
            {
                lspf[i]=lspv[i]? 0.02:-0.02; /* 2.4.3.3.1-1 */
            }
    }
}

/**
 * Converts codebook transmission codes to GAIN and INDEX
 * (and cbseed for rate 1/4)
 * TIA/EIA/IS-733 2.4.6.2
 */
void qcelp_decode_params(const QCELPFrame *frame, int *g0, uint16_t *cbseed,
     float *gain, int *index)
{
    int           i, gs[16], g1[16], predictor;
    const uint8_t *cbgain, *cbsign, *cindex, *data;
    float         ga[16];

    /* FIXME need to get rid of g0, sanity checks should be done here */
    /* WIP this is almost verbatim from spec, seeking workability first */

    cbsign=frame->data+QCELP_CBSIGN0_POS;
    cbgain=frame->data+QCELP_CBGAIN0_POS;
    cindex=frame->data+QCELP_CINDEX0_POS;

    switch(frame->rate)
    {
        case RATE_FULL:
        case RATE_HALF:
            for(i=0; i<16; i++)
            {
                if(frame->rate == RATE_HALF && i>=4) break;

                gs[i]=QCELP_CBSIGN2GS(cbsign[i]);
                g0[i]=QCELP_CBGAIN2G0(cbgain[i]);

                /* FIXME this needs to be further examinated */
                if(frame->rate == RATE_FULL && i > 0 && !((i+1) & 3))
                    predictor=av_clip(6, 38, (g1[i-1]+g1[i-2]+g1[i-3])/3);
                else
                    predictor=0;

                g1[i]=g0[i]+predictor;
                ga[i]=qcelp_g12ga[g1[i]];

                gain[i]=ga[i]*gs[i];
                index[i]=(gs[i] > 0)? cindex[i]:(cindex[i]-89) & 127;/* FIXME */
            }

            break;
        case RATE_QUARTER:
            for(i=0; i<5; i++)
            {
                g0[i]=g1[i]=QCELP_CBGAIN2G0(cbgain[i]);
                gs[i]=1;
                ga[i]=qcelp_g12ga[g1[i]];
            }
            /**
             * 5->8 Interpolation to 'Provide smoothing of the energy
             * of the unvoiced excitation' 2.4.6.2
             */
            gain[0]=    ga[0];
            gain[1]=0.6*ga[0]+0.4*ga[1];
            gain[2]=    ga[1];
            gain[3]=0.2*ga[1]+0.8*ga[2];
            gain[4]=0.8*ga[2]+0.2*ga[3];
            gain[5]=    ga[3];
            gain[7]=0.4*ga[3]+0.6*ga[4];
            gain[7]=    ga[4];
            /**
             * Build random* seed needed to make Cdn
             */
            data=frame->data;
            *cbseed=(0x0003 & data[QCELP_LSPV0_POS+4])<<14 |
                    (0x003C & data[QCELP_LSPV0_POS+3])<< 8 |
                    (0x0060 & data[QCELP_LSPV0_POS+2])<< 1 |
                    (0x0007 & data[QCELP_LSPV0_POS+1])<< 3 |
                    (0x0038 & data[QCELP_LSPV0_POS  ])>> 3 ;
            break;
        case RATE_OCTAVE:
            switch(cbgain[0])
            {
                case 0: gain[0]=-4; break;
                case 1: gain[0]=-2; break;
                case 2: gain[0]= 0; break;
                case 3: gain[0]= 2; break;
                default:; /* shouldn't happen.. must propagate some error */
            }
            gs[0]=1;
            /* WIP finish rate 1/8 calculations, spec is kind of fuzzy here */
    }
}

/**
 * Computes the scaled codebook vector Cdn From INDEX and GAIN
 * For all rates
 *
 * FIXME:
 * - Needs outbound reading checks and error propagation if weirdness
 *   is detected :-).
 */
static int qcelp_compute_svector(qcelp_packet_rate rate, const float *gain,
           const int *index, uint16_t cbseed, float *cdn_vector)
{
    int      i,j;
    uint16_t new_cbseed;
    float    rnd[160];

    switch(rate)
    {
        case RATE_FULL:
             for(i=0; i<160; i++)
                cdn_vector[i]=
                gain[i/10]*qcelp_fullrate_ccodebook[(i+1-index[i/10]) & 127];
             break;
        case RATE_HALF:
             for(i=0; i<160; i++)
                cdn_vector[i]=
                gain[i/40]*qcelp_halfrate_ccodebook[(i+1-index[i/40]) & 127];
             break;
        case RATE_QUARTER:
            for(i=0; i<160; i++)
            {
                new_cbseed=(521*cbseed+259) & 65535;
                cbseed=rnd[i]=
                QCELP_SQRT1887*(((new_cbseed+32768) & 65535)-32768)/32768.0;

                /* FIR filter */
                cdn_vector[i]=qcelp_rnd_fir_coefs[1]*rnd[i];
                for(j=1; j<22 && !(i-j+1); j++)
                {
                    cdn_vector[i]+=qcelp_rnd_fir_coefs[j]*rnd[i-j];
                }
                /* final scaling */
                cdn_vector[i]*=gain[i/20];
            }
            break;
        case RATE_OCTAVE:
            for(i=0; i<160; i++)
            {
                new_cbseed=(521*cbseed+259) & 65535;
                cbseed=rnd[i]=
                QCELP_SQRT1887*(((new_cbseed+32768) & 65535)-32768)/32768.0;

                cdn_vector[i]=gain[0]*rnd[i];
            }
    }

    return 1;
}

/**
 * Computes hammsinc(x), this will probably be replaced
 * by a lookup table
 */
static float qcelp_hammsinc(int i)
{
    return (sin(M_PI*i)/(M_PI*i))*(0.5+0.46*cos(M_PI*i/4));
}

/**
 * pitch filters cdn_vector, returns 0 if everything goes well, otherwise
 * it returns the index of the failing-to-be-pitched element or -1 if
 * an invalid (140.5, 141.5, 142.5) fractional lag is found.
 *
 * For details see 2.4.5.2
 */
static int qcelp_do_pitchfilter(QCELPFrame *frame, float *cdn_vector)
{
    int     i,j;
    uint8_t *pgain, *plag, *pfrac;
    float   gain, lag;

    switch(frame->rate)
    {
        case RATE_FULL:
        case RATE_HALF:
            pgain=frame->data+QCELP_PGAIN0_POS;
            plag =frame->data+QCELP_PLAG0_POS;
            pfrac=frame->data+QCELP_PFRAC0_POS;

            for(i=0; i<160; i++)
            {
                gain=plag[i/40]? (pgain[i/40]+1)/4.0 : 0.0;
                lag =plag[i/40]+16+0.5*pfrac[i/40];

                if(lag == 140.5 || lag == 141.5 || lag == 142.5)
                    return -1;

                if(pfrac[i/40]==1) /* if lag is a fractional lag... */
                {
                    /* see equation 2.4.5.2-2 */
                    for(j=-4; j<4; j++)
                    {
                        if(i + j - lag + 0.5 < 0)
                            break; /*XXX may be unneded */

                        cdn_vector[i]+=
                        gain*qcelp_hammsinc(j+0.5)*
                        cdn_vector[i + j - (int)(lag + 0.5)];
                    }
                    if(j<4) break;

                }else if(i >= lag) /*XXX may be unneded*/
                    cdn_vector[i]+=gain*cdn_vector[i - (int)lag];
                else
                    break; /* shouldn't happen */
            }
            if(i < 160) return i;

            break;
        case RATE_QUARTER:
        case RATE_OCTAVE:
            /* ? */
        break;
    }

    return 0;
}

static int qcelp_decode_frame(AVCodecContext *avctx, void *data,
           int *data_size, uint8_t *buf, int buf_size)
{
    QCELPContext *q    = avctx->priv_data;
    const QCELPBitmap *order = NULL;
    int16_t  *outbuffer = data, cbseed;
    int      n, is_ifq = 0, is_codecframe_fmt = 0;
    uint16_t first16 = 0;
    float    qtzd_lspf[10], gain[16], cdn_vector[160];
    int      g0[16], index[16];
    uint8_t  claimed_rate;

    init_get_bits(&q->gb, buf, buf_size*8);

    /**
     * figure out frame's rate by its size, set up a few utility vars
     * and point 'order' to the rate's reference _slice_ inside the
     * big REFERENCE_FRAME array.
     */

    switch(buf_size)
    {
        case 35:
        case 34:
            q->frame->rate = RATE_FULL;
            q->frame->bits = qcelp_bits_per_rate[RATE_FULL];
            order = QCELP_REFERENCE_FRAME + QCELP_FULLPKT_REFERENCE_POS;
            if(buf_size == 35) is_codecframe_fmt=1;
            break;
        case 17:
        case 16:
            q->frame->rate = RATE_HALF;
            q->frame->bits = qcelp_bits_per_rate[RATE_HALF];
            order = QCELP_REFERENCE_FRAME + QCELP_HALFPKT_REFERENCE_POS;
            if(buf_size == 17) is_codecframe_fmt=1;
            break;
        case  8:
        case  7:
            q->frame->rate = RATE_QUARTER;
            q->frame->bits = qcelp_bits_per_rate[RATE_QUARTER];
            order = QCELP_REFERENCE_FRAME + QCELP_4THRPKT_REFERENCE_POS;
            if(buf_size == 8) is_codecframe_fmt=1;
            break;
        case  4:
        case  3:
            q->frame->rate = RATE_OCTAVE;
            q->frame->bits = qcelp_bits_per_rate[RATE_OCTAVE];
            order = QCELP_REFERENCE_FRAME + QCELP_8THRPKT_REFERENCE_POS;
            if(buf_size == 9) is_codecframe_fmt=1;
            break;
        case  0: /* FIXME */
        case  1:
            q->frame->rate = BLANK;
            q->frame->bits = 0;
            if(buf_size == 1) is_codecframe_fmt=1;
            break;
        default:
            q->frame->rate = RATE_UNKNOWN;
            q->frame->bits = 0;
            /*
            printf("UNKNOWN PACKET RATE\n");
            */
    }

    if(is_codecframe_fmt)
    {
        claimed_rate=get_bits(&q->gb, 8);

        if((claimed_rate ==  0 && q->frame->rate != BLANK       ) ||
           (claimed_rate ==  1 && q->frame->rate != RATE_OCTAVE ) ||
           (claimed_rate ==  2 && q->frame->rate != RATE_QUARTER) ||
           (claimed_rate ==  3 && q->frame->rate != RATE_HALF   ) ||
           (claimed_rate ==  4 && q->frame->rate != RATE_FULL   ))
        {
           av_log(NULL, AV_LOG_ERROR,
           "Claimed rate and buffer size missmatch\n");
           is_ifq=1;
        }
    }

    av_log(NULL, AV_LOG_INFO, "Rate %d Size %d\n", q->frame->rate, buf_size);

    /**
     * reordering loop
     */
    memset(q->frame->data, 0, 76);
    for(n=0; n < q->frame->bits; n++)
    {
        q->frame->data[ order[n].index ] |=
        get_bits1(&q->gb)<<order[n].bitpos;

        /* FIXME Should rework this a bit */
        if(n<20)
        {
            if(n>3)  /* this is the random seed for rate 1/8 frames */
                cbseed |= q->frame->data[ order[n].index ]>>n;
            if(n<16) /* this is for a rate 1/8 only sanity check */
                first16 |= q->frame->data[ order[n].index ]>>n;
        }

    }

    /* skip padding bits */
    skip_bits(&q->gb, 8*(buf_size - is_codecframe_fmt) - q->frame->bits);

    /**
     * check for erasures/blanks on rates 1, 1/4 and 1/8
     */

    if(q->frame->rate != RATE_HALF && q->frame->data[QCELP_RSRVD_POS])
    {
        av_log(NULL, AV_LOG_ERROR, "Wrong data in reserved frame area:%d\n",
               q->frame->data[QCELP_RSRVD_POS]);
        is_ifq=1;
    }

    if(q->frame->rate == RATE_OCTAVE && first16==0xFFFF)
    {
        av_log(NULL, AV_LOG_ERROR,
        "Wrong frame data, rate 1/8 and first 16 bits are on\n");
        is_ifq=1;
    }

    /**
     * Preliminary decoding of frame's transmission codes
     */

    qcelp_decode_lspf(q->frame, qtzd_lspf);
    qcelp_decode_params(q->frame, g0, &cbseed, gain, index);

    /**
     * Check for badly received packets
     * TIA/EIA/IS-733 2.4.8.7.3
     */

    if(q->frame->rate != RATE_OCTAVE)
    {

        /* check for outbound LSP freqs and codebook gain params */
        if(q->frame->rate != RATE_QUARTER)
        {
            if(qtzd_lspf[9] <= .66 || qtzd_lspf[9] >= .985)
            {
                av_log(NULL, AV_LOG_ERROR,
                "Wrong frame data, 9th LSPF=%4f outside [.66,.985]\n", qtzd_lspf[9]);
                is_ifq=1; /* FIXME 'erase packet'==ifq? */
            }

            for(n=4; !is_ifq && n<10; n++)
            {
                if(FFABS(qtzd_lspf[n]-qtzd_lspf[n-4]) < .0931)
                {
                    av_log(NULL, AV_LOG_ERROR,
                    "Wrong frame data, outbound LSP freqs\n");
                    is_ifq=1;
                }
            }
        }else
        {
            if(qtzd_lspf[9] <= .70 || qtzd_lspf[9] >=  .97)
                is_ifq=1;

            for(n=3; !is_ifq && n<10; n++)
            {
                if(FFABS(qtzd_lspf[n]-qtzd_lspf[n-2]) < .08)
                    is_ifq=1;
            }
            /* codebook gain sanity check */
            /* FIXME This should be implemented into qcelp_decode_params() */
            for(n=0; !is_ifq && n<4; n++)
            {
                if(FFABS(g0[n+1]-g0[n]) > 40) is_ifq=1;
                /* FIXME: spec with typing errors here? */
                if(n<3 && FFABS(g0[n+2] - 2*g0[n+1] + g0[n]) > 48) is_ifq=1;
            }

        }
    }

    /**
     * decode loop glue code. WIP - mean it, WIP. :-)
     */

    if(!is_ifq)
    {
        qcelp_compute_svector(q->frame->rate, gain, index, cbseed, cdn_vector);
        if((is_ifq = qcelp_do_pitchfilter(q->frame, cdn_vector)))
        {
            av_log(NULL, AV_LOG_ERROR, "Error can't pitch cdn_vector[%d]\n",
            is_ifq);
            is_ifq=1;
        }
    }

    if(is_ifq)
    {
        /**
         * Insufficient frame quality (erasure) decoding
         */
        q->ifq_count++;

    }

    *data_size=buf_size;
    return *data_size;
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
