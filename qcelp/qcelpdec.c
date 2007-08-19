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
    int     bits;
} QCELPFrame;

typedef struct {
    GetBitContext gb;
    QCELPFrame    *frame;
    uint8_t       erasure_count;
    uint8_t       ifq_count;
    float         prev_lspf[10];
    float         pitchf_mem[144];
    float         pitchp_mem[144];
    float         hammsinc_table[8];
    int           frame_num;
} QCELPContext;

static int qcelp_decode_init(AVCodecContext *avctx);
static int qcelp_decode_frame(AVCodecContext *avctx, void *data,
           int *data_size, uint8_t *buf, int buf_size);
static int qcelp_decode_close(AVCodecContext *avctx);


static float qcelp_hammsinc(float i)
{
    return (sin(M_PI*i)/(M_PI*i))*(0.5+0.46*cos(M_PI*i/4));
}

static void qcelp_update_pitchf_mem(float *pitchf_mem, float last)
{
    float tmp[144];

    memcpy(tmp, pitchf_mem+1, 143*sizeof(float));
    pitchf_mem[143]=last;
    memcpy(pitchf_mem, tmp, 143*sizeof(float));
}

static int qcelp_decode_init(AVCodecContext *avctx)
{
    QCELPContext *q = (QCELPContext *) avctx->priv_data;
    int i;

    avctx->sample_rate = 8000;
    avctx->channels = 1;

    q->frame = av_mallocz(sizeof(QCELPFrame));


    if(q->frame == NULL)
        return -1;

    q->frame_num=0;

    memset(q->prev_lspf , 0, sizeof(q->prev_lspf ));
    memset(q->pitchf_mem, 0, sizeof(q->pitchf_mem));
    memset(q->pitchp_mem, 0, sizeof(q->pitchp_mem));

    /**
     * Fill hammsinc table
     */

    for(i=0; i<8; i++)
        q->hammsinc_table[i]=qcelp_hammsinc(i-3.5);

    return 0;
}

static int qcelp_decode_close(AVCodecContext *avctx)
{
    QCELPContext *q = avctx->priv_data;

    av_free(q->frame);

    return 0;
}

/**
 * Decodes the 10 quantized LSP frequencies from the LSPV/LSP
 * transsmision codes of any frame rate.
 *
 * For details see TIA/EIA/IS-733 2.4.3.2.6.2-2
 *
 * WIP: implement I_F_Q handling?
 */
void qcelp_decode_lspf(AVCodecContext *avctx, const QCELPFrame *frame,
     float *lspf)
{
    const uint8_t *lspv;
    int i;

    switch(frame->rate)
    {
        case RATE_FULL:
        case RATE_HALF:
        case RATE_QUARTER:
            lspv=frame->data+QCELP_LSPV0_POS;

            av_log(avctx, AV_LOG_DEBUG,
                   "-------- Decoded values from frame --------\n");
            av_log(avctx, AV_LOG_DEBUG, "[LSPV] %5d %5d %5d %5d %5d\n",
                   lspv[0],lspv[1],lspv[2],lspv[3],lspv[4]);

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
void qcelp_decode_params(AVCodecContext *avctx, const QCELPFrame *frame,
     int *g0, uint16_t *cbseed, float *gain, int *index)
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
                    predictor=av_clip(floor((g1[i-1]+g1[i-2]+g1[i-3])/3), 6,
                                      38);
                else
                    predictor=0;

                g1[i]=g0[i]+predictor;

                /* *
                 * FIXME
                 *
                 * This shouldn't but does happen quite a few
                 * times during decoding, my guess would be an
                 * error at predictor computation.
                 */

                if(g1[i]<0 || g1[i]>60)
                {
                    av_log(avctx, AV_LOG_WARNING,
                           "Gain Ga %d out of range for CBGAIN number %d\n",
                           g1[i], i);
                    g1[i]=av_clip(g1[i], 0, 60);
                }

                ga[i]=qcelp_g12ga[g1[i]];

                gain[i]=ga[i]*gs[i];
                index[i]=(gs[i] == 1)? cindex[i]:(cindex[i]-89) & 127;
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

    av_log(NULL, AV_LOG_ERROR,
           "[ GAIN ] %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n",
           gain[0],gain[1],gain[2],gain[3],gain[4],gain[5],gain[6],gain[7],
           gain[8],gain[9],gain[10],gain[11],gain[12],gain[13],gain[14],
           gain[15],gain[16]);
    av_log(NULL, AV_LOG_ERROR,
           "[ INDEX] %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
           index[0],index[1],index[2],index[3],index[4],index[5],index[6],
           index[7],index[8],index[9],index[10],index[11],index[12],index[13],
           index[14], index[15],index[16]);
    av_log(NULL, AV_LOG_ERROR, "[  CDN]");

    switch(rate)
    {
        case RATE_FULL:
             for(i=0; i<160; i++)
             {
                cdn_vector[i]=
                gain[i/10]*qcelp_fullrate_ccodebook[(i+1-index[i/10]) & 127];
                av_log(NULL, AV_LOG_ERROR, " %f", cdn_vector[i]);
             }
             av_log(NULL, AV_LOG_ERROR, "\n");
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
 * Computes energy of the subframeno-ith subvector, using equations
 * 2.4.8.3-2 and 2.4.3.8-3
 */
static float qcelp_compute_subframe_energy(const float *vector, int subframeno)
{
    int   i;
    float energy=0;

    vector+=40*subframeno;

    for(i=0; i<40; i++)
        energy+=vector[i]*vector[i];

    return energy;
}

static void qcelp_get_gain_scalefactors(const float *in, const float *out,
            float *scalefactors)
{
    int i;

    for(i=0; i<4; i++)
          scalefactors[i]=sqrt(qcelp_compute_subframe_energy(in , i)/
                               qcelp_compute_subframe_energy(out, i));
}

static void qcelp_apply_gain_ctrl(int do_iirf, const float *in, float *out)
{
    int i;
    float scalefactors[4];

    qcelp_get_gain_scalefactors(in, out, scalefactors);

    /* 2.4.8.6-6 */

    if(do_iirf)
    {
        scalefactors[0]*=0.0625;

        for(i=1;i<4;i++)
            scalefactors[i]=0.9375*scalefactors[i-1]+0.0625*scalefactors[i];
    }

    for(i=0; i<160; i++)
        out[i]=scalefactors[i/40]*out[i];

}

/**
 * pitch filters or pre-filters pv, returns 0 if everything goes
 * well, otherwise it returns the index of the failing-to-be-pitched
 * element and -1 if an invalid (140.5, 141.5, 142.5) lag is found or
 * an invalid operation mode is requested.
 *
 * This function implements both, the pitch filter and the pitch pre-filter
 * whose results gets stored in pv.
 *
 * For details see 2.4.5.2
 *
 * WIP (but should work)
 *
 * @param step Mode, 1 for pitch filter or 2 for pitch pre-filter
 *
 */
static int qcelp_do_pitchfilter(QCELPFrame *frame, float *pitch_mem, int step,
           float *pv, float *hammsinc_table)
{
    int     i, j, tmp;
    uint8_t *pgain, *plag, *pfrac;
    float   gain[4], lag[4];

    if(step != 1 && step != 2)
        return -1;

    switch(frame->rate)
    {
        case RATE_FULL:
        case RATE_HALF:

            pgain=frame->data+QCELP_PGAIN0_POS;
            plag =frame->data+QCELP_PLAG0_POS;
            pfrac=frame->data+QCELP_PFRAC0_POS;

            /**
             * Compute Gain & Lag
             */

            for(i=0; i<4; i++)
            {
                gain[i]=plag[i]? (pgain[i]+1)/4.0 : 0.0;

                if(step == 2) /* become pitch pre-filter */
                    gain[i]=0.5*FFMIN(gain[i],1.0);

                lag[i]  =plag[i]+16;

                if(pfrac[i])
                    lag[i]+=0.5;

                if(lag[i] == 140.5 || lag[i] == 141.5 || lag[i] == 142.5)
                    return -1;
            }

            /**
             * Apply filter
             */

            for(i=0; i<160; i++)
            {
                if(pfrac[i/40]) /* if is a fractional lag... */
                {
                    for(j=-4; j<4; j++)
                    {
                        tmp = i+j+0.5-lag[i/40];

                        if(tmp < 0)
                            pv[i]+=gain[i/40]*hammsinc_table[j+4]
                                   * pitch_mem[144+tmp];
                        else
                            pv[i]+=gain[i/40]*hammsinc_table[j+4]
                                   * pv [tmp];
                    }

                }else
                {
                    tmp=i-lag[i/40];

                    if(tmp < 0)
                        pv[i]+=gain[i/40]*pitch_mem[144+tmp];
                    else
                        pv[i]+=gain[i/40]*pv[i - lrintf(lag[i/40])];
                }

                qcelp_update_pitchf_mem(pitch_mem, pv[i]);
            }

            break;
        case RATE_QUARTER:
        case RATE_OCTAVE:
            break;
    }

    return 0;
}

/**
 * Computes interpolated lsp frequencies for a given rate & pitch subframe
 *
 * For details see 2.4.3.3.4
 *
 * @param rate Frame rate
 * @param prev_lspf Previous frame LSP freqs vector
 * @param curr_lspf Current frame LSP freqs vector
 * @param interpolated_lspf Float vector to put the resulting LSP freqs
 * @param frame_num Frame number in decoded stream
 *
 */
void qcelp_do_interpolate_lspf(qcelp_packet_rate rate, float *prev_lspf,
     float *curr_lspf, float *interpolated_lspf, int sample_num, int frame_num)
{
    int   i;
    float curr_weight, prev_weight;

    switch(rate)
    {
        case RATE_FULL:
        case RATE_HALF:
        case RATE_QUARTER:

                if(!frame_num)
                {
                    curr_weight=1.0;
                    prev_weight=0.0;
                }else
                {
                    switch(sample_num)
                    {
                        case 0:
                            curr_weight=0.25;
                            prev_weight=0.75;
                            break;
                        case 39:
                            curr_weight=0.5;
                            prev_weight=0.5;
                            break;
                        case 79:
                            curr_weight=0.75;
                            prev_weight=0.25;
                            break;
                        case 119:
                            curr_weight=1.0;
                            prev_weight=0;
                    }
                }

            for(i=0;i<10;i++)
                interpolated_lspf[i]=prev_weight*prev_lspf[i]+
                                     curr_weight*curr_lspf[i];
            break;
        case RATE_OCTAVE:

            curr_weight=0.625;
            prev_weight=0.375;

            for(i=0;i<10;i++)
                interpolated_lspf[i]=prev_weight*prev_lspf[i]+
                                     curr_weight*curr_lspf[i];
            break;
        case I_F_Q:
            memcpy(interpolated_lspf, prev_lspf, 10*sizeof(float));
    }
}

/**
 * Linear convolution of two vectors, with max resultant
 * vector dim = 12 -- just what we need. Result gets stored
 * in vector_a so it must have enough room to hold d1+d2-1.
 *
 * WIP this is a heavily suboptimal implementation
 *
 * @param d1 real dimension of vector_a prior convolution
 * @param d2 dimension of vector_b
 *
 */
static void qcelp_convolve(float *vector_a, const float *vector_b, int d1,
            int d2)
{
    float copy[12];
    int i,j;

    memcpy(copy, vector_a, sizeof(copy));

    for(i=0;i<(d1+d2-1);i++)
    {
        vector_a[i]=0.0;
        for(j=0;j<=i;j++)
            vector_a[i]+=(((i-j)>=d1 || (i-j)<0)? 0:copy[i-j])*(j>=d2? 0:vector_b[j]);
    }

}

/**
 * 2.4.3.3.5-1/2
 */
static void qcelp_lsp2poly(float *lspf, float *pa, float *qa)
{
    int i;
    float vector_a[12];
    float vector_b[3];
    int   a_limit[]={2,4,6,8,10};

    vector_b[0]=1;
    vector_b[2]=1;

    /**
     * Compute Pa coefs
     */

    vector_a[0]=1.0;
    vector_a[1]=1.0;

    for(i=0; i<5; i++)
    {
        vector_b[1]=-2*cos(M_PI*lspf[2*i]);
        qcelp_convolve(vector_a, vector_b, a_limit[i], 3);
    }

    for(i=0;i<5;i++)
        pa[i]=vector_a[i+1];

    /**
     * Compute Qa coefs
     */

    vector_a[0]= 1.0;
    vector_a[1]=-1.0;

    for(i=0; i<5; i++)
    {
        vector_b[1]=-2*cos(M_PI*lspf[2*i+1]);
        qcelp_convolve(vector_a, vector_b, a_limit[i], 3);
    }

    for(i=0;i<5;i++)
        qa[i]=vector_a[i+1];

}

/**
 * 2.4.3.3.5
 */
static void qcelp_lsp2lpc(AVCodecContext *avctx, float *lspf, float *lpc)
{
    float pa[5],qa[5];
    int   i;

    qcelp_lsp2poly(lspf, pa, qa);

    for(i=0; i< 5; i++)
            lpc[i]=-(pa[i]+qa[i])/2.0;
    for(i=5; i<10; i++)
            lpc[i]=-(pa[9-i]-qa[9-i])/2.0;

    av_log(avctx, AV_LOG_DEBUG, "[LPC ] %f %f %f %f %f %f %f %f %f %f\n",
           lpc[0], lpc[1], lpc[2], lpc[3], lpc[4],
           lpc[5], lpc[6], lpc[7], lpc[8], lpc[9]);
}

/**
 * 2.4.3.1
 *
 * This is the 10th order predictor error filter -- the reciprocal
 * of the formant synthesis filter.
 */
static float qcelp_prede_filter(float *lpc, float z)
{
    int   i;
    float tmp=0.0;

    for(i=0; i<10; i++)
       tmp+=lpc[i]*1.0/z;

    return(1.0-tmp);
}

/**
 * 2.4.8.6-2
 * Used in the adaptive postfilter at sample generation stage.
 */
static float qcelp_detilt(float z)
{
    if(z)
        return (1.0/(1.0 + 0.3 / z));
    else
        return 0;
}

static int qcelp_decode_frame(AVCodecContext *avctx, void *data,
           int *data_size, uint8_t *buf, int buf_size)
{
    QCELPContext *q    = avctx->priv_data;
    const QCELPBitmap *order = NULL;
    int16_t  *outbuffer = data, cbseed;
    int      i, n, is_ifq = 0, is_codecframe_fmt = 0;
    uint16_t first16 = 0;
    float    qtzd_lspf[10], gain[16], cdn_vector[160], ppf_vector[160], lpc[10];
    float    interpolated_lspf[10];
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
            is_codecframe_fmt=1;
        case 34:
            q->frame->rate = RATE_FULL;
            q->frame->bits = qcelp_bits_per_rate[RATE_FULL];
            order = QCELP_REFERENCE_FRAME + QCELP_FULLPKT_REFERENCE_POS;
            break;
        case 17:
            is_codecframe_fmt=1;
        case 16:
            q->frame->rate = RATE_HALF;
            q->frame->bits = qcelp_bits_per_rate[RATE_HALF];
            order = QCELP_REFERENCE_FRAME + QCELP_HALFPKT_REFERENCE_POS;
            break;
        case  8:
            is_codecframe_fmt=1;
        case  7:
            q->frame->rate = RATE_QUARTER;
            q->frame->bits = qcelp_bits_per_rate[RATE_QUARTER];
            order = QCELP_REFERENCE_FRAME + QCELP_4THRPKT_REFERENCE_POS;
            break;
        case  4:
            is_codecframe_fmt=1;
        case  3:
            q->frame->rate = RATE_OCTAVE;
            q->frame->bits = qcelp_bits_per_rate[RATE_OCTAVE];
            order = QCELP_REFERENCE_FRAME + QCELP_8THRPKT_REFERENCE_POS;
            break;
        case  1:
            is_codecframe_fmt=1;
        case  0:
            q->frame->rate = BLANK;
            q->frame->bits = 0;
            order = NULL;
            break;
        default:
            q->frame->rate = RATE_UNKNOWN;
            q->frame->bits = 0;
            av_log(avctx, AV_LOG_ERROR, "UNKNOWN PACKET RATE\n");
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
           av_log(avctx, AV_LOG_ERROR,
                  "Claimed rate and buffer size missmatch\n");
           is_ifq=1;
        }
    }

    av_log(avctx, AV_LOG_DEBUG, "Rate %d Size %d\n",
           q->frame->rate, q->frame->bits);

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

    /* skip padding byte if codec_frame_fmt */

    skip_bits(&q->gb, 8*(buf_size - is_codecframe_fmt) - q->frame->bits);

    /**
     * check for erasures/blanks on rates 1, 1/4 and 1/8
     */

    if(q->frame->rate != RATE_HALF && q->frame->data[QCELP_RSRVD_POS])
    {
        av_log(avctx, AV_LOG_ERROR, "Wrong data in reserved frame area:%d\n",
               q->frame->data[QCELP_RSRVD_POS]);
        is_ifq=1;
    }

    if(q->frame->rate == RATE_OCTAVE && first16==0xFFFF)
    {
        av_log(avctx, AV_LOG_ERROR,
               "Wrong frame data, rate 1/8 and first 16 bits are on\n");
        is_ifq=1;
    }

    /**
     * Preliminary decoding of frame's transmission codes
     */

    qcelp_decode_lspf(avctx, q->frame, qtzd_lspf);
    qcelp_decode_params(avctx, q->frame, g0, &cbseed, gain, index);

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
                av_log(avctx, AV_LOG_ERROR,
                       "IFQ: 9th LSPF=%4f outside [.66,.985]\n", qtzd_lspf[9]);
                is_ifq=1;
            }

            for(n=4; !is_ifq && n<10; n++)
            {
                if(FFABS(qtzd_lspf[n]-qtzd_lspf[n-4]) < .0931)
                {
                    av_log(avctx, AV_LOG_ERROR, "Wrong data, outbound LSPFs\n");
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

            /* FIXME This should be implemented into qcelp_decode_params() */

            for(n=0; !is_ifq && n<4; n++)
            {
                if(FFABS(g0[n+1]-g0[n]) > 40) is_ifq=1;
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

        av_log(avctx, AV_LOG_DEBUG, "\n-------- Pre pitch filters --------\n");
        av_log(avctx, AV_LOG_DEBUG, "[CDN_VECTOR]:\n");
        for(i=0; i<160; i++)
        {
            av_log(avctx, AV_LOG_DEBUG, " %f", cdn_vector[i]);
        }

        /* pitch filter */

        if((is_ifq = qcelp_do_pitchfilter(q->frame, q->pitchf_mem,
                                          1, cdn_vector, q->hammsinc_table)))
        {
            av_log(avctx, AV_LOG_ERROR,
                   "Error can't pitchfilter cdn_vector[%d]\n",
                   is_ifq);
            is_ifq=1;
        }

        memcpy(ppf_vector, cdn_vector, 160*sizeof(float));

        /* pitch pre-filter */

        if((is_ifq = qcelp_do_pitchfilter(q->frame, q->pitchp_mem,
                                          2, ppf_vector, q->hammsinc_table)))
        {
            av_log(avctx, AV_LOG_ERROR,
                   "Error can't pitch-prefilter ppf_vector[%d]\n",
                   is_ifq);
            is_ifq=1;
        }
    }

    av_log(avctx, AV_LOG_DEBUG, "-------- Post pitch filters --------\n");
    av_log(avctx, AV_LOG_DEBUG, "[CDN_VECTOR]:\n");
    for(i=0; i<160; i++)
    {
        av_log(avctx, AV_LOG_DEBUG, " %f", cdn_vector[i]);
    }

    av_log(avctx, AV_LOG_DEBUG, "[PPF_VECTOR]:\n");
    for(i=0; i<160; i++)
    {
        av_log(avctx, AV_LOG_DEBUG, " %f", ppf_vector[i]);
    }
    av_log(avctx, AV_LOG_DEBUG, "\n");

    /* pitch gain control */

    qcelp_apply_gain_ctrl(0, cdn_vector, ppf_vector);

    av_log(avctx, AV_LOG_DEBUG, "-------- Post Gain control --------\n");
    av_log(avctx, AV_LOG_DEBUG, "[CDN_VECTOR]:\n");
    for(i=0; i<160; i++)
    {
        av_log(avctx, AV_LOG_DEBUG, " %f", cdn_vector[i]);
    }
    av_log(avctx, AV_LOG_DEBUG, "\n");

    av_log(avctx, AV_LOG_DEBUG, "[PPF_VECTOR]:\n");
    for(i=0; i<160; i++)
    {
        av_log(avctx, AV_LOG_DEBUG, " %f", ppf_vector[i]);
    }
    av_log(avctx, AV_LOG_DEBUG, "\n");

    /* Apply formant synthesis filter over the pitch pre-filter output. */

    av_log(avctx, AV_LOG_DEBUG, "-------- Output --------\n");
    for(i=0; i<160; i++)
    {
        /* interpolate lsp freqs */

        if(i == 0 || i == 39  || i == 79 || i == 119)
        {
            qcelp_do_interpolate_lspf(q->frame->rate, q->prev_lspf, qtzd_lspf,
                                      interpolated_lspf, i, q->frame_num);
            qcelp_lsp2lpc(avctx, interpolated_lspf, lpc);
        }

        /* formant */

        ppf_vector[i]=1.0/qcelp_prede_filter(lpc, ppf_vector[i]);

        /* adaptive postfilter */

        ppf_vector[i]=qcelp_detilt(ppf_vector[i])*
                      qcelp_prede_filter(lpc, ppf_vector[i]/0.625)/
                      qcelp_prede_filter(lpc, ppf_vector[i]/0.775);

        av_log(avctx, AV_LOG_DEBUG, " %f/", ppf_vector[i]);

        /* output stage */

        outbuffer[i]=av_clip_int16(lrintf(4*ppf_vector[i]));
        av_log(avctx, AV_LOG_DEBUG, "%d", outbuffer[i]);

    }
    av_log(avctx, AV_LOG_DEBUG, "\n");

    if(is_ifq)
    {
        /**
         * Insufficient frame quality (erasure) decoding
         */

        q->ifq_count++;

    }

    /**
     * Copy current lspf freqs over to prev_lspf
     */

    memcpy(q->prev_lspf, qtzd_lspf, 10*sizeof(float));

    q->frame_num++;
    *data_size=160;

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
