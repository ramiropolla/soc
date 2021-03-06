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
 * @file
 * QCELP decoder
 * @author Reynaldo H. Verdejo Pinochet
 * @remark FFmpeg merging spearheaded by Kenan Gillet
 */

#include <math.h>
#include <stddef.h>

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"

#include "qcelpdata.h"

typedef struct
{
    GetBitContext     gb;
    qcelp_packet_rate rate;
    uint8_t           data[76];       /*!< data from a _parsed_ frame */
    uint8_t           erasure_count;
    uint8_t           ifq_count;
    float             prev_lspf[10];
    float             pitchf_mem[150];
    float             pitchp_mem[150];
    float             formant_mem[10];
    float             last_codebook_gain;
    int               frame_num;
    int               bits;
} QCELPContext;

static void qcelp_update_filter_mem(float *pitchf_mem, const float *last)
{
    memmove(pitchf_mem, pitchf_mem+40, 110*sizeof(float));
    memmove(pitchf_mem+110, last, 40*sizeof(float));
}

static int qcelp_decode_init(AVCodecContext *avctx)
{
    return 0;
}

/**
 * Decodes the 10 quantized LSP frequencies from the LSPV/LSP
 * transmission codes of any framerate.
 *
 * TIA/EIA/IS-733 2.4.3.2.6.2-2
 */
static void qcelp_decode_lspf(const QCELPContext *q, float *lspf, int *is_ifq)
{
    const uint8_t *lspv;
    int i;

    if(q->rate == RATE_OCTAVE)
    {
        lspv=q->data+QCELP_LSP0_POS;
        for(i=0; i<10; i++)
        {
            lspf[i]=lspv[i]? 0.02:-0.02; /* 2.4.3.3.1-1 */
        }
    }else
    {
        lspv=q->data+QCELP_LSPV0_POS;

        lspf[0]=        qcelp_lspvq1[lspv[0]].x / 10000.0;
        lspf[1]=lspf[0]+qcelp_lspvq1[lspv[0]].y / 10000.0;
        lspf[2]=lspf[1]+qcelp_lspvq2[lspv[1]].x / 10000.0;
        lspf[3]=lspf[2]+qcelp_lspvq2[lspv[1]].y / 10000.0;
        lspf[4]=lspf[3]+qcelp_lspvq3[lspv[2]].x / 10000.0;
        lspf[5]=lspf[4]+qcelp_lspvq3[lspv[2]].y / 10000.0;
        lspf[6]=lspf[5]+qcelp_lspvq4[lspv[3]].x / 10000.0;
        lspf[7]=lspf[6]+qcelp_lspvq4[lspv[3]].y / 10000.0;
        lspf[8]=lspf[7]+qcelp_lspvq5[lspv[4]].x / 10000.0;
        lspf[9]=lspf[8]+qcelp_lspvq5[lspv[4]].y / 10000.0;

        // Check for badly received packets TIA/EIA/IS-733 2.4.8.7.3

        if(!*is_ifq)
        {
            if(q->rate != RATE_QUARTER)
            {
                if(lspf[9] <= .66 || lspf[9] >= .985) *is_ifq=1;

                for(i=4; !*is_ifq && i<10; i++)
                    if(FFABS(lspf[i] - lspf[i-4]) < .0931) *is_ifq=1;
            }else
            {
                if(lspf[9] <= .70 || lspf[9] >=  .97) *is_ifq=1;

                for(i=3; !*is_ifq && i<10; i++)
                    if(FFABS(lspf[i] - lspf[i-2]) < .08) *is_ifq=1;
            }

        }
    }
}

/**
 * Converts codebook transmission codes to GAIN and INDEX
 * (and cbseed for rate 1/4).
 *
 * TIA/EIA/IS-733 2.4.6.2
 */
static void qcelp_decode_params(AVCodecContext *avctx, uint16_t *cbseed,
                                float *gain, int *index, int *is_ifq)
{
    int                 i, gs[16], g0[16], g1[16], predictor;
    const uint8_t       *cbgain, *cbsign, *cindex, *data;
    float               ga[16], gain_memory;
    QCELPContext  *q = avctx->priv_data;

    cbsign=q->data+QCELP_CBSIGN0_POS;
    cbgain=q->data+QCELP_CBGAIN0_POS;
    cindex=q->data+QCELP_CINDEX0_POS;

    switch(q->rate)
    {
        case RATE_FULL:
        case RATE_HALF:
            for(i=0; i<16; i++)
            {
                if(q->rate == RATE_HALF && i>=4) break;

                gs[i]=cbsign[i]? -1:1;
                g0[i]=4*cbgain[i];

                /*
                 * TIA/EIA/IS-733 Spec has errors on the predictor determination
                 * formula at equation 2.4.6.1-4 -- The predictor there needs 6
                 * to be subtracted from it to give RI-compliant results. The
                 * problem is it ignores the fact that codebook subframes 4, 8,
                 * 12 and 16 on a FULL_RATE frame use a different quantizer
                 * table.
                 */

                if(q->rate == RATE_FULL && i > 0 && !((i+1) & 3))
                    predictor=av_clip(floor((g1[i-1]+g1[i-2]+g1[i-3])/3.0), 6,
                              38)-6;
                else
                    predictor=0;

                g1[i]=g0[i]+predictor;

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

            q->last_codebook_gain=gain[i-1];

            break;
        case RATE_QUARTER:
            for(i=0; i<5; i++)
            {
                g0[i]=g1[i]=4*cbgain[i];

                if(!*is_ifq)
                {
                    if(i>0 && FFABS(g0[i] - g0[i-1]) > 40) *is_ifq=1;
                    if(i<3 && FFABS(g0[i+2] - 2*g0[i+1] + g0[i]) > 48) *is_ifq=1;
                }
                ga[i]=qcelp_g12ga[g1[i]];
            }

            q->last_codebook_gain=gain[4];

            /*
             * 5->8 interpolation to 'Provide smoothing of the energy
             * of the unvoiced excitation' TIA/EIA/IS-733 2.4.6.2
             */

            gain[0]=    ga[0];
            gain[1]=0.6*ga[0]+0.4*ga[1];
            gain[2]=    ga[1];
            gain[3]=0.2*ga[1]+0.8*ga[2];
            gain[4]=0.8*ga[2]+0.2*ga[3];
            gain[5]=    ga[3];
            gain[7]=0.4*ga[3]+0.6*ga[4];
            gain[7]=    ga[4];

            // Build random* seed needed to make Cdn

            data=q->data;
            *cbseed=(0x0003 & data[QCELP_LSPV0_POS+4])<<14 |
                    (0x003C & data[QCELP_LSPV0_POS+3])<< 8 |
                    (0x0060 & data[QCELP_LSPV0_POS+2])<< 1 |
                    (0x0007 & data[QCELP_LSPV0_POS+1])<< 3 |
                    (0x0038 & data[QCELP_LSPV0_POS  ])>> 3 ;
            break;
        case RATE_OCTAVE:
            switch(cbgain[0])
            {
                case 0: ga[0]=-4; break;
                case 1: ga[0]=-2; break;
                case 2: ga[0]= 0; break;
                case 3: ga[0]= 2; break;
            }

            gain_memory=FFABS(q->last_codebook_gain);

            gain[0]=0.5*gain_memory + 0.5*ga[0];

            q->last_codebook_gain=gain[0];

            /*
             * This interpolation is done to produce smoother
             * background noise. TIA/EIA/IS-733 2.4.6.2.3
             */

            for(i=0; i<8; i++)
                gain[i]=(0.875-0.125*i)*gain_memory+(0.125+0.125*i)*gain[0];
    }
}

/**
 * Computes the scaled codebook vector Cdn From INDEX and GAIN
 * for all rates.
 *
 * @param rate rate of the current frame/packet
 * @param gain array holding the 4 pitch subframe gain values
 * @param index array holding the 4 pitch subframe index values
 * @param cbseed seed needed for scaled codebook vector generation on rates
 * other than RATE_FULL or RATE_HALF
 * @param cnd_vector array for the generated scaled codebook vector
 */
static int qcelp_compute_svector(qcelp_packet_rate rate, const float *gain,
                                 const int *index, uint16_t cbseed,
                                 float *cdn_vector)
{
    int      i,j;
    uint16_t new_cbseed;
    float    rnd[160];


   /*
    * The specification misses some information here.
    *
    * TIA/EIA/IS-733 has an omission on the codebook index determination
    * formula for RATE_FULL and RATE_HALF frames at section 2.4.8.1.1. It says
    * you have to subtract the decoded index parameter from the given scaled
    * codebook vector index 'n' to get the desired circular codebook index, but
    * it does not mention that you have to clamp 'n' to [0-9] in order to get
    * RI-compliant results.
    *
    * The reason for this mistake seems to be the fact they forget to mention
    * you have to do these calculations per codebook subframe and adjust given
    * equation values accordingly.
    */

    j=0;

    switch(rate)
    {
        case RATE_FULL:

            for(i=0; i<160; i++)
            {
                cdn_vector[i]=
                gain[i/10]*qcelp_fullrate_ccodebook[(j-index[i/10]) & 127];

                j=j<9? j+1:0;  // See comment above
            }
            break;
        case RATE_HALF:

            for(i=0; i<160; i++)
            {
                cdn_vector[i]=
                gain[i/40]*qcelp_halfrate_ccodebook[(j-index[i/40]) & 127];

                j=j<9? j+1:0;  // See comment above
            }
            break;
        case RATE_QUARTER:
            for(i=0; i<160; i++)
            {
                new_cbseed=(521*cbseed+259) & 65535;
                cbseed=rnd[i]=
                QCELP_SQRT1887*(((new_cbseed+32768) & 65535)-32768)/32768.0;

                // FIR filter

                cdn_vector[i]=qcelp_rnd_fir_coefs[1]*rnd[i];
                for(j=1; j<22 && !(i-j+1); j++)
                {
                    cdn_vector[i]+=qcelp_rnd_fir_coefs[j]*rnd[i-j];
                }

                // final scaling

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
 * Computes energy of the subframeno-ith subvector. These values are
 * used to generate the scalefactors for the gain control stages.
 *
 * @param vector vector from which to measure the subframe's energy
 * @param subframeno size 40 subframe number that should be measured
 *
 * TIA/EIA/IS-733 2.4.8.3-2/3
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

/**
 * Computes scalefactors needed to gain-control 'in' and 'out' vectors.
 *
 * @param scalefactors array to place the resulting four scalefactors
 */
static void qcelp_get_gain_scalefactors(const float *in, const float *out,
                                        float *scalefactors)
{
    int i;

    for(i=0; i<4; i++)
          scalefactors[i]=sqrt(qcelp_compute_subframe_energy(in , i)/
                               qcelp_compute_subframe_energy(out, i));
}

/**
 * Generic gain control stage to implement TIA/EIA/IS-733 2.4.8.6-6 and
 * FIXME_MISSINGSPECSECTION.
 *
 * @param do_iirf whether or not to apply hardcoded coefficient infinite
 * impulse response filter
 * @param in vector to control gain off
 * @param out gain-controlled output vector
 */
static void qcelp_apply_gain_ctrl(int do_iirf, const float *in, float *out)
{
    int i;
    float scalefactors[4];

    qcelp_get_gain_scalefactors(in, out, scalefactors);

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
 * Pitch filters or pre-filters pv, returns 0 if everything goes
 * well, otherwise it returns the index of the failing-to-be-pitched
 * element or -1 if an invalid (140.5, 141.5, 142.5, 243.5) lag is found.
 *
 * This function implements both, the pitch filter and the pitch pre-filter
 * whose results gets stored in pv.
 *
 * TIA/EIA/IS-733 2.4.5.2
 *
 * @param step Mode, 1 for pitch filter or 2 for pitch pre-filter
 */
static int qcelp_do_pitchfilter(QCELPContext *q, int step, float *pv)
{
    int     i, j, k, tmp;
    uint8_t *pgain, *plag, *pfrac;
    float   gain[4], lag[4], hamm_tmp;

    assert(step == 1 || step == 2);

    switch(q->rate)
    {
        case RATE_FULL:
        case RATE_HALF:

            pgain=q->data+QCELP_PGAIN0_POS;
            plag =q->data+QCELP_PLAG0_POS;
            pfrac=q->data+QCELP_PFRAC0_POS;

            // Compute gain & lag for the whole frame

            for(i=0; i<4; i++)
            {
                gain[i]=plag[i]? (pgain[i]+1)/4.0 : 0.0;

                if(step == 2) // Become pitch pre-filter
                    gain[i]=0.5*FFMIN(gain[i],1.0);

                lag[i]=plag[i]+16;

                if(pfrac[i])
                {
                    if(lag[i] >= 140) return -1; // WIP: hook to IFQ decoding
                    lag[i]+=0.5;
                }

            }

            /*
             * Apply filter.
             *
             * TIA/EIA/IS-733 2.4.5.2-2/3 equations are not clear enough but
             * we know this filter has to be applied in pitch-subframe steps.
             */

            k=0;
            for(i=0; i<160; i++)
            {
                if(pfrac[i/40]) // If is a fractional lag...
                {
                    hamm_tmp=0.0;

                    for(j=-4; j<4; j++)
                    {
                        tmp = k+j+0.5-lag[i/40];

                        if(tmp < 0)
                            hamm_tmp+=qcelp_hammsinc_table[j+4]
                                   * q->pitchf_mem[150+tmp];
                        else
                            hamm_tmp+=qcelp_hammsinc_table[j+4]
                                   * pv [tmp];
                    }

                    pv[i]+=gain[i/40]*hamm_tmp;

                }else
                {
                    tmp=k-lag[i/40];

                    if(tmp < 0)
                        pv[i]+=lrintf(gain[i/40]*q->pitchf_mem[150+tmp]);
                    else
                        pv[i]+=lrintf(gain[i/40]*pv[i - lrintf(lag[i/40])]);
                }

                // Done with the pitch subframe -- update filter memory.

                if(k==39)
                {
                    qcelp_update_filter_mem(q->pitchf_mem, &pv[i-k]);
                }

                k=(k<39)? k+1:0;
            }

            break;
        case RATE_QUARTER:
        case RATE_OCTAVE:
            break;
    }

    return 0;
}

/**
 * Computes interpolated LSP frequencies for a given rate & pitch subframe.
 *
 * TIA/EIA/IS-733 2.4.3.3.4
 *
 * @param rate framerate
 * @param prev_lspf LSP frequencies vector of the previous frame
 * @param curr_lspf LSP frequencies vector of the  current frame
 * @param interpolated_lspf float vector for the resulting LSP frequencies
 * @param frame_num frame number in decoded stream
 */
void qcelp_do_interpolate_lspf(qcelp_packet_rate rate, float *prev_lspf,
                               float *curr_lspf, float *interpolated_lspf,
                               int sample_num, int frame_num)
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
                        case 40:
                            curr_weight=0.5;
                            prev_weight=0.5;
                            break;
                        case 80:
                            curr_weight=0.75;
                            prev_weight=0.25;
                            break;
                        case 120:
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
 * vector dim = 12 -- just what we need. The result gets stored
 * in v1 so it must have enough room to hold d1 + d2 - 1.
 *
 * WIP This is a heavily suboptimal implementation.
 *
 * @param d1 real dimension of v1 prior to convolution
 * @param d2 dimension of v2
 */
static void qcelp_convolve(float *v1, const float *v2, int d1, int d2)
{
    float copy[12];
    int   i,j,dim;

    memcpy(copy, v1, sizeof(copy));
    dim=d1+d2-1;

    for(i=0; i<dim; i++)
    {
        v1[i]=0.0;
        for(j=0;j<=i;j++)
            v1[i]+=(((i-j)>=d1 || (i-j)<0)? 0:copy[i-j])*(j>=d2? 0:v2[j]);
    }

}

/**
 * Computes the Pa and Qa coefficients needed for LSP to LPC conversion.
 *
 * TIA/EIA/IS-733 2.4.3.3.5-1/2
 */
static void qcelp_lsp2poly(float *lspf, float *pa, float *qa)
{
    int i;
    float v1[12];
    float v2[3];
    int   limit[]={2,4,6,8,10};

    v2[0]=1;
    v2[2]=1;

    // Compute Pa coefficients

    v1[0]=1.0;
    v1[1]=1.0;

    for(i=0; i<5; i++)
    {
        v2[1]=-2*cos(M_PI*lspf[2*i]);
        qcelp_convolve(v1, v2, limit[i], 3);
    }

    for(i=0;i<5;i++)
        pa[i]=v1[i+1];

    // Compute Qa coefficients

    v1[0]= 1.0;
    v1[1]=-1.0;

    for(i=0; i<5; i++)
    {
        v2[1]=-2*cos(M_PI*lspf[2*i+1]);
        qcelp_convolve(v1, v2, limit[i], 3);
    }

    for(i=0;i<5;i++)
        qa[i]=v1[i+1];

}

/**
 * Reconstructs LPC coefficients from the line spectral pairs frequencies.
 *
 * TIA/EIA/IS-733 2.4.3.3.5
 */
static void qcelp_lsp2lpc(float *lspf, float *lpc)
{
    float pa[5],qa[5];
    int   i;

    qcelp_lsp2poly(lspf, pa, qa);

    for(i=0; i< 5; i++)
            lpc[i]=-(pa[i]+qa[i])/2.0;
    for(i=5; i<10; i++)
            lpc[i]=-(pa[9-i]-qa[9-i])/2.0;
}

/**
 * Formant synthesis filter
 *
 * TIA/EIA/IS-733 2.4.3.1 (NOOOOT)
 */
static void qcelp_do_formant(float *in, float *out, float *lpc_coefs,
                             float *memory)
{
    float tmp[50];
    int i,j;

    // Copy over previous ten generated samples

    memcpy(tmp, memory, 10*sizeof(float));
    memcpy(tmp+10, in, 40*sizeof(float));

    for(i=10;i<50;i++)
    {
        for(j=1;j<11;j++)
        {
            tmp[i]+=tmp[i-j]*lpc_coefs[j-1];
        }
    }

    // Update memory for next pitch subframe

    memcpy(memory, tmp+40, 10*sizeof(float));

    // Write filtered samples to *out

    memcpy(out, tmp+10, 40*sizeof(float));
}

/**
 * Detilt used in the adaptive postfilter after the formant synthesis
 * filter.
 *
 * TIA/EIA/IS-733 2.4.8.6-2
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
    int16_t  *outbuffer = data;
    int      i, n, is_ifq = 0, is_codecframe_fmt = 0;
    uint16_t first16 = 0, cbseed = 0;
    float    qtzd_lspf[10], gain[16], cdn_vector[160], ppf_vector[160], lpc[10];
    float    interpolated_lspf[10];
    int      index[16];
    uint8_t  claimed_rate;

    init_get_bits(&q->gb, buf, buf_size*8);

    /*
     * Figure out framerate by its size, set up a few utility variables
     * and point 'order' to the rate's reference _slice_ inside the
     * big REFERENCE_FRAME array.
     */

    switch(buf_size)
    {
        case 35:
            is_codecframe_fmt=1;
        case 34:
            q->rate = RATE_FULL;
            q->bits = qcelp_bits_per_rate[RATE_FULL];
            order = QCELP_REFERENCE_FRAME + QCELP_FULLPKT_REFERENCE_POS;
            break;
        case 17:
            is_codecframe_fmt=1;
        case 16:
            q->rate = RATE_HALF;
            q->bits = qcelp_bits_per_rate[RATE_HALF];
            order = QCELP_REFERENCE_FRAME + QCELP_HALFPKT_REFERENCE_POS;
            break;
        case  8:
            is_codecframe_fmt=1;
        case  7:
            q->rate = RATE_QUARTER;
            q->bits = qcelp_bits_per_rate[RATE_QUARTER];
            order = QCELP_REFERENCE_FRAME + QCELP_4THRPKT_REFERENCE_POS;
            break;
        case  4:
            is_codecframe_fmt=1;
        case  3:
            q->rate = RATE_OCTAVE;
            q->bits = qcelp_bits_per_rate[RATE_OCTAVE];
            order = QCELP_REFERENCE_FRAME + QCELP_8THRPKT_REFERENCE_POS;

            /*
             * We assume it is IFQ unless we find at least one '0'
             * in the first 16 bits, this check is performed as part of
             * the reordering loop that follows.
             */

            is_ifq = 1;
            break;
        case  1:
            is_codecframe_fmt=1;
        case  0:
            q->rate = BLANK;
            q->bits = 0;
            order = NULL;
            break;
        default:
            av_log(avctx, AV_LOG_ERROR, "Error decoding frame"
                   " -- Unknown framerate, unsupported size: %d\n",
                   buf_size);
            return -1;
    }

    if(is_codecframe_fmt)
    {
        claimed_rate=get_bits(&q->gb, 8);

        if((claimed_rate ==  0 && q->rate != BLANK       ) ||
           (claimed_rate ==  1 && q->rate != RATE_OCTAVE ) ||
           (claimed_rate ==  2 && q->rate != RATE_QUARTER) ||
           (claimed_rate ==  3 && q->rate != RATE_HALF   ) ||
           (claimed_rate ==  4 && q->rate != RATE_FULL   ))
        {
           av_log(avctx, AV_LOG_WARNING,
                  "Claimed rate and buffer size mismatch\n");
        }
    }

    // Data reordering loop

    memset(q->data, 0, 76);
    for(n=0; n < q->bits; n++)
    {
        q->data[ order[n].index ] |=
        get_bits1(&q->gb)<<order[n].bitpos;

        if(q->rate == RATE_OCTAVE)
        {
            if(n>3)  // Random seed
                cbseed  |= (uint16_t)q->data[ order[n].index ]<<(n-4);
            if(n<16 && is_ifq && !q->data[ order[n].index ]) is_ifq = 0;
        }

    }

    // Check for erasures/blanks on rates 1, 1/4 and 1/8

    if(q->rate != RATE_HALF && q->data[QCELP_RSRVD_POS])
    {
        av_log(avctx, AV_LOG_ERROR, "Wrong data in reserved frame area:%d\n",
               q->data[QCELP_RSRVD_POS]);
        is_ifq=1;
    }

    if(q->rate == RATE_OCTAVE && first16==0xFFFF)
    {
        av_log(avctx, AV_LOG_ERROR,
               "Wrong frame data, rate 1/8 and first 16 bits are on\n");
        is_ifq=1;
    }

    // Preliminary decoding of frame's transmission codes

    qcelp_decode_lspf(q, qtzd_lspf, &is_ifq);
    qcelp_decode_params(avctx, &cbseed, gain, index, &is_ifq);

    // Decode loop glue code. WIP - mean it, WIP. :-)

    if(!is_ifq)
    {
        qcelp_compute_svector(q->rate, gain, index, cbseed, cdn_vector);

        // Pitch filter

        if((is_ifq = qcelp_do_pitchfilter(q, 1, cdn_vector)))
        {
            av_log(avctx, AV_LOG_WARNING,
                   "Error can't pitchfilter cdn_vector[%d]\n", is_ifq);
            is_ifq=1;
        }

        memcpy(ppf_vector, cdn_vector, 160*sizeof(float));

        /*
         * Pitch pre-filter
         *
         * Making this runtime-selectable might be a good speed-wise compromise
         * at the cost of a small degradation in perceived output quality.
         */

        if((is_ifq = qcelp_do_pitchfilter(q, 2, ppf_vector)))
        {
            av_log(avctx, AV_LOG_WARNING,
                   "Error: Cannot pitch-prefilter ppf_vector[%d]\n", is_ifq);
            is_ifq=1;
        }
    }

    // Pitch gain control

    qcelp_apply_gain_ctrl(0, cdn_vector, ppf_vector);

    // Interpolate LSP frequencies and apply formant synthesis filter.

    for(i=0; i<4; i++)
    {
        qcelp_do_interpolate_lspf(q->rate, q->prev_lspf, qtzd_lspf,
                                  interpolated_lspf, i*40, q->frame_num);

        qcelp_lsp2lpc(interpolated_lspf, lpc);

        qcelp_do_formant(ppf_vector+i*40, cdn_vector+i*40, lpc, q->formant_mem);

        // WIP Adaptive postfilter should be here
    }

    // WIP Final gain control stage should be here

    for(i=0; i<160; i++)
    {
        outbuffer[i]=av_clip_int16(lrintf(4*cdn_vector[i]));
    }

    if(is_ifq)
    {
        av_log(avctx, AV_LOG_WARNING, "IFQ Frame %d\n",
               q->frame_num);
        q->ifq_count++;

    }

    // Copy current lspf frequencies over to prev_lspf

    memcpy(q->prev_lspf, qtzd_lspf, sizeof(q->prev_lspf));

    q->frame_num++;
    *data_size=160*2;

    return *data_size;
}

AVCodec qcelp_decoder =
{
    .name   = "qcelp",
    .type   = CODEC_TYPE_AUDIO,
    .id     = CODEC_ID_QCELP,
    .init   = qcelp_decode_init,
    .decode = qcelp_decode_frame,
    .priv_data_size = sizeof(QCELPContext),
};
