/*
 * AAC encoder psychoacoustic model
 * Copyright (C) 2008 Konstantin Shishkov
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
 * @file aacpsy.c
 * AAC encoder psychoacoustic model
 */

#include "avcodec.h"
#include "dsputil.h"
#include "aacpsy.h"

//borrowed from aac.c
static float pow2sf_tab[340];


#define SCALE_ONE_POS   140
#define SCALE_MAX_POS   255
#define SCALE_MAX_DIFF   60


/**
 * Convert coefficients to integers.
 * @return sum of coefficients
 * @see 3GPP TS26.403 5.6.2
 */
static inline int convert_coeffs(float *in, int *out, int size, int scale_idx)
{
    int i, sign, sum = 0;
    for(i = 0; i < size; i++){
        sign = in[i] > 0.0;
        out[i] = (int)(pow(FFABS(in[i]) * pow2sf_tab[200 - scale_idx + SCALE_ONE_POS], 0.75) + 0.4054);
        if(out[i] > 8191) out[i] = 8191;
        sum += out[i];
        if(sign) out[i] = -out[i];
    }
    return sum;
}

static inline float unquant(int q, int scale_idx){
    return (FFABS(q) * cbrt(q*1.0)) * pow2sf_tab[200 + scale_idx - SCALE_ONE_POS];
}
static inline float calc_distortion(float *c, int size, int scale_idx)
{
    int i;
    int q;
    float coef, unquant, sum = 0.0f;
    for(i = 0; i < size; i++){
        coef = FFABS(c[i]);
        q = (int)(pow(FFABS(coef) * pow2sf_tab[200 - scale_idx + SCALE_ONE_POS], 0.75) + 0.4054);
        q = av_clip(q, 0, 8191);
        unquant = (q * cbrt(q)) * pow2sf_tab[200 + scale_idx - SCALE_ONE_POS];
        sum += (coef - unquant) * (coef - unquant);
    }
    return sum;
}

/**
 * Produce integer coefficients from scalefactors provided by model.
 */
static void psy_create_output(AACPsyContext *apc, ChannelElement *cpe, int search_pulses)
{
    int i, w, w2, g, ch;
    int start, sum, maxsfb, cmaxsfb;
    int pulses, poff[4], pamp[4];

    for(ch = 0; ch < apc->avctx->channels; ch++){
        start = 0;
        maxsfb = 0;
        cpe->ch[ch].pulse.present = 0;
        for(w = 0; w < cpe->ch[ch].ics.num_windows; w++){
            for(g = 0; g < cpe->ch[ch].ics.num_swb; g++){
                sum = 0;
                //apply M/S
                if(!ch && cpe->ms.mask[w][g]){
                    for(i = 0; i < cpe->ch[ch].ics.swb_sizes[g]; i++){
                        cpe->ch[0].coeffs[start+i] = (cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i]) / 2.0;
                        cpe->ch[1].coeffs[start+i] =  cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                    }
                }
                if(!cpe->ch[ch].zeroes[w][g])
                    sum = convert_coeffs(cpe->ch[ch].coeffs + start, cpe->ch[ch].icoefs + start, cpe->ch[ch].ics.swb_sizes[g], cpe->ch[ch].sf_idx[w][g]);
                cpe->ch[ch].zeroes[w][g] = !sum;
                //try finding pulses
                if(search_pulses && cpe->ch[ch].ics.num_windows == 1 && !cpe->ch[ch].pulse.present){
                    pulses = 0;
                    memset(poff,0,sizeof(poff));
                    memset(pamp,0,sizeof(pamp));
                    for(i = 0; i < cpe->ch[ch].ics.swb_sizes[g]; i++){
                        if(pulses > 4 || (pulses && i > cpe->ch[ch].pulse.offset[pulses-1] - 31)) break;
                        if(FFABS(cpe->ch[ch].icoefs[start+i]) > 4 && pulses < 4){
                            poff[pulses] = i;
                            pamp[pulses] = FFMIN(FFABS(cpe->ch[ch].icoefs[start+i]) - 1, 15);
                            pulses++;
                        }
                    }
                    if(pulses){
                        cpe->ch[ch].pulse.present = 1;
                        cpe->ch[ch].pulse.start = g;
                        cpe->ch[ch].pulse.num_pulse = pulses;
                        for(i = 0; i < pulses; i++){
                            cpe->ch[ch].pulse.amp[i] = pamp[i];
                            cpe->ch[ch].pulse.offset[i] = i ? poff[i] - poff[i-1] : poff[0];

                            if(cpe->ch[ch].icoefs[start+poff[i]] > 0)
                                cpe->ch[ch].icoefs[start+poff[i]] -= pamp[i];
                            else
                                cpe->ch[ch].icoefs[start+poff[i]] += pamp[i];
                        }
                    }
                }
                start += cpe->ch[ch].ics.swb_sizes[g];
            }
            for(cmaxsfb = cpe->ch[ch].ics.num_swb; cmaxsfb > 0 && cpe->ch[ch].zeroes[0][cmaxsfb-1]; cmaxsfb--);
            maxsfb = FFMAX(maxsfb, cmaxsfb);
        }
        cpe->ch[ch].ics.max_sfb = maxsfb;

        //adjust zero bands for window groups
        for(w = 0; w < cpe->ch[ch].ics.num_windows; w++){
            if(cpe->ch[ch].ics.group_len[w]) continue;
            for(g = 0; g < cpe->ch[ch].ics.max_sfb; g++){
                i = 1;
                w2 = w;
                do{
                    if(!cpe->ch[ch].zeroes[w2][g]){
                        i = 0;
                        break;
                    }
                    w2++;
                }while(w2 < cpe->ch[ch].ics.num_windows && cpe->ch[ch].ics.group_len[w2]);
                cpe->ch[ch].zeroes[w][g] = i;
            }
        }
    }

    if(apc->avctx->channels > 1 && cpe->common_window){
        int msc = 0;
        cpe->ch[0].ics.max_sfb = FFMAX(cpe->ch[0].ics.max_sfb, cpe->ch[1].ics.max_sfb);
        cpe->ch[1].ics.max_sfb = cpe->ch[0].ics.max_sfb;
        for(w = 0; w < cpe->ch[0].ics.num_windows; w++)
            for(i = 0; i < cpe->ch[0].ics.max_sfb; i++)
                if(cpe->ms.mask[w][i]) msc++;
        if(msc == 0 || cpe->ch[0].ics.max_sfb == 0) cpe->ms.present = 0;
        else cpe->ms.present = msc < cpe->ch[0].ics.max_sfb ? 1 : 2;
    }
}

static void psy_null_window(AACPsyContext *apc, int16_t *audio, int channel, ChannelElement *cpe)
{
    int ch;

    for(ch = 0; ch < apc->avctx->channels; ch++){
        cpe->ch[ch].ics.window_sequence = ONLY_LONG_SEQUENCE;
        cpe->ch[ch].ics.use_kb_window[0] = 1;
        cpe->ch[ch].ics.num_windows = 1;
        cpe->ch[ch].ics.swb_sizes = apc->bands1024;
        cpe->ch[ch].ics.num_swb = apc->num_bands1024;
        cpe->ch[ch].ics.group_len[0] = 0;
    }
    cpe->common_window = cpe->ch[0].ics.use_kb_window[0] == cpe->ch[1].ics.use_kb_window[0];
}

static void psy_null_process(AACPsyContext *apc, int16_t *audio, int channel, ChannelElement *cpe)
{
    int start;
    int ch, g, i;
    int minscale;

    for(ch = 0; ch < apc->avctx->channels; ch++){
        start = 0;
        for(g = 0; g < apc->num_bands1024; g++){
            float energy = 0.0f, ffac = 0.0f, thr, dist;

            for(i = 0; i < apc->bands1024[g]; i++){
                energy += cpe->ch[ch].coeffs[start+i]*cpe->ch[ch].coeffs[start+i];
                ffac += sqrt(FFABS(cpe->ch[ch].coeffs[start+i]));
            }
            thr = energy * 0.001258925f;
            cpe->ch[ch].sf_idx[ch][g] = 136;
            cpe->ch[ch].zeroes[ch][g] = (energy == 0.0);
            if(cpe->ch[ch].zeroes[ch][g]) continue;
            minscale = (int)(2.66667 * (log2(6.75*thr) - log2(ffac)));
            cpe->ch[ch].sf_idx[ch][g] = SCALE_ONE_POS - minscale;
            while(cpe->ch[ch].sf_idx[ch][g] > 3){
                dist = calc_distortion(cpe->ch[ch].coeffs + start, apc->bands1024[g], cpe->ch[ch].sf_idx[ch][g]);
                if(dist < thr) break;
                cpe->ch[ch].sf_idx[ch][g] -= 3;
            }
        }
    }
    for(ch = 0; ch < apc->avctx->channels; ch++){
        minscale = 255;
        for(g = 0; g < apc->num_bands1024; g++)
            if(!cpe->ch[ch].zeroes[0][g])
                minscale = FFMIN(minscale, cpe->ch[ch].sf_idx[0][g]);
        cpe->ch[ch].gain = minscale;
        for(g = 0; g < apc->num_bands1024; g++)
            if(!cpe->ch[ch].zeroes[0][g])
                cpe->ch[ch].sf_idx[0][g] = FFMIN(minscale + SCALE_MAX_DIFF, cpe->ch[ch].sf_idx[0][g]);
    }
    psy_create_output(apc, cpe, 1);
}

static void psy_null8_window(AACPsyContext *apc, int16_t *audio, int channel, ChannelElement *cpe)
{
    int ch, i;

    for(ch = 0; ch < apc->avctx->channels; ch++){
        cpe->ch[ch].ics.window_sequence = EIGHT_SHORT_SEQUENCE;
        cpe->ch[ch].ics.use_kb_window[0] = 1;
        cpe->ch[ch].ics.num_windows = 8;
        cpe->ch[ch].ics.swb_sizes = apc->bands128;
        cpe->ch[ch].ics.num_swb = apc->num_bands128;
        for(i = 0; i < cpe->ch[ch].ics.num_windows; i++)
            cpe->ch[ch].ics.group_len[i] = i & 1;
    }
    cpe->common_window = cpe->ch[0].ics.use_kb_window[0] == cpe->ch[1].ics.use_kb_window[0];
}

static void psy_null8_process(AACPsyContext *apc, int16_t *audio, int channel, ChannelElement *cpe)
{
    int start;
    int w, ch, g, i;

    //detect M/S
    if(apc->avctx->channels > 1 && cpe->common_window){
        start = 0;
        for(w = 0; w < cpe->ch[0].ics.num_windows; w++){
            for(g = 0; g < cpe->ch[0].ics.num_swb; g++){
                float diff = 0.0f;

                for(i = 0; i < cpe->ch[0].ics.swb_sizes[g]; i++)
                    diff += fabs(cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i]);
                cpe->ms.mask[w][g] = diff == 0.0;
            }
        }
    }
    for(ch = 0; ch < apc->avctx->channels; ch++){
        cpe->ch[ch].gain = SCALE_ONE_POS;
        for(w = 0; w < cpe->ch[ch].ics.num_windows; w++){
            for(g = 0; g < cpe->ch[ch].ics.num_swb; g++){
                cpe->ch[ch].sf_idx[w][g] = SCALE_ONE_POS;
                cpe->ch[ch].zeroes[w][g] = 0;
            }
        }
    }
    psy_create_output(apc, cpe, 0);
}

/**
 * constants for 3GPP AAC psychoacoustic model
 * @{
 */
#define PSY_3GPP_C1 3.0f                    // log2(8.0)
#define PSY_3GPP_C2 1.32192809488736234787f // log2(2.5)
#define PSY_3GPP_C3 0.55935730170421255071f // 1 - C2/C1

#define PSY_3GPP_SPREAD_LOW  1.5f // spreading factor for ascending threshold spreading  (15 dB/Bark)
#define PSY_3GPP_SPREAD_HI   3.0f // spreading factor for descending threshold spreading (30 dB/Bark)

#define PSY_3GPP_RPEMIN      0.01f
#define PSY_3GPP_RPELEV      2.0f
/**
 * @}
 */

/**
 * information for single band used by 3GPP TS26.403-inspired psychoacoustic model
 */
typedef struct Psy3gppBand{
    float energy;    ///< band energy
    float ffac;      ///< form factor
    float thr;       ///< energy threshold
    float pe;        ///< perceptual entropy
    float a;         ///< constant part in perceptual entropy
    float b;         ///< variable part in perceptual entropy
    float nl;        ///< predicted number of lines left after quantization
    float min_snr;   ///< minimal SNR
    float thr_quiet; ///< threshold in quiet
}Psy3gppBand;

/**
 * 3GPP TS26.403-inspired psychoacoustic model specific data
 */
typedef struct Psy3gppContext{
    float       barks [1024];
    float       bark_l[64];
    float       bark_s[16];
    float       s_low [1024];
    float       s_hi  [1024];
    Psy3gppBand band[2][128];
    Psy3gppBand prev_band[2][128];
    int         reservoir;
    int         avg_bits;
    float       a[2];
    float       b[2];
    float       thr[2];
    float       ath[64];
}Psy3gppContext;

/**
 * Calculate Bark value for given line.
 */
static inline float calc_bark(float f)
{
    return 13.3f * atanf(0.00076f * f) + 3.5f * atanf((f / 7500.0f) * (f / 7500.0f));
}

#define ATH_ADD 4
/**
 * Calculate ATH value for given frequency.
 * Borrowed from Lame.
 */
static inline float ath(float f, float add)
{
    f /= 1000.0f;
    return   3.64 * pow(f, -0.8)
            - 6.8  * exp(-0.6  * (f - 3.4) * (f - 3.4))
            + 6.0  * exp(-0.15 * (f - 8.7) * (f - 8.7))
            + (0.6 + 0.04 * add) * 0.001 * f * f * f * f;
}

static int psy_3gpp_init(AACPsyContext *apc)
{
    Psy3gppContext *pctx;
    int i, g, start;
    float prev, minscale, minath;
    apc->model_priv_data = av_mallocz(sizeof(Psy3gppContext));
    pctx = (Psy3gppContext*) apc->model_priv_data;

    for(i = 0; i < 1024; i++)
        pctx->barks[i] = calc_bark(i * apc->avctx->sample_rate / 2048.0);
    i = 0;
    prev = 0.0;
    for(g = 0; g < apc->num_bands1024; g++){
        i += apc->bands1024[g];
        pctx->bark_l[g] = (pctx->barks[i - 1] + prev) / 2.0;
        prev = pctx->barks[i - 1];
    }
    for(g = 0; g < apc->num_bands1024 - 1; g++){
        pctx->s_low[g] = pow(10.0, -(pctx->bark_l[g+1] - pctx->bark_l[g]) * PSY_3GPP_SPREAD_LOW);
        pctx->s_hi [g] = pow(10.0, -(pctx->bark_l[g+1] - pctx->bark_l[g]) * PSY_3GPP_SPREAD_HI);
    }
    start = 0;
    minath = ath(3410, ATH_ADD);
    for(g = 0; g < apc->num_bands1024; g++){
        minscale = ath(apc->avctx->sample_rate * start / 1024.0, ATH_ADD);
        for(i = 1; i < apc->bands1024[g]; i++){
            minscale = fminf(minscale, ath(apc->avctx->sample_rate * (start + i) / 1024.0 / 2.0, ATH_ADD));
        }
        pctx->ath[g] = minscale - minath;
        start += apc->bands1024[g];
    }

    pctx->avg_bits = apc->avctx->bit_rate * 1024 / apc->avctx->sample_rate;
    return 0;
}

/**
 * Tell encoder which window types to use.
 * @see 3GPP TS26.403 5.4.1
 */
static void psy_3gpp_window(AACPsyContext *apc, int16_t *audio, int channel, ChannelElement *cpe)
{
    int ch;

//XXX: stub, because encoder does not support long to short window transition yet :(
    for(ch = 0; ch < apc->avctx->channels; ch++){
        cpe->ch[ch].ics.window_sequence = ONLY_LONG_SEQUENCE;
        cpe->ch[ch].ics.use_kb_window[0] = 1;
        cpe->ch[ch].ics.num_windows = 1;
        cpe->ch[ch].ics.swb_sizes = apc->bands1024;
        cpe->ch[ch].ics.num_swb = apc->num_bands1024;
        cpe->ch[ch].ics.group_len[0] = 0;
    }
    cpe->common_window = cpe->ch[0].ics.use_kb_window[0] == cpe->ch[1].ics.use_kb_window[0];
}

/**
 * Modify threshold by adding some value in loudness domain.
 * @see 3GPP TS26.403 5.6.1.1.1
 */
static inline float modify_thr(float thr, float r){
    float t;
    t = pow(thr, 0.25) + r;
    return t*t*t*t;
}

/**
 * Determine scalefactors and prepare coefficients for encoding.
 * @see 3GPP TS26.403 5.4
 */
static void psy_3gpp_process(AACPsyContext *apc, int16_t *audio, int channel, ChannelElement *cpe)
{
    int start;
    int ch, g, i;
    int prev_scale;
    Psy3gppContext *pctx = (Psy3gppContext*) apc->model_priv_data;
    float stereo_att, pe_target;
    int bits_avail;

    //calculate and apply stereo attenuation factor - 5.2
    if(apc->avctx->channels > 1){
        float l, r;
        stereo_att = 1.0 / 2.0; //XXX: find some way to determine it
        for(i = 0; i < 1024; i++){
            l = cpe->ch[0].coeffs[i];
            r = cpe->ch[1].coeffs[i];
            cpe->ch[0].coeffs[i] = (0.5 + stereo_att) * l + (0.5 - stereo_att) * r;
            cpe->ch[1].coeffs[i] = (0.5 - stereo_att) * l + (0.5 + stereo_att) * r;
        }
    }

    //calculate energies, initial thresholds and related values - 5.4.2
    memset(pctx->band, 0, sizeof(pctx->band));
    for(ch = 0; ch < apc->avctx->channels; ch++){
        start = 0;
        cpe->ch[ch].gain = 0;
        for(g = 0; g < apc->num_bands1024; g++){
            for(i = 0; i < apc->bands1024[g]; i++)
                pctx->band[ch][g].energy +=  cpe->ch[ch].coeffs[start+i] *  cpe->ch[ch].coeffs[start+i];
            pctx->band[ch][g].energy *= 1048576.0;
            pctx->band[ch][g].thr = pctx->band[ch][g].energy * 0.001258925f;
            start += apc->bands1024[g];
            if(pctx->band[ch][g].energy != 0.0){
                float ffac = 0.0;

                for(i = 0; i < apc->bands1024[g]; i++)
                    ffac += sqrt(FFABS(cpe->ch[ch].coeffs[start+i]));
                pctx->band[ch][g].ffac = ffac * 32.0;

                pctx->band[ch][g].nl = pctx->band[ch][g].ffac / pow(pctx->band[ch][g].energy/apc->bands1024[g], 0.25);
                if(pctx->band[ch][g].energy / pctx->band[ch][g].thr >= 8.0){
                    pctx->band[ch][g].a = pctx->band[ch][g].nl * log2(pctx->band[ch][g].energy);
                    pctx->band[ch][g].b = pctx->band[ch][g].nl;
                }else{
                    pctx->band[ch][g].a = pctx->band[ch][g].nl * (PSY_3GPP_C2 + PSY_3GPP_C3 * log2(pctx->band[ch][g].energy));
                    pctx->band[ch][g].b = pctx->band[ch][g].nl * PSY_3GPP_C3;
                }
                pctx->band[ch][g].pe = pctx->band[ch][g].a - pctx->band[ch][g].b * log2(pctx->band[ch][g].thr);
                pctx->band[ch][g].min_snr = 1.0 / (pow(2.0, pctx->band[ch][g].pe / apc->bands1024[g]) - 1.5);
                if(pctx->band[ch][g].min_snr < 1.26f)     pctx->band[ch][g].min_snr = 1.26f;
                if(pctx->band[ch][g].min_snr > 316.2277f) pctx->band[ch][g].min_snr = 316.2277f;

                cpe->ch[ch].zeroes[0][g] = 0;
            }else{
                cpe->ch[ch].zeroes[0][g] = 1;
            }
            pctx->a[ch]   += pctx->band[ch][g].a;
            pctx->b[ch]   += pctx->band[ch][g].b;
            pctx->thr[ch] += pctx->band[ch][g].thr;
        }
        pctx->a[ch]   /= 1024.0f;
        pctx->b[ch]   /= 1024.0f;
        pctx->thr[ch] /= 1024.0f;
    }

    //modify thresholds - spread, threshold in quiet - 5.4.3
    for(ch = 0; ch < apc->avctx->channels; ch++){
         for(g = 1; g < apc->num_bands1024; g++)
             pctx->band[ch][g].thr = FFMAX(pctx->band[ch][g].thr, pctx->band[ch][g-1].thr * pctx->s_low[g-1]);
         for(g = apc->num_bands1024 - 2; g >= 0; g--)
             pctx->band[ch][g].thr = FFMAX(pctx->band[ch][g].thr, pctx->band[ch][g+1].thr * pctx->s_hi[g+1]);
         for(g = 0; g < apc->num_bands1024; g++){
             pctx->band[ch][g].thr_quiet = FFMAX(pctx->band[ch][g].thr, pctx->ath[g] *0.5);
             pctx->band[ch][g].thr = fmaxf(PSY_3GPP_RPEMIN*pctx->band[ch][g].thr_quiet, fminf(pctx->band[ch][g].thr_quiet, PSY_3GPP_RPELEV*pctx->prev_band[ch][g].thr_quiet));
         }
    }

    // M/S detection - 5.5.2
    if(apc->avctx->channels > 1 && cpe->common_window){
        start = 0;
        for(g = 0; g < cpe->ch[0].ics.num_swb; g++){
            double en_m = 0.0, en_s = 0.0, l1;
            float m, s;

            cpe->ms.mask[0][g] = 0;
            if(pctx->band[0][g].energy + pctx->band[1][g].energy == 0.0)
                continue;
            for(i = 0; i < cpe->ch[0].ics.swb_sizes[g]; i++){
                m = (cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i]) / 2.0;
                s = (cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i]) / 2.0;
                en_m += m*m;
                en_s += s*s;
            }
            l1 = FFMIN(pctx->band[0][g].thr, pctx->band[1][g].thr);
            l1 = l1*l1 / (en_m + en_s);
            if(l1 >= (pctx->band[0][g].thr * pctx->band[1][g].thr / (pctx->band[0][g].energy + pctx->band[1][g].energy)))
                cpe->ms.mask[0][g] = 1;
        }
    }

    //bitrate reduction - 5.6.1
    //TODO: add more that first step estimation
    pctx->reservoir += pctx->avg_bits - apc->avctx->frame_bits;
    bits_avail = pctx->avg_bits + pctx->reservoir;
    bits_avail = FFMIN(bits_avail, pctx->avg_bits * 1.5);
    pe_target = 1.18f * bits_avail / apc->avctx->channels / 1024.0f;
    for(ch = 0; ch < apc->avctx->channels; ch++){
        float t0, pe, r;
        if(pctx->b[ch] == 0.0f) continue;
        for(i = 0; i < 2; i++){
            pe = pctx->a[ch] - pctx->b[ch] * 4.0f * log2(pow(pctx->thr[ch]/cpe->ch[ch].ics.num_swb, 0.25));
            t0 = pow(2.0, (pctx->a[ch] - pe)        / (4.0 * pctx->b[ch]));
            r  = pow(2.0, (pctx->a[ch] - pe_target) / (4.0 * pctx->b[ch])) - t0;

            //add correction factor to thresholds
            for(g = 0; g < apc->num_bands1024; g++)
                pctx->band[ch][g].thr = modify_thr(pctx->band[ch][g].thr, r);
        }
    }

    //determine scalefactors - 5.6.2
    //TODO: quantization optimization, scalefactor difference reduction
    for(ch = 0; ch < apc->avctx->channels; ch++){
        int min_scale = 256;
        prev_scale = -1;
        cpe->ch[ch].gain = 0;
        for(g = 0; g < apc->num_bands1024; g++){
            if(cpe->ch[ch].zeroes[0][g]) continue;
            //spec gives constant for lg() but we scaled it for log2()
            cpe->ch[ch].sf_idx[0][g] = (int)(2.66667 * (log2(6.75*pctx->band[ch][g].thr) - log2(pctx->band[ch][g].ffac)));
            cpe->ch[ch].sf_idx[0][g] = av_clip(cpe->ch[ch].sf_idx[0][g], 0, 255);
            if(prev_scale != -1)
                cpe->ch[ch].sf_idx[0][g] = av_clip(cpe->ch[ch].sf_idx[0][g], prev_scale - SCALE_MAX_DIFF, prev_scale + SCALE_MAX_DIFF);
            prev_scale = cpe->ch[ch].sf_idx[0][g];
        }
        //limit scalefactors
        for(g = 0; g < apc->num_bands1024; g++){
            if(cpe->ch[ch].zeroes[0][g]) continue;
            min_scale = FFMIN(min_scale, cpe->ch[ch].sf_idx[0][g]);
        }
        for(g = 0; g < apc->num_bands1024; g++){
            if(cpe->ch[ch].zeroes[0][g]) continue;
            cpe->ch[ch].sf_idx[0][g] = FFMIN(cpe->ch[ch].sf_idx[0][g], min_scale + SCALE_MAX_DIFF);
            if(!cpe->ch[ch].gain) cpe->ch[ch].gain = cpe->ch[ch].sf_idx[0][g];
        }
    }

    memcpy(pctx->prev_band, pctx->band, sizeof(pctx->band));
    psy_create_output(apc, cpe, 0);
}

static void psy_3gpp_end(AACPsyContext *apc)
{
    av_freep(&apc->model_priv_data);
}

static const AACPsyModel psy_models[AAC_NB_PSY_MODELS] =
{
    {
       "Null model",
        NULL,
        psy_null_window,
        psy_null_process,
        NULL,
    },
    {
       "Null model - short windows",
        NULL,
        psy_null8_window,
        psy_null8_process,
        NULL,
    },
    {
       "3GPP TS 26.403-inspired model",
        psy_3gpp_init,
        psy_3gpp_window,
        psy_3gpp_process,
        psy_3gpp_end,
    },
};

int ff_aac_psy_init(AACPsyContext *ctx, AVCodecContext *avctx, int model, int flags,
                    const uint8_t *bands1024, int num_bands1024,
                    const uint8_t *bands128,  int num_bands128)
{
    int i;

    if(model >= AAC_NB_PSY_MODELS || !psy_models[model].window || !psy_models[model].process){
         av_log(avctx, AV_LOG_ERROR, "Invalid psy model\n");
         return -1;
    }

    for (i = 0; i < 340; i++)
        pow2sf_tab[i] = pow(2, (i - 200)/4.);

    ctx->avctx = avctx;
    ctx->bands1024 = bands1024;
    ctx->num_bands1024 = num_bands1024;
    ctx->bands128 = bands128;
    ctx->num_bands128 = num_bands128;
    dsputil_init(&ctx->dsp, avctx);
    ctx->model = &psy_models[model];

    if(ctx->model->init)
        return ctx->model->init(ctx);
    return 0;
}

void ff_aac_psy_suggest_window(AACPsyContext *ctx, int16_t *audio, int channel, ChannelElement *cpe)
{
    ctx->model->window(ctx, audio, channel, cpe);
}

void ff_aac_psy_analyze(AACPsyContext *ctx, int16_t *audio, int channel, ChannelElement *cpe)
{
    ctx->model->process(ctx, audio, channel, cpe);
}

void ff_aac_psy_end(AACPsyContext *ctx)
{
    if(ctx->model->end)
        return ctx->model->end(ctx);
}
