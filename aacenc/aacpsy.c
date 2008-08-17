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
#include "aacpsy.h"
#include "aactab.h"

/***********************************
 *              TODOs:
 * General:
 * better audio preprocessing (add DC highpass filter?)
 * more psy models
 * maybe improve coefficient quantization function in some way
 *
 * 3GPP-based psy model:
 * thresholds linearization after their modifications for attaining given bitrate
 * try other bitrate controlling mechanism (maybe use ratecontrol.c?)
 * control quality for quality-based output
 **********************************/

/**
 * Quantize one coefficient.
 * @return absolute value of the quantized coefficient
 * @see 3GPP TS26.403 5.6.2 "Scalefactor determination"
 */
static av_always_inline int quant(float coef, const float Q)
{
    return av_clip((int)(pow(fabsf(coef) * Q, 0.75) + 0.4054), 0, 8191);
}

/**
 * Convert coefficients to integers.
 * @return sum of coefficient absolute values
 */
static inline int quantize_coeffs(float *in, int *out, int size, int scale_idx)
{
    int i, sign, sum = 0;
    const float Q = ff_aac_pow2sf_tab[200 - scale_idx + SCALE_ONE_POS - SCALE_DIV_512];
    for(i = 0; i < size; i++){
        sign = in[i] > 0.0;
        out[i] = quant(in[i], Q);
        sum += out[i];
        if(sign) out[i] = -out[i];
    }
    return sum;
}

static inline float calc_distortion(float *c, int size, int scale_idx)
{
    int i;
    int q;
    float coef, unquant, sum = 0.0f;
    const float Q  = ff_aac_pow2sf_tab[200 - scale_idx + SCALE_ONE_POS - SCALE_DIV_512];
    const float IQ = ff_aac_pow2sf_tab[200 + scale_idx - SCALE_ONE_POS + SCALE_DIV_512];
    for(i = 0; i < size; i++){
        coef = fabs(c[i]);
        q = quant(c[i], Q);
        unquant = (q * cbrt(q)) * IQ;
        sum += (coef - unquant) * (coef - unquant);
    }
    return sum;
}

/**
 * Produce integer coefficients from scalefactors provided by the model.
 */
static void psy_create_output(AACPsyContext *apc, ChannelElement *cpe, int chans)
{
    int i, w, w2, wg, g, ch;
    int start, sum, maxsfb, cmaxsfb;

    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        start = 0;
        maxsfb = 0;
        cpe->ch[ch].pulse.num_pulse = 0;
        for(w = 0; w < ics->num_windows*16; w += 16){
            for(g = 0; g < ics->num_swb; g++){
                sum = 0;
                //apply M/S
                if(!ch && cpe->ms_mask[w + g]){
                    for(i = 0; i < ics->swb_sizes[g]; i++){
                        cpe->ch[0].coeffs[start+i] = (cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i]) / 2.0;
                        cpe->ch[1].coeffs[start+i] =  cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                    }
                }
                if(!cpe->ch[ch].zeroes[w + g])
                    sum = quantize_coeffs(cpe->ch[ch].coeffs + start,
                                          cpe->ch[ch].icoefs + start,
                                          ics->swb_sizes[g],
                                          cpe->ch[ch].sf_idx[w + g]);
                else
                    memset(cpe->ch[ch].icoefs + start, 0, ics->swb_sizes[g] * sizeof(cpe->ch[0].icoefs[0]));
                cpe->ch[ch].zeroes[w + g] = !sum;
                start += ics->swb_sizes[g];
            }
            for(cmaxsfb = ics->num_swb; cmaxsfb > 0 && cpe->ch[ch].zeroes[w+cmaxsfb-1]; cmaxsfb--);
            maxsfb = FFMAX(maxsfb, cmaxsfb);
        }
        ics->max_sfb = maxsfb;

        //adjust zero bands for window groups
        w = 0;
        for(wg = 0; wg < ics->num_window_groups; wg++){
            for(g = 0; g < ics->max_sfb; g++){
                i = 1;
                for(w2 = 0; w2 < ics->group_len[wg]*16; w2 += 16){
                    if(!cpe->ch[ch].zeroes[w + w2 + g]){
                        i = 0;
                        break;
                    }
                }
                cpe->ch[ch].zeroes[w + g] = i;
            }
            w += ics->group_len[wg] * 16;
        }
    }

    if(chans > 1 && cpe->common_window){
        IndividualChannelStream *ics0 = &cpe->ch[0].ics;
        IndividualChannelStream *ics1 = &cpe->ch[1].ics;
        int msc = 0;
        ics0->max_sfb = FFMAX(ics0->max_sfb, ics1->max_sfb);
        ics1->max_sfb = ics0->max_sfb;
        for(w = 0; w < ics0->num_windows*16; w += 16)
            for(i = 0; i < ics0->max_sfb; i++)
                if(cpe->ms_mask[w+i]) msc++;
        if(msc == 0 || ics0->max_sfb == 0) cpe->ms_mode = 0;
        else cpe->ms_mode = msc < ics0->max_sfb ? 1 : 2;
    }
}

static void psy_test_window(AACPsyContext *apc, int16_t *audio, int16_t *la,
                            int tag, int type, ChannelElement *cpe)
{
    int ch, i;
    int chans = type == TYPE_CPE ? 2 : 1;

    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        int prev_seq = ics->window_sequence[1];
        ics->use_kb_window[1] = ics->use_kb_window[0];
        ics->window_sequence[1] = ics->window_sequence[0];
        switch(ics->window_sequence[0]){
        case ONLY_LONG_SEQUENCE:
            if(prev_seq == ONLY_LONG_SEQUENCE)
                ics->window_sequence[0] = LONG_START_SEQUENCE;
            break;
        case LONG_START_SEQUENCE:
            ics->window_sequence[0] = EIGHT_SHORT_SEQUENCE;
            break;
        case EIGHT_SHORT_SEQUENCE:
            if(prev_seq == EIGHT_SHORT_SEQUENCE)
                ics->window_sequence[0] = LONG_STOP_SEQUENCE;
            break;
        case LONG_STOP_SEQUENCE:
            ics->window_sequence[0] = ONLY_LONG_SEQUENCE;
            break;
        }

        if(ics->window_sequence[0] != EIGHT_SHORT_SEQUENCE){
            ics->use_kb_window[0] = 1;
            ics->num_windows = 1;
            ics->swb_sizes = apc->bands1024;
            ics->num_swb = apc->num_bands1024;
            ics->num_window_groups = 1;
            ics->group_len[0] = 1;
        }else{
            ics->use_kb_window[0] = 1;
            ics->num_windows = 8;
            ics->swb_sizes = apc->bands128;
            ics->num_swb = apc->num_bands128;
            ics->num_window_groups = 4;
            for(i = 0; i < 4; i++)
                ics->group_len[i] = 2;
        }
    }
    cpe->common_window = cpe->ch[0].ics.use_kb_window[0] == cpe->ch[1].ics.use_kb_window[0];
}

static void psy_test_process(AACPsyContext *apc, int tag, int type, ChannelElement *cpe)
{
    int start;
    int w, ch, g, i;
    int chans = type == TYPE_CPE ? 2 : 1;

    //detect M/S
    if(chans > 1 && cpe->common_window){
        start = 0;
        for(w = 0; w < cpe->ch[0].ics.num_windows*16; w += 16){
            for(g = 0; g < cpe->ch[0].ics.num_swb; g++){
                float diff = 0.0f;

                for(i = 0; i < cpe->ch[0].ics.swb_sizes[g]; i++)
                    diff += fabs(cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i]);
                cpe->ms_mask[w + g] = diff == 0.0;
            }
        }
    }
    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        for(w = 0; w < ics->num_windows*16; w += 16){
            for(g = 0; g < ics->num_swb; g++){
                cpe->ch[ch].sf_idx[w + g] = SCALE_ONE_POS;
                cpe->ch[ch].zeroes[w + g] = 0;
            }
        }
    }
    psy_create_output(apc, cpe, chans);
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
 * single/pair channel context for psychoacoustic model
 */
typedef struct Psy3gppChannel{
    float       a[2];                       ///< parameter used for perceptual entropy - constant part
    float       b[2];                       ///< parameter used for perceptual entropy - variable part
    float       pe[2];                      ///< channel perceptual entropy
    float       thr[2];                     ///< channel thresholds sum
    Psy3gppBand band[2][128];               ///< bands information
    Psy3gppBand prev_band[2][128];          ///< bands information from the previous frame

    float       win_energy[2];              ///< sliding average of channel energy
    float       iir_state[2][2];            ///< hi-pass IIR filter state
    uint8_t     next_grouping[2];           ///< stored grouping scheme for the next frame (in case of 8 short window sequence)
    enum WindowSequence next_window_seq[2]; ///< window sequence to be used in the next frame
}Psy3gppChannel;

/**
 * psychoacoustic model frame type-dependent coefficients
 */
typedef struct Psy3gppCoeffs{
    float ath       [64]; ///< absolute threshold of hearing per bands
    float barks     [64]; ///< Bark value for each spectral band in long frame
    float spread_low[64]; ///< spreading factor for low-to-high threshold spreading in long frame
    float spread_hi [64]; ///< spreading factor for high-to-low threshold spreading in long frame
}Psy3gppCoeffs;

/**
 * 3GPP TS26.403-inspired psychoacoustic model specific data
 */
typedef struct Psy3gppContext{
    Psy3gppCoeffs psy_coef[2];
    int         reservoir;    ///< bit reservoir fullness
    int         avg_bits;     ///< average frame size of bits for CBR
    Psy3gppChannel *ch;
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

static av_cold int psy_3gpp_init(AACPsyContext *apc, int elements)
{
    Psy3gppContext *pctx;
    float barks[1024];
    int i, j, g, start;
    float prev, minscale, minath;
    apc->model_priv_data = av_mallocz(sizeof(Psy3gppContext));
    pctx = (Psy3gppContext*) apc->model_priv_data;

    for(i = 0; i < 1024; i++)
        barks[i] = calc_bark(i * apc->avctx->sample_rate / 2048.0);
    minath = ath(3410, ATH_ADD);
    for(j = 0; j < 2; j++){
        Psy3gppCoeffs *coeffs = &pctx->psy_coef[j];
        int bands = j ? apc->num_bands128 : apc->num_bands1024;
        i = 0;
        prev = 0.0;
        for(g = 0; g < bands; g++){
            i += j ? apc->bands128[g] : apc->bands1024[g];
            coeffs->barks[g] = (barks[i - 1] + prev) / 2.0;
            prev = barks[i - 1];
        }
        for(g = 0; g < bands - 1; g++){
            coeffs->spread_low[g] = pow(10.0, -(coeffs->barks[g+1] - coeffs->barks[g]) * PSY_3GPP_SPREAD_LOW);
            coeffs->spread_hi [g] = pow(10.0, -(coeffs->barks[g+1] - coeffs->barks[g]) * PSY_3GPP_SPREAD_HI);
        }
        start = 0;
        for(g = 0; g < bands; g++){
            int size = j ? apc->bands128[g] : apc->bands1024[g];
            minscale = ath(apc->avctx->sample_rate * start / 1024.0, ATH_ADD);
            for(i = 1; i < size; i++){
                minscale = fminf(minscale, ath(apc->avctx->sample_rate * (start + i) / 1024.0 / 2.0, ATH_ADD));
            }
            coeffs->ath[g] = minscale - minath;
            start += size;
        }
    }

    pctx->avg_bits = apc->avctx->bit_rate * 1024 / apc->avctx->sample_rate;
    pctx->ch = av_mallocz(sizeof(Psy3gppChannel) * elements);
    return 0;
}

/**
 * IIR filter used in block switching decision
 */
static float iir_filter(int in, float state[2])
{
    float ret;

    ret = 0.7548f * (in - state[0]) + 0.5095f * state[1];
    state[0] = in;
    state[1] = ret;
    return ret;
}

/**
 * window grouping information stored as bits (0 - new group, 1 - group continues)
 */
static const uint8_t window_grouping[9] = {
    0xB6, 0x6C, 0xD8, 0xB2, 0x66, 0xC6, 0x96, 0x36, 0x36
};

/**
 * Tell encoder which window types to use.
 * @see 3GPP TS26.403 5.4.1 "Blockswitching"
 */
static void psy_3gpp_window(AACPsyContext *apc, int16_t *audio, int16_t *la,
                            int tag, int type, ChannelElement *cpe)
{
    int ch;
    int chans = type == TYPE_CPE ? 2 : 1;
    int i, j;
    int br = apc->avctx->bit_rate / apc->avctx->channels;
    int attack_ratio = (br <= 16000 + 8000*chans) ? 18 : 10;
    Psy3gppContext *pctx = (Psy3gppContext*) apc->model_priv_data;
    Psy3gppChannel *pch = &pctx->ch[tag];
    uint8_t grouping[2];
    enum WindowSequence win[2];
    IndividualChannelStream *ics0 = &cpe->ch[0].ics, *ics1 = &cpe->ch[1].ics;

    if(la && !(apc->flags & PSY_MODEL_NO_SWITCH)){
        float s[8], v;
        for(ch = 0; ch < chans; ch++){
            enum WindowSequence last_window_sequence = cpe->ch[ch].ics.window_sequence[0];
            int switch_to_eight = 0;
            float sum = 0.0, sum2 = 0.0;
            int attack_n = 0;
            for(i = 0; i < 8; i++){
                for(j = 0; j < 128; j++){
                    v = iir_filter(audio[(i*128+j)*apc->avctx->channels+ch], pch->iir_state[ch]);
                    sum += v*v;
                }
                s[i] = sum;
                sum2 += sum;
            }
            for(i = 0; i < 8; i++){
                if(s[i] > pch->win_energy[ch] * attack_ratio){
                    attack_n = i + 1;
                    switch_to_eight = 1;
                    break;
                }
            }
            pch->win_energy[ch] = pch->win_energy[ch]*7/8 + sum2/64;

            switch(last_window_sequence){
            case ONLY_LONG_SEQUENCE:
                win[ch] = switch_to_eight ? LONG_START_SEQUENCE : ONLY_LONG_SEQUENCE;
                grouping[ch] = 0;
                break;
            case LONG_START_SEQUENCE:
                win[ch] = EIGHT_SHORT_SEQUENCE;
                grouping[ch] = pch->next_grouping[ch];
                break;
            case LONG_STOP_SEQUENCE:
                win[ch] = ONLY_LONG_SEQUENCE;
                grouping[ch] = 0;
                break;
            case EIGHT_SHORT_SEQUENCE:
                win[ch] = switch_to_eight ? EIGHT_SHORT_SEQUENCE : LONG_STOP_SEQUENCE;
                grouping[ch] = switch_to_eight ? pch->next_grouping[ch] : 0;
                break;
            }
            pch->next_grouping[ch] = window_grouping[attack_n];
        }
    }else{
        for(ch = 0; ch < chans; ch++){
            IndividualChannelStream *ics = &cpe->ch[ch].ics;
            win[ch] = (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE)
                      ? EIGHT_SHORT_SEQUENCE
                      : ONLY_LONG_SEQUENCE;
            grouping[ch] = (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE) ? window_grouping[0] : 0;
        }
    }

    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        ics->window_sequence[0] = win[ch];
        ics->use_kb_window[0] = 1;
        if(win[ch] != EIGHT_SHORT_SEQUENCE){
            ics->num_windows = 1;
            ics->swb_sizes = apc->bands1024;
            ics->num_swb = apc->num_bands1024;
            ics->num_window_groups = 1;
            ics->group_len[0] = 1;
        }else{
            ics->num_windows = 8;
            ics->swb_sizes = apc->bands128;
            ics->num_swb = apc->num_bands128;
            ics->num_window_groups = 0;
            ics->group_len[0] = 1;
            for(i = 0; i < 8; i++){
                if((grouping[ch] >> i) & 1){
                    ics->group_len[ics->num_window_groups - 1]++;
                }else{
                    ics->num_window_groups++;
                    ics->group_len[ics->num_window_groups - 1] = 1;
                }
            }
        }
    }
    cpe->common_window = 0;
    if(chans > 1
            && ics0->window_sequence[0] == ics1->window_sequence[0]
            && ics0->use_kb_window[0]   == ics1->use_kb_window[0]
            && !(ics0->window_sequence[0] == EIGHT_SHORT_SEQUENCE && grouping[0] != grouping[1]))
        cpe->common_window = 1;
    if(PSY_MODEL_MODE(apc->flags) > PSY_MODE_QUALITY){
        av_log(apc->avctx, AV_LOG_ERROR, "Unknown mode %d, defaulting to CBR\n", PSY_MODEL_MODE(apc->flags));
    }
}

/**
 * Modify threshold by adding some value in loudness domain.
 * @see 3GPP TS26.403 5.6.1.1.1 "Addition of noise with equal loudness"
 */
static inline float modify_thr(float thr, float r){
    float t;
    t = pow(thr, 0.25) + r;
    return (t*t)*(t*t);
}

/**
 * Calculate perceptual entropy and its corresponding values for one band.
 * @see 3GPP TS26.403 5.6.1.3 "Calculation of the reduction value"
 */
static void calc_pe(Psy3gppBand *band, int band_width)
{
    if(band->energy <= band->thr){
        band->a  = 0.0f;
        band->b  = 0.0f;
        band->nl = 0.0f;
        return;
    }
    band->nl = band->ffac / pow(band->energy/band_width, 0.25);
    if(band->energy >= band->thr * 8.0){
        band->a = band->nl * log2(band->energy);
        band->b = band->nl;
    }else{
        band->a = band->nl * (PSY_3GPP_C2 + PSY_3GPP_C3 * log2(band->energy));
        band->b = band->nl * PSY_3GPP_C3;
    }
    band->pe = band->a - band->b * log2(band->thr);
    band->min_snr = 1.0 / (pow(2.0, band->pe / band_width) - 1.5);
    band->min_snr = av_clipf(band->min_snr, 1.26f, 316.2277f);
}

/**
 * Determine scalefactor from band threshold and form factor.
 * @see 3GPP TS26.403 5.4 5.6.2 "Scalefactor determination"
 */
static inline int determine_scalefactor(Psy3gppBand *band)
{
    //spec gives constant for lg() but we scaled it for log2()
    return (int)(2.66667 * log2(6.75*band->thr/band->ffac));
}

/**
 * Determine scalefactors and prepare coefficients for encoding.
 * @see 3GPP TS26.403 5.4 "Psychoacoustic model"
 */
static void psy_3gpp_process(AACPsyContext *apc, int tag, int type, ChannelElement *cpe)
{
    int start;
    int ch, w, wg, g, i;
    int prev_scale;
    Psy3gppContext *pctx = (Psy3gppContext*) apc->model_priv_data;
    float pe_target;
    int bits_avail;
    int chans = type == TYPE_CPE ? 2 : 1;
    Psy3gppChannel *pch = &pctx->ch[tag];

    //calculate energies, initial thresholds and related values - 5.4.2 "Threshold Calculation"
    memset(pch->band, 0, sizeof(pch->band));
    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        start = 0;
        for(w = 0; w < ics->num_windows*16; w += 16){
            for(g = 0; g < ics->num_swb; g++){
                Psy3gppBand *band = &pch->band[ch][w+g];
                for(i = 0; i < ics->swb_sizes[g]; i++)
                    band->energy +=  cpe->ch[ch].coeffs[start+i] * cpe->ch[ch].coeffs[start+i];
                band->energy *= 1.0f / (512*512);
                band->thr = band->energy * 0.001258925f;
                start += ics->swb_sizes[g];
                if(band->energy != 0.0){
                    float ffac = 0.0;

                    for(i = 0; i < ics->swb_sizes[g]; i++)
                        ffac += sqrt(FFABS(cpe->ch[ch].coeffs[start+i]));
                    band->ffac = ffac / sqrt(512.0);
                }
            }
        }
    }

    //modify thresholds - spread, threshold in quiet - 5.4.3 "Spreaded Energy Calculation"
    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        Psy3gppCoeffs *coeffs = &pctx->psy_coef[ics->num_windows == 8];
        for(w = 0; w < ics->num_windows*16; w += 16){
            Psy3gppBand *band = &pch->band[ch][w];
            for(g = 1; g < ics->num_swb; g++){
                band[g].thr = FFMAX(band[g].thr, band[g-1].thr * coeffs->spread_low[g-1]);
            }
            for(g = ics->num_swb - 2; g >= 0; g--){
                band[g].thr = FFMAX(band[g].thr, band[g+1].thr * coeffs->spread_hi [g+1]);
            }
            for(g = 0; g < ics->num_swb; g++){
                band[g].thr_quiet = FFMAX(band[g].thr, coeffs->ath[g]);
                band[g].thr_quiet = fmaxf(PSY_3GPP_RPEMIN*band[g].thr_quiet,
                                          fminf(band[g].thr_quiet,
                                                PSY_3GPP_RPELEV*pch->prev_band[ch][w+g].thr_quiet));
                band[g].thr = FFMAX(band[g].thr, band[g].thr_quiet * 0.25);
            }
        }
    }

    // M/S detection - 5.5.2 "Mid/Side Stereo"
    if(chans > 1 && cpe->common_window){
        start = 0;
        for(w = 0; w < cpe->ch[0].ics.num_windows*16; w += 16){
            for(g = 0; g < cpe->ch[0].ics.num_swb; g++){
                Psy3gppBand *band0 = &pch->band[0][w+g];
                Psy3gppBand *band1 = &pch->band[1][w+g];
                double en_m = 0.0, en_s = 0.0, ff_m = 0.0, ff_s = 0.0, minthr;
                float m, s;

                cpe->ms_mask[w+g] = 0;
                if(band0->energy == 0.0 || band1->energy == 0.0)
                    continue;
                for(i = 0; i < cpe->ch[0].ics.swb_sizes[g]; i++){
                    m = cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i];
                    s = cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                    en_m += m*m;
                    en_s += s*s;
                }
                en_m *= 1.0f / (512*512*4);
                en_s *= 1.0f / (512*512*4);
                minthr = FFMIN(band0->thr, band1->thr);
                if(minthr * minthr * band0->energy * band1->energy >= band0->thr * band1->thr * en_m * en_s){
                    cpe->ms_mask[w+g] = 1;
                    band0->energy = en_m;
                    band1->energy = en_s;
                    band0->thr = en_m * 0.001258925f;
                    band1->thr = en_s * 0.001258925f;
                    for(i = 0; i < cpe->ch[0].ics.swb_sizes[g]; i++){
                        m = cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i];
                        s = cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                        ff_m += sqrt(fabs(m));
                        ff_s += sqrt(fabs(s));
                    }
                    band0->ffac = ff_m * (1.0f / 32.0f); // sqrt(512)*sqrt(2)
                    band1->ffac = ff_s * (1.0f / 32.0f);
                }
            }
        }
    }

    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        pch->a[ch] = pch->b[ch] = pch->pe[ch] = pch->thr[ch] = 0.0f;
        for(w = 0; w < ics->num_windows*16; w += 16){
            for(g = 0; g < ics->num_swb; g++){
                Psy3gppBand *band = &pch->band[ch][w+g];
                if(band->energy != 0.0)
                    calc_pe(band, ics->swb_sizes[g]);
                if(band->thr < band->energy){
                    pch->a[ch]   += band->a;
                    pch->b[ch]   += band->b;
                    pch->pe[ch]  += band->pe;
                    pch->thr[ch] += band->thr;
                }
            }
        }
    }

    switch(PSY_MODEL_MODE(apc->flags)){
    case PSY_MODE_CBR:
    case PSY_MODE_ABR:
        //bitrate reduction - 5.6.1 "Reduction of psychoacoustic requirements"
        if(PSY_MODEL_MODE(apc->flags) != PSY_MODE_ABR){
            pctx->reservoir += pctx->avg_bits - apc->avctx->frame_bits;
            bits_avail = pctx->avg_bits + pctx->reservoir;
            bits_avail = FFMIN(bits_avail, pctx->avg_bits * 1.5);
            pe_target = 1.18f * bits_avail / apc->avctx->channels * chans;
        }else{
            pe_target = pctx->avg_bits / apc->avctx->channels * chans;
        }
        for(i = 0; i < 2; i++){
            float t0, pe, r, a0 = 0.0f, pe0 = 0.0f, b0 = 0.0f;
            for(ch = 0; ch < chans; ch++){
                a0  += pch->a[ch];
                b0  += pch->b[ch];
                pe0 += pch->pe[ch];
            }
            if(pe0 == 0.0f) break;
            t0 = pow(2.0, (a0 - pe0)       / (4.0 * b0));
            r  = pow(2.0, (a0 - pe_target) / (4.0 * b0)) - t0;

            //add correction factor to thresholds and recalculate perceptual entropy
            for(ch = 0; ch < chans; ch++){
                IndividualChannelStream *ics = &cpe->ch[ch].ics;
                pch->a[ch] = pch->b[ch] = pch->pe[ch] = pch->thr[ch] = 0.0;
                pe = 0.0f;
                for(w = 0; w < ics->num_windows*16; w += 16){
                    for(g = 0; g < ics->num_swb; g++){
                        Psy3gppBand *band = &pch->band[ch][w+g];
                        band->thr = modify_thr(band->thr, r);
                        calc_pe(band, ics->swb_sizes[g]);
                        if(band->thr < band->energy){
                            pch->a[ch]   += band->a;
                            pch->b[ch]   += band->b;
                            pch->pe[ch]  += band->pe;
                            pch->thr[ch] += band->thr;
                        }
                    }
                }
            }
        }

        //determine scalefactors - 5.6.2 "Scalefactor determination"
        for(ch = 0; ch < chans; ch++){
            IndividualChannelStream *ics = &cpe->ch[ch].ics;
            for(w = 0; w < ics->num_windows*16; w += 16){
                for(g = 0; g < ics->num_swb; g++){
                    Psy3gppBand *band = &pch->band[ch][w+g];
                    cpe->ch[ch].zeroes[w+g] = band->thr >= band->energy;
                    if(cpe->ch[ch].zeroes[w+g]) continue;
                    cpe->ch[ch].sf_idx[w+g] = determine_scalefactor(band);
                }
            }
        }
        break;
    case PSY_MODE_QUALITY:
        for(ch = 0; ch < chans; ch++){
            IndividualChannelStream *ics = &cpe->ch[ch].ics;
            start = 0;
            for(w = 0; w < ics->num_windows*16; w += 16){
                for(g = 0; g < ics->num_swb; g++){
                    Psy3gppBand *band = &pch->band[ch][w+g];
                    if(band->thr >= band->energy){
                        cpe->ch[ch].sf_idx[w+g] = 0;
                        cpe->ch[ch].zeroes[w+g] = 1;
                    }else{
                        cpe->ch[ch].zeroes[w+g] = 0;
                        cpe->ch[ch].sf_idx[w+g] = determine_scalefactor(band);
                        while(cpe->ch[ch].sf_idx[w+g] > 3){
                            float dist = calc_distortion(cpe->ch[ch].coeffs + start,
                                                         ics->swb_sizes[g],
                                                         SCALE_ONE_POS + cpe->ch[ch].sf_idx[w+g]);
                            if(dist < band->thr) break;
                            cpe->ch[ch].sf_idx[w+g] -= 3;
                        }
                    }
                    start += ics->swb_sizes[g];
                }
            }
        }
        break;
    }

    //limit scalefactors
    for(ch = 0; ch < chans; ch++){
        int min_scale = 256;
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        for(w = 0; w < ics->num_windows*16; w += 16)
            for(g = 0; g < ics->num_swb; g++){
                if(cpe->ch[ch].zeroes[w + g]) continue;
                min_scale = FFMIN(min_scale, cpe->ch[ch].sf_idx[w + g]);
            }
        for(w = 0; w < ics->num_windows*16; w += 16)
            for(g = 0; g < ics->num_swb; g++){
                if(cpe->ch[ch].zeroes[w + g]) continue;
                cpe->ch[ch].sf_idx[w + g] = FFMIN(cpe->ch[ch].sf_idx[w + g], min_scale + SCALE_MAX_DIFF);
            }
        for(w = 0; w < ics->num_windows*16; w += 16)
            for(g = 0; g < ics->num_swb; g++){
                if(cpe->ch[ch].zeroes[w + g])
                    cpe->ch[ch].sf_idx[w + g] = 256;
                else
                    cpe->ch[ch].sf_idx[w + g] = av_clip(SCALE_ONE_POS + cpe->ch[ch].sf_idx[w + g],
                                                        0,
                                                        SCALE_MAX_POS);
            }

        //adjust scalefactors for window groups
        w = 0;
        for(wg = 0; wg < ics->num_window_groups; wg++){
            int min_scale = 256;

            for(g = 0; g < ics->num_swb; g++){
                for(i = w; i < w + ics->group_len[wg]*16; i += 16){
                    if(cpe->ch[ch].zeroes[i + g]) continue;
                    min_scale = FFMIN(min_scale, cpe->ch[ch].sf_idx[i + g]);
                }
                for(i = w; i < w + ics->group_len[wg]*16; i += 16)
                    cpe->ch[ch].sf_idx[i + g] = min_scale;
            }
            w += ics->group_len[wg] * 16;
        }
    }

    memcpy(pch->prev_band, pch->band, sizeof(pch->band));
    psy_create_output(apc, cpe, chans);
}

static av_cold void psy_3gpp_end(AACPsyContext *apc)
{
    Psy3gppContext *pctx = (Psy3gppContext*) apc->model_priv_data;
    av_freep(&pctx->ch);
    av_freep(&apc->model_priv_data);
}

static const AACPsyModel psy_models[AAC_NB_PSY_MODELS] =
{
    {
       "Test model",
        NULL,
        psy_test_window,
        psy_test_process,
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

int av_cold ff_aac_psy_init(AACPsyContext *ctx, AVCodecContext *avctx,
                            enum AACPsyModelType model, int elements, int flags,
                            const uint8_t *bands1024, int num_bands1024,
                            const uint8_t *bands128,  int num_bands128)
{
    int i;

    if(model < 0 || model >= AAC_NB_PSY_MODELS){
         av_log(avctx, AV_LOG_ERROR, "Invalid psy model\n");
         return -1;
    }

#ifndef CONFIG_HARDCODED_TABLES
   for (i = 0; i < 316; i++)
        ff_aac_pow2sf_tab[i] = pow(2, (i - 200)/4.);
#endif /* CONFIG_HARDCODED_TABLES */

    ctx->avctx = avctx;
    ctx->flags = flags;
    ctx->bands1024 = bands1024;
    ctx->num_bands1024 = num_bands1024;
    ctx->bands128 = bands128;
    ctx->num_bands128 = num_bands128;
    ctx->model = &psy_models[model];

    if(ctx->flags & PSY_MODEL_NO_ST_ATT || PSY_MODEL_MODE(ctx->flags) == PSY_MODE_QUALITY){
        ctx->flags |= PSY_MODEL_NO_ST_ATT;
        ctx->stereo_att = 0.5f;
    }else{
        ctx->stereo_att = av_clipf(avctx->bit_rate / elements / 192000.0, 0.0f, 0.5f);
    }
    if(ctx->flags & PSY_MODEL_NO_LOWPASS || PSY_MODEL_MODE(ctx->flags) == PSY_MODE_QUALITY){
        ctx->flags |= PSY_MODEL_NO_LOWPASS;
    }else{
        int cutoff;
        cutoff = avctx->bit_rate / elements / 8;
        if(ff_lowpass_filter_init_coeffs(&ctx->lp_coeffs, avctx->sample_rate/2, cutoff) < 0){
            ctx->flags |= PSY_MODEL_NO_LOWPASS;
        }else{
            ctx->lp_state = av_mallocz(sizeof(LPFilterState) * elements * 2);
        }
    }
    if(ctx->model->init)
        return ctx->model->init(ctx, elements);
    return 0;
}

void ff_aac_psy_suggest_window(AACPsyContext *ctx, int16_t *audio, int16_t *la,
                               int tag, int type, ChannelElement *cpe)
{
    ctx->model->window(ctx, audio, la, tag, type, cpe);
}

void ff_aac_psy_analyze(AACPsyContext *ctx, int tag, int type, ChannelElement *cpe)
{
    ctx->model->process(ctx, tag, type, cpe);
}

void av_cold ff_aac_psy_end(AACPsyContext *ctx)
{
    av_freep(&ctx->lp_state);
    if(ctx->model->end)
        return ctx->model->end(ctx);
}

void ff_aac_psy_preprocess(AACPsyContext *ctx, int16_t *audio, int16_t *dest, int tag, int type)
{
    int chans = type == TYPE_CPE ? 2 : 1;
    const int chstride = ctx->avctx->channels;
    int i, ch;
    float t[2];

    if(chans == 1){
        for(ch = 0; ch < chans; ch++){
            for(i = 0; i < 1024; i++){
                dest[i * chstride + ch] = audio[i * chstride + ch];
            }
        }
    }else{
        for(i = 0; i < 1024; i++){
            if(ctx->flags & PSY_MODEL_NO_ST_ATT){
                for(ch = 0; ch < 2; ch++)
                    t[ch] = audio[i * chstride + ch];
            }else{
                t[0] = audio[i * chstride + 0] * (0.5 + ctx->stereo_att) + audio[i * chstride + 1] * (0.5 - ctx->stereo_att);
                t[1] = audio[i * chstride + 0] * (0.5 - ctx->stereo_att) + audio[i * chstride + 1] * (0.5 + ctx->stereo_att);
            }
            if(!(ctx->flags & PSY_MODEL_NO_LOWPASS)){
                LPFilterState *is = (LPFilterState*)ctx->lp_state + tag*2;
                for(ch = 0; ch < 2; ch++)
                    t[ch] = ff_lowpass_filter(&ctx->lp_coeffs, is + ch, t[ch]);
            }
            for(ch = 0; ch < 2; ch++)
                dest[i * chstride + ch] = av_clip_int16(t[ch]);
        }
    }
}

