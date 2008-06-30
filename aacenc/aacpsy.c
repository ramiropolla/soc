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

static void psy_null_window(AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe)
{
    int ch;

    for(ch = 0; ch < apc->avctx->channels; ch++){
        cpe->ch[ch].ics.window_sequence = ONLY_LONG_SEQUENCE;
        cpe->ch[ch].ics.window_shape = 1;
        cpe->ch[ch].ics.num_windows = 1;
        cpe->ch[ch].ics.swb_sizes = apc->bands1024;
        cpe->ch[ch].ics.num_swb = apc->num_bands1024;
        cpe->ch[ch].ics.group_len[0] = 0;
    }
    cpe->common_window = cpe->ch[0].ics.window_shape == cpe->ch[1].ics.window_shape;
}

static void psy_null_process(AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe)
{
    int start, sum, maxsfb;
    int ch, g, i;
    int pulses, poff[4], pamp[4];

    //detect M/S
    if(apc->avctx->channels > 1 && cpe->common_window){
        start = 0;
        for(g = 0; g < apc->num_bands1024; g++){
            float diff = 0.0f;

            for(i = 0; i < apc->bands1024[g]; i++)
                diff += fabs(cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i]);
            cpe->ms.mask[0][g] = diff == 0.0;
        }
    }
    for(ch = 0; ch < apc->avctx->channels; ch++){
        start = 0;
        cpe->ch[ch].gain = SCALE_ONE_POS;
        cpe->ch[ch].pulse.present = 0;
        for(g = 0; g < apc->num_bands1024; g++){
            cpe->ch[ch].sf_idx[0][g] = SCALE_ONE_POS;
            //apply M/S
            if(!ch && cpe->ms.mask[0][g]){
                for(i = 0; i < apc->bands1024[g]; i++){
                    cpe->ch[0].coeffs[start+i] = (cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i]) / 2.0;
                    cpe->ch[1].coeffs[start+i] =  cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                }
            }
            sum = convert_coeffs(cpe->ch[ch].coeffs + start, cpe->ch[ch].icoefs + start, apc->bands1024[g], cpe->ch[ch].sf_idx[0][g]);
            cpe->ch[ch].zeroes[0][g] = !sum;
            //try finding pulses
            if(!cpe->ch[ch].pulse.present){
                pulses = 0;
                memset(poff,0,sizeof(poff));
                memset(pamp,0,sizeof(pamp));
                for(i = 0; i < apc->bands1024[g]; i++){
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
                    cpe->ch[ch].pulse.num_pulse_minus1 = pulses - 1;
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
            start += apc->bands1024[g];
        }
        for(maxsfb = apc->num_bands1024; maxsfb > 0 && cpe->ch[ch].zeroes[0][maxsfb-1]; maxsfb--);
        cpe->ch[ch].ics.max_sfb = maxsfb;
    }
    if(apc->avctx->channels > 1 && cpe->common_window){
        int msc = 0;
        cpe->ch[0].ics.max_sfb = FFMAX(cpe->ch[0].ics.max_sfb, cpe->ch[1].ics.max_sfb);
        cpe->ch[1].ics.max_sfb = cpe->ch[0].ics.max_sfb;
        for(i = 0; i < cpe->ch[0].ics.max_sfb; i++)
            if(cpe->ms.mask[0][i]) msc++;
        if(msc == 0 || cpe->ch[0].ics.max_sfb == 0) cpe->ms.present = 0;
        else cpe->ms.present = msc < cpe->ch[0].ics.max_sfb ? 1 : 2;
    }
}

static void psy_null8_window(AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe)
{
    int ch, i;

    for(ch = 0; ch < apc->avctx->channels; ch++){
        cpe->ch[ch].ics.window_sequence = EIGHT_SHORT_SEQUENCE;
        cpe->ch[ch].ics.window_shape = 1;
        cpe->ch[ch].ics.num_windows = 8;
        cpe->ch[ch].ics.swb_sizes = apc->bands128;
        cpe->ch[ch].ics.num_swb = apc->num_bands128;
        for(i = 0; i < cpe->ch[ch].ics.num_windows; i++)
            cpe->ch[ch].ics.group_len[i] = i & 1;
    }
    cpe->common_window = cpe->ch[0].ics.window_shape == cpe->ch[1].ics.window_shape;
}

static void psy_null8_process(AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe)
{
    int start, sum, cmaxsfb, maxsfb;
    int w, w2, ch, g, i;

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
        start = 0;
        cpe->ch[ch].gain = SCALE_ONE_POS;
        maxsfb = 0;
        for(w = 0; w < cpe->ch[ch].ics.num_windows; w++){
            for(g = 0; g < cpe->ch[ch].ics.num_swb; g++){
                sum = 0;
                cpe->ch[ch].sf_idx[w][g] = SCALE_ONE_POS;
                //apply M/S
                if(!ch && cpe->ms.mask[w][g]){
                    for(i = 0; i < cpe->ch[ch].ics.swb_sizes[g]; i++){
                        cpe->ch[0].coeffs[start+i] = (cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i]) / 2.0;
                        cpe->ch[1].coeffs[start+i] =  cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                    }
                }
                sum = convert_coeffs(cpe->ch[ch].coeffs + start, cpe->ch[ch].icoefs + start, cpe->ch[ch].ics.swb_sizes[g], cpe->ch[ch].sf_idx[w][g]);
                cpe->ch[ch].zeroes[w][g] = !sum;
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

void ff_aac_psy_suggest_window(AACPsyContext *ctx, int16_t *audio, int channel, cpe_struct *cpe)
{
    ctx->model->window(ctx, audio, channel, cpe);
}

void ff_aac_psy_analyze(AACPsyContext *ctx, int16_t *audio, int channel, cpe_struct *cpe)
{
    ctx->model->process(ctx, audio, channel, cpe);
}

void ff_aac_psy_end(AACPsyContext *ctx)
{
    if(ctx->model->end)
        return ctx->model->end(ctx);
}
