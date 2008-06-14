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
static float pow2sf_tab[316];


#define SCALE_ONE_POS   140
#define SCALE_MAX_POS   255
#define SCALE_MAX_DIFF   60

static void psy_null_window(AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe)
{
    int ch;

    for(ch = 0; ch < apc->avctx->channels; ch++){
        cpe->ch[ch].ics.window_sequence = 0;
        cpe->ch[ch].ics.window_shape = 1;
    }
    cpe->common_window = cpe->ch[0].ics.window_shape == cpe->ch[1].ics.window_shape;
}

static void psy_null_process(AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe)
{
    int start, sum, maxsfb;
    int ch, g, i;

    //detect M/S
    if(apc->avctx->channels > 1 && cpe->common_window){
        start = 0;
        for(g = 0; g < apc->num_bands; g++){
            float diff = 0.0f;

            for(i = 0; i < apc->bands[g]; i++)
                diff += fabs(cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i]);
            cpe->ms.mask[0][g] = diff == 0.0;
        }
    }
    for(ch = 0; ch < apc->avctx->channels; ch++){
        start = 0;
        cpe->ch[ch].gain = SCALE_ONE_POS;
        for(g = 0; g < apc->num_bands; g++){
            sum = 0;
            cpe->ch[ch].sf_idx[g] = SCALE_ONE_POS;
            //apply M/S
            if(!ch && cpe->ms.mask[0][g]){
                for(i = 0; i < apc->bands[g]; i++){
                    cpe->ch[0].coeffs[start+i] = (cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i]) / 2.0;
                    cpe->ch[1].coeffs[start+i] =  cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                }
            }
            for(i = 0; i < apc->bands[g]; i++){
                cpe->ch[ch].icoefs[start+i] = av_clip((int)(roundf(cpe->ch[ch].coeffs[start+i] / pow2sf_tab[cpe->ch[ch].sf_idx[g]+60])), -8191, 8191);
                sum += !!cpe->ch[ch].icoefs[start+i];
            }
            cpe->ch[ch].zeroes[g] = !sum;
            start += apc->bands[g];
        }
        for(maxsfb = apc->num_bands; maxsfb > 0 && cpe->ch[ch].zeroes[maxsfb-1]; maxsfb--);
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

static const AACPsyModel psy_models[AAC_NB_PSY_MODELS] =
{
    {
       "Null model",
        NULL,
        psy_null_window,
        psy_null_process,
        NULL,
    },
};

int ff_aac_psy_init(AACPsyContext *ctx, AVCodecContext *avctx, int model, int flags,
                    const uint8_t *bands, int num_bands)
{
    int i;

    if(model >= AAC_NB_PSY_MODELS || !psy_models[model].window || !psy_models[model].process){
         av_log(avctx, AV_LOG_ERROR, "Invalid psy model\n");
         return -1;
    }

    for (i = 0; i < 316; i++)
        pow2sf_tab[i] = pow(2, (i - 200)/4.);

    ctx->avctx = avctx;
    ctx->bands = bands;
    ctx->num_bands = num_bands;
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
