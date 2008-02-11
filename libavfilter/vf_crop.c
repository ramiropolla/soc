/*
 * video crop filter
 * copyright (c) 2007 Bobby Bingham
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

#include <stdio.h>

#include "avfilter.h"

typedef struct
{
    int  x,  y,  w,  h;
    int cx, cy, cw, ch;

    int bpp;                //< bytes per pixel
    int hsub, vsub;         //< chroma subsampling
} CropContext;

static int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    CropContext *crop = ctx->priv;

    /* default parameters */
    crop->x = 0;
    crop->y = 0;
    crop->w = -1;
    crop->h = -1;

    if(args)
        sscanf(args, "%d:%d:%d:%d", &crop->x, &crop->y, &crop->w, &crop->h);

    return 0;
}

static int config_input(AVFilterLink *link)
{
    CropContext *crop = link->dst->priv;

    crop->cx = FFMIN(crop->x, link->w - 1);
    crop->cy = FFMIN(crop->y, link->h - 1);
    crop->cw = FFMIN(crop->w, link->w - crop->cx);
    crop->ch = FFMIN(crop->h, link->h - crop->cy);

    if(crop->cw <= 0) crop->cw = link->w - crop->cx;
    if(crop->ch <= 0) crop->ch = link->h - crop->cy;

    switch(link->format) {
    case PIX_FMT_RGB32:
    case PIX_FMT_BGR32:
        crop->bpp = 4;
        break;
    case PIX_FMT_RGB24:
    case PIX_FMT_BGR24:
        crop->bpp = 3;
        break;
    case PIX_FMT_RGB565:
    case PIX_FMT_RGB555:
    case PIX_FMT_BGR565:
    case PIX_FMT_BGR555:
    case PIX_FMT_GRAY16BE:
    case PIX_FMT_GRAY16LE:
        crop->bpp = 2;
        break;
    default:
        crop->bpp = 1;
    }

    avcodec_get_chroma_sub_sample(link->format, &crop->hsub, &crop->vsub);
    crop->cx &= ~((1 << crop->hsub) - 1);
    crop->cw &= ~((1 << crop->hsub) - 1);
    crop->ch &= ~((1 << crop->vsub) - 1);

    return 0;
}

static int config_output(AVFilterLink *link)
{
    CropContext *crop = link->src->priv;

    link->w = crop->cw;
    link->h = crop->ch;

    return 0;
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    CropContext *crop = link->dst->priv;
    AVFilterPicRef *ref2 = avfilter_ref_pic(picref, ~0);
    int i;

    ref2->w = crop->cw;
    ref2->h = crop->ch;
    ref2->data[0] += crop->cy * ref2->linesize[0];
    ref2->data[0] += crop->cx * crop->bpp;
    for(i = 1; i < 4; i ++) {
        if(ref2->data[i]) {
            ref2->data[i] += (crop->cy >> crop->vsub) * ref2->linesize[i];
            ref2->data[i] +=  crop->cx >> crop->hsub;
        }
    }

    avfilter_start_frame(link->dst->outputs[0], ref2);
}

static void draw_slice(AVFilterLink *link, int y, int h)
{
    AVFilterContext *ctx = link->dst;
    CropContext *crop = ctx->priv;

    int top = y;
    int height = h;

    if(y >= crop->cy + crop->ch || y + h <= crop->cy) return;

    if(top < crop->cy) {
        height -=  crop->cy - top;
        top     =  crop->cy;
    }
    if(top + height > crop->cy + crop->ch)
        height = crop->cy + crop->ch - top;

    avfilter_draw_slice(ctx->outputs[0], top - crop->cy, height);
}

AVFilter avfilter_vf_crop =
{
    .name      = "crop",
    .priv_size = sizeof(CropContext),

    .init      = init,

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .config_props    = config_input, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .config_props    = config_output, },
                                  { .name = NULL}},
};

