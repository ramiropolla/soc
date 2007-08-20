/*
 * Video scaling/colorspace conversion
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
#include "swscale.h"

typedef struct
{
    struct SwsContext *sws;     ///< software scaler context

    /**
     * New dimensions. Special values are:
     *   0 = original width/height
     *  -1 = keep original aspect
     */
    int w, h;
} ScaleContext;

static int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    ScaleContext *scale = ctx->priv;

    /* default to no scaling */
    scale->w =
    scale->h = 0;

    if(args)
        sscanf(args, "%d:%d", &scale->w, &scale->h);

    /* sanity check parms */
    if(scale->w <  -1 || scale->h <  -1)
        return -1;
    if(scale->w == -1 && scale->h == -1)
        scale->w =
        scale->h = 0;

    return 0;
}

static void uninit(AVFilterContext *ctx)
{
    ScaleContext *scale = ctx->priv;
    if(scale->sws)
        sws_freeContext(scale->sws);
}

static int *query_formats(AVFilterLink *link)
{
    return avfilter_make_format_list(31,
                PIX_FMT_YUV444P,  PIX_FMT_YUV422P,  PIX_FMT_YUV420P,
                PIX_FMT_YUV411P,  PIX_FMT_YUV410P,
                PIX_FMT_YUYV422,  PIX_FMT_UYVY422,  PIX_FMT_UYYVYY411,
                PIX_FMT_YUVJ444P, PIX_FMT_YUVJ422P, PIX_FMT_YUVJ420P,
                PIX_FMT_YUV440P,  PIX_FMT_YUVJ440P,
                PIX_FMT_RGB32,    PIX_FMT_BGR32,
                PIX_FMT_RGB32_1,  PIX_FMT_BGR32_1,
                PIX_FMT_RGB24,    PIX_FMT_BGR24,
                PIX_FMT_RGB565,   PIX_FMT_BGR565,
                PIX_FMT_RGB555,   PIX_FMT_BGR555,
                PIX_FMT_RGB8,     PIX_FMT_BGR8,
                PIX_FMT_RGB4_BYTE,PIX_FMT_BGR4_BYTE,
                PIX_FMT_GRAY16BE, PIX_FMT_GRAY16LE,
                PIX_FMT_GRAY8,    PIX_FMT_PAL8);
}

static int config_props(AVFilterLink *link)
{
    ScaleContext *scale = link->src->priv;
    int w, h;

    if(scale->sws)
        sws_freeContext(scale->sws);

    w = scale->w;
    h = scale->h;
    if(!w)      w = link->src->inputs[0]->w;
    if(!h)      h = link->src->inputs[0]->h;
    if(w == -1) w = scale->h*link->src->inputs[0]->w/link->src->inputs[0]->h;
    if(h == -1) h = scale->w*link->src->inputs[0]->h/link->src->inputs[0]->w;

    /* TODO: make algorithm configurable */
    scale->sws = sws_getContext(link->src->inputs[0]->w,
                                link->src->inputs[0]->h,
                                link->src->inputs[0]->format,
                                w, h, link->format, SWS_BILINEAR,
                                NULL, NULL, NULL);

    link->w = w;
    link->h = h;

    return !scale->sws;
}

/* TODO: figure out the swscale API well enough to scale slice at a time */
static void end_frame(AVFilterLink *link)
{
    ScaleContext *scale = link->dst->priv;

    sws_scale(scale->sws, link->cur_pic->data, link->cur_pic->linesize, 0,
              link->cur_pic->h, link->dst->outputs[0]->outpic->data,
              link->dst->outputs[0]->outpic->linesize);
    avfilter_draw_slice(link->dst->outputs[0], 0, link->dst->outputs[0]->h);
    avfilter_end_frame(link->dst->outputs[0]);

    avfilter_unref_pic(link->cur_pic);
    avfilter_unref_pic(link->dst->outputs[0]->outpic);
}

static void draw_slice(AVFilterLink *link, int y, int h)
{
}

AVFilter avfilter_vf_scale =
{
    .name      = "scale",
    .author    = "Bobby Bingham",

    .init      = init,
    .uninit    = uninit,

    .priv_size = sizeof(ScaleContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .draw_slice      = draw_slice,
                                    .end_frame       = end_frame,
                                    .query_formats   = query_formats,
                                    .min_perms       = AV_PERM_READ, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .query_formats   = query_formats,
                                    .config_props    = config_props, },
                                  { .name = NULL}},
};

