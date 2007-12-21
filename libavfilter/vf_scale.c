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

static int query_formats(AVFilterContext *ctx)
{
    AVFilterFormats *formats;

    if(ctx->inputs[0]) {
        formats = avfilter_all_colorspaces();
        avfilter_formats_ref(formats, &ctx->inputs[0]->out_formats);
    }
    if(ctx->outputs[0]) {
        formats = avfilter_all_colorspaces();
        avfilter_formats_ref(formats, &ctx->outputs[0]->in_formats);
    }

    return 0;
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

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    AVFilterLink *out = link->dst->outputs[0];

    out->outpic      = avfilter_get_video_buffer(out, AV_PERM_WRITE);
    out->outpic->pts = picref->pts;

    out->outpic->pixel_aspect.num = picref->pixel_aspect.num * out->h * link->w;
    out->outpic->pixel_aspect.den = picref->pixel_aspect.den * out->w * link->h;
    av_reduce(&out->outpic->pixel_aspect.num, &out->outpic->pixel_aspect.den,
               out->outpic->pixel_aspect.num,  out->outpic->pixel_aspect.den,
         FFMAX(out->outpic->pixel_aspect.num,  out->outpic->pixel_aspect.den));

    avfilter_start_frame(out, avfilter_ref_pic(out->outpic, ~0));
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

AVFilter avfilter_vf_scale =
{
    .name      = "scale",
    .author    = "Bobby Bingham",

    .init      = init,
    .uninit    = uninit,

    .query_formats = query_formats,

    .priv_size = sizeof(ScaleContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .start_frame     = start_frame,
                                    .end_frame       = end_frame,
                                    .min_perms       = AV_PERM_READ, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .config_props    = config_props, },
                                  { .name = NULL}},
};

