/*
 * Video framerate modification filter
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

/* TODO: improve handling of non-continuous timestamps (mpeg, seeking, etc) */

#include <stdio.h>

#include "avfilter.h"

typedef struct {
    uint64_t timebase;
    uint64_t pts;
    AVFilterPicRef *pic;
} FPSContext;

static int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    FPSContext *fps = ctx->priv;
    int framerate;

    /* TODO: support framerates specified as decimals or fractions */
    if(args && sscanf(args, "%d", &framerate))
        fps->timebase = 1000 / framerate;
    else
        /* default to 25 fps */
        fps->timebase = 1000 / 25;

    return 0;
}

static void uninit(AVFilterContext *ctx)
{
    FPSContext *fps = ctx->priv;
    if(fps->pic) avfilter_unref_pic(fps->pic);
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

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    FPSContext *fps = link->dst->priv;
    if(fps->pic) avfilter_unref_pic(fps->pic);
    fps->pic = picref;
}

static void end_frame(AVFilterLink *link)
{
}

static void draw_slice(AVFilterLink *link, int y, int h)
{
}

static int request_frame(AVFilterLink *link)
{
    FPSContext *fps = link->src->priv;

    while(!fps->pic || fps->pic->pts < fps->pts)
        if(avfilter_request_frame(link->src->inputs[0]))
            return -1;

    avfilter_start_frame(link, avfilter_ref_pic(fps->pic, ~AV_PERM_WRITE));
    avfilter_draw_slice (link, 0, fps->pic->h);
    avfilter_end_frame  (link);

    fps->pts += fps->timebase;

    return 0;
}

AVFilter vf_fps =
{
    .name      = "fps",
    .author    = "Bobby Bingham",

    .init      = init,
    .uninit    = uninit,

    .priv_size = sizeof(FPSContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .query_formats   = query_formats,
                                    .end_frame       = end_frame, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .request_frame   = request_frame, },
                                  { .name = NULL}},
};

