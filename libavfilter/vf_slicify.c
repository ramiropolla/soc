/*
 * Video slicing filter
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

#include <string.h>
#include <stdio.h>

#include "avfilter.h"

typedef struct {
    int h;          ///< output slice height
    int vshift;     ///< chroma subsampling shift
} SliceContext;

static int init(AVFilterContext *ctx, const char *args, const void *opaque)
{
    SliceContext *slice = ctx->priv;

    /* set slice height */
    slice->h = 16;
    if(args)
        sscanf(args, "%d", &slice->h);

    return 0;
}

static int *query_formats(AVFilterLink *link)
{
    return avfilter_make_format_list(1, PIX_FMT_RGB24);
}

static int config_props(AVFilterLink *link)
{
    SliceContext *slice = link->dst->priv;
    int tmp;

    avcodec_get_chroma_sub_sample(link->format, &tmp, &slice->vshift);

    /* ensure that slices play nice with chroma subsampling, and enforce
     * a reasonable minimum size for the slices */
    slice->h = FFMAX(8, slice->h & (-1 << slice->vshift));
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    avfilter_default_start_frame(link, picref);
    avfilter_start_frame(link->dst->outputs[0], avfilter_ref_pic(picref, ~0));
}

static void end_frame(AVFilterLink *link)
{
    avfilter_default_end_frame(link);
    avfilter_end_frame(link->dst->outputs[0]);
}

static void draw_slice(AVFilterLink *link, uint8_t *data[4], int y, int h)
{
    SliceContext *slice = link->dst->priv;
    uint8_t *src[4];
    int y2, i;

    memcpy(src, data, sizeof(src));

    for(y2 = y; y2 + slice->h <= y + h; y2 += slice->h) {
        avfilter_draw_slice(link->dst->outputs[0], src, y2, slice->h);
        src[0] += link->cur_pic->linesize[0] * slice->h;
        src[3] += link->cur_pic->linesize[3] * slice->h;
        /* TODO: make sure this works once other filters support YUV too */
        for(i = 1; i < 3; i ++)
            src[i] += link->cur_pic->linesize[i] * (slice->h >> slice->vshift);
    }

    if(y2 < y + h)
        avfilter_draw_slice(link->dst->outputs[0], src, y2, y + h - y2);
}

static void request_frame(AVFilterLink *link)
{
    avfilter_request_frame(link->src->inputs[0]);
}

AVFilter vf_slicify =
{
    .name      = "slicify",
    .author    = "Bobby Bingham",

    .init      = init,

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .query_formats   = query_formats,
                                    .config_props    = config_props,
                                    .end_frame       = end_frame, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .request_frame   = request_frame, },
                                  { .name = NULL}},
};

