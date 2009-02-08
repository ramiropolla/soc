/*
 * video slicing filter
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

typedef struct {
    int h;          ///< output slice height
    int vshift;     ///< chroma subsampling shift
} SliceContext;

static av_cold int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    SliceContext *slice = ctx->priv;

    /* set slice height */
    slice->h = 16;
    if(args)
        sscanf(args, "%d", &slice->h);

    return 0;
}

static int config_props(AVFilterLink *link)
{
    SliceContext *slice = link->dst->priv;
    int tmp;

    avcodec_get_chroma_sub_sample(link->format, &tmp, &slice->vshift);

    /* ensure that slices play nice with chroma subsampling, and enforce
     * a reasonable minimum size for the slices */
    slice->h = FFMAX(8, slice->h & (-1 << slice->vshift));

    return 0;
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    avfilter_start_frame(link->dst->outputs[0], picref);
}

static void end_frame(AVFilterLink *link)
{
    avfilter_end_frame(link->dst->outputs[0]);
}

static void draw_slice(AVFilterLink *link, int y, int h)
{
    SliceContext *slice = link->dst->priv;
    int y2;

    for(y2 = y; y2 + slice->h <= y + h; y2 += slice->h) {
        avfilter_draw_slice(link->dst->outputs[0], y2, slice->h);
    }

    if(y2 < y + h)
        avfilter_draw_slice(link->dst->outputs[0], y2, y + h - y2);
}

AVFilter avfilter_vf_slicify =
{
    .name      = "slicify",

    .init      = init,

    .priv_size = sizeof(SliceContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .config_props    = config_props,
                                    .end_frame       = end_frame, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO, },
                                  { .name = NULL}},
};

