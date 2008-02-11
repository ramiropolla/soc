/*
 * vertical flip filter
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

#include "avfilter.h"

typedef struct
{
    int vsub;   //< chroma subsampling
} FlipContext;

static int config_input(AVFilterLink *link)
{
    FlipContext *flip = link->dst->priv;
    int tmp;

    avcodec_get_chroma_sub_sample(link->format, &tmp, &flip->vsub);

    return 0;
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    FlipContext *flip = link->dst->priv;
    AVFilterPicRef *ref2 = avfilter_ref_pic(picref, ~0);
    int i;

    ref2->data[0] += (ref2->h-1) * ref2->linesize[0];
    ref2->linesize[0] = -ref2->linesize[0];
    for(i = 1; i < 4; i ++) {
        if(ref2->data[i]) {
            ref2->data[i] += ((ref2->h >> flip->vsub)-1) * ref2->linesize[i];
            ref2->linesize[i] = -ref2->linesize[i];
        }
    }

    avfilter_start_frame(link->dst->outputs[0], ref2);
}

static void draw_slice(AVFilterLink *link, int y, int h)
{
    AVFilterContext *ctx = link->dst;
    AVFilterPicRef *pic = link->cur_pic;

    avfilter_draw_slice(ctx->outputs[0], pic->h - (y+h), h);
}

AVFilter avfilter_vf_vflip =
{
    .name      = "vflip",
    .priv_size = sizeof(FlipContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .config_props    = config_input, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO, },
                                  { .name = NULL}},
};

