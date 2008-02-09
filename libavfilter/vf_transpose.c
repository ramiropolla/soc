/*
 * transpose (line => column) video filter
 * Copyright (c) 2008 Vitor Sessak
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
 * @file vf_transpose.c
 * Transposition filter
 *
 * @todo handle packed pixel formats
 */

#include "avfilter.h"

typedef struct
{
    int hsub, vsub;
} TransContext;

static int config_props_input(AVFilterLink *link)
{
    TransContext *trans = link->dst->priv;

    avcodec_get_chroma_sub_sample(link->format, &trans->hsub, &trans->vsub);

    return 0;
}

static int query_formats(AVFilterContext *ctx)
{
    avfilter_set_common_formats(ctx,
        avfilter_make_format_list(10,
                PIX_FMT_YUV444P,  PIX_FMT_YUV422P,  PIX_FMT_YUV420P,
                PIX_FMT_YUV411P,  PIX_FMT_YUV410P,
                PIX_FMT_YUVJ444P, PIX_FMT_YUVJ422P, PIX_FMT_YUVJ420P,
                PIX_FMT_YUV440P,  PIX_FMT_YUVJ440P));
    return 0;
}

static int config_props_output(AVFilterLink *link)
{
    link->w = link->src->inputs[0]->h;
    link->h = link->src->inputs[0]->w;

    return 0;
}

static void draw_slice(AVFilterLink *link, int y, int h)
{
    TransContext *trans = link->dst->priv;
    AVFilterPicRef *in  = link->cur_pic;
    AVFilterPicRef *out = link->dst->outputs[0]->outpic;
    int i, j, plane;

    /* luma plane */
    for(i = y; i < h; i ++)
        for(j = 0; j < link->w; j ++)
            *(out->data[0] +   j *out->linesize[0] + i) =
                *(in->data[0]+ i * in->linesize[0] + j);

    /* chroma planes */
    for(plane = 1; plane < 3; plane ++) {
        for(i = y >> trans->vsub; i < h >> trans->vsub; i++) {
            for(j = 0; j < link->w >> trans->hsub; j++)
                *(out->data[plane] +   j *out->linesize[plane] + i) =
                    *(in->data[plane]+ i * in->linesize[plane] + j);
        }
    }

    avfilter_draw_slice(link->dst->outputs[0], y, h);
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    AVFilterLink *out = link->dst->outputs[0];

    out->outpic      = avfilter_get_video_buffer(out, AV_PERM_WRITE);
    out->outpic->pts = picref->pts;

    out->outpic->pixel_aspect.num = picref->pixel_aspect.den;
    out->outpic->pixel_aspect.den = picref->pixel_aspect.num;

    avfilter_start_frame(out, avfilter_ref_pic(out->outpic, ~0));
}

AVFilter avfilter_vf_transpose =
{
    .name      = "transpose",
    .author    = "Vitor Sessak",

    .priv_size = sizeof(TransContext),

    .query_formats = query_formats,

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .config_props    = config_props_input,
                                    .min_perms       = AV_PERM_READ, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .config_props    = config_props_output,
                                    .type            = AV_PAD_VIDEO, },
                                  { .name = NULL}},
};

