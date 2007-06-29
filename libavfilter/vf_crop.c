/*
 * Video crop filter
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

typedef struct
{
    int x, y, w, h;
} CropContext;

static int init(AVFilterContext *ctx)
{
    CropContext *crop = ctx->priv;

    if(ctx->inputs[0]->format != PIX_FMT_RGB24) {
        av_log(ctx, AV_LOG_FATAL, "unsupported input format\n");
        return -1;
    }

    crop->x = 20;
    crop->y = 15;
    crop->w = 467;
    crop->h = 45;

    return 0;
}

static int set_video_props(AVFilterLink *link)
{
    CropContext *crop = link->src->priv;

    link->w = crop->w;
    link->h = crop->h;
    link->format = link->src->inputs[0]->format;

    return 0;
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    CropContext *crop = link->dst->priv;
    AVFilterPicRef *ref2 = avfilter_ref_pic(picref);

    ref2->w = crop->w;
    ref2->h = crop->h;
    ref2->data[0] += crop->y * ref2->pic->linesize[0];
    ref2->data[0] += 3 * crop->x;

    av_log(link->dst, AV_LOG_INFO, "start_frame()\n");
    avfilter_default_start_frame(link, picref);

    avfilter_start_frame(link->dst->outputs[0], ref2);
}

static void end_frame(AVFilterLink *link)
{
    avfilter_default_end_frame(link);

    av_log(link->dst, AV_LOG_INFO, "end_frame()\n");
    avfilter_end_frame(link->dst->outputs[0]);
}

static void draw_slice(AVFilterLink *link, uint8_t *data[4], int y, int h)
{
    AVFilterContext *ctx = link->dst;
    AVFilterPic *pic = link->cur_pic->pic;
    CropContext *crop = ctx->priv;

    uint8_t *src[4];
    int top = y;
    int height = h;

    av_log(link->dst, AV_LOG_INFO, "draw_slice()\n");

    if(y >= crop->y + crop->h || y + h <= crop->y) return;

    memcpy(src, data, sizeof(uint8_t *) * 4);

    if(top < crop->y) {
        height -=  crop->y - top;
        src[0] += (crop->y - top) * pic->linesize[0];
        top     =  crop->y;
    }
    if(top + height > crop->y + crop->h)
        height = crop->y + crop->h - top;
    src[0] += 3 * crop->x;

    avfilter_draw_slice(ctx->outputs[0], src, top - crop->y, height);
}

/* XXX: maybe make the default implementation do this? */
static void request_frame(AVFilterLink *link)
{
    avfilter_request_frame(link->src->inputs[0]);
}

AVFilter vf_crop =
{
    .name      = "crop",
    .author    = "Bobby Bingham",
    .priv_size = sizeof(CropContext),

    .init      = init,

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .end_frame       = end_frame, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .request_frame   = request_frame,
                                    .set_video_props = set_video_props},
                                  { .name = NULL}},
};

