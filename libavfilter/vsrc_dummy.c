/*
 * Dummy video source filter for testing
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

#define RED     0x20
#define GREEN   0x60
#define BLUE    0xC0

typedef struct {
    int64_t pts;
} DummyContext;

static int *query_formats(AVFilterLink *link)
{
    return avfilter_make_format_list(1, PIX_FMT_RGB24);
}

static int config_props(AVFilterLink *link)
{
    link->w = 640;
    link->h = 480;

    return 0;
}

static int request_frame(AVFilterLink *link)
{
    DummyContext *ctx = link->src->priv;
    AVFilterPicRef *pic;

    int x, y;
    uint8_t *row, *cur;

    pic = avfilter_get_video_buffer(link, AV_PERM_WRITE);
    pic->pts  =
    ctx->pts += 30;
    avfilter_start_frame(link, avfilter_ref_pic(pic,
                                                link->dst->outputs[0]->dst,
                                                ~0));

    row = pic->data[0];
    for(y = 0; y < pic->h; y ++) {
        cur = row;
        for(x = 0; x < pic->w; x ++) {
            *cur ++ = BLUE;
            *cur ++ = GREEN;
            *cur ++ = RED;
        }
        row += pic->linesize[0];
    }

    avfilter_draw_slice(link, 0, pic->h);

    avfilter_end_frame(link);
    avfilter_unref_pic(pic);

    return 0;
}

AVFilter vsrc_dummy =
{
    .name      = "dummy",
    .author    = "Bobby Bingham",
    .priv_size = sizeof(DummyContext),

    .inputs    = (AVFilterPad[]) {{ .name = NULL }},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .request_frame   = request_frame,
                                    .query_formats   = query_formats,
                                    .config_props    = config_props, },
                                  { .name = NULL}},
};

