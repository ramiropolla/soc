/*
 * RGB <-> BGR conversion filter
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

static int *query_in_formats(AVFilterLink *link)
{
    return avfilter_make_format_list(2, PIX_FMT_RGB24, PIX_FMT_BGR24);
}

static int *query_out_formats(AVFilterLink *link)
{
    int format;

    if(link->src->inputs[0]->format == PIX_FMT_RGB24)
        format = PIX_FMT_BGR24;
    else
        format = PIX_FMT_RGB24;

    return avfilter_make_format_list(1, format);
}

static void draw_slice(AVFilterLink *link, uint8_t *data[4], int y, int h)
{
    AVFilterPicRef *outpic = link->dst->outputs[0]->outpic;
    uint8_t *out[4];
    uint8_t *row[2], *cur[2];
    int i, j;

    row[0] = data[0];
    row[1] = out[0] = &outpic->data[0][y * outpic->linesize[0]];
    out[1] = out[2] = out[3] = 0;
    for(i = 0; i < h; i ++) {
        cur[0] = row[0];
        cur[1] = row[1];
        for(j = 0; j < link->cur_pic->w; j ++) {
            cur[1][0] = cur[0][2];
            cur[1][1] = cur[0][1];
            cur[1][2] = cur[0][0];

            cur[0] += 3;
            cur[1] += 3;
        }
        row[0] += link->cur_pic->linesize[0];
        row[1] += outpic->       linesize[0];
    }
    avfilter_draw_slice(link->dst->outputs[0], out, y, h);
}

AVFilter vf_rgb2bgr =
{
    .name      = "rgb2bgr",
    .author    = "Bobby Bingham",

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .draw_slice      = draw_slice,
                                    .query_formats   = query_in_formats, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .query_formats   = query_out_formats, },
                                  { .name = NULL}},
};

