/*
 * PPM file video source filter for testing
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
    int w, h;
    FILE *in;
    AVFilterPicRef *pic;
} PPMContext;

static int init(AVFilterContext *ctx, const char *args)
{
    PPMContext *ppm = ctx->priv;
    FILE *in;
    int max, x, y;

    if(!args) return -1;
    if(!(ppm->in = fopen(args, "r"))) return -1;

    if(fscanf(ppm->in, "P6 %d %d %d ", &ppm->w, &ppm->h, &max) < 3)
        goto fail;
    if(ppm->w < 1 || ppm->h < 1 || max != 255)
        goto fail;

    return 0;

fail:
    fclose(ppm->in);
    ppm->in = NULL;
    return -1;
}

static int *query_formats(AVFilterLink *link)
{
    return avfilter_make_format_list(1, PIX_FMT_BGR24);
}

static int config_props(AVFilterLink *link)
{
    PPMContext *ppm = link->src->priv;

    link->w = ppm->w;
    link->h = ppm->h;

    return 0;
}

static void request_frame(AVFilterLink *link)
{
    PPMContext *ppm = link->src->priv;
    AVFilterPicRef *out;

    int x, y;
    uint8_t *row, *cur;

    if(!ppm->pic) {
        ppm->pic = avfilter_get_video_buffer(link,
                        AV_PERM_WRITE | AV_PERM_PRESERVE | AV_PERM_REUSE);

        row = ppm->pic->data[0];
        for(y = 0; y < ppm->h; y ++) {
            fread(row, 3, ppm->w, ppm->in);
            row += ppm->pic->linesize[0];
        }

        fclose(ppm->in);
        ppm->in = NULL;
    }

    out = avfilter_ref_pic(ppm->pic, ~AV_PERM_WRITE);
    avfilter_start_frame(link, out);
    avfilter_draw_slice(link, out->data, 0, out->h);
    avfilter_end_frame(link);
}

static void uninit(AVFilterContext *ctx)
{
    PPMContext *ppm = ctx->priv;

    if(ppm->pic)
        avfilter_unref_pic(ppm->pic);
    if(ppm->in)
        fclose(ppm->in);
}

AVFilter vsrc_ppm =
{
    .name      = "ppm",
    .author    = "Bobby Bingham",
    .priv_size = sizeof(PPMContext),

    .init      = init,
    .uninit    = uninit,

    .inputs    = (AVFilterPad[]) {{ .name = NULL }},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = AV_PAD_VIDEO,
                                    .request_frame   = request_frame,
                                    .query_formats   = query_formats,
                                    .config_props    = config_props, },
                                  { .name = NULL}},
};

