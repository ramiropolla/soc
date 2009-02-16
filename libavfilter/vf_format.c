/*
 * video (no)format filter
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
    /** nonzero for each format included in the list */
    uint8_t formats[PIX_FMT_NB];
} FormatContext;

#define FMT_NAME_MAXSIZE 32

static av_cold int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    FormatContext *format = ctx->priv;
    const char *cur, *sep;
    char name[FMT_NAME_MAXSIZE];
    int fmt, len;

    /* parse the list of formats */
    for(cur = args; cur; cur = sep ? sep+1 : NULL) {
        if(!(sep = strchr(cur, ':')))
            len = strlen(cur);
        else
            len = sep - cur;
        if(len >= FMT_NAME_MAXSIZE) {
            av_log(ctx, AV_LOG_ERROR, "Format name too long\n");
            return -1;
        }

        memcpy(name, cur, len);
        name[len] = 0;
        fmt = avcodec_get_pix_fmt(name);

        if(fmt == PIX_FMT_NONE) {
            av_log(ctx, AV_LOG_ERROR, "Unknown pixel format: %s\n", name);
            return -1;
        }

        format->formats[fmt] = 1;
    }

    return 0;
}

static AVFilterFormats *make_format_list(FormatContext *format, uint8_t val)
{
    AVFilterFormats *ret;
    int i;

    ret = av_mallocz(sizeof(AVFilterFormats));
    ret->formats = av_malloc(sizeof(int) * PIX_FMT_NB);

    for(i = 0; i < PIX_FMT_NB; i ++)
        if(format->formats[i] == val)
            ret->formats[ret->format_count ++] = i;

    return ret;
}

static int query_formats_noformat(AVFilterContext *ctx)
{
    avfilter_set_common_formats(ctx, make_format_list(ctx->priv, 0));
    return 0;
}

static int query_formats_format(AVFilterContext *ctx)
{
    avfilter_set_common_formats(ctx, make_format_list(ctx->priv, 1));
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
    avfilter_draw_slice(link->dst->outputs[0], y, h);
}

AVFilter avfilter_vf_noformat =
{
    .name      = "noformat",

    .init      = init,

    .query_formats = query_formats_noformat,

    .priv_size = sizeof(FormatContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .end_frame       = end_frame, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO },
                                  { .name = NULL}},
};

AVFilter avfilter_vf_format =
{
    .name      = "format",

    .init      = init,

    .query_formats = query_formats_format,

    .priv_size = sizeof(FormatContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .draw_slice      = draw_slice,
                                    .end_frame       = end_frame, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO },
                                  { .name = NULL}},
};

