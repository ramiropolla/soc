/*
 * Box drawing filter. Also a nice template for a filter that needs to write
 * in the input frame.
 *
 * Copyright (c) 2008 Affine Systems, Inc (Michael Sullivan, Bobby Impollonia)
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
#include <ctype.h>

#include "avfilter.h"

typedef struct
{
    const char *name;
    unsigned char y;
    unsigned char cb;
    unsigned char cr;
} box_color;

#define NUM_COLORS 8
static box_color colors[] =
{
    {.name = "blue",   .y =  40, .cr = 109, .cb = 240},
    {.name = "red",    .y =  81, .cr = 240, .cb =  90},
    {.name = "black",  .y =  16, .cr = 128, .cb = 128},
    {.name = "white",  .y = 255, .cr = 128, .cb = 128},
    {.name = "green",  .y = 144, .cr =  34, .cb =  53},
    {.name = "yellow", .y = 210, .cr = 146, .cb =  16},
    {.name = "gray",   .y = 128, .cr = 128, .cb = 128},
    {.name = "grey",   .y = 128, .cr = 128, .cb = 128}
};

typedef struct
{
    int x, y, w, h;
    box_color color;
    int vsub, hsub;   //< chroma subsampling
} BoxContext;

static void box_strlwr(char *str)
{
    int i = 0;
    while(str[i] != '\0') {
        if(isupper(str[i]))
            str[i] = tolower(str[i]);

        i++;
    }
}

static box_color get_color(char *color)
{
    box_color tmp;
    int i;
    box_strlwr(color);
    for(i = NUM_COLORS - 1; i >= 0; i--) {
        tmp = colors[i];
        if(!strcmp(color, colors[i].name))
            break;
    }
    return tmp;
}

static av_cold int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    BoxContext *context= ctx->priv;
    char tmp[1024];

    if(!args || strlen(args) > 1024) {
        av_log(ctx, AV_LOG_ERROR, "Invalid arguments!\n");
        return -1;
    }

    sscanf(args, "%d:%d:%d:%d:%s", &(context->x), &(context->y),
           &(context->w), &(context->h), tmp);

    context->color = get_color(tmp);
    return 0;
}

static int query_formats(AVFilterContext *ctx)
{
    avfilter_set_common_formats(ctx,
                                avfilter_make_format_list(10,
                PIX_FMT_YUV444P,  PIX_FMT_YUV422P,  PIX_FMT_YUV420P,
                PIX_FMT_YUV411P,  PIX_FMT_YUV410P,
                PIX_FMT_YUVJ444P, PIX_FMT_YUVJ422P, PIX_FMT_YUVJ420P,
                PIX_FMT_YUV440P,  PIX_FMT_YUVJ440P)
        );
    return 0;
}

static int config_input(AVFilterLink *link)
{
    BoxContext *context = link->dst->priv;

    avcodec_get_chroma_sub_sample(link->format,
                                  &context->hsub, &context->vsub);

    return 0;
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    avfilter_start_frame(link->dst->outputs[0], picref);
}


static void draw_box(AVFilterPicRef *pic, BoxContext* context, box_color color)
{
    int x, y;
    int channel;
    unsigned char *row[4];
    int xb = context->x;
    int yb = context->y;

    for (y = yb; (y < yb + context->h) && y < pic->h; y++) {
        row[0] = pic->data[0] + y * pic->linesize[0];

        for (channel = 1; channel < 3; channel++)
            row[channel] = pic->data[channel] +
                pic->linesize[channel] * (y>> context->vsub);

        for (x = xb; (x < xb + context->w) && x < pic->w; x++)
            if((y - yb < 3) || (yb + context->h - y < 4) ||
               (x - xb < 3) || (xb + context->w - x < 4)) {
                row[0][x] = color.y;
                row[1][x >> context->hsub] = color.cb;
                row[2][x >> context->hsub] = color.cr;
            }
    }
}

static void end_frame(AVFilterLink *link)
{
    BoxContext *context = link->dst->priv;
    AVFilterLink *output = link->dst->outputs[0];
    AVFilterPicRef *pic = link->cur_pic;

    draw_box(pic,context,context->color);

    avfilter_draw_slice(output, 0, pic->h);
    avfilter_end_frame(output);
}

AVFilter avfilter_vf_drawbox=
{
    .name      = "drawbox",
    .priv_size = sizeof(BoxContext),
    .init      = init,

    .query_formats   = query_formats,
    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .end_frame       = end_frame,
                                    .config_props    = config_input,
                                    .min_perms       = AV_PERM_WRITE |
                                                       AV_PERM_READ,
                                    .rej_perms       = AV_PERM_REUSE |
                                                       AV_PERM_REUSE2},
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO, },
                                  { .name = NULL}},
};
