/*
 * filter to overlay one video on top of another
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
    int x, y;                   //< position of subpicture

    /** pics[0][0..1] are pictures for the main image.
     *  pics[1][0..1] are pictures for the sub image */
    AVFilterPicRef *pics[2][2];

    int bpp;                    //< bytes per pixel
    int hsub, vsub;             //< chroma subsampling
} OverlayContext;

static int init(AVFilterContext *ctx, const char *args, void *opaque)
{
    OverlayContext *over = ctx->priv;

    if(!args || sscanf(args, "%d:%d", &over->x, &over->y) != 2) {
        over->x =
        over->y = 0;
    }

    return 0;
}

static void uninit(AVFilterContext *ctx)
{
    OverlayContext *over = ctx->priv;
    int i, j;

    for(i = 0; i < 2; i ++)
        for(j = 0; j < 2; j ++)
            if(over->pics[i][j])
                avfilter_unref_pic(over->pics[i][j]);
}

static int config_input_main(AVFilterLink *link)
{
    OverlayContext *over = link->dst->priv;

    switch(link->format) {
    case PIX_FMT_RGB32:
    case PIX_FMT_BGR32:
        over->bpp = 4;
        break;
    case PIX_FMT_RGB24:
    case PIX_FMT_BGR24:
        over->bpp = 3;
        break;
    case PIX_FMT_RGB565:
    case PIX_FMT_RGB555:
    case PIX_FMT_BGR565:
    case PIX_FMT_BGR555:
    case PIX_FMT_GRAY16BE:
    case PIX_FMT_GRAY16LE:
        over->bpp = 2;
        break;
    default:
        over->bpp = 1;
    }

    avcodec_get_chroma_sub_sample(link->format, &over->hsub, &over->vsub);

    return 0;
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    OverlayContext *over = link->dst->priv;
    if(over->pics[link->dstpad][0])
        avfilter_unref_pic(over->pics[link->dstpad][0]);
    over->pics[link->dstpad][0] = over->pics[link->dstpad][1];
    over->pics[link->dstpad][1] = picref;
}

static void end_frame(AVFilterLink *link)
{
}

static int lower_timestamp(OverlayContext *over)
{
    if(!over->pics[0][0] &&
       !over->pics[1][0]) return 2;
    if(!over->pics[0][1]) return 0;
    if(!over->pics[1][1]) return 1;

    if(over->pics[0][1]->pts == over->pics[1][1]->pts) return 2;
    return (over->pics[0][1]->pts > over->pics[1][1]->pts);
}

static void copy_image(AVFilterPicRef *dst, int x, int y,
                       AVFilterPicRef *src, int w, int h,
                       int bpp, int hsub, int vsub)
{
    AVPicture pic;
    int i;

    memcpy(&pic, &dst->data, sizeof(AVPicture));
    pic.data[0] += x * bpp;
    pic.data[0] += y * pic.linesize[0];
    for(i = 1; i < 4; i ++) {
        if(pic.data[i]) {
            pic.data[i] +=  x >> hsub;
            pic.data[i] += (y >> vsub) * pic.linesize[i];
        }
    }

    av_picture_copy(&pic, (AVPicture *)src->data, dst->pic->format, w, h);
}

static int request_frame(AVFilterLink *link)
{
    AVFilterPicRef *pic;
    OverlayContext *over = link->src->priv;
    int idx;
    int x, y, w, h;

    /* the first time through, we need to pull a couple frames */
    if(!over->pics[0][1] && !over->pics[1][1] &&
       (avfilter_request_frame(link->src->inputs[0]) ||
        avfilter_request_frame(link->src->inputs[1])))
        return -1;

    /* then we pull a frame from the stream which currently has a lower pts */
    if((idx = lower_timestamp(over)) == 2) {
        if(avfilter_request_frame(link->src->inputs[0]) ||
           avfilter_request_frame(link->src->inputs[1]))
            return -1;
    } else
        if(avfilter_request_frame(link->src->inputs[idx]))
            return -1;

    /* we draw the output frame */
    pic = avfilter_get_video_buffer(link, AV_PERM_WRITE);
    if(over->pics[0][0]) {
        pic->pixel_aspect = over->pics[0][0]->pixel_aspect;
        copy_image(pic, 0, 0, over->pics[0][0], link->w, link->h,
                   over->bpp, over->hsub, over->vsub);
    }
    x = FFMIN(over->x, link->w-1);
    y = FFMIN(over->y, link->h-1);
    w = FFMIN(link->w-x, over->pics[1][0]->w);
    h = FFMIN(link->h-y, over->pics[1][0]->h);
    if(over->pics[1][0])
        copy_image(pic, x, y, over->pics[1][0], w, h,
                   over->bpp, over->hsub, over->vsub);

    /* we give the output frame the higher of the two current pts values */
    pic->pts = FFMAX(over->pics[0][0]->pts, over->pics[1][0]->pts);

    /* and send it to the next filter */
    avfilter_start_frame(link, avfilter_ref_pic(pic, ~0));
    avfilter_draw_slice (link, 0, pic->h);
    avfilter_end_frame  (link);
    avfilter_unref_pic(pic);

    return 0;
}

AVFilter avfilter_vf_overlay =
{
    .name      = "overlay",

    .init      = init,
    .uninit    = uninit,

    .priv_size = sizeof(OverlayContext),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .config_props    = config_input_main,
                                    .end_frame       = end_frame,
                                    .min_perms       = AV_PERM_READ,
                                    .rej_perms       = AV_PERM_REUSE2, },
                                  { .name            = "sub",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .start_frame     = start_frame,
                                    .end_frame       = end_frame,
                                    .min_perms       = AV_PERM_READ,
                                    .rej_perms       = AV_PERM_REUSE2, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_VIDEO,
                                    .request_frame   = request_frame, },
                                  { .name = NULL}},
};

