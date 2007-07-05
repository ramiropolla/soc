/*
 * SDL video output filter
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

/* TODO: right now only rgb24 supported for testing.  add more formats */

#include <SDL/SDL.h>

#include "avfilter.h"

typedef struct
{
    SDL_Surface *surface;
} SDLContext;

static int init(AVFilterContext *ctx, const char *args)
{
    SDLContext *sdl = ctx->priv;

    if(SDL_Init(SDL_INIT_VIDEO | SDL_INIT_NOPARACHUTE)) {
        av_log(ctx, AV_LOG_FATAL, "unable to initialize SDL: %s\n",
               SDL_GetError());
        return -1;
    }

    return 0;
}

static int *query_formats(AVFilterLink *link)
{
    return avfilter_make_format_list(1, PIX_FMT_RGB24);
}

static int config_props(AVFilterLink *link)
{
    SDLContext *sdl = link->dst->priv;

    /* SDL docs claim that calling this a second time to resize the
     * surface will automatically free the old surface for us */
    sdl->surface = SDL_SetVideoMode(link->w, link->h, 24,
                                    SDL_HWSURFACE | SDL_DOUBLEBUF);

    return !sdl->surface;
}

static void uninit(AVFilterContext *ctx)
{
    SDL_QuitSubSystem(SDL_INIT_VIDEO);
}

static void start_frame(AVFilterLink *link, AVFilterPicRef *picref)
{
    SDLContext *sdl = link->dst->priv;
    avfilter_default_start_frame(link, picref);
    SDL_LockSurface(sdl->surface);

    av_log(link->dst, AV_LOG_INFO, "start_frame()\n");
}

static void draw_slice(AVFilterLink *link, uint8_t *data[4], int y, int h)
{
    AVFilterPicRef *picref = link->cur_pic;
    SDLContext *sdl = link->dst->priv;

    uint8_t *src = data[0];
    uint8_t *dst = &((uint8_t *)sdl->surface->pixels)[y * sdl->surface->pitch];
    int copysize = 3 * link->w;

    int i;

    for(i = 0; i < h; i ++) {
        memcpy(dst, src, copysize);
        src += picref->linesize[0];
        dst += sdl->surface->pitch;
    }

    av_log(link->dst, AV_LOG_INFO, "draw_slice()\n");
}

static void end_frame(AVFilterLink *link)
{
    SDLContext *sdl = link->dst->priv;
    SDL_UnlockSurface(sdl->surface);
    avfilter_default_end_frame(link);

    av_log(link->dst, AV_LOG_INFO, "end_frame()\n");
}

/* XXX: this is a hack.  should provide a proper vout interface */
void sdl_display(AVFilterContext *ctx)
{
    SDLContext *sdl = ctx->priv;

    SDL_Flip(sdl->surface);
    avfilter_request_frame(ctx->inputs[0]);
}

AVFilter vo_sdl =
{
    .name      = "sdl",
    .author    = "Bobby Bingham",
    .priv_size = sizeof(SDLContext),

    .init      = init,
    .uninit    = uninit,

    .inputs    = (AVFilterPad[]) {{ .name          = "default",
                                    .type          = AV_PAD_VIDEO,
                                    .start_frame   = start_frame,
                                    .end_frame     = end_frame,
                                    .draw_slice    = draw_slice,
                                    .query_formats = query_formats,
                                    .config_props  = config_props, },
                                  { .name = NULL}},
    .outputs   = (AVFilterPad[]) {{ .name = NULL}},
};

