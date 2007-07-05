/*
 * Video filter test program
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

#include <unistd.h>

#include "avfilter.h"

void sdl_display(AVFilterContext *ctx);

int main()
{
    int i;
    AVFilterContext *src, *crop, *out;

    avfilter_init();

    src  = avfilter_create_by_name("dummy");
    crop = avfilter_create_by_name("crop");
    out  = avfilter_create_by_name("sdl");

    avfilter_init_filter(src,  NULL);
    avfilter_init_filter(crop, "20:40:320:240");
    avfilter_init_filter(out,  NULL);

    avfilter_link(src,  0, crop, 0);
    avfilter_link(crop, 0, out,  0);

    for(i = 0; i < 10; i ++) {
        sdl_display(out);
        sleep(1);
    }

    avfilter_destroy(src);
    avfilter_destroy(crop);
    avfilter_destroy(out);

    return 0;
}

