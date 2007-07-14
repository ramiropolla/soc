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
#include "avfiltergraph.h"

int64_t sdl_display(AVFilterContext *ctx);

int main(int argc, char **argv)
{
    int i;
    int ret = -1;
    int64_t pts = 0, newpts;
    AVFilterGraph   *graph;
    AVFilterContext *out;

    if(argc < 3) {
        av_log(NULL, AV_LOG_ERROR, "require at least two filters\n");
        return -1;
    }

    avfilter_init();
    graph = avfilter_create_graph();
    if(avfilter_graph_load_chain(graph, argc - 1, argv + 1, NULL, &out) < 0)
        goto done;

    while(pts < 5000) {
        newpts = sdl_display(out);
        usleep(newpts - pts);
        pts = newpts;
    }

    ret = 0;

done:
    avfilter_destroy_graph(graph);

    return ret;
}

