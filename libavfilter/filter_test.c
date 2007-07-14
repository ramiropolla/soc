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
#include <string.h>

#include "avfilter.h"

int64_t sdl_display(AVFilterContext *ctx);

AVFilterContext *create_filter(char *argv)
{
    AVFilterContext *ret;
    char *name, *args;

    name = argv;
    if((args = strchr(argv, '='))) {
        /* ensure we at least have a name */
        if(args == argv)
            return NULL;

        *args ++ = 0;
    }

    av_log(NULL, AV_LOG_INFO, "creating filter \"%s\" with args \"%s\"\n",
           name, args ? args : "(none)");

    if((ret = avfilter_create_by_name(name, NULL))) {
        if(avfilter_init_filter(ret, args)) {
            av_log(NULL, AV_LOG_ERROR, "error initializing filter!\n");
            avfilter_destroy(ret);
            ret = NULL;
        }
    } else av_log(NULL, AV_LOG_ERROR, "error creating filter!\n");

    return ret;
}

int main(int argc, char **argv)
{
    int i;
    int ret = -1;
    int64_t pts = 0, newpts;
    AVFilterGraph   *graph;
    AVFilterContext *filters[2];

    if(argc < 3) {
        av_log(NULL, AV_LOG_ERROR, "require at least two filters\n");
        return -1;
    }

    avfilter_init();
    graph = avfilter_create_graph();

    for(i = 0; i < argc-1; i ++) {
        if(!(filters[1] = create_filter(argv[i+1])))
            goto done;
        avfilter_graph_add_filter(graph, filters[1]);
        if(i && avfilter_link(filters[0], 0, filters[1], 0)) {
            av_log(NULL, AV_LOG_ERROR, "error linking filters!\n");
            goto done;
        }
        filters[0] = filters[1];
    }

    while(pts < 5000) {
        newpts = sdl_display(filters[1]);
        usleep(newpts - pts);
        pts = newpts;
    }

    ret = 0;

done:
    avfilter_destroy_graph(graph);

    return ret;
}

