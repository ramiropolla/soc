/*
 * Loading filter graph descriptions from files
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
#include "avfiltergraph.h"

#define LINESIZE    240             ///< maximum length of an input line

AVFilterGraphDesc *avfilter_graph_load_desc(const char *filename)
{
    AVFilterGraphDesc       *ret    = NULL;
    AVFilterGraphDescParser *parser = NULL;

    char line[LINESIZE];
    FILE *in = NULL;

    /* TODO: maybe allow searching in a predefined set of directories to
     * allow users to build up libraries of useful graphs? */
    if(!(in = fopen(filename, "r")))
        goto fail;

    while(fgets(line, LINESIZE, in))
        if(avfilter_graph_parse_desc(&ret, &parser, line) < 0)
            goto fail;

    fclose(in);
    av_free(parser);
    return ret;

fail:
    av_free(ret);
    av_free(parser);
    if(in) fclose(in);

    return NULL;
}

