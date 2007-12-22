/*
 * Filter graph descriptions
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

#include <ctype.h>
#include <string.h>

#include "avfilter.h"
#include "avfiltergraph.h"

#define LINESIZE    240             ///< maximum length of an input line

/** a comment is a line which is empty, or starts with whitespace, ';' or '#' */
static inline int is_line_comment(char *line)
{
    return line[0] == 0     ||
           isspace(line[0]) ||
           line[0] == ';'   ||
           line[0] == '#';
}

static AVFilterGraphDescSection parse_section_name(char *line)
{
    struct {
        char *str;
        int section;
    } *sec, sections[] = { { "[filters]", SEC_FILTERS },
                           { "[links]",   SEC_LINKS   },
                           { "[inputs]",  SEC_INPUTS  },
                           { "[outputs]", SEC_OUTPUTS },
                           { NULL, 0 } };

    for(sec = sections; sec->str; sec ++)
        if(!strncmp(line, sec->str, strlen(sec->str)))
            return sec->section;

    av_log(NULL, AV_LOG_ERROR, "unknown section name in graph description\n");
    return SEC_NONE;
}

static AVFilterGraphDescFilter *parse_filter(char *line)
{
    AVFilterGraphDescFilter *ret = av_mallocz(sizeof(AVFilterGraphDescFilter));
    char *tok;

    if(!(tok = strchr(line, '='))) {
        av_log(NULL, AV_LOG_ERROR, "filter line missing type of filter");
        goto fail;
    }
    *tok = '\0';
    ret->name = av_strdup(line);
    line = tok+1;

    if((tok = strchr(line, '='))) {
        *tok ++ = '\0';
        ret->args = av_strdup(tok);
    }
    ret->filter = av_strdup(line);

    return ret;

fail:
    av_free(ret->name);
    av_free(ret->filter);
    av_free(ret->args);
    av_free(ret);
    return NULL;
}

/* TODO: allow referencing pad names, not just indices */
static AVFilterGraphDescLink *parse_link(char *line)
{
    AVFilterGraphDescLink *ret = av_mallocz(sizeof(AVFilterGraphDescLink));
    ret->src = av_malloc(32);
    ret->dst = av_malloc(32);

    if(sscanf(line, "%31[a-zA-Z0-9]:%u=%31[a-zA-Z0-9]:%u",
              ret->src, &ret->srcpad, ret->dst, &ret->dstpad) < 4) {
        av_free(ret->src);
        av_free(ret->dst);
        av_free(ret);
        return NULL;
    }

    return ret;
}

/* TODO: allow referencing pad names, not just indices */
static AVFilterGraphDescExport *parse_export(char *line)
{
    AVFilterGraphDescExport *ret = av_mallocz(sizeof(AVFilterGraphDescLink));
    ret->name   = av_malloc(32);
    ret->filter = av_malloc(32);

    if(sscanf(line, "%31[a-zA-Z0-9]=%31[a-zA-Z0-9]:%u",
              ret->name, ret->filter, &ret->pad) < 3) {
        av_free(ret->name);
        av_free(ret->filter);
        av_free(ret);
        return NULL;
    }

    return ret;
}

int avfilter_graph_parse_desc(AVFilterGraphDesc **desc,
                              AVFilterGraphDescParser **parser,
                              char *line)
{
    void *next;
    int len;

    if(!*desc)
        if(!(*desc = av_mallocz(sizeof(AVFilterGraphDesc))))
            return AVERROR_NOMEM;

    if(!*parser) {
        if(!(*parser = av_mallocz(sizeof(AVFilterGraphDescParser))))
            return AVERROR_NOMEM;

        (*parser)->filterp = &(*desc)->filters;
        (*parser)->linkp   = &(*desc)->links;
        (*parser)->inputp  = &(*desc)->inputs;
        (*parser)->outputp = &(*desc)->outputs;
    }

    /* ignore comments */
    if(is_line_comment(line)) return 0;

    /* check if a new section is starting */
    if(line[0] == '[') {
        if(((*parser)->section = parse_section_name(line)) == SEC_NONE)
            return AVERROR_INVALIDDATA;
        return 0;
    }

    /* remove any trailing newline characters */
    for(len = strlen(line); len && (line[len-1]=='\n'||line[len-1]=='\r');)
        line[--len] = '\0';

    switch((*parser)->section) {
    case SEC_FILTERS:
        if(!(next = parse_filter(line)))
            return AVERROR_INVALIDDATA;
        *(*parser)->filterp = next;
        (*parser)->filterp  = &(*(*parser)->filterp)->next;
        break;
    case SEC_LINKS:
        if(!(next = parse_link(line)))
            return AVERROR_INVALIDDATA;
        *(*parser)->linkp = next;
        (*parser)->linkp  = &(*(*parser)->linkp)->next;
        break;
    case SEC_INPUTS:
        if(!(next = parse_export(line)))
            return AVERROR_INVALIDDATA;
        *(*parser)->inputp = next;
        (*parser)->inputp  = &(*(*parser)->inputp)->next;
        break;
    case SEC_OUTPUTS:
        if(!(next = parse_export(line)))
            return AVERROR_INVALIDDATA;
        *(*parser)->outputp = next;
        (*parser)->outputp  = &(*(*parser)->outputp)->next;
        break;
    default:
        return AVERROR_INVALIDDATA;
    }

    return 0;
}

void avfilter_graph_free_desc(AVFilterGraphDesc *desc)
{
    void *next;

    while(desc->filters) {
        next = desc->filters->next;
        av_free(desc->filters->name);
        av_free(desc->filters->filter);
        av_free(desc->filters->args);
        av_free(desc->filters);
        desc->filters = next;
    }

    while(desc->links) {
        next = desc->links->next;
        av_free(desc->links->src);
        av_free(desc->links->dst);
        av_free(desc->links);
        desc->links = next;
    }

    while(desc->inputs) {
        next = desc->inputs->next;
        av_free(desc->inputs->filter);
        av_free(desc->inputs);
        desc->inputs = next;
    }

    while(desc->outputs) {
        next = desc->outputs->next;
        av_free(desc->outputs->filter);
        av_free(desc->outputs);
        desc->outputs = next;
    }
}

