/*
 * filter graph descriptions
 * copyright (c) 2008 Vitor Sessak
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


/**
 * For use in av_log
 */
static const char *log_name(void *p)
{
    return "Filter parser";
}

static const AVClass filter_parser_class = {
    "Filter parser",
    log_name
};

static const AVClass *log_ctx = &filter_parser_class;

static void consume_whitespace(const char **buf)
{
    *buf += strspn(*buf, " \n\t");
}

/**
 * get the next non-whitespace char
 */
static char consume_char(const char **buf)
{
    char out;
    consume_whitespace(buf);

    out = **buf;

    if (out)
        (*buf)++;

    return out;
}

/**
 * remove the quotation marks from a string. Ex: "aaa'bb'cc" -> "aaabbcc"
 */
static void unquote(char *str)
{
    char *p1, *p2;
    p1=p2=str;
    while (*p1 != 0) {
        if (*p1 != '\'')
            *p2++ = *p1;
        p1++;
    }

    *p2 = 0;
}

/**
 * Consumes a string from *buf.
 * @return a copy of the consumed string, which should be free'd after use
 */
static char *consume_string(const char **buf)
{
    const char *start;
    char *ret;
    int size;

    consume_whitespace(buf);

    if (!(**buf))
        return av_mallocz(1);

    start = *buf;

    *buf += strcspn(*buf, " ()=,'");

    if (**buf == '\'') {
        char *p = strchr(*buf + 1, '\'');
        if (p)
            *buf = p + 1;
        else
            *buf += strlen(*buf); // Move the pointer to the null end byte
    }

    size = *buf - start + 1;
    ret = av_malloc(size);
    memcpy(ret, start, size - 1);
    ret[size-1] = 0;

    unquote(ret);

    return ret;
}

/**
 * Parse "(linkname)"
 * @arg name a pointer (that need to be free'd after use) to the name between
 *           parenthesis
 */
static void parse_link_name(const char **buf, char **name)
{
    consume_char(buf);

    *name = consume_string(buf);

    if (!*name[0])
        goto fail;

    if (consume_char(buf) != ')')
        goto fail;

    return;
 fail:
    av_freep(name);
    av_log(&log_ctx, AV_LOG_ERROR, "Could not parse link name!\n");
}

/**
 * Parse "filter=params"
 * @arg name a pointer (that need to be free'd after use) to the name of the
 *           filter
 * @arg ars  a pointer (that need to be free'd after use) to the args of the
 *           filter
 */
static void parse_filter(const char **buf, char **name, char **opts)
{
    *name = consume_string(buf);

    if (**buf == '=') {
        consume_char(buf);
        *opts = consume_string(buf);
    } else {
        *opts = NULL;
    }

}

enum LinkType {
    LinkTypeIn,
    LinkTypeOut,
};

/**
 * A linked-list of the inputs/outputs of the filter chain.
 */
typedef struct AVFilterInOut {
    enum LinkType type;
    char *name;
    int instance;
    int pad_idx;

    struct AVFilterInOut *next;
} AVFilterInOut;

static void free_inout(AVFilterInOut *head)
{
    while (head) {
        AVFilterInOut *next;
        next = head->next;
        av_free(head);
        head = next;
    }
}

/**
 * Parse "(a1)(link2) ... (etc)"
 */
static int parse_inouts(const char **buf, AVFilterInOut **inout, int firstpad,
                        enum LinkType type, int instance)
{
    int pad = firstpad;
    while (**buf == '(') {
        AVFilterInOut *inoutn = av_malloc(sizeof(AVFilterInOut));
        parse_link_name(buf, &inoutn->name);
        inoutn->type = type;
        inoutn->instance = instance;
        inoutn->pad_idx = pad++;
        inoutn->next = *inout;
        *inout = inoutn;
    }
    return pad;
}

static AVFilterGraphDesc *parse_chain(const char *filters, int has_in)
{
    AVFilterGraphDesc        *ret;
    AVFilterGraphDescFilter **filterp, *filtern;
    AVFilterGraphDescLink   **linkp,   *linkn;
    AVFilterInOut           *inout=NULL;
    AVFilterInOut           *head;

    int index = 0;
    char chr = 0;
    int pad = 0;
    int has_out = 0;

    consume_whitespace(&filters);

    if(!(ret = av_mallocz(sizeof(AVFilterGraphDesc))))
        return NULL;

    filterp = &ret->filters;
    linkp   = &ret->links;

    do {
        if(chr == ',') {
            linkn = av_mallocz(sizeof(AVFilterGraphDescLink));
            linkn->src = index-1;
            linkn->srcpad = pad;
            linkn->dst = index;
            linkn->dstpad = 0;

            *linkp = linkn;
            linkp = &linkn->next;
        }
        pad = parse_inouts(&filters, &inout, chr == ',' || (!has_in),
                           LinkTypeIn, index);

        filtern = av_mallocz(sizeof(AVFilterGraphDescFilter));
        filtern->index = index;
        parse_filter(&filters, &filtern->filter, &filtern->args);
        *filterp = filtern;
        filterp = &filtern->next;

        pad = parse_inouts(&filters, &inout, 0,
                           LinkTypeOut, index);
        chr = consume_char(&filters);
        index++;
    } while (chr == ',' || chr == ';');

    head = inout;
    for (; inout != NULL; inout = inout->next) {
        if (inout->instance == -1)
            continue; // Already processed

        if (!strcmp(inout->name, "in")) {
            if (!has_in)
                goto fail;
            ret->inputs = av_mallocz(sizeof(AVFilterGraphDescExport));
            ret->inputs->filter = inout->instance;
            ret->inputs->pad = inout->pad_idx;
        } else if (!strcmp(inout->name, "out")) {
            has_out = 1;
            ret->outputs = av_mallocz(sizeof(AVFilterGraphDescExport));
            ret->outputs->filter = inout->instance;
            ret->outputs->pad = inout->pad_idx;
        } else {
            AVFilterInOut *p, *src, *dst;
            for (p = inout->next;
                 p && strcmp(p->name,inout->name); p = p->next);

            if (!p) {
                av_log(&log_ctx, AV_LOG_ERROR, "Unmatched link: %s.\n",
                       inout->name);
                goto fail;
            }

            if (p->type == LinkTypeIn && inout->type == LinkTypeOut) {
                src = inout;
                dst = p;
            } else if (p->type == LinkTypeOut && inout->type == LinkTypeIn) {
                src = p;
                dst = inout;
            } else {
                av_log(&log_ctx, AV_LOG_ERROR, "Two links named '%s' are either both input or both output\n",
                       inout->name);
                goto fail;
            }
            linkn = av_mallocz(sizeof(AVFilterGraphDescLink));

            linkn->src = src->instance;
            linkn->srcpad = src->pad_idx;
            linkn->dst = dst->instance;
            linkn->dstpad = dst->pad_idx;

            *linkp = linkn;
            linkp = &linkn->next;

            src->instance = -1;
            dst->instance = -1;
        }
    }

    free_inout(head);

    if (!has_in) {
        ret->inputs = av_mallocz(sizeof(AVFilterGraphDescExport));
        ret->inputs->filter = 0;
    }
    if (!has_out) {
        ret->outputs = av_mallocz(sizeof(AVFilterGraphDescExport));
        ret->outputs->filter = index-1;
    }

    return ret;

 fail:
    free_inout(head);

    avfilter_graph_free_desc(ret);
    return NULL;
}

/**
 * Parse a string describing a filter graph.
 */
AVFilterGraphDesc *avfilter_graph_parse_chain(const char *filters)
{
    AVFilterGraphDesc *ret;

    /* Try first to parse supposing there is no (in) element */
    if ((ret = parse_chain(filters, 0)))
        return ret;

    /* Parse supposing there is an (in) element */
    return parse_chain(filters, 1);
}

/**
 * Free a graph description.
 */
void avfilter_graph_free_desc(AVFilterGraphDesc *desc)
{
    void *next;

    while(desc->filters) {
        next = desc->filters->next;
        av_free(desc->filters->filter);
        av_free(desc->filters->args);
        av_free(desc->filters);
        desc->filters = next;
    }

    while(desc->links) {
        next = desc->links->next;
        av_free(desc->links);
        desc->links = next;
    }

    while(desc->inputs) {
        next = desc->inputs->next;
        av_free(desc->inputs);
        desc->inputs = next;
    }

    while(desc->outputs) {
        next = desc->outputs->next;
        av_free(desc->outputs);
        desc->outputs = next;
    }
}

