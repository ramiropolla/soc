/*
 * PLS playlist demuxer
 * Copyright (c) 2009 Geza Kovacs
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

/** @file libavformat/pls.c
 *  @author Geza Kovacs ( gkovacs mit edu )
 *
 *  @brief PLS playlist demuxer
 */

#include "concatgen.h"

/* The ffmpeg codecs we support, and the IDs they have in the file */
static const AVCodecTag codec_pls_tags[] = {
    { 0, 0 },
};

static int pls_probe(AVProbeData *p)
{
    if (!strncmp(p->buf, "[playli", 7))
        return AVPROBE_SCORE_MAX;
    else
        return 0;
}

static int pls_list_files(ByteIOContext *b, PlaylistContext *ctx, const char *filename)
{
    int i, j, k, c;
    unsigned int buflen;
    char state;
    char **flist;
    char buf[1024];
    char buf_tag[5] = {0};
    const char match_tag[] = "\nFile";
    flist = NULL;
    state = buflen = i = j = 0;
    while ((c = url_fgetc(b))) {
        if (c == EOF)
            break;
        if (state == 0) {
            memmove(buf_tag, buf_tag+1, 4);
            buf_tag[4] = c;
            if (!memcmp(buf_tag, match_tag, 5))
                state = 1;
        } else if (state == 1) {
            if (c == '=')
                state = 2;
            else if (c == '#')
                state = 0;
        } else {
            if (c == '\n' || c == '#') {
                termfn:
                buf[i++] = 0;
                flist = av_fast_realloc(flist, &buflen, sizeof(*flist) * (j+2));
                flist[j] = av_malloc(i);
                av_strlcpy(flist[j++], buf, i);
                i = 0;
                state = 0;
                buf_tag[sizeof(buf_tag)-1] = c;
                continue;
            } else {
                buf[i++] = c;
                if (i >= sizeof(buf)-1)
                    goto termfn;
            }
        }
    }
    if (!flist) // no files have been found
        return AVERROR_EOF;
    flist[j] = 0;
    ff_playlist_relative_paths(flist, j, dirname(filename));
    for (k = 0; k < j; ++k)
        ff_playlist_add_path(ctx, flist[k]);
    av_free(flist);
    return 0;
}

static int pls_read_header(AVFormatContext *s,
                           AVFormatParameters *ap)
{
    PlaylistContext *ctx = av_mallocz(sizeof(*ctx));
    if (pls_list_files(s->pb, ctx, s->filename)) {
        fprintf(stderr, "no playlist items found in %s\n", s->filename);
        return AVERROR_EOF;
    }
    s->priv_data = ctx;
    ff_playlist_populate_context(ctx, ctx->pe_curidx);
    ff_playlist_set_streams(s);
    return 0;
}

AVInputFormat pls_demuxer = {
    "pls",
    NULL_IF_CONFIG_SMALL("CONCAT PLS format"),
    sizeof(PlaylistContext),
    pls_probe,
    pls_read_header,
    ff_concatgen_read_packet,
    ff_concatgen_read_close,
    ff_concatgen_read_seek,
    ff_concatgen_read_timestamp,
    0, //flags
    NULL, //extensions
    0, //value
    ff_concatgen_read_play,
    ff_concatgen_read_pause,
    (const AVCodecTag* const []){codec_pls_tags, 0},
};
