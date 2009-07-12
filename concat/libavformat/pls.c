/*
 * PLS muxer and demuxer
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

#include "playlist.h"

static int ff_concatgen_read_packet(AVFormatContext *s, AVPacket *pkt);

static int ff_concatgen_read_seek(AVFormatContext *s, int stream_index, int64_t pts, int flags);

static int ff_concatgen_read_timestamp(AVFormatContext *s, int stream_index, int64_t *pos, int64_t pos_limit);

static int ff_concatgen_read_close(AVFormatContext *s);

static int ff_concatgen_read_play(AVFormatContext *s);

static int ff_concatgen_read_pause(AVFormatContext *s);

/* The ffmpeg codecs we support, and the IDs they have in the file */
static const AVCodecTag codec_pls_tags[] = {
    { 0, 0 },
};

static int pls_probe(AVProbeData *p)
{
    if (p->buf != 0) {
        if (!strncmp(p->buf, "[playli", 7))
            return AVPROBE_SCORE_MAX;
        else
            return 0;
    }
    if (match_ext(p->filename, "pls"))
        return AVPROBE_SCORE_MAX/2;
    else
        return 0;
}

static int pls_list_files(ByteIOContext *b, PlaylistContext *ctx)
{
    int i, j, c;
    unsigned int buflen;
    char state;
    char **flist;
    char buf[1024];
    char s[5];
    char t[] = "\nFile";
    state = flist = buflen = i = j = 0;
    while ((c = url_fgetc(b))) {
        if (c == EOF)
            break;
        if (state == 0) {
            s[0] = s[1];
            s[1] = s[2];
            s[2] = s[3];
            s[3] = s[4];
            s[4] = c;
            if (s[0] == t[0] && s[1] == t[1] && s[2] == t[2] && s[3] == t[3] && s[4] == t[4])
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
                s[sizeof(s)-1] = c;
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
    ctx->pelist_size = j;
    ff_playlist_relative_paths(flist, ctx->workingdir);
    ctx->pelist = av_malloc(ctx->pelist_size * sizeof(*(ctx->pelist)));
    memset(ctx->pelist, 0, ctx->pelist_size * sizeof(*(ctx->pelist)));
    for (i = 0; i < ctx->pelist_size; ++i) {
        ctx->pelist[i] = av_malloc(sizeof(*(ctx->pelist[i])));
        ctx->pelist[i]->filename = flist[i];
    }
    av_free(flist);
    return 0;
}

static int pls_read_header(AVFormatContext *s,
                           AVFormatParameters *ap)
{
    int i;
    PlaylistContext *ctx = ff_playlist_alloc_context(s->filename);
    if (pls_list_files(s->pb, ctx)) {
        fprintf(stderr, "no playlist items found in %s\n", s->filename);
        return AVERROR_EOF;
    }
    s->priv_data = ctx;
    for (i = 0; i < ctx->pe_curidxs_size; ++i) {
        ff_playlist_populate_context(ctx, s, i);
    }
    return 0;
}

#if CONFIG_PLS_DEMUXER
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
    NULL, //flags
    NULL, //extensions
    NULL, //value
    ff_concatgen_read_play,
    ff_concatgen_read_pause,
    (const AVCodecTag* const []){codec_pls_tags, 0},
    ff_concatgen_read_seek, //m3u_read_seek2
    NULL, //metadata_conv
    NULL, //next
};
#endif
