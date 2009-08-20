/*
 * XSPF playlist demuxer
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
 *  @brief XSPF playlist demuxer
 */

#include "concatgen.h"

/* The ffmpeg codecs we support, and the IDs they have in the file */
static const AVCodecTag codec_xspf_tags[] = {
    { 0, 0 },
};

static int xspf_probe(AVProbeData *p)
{
    int i, len;
    char fxml;
    char ftag;
    unsigned char c;
    char t[] = "<?xml";
    char s[5] = {0};
    char v[] = "<playlist";
    char u[9] = {0};
    fxml = ftag = 0;
    len = p->buf_size;
    for (i = 0; i < len; i++) {
        if (ftag)
            break;
        c = p->buf[i];
        if (!fxml) {
            s[0] = s[1];
            s[1] = s[2];
            s[2] = s[3];
            s[3] = s[4];
            s[4] = c;
            if (!memcmp(s, t, 5))
                fxml = 1;
        }
        if (!ftag) {
            u[0] = u[1];
            u[1] = u[2];
            u[2] = u[3];
            u[3] = u[4];
            u[4] = u[5];
            u[5] = u[6];
            u[6] = u[7];
            u[7] = u[8];
            u[8] = c;
            if (!memcmp(u, v, 9))
                ftag = 1;
        }
    }
    if (fxml && ftag)
        return AVPROBE_SCORE_MAX;
    else if (fxml || ftag)
        return AVPROBE_SCORE_MAX / 2;
    else
        return 0;
}

static int xspf_list_files(ByteIOContext *b, PlaylistContext *ctx, const char *filename)
{
    int i, j, k, c;
    unsigned int buflen;
    char state;
    char **flist;
    char buf[1024];
    char s[10] = {0};
    char t[] = "<location>";
    flist = NULL;
    state = buflen = i = j = 0;
    while ((c = url_fgetc(b))) {
        if (c == EOF)
            break;
        if (state == 0) {
            s[0] = s[1];
            s[1] = s[2];
            s[2] = s[3];
            s[3] = s[4];
            s[4] = s[5];
            s[5] = s[6];
            s[6] = s[7];
            s[7] = s[8];
            s[8] = s[9];
            s[9] = c;
            if (!memcmp(s, t, 10))
                state = 1;
        } else {
            if (c == '<') {
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
    ff_playlist_relative_paths(flist, j, dirname(filename));
    for (k = 0; k < j; ++k)
        ff_playlist_add_path(ctx, flist[k]);
    av_free(flist);
    return 0;
}

static int xspf_read_header(AVFormatContext *s,
                            AVFormatParameters *ap)
{
    PlaylistContext *ctx = av_mallocz(sizeof(*ctx));
    if (xspf_list_files(s->pb, ctx, s->filename)) {
        fprintf(stderr, "no playlist items found in %s\n", s->filename);
        return AVERROR_EOF;
    }
    s->priv_data = ctx;
    ff_playlist_populate_context(ctx, ctx->pe_curidx);
    ff_playlist_set_streams(s);
    return 0;
}

AVInputFormat xspf_demuxer = {
    "xspf",
    NULL_IF_CONFIG_SMALL("CONCAT XSPF format"),
    sizeof(PlaylistContext),
    xspf_probe,
    xspf_read_header,
    ff_concatgen_read_packet,
    ff_concatgen_read_close,
    ff_concatgen_read_seek,
    ff_concatgen_read_timestamp,
    0, //flags
    NULL, //extensions
    0, //value
    ff_concatgen_read_play,
    ff_concatgen_read_pause,
    (const AVCodecTag* const []){codec_xspf_tags, 0},
};
