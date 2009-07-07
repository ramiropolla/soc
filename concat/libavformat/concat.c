/*
 * Standard playlist/concatenation demuxer
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

#include "concat.h"

/* The ffmpeg codecs we support, and the IDs they have in the file */
static const AVCodecTag codec_concat_tags[] = {
    { 0, 0 },
};

static int concat_probe(AVProbeData *p)
{
    // concat demuxer should only be manually constructed in ffmpeg
    return 0;
}

static int concat_read_header(AVFormatContext *s,
                       AVFormatParameters *ap)
{
    // PlaylistD should be constructed externally
    return 0;
}

void concat_init_demuxer(AVInputFormat *cdm)
{
    cdm->name = "concat";
    cdm->long_name = NULL_IF_CONFIG_SMALL("CONCAT format");
    cdm->priv_data_size = sizeof(PlaylistContext);
    cdm->read_probe = concat_probe;
    cdm->read_header = concat_read_header;
    cdm->read_packet = ff_concatgen_read_packet;
    cdm->read_close = ff_concatgen_read_close;
    cdm->read_seek = ff_concatgen_read_seek;
    cdm->read_timestamp = ff_concatgen_read_timestamp;
    cdm->flags = NULL;
    cdm->extensions = NULL;
    cdm->value = NULL;
    cdm->read_play = ff_concatgen_read_play;
    cdm->read_pause = ff_concatgen_read_pause;
    cdm->codec_tag = codec_concat_tags;
    cdm->read_seek2 = ff_concatgen_read_seek;
    cdm->metadata_conv = NULL;
    cdm->next = NULL;
}

#if CONFIG_CONCAT_DEMUXER
AVInputFormat concat_demuxer = {
    "concat",
    NULL_IF_CONFIG_SMALL("CONCAT format"),
    0,
    concat_probe,
    concat_read_header,
    ff_concatgen_read_packet,
    ff_concatgen_read_close,
    ff_concatgen_read_seek,
    ff_concatgen_read_timestamp,
    NULL, //flags
    NULL, //extensions
    NULL, //value
    ff_concatgen_read_play,
    ff_concatgen_read_pause,
    (const AVCodecTag* const []){codec_concat_tags, 0},
    NULL, //m3u_read_seek2
    NULL, //metadata_conv
    NULL, //next
};
#endif
