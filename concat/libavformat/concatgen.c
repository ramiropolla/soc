/*
 * Generic functions used by playlist/concatenation demuxers
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

/** @file libavformat/concatgen.c
 *  @author Geza Kovacs ( gkovacs mit edu )
 *
 *  @brief Generic functions used by playlist/concatenation demuxers
 *
 *  @details These functions are used to read packets and seek streams
 *  for concat-type demuxers, abstracting away the playlist element switching
 *  process.
 */

#include "concatgen.h"

int ff_concatgen_read_packet(AVFormatContext *s,
                             AVPacket *pkt)
{
    int ret, i, stream_index;
    PlaylistContext *ctx;
    AVFormatContext *ic;
    char have_switched_streams = 0;
    ctx = s->priv_data;
    stream_index = 0;
//    pkt = NULL;
    for (;;) {
        ic = ctx->icl[ctx->pe_curidx];
//        ff_playlist_set_streams(s);
        ret = ic->iformat->read_packet(ic, pkt);
//        ff_playlist_set_streams(s);
        s->cur_st = ic->cur_st;
        if (ret >= 0) {
            if (pkt) {
                pkt->stream = ic->streams[pkt->stream_index];
            stream_index = pkt->stream_index;
            pkt->index_offset = ff_playlist_streams_offset_from_playidx(ctx, ctx->pe_curidx);
            pkt->stream_index += pkt->index_offset;
//                ff_playlist_set_streams(s);
                if (!ic->streams[stream_index]->codec->has_b_frames) {
                    pkt->dts += av_rescale_q(ff_playlist_time_offset(ctx->durations, ctx->pe_curidx),
                                             AV_TIME_BASE_Q,
                                             ic->streams[stream_index]->time_base);
                    pkt->pts = pkt->dts + 1;
                }
            }
            break;
        } else {
            if (!have_switched_streams &&
                ctx->pe_curidx < ctx->pelist_size - 1 &&
                ic->cur_st) {
            // TODO switch from AVERROR_EOF to AVERROR_EOS
            // -32 AVERROR_EOF for avi, -51 for ogg
                av_log(ic, AV_LOG_DEBUG, "Switching stream %d to %d\n", stream_index, ctx->pe_curidx+1);
                ctx->durations[ctx->pe_curidx] = ic->duration;
                ctx->pe_curidx = ff_playlist_stream_index_from_time(ctx,
                                                                    ff_playlist_time_offset(ctx->durations, ctx->pe_curidx),
                                                                    NULL);
                ff_playlist_populate_context(ctx, ctx->pe_curidx);
                ff_playlist_set_streams(s);
                // have_switched_streams is set to avoid infinite loop
                have_switched_streams = 1;
                // duration is updated in case it's checked by a parent demuxer (chained concat demuxers)
                s->duration = 0;
                for (i = 0; i < ctx->pe_curidx; ++i)
                    s->duration += ctx->durations[i];
                continue;
            } else {
                av_log(ic, AV_LOG_ERROR, "Packet read error %d\n", ret);
                break;
            }
        }
    }
    return ret;
}

int ff_concatgen_read_seek(AVFormatContext *s,
                           int stream_index,
                           int64_t pts,
                           int flags)
{
    int i;
    int64_t localpts_avtimebase, localpts, pts_avtimebase;
    PlaylistContext *ctx;
    AVFormatContext *ic;
    ctx = s->priv_data;
    ic = ctx->icl[ctx->pe_curidx];
    ctx->durations[ctx->pe_curidx] = ic->duration;
    pts_avtimebase = av_rescale_q(pts,
                                  ic->streams[stream_index]->time_base,
                                  AV_TIME_BASE_Q);
    ctx->pe_curidx = ff_playlist_stream_index_from_time(ctx,
                                                        pts_avtimebase,
                                                        &localpts_avtimebase);
    ff_playlist_populate_context(ctx, ctx->pe_curidx);
    ff_playlist_set_streams(s);
    ic = ctx->icl[ctx->pe_curidx];
    localpts = av_rescale_q(localpts_avtimebase,
                            AV_TIME_BASE_Q,
                            ic->streams[stream_index]->time_base);
    s->duration = 0;
    for (i = 0; i < ctx->pe_curidx; ++i)
        s->duration += ctx->durations[i];
    return ic->iformat->read_seek(ic, stream_index, localpts, flags);
}

int64_t ff_concatgen_read_timestamp(AVFormatContext *s,
                                    int stream_index,
                                    int64_t *pos,
                                    int64_t pos_limit)
{
    PlaylistContext *ctx;
    AVFormatContext *ic;
    ctx = s->priv_data;
    ic = ctx->icl[ctx->pe_curidx];
    if (ic->iformat->read_timestamp)
        return ic->iformat->read_timestamp(ic, stream_index, pos, pos_limit);
    return 0;
}

int ff_concatgen_read_close(AVFormatContext *s)
{
    PlaylistContext *ctx;
    AVFormatContext *ic;
    ctx = s->priv_data;
    ic = ctx->icl[ctx->pe_curidx];
    if (ic->iformat->read_close)
        return ic->iformat->read_close(ic);
    return 0;
}

int ff_concatgen_read_play(AVFormatContext *s)
{
    PlaylistContext *ctx;
    AVFormatContext *ic;
    ctx = s->priv_data;
    ic = ctx->icl[ctx->pe_curidx];
    return av_read_play(ic);
}

int ff_concatgen_read_pause(AVFormatContext *s)
{
    PlaylistContext *ctx;
    AVFormatContext *ic;
    ctx = s->priv_data;
    ic = ctx->icl[ctx->pe_curidx];
    return av_read_pause(ic);
}
