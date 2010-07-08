/*
 * General components used by playlist formats
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

/** @file
 *  @author Geza Kovacs ( gkovacs mit edu )
 *
 *  @brief General components used by playlist formats
 *
 *  @details These functions are used to initialize and manipulate playlists
 *  (AVPlaylistContext) and AVFormatContext for each playlist item.
 *  This abstraction is used to implement file concatenation and
 *  support for playlist formats.
 */

#include "avplaylist.h"
#include "riff.h"
#include "libavutil/avstring.h"
#include "internal.h"
#include "playlist.h"

AVPlaylistContext *av_playlist_alloc(void)
{
    return av_mallocz(sizeof(AVPlaylistContext));
}

int av_playlist_insert_item(AVPlaylistContext *ctx, const char *itempath, int pos)
{
    int i, itempath_len;
    int64_t *durations_tmp, durations_offset;
    unsigned int *nb_streams_list_tmp, nb_streams_offset;
    AVFormatContext **formatcontext_list_tmp;
    char **flist_tmp;
    AVFormatContext *ic;
    ++ctx->pelist_size;
    flist_tmp = av_realloc(ctx->flist, sizeof(*(ctx->flist)) * ctx->pelist_size);
    durations_tmp = av_realloc(ctx->durations,
                               sizeof(*(ctx->durations)) * ctx->pelist_size);
    nb_streams_list_tmp = av_realloc(ctx->nb_streams_list,
                                     sizeof(*(ctx->nb_streams_list)) * ctx->pelist_size);
    formatcontext_list_tmp = av_realloc(ctx->formatcontext_list,
                                        sizeof(*(ctx->formatcontext_list)) * ctx->pelist_size);
    if (!flist_tmp ||
        !durations_tmp ||
        !nb_streams_list_tmp ||
        !formatcontext_list_tmp) {
        av_log(ctx, AV_LOG_ERROR,
               "av_realloc error in av_playlist_insert_item\n");
        return AVERROR_NOMEM;
    } else {
        ctx->flist              = flist_tmp;
        ctx->durations          = durations_tmp;
        ctx->nb_streams_list    = nb_streams_list_tmp;
        ctx->formatcontext_list = formatcontext_list_tmp;
    }
    ic = ff_playlist_alloc_formatcontext(itempath);
    if (!ic) {
        av_log(ctx, AV_LOG_ERROR,
               "failed to allocate and open %s in av_playlist_insert_item\n",
               itempath);
        return AVERROR_NOMEM;
    }
    if (pos > 0) {
        durations_offset  = ic->duration + ctx->durations[pos - 1];
        nb_streams_offset = ic->nb_streams + ctx->nb_streams_list[pos - 1];
    } else {
        durations_offset  = ic->duration;
        nb_streams_offset = ic->nb_streams;
    }
    av_close_input_file(ic);
    for (i = ctx->pelist_size; i > pos; --i) {
        ctx->flist[i] = ctx->flist[i - 1];
        ctx->durations[i] = ctx->durations[i - 1] + durations_offset;
        ctx->nb_streams_list[i] = ctx->nb_streams_list[i - 1] + nb_streams_offset;
        ctx->formatcontext_list[i] = ctx->formatcontext_list[i - 1];
    }
    ctx->formatcontext_list[pos] = NULL;
    ctx->durations[pos]          = durations_offset;
    ctx->nb_streams_list[pos]    = nb_streams_offset;
    itempath_len                 = strlen(itempath);
    ctx->flist[pos]              = av_malloc(itempath_len + 1);
    if (!ctx->flist[pos]) {
        av_log(ctx, AV_LOG_ERROR,
               "av_malloc error in av_playlist_insert_item\n");
        return AVERROR_NOMEM;
    }
    av_strlcpy(ctx->flist[pos], itempath, itempath_len + 1);
    // duration is updated in case it's checked by a parent demuxer (chained concat demuxers)
    if (ctx->master_formatcontext)
        ctx->master_formatcontext->duration = ctx->durations[ctx->pelist_size - 1];
    return 0;
}

int av_playlist_insert_playlist(AVPlaylistContext *ctx, AVPlaylistContext *insert_ctx, int pos)
{
    int i, err;
    for (i = 0; i < insert_ctx->pelist_size; ++i) {
        err = av_playlist_insert_item(ctx, insert_ctx->flist[i], pos + i);
        if (err) {
            av_log(ctx, AV_LOG_ERROR,
                   "failed to insert item %d to new position %d in av_playlist_insert_playlist\n",
                   i, pos + i);
            return err;
        }
    }
    return 0;
}

int av_playlist_remove_item(AVPlaylistContext *ctx, int pos)
{
    int i;
    int64_t *durations_tmp, durations_offset;
    unsigned int *nb_streams_list_tmp, nb_streams_offset;
    AVFormatContext **formatcontext_list_tmp;
    char **flist_tmp;
    if (pos >= ctx->pelist_size || !ctx->flist || !ctx->durations || !ctx->nb_streams_list) {
        av_log(ctx, AV_LOG_ERROR,
               "failed to remove item %d which is not present in playlist\n",
               pos);
        return AVERROR_INVALIDDATA;
    }
    --ctx->pelist_size;
    av_free(ctx->flist[pos]);
    if (pos > 0) {
        durations_offset  = ctx->durations[pos] - ctx->durations[pos - 1];
        nb_streams_offset = ctx->nb_streams_list[pos] - ctx->nb_streams_list[pos - 1];
    } else {
        durations_offset  = ctx->durations[pos];
        nb_streams_offset = ctx->nb_streams_list[pos];
    }
    if (ctx->formatcontext_list && ctx->formatcontext_list[pos]) {
        av_close_input_file(ctx->formatcontext_list[pos]);
        av_close_input_stream(ctx->formatcontext_list[pos]);
        av_free(ctx->formatcontext_list[pos]);
        ctx->formatcontext_list[pos] = NULL;
    }
    for (i = pos; i < ctx->pelist_size; ++i) {
        ctx->flist[i] = ctx->flist[i + 1];
        ctx->durations[i] = ctx->durations[i + 1] - durations_offset;
        ctx->nb_streams_list[i] = ctx->nb_streams_list[i + 1] - nb_streams_offset;
        ctx->formatcontext_list[i] = ctx->formatcontext_list[i + 1];
    }
    flist_tmp = av_realloc(ctx->flist, sizeof(*(ctx->flist)) * ctx->pelist_size);
    durations_tmp = av_realloc(ctx->durations,
                               sizeof(*(ctx->durations)) * ctx->pelist_size);
    nb_streams_list_tmp = av_realloc(ctx->nb_streams_list,
                                     sizeof(*(ctx->nb_streams_list)) * ctx->pelist_size);
    formatcontext_list_tmp = av_realloc(ctx->formatcontext_list,
                                        sizeof(*(ctx->formatcontext_list)) * ctx->pelist_size);
    if (!flist_tmp ||
        !durations_tmp ||
        !nb_streams_list_tmp ||
        !formatcontext_list_tmp) {
        av_log(ctx, AV_LOG_ERROR,
               "av_realloc error in av_playlist_remove_item\n");
        return AVERROR_NOMEM;
    } else {
        ctx->flist              = flist_tmp;
        ctx->durations          = durations_tmp;
        ctx->nb_streams_list    = nb_streams_list_tmp;
        ctx->formatcontext_list = formatcontext_list_tmp;
    }
    // duration is updated in case it's checked by a parent demuxer (chained concat demuxers)
    if (ctx->master_formatcontext)
        ctx->master_formatcontext->duration = ctx->durations[ctx->pelist_size - 1];
    return 0;
}

int av_playlist_close(AVPlaylistContext *ctx)
{
    int err;
    while (ctx->pelist_size > 0) {
        err = av_playlist_remove_item(ctx, ctx->pelist_size-1);
        if (err) {
            av_log(ctx, AV_LOG_ERROR,
                   "failed to remove item %d from playlist",
                   ctx->pelist_size-1);
            return err;
        }
    }
    av_free(ctx->flist);
    av_free(ctx->durations);
    av_free(ctx->nb_streams_list);
    av_free(ctx->formatcontext_list);
    av_free(ctx);
    return 0;
}
