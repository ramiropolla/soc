/*
 * Internal functions used to manipulate playlists
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
 *  @brief Internal functions used to manipulate playlists
 *
 *  @details These functions are used internally for manipulating playlists.
 *  The public playlist API can be found in avplaylist.h
 */

#ifndef AVFORMAT_PLAYLIST_H
#define AVFORMAT_PLAYLIST_H

#include "avplaylist.h"

/** @brief Allocates AVFormatContext, then opens file, and probes and opens streams.
 *  @param filename Null-terminated string of path to file to open.
 *  @return Returns an allocated AVFormatContext upon success, or NULL upon failure.
 */
AVFormatContext *ff_playlist_alloc_formatcontext(char *filename);

/** @brief Allocates a new AVFormatContext for a concat-type demuxer.
 *  @return Returns NULL if failed, or AVFormatContext if succeeded.
 */
AVFormatContext *ff_playlist_alloc_concat_formatcontext(void);

/** @brief Sets the master concat-type demuxer's streams to those of its currently opened playlist element.
 *  Does nothing if using a standalone playlist (master_formatcontext is NULL).
 *  @param ctx AVPlaylistContext within which the list of playlist elements and durations are stored.
 *  @return Returns 0 upon success, or negative upon failure.
 */
int ff_playlist_set_streams(AVPlaylistContext *ctx);

/** @brief Converts a list of mixed relative or absolute paths into all absolute paths.
 *  @param flist List of null-terminated strings of relative or absolute paths.
 *  @param len Number of paths in flist.
 *  @param workingdir Path that strings in flist are relative to.
 */
void ff_playlist_relative_paths(char **flist,
                                int len,
                                const char *workingdir);

/** @brief Splits a character-delimited string into a list of strings.
 *  @param s The input character-delimited string ("one,two,three").
 *  @param sep The delimiter character (',').
 *  @param flist_ptr Pointer to string list which will be allocated by function.
 *  @param len_ptr Number of segments the string was split into.
 *  @return Returns 0 upon success, or negative upon failure.
 */
int ff_playlist_split_encodedstring(const char *s,
                                    const char sep,
                                    char ***flist_ptr,
                                    int *len_ptr);

/** @brief Calculates the index of the playlist item which would contain the timestamp specified in AV_TIME_BASE units.
 *  @param ctx AVPlaylistContext within which the list of playlist elements and durations are stored.
 *  @param pts Timestamp in AV_TIME_BASE.
 *  @param localpts Time in the local demuxer's timeframe in AV_TIME_BASE units; if null, not calculated.
 *  @return Returns the index of the stream which covers the specified time range.
 */
int ff_playlist_stream_index_from_time(AVPlaylistContext *ctx,
                                       int64_t pts,
                                       int64_t *localpts);

/** @brief Calculates the local stream index which corresponds to a global stream index.
 *  @param ctx AVPlaylistContext within which the list of playlist elements and durations are stored.
 *  @param stream_index Global stream index, the index of the stream within the playlist demuxer.
 *  @return Returns the local stream index, the index of the stream within the child demuxer.
 */
int ff_playlist_localstidx_from_streamidx(AVPlaylistContext *ctx, int stream_index);

#endif /* AVFORMAT_PLAYLIST_H */

