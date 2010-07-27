/*
 * MMS protocol over HTTP
 * Copyright (c) 2010 Zhentan Feng <spyfeng at gmail dot com>
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
#ifndef AVFORMAT_MMS_H
#define AVFORMAT_MMS_H

#include "avformat.h"

typedef struct {
    int id;
}MMSStream;

typedef struct {
    int outgoing_packet_seq;             ///< Outgoing packet sequence number.
    char path[256];                      ///< Path of the resource being asked for.
    char host[128];                      ///< Host of the resources.

    URLContext *mms_hd;                  ///< TCP connection handle
    MMSStream streams[MAX_STREAMS];

    /** Buffer for outgoing packets. */
    /*@{*/
    uint8_t *write_out_ptr;              ///< Pointer for writting the buffer.
    uint8_t out_buffer[512];             ///< Buffer for outgoing packet.
    /*@}*/

    /** Buffer for incoming packets. */
    /*@{*/
    uint8_t in_buffer[8192];             ///< Buffer for incoming packets.
    uint8_t *read_in_ptr;                ///< Pointer for reading from incoming buffer.
    int remaining_in_len;                ///< Reading length from incoming buffer.
    /*@}*/

    int incoming_packet_seq;             ///< Incoming packet sequence number.
    int incoming_flags;                  ///< Incoming packet flags.

    int packet_id;                       ///< Identifier for packets in the current stream.
    unsigned int header_packet_id;       ///< default is 2.

    /** Internal handling of the ASF header */
    /*@{*/
    uint8_t *asf_header;                 ///< Stored ASF header.
    int asf_header_size;                 ///< Size of stored ASF header.
    int header_parsed;                   ///< The header has been received and parsed.
    int asf_packet_len;
    int asf_header_read_size;
    /*@}*/

    int stream_num;                      ///< stream numbers.
    int is_playing;
} MMSContext;
#endif
