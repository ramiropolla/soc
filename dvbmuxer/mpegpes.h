/*
 * Copyright (c) 2000-2002 Fabrice Bellard
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

/**
 * @file mpegpes.h
 * MPEG PES packetizer API header
 */

#ifndef FFMPEG_MPEGPES_H
#define FFMPEG_MPEGPES_H

#include "avformat.h"
#include "libavutil/fifo.h"

#define PES_FMT_MPEG2 0x01
#define PES_FMT_VCD   0x02
// formats below are all mpeg-2
#define PES_FMT_SVCD  0x05
#define PES_FMT_DVD   0x09
#define PES_FMT_TS    0x11

/**
 * PES packet description
 */
typedef struct PacketDesc {
    int64_t pts;
    int64_t dts;
    int size;
    int unwritten_size;
    int flags;
    struct PacketDesc *next;
} PacketDesc;

/**
 * PES stream structure
 */
typedef struct {
    AVFifoBuffer fifo;
    uint8_t id;
    int max_buffer_size; /**< in bytes */
    int buffer_index;
    PacketDesc *predecode_packet;
    PacketDesc *premux_packet;
    PacketDesc **next_packet;
    int packet_number;
    uint8_t lpcm_header[3];
    int lpcm_align;
    int bytes_to_iframe;
    int align_iframe;
    int64_t vobu_start_pts;
    int format; ///< mux format (ts, svcd ...)
} StreamInfo;

#define PRIVATE_STREAM_1   0x1bd
#define PADDING_STREAM     0x1be

/**
 * Initialization of PES muxer.
 * @param[in] ctx AVFormatContext
 * @return        Negative value on error, zero on success.
 */
int ff_pes_muxer_init(AVFormatContext *ctx);

/**
 * Finalization of PES muxer.
 * @param [in] ctx AVFormatContext
 */
void ff_pes_muxer_end(AVFormatContext *ctx);

/**
 * Write packet into PES FIFO.
 * @param [in] ctx  AVFormatContext
 * @param [in] pkt  Packet to write.
 */
void ff_pes_write_packet(AVFormatContext *ctx, AVPacket *pkt);

/**
 * Find the best stream to mux into the PES stream.
 * @param[in]  ctx         AVFormatContext
 * @param[in]  packet_size PES stream packet size
 * @param[in]  flush       Flush after every single subtitle packet.
 * @param[in]  scr         Current clock reference
 * @param[out] scr         Updated clock reference, bumped if needed
 * @param[out] best_i      Index of the stream to be muxed
 * @return                 Negative on error, zero if not found, 1 if found.
 */
int ff_pes_find_beststream(AVFormatContext *ctx, int packet_size, int flush, int64_t *scr, int *best_i);

/**
 * Get total number of frames of PES stream to be muxed considering len
 * @param[in] ctx    AVFormatContext
 * @param[in] stream PES stream
 * @param[in] len    Bytes available in next PES packet
 * @return           Number of frames muxed.
 */
int ff_pes_get_nb_frames(AVFormatContext *ctx, StreamInfo *stream, int len);

/**
 * Caculate next PES packet informations
 * @param[in]  stream            PES stream
 * @param[out] packet_size       PES packet size
 * @param[out] header_len        PES header length
 * @param[in]  pts               Current pts
 * @param[out] pts               AV_NOPTS_VALUE if reset
 * @param[in]  dts               Current dts
 * @param[out] dts               AV_NOPTS_VALUE if reset
 * @param[out] payload_size      PES payload size
 * @param[out] startcode         PES startcode
 * @param[out] stuffing_size     PES stuffing size
 * @param[out] trailer_size      Unwritten trailer size
 * @param[out] pad_packet_bytes  Packet padding size
 */
void ff_pes_cal_header(StreamInfo *stream,
          int *packet_size,  int *header_len, int64_t *pts, int64_t *dts,
          int *payload_size, int *startcode, int *stuffing_size,
          int *trailer_size, int *pad_packet_bytes);

/**
 * Write PES data from PES Stream into supplied buffer.
 * @param [in] ctx            AVFormatContext
 * @param [in] stream_index   Stream index to write from
 * @param [in] buf            Buffer to write to
 * @param [in] pts            Packet pts
 * @param [in] dts            Packet dts
 * @param [in] start_code     PES packet startcode
 * @param [in] header_len     PES header size
 * @param [in] packet_size    Total packet size
 * @param [in] payload_size   Packet payload size
 * @param [in] stuffing_size  Packet stuffing size
 * @return                    Bytes written to buffer
 */
int ff_pes_write_buf(AVFormatContext *ctx, int stream_index, uint8_t *buf,
          int64_t pts, int64_t dts, int startcode,
          int header_len, int packet_size, int payload_size, int stuffing_size);

/**
 * Remove decoded packets of each stream.
 * @param[in] ctx  the AVFormatContext
 * @param[in] scr  System Clock Reference of PES stream
 * @return         Negative value on error, 0 on success.
 */
int ff_pes_remove_decoded_packets(AVFormatContext *ctx, int64_t scr);

#endif /* FFMPEG_MPEGPES_H */
