/*
 * Copyright (c) 2000-2002 Fabrice Bellard
 * Copyright (c) 2007 Xiaohui Sun <sunxiaohui@dsp.ac.cn>
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
 * @file mpeg_pes.h
 * MPEG PES packetizer API header
 */

#ifndef AVFORMAT_MPEG_PES_H
#define AVFORMAT_MPEG_PES_H

#include "avformat.h"
#include "fifo.h"

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
 * muxer type for PES
 */
typedef enum {
    PESMUXER_PS,
    PESMUXER_TS,
    PESMUXER_PES
} PESMuxerType;

/**
 * PES context
 */
typedef struct {
    PESMuxerType muxer_type;
    int packet_number;
} PESContext;

/**
 * PES stream structure
 */
typedef struct {
    AVFifoBuffer fifo;
    int max_buffer_size; /**< in bytes */
    int buffer_index;
    PacketDesc *predecode_packet;
    PacketDesc *premux_packet;
    PacketDesc **next_packet;
} PESStream;


#define AUDIO_ID 0xc0
#define VIDEO_ID 0xe0
#define AC3_ID   0x80
#define DTS_ID   0x8a
#define LPCM_ID  0xa0
#define SUB_ID   0x20

#define PROGRAM_STREAM_MAP 0x1bc
#define PRIVATE_STREAM_1   0x1bd
#define PADDING_STREAM     0x1be
#define PRIVATE_STREAM_2   0x1bf



/**
 * Initialization of PES muxer.
 * @param[in] ctx the AVFormatContext which contains streams
 * @return  On error a negative value is returned, on success zero.
 */
int ff_pes_muxer_init(AVFormatContext *ctx);

/**
 * Finalization of PES muxer.
 * @param [in] ctx the AVFormatContext which contains streams.
 * @return  NULL
 */
void ff_pes_muxer_end(AVFormatContext *ctx);

/**
 * Write packet into PES FIFO.
 * @param [in] ctx  the AVFormatContext which contains streams.
 * @param [in] pkt  the packet to write.
 * @return  NULL
 */
void ff_pes_write_packet(AVFormatContext *ctx, AVPacket *pkt);

/**
 * Find the stream to mux into the PES stream.
 * @param[in] ctx          the AVFormatContext
 * @param[in] packet_size  PES stream packet size
 * @param[in] flush        Flush after every single subtitle packet.
 * @param[out] best_i      index of stream to be muxed
 * @return  On error a negative or zero value is returned, on success 1 is returned.
 */
int ff_pes_find_beststream(AVFormatContext *ctx, int packet_size, int flush, int64_t *scr, int* best_i);

/**
 * Get total number of frames that have been muxed.
 * @param[in] ctx    the AVFormatContext
 * @param[in] stream the PES stream
 * @param[in] len    PES packet size
 * @return  the frame number to be muxed
 */
int ff_pes_get_nb_frames(AVFormatContext *ctx, PESStream *stream, int len);

/**
 * Mux streams into a PES packet.
 * @param [in]      ctx            the AVFormatContext which contains streams
 * @param [in]      stream_index   the stream index to write
 * @param [in]      pes_buffer     PES payload data
 * @param [in]      pts            packet presentation timestamp
 * @param [in]      dts            packet decoding timestamp
 * @param [in]      id             stream ID
 * @param [in]      start_code     PES packet start code
 * @param [in]      header_len     PES header size
 * @param [in]      packet_size    total packet size
 * @param [in]      payload_size   packet payload size
 * @param [in]      stuffing_size  packet stuffing size
 * @return   bytes written to PES stream.
 */
int ff_pes_muxer_write(AVFormatContext *ctx, int stream_index, uint8_t *pes_buffer,
          int64_t pts,int64_t dts, int  id, int startcode,
          uint8_t* pes_content, int pes_content_len,
          int header_len, int packet_size, int payload_size, int stuffing_size);

/**
 * Remove decoded packets of each stream.
 * @param[in] ctx  the AVFormatContext
 * @param[in] scr  System Clock Reference of PES stream
 * @return  On error a negative or zero value is returned, on success 1 is returned.
 */
int ff_pes_remove_decoded_packets(AVFormatContext *ctx, int64_t scr);

#endif/* AVFORMAT_MPEG_PES_H */
