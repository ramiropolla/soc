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
 * @file mpegpes.h
 * MPEG PES packetizer API header
 */

#ifndef AVFORMAT_MPEG_PES_H
#define AVFORMAT_MPEG_PES_H

#include "avformat.h"
#include "fifo.h"

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
    int format;
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
} StreamInfo;


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
 * @return  the number of frames have been muxed.
 */
int ff_pes_get_nb_frames(AVFormatContext *ctx, StreamInfo *stream, int len);


/**
 * Caculate the PES header size
 * @param[in] id                stream id
 * @param[in] stream            pes stream
 * @param[in] packet_size       pes packet size
 * @param[in] header_len        pes header length
 * @param[in] pts               current pts
 * @param[in] dts               current dts
 * @param[in] payload_size      pes payload size
 * @param[in] startcode         pes startcode
 * @param[in] stuffing_size     pes stuffing size
 * @param[in] trailer_size      unwritten trailer size
 * @param[in] pad_packet_bytes  padding size for packet
 */
void ff_pes_cal_header(int id, StreamInfo *stream,
          int *packet_size,  int *header_len, int64_t *pts,int64_t *dts,
          int *payload_size, int *startcode, int *stuffing_size,
          int *trailer_size, int *pad_packet_bytes);

/**
 * Mux one stream into PES stream.
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
