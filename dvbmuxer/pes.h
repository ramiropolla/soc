/*
 * Copyright (c) 2000, 2001, 2002 Fabrice Bellard.
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
 * @file pes.h
 * PES packetizer api header.
 */

#ifndef PES_H
#define PES_H

#include "avformat.h"
#include "fifo.h"

typedef struct PacketDesc {
    int64_t pts;
    int64_t dts;
    int size;
    int unwritten_size;
    int flags;
    struct PacketDesc *next;
} PacketDesc;

typedef struct {
    int is_ps;  /*< whether it is a Program Stream */
    int packet_number;
} PESContext;

typedef struct {
    AVFifoBuffer fifo;
    int max_buffer_size; /*< in bytes */
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

/* mpeg2 */
#define PROGRAM_STREAM_MAP 0x1bc
#define PRIVATE_STREAM_1   0x1bd
#define PADDING_STREAM     0x1be
#define PRIVATE_STREAM_2   0x1bf


static const int lpcm_freq_tab[4] = { 48000, 96000, 44100, 32000 };


/*
 * Initialization of PES mux.
 * @param[in] ctx the AVFormatContext which contains streams
 * @return  On error a negative value is returned, on success zero.
 */
int ff_pes_mux_init(AVFormatContext *ctx);

/*
 * Finalization of PES mux.
 * @param [in] ctx the AVFormatContext which contains streams.
 * @return  NULL
 */
void ff_pes_mux_end(AVFormatContext *ctx);

/*
 * Write packet into PES fifo.
 * @param [in] ctx  the AVFormatContext which contains streams.
 * @param [in] pkt  the packet to write.
 * @return  NULL
 */
void ff_pes_write_packet(AVFormatContext *ctx, AVPacket *pkt);

/*
 * Find the most fit stream to be muxed.
 * @param[in] ctx   the AVFormatContext
 * @param[in] packet_size  the packet size of PES stream
 * @param[in] flush   whether we flush after every single subtitle packet for subtitle
 * @param[out] best_i       the best fit stream index
 * @return  On error a negative or zero value is returned, on success 1 is returned
 */
int ff_pes_find_beststream(AVFormatContext *ctx, int packet_size, int flush, int64_t scr, int* best_i);

/*
 * Get how many frames is muxed.
 * @param[in] ctx    the AVFormatContext
 * @param[in] stream the PES stream
 * @param[in] len    PES packet size
 * @return  the frame number to be muxed
 */
int ff_get_nb_frames(AVFormatContext *ctx, PESStream *stream, int len);

/*
 * Mux streams into a PES packet.
 * @param [in]      ctx            the AVFormatContext which contains streams
 * @param [in]      stream_index   the stream index to write
 * @param [in]      pts            packet presentation time stamp
 * @param [in]      dts            packet decoding time stamp
 * @param [in]      id             stream id
 * @param [in]      start_code     PES packet start code
 * @param [in]      header_len     PES header size
 * @param [in]      packet_size    the total packet size
 * @param [in]      payload_size   the payload size of the packet
 * @param [in]      stuffing_size  the stuffing size of the packet
 * @param [in]      trailer_size   the trailer size of the packet
 * @return   bytes wirtten to PES stream.
 */
int ff_pes_mux_write(AVFormatContext *ctx, int stream_index,
          int64_t pts,int64_t dts, int  id, int startcode,
          uint8_t* pes_content, int pes_content_len,
          int header_len, int packet_size, int payload_size, int stuffing_size, int tailer_size);

/*
 * Remove decoded packets of each stream.
 * @param[in] ctx  the AVFormatContext
 * @param[in] scr  System Clock Reference of PES stream.
 * @return  On error a negative or zero value is returned, on success 1 is returned
 */
int ff_pes_remove_decoded_packets(AVFormatContext *ctx, int64_t scr);


#endif/* PES_H */
