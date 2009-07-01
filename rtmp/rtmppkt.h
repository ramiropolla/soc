/*
 * RTMP packet utilities
 * Copyright (c) 2009 Kostya Shishkov
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

#ifndef AVFORMAT_RTMPPKT_H
#define AVFORMAT_RTMPPKT_H

#include "avformat.h"

/* maximum possible number of different RTMP channels */
#define RTMP_CHANNELS 64

enum RTMPChannel {
    RTMP_NETWORK_CHANNEL = 2,   ///< channel for network-related messages (bandwidth report, ping, etc)
    RTMP_SYSTEM_CHANNEL,        ///< channel for sending server control messages
    RTMP_VIDEO_CHANNEL = 8,     ///< channel for video data
    RTMP_AUDIO_CHANNEL,         ///< channel for audio data
};

typedef enum RTMPPacketType {
    RTMP_PT_CHUNK_SIZE   =  1,  ///< chunk size change
    RTMP_PT_BYTES_READ   =  3,  ///< number of bytes read
    RTMP_PT_PING,               ///< ping
    RTMP_PT_SERVER_BW,          ///< server bandwidth
    RTMP_PT_CLIENT_BW,          ///< client bandwidth
    RTMP_PT_AUDIO        =  8,  ///< audio packet
    RTMP_PT_VIDEO,              ///< video packet
    RTMP_PT_FLEX_STREAM  = 15,  ///< Flex shared stream
    RTMP_PT_FLEX_OBJECT,        ///< Flex shared object
    RTMP_PT_FLEX_MESSAGE,       ///< Flex shared message
    RTMP_PT_NOTIFY,             ///< some notification
    RTMP_PT_SHARED_OBJ,         ///< shared object
    RTMP_PT_INVOKE,             ///< invoke some stream action
    RTMP_PT_METADATA     = 22,  ///< FLV metadata
} RTMPPacketType;

enum RTMPPacketSize {
    RTMP_PS_TWELVEBYTES = 0,
    RTMP_PS_EIGHTBYTES,
    RTMP_PS_FOURBYTES,
    RTMP_PS_ONEBYTE
};

typedef enum AMFType {
    AMF_NUMBER = 0,
    AMF_BOOLEAN,
    AMF_STRING,
    AMF_OBJECT,
    AMF_MOVIE,
    AMF_NULL,
    AMF_UNDEFINED,
    AMF_REFERENCE,
    AMF_ECMA_ARRAY,
    AMF_OBJECT_END,
    AMF_STRICT_ARRAY,
    AMF_DATE,
    AMF_LONG_STRING,
    AMF_UNSUPPORTED,
    AMD_RECORD_SET,
    AMF_XML_OBJECT,
    AMF_TYPED_OBJECT,

    AMF_STRING_IN_OBJECT = 99,
} AMFType;

/**
 * structure for holding RTMP packets
 */
typedef struct RTMPPacket {
    uint8_t        channel_id; ///< RTMP channel ID
    RTMPPacketType type;       ///< packet type
    int            timestamp;  ///< packet timestamp
    int            extra;      ///< additional data
    uint8_t        *data;      ///< packet payload
    int            data_size;  ///< packet payload size
} RTMPPacket;

int rtmp_packet_create(RTMPPacket *pkt, int channel_id, RTMPPacketType type,
                       int timestamp, int size);

void rtmp_packet_destroy(RTMPPacket *pkt);

int rtmp_packet_read(URLContext *h, RTMPPacket *p,
                     int chunk_size, RTMPPacket *prev_pkt);

int rtmp_packet_write(URLContext *h, RTMPPacket *p,
                      int chunk_size, RTMPPacket *prev_pkt);

int rtmp_amf_skip_data(const uint8_t *data);

int rtmp_amf_find_field(const uint8_t *data, const uint8_t *name,
                        uint8_t *dst, int dst_size);

void rtmp_amf_write_tag(uint8_t **dst, AMFType type, const void *data);

void rtmp_packet_inspect(RTMPPacket *pkt);

#endif /* AVFORMAT_RTMPPKT_H */
