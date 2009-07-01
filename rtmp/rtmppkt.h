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

/**
 * channels used to for RTMP packets with different purposes (i.e. data, network
 * control, remote procedure calls, etc.)
 */
enum RTMPChannel {
    RTMP_NETWORK_CHANNEL = 2,   ///< channel for network-related messages (bandwidth report, ping, etc)
    RTMP_SYSTEM_CHANNEL,        ///< channel for sending server control messages
    RTMP_VIDEO_CHANNEL = 8,     ///< channel for video data
    RTMP_AUDIO_CHANNEL,         ///< channel for audio data
};

/**
 * known RTMP packet types
 */
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

/**
 * possible RTMP packet header sizes
 */
enum RTMPPacketSize {
    RTMP_PS_TWELVEBYTES = 0, ///< packet has 12-byte header
    RTMP_PS_EIGHTBYTES,      ///< packet has 8-byte header
    RTMP_PS_FOURBYTES,       ///< packet has 4-byte header
    RTMP_PS_ONEBYTE          ///< packet is really a next chunk of a packet
};

/**
 * AMF types used in RTMP packets
 */
typedef enum AMFType {
    AMF_NUMBER = 0,   ///< number (double precision)
    AMF_BOOLEAN,      ///< boolean value
    AMF_STRING,       ///< Pascal-style string with length < 65536
    AMF_OBJECT,       ///< AMF object, contains property names and values
    AMF_MOVIE,        ///< Flash object
    AMF_NULL,         ///< NULL value
    AMF_UNDEFINED,    ///< undefined (return?) value
    AMF_REFERENCE,    ///< reference
    AMF_ECMA_ARRAY,   ///< ECMA array, almost like AMF object but has number of entries
    AMF_OBJECT_END,   ///< marker for end of AMF object or ECMA array
    AMF_STRICT_ARRAY, ///< strict array
    AMF_DATE,         ///< date
    AMF_LONG_STRING,  ///< Pascal-style string with possible length up to 4GB
    AMF_UNSUPPORTED,  ///< unsipported feature indicator
    AMD_RECORD_SET,   ///< record set
    AMF_XML_OBJECT,   ///< XML object
    AMF_TYPED_OBJECT, ///< typed object

    AMF_STRING_IN_OBJECT = 99, ///< internal type used for AMF object field names
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

/**
 * Creates new RTMP packet with given attributes.
 *
 * @param pkt        packet
 * @param channel_id packet channel ID
 * @param type       packet type
 * @param timestamp  packet timestamp
 * @param size       packet size
 * @return zero on success, -1 otherwise
 */
int rtmp_packet_create(RTMPPacket *pkt, int channel_id, RTMPPacketType type,
                       int timestamp, int size);

/**
 * Frees RTMP packet.
 *
 * @param pkt packet
 */
void rtmp_packet_destroy(RTMPPacket *pkt);

/**
 * Reads RTMP packet.
 *
 * @param h          reader context
 * @param p          packet
 * @param chunk_size current chunk size
 * @param prev_pkt   previously read packet headers for all channels
 *                   (may be needed for restoring incomplete packet header)
 * @return zero on success, -1 otherwise
 */
int rtmp_packet_read(URLContext *h, RTMPPacket *p,
                     int chunk_size, RTMPPacket *prev_pkt);

/**
 * Sends RTMP packet.
 *
 * @param h          reader context
 * @param p          packet to send
 * @param chunk_size current chunk size
 * @param prev_pkt   previously sent packet headers for all channels
 *                   (may be used for packet header compressing)
 * @return zero on success, -1 otherwise
 */
int rtmp_packet_write(URLContext *h, RTMPPacket *p,
                      int chunk_size, RTMPPacket *prev_pkt);

/**
 * Calculates number of bytes needed to skip first AMF entry in data.
 *
 * @param data input data
 * @return number of bytes used by first AMF entry
 */
int rtmp_amf_skip_data(const uint8_t *data);

/**
 * Retrieves value of given AMF object field in string form.
 *
 * @param data     AMF object data
 * @param name     name of field to retrieve
 * @param dst      buffer for storing result
 * @param dst_size output buffer size
 * @return 0 if search and retrieval succeeded, -1 otherwise
 */
int rtmp_amf_find_field(const uint8_t *data, const uint8_t *name,
                        uint8_t *dst, int dst_size);

/**
 * Write AMF tag to buffer.
 *
 * @param dst  pointer to the input buffer (will be modified)
 * @param type tag type
 * @param data optional tag value
 */
void rtmp_amf_write_tag(uint8_t **dst, AMFType type, const void *data);

void rtmp_packet_inspect(RTMPPacket *pkt);

#endif /* AVFORMAT_RTMPPKT_H */
