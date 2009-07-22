/*
 * RTMP input format
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

/* needed for gethostname() */
#define _XOPEN_SOURCE 600

#include "libavcodec/bytestream.h"
#include "libavutil/avstring.h"
#include "avformat.h"

#include "rtmppkt.h"

void ff_amf_write_tag(uint8_t **dst, AMFType type, const void *data)
{
    if (type != AMF_OBJECT_END && type != AMF_STRING_IN_OBJECT)
        bytestream_put_byte(dst, type);
    switch(type){
    case AMF_NUMBER:
        bytestream_put_be64(dst, av_dbl2int(*(const double*)data));
        break;
    case AMF_BOOLEAN:
        bytestream_put_byte(dst, *(const uint8_t*)data);
        break;
    case AMF_STRING:
    case AMF_STRING_IN_OBJECT:
        bytestream_put_be16(dst, strlen(data));
        bytestream_put_buffer(dst, data, strlen(data));
        break;
    case AMF_OBJECT_END:
        bytestream_put_be24(dst, AMF_OBJECT_END);
        break;
    case AMF_STRICT_ARRAY:
        bytestream_put_be32(dst, *(const uint32_t*)data);
        break;
    }
}

int ff_rtmp_packet_read(URLContext *h, RTMPPacket *p,
                        int chunk_size, RTMPPacket *prev_pkt)
{
    uint8_t hdr, t, buf[16];
    int channel_id, timestamp, data_size, offset = 0, extra = 0;
    uint8_t type;

    if (url_read(h, &hdr, 1) != 1) {
        return AVERROR(EIO);
    }
    channel_id = hdr & 0x3F;

    hdr >>= 6;
    if (hdr == RTMP_PS_ONEBYTE) {
        //todo
        return -1;
    } else {
        if (url_read_complete(h, buf, 3) != 3)
            return AVERROR(EIO);
        timestamp = AV_RB24(buf);
        if (hdr != RTMP_PS_FOURBYTES) {
            if (url_read_complete(h, buf, 3) != 3)
                return AVERROR(EIO);
            data_size = AV_RB24(buf);
            if (url_read_complete(h, &type, 1) != 1)
                return AVERROR(EIO);
            if (hdr == RTMP_PS_TWELVEBYTES) {
                if (url_read_complete(h, buf, 4) != 4)
                    return AVERROR(EIO);
                extra = AV_RL32(buf);
            } else {
                extra = prev_pkt[channel_id].extra;
            }
        } else {
            data_size = prev_pkt[channel_id].data_size;
            type      = prev_pkt[channel_id].type;
            extra     = prev_pkt[channel_id].extra;
        }
    }
    if (ff_rtmp_packet_create(p, channel_id, type, timestamp, data_size))
        return -1;
    p->extra = extra;
    // save history
    prev_pkt[channel_id].channel_id = channel_id;
    prev_pkt[channel_id].type       = type;
    prev_pkt[channel_id].data_size  = data_size;
    prev_pkt[channel_id].timestamp  = timestamp;
    prev_pkt[channel_id].extra      = extra;
    while (data_size > 0) {
        int toread = FFMIN(data_size, chunk_size);
        int r;
        if ((r = url_read_complete(h, p->data + offset, toread)) != toread) {
            ff_rtmp_packet_destroy(p);
            return AVERROR(EIO);
        }
        data_size -= chunk_size;
        offset    += chunk_size;
        if (data_size > 0) {
            url_read_complete(h, &t, 1); //marker
            if (t != (0xC0 + channel_id))
                return -1;
        }
    }
    return 0;
}

int ff_rtmp_packet_write(URLContext *h, RTMPPacket *pkt,
                         int chunk_size, RTMPPacket *prev_pkt)
{
    uint8_t pkt_hdr[16], *p = pkt_hdr;
    int mode = RTMP_PS_TWELVEBYTES;
    int off = 0;

    //TODO: header compression
    bytestream_put_byte(&p, pkt->channel_id | (mode << 6));
    if (mode != RTMP_PS_ONEBYTE) {
        bytestream_put_be24(&p, pkt->timestamp);
        if (mode != RTMP_PS_FOURBYTES) {
            bytestream_put_be24(&p, pkt->data_size);
            bytestream_put_byte(&p, pkt->type);
            if (mode == RTMP_PS_TWELVEBYTES)
                bytestream_put_le32(&p, pkt->extra);
        }
    }
    url_write(h, pkt_hdr, p-pkt_hdr);
    while (off < pkt->data_size) {
        int towrite = FFMIN(chunk_size, pkt->data_size - off);
        url_write(h, pkt->data + off, towrite);
        off += towrite;
        if (off < pkt->data_size) {
            uint8_t marker = 0xC0 | pkt->channel_id;
            url_write(h, &marker, 1);
        }
    }
    return 0;
}

int ff_rtmp_packet_create(RTMPPacket *pkt, int channel_id, RTMPPacketType type,
                          int timestamp, int size)
{
    pkt->data = av_malloc(size);
    if (!pkt->data)
        return AVERROR(ENOMEM);
    pkt->data_size  = size;
    pkt->channel_id = channel_id;
    pkt->type       = type;
    pkt->timestamp  = timestamp;
    pkt->extra      = 0;

    return 0;
}

void ff_rtmp_packet_destroy(RTMPPacket *pkt)
{
    if (!pkt)
        return;
    av_freep(&pkt->data);
    pkt->data_size = 0;
}

int ff_amf_skip_data(const uint8_t *data)
{
    const uint8_t *base = data;

    switch (*data++) {
    case AMF_NUMBER:      return 9;
    case AMF_BOOLEAN:     return 2;
    case AMF_STRING:      return 3 + AV_RB16(data);
    case AMF_LONG_STRING: return 5 + AV_RB32(data);
    case AMF_NULL:        return 1;
    case AMF_ECMA_ARRAY:
        data += 4;
    case AMF_OBJECT:
        for (;;) {
            int size = bytestream_get_be16(&data);
            if (!size) {
                data++;
                break;
            }
            data += size;
            data += ff_amf_skip_data(data);
        }
        return data - base;
    case AMF_OBJECT_END:  return 1;
    default:              return -1;
    }
}

int ff_amf_find_field(const uint8_t *data, const uint8_t *name,
                      uint8_t *dst, int dst_size)
{
    int namelen = strlen(name);
    int len;

    if (*data++ != AMF_OBJECT)
        return -1;
    for (;;) {
        int size = bytestream_get_be16(&data);
        if (!size)
            break;
        data += size;
        if (size == namelen && !memcmp(data-size, name, namelen)) {
            switch (*data++) {
            case AMF_NUMBER:
                snprintf(dst, dst_size, "%g", av_int2dbl(AV_RB64(data)));
                return 0;
            case AMF_BOOLEAN:
                snprintf(dst, dst_size, "%s", *data ? "true" : "false");
                return 0;
            case AMF_STRING:
                len = bytestream_get_be16(&data);
                av_strlcpy(dst, data, FFMIN(len+1, dst_size));
                return 0;
            default:
                return -1;
            }
        }
        len = ff_amf_skip_data(data);
        if (len < 0)
            return -1;
        data += len;
    }
    return -1;
}
