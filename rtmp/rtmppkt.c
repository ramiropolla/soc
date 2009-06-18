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

static const int8_t amf_tag_sizes[17] = {
    8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 0, 0, 0, 0
};

int rtmp_amf_tag_size(int type, const void *data)
{
    if (type == AMF_STRING_IN_OBJECT || type == AMF_STRING)
        return strlen(data) + 2;
    if (type > AMF_TYPED_OBJECT)
        return -1;
    if (type == AMF_LONG_STRING)
        return strlen(data) + 4;
    return amf_tag_sizes[type];
}

void rtmp_amf_write_tag(uint8_t **dst, AMFType type, const void *data)
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

int rtmp_packet_read(AVFormatContext *ctx, URLContext *h, RTMPPacket *p,
                     RTMPPacketHistory *hist)
{
    uint8_t hdr, t, buf[16];
    int stream_id, timestamp, data_size, offset = 0;
    uint8_t type;

    if (url_read(h, &hdr, 1) != 1) {
        av_log(ctx, AV_LOG_ERROR, "Cannot read packet header\n");
        return -1;
    }
    stream_id = hdr & 0x3F;

    hdr >>= 6;
    if (hdr == RTMP_PS_ONEBYTE) {
        av_log(ctx, AV_LOG_ERROR, "Onebyte header!\n\n\n");
        //todo
        return -1;
    } else {
        if (url_read_complete(h, buf, 3) != 3) {
            av_log(ctx, AV_LOG_ERROR, "reading timestamp failed\n");
            return -1;
        }
        timestamp = AV_RB24(buf);
        if (hdr != RTMP_PS_FOURBYTES) {
            if (url_read_complete(h, buf, 3) != 3) {
                av_log(ctx, AV_LOG_ERROR, "reading size failed\n");
                return -1;
            }
            data_size = AV_RB24(buf);
            if (url_read_complete(h, &type, 1) != 1) {
                av_log(ctx, AV_LOG_ERROR, "reading type failed\n");
                return -1;
            }
            if (hdr == RTMP_PS_TWELVEBYTES)
                if (url_read(h, buf, 4) != 4) {
                    av_log(ctx, AV_LOG_ERROR, "reading timestamp2 failed\n");
                    return -1;
                }
        } else {
            av_log(ctx, AV_LOG_ERROR, "Fourbytes!\n\n\n");
            //todo
            return -1;
        }
    }
    rtmp_packet_create(p, stream_id, type, type, data_size);
    if (hdr == RTMP_PS_TWELVEBYTES)
        p->extra = AV_RL32(buf);
    while (data_size > 0) {
        int toread = FFMIN(data_size, hist->chunk_size[stream_id]);
        int r;
        if ((r = url_read_complete(h, p->data + offset, toread)) != toread) {
            av_log(ctx, AV_LOG_ERROR, "Need %d read %d\n", toread, r);
            rtmp_packet_destroy(p);
            return -1;
        }
        data_size -= hist->chunk_size[stream_id];
        offset    += hist->chunk_size[stream_id];
        if (data_size > 0) {
            url_read_complete(h, &t, 1); //marker
            if (t != (0xC0 + stream_id)) {
                av_log(ctx, AV_LOG_ERROR, "Expected marker %02X, got %02X\n", t, 0xC0 + stream_id);
                return -1;
            }
        }
    }
    return 0;
}

int rtmp_packet_write(AVFormatContext *ctx, URLContext *h, RTMPPacket *pkt,
                      RTMPPacketHistory *hist)
{
    uint8_t pkt_hdr[16], *p = pkt_hdr;
    int mode = RTMP_PS_TWELVEBYTES;
    int chunk_size = hist->chunk_size[pkt->stream_id];
    int off = 0;

    if (pkt->type != RTMP_PT_INVOKE)
        mode = RTMP_PS_EIGHTBYTES;
    bytestream_put_byte(&p, pkt->stream_id | (mode << 6));
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
            uint8_t marker = 0xC0 | pkt->stream_id;
            url_write(h, &marker, 1);
        }
    }
    return 0;
}

int rtmp_packet_create(RTMPPacket *pkt, int stream_id, RTMPPacketType type,
                       int timestamp, int size)
{
    pkt->data = av_malloc(size);
    if (!pkt->data)
        return -1;
    pkt->data_size = size;
    pkt->stream_id = stream_id;
    pkt->type      = type;
    pkt->timestamp = timestamp;
    pkt->extra     = 0;

    return 0;
}

void rtmp_packet_destroy(RTMPPacket *pkt)
{
    if (!pkt)
        return;
    av_freep(&pkt->data);
    pkt->data_size = 0;
}


static void parse_amf(const uint8_t *data, int size)
{
    const uint8_t *ptr = data, *end = data+size;
    int object = 0,array=0,i,len;
    char buf[65536];
    while(ptr < end){
        if(object){
            len=bytestream_get_be16(&ptr);
            if(len){
            for(i=0;i<object;i++)av_log(NULL,0,"    ");
            bytestream_get_buffer(&ptr, buf, len);
            buf[len] = 0;
            av_log(NULL,0,"%s: ",buf);
            }
        }
        switch (*ptr++) {
        case AMF_NUMBER:
            {
                double d;
                d = av_int2dbl(bytestream_get_be64(&ptr));
                av_log(NULL,0,"Number %g\n",d);
            }
            break;
        case AMF_BOOLEAN:
            av_log(NULL,0,"Boolean %d\n",*ptr++);
            break;
        case AMF_STRING:
            len=bytestream_get_be16(&ptr);
            bytestream_get_buffer(&ptr, buf, len);
            buf[len] = 0;
            av_log(NULL,0,"String '%s'\n",buf);
            break;
        case AMF_NULL:
            av_log(NULL,0,"NULL\n");
            break;
        case AMF_OBJECT:
            for(i=0;i<object;i++)av_log(NULL,0,"    ");
            av_log(NULL,0,"Object{\n");
            object++;
            break;
        case AMF_ECMA_ARRAY:
            array = bytestream_get_be32(&ptr);
            av_log(NULL,0,"Array of %d elements: [\n",array);
            object++;
            break;
        case AMF_OBJECT_END:
            object--;
            for(i=0;i<object;i++)av_log(NULL,0,"    ");
            if(array){
            array = 0;
            av_log(NULL,0,"]\n");
            }else
            av_log(NULL,0,"}\n");
            break;
        default:
            av_log(NULL,0,"Type %02X\n",ptr[-1]);
            return;
        }
    }
}

void rtmp_packet_inspect(AVFormatContext *ctx, RTMPPacket *pkt)
{
    av_log(NULL,0,"Packet on ");
    switch (pkt->stream_id) {
    case RTMP_NETWORK_CHANNEL: av_log(NULL,0,"network channel");break;
    case RTMP_SYSTEM_CHANNEL:  av_log(NULL,0,"system channel");break;
    case RTMP_VIDEO_CHANNEL:   av_log(NULL,0,"video channel");break;
    case RTMP_AUDIO_CHANNEL:   av_log(NULL,0,"audio channel");break;
    default:                   av_log(NULL,0,"channel %d",pkt->stream_id);
    }
    av_log(NULL,0," type ");
    switch (pkt->type) {
    case RTMP_PT_CHUNK_SIZE:   av_log(NULL,0,"chunk size %d",AV_RB32(pkt->data));break;
    case RTMP_PT_BYTES_READ:   av_log(NULL,0,"bytes read");break;
    case RTMP_PT_PING:         av_log(NULL,0,"ping type = %d", AV_RB16(pkt->data));break;
    case RTMP_PT_SERVER_BW:    av_log(NULL,0,"server BW=%d",AV_RB32(pkt->data));break;
    case RTMP_PT_CLIENT_BW:    av_log(NULL,0,"client BW=%d",AV_RB32(pkt->data));break;
    case RTMP_PT_AUDIO:        av_log(NULL,0,"audio");break;
    case RTMP_PT_VIDEO:        av_log(NULL,0,"video");break;
    case RTMP_PT_NOTIFY:       av_log(NULL,0,"notify");break;
    case RTMP_PT_INVOKE:       av_log(NULL,0,"invoke");break;
    case RTMP_PT_METADATA:     av_log(NULL,0,"metadata");break;
    default:                   av_log(NULL,0,"%X",pkt->type);
    }
    av_log(NULL,0," ts %d/%d size %d\n", pkt->timestamp, pkt->extra, pkt->data_size);
    if (pkt->type == RTMP_PT_INVOKE || pkt->type == RTMP_PT_NOTIFY
     || pkt->type == RTMP_PT_METADATA)
        parse_amf(pkt->data, pkt->data_size);
}
