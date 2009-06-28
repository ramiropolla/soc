/*
 * IEC958 muxer
 * Copyright (c) 2009 Bartlomiej Wolowiec
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

#include "avformat.h"
#include "libavcodec/ac3.h"
#include "libavcodec/dca.h"

#define SYNCWORD1 0xF872
#define SYNCWORD2 0x4E1F
#define BURST_HEADER_SIZE 0x8

#define IEC958_AC3                0x01
#define IEC958_MPEG1_LAYER1       0x04
#define IEC958_MPEG1_LAYER23      0x05
//#define IEC958_MPEG2_NO_EXT       0x05 /* No extension */
//#define IEC958_MPEG2_EXT          0x06 /* With extension */
//#define IEC958_MPEG2_AAC          0x07
//#define IEC958_MPEG2_LAYER1_LSW   0x08 /* Low Sampling Frequency */
//#define IEC958_MPEG2_LAYER2_LSW   0x09 /* Low Sampling Frequency */
//#define IEC958_MPEG2_LAYER3_LSW   0x0A /* Low Sampling Frequency */
#define IEC958_DTS1               0x0B
#define IEC958_DTS2               0x0C
#define IEC958_DTS3               0x0D
//#define IEC958_EAC3               0x15

typedef struct IEC958Context{
    int data_type;
    int pkt_size;
    int pkt_offset; // bytes
    int (*header_info)(AVFormatContext *s, AVPacket *pkt);
} IEC958Context;

static int spdif_header_ac3(AVFormatContext *s, AVPacket *pkt){
    IEC958Context *ctx = s->priv_data;
    int bitstream_mode = pkt->data[6] & 0x7;

    ctx->data_type = IEC958_AC3 | (bitstream_mode << 8);
    ctx->pkt_offset = AC3_FRAME_SIZE<<2;
    return 0;
}

static int spdif_header_dts(AVFormatContext *s, AVPacket *pkt){
    IEC958Context *ctx = s->priv_data;
    uint32_t syncword_dts = (pkt->data[0] << 24) | (pkt->data[1]<<16) | (pkt->data[2]<<8) | pkt->data[3];
    int samples;

    if(syncword_dts != DCA_MARKER_RAW_BE){
        av_log(NULL, AV_LOG_ERROR, "bad DTS syncword\n");
        return -1;
    }
    samples = ((((pkt->data[4] & 0x01) << 6) | (pkt->data[5] >> 2)) + 1) << 5; // :)
    av_log(NULL, AV_LOG_DEBUG, "samples=%i\n", samples);
    switch(samples){
        case 512:
            ctx->data_type = IEC958_DTS1;
            break;
        case 1024:
            ctx->data_type = IEC958_DTS2;
            break;
        case 2048:
            ctx->data_type = IEC958_DTS3;
            break;
        default:
            av_log(NULL, AV_LOG_ERROR, "%i samples in DTS frame not supported\n", samples);
            return -1;
    }
    ctx->pkt_offset = samples<<2;

    return 0;
}

static int spdif_header_mpeg1_layer1(AVFormatContext *s, AVPacket *pkt){
    IEC958Context *ctx = s->priv_data;

    ctx->data_type = IEC958_MPEG1_LAYER1;
    ctx->pkt_offset = 384<<2; //TODO
    return 0;
}

static int spdif_header_mpeg1_layer23(AVFormatContext *s, AVPacket *pkt){
    IEC958Context *ctx = s->priv_data;

    ctx->data_type = IEC958_MPEG1_LAYER23;
    // TODO Data type dependant info (normal/karaoke, dynamic range control)
    ctx->pkt_offset = 1152<<2; //TODO
    return 0;
}

static int spdif_write_header(AVFormatContext *s){
    IEC958Context *ctx = s->priv_data;

    switch(s->streams[0]->codec->codec_id){
        case CODEC_ID_AC3:
            ctx->header_info = spdif_header_ac3;
            break;
        case CODEC_ID_MP1:
            ctx->header_info = spdif_header_mpeg1_layer1;
            break;
        case CODEC_ID_MP2:
        case CODEC_ID_MP3:
            ctx->header_info = spdif_header_mpeg1_layer23;
            break;
        case CODEC_ID_DTS:
            ctx->header_info = spdif_header_dts;
            break;

        default:
            av_log(NULL, AV_LOG_ERROR, "codec not supported\n");
            return -1;
    }
    put_le16(s->pb, 0);
    put_le16(s->pb, 0);
    put_le16(s->pb, 0);
    put_le16(s->pb, 0);
    return 0;
}

static int spdif_write_packet(struct AVFormatContext *s, AVPacket *pkt){
    IEC958Context *ctx = s->priv_data;
    uint16_t *data = (uint16_t *)pkt->data;
    int i;

    ctx->pkt_size = pkt->size << 3;

    (*ctx->header_info)(s, pkt);

    put_le16(s->pb, SYNCWORD1);      //Pa
    put_le16(s->pb, SYNCWORD2);      //Pb
    put_le16(s->pb, ctx->data_type); //Pc
    put_le16(s->pb, ctx->pkt_size);  //Pd

    //put_buffer(s->pb, pkt->data, pkt->size);
    //XXX memcpy... ?
    for(i=0; i<pkt->size>>1; i++)
        put_be16(s->pb, data[i]); //XXX be?

    if(pkt->size&1)
        put_be16(s->pb, pkt->data[pkt->size-1]); //XXX be?

    i=(ctx->pkt_offset - BURST_HEADER_SIZE - pkt->size) >> 1;
    if(i < 0){
        av_log(NULL, AV_LOG_ERROR, "bitrate is too high\n");
        return -1;
    }

    for(; i>0; i--)
        put_le16(s->pb, 0);

    av_log(NULL, AV_LOG_DEBUG, "type=%x len=%i pkt_offset=%i\n", ctx->data_type, pkt->size, ctx->pkt_offset);

    put_flush_packet(s->pb);
    return 0;
}

AVOutputFormat spdif_muxer = {
    "spdif",
    NULL_IF_CONFIG_SMALL("IEC958 (IEC-61937)"),
    NULL,
    "spdif",
    sizeof(IEC958Context),
    CODEC_ID_AC3,
    CODEC_ID_NONE,
    spdif_write_header,
    spdif_write_packet,
    NULL,
    //.flags=
};

