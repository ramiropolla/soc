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
#include "libavcodec/aac_parser.h"

#define SYNCWORD1 0xF872
#define SYNCWORD2 0x4E1F
#define BURST_HEADER_SIZE 0x8

#define IEC958_AC3                0x01
#define IEC958_MPEG1_LAYER1       0x04
#define IEC958_MPEG1_LAYER23      0x05
//#define IEC958_MPEG2_EXT          0x06 /* With extension */
#define IEC958_MPEG2_AAC          0x07
#define IEC958_MPEG2_LAYER1_LSF   0x08  /* Low Sampling Frequency */
#define IEC958_MPEG2_LAYER2_LSF   0x09  /* Low Sampling Frequency */
#define IEC958_MPEG2_LAYER3_LSF   0x0A  /* Low Sampling Frequency */
#define IEC958_DTS1               0x0B
#define IEC958_DTS2               0x0C
#define IEC958_DTS3               0x0D
#define IEC958_MPEG2_AAC_LSF_2048 0x13
#define IEC958_MPEG2_AAC_LSF_4096 (0x13|0x20)
//#define IEC958_EAC3               0x15

typedef struct IEC958Context {
    int data_type;              ///< Burst info
    int pkt_size;               ///< Length code (number of bits or bytes - according to data_type)
    int pkt_offset;             ///< Repetition period of a data burst in bytes
    int (*header_info) (AVFormatContext *s, AVPacket *pkt);
} IEC958Context;

static int spdif_header_ac3(AVFormatContext *s, AVPacket *pkt)
{
    IEC958Context *ctx = s->priv_data;
    int bitstream_mode = pkt->data[6] & 0x7;

    ctx->data_type = IEC958_AC3 | (bitstream_mode << 8);
    ctx->pkt_offset = AC3_FRAME_SIZE << 2;
    return 0;
}

static int spdif_header_dts(AVFormatContext *s, AVPacket *pkt)
{
    IEC958Context *ctx = s->priv_data;
    uint32_t syncword_dts =
        (pkt->data[0] << 24) | (pkt->data[1] << 16) | (pkt-> data[2] << 8) | pkt->data[3];
    int samples;

    switch (syncword_dts) {
    case DCA_MARKER_RAW_BE:
        samples =
            ((((pkt->data[4] & 0x01) << 6) | (pkt->data[5] >> 2)) + 1) << 5;
        break;
    case DCA_MARKER_RAW_LE:
        samples =
            ((((pkt->data[5] & 0x01) << 6) | (pkt->data[4] >> 2)) + 1) << 5;
        break;
    case DCA_MARKER_14B_BE:
        samples =
            ((((pkt->data[5] & 0x07) << 4) | (pkt->data[6] & 0x3f)) >> 2) << 5;
        break;
    case DCA_MARKER_14B_LE:
        samples =
            ((((pkt->data[4] & 0x07) << 4) | (pkt->data[7] & 0x3f)) >> 2) << 5;
        break;
    default:
        av_log(s, AV_LOG_ERROR, "bad DTS syncword\n");
        return -1;
    }
    av_log(s, AV_LOG_DEBUG, "samples=%i\n", samples);
    switch (samples) {
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
        av_log(s, AV_LOG_ERROR, "%i samples in DTS frame not supported\n",
               samples);
        return -1;
    }
    ctx->pkt_offset = samples << 2;

    return 0;
}

static const uint8_t mpeg_data_type[2][3] = {
    //     LAYER1                      LAYER2                  LAYER3
    { IEC958_MPEG2_LAYER1_LSF, IEC958_MPEG2_LAYER2_LSF, IEC958_MPEG2_LAYER3_LSF },  //MPEG2 LSF
    { IEC958_MPEG1_LAYER1,     IEC958_MPEG1_LAYER23,    IEC958_MPEG1_LAYER23 },     //MPEG1
};

static const uint16_t mpeg_pkt_offset[2][3] = {
    //LAYER1  LAYER2  LAYER3
    { 768,    2304,   1152 }, // MPEG2 LSF
    { 384,    1152,   1152 }, // MPEG1
};

static int spdif_header_mpeg(AVFormatContext *s, AVPacket *pkt)
{
    IEC958Context *ctx = s->priv_data;
    int lsf = (pkt->data[1] >> 3) & 1;
    int layer = 3 - ((pkt->data[1] >> 1) & 3);

    av_log(s, AV_LOG_DEBUG, "lsf: %i layer: %i\n", lsf, layer);
    ctx->data_type = mpeg_data_type[lsf][layer];
    ctx->pkt_offset = mpeg_pkt_offset[lsf][layer] << 2;
    // TODO Data type dependant info (normal/karaoke, dynamic range control)
    return 0;
}

static int spdif_header_aac(AVFormatContext *s, AVPacket *pkt)
{
    IEC958Context *ctx = s->priv_data;
    AACADTSHeaderInfo hdr;
    GetBitContext gbc;
    int ret;

    init_get_bits(&gbc, pkt->data, AAC_ADTS_HEADER_SIZE * 8);
    ret = ff_aac_parse_header(&gbc, &hdr);

    ctx->pkt_offset = hdr.samples << 2;
    switch (hdr.num_aac_frames) {
    case 1:
        ctx->data_type = IEC958_MPEG2_AAC;
        break;
    case 2:
        ctx->data_type = IEC958_MPEG2_AAC_LSF_2048;
        break;
    case 4:
        ctx->data_type = IEC958_MPEG2_AAC_LSF_4096;
        break;
    default:
        av_log(s, AV_LOG_ERROR, "%i samples in AAC frame not supported\n",
               hdr.samples);
        return -1;
    }
    //TODO Data type dependent info (LC profile/SBR)
    return 0;
}

static int spdif_write_header(AVFormatContext *s)
{
    IEC958Context *ctx = s->priv_data;

    switch (s->streams[0]->codec->codec_id) {
    case CODEC_ID_AC3:
        ctx->header_info = spdif_header_ac3;
        break;
    case CODEC_ID_MP1:
    case CODEC_ID_MP2:
    case CODEC_ID_MP3:
        ctx->header_info = spdif_header_mpeg;
        break;
    case CODEC_ID_DTS:
        ctx->header_info = spdif_header_dts;
        break;
    case CODEC_ID_AAC:
        ctx->header_info = spdif_header_aac;
        break;
    default:
        av_log(s, AV_LOG_ERROR, "codec not supported\n");
        return -1;
    }
    put_le16(s->pb, 0);
    put_le16(s->pb, 0);
    put_le16(s->pb, 0);
    put_le16(s->pb, 0);
    return 0;
}

static int spdif_write_packet(struct AVFormatContext *s, AVPacket *pkt)
{
    IEC958Context *ctx = s->priv_data;
    uint16_t *data = (uint16_t *) pkt->data;
    int i;

    ctx->pkt_size = ((pkt->size + 1) >> 1) << 4;
    (*ctx->header_info) (s, pkt);

    put_le16(s->pb, SYNCWORD1);      //Pa
    put_le16(s->pb, SYNCWORD2);      //Pb
    put_le16(s->pb, ctx->data_type); //Pc
    put_le16(s->pb, ctx->pkt_size);  //Pd

    //XXX memcpy... ?
    for (i = 0; i < pkt->size >> 1; i++)
        put_be16(s->pb, data[i]);

    if (pkt->size & 1)
        put_be16(s->pb, pkt->data[pkt->size - 1]);

    i = (ctx->pkt_offset - BURST_HEADER_SIZE - pkt->size) >> 1;
    if (i < 0) {
        av_log(s, AV_LOG_ERROR, "bitrate is too high\n");
        return -1;
    }

    for (; i > 0; i--)
        put_le16(s->pb, 0);

    av_log(s, AV_LOG_DEBUG, "type=%x len=%i pkt_offset=%i\n",
           ctx->data_type, pkt->size, ctx->pkt_offset);

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
