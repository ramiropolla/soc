/*
 * VC1 Test Bitstreams Format Demuxer
 * Copyright (c) 2006 Konstantin Shishkov
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/**
 * @file vc1test.c
 * VC1 test bitstreams file demuxer
 * by Konstantin Shishkov
 * Format specified in SMPTE standard 421 Annex L
 */

#include "avformat.h"

typedef struct VC1testDemuxContext {
    int frame_count;
    int current_frame;
    uint8_t struct_c[4];
    uint32_t width, height;
    int fps;
    int stream_index;
} VC1testDemuxContext;

static int vc1t_probe(AVProbeData *p)
{
    if (p->buf_size < 9*4)
        return 0;

    if (p->buf[3] != 0xC5)
        return 0;

    if (LE_32(&p->buf[4]) != 4)
        return 0;

    return AVPROBE_SCORE_MAX;
}

static int vc1t_read_header(AVFormatContext *s,
                           AVFormatParameters *ap)
{
    VC1testDemuxContext *vcd = (VC1testDemuxContext *)s->priv_data;
    ByteIOContext *pb = &s->pb;
    AVStream *st;
    uint32_t value;

    value = get_le32(pb);

    if((value >> 24) != 0xC5){
        av_log(s, AV_LOG_ERROR, "Not a VC1 Simple/Main profile bitstream!\n");
        return -1;
    }
    vcd->frame_count = value & 0xFFFFFF;
    vcd->current_frame = 0;

    value = get_le32(pb);
    if(value != 4){
        av_log(s, AV_LOG_ERROR, "Not a VC1 Simple/Main profile bitstream!\n");
        return -1;
    }

    get_buffer(pb, vcd->struct_c, 4);

    //profile = vcd->struct_c
    vcd->height = get_le32(pb);
    vcd->width = get_le32(pb);

    value = get_le32(pb);
    if(value != 0xC){
        av_log(s, AV_LOG_ERROR, "Not a VC1 Simple/Main profile bitstream!\n");
        return -1;
    }

    value = get_le32(pb);
    value = get_le32(pb);
    vcd->fps = get_le32(pb);

    /* init video codec */
    st = av_new_stream(s, 0);
    if (!st)
        return -1;

    st->codec->width = vcd->width;
    st->codec->height = vcd->height;
    st->codec->codec_type = CODEC_TYPE_VIDEO;
    st->codec->codec_id = CODEC_ID_WMV3;

    st->codec->extradata = av_malloc(4);
    st->codec->extradata_size = 4;
    memcpy(st->codec->extradata, vcd->struct_c, 4);
    av_set_pts_info(st, 33, 1, vcd->fps);
    vcd->stream_index = st->index;

    return 0;
}

static int vc1t_read_packet(AVFormatContext *s,
                           AVPacket *pkt)
{
    VC1testDemuxContext *vcd = (VC1testDemuxContext *)s->priv_data;
    ByteIOContext *pb = &s->pb;
    int ret = 0;
    uint32_t value;
    int frame_size;

    if (vcd->current_frame >= vcd->frame_count)
        return AVERROR_IO;

    value = get_le32(pb);
    frame_size = value & 0xFFFFFF;
    value = get_le32(pb);
    if (av_new_packet(pkt, frame_size))
        return AVERROR_NOMEM;

    ret = get_buffer(pb, pkt->data, frame_size);

    if (ret != frame_size) {
        av_free_packet(pkt);
        ret = AVERROR_IO;
    }
    pkt->stream_index = vcd->stream_index;
    vcd->current_frame++;

    return ret;
}

static int vc1t_read_close(AVFormatContext *s)
{
//    VC1testDemuxContext *vcd = (VC1testDemuxContext *)s->priv_data;

    return 0;
}

static AVInputFormat vc1t_iformat = {
    "vc1test",
    "VC1 bitstream format",
    sizeof(VC1testDemuxContext),
    vc1t_probe,
    vc1t_read_header,
    vc1t_read_packet,
    vc1t_read_close,
};

int vc1t_init(void)
{
    av_register_input_format(&vc1t_iformat);
    return 0;
}
