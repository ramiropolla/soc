/*
 * MPEG PES muxer
 * Copyright (c) 2000-2002 Fabrice Bellard
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

#include "mpegpes.h"
#include "mpeg.h"
#include "libavcodec/bytestream.h"

int ff_pes_muxer_init(AVFormatContext *ctx)
{
    int i;

    for(i=0;i<ctx->nb_streams;i++) {
        AVStream *st = ctx->streams[i];
        StreamInfo *stream = st->priv_data;
        av_set_pts_info(st, 64, 1, 90000);

        switch(st->codec->codec_type) {
        case CODEC_TYPE_AUDIO:
            /* This value HAS to be used for VCD (see VCD standard, p. IV-7).
               Right now it is also used for everything else.*/
            stream->max_buffer_size = 4 * 1024;
            break;
        case CODEC_TYPE_VIDEO:
            if (st->codec->rc_buffer_size)
                stream->max_buffer_size = 6*1024 + st->codec->rc_buffer_size/8;
            else
                stream->max_buffer_size = 230*1024; //FIXME this is probably too small as default
#if 0
                /* see VCD standard, p. IV-7*/
                stream->max_buffer_size = 46 * 1024;
            else
                /* This value HAS to be used for SVCD (see SVCD standard, p. 26 V.2.3.2).
                   Right now it is also used for everything else.*/
                stream->max_buffer_size = 230 * 1024;
#endif
            break;
        case CODEC_TYPE_SUBTITLE:
            stream->max_buffer_size = 16 * 1024;
            break;
        default:
            return -1;
        }
        av_fifo_init(&stream->fifo, 16);
    }
    return 0;
}

static inline void put_timestamp(uint8_t** p, int id, int64_t timestamp)
{
    bytestream_put_byte(p,
        (id << 4) |
        (((timestamp >> 30) & 0x07) << 1) |
        1);
    bytestream_put_be16(p, (uint16_t)((((timestamp >> 15) & 0x7fff) << 1) | 1));
    bytestream_put_be16(p, (uint16_t)((((timestamp) & 0x7fff) << 1) | 1));
}

static int get_nb_frames(AVFormatContext *ctx, StreamInfo *stream, int len){
    int nb_frames=0;
    PacketDesc *pkt_desc= stream->premux_packet;

    while(len>0){
        if(pkt_desc->size == pkt_desc->unwritten_size)
            nb_frames++;
        len -= pkt_desc->unwritten_size;
        pkt_desc= pkt_desc->next;
    }

    return nb_frames;
}

int ff_pes_write_buf(AVFormatContext *ctx, int stream_index, uint8_t *buf,
                     int64_t *pts, int64_t *dts,
                     int trailer_size, int *packet_size, int *pad_packet_bytes,
                     int *payload_size, int *stuffing_size)
{
    StreamInfo *stream = ctx->streams[stream_index]->priv_data;
    int startcode, i, header_len;
    int pes_flags = 0;
    uint8_t *p = buf;
    int nb_frames;

    /* packet header size */
    *packet_size -= 6;

    /* packet header */
    if (stream->format & PES_FMT_MPEG2) {
        header_len = 3;
        if (stream->format != PES_FMT_TS && stream->packet_number==0)
            header_len += 3; /* PES extension */
        header_len += 1; /* obligatory stuffing byte */
    } else {
        header_len = 0;
    }

    if (*pts != AV_NOPTS_VALUE) {
        if (*dts != *pts)
            header_len += 5 + 5;
        else
            header_len += 5;
    } else {
        if (!(stream->format & PES_FMT_MPEG2))
            header_len++;
    }

    *payload_size = *packet_size - header_len;
    if (stream->id < 0xc0) {
        startcode = PRIVATE_STREAM_1;
        *payload_size -= 1;
        if (stream->id >= 0x40) {
            *payload_size -= 3;
            if (stream->id >= 0xa0)
                *payload_size -= 3;
        }
    } else {
        startcode = 0x100 + stream->id;
    }

    *stuffing_size = *payload_size - av_fifo_size(&stream->fifo);

    // first byte does not fit -> reset pts/dts + stuffing
    if(*payload_size <= trailer_size && *pts != AV_NOPTS_VALUE){
        int timestamp_len=0;
        if(*dts != *pts)
            timestamp_len += 5;
        if(*pts != AV_NOPTS_VALUE)
            timestamp_len += stream->format & PES_FMT_MPEG2 ? 5 : 4;
        *pts=*dts= AV_NOPTS_VALUE;
        header_len -= timestamp_len;
        if (stream->format == PES_FMT_DVD && stream->align_iframe) {
            *pad_packet_bytes += timestamp_len;
            *packet_size -= timestamp_len;
        } else {
            *payload_size += timestamp_len;
        }
        *stuffing_size += timestamp_len;
        if(*payload_size > trailer_size)
            *stuffing_size += *payload_size - trailer_size;
    }

    if (stream->format != PES_FMT_TS) {
        if (*pad_packet_bytes > 0 && *pad_packet_bytes <= 7) { // can't use padding, so use stuffing
            *packet_size += *pad_packet_bytes;
            *payload_size += *pad_packet_bytes; // undo the previous adjustment
            if (*stuffing_size < 0) {
                *stuffing_size = *pad_packet_bytes;
            } else {
                *stuffing_size += *pad_packet_bytes;
            }
            *pad_packet_bytes = 0;
        }
    }

    if (*stuffing_size < 0)
        *stuffing_size = 0;
    if (*stuffing_size > 16) {    /*<=16 for MPEG-1, <=32 for MPEG-2*/
        *pad_packet_bytes += *stuffing_size;
        *packet_size -= *stuffing_size;
        *payload_size -= *stuffing_size;
        *stuffing_size = 0;
    }

    nb_frames= get_nb_frames(ctx, stream, *payload_size - *stuffing_size);

    bytestream_put_be32(&p, startcode);

    bytestream_put_be16(&p, *packet_size);

    if (!(stream->format & PES_FMT_MPEG2))
        for(i=0;i<*stuffing_size;i++)
            bytestream_put_byte(&p, 0xff);

    if (stream->format & PES_FMT_MPEG2) {
        bytestream_put_byte(&p, 0x80); /* mpeg2 id */

        pes_flags=0;

        if (*pts != AV_NOPTS_VALUE) {
            pes_flags |= 0x80;
            if (*dts != *pts)
                pes_flags |= 0x40;
        }

        /* Both the MPEG-2 and the SVCD standards demand that the
           P-STD_buffer_size field be included in the first packet of
           every stream. (see SVCD standard p. 26 V.2.3.1 and V.2.3.2
           and MPEG-2 standard 2.7.7) */
        if (stream->format != PES_FMT_TS && stream->packet_number == 0)
            pes_flags |= 0x01;

        bytestream_put_byte(&p, pes_flags); /* flags */
        bytestream_put_byte(&p, header_len - 3 + *stuffing_size);

        if (pes_flags & 0x80)  /*write pts*/
            put_timestamp(&p, (pes_flags & 0x40) ? 0x03 : 0x02, *pts);
        if (pes_flags & 0x40)  /*write dts*/
            put_timestamp(&p, 0x01, *dts);

        if (pes_flags & 0x01) {  /*write pes extension*/
            bytestream_put_byte(&p, 0x10); /* flags */

            /* P-STD buffer info */
            if (stream->id == AUDIO_ID)
                bytestream_put_be16(&p, 0x4000 | stream->max_buffer_size/128);
            else
                bytestream_put_be16(&p, 0x6000 | stream->max_buffer_size/1024);
        }

    } else {
        if (*pts != AV_NOPTS_VALUE) {
            if (*dts != *pts) {
                put_timestamp(&p, 0x03, *pts);
                put_timestamp(&p, 0x01, *dts);
            } else {
                put_timestamp(&p, 0x02, *pts);
            }
        } else {
            bytestream_put_byte(&p, 0x0f);
        }
    }

    if (stream->format & PES_FMT_MPEG2) {
        /* special stuffing byte that is always written
           to prevent accidental generation of start codes. */
        bytestream_put_byte(&p, 0xff);

        for(i=0;i<*stuffing_size;i++)
            bytestream_put_byte(&p, 0xff);
    }

    if (stream->format != PES_FMT_TS && startcode == PRIVATE_STREAM_1) {
        bytestream_put_byte(&p, stream->id);
        if (stream->id >= 0xa0) {
            /* LPCM (XXX: check nb_frames) */
            bytestream_put_byte(&p, 7);
            bytestream_put_be16(&p, 4); /* skip 3 header bytes */
            bytestream_put_byte(&p, stream->lpcm_header[0]);
            bytestream_put_byte(&p, stream->lpcm_header[1]);
            bytestream_put_byte(&p, stream->lpcm_header[2]);
        } else if (stream->id >= 0x40) {
            /* AC3 */
            bytestream_put_byte(&p, nb_frames);
            bytestream_put_be16(&p, trailer_size+1);
        }
    }

    /* output data */
    assert(*payload_size - *stuffing_size <= av_fifo_size(&stream->fifo));
    if(av_fifo_read(&stream->fifo, p, *payload_size - *stuffing_size) < 0)
        return -1;
    return p - buf + *payload_size - *stuffing_size;
}

int ff_pes_remove_decoded_packets(AVFormatContext *ctx, int64_t scr)
{
    int i;

    for(i=0; i<ctx->nb_streams; i++){
        AVStream *st = ctx->streams[i];
        StreamInfo *stream = st->priv_data;
        PacketDesc *pkt_desc;

        while((pkt_desc= stream->predecode_packet)
              && scr > pkt_desc->dts){ //FIXME > vs >=
            if(stream->buffer_index < pkt_desc->size ||
               stream->predecode_packet == stream->premux_packet){
                av_log(ctx, AV_LOG_ERROR,
                       "buffer underflow i=%d bufi=%d size=%d\n",
                       i, stream->buffer_index, pkt_desc->size);
                break;
            }
            stream->buffer_index -= pkt_desc->size;

            stream->predecode_packet= pkt_desc->next;
            av_freep(&pkt_desc);
        }
    }

    return 0;
}

void ff_pes_write_packet(AVFormatContext *ctx, AVPacket *pkt, int packet_number)
{
    int stream_index= pkt->stream_index;
    int size= pkt->size;
    uint8_t *buf= pkt->data;
    AVStream *st= ctx->streams[stream_index];
    StreamInfo *stream= st->priv_data;
    int64_t pts, dts;
    PacketDesc *pkt_desc;
    const int preload= av_rescale(ctx->preload, 90000, AV_TIME_BASE);
    const int is_iframe = st->codec->codec_type == CODEC_TYPE_VIDEO && (pkt->flags & PKT_FLAG_KEY);

    pts= pkt->pts;
    dts= pkt->dts;

    if(pts != AV_NOPTS_VALUE) pts += preload;
    if(dts != AV_NOPTS_VALUE) dts += preload;

//av_log(ctx, AV_LOG_DEBUG, "dts:%f pts:%f flags:%d stream:%d nopts:%d\n", dts/90000.0, pts/90000.0, pkt->flags, pkt->stream_index, pts != AV_NOPTS_VALUE);
    if (!stream->premux_packet)
        stream->next_packet = &stream->premux_packet;
    *stream->next_packet=
    pkt_desc= av_mallocz(sizeof(PacketDesc));
    pkt_desc->pts= pts;
    pkt_desc->dts= dts;
    pkt_desc->unwritten_size=
    pkt_desc->size= size;
    if(!stream->predecode_packet)
        stream->predecode_packet= pkt_desc;
    stream->next_packet= &pkt_desc->next;

    av_fifo_realloc(&stream->fifo, av_fifo_size(&stream->fifo) + size);

    if (stream->format == PES_FMT_DVD){
        if (is_iframe && (packet_number == 0 || (pts - stream->vobu_start_pts >= 36000))) { // min VOBU length 0.4 seconds (mpucoder)
            stream->bytes_to_iframe = av_fifo_size(&stream->fifo);
            stream->align_iframe = 1;
            stream->vobu_start_pts = pts;
        }
    }

    av_fifo_generic_write(&stream->fifo, buf, size, NULL);
}

int ff_pes_output_packet(AVFormatContext *ctx, int packet_size, int64_t *scr,
                         int *best_i, int flush, int (*flush_packet)())
{
    AVStream *st;
    StreamInfo *stream;
    int i, avail_space=0, trailer_size;
    int best_score= INT_MIN;
    int ignore_constraints=0;
    PacketDesc *timestamp_packet;
    const int64_t max_delay= av_rescale(ctx->max_delay, 90000, AV_TIME_BASE);

retry:
    for(i=0; i<ctx->nb_streams; i++){
        AVStream *st= ctx->streams[i];
        StreamInfo *stream= st->priv_data;
        const int avail_data= av_fifo_size(&stream->fifo);
        const int space= stream->max_buffer_size - stream->buffer_index;
        int rel_space= 1024*space / stream->max_buffer_size;
        PacketDesc *next_pkt= stream->premux_packet;

        /* for subtitle, a single PES packet must be generated,
           so we flush after every single subtitle packet */
        if(packet_size > avail_data && !flush
           && st->codec->codec_type != CODEC_TYPE_SUBTITLE)
            return 0;
        if(avail_data==0)
            continue;
        assert(avail_data>0);

        if(space < packet_size && !ignore_constraints)
            continue;

        if(next_pkt && next_pkt->dts - *scr > max_delay)
            continue;

        if(rel_space > best_score){
            best_score= rel_space;
            *best_i = i;
            avail_space= space;
        }
    }

    if(*best_i < 0){
        int64_t best_dts= INT64_MAX;

        for(i=0; i<ctx->nb_streams; i++){
            AVStream *st = ctx->streams[i];
            StreamInfo *stream = st->priv_data;
            PacketDesc *pkt_desc= stream->predecode_packet;
            if(pkt_desc && pkt_desc->dts < best_dts)
                best_dts= pkt_desc->dts;
        }

#if 0
        av_log(ctx, AV_LOG_DEBUG, "bumping scr, scr:%f, dts:%f\n",
               scr/90000.0, best_dts/90000.0);
#endif
        if(best_dts == INT64_MAX)
            return 0;

        if(*scr >= best_dts+1 && !ignore_constraints){
            av_log(ctx, AV_LOG_ERROR, "packet too large, ignoring buffer limits to mux it\n");
            ignore_constraints= 1;
        }
        *scr= FFMAX(best_dts+1, *scr);
        if(ff_pes_remove_decoded_packets(ctx, *scr) < 0)
            return -1;
        goto retry;
    }

    assert(*best_i >= 0);

    st = ctx->streams[*best_i];
    stream = st->priv_data;

    assert(av_fifo_size(&stream->fifo) > 0);

    assert(avail_space >= packet_size || ignore_constraints);

    timestamp_packet= stream->premux_packet;
    if(timestamp_packet->unwritten_size == timestamp_packet->size){
        trailer_size= 0;
    }else{
        trailer_size= timestamp_packet->unwritten_size;
        timestamp_packet= timestamp_packet->next;
    }

    if(timestamp_packet){
//av_log(ctx, AV_LOG_DEBUG, "dts:%f pts:%f pcr:%f stream:%d\n", timestamp_packet->dts/90000.0, timestamp_packet->pts/90000.0, pcr/90000.0, best_i);
        return flush_packet(ctx, *best_i, timestamp_packet->pts, timestamp_packet->dts, *scr, trailer_size);
    }else{
        assert(av_fifo_size(&stream->fifo) == trailer_size);
        return flush_packet(ctx, *best_i, AV_NOPTS_VALUE, AV_NOPTS_VALUE, *scr, trailer_size);
    }
}


void ff_pes_muxer_end(AVFormatContext *ctx)
{
    StreamInfo *stream;
    int i;

    for(i=0;i<ctx->nb_streams;i++) {
        stream = ctx->streams[i]->priv_data;

        assert(av_fifo_size(&stream->fifo) == 0);
        av_fifo_free(&stream->fifo);
    }
}
