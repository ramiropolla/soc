/*
 * MPEG PES muxer
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

#include "mpegpes.h"
#include "bytestream.h"

/**
 * Initialization of PES muxer.
 * @param[in] ctx the AVFormatContext which contains streams
 * @return  On error a negative value is returned, on success zero.
 */
int ff_pes_muxer_init(AVFormatContext *ctx)
{
    AVStream *st;
    PESStream *stream;
    int i;

    for(i=0;i<ctx->nb_streams;i++) {
        st = ctx->streams[i];
        stream = st->priv_data;
        av_set_pts_info(st, 64, 1, 90000);

        switch(st->codec->codec_type) {
        case CODEC_TYPE_AUDIO:
            /* This value HAS to be used for VCD (see VCD standard, p. IV-7).
               Right now it is also used for everything else.*/
            stream->max_buffer_size = 4 * 1024;
            break;
        case CODEC_TYPE_VIDEO:
#if 0
                /* see VCD standard, p. IV-7*/
                stream->max_buffer_size = 46 * 1024;
            else
                /* This value HAS to be used for SVCD (see SVCD standard, p. 26 V.2.3.2).
                   Right now it is also used for everything else.*/
                stream->max_buffer_size = 230 * 1024;
#endif
            if (st->codec->rc_buffer_size)
                stream->max_buffer_size = 6*1024 + st->codec->rc_buffer_size/8;
            else
                stream->max_buffer_size = 230*1024; //FIXME this is probably too small as default
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

static inline void insert_timestamp(uint8_t** p, int id, int64_t timestamp)
{
    bytestream_put_byte(p,
        (id << 4) |
        (((timestamp >> 30) & 0x07) << 1) |
        1);
    bytestream_put_be16(p, (uint16_t)((((timestamp >> 15) & 0x7fff) << 1) | 1));
    bytestream_put_be16(p, (uint16_t)((((timestamp) & 0x7fff) << 1) | 1));
}

/**
 * Get total number of frames that have been muxed.
 * @param[in] ctx    the AVFormatContext
 * @param[in] stream the PES stream
 * @param[in] len    PES packet size
 * @return  the number of frames have been muxed.
 */
int ff_pes_get_nb_frames(AVFormatContext *ctx, PESStream *stream, int len){
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
int ff_pes_muxer_write(AVFormatContext *ctx, int stream_index, uint8_t* pes_buffer,
    int64_t pts,int64_t dts, int  id, int startcode,
    uint8_t* pes_content, int pes_content_len,
    int header_len, int packet_size, int payload_size, int stuffing_size)
{
    PESStream *stream = ctx->streams[stream_index]->priv_data;
    PESContext *context = ctx->priv_data;
    int pes_flags, i;
    int data_size = payload_size - stuffing_size;
    uint8_t *q = pes_buffer;

    bytestream_put_be32(&q, startcode);

    bytestream_put_be16(&q, packet_size);
    bytestream_put_byte(&q, 0x80); /* mpeg2 id */

    pes_flags=0;

    if (pts != AV_NOPTS_VALUE) {
        pes_flags |= 0x80;
        if (dts != pts)
           pes_flags |= 0x40;
    }

    /* Both the MPEG-2 and the SVCD standards demand that the
       P-STD_buffer_size field be included in the first packet of
       every stream. (see SVCD standard p. 26 V.2.3.1 and V.2.3.2
       and MPEG-2 standard 2.7.7) */
    if (context->packet_number == 0 && context->muxer_type == PESMUXER_PS)
        pes_flags |= 0x01;

    bytestream_put_byte(&q, pes_flags); /* flags */
    bytestream_put_byte(&q, header_len - 3 + stuffing_size);

    if (pes_flags & 0x80)  /*write pts*/
        insert_timestamp(&q, (pes_flags & 0x40) ? 0x03 : 0x02, pts);
    if (pes_flags & 0x40)  /*write dts*/
        insert_timestamp(&q, 0x01, dts);

    if (pes_flags & 0x01) {  /*write pes extension*/
        bytestream_put_byte(&q, 0x10); /* flags */

        /* P-STD buffer info */
        if (id == AUDIO_ID)
            bytestream_put_be16(&q, 0x4000 | stream->max_buffer_size/128);
        else
            bytestream_put_be16(&q, 0x6000 | stream->max_buffer_size/1024);
    }

    /* special stuffing byte that is always written
       to prevent accidental generation of startcodes. */
    bytestream_put_byte(&q, 0xff);

    for(i=0;i<stuffing_size;i++)
        bytestream_put_byte(&q, 0xff);

    if(pes_content != NULL)
        bytestream_put_buffer(&q, pes_content, pes_content_len);

    /* output data */
    if(av_fifo_read(&stream->fifo, q, data_size) < 0)
        return -1;
    return (q - pes_buffer + data_size);
}

/**
 * Remove decoded packets of each stream.
 * @param[in] ctx  the AVFormatContext
 * @param[in] scr  System Clock Reference of PES stream
 * @return  On error a negative or zero value is returned, on success 1 is returned.
 */
int ff_pes_remove_decoded_packets(AVFormatContext *ctx, int64_t scr)
{
    int i;

    for(i=0; i<ctx->nb_streams; i++){
        AVStream *st = ctx->streams[i];
        PESStream  *stream = st->priv_data;
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

/**
 * Find the stream to mux into the PES stream.
 * @param[in] ctx          the AVFormatContext
 * @param[in] packet_size  PES stream packet size
 * @param[in] flush        Flush after every single subtitle packet.
 * @param[out] best_i      index of stream to be muxed
 * @return  On error a negative or zero value is returned, on success 1 is returned.
 */
int ff_pes_find_beststream(AVFormatContext *ctx, int packet_size, int flush, int64_t *scr, int* best_i)
{
    int i, avail_space;
    int best_score= INT_MIN;
    int ignore_constraints=0;
    const int64_t max_delay= av_rescale(ctx->max_delay, 90000, AV_TIME_BASE);

retry:
    for(i=0; i<ctx->nb_streams; i++){
        AVStream *st = ctx->streams[i];
        PESStream*stream = st->priv_data;
        const int avail_data=  av_fifo_size(&stream->fifo);
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
            PESStream *stream = st->priv_data;
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
    assert(avail_space >= packet_size || ignore_constraints);
    return 1;
}

/**
 * Write packet into PES FIFO.
 * @param [in] ctx  the AVFormatContext which contains streams.
 * @param [in] pkt  the packet to write.
 * @return  NULL
 */
void ff_pes_write_packet(AVFormatContext *ctx, AVPacket *pkt)
{
    int stream_index= pkt->stream_index;
    int size= pkt->size;
    uint8_t *buf= pkt->data;
    AVStream *st = ctx->streams[stream_index];
    PESStream *stream = st->priv_data;
    int64_t pts, dts;
    PacketDesc *pkt_desc;
    const int preload= av_rescale(ctx->preload, 90000, AV_TIME_BASE);

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

    av_fifo_realloc(&stream->fifo, av_fifo_size(&stream->fifo) + size + 1);
    av_fifo_write(&stream->fifo, buf, size);
}

/**
 * Finalization of PES muxer.
 * @param [in] ctx the AVFormatContext which contains streams.
 * @return  NULL
 */
void ff_pes_muxer_end(AVFormatContext *ctx)
{
    PESStream *stream;
    int i;

    for(i=0;i<ctx->nb_streams;i++) {
        stream = ctx->streams[i]->priv_data;

        assert(av_fifo_size(&stream->fifo) == 0);
        av_fifo_free(&stream->fifo);
    }
}
