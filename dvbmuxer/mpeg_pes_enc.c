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

#include "mpeg_pes.h"

int ff_pes_muxer_init(AVFormatContext *ctx)
{
    AVStream *st;
    PESStream *stream;
    int i;

    for(i=0;i<ctx->nb_streams;i++) {
        st = ctx->streams[i];
        stream = (PESStream*)st->priv_data;
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

void ff_put_timestamp(ByteIOContext *pb, int id, int64_t timestamp)
{
    put_byte(pb,
             (id << 4) |
             (((timestamp >> 30) & 0x07) << 1) |
             1);
    put_be16(pb, (uint16_t)((((timestamp >> 15) & 0x7fff) << 1) | 1));
    put_be16(pb, (uint16_t)((((timestamp) & 0x7fff) << 1) | 1));
}

int ff_get_nb_frames(AVFormatContext *ctx, PESStream *stream, int len){
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

int ff_pes_muxer_write(AVFormatContext *ctx, int stream_index,
    int64_t pts,int64_t dts, int  id, int startcode,
    uint8_t* pes_content, int pes_content_len,
    int header_len, int packet_size, int payload_size, int stuffing_size)
{
    PESStream *stream = ctx->streams[stream_index]->priv_data;
    PESContext *context = ctx->priv_data;
    int pes_flags, i;
    int data_size = payload_size - stuffing_size;

        put_be32(&ctx->pb, startcode);

        put_be16(&ctx->pb, packet_size);
    put_byte(&ctx->pb, 0x80); /* mpeg2 id */

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

            put_byte(&ctx->pb, pes_flags); /* flags */
            put_byte(&ctx->pb, header_len - 3 + stuffing_size);

            if (pes_flags & 0x80)  /*write pts*/
                ff_put_timestamp(&ctx->pb, (pes_flags & 0x40) ? 0x03 : 0x02, pts);
            if (pes_flags & 0x40)  /*write dts*/
                ff_put_timestamp(&ctx->pb, 0x01, dts);

            if (pes_flags & 0x01) {  /*write pes extension*/
                put_byte(&ctx->pb, 0x10); /* flags */

                /* P-STD buffer info */
                if (id == AUDIO_ID)
                    put_be16(&ctx->pb, 0x4000 | stream->max_buffer_size/128);
                else
                    put_be16(&ctx->pb, 0x6000 | stream->max_buffer_size/1024);
    }

            /* special stuffing byte that is always written
               to prevent accidental generation of startcodes. */
             put_byte(&ctx->pb, 0xff);

             for(i=0;i<stuffing_size;i++)
                 put_byte(&ctx->pb, 0xff);

             put_buffer(&ctx->pb, pes_content, pes_content_len);

    /* output data */
    if(av_fifo_generic_read(&stream->fifo, data_size, &put_buffer, &ctx->pb) < 0)
        return -1;
    return data_size;
}

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


int ff_pes_find_beststream(AVFormatContext *ctx, int packet_size, int flush, int64_t scr, int* best_i)
{
    int best_score = INT_MIN;
    int i, avail_space = 0;
    int ignore_constraints = 0;
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

        if(next_pkt && next_pkt->dts - scr > max_delay)
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

        if(scr >= best_dts+1 && !ignore_constraints){
            av_log(ctx, AV_LOG_ERROR, "packet too large, ignoring buffer limits to mux it\n");
            ignore_constraints= 1;
        }
        scr= FFMAX(best_dts+1, scr);
        if(ff_pes_remove_decoded_packets(ctx, scr) < 0)
            return -1;
        goto retry;
    }
    assert(avail_space >= packet_size || ignore_constraints);
    return 1;
}


void ff_pes_write_packet(AVFormatContext *ctx, AVPacket *pkt)
{
    int stream_index= pkt->stream_index;
    AVStream *st = ctx->streams[stream_index];
    PESStream *stream = st->priv_data;
    PacketDesc *pkt_desc;
    int size= pkt->size;
    uint8_t *buf= pkt->data;
    int64_t pts, dts;
    const int preload= av_rescale(ctx->preload, 90000, AV_TIME_BASE);

    pts= pkt->pts;
    dts= pkt->dts;

    if(pts != AV_NOPTS_VALUE) pts += preload;
    if(dts != AV_NOPTS_VALUE) dts += preload;

    if (!stream->premux_packet)
        stream->next_packet = &stream->premux_packet;
    *stream->next_packet=
    pkt_desc= av_mallocz(sizeof(PacketDesc));
    pkt_desc->pts= pkt->pts;
    pkt_desc->dts= pkt->dts;
    pkt_desc->unwritten_size=
    pkt_desc->size= size;
    if(!stream->predecode_packet)
        stream->predecode_packet= pkt_desc;
    stream->next_packet= &pkt_desc->next;

    av_fifo_realloc(&stream->fifo, av_fifo_size(&stream->fifo) + size + 1);
    av_fifo_write(&stream->fifo, buf, size);
}


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
