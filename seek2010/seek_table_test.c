/*
 * Copyright (c) 2003 Fabrice Bellard
 * Copyright (c) 2007 Michael Niedermayer
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

/** @file
 *  Tests sample accuracy of seek by table
 */

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "libavutil/common.h"
#include "libavformat/avformat.h"

#undef exit
#undef printf
#undef fprintf

static char buffer[20];

static const char *ret_str(int v)
{
    switch (v) {
    case AVERROR_EOF:     return "-EOF";
    case AVERROR(EIO):    return "-EIO";
    case AVERROR(ENOMEM): return "-ENOMEM";
    case AVERROR(EINVAL): return "-EINVAL";
    default:
        snprintf(buffer, sizeof(buffer), "%2d", v);
        return buffer;
    }
}

static void ts_str(char buffer[60], int64_t ts, AVRational base)
{
    double tsval;
    if (ts == AV_NOPTS_VALUE) {
        strcpy(buffer, " NOPTS   ");
        return;
    }
    tsval = ts * av_q2d(base);
    snprintf(buffer, 60, "%9f", tsval);
}

int main(int argc, char **argv)
{
    const char *filename;
    AVFormatContext *ic = NULL;
    AVCodecContext *avctx;
    int i, ret, stream_id;
    int64_t timestamp;
    AVFormatParameters params, *ap= &params;
    int16_t *ref_buf;
    int16_t *test_buf;

    int offset = 0;
    int ref_len = 0;
    int data_size;
    AVStream* st;
    AVCodec *codec;
    int len;
    int channels;
    int done_size;
    int j;
    int64_t error;
    int mismatches;
    int64_t start_ts;
    int new_test;
    int done_samples;
    int offset_samples;
    int64_t start_pos;
#define REF_BUF_NB_SAMPLES (44100 * 120)

    ref_buf  = (int16_t*)malloc(REF_BUF_NB_SAMPLES * 2 * sizeof(*ref_buf));
    test_buf = (int16_t*)malloc(REF_BUF_NB_SAMPLES * 2 * sizeof(*test_buf));

    memset(ap, 0, sizeof(params));
    ap->channels    = 2;
    ap->sample_rate = 44100;

    av_log_set_level(99);
    /* initialize libavcodec, and register all codecs and formats */
    av_register_all();

    if (argc != 2) {
        printf("usage: %s input_file\n"
               "\n", argv[0]);
        exit(1);
    }

    filename = argv[1];

    ret = av_open_input_file(&ic, filename, NULL, 0, ap);
    if (ret < 0) {
        fprintf(stderr, "cannot open %s\n", filename);
        exit(1);
    }

    ret = av_find_stream_info(ic);
    if (ret < 0) {
        fprintf(stderr, "%s: could not find codec parameters\n", filename);
        exit(1);
    }

    for(i = 0; i < ic->nb_streams; i++) {
        avctx = ic->streams[i]->codec;
        if (avctx->codec_type == AVMEDIA_TYPE_AUDIO) {
            if (avctx->channels > 0) {
                channels                =
                avctx->request_channels = FFMIN(2, avctx->channels);
            } else {
                channels                =
                avctx->request_channels = 2;
            }
            st = ic->streams[i];
            start_ts = st->start_time != AV_NOPTS_VALUE ?
                       st->start_time :
                      (st->first_dts  != AV_NOPTS_VALUE ?
                       st->first_dts  : 0);

        }

        codec = avcodec_find_decoder(avctx->codec_id);
        if (!codec || avcodec_open(avctx, codec) < 0)
            return -1;
    }
    offset = 0;
    printf("generating ref\n");
    for(;;) {
        AVPacket pkt, pkt_tmp;
        AVStream* st;
        char ts_buf[60];
        memset(&pkt, 0, sizeof(pkt));
        if(ret >= 0){
            ret= av_read_frame(ic, &pkt);
            pkt_tmp = pkt;
            if(ret >= 0 && pkt.size > 0){
                char dts_buf[60];
                st= ic->streams[pkt.stream_index];

                while(pkt_tmp.size > 0) {
                    data_size = sizeof(*ref_buf) * (2 * REF_BUF_NB_SAMPLES - offset);
                    len = avcodec_decode_audio3(st->codec, ref_buf + offset / sizeof(*ref_buf),
                                                &data_size, &pkt);
                    if (len < 0) {
                        pkt_tmp.size = 0;
                        break;
                    }

                    pkt_tmp.data += len;
                    pkt_tmp.size -= len;
                    if (data_size <= 0)
                        continue;

                    printf("offset %d\n",offset);
                    offset += data_size;
                    ts_str(dts_buf, pkt.dts, st->time_base);
                    ts_str(ts_buf,  pkt.pts, st->time_base);
                    printf("ref ret:%-10s st:%2d flags:%d dts:%s pts:%s pos:%7" PRId64 " size:%6d", ret_str(ret), pkt.stream_index, pkt.flags, dts_buf, ts_buf, pkt.pos, pkt.size);

                    if(offset >= REF_BUF_NB_SAMPLES) {
                        if(pkt.data)
                            av_free_packet(&pkt);
                        goto after_ref;
                    }
                }
                if(pkt.data)
                    av_free_packet(&pkt);
            } else
                printf("ref ret:%s", ret_str(ret)); // necessary to avoid trailing whitespace
            printf("\n");
        }
    }
after_ref:
    ret     = 0;
    ref_len = offset;
    printf("testing seek\n");

    av_build_index(ic, 0);

    timestamp = 0;
    if (ret < 0)
        goto new_seek;

    for(i = 0; ; i++) {
        AVPacket pkt, pkt_tmp;
        char ts_buf[60];

        new_test  = 1;
        done_size = 0;
        if (ret < 0)
            goto new_seek;
        for(;;) {
            memset(&pkt, 0, sizeof(pkt));
            if(ret >= 0){
                start_pos = url_ftell(ic->pb);
                ret= av_read_frame(ic, &pkt);
                pkt_tmp = pkt;
                if(ret >= 0 && pkt.size > 0){
                    char dts_buf[60];
                    st= ic->streams[pkt.stream_index];

                    while(pkt_tmp.size > 0) {
                        if (new_test) {
                            new_test = 0;
                            offset   = ((pkt.dts - start_ts)   *
                                        avctx->sample_rate     *
                                        channels               *
                                        av_q2d(st->time_base)) * sizeof(*test_buf);
                            if(offset < 0 || offset >= REF_BUF_NB_SAMPLES) {
                                printf("ERROR: negative offset\n");
                                goto new_seek;
                            }
                            printf("testing offset %d, ts %f (pts %f), start_ts %f,"
                                   "pos %lld\n",
                                   offset,
                                   pkt.dts  * av_q2d(st->time_base),
                                   pkt.pts  * av_q2d(st->time_base),
                                   start_ts * av_q2d(st->time_base),
                                   start_pos);
                            memset(test_buf, 0xFF, sizeof(*test_buf) * 2 * REF_BUF_NB_SAMPLES);
                        }

                        data_size = sizeof(*test_buf) * (2 * REF_BUF_NB_SAMPLES - offset);
                        len = avcodec_decode_audio3(st->codec,
                                                    test_buf + offset / sizeof(*test_buf),
                                                    &data_size, &pkt);
                        if (len < 0) {
                            pkt_tmp.size = 0;
                            break;
                        }

                        pkt_tmp.data += len;
                        pkt_tmp.size -= len;
                        if (data_size <= 0)
                            continue;

                        offset    += data_size;
                        done_size += data_size;
                        ts_str(dts_buf, pkt.dts, st->time_base);
                        ts_str(ts_buf,  pkt.pts, st->time_base);
                        printf("ref ret:%-10s st:%2d flags:%d dts:%s pts:%s pos:%7"
                               PRId64 " size:%6d",
                               ret_str(ret), pkt.stream_index, pkt.flags,
                               dts_buf, ts_buf, pkt.pos, pkt.size);

#define TEST_BUF_READ_MAX (44100 * 1)
                        if(done_size >= TEST_BUF_READ_MAX) {
                            if(pkt.data)
                                av_free_packet(&pkt);
                            goto single_test_done;
                        }
                    }
                    if(pkt.data)
                        av_free_packet(&pkt);
                } else
                    printf("ref ret:%s", ret_str(ret)); // necessary to avoid trailing whitespace
                printf("\n");
            }
        }
    single_test_done:
        error          = 0;
        mismatches     = 0;
        done_samples   = done_size / sizeof(*test_buf);
        offset_samples = offset    / sizeof(*test_buf);
        printf("compare start %d\n", offset - done_size);
        for (j = 0; j < done_samples; j++) {
            error += (int) fabs(test_buf[offset_samples - done_samples + j] -
                                ref_buf [offset_samples - done_samples + j]);
            if (test_buf[offset_samples - done_samples + j] -
                ref_buf [offset_samples - done_samples + j]) {
                mismatches++;
            }
            if(error)
                printf("%8d %8d  %d  %d\n", test_buf[offset_samples - done_samples + j],
                       ref_buf [offset_samples - done_samples + j], i, j);
        }

        if(error)
            printf("error of %lld found with %d of %d samples mismatching\n",
                   error, mismatches, done_samples);
        else
            printf("No errors\n");

    new_seek:
        if(i>25) break;

        stream_id = (i >> 1) % (ic->nb_streams + 1) - 1;
        timestamp = (i * 19362894167LL) % (20 * AV_TIME_BASE) - AV_TIME_BASE;
        if(stream_id>=0){
            st= ic->streams[stream_id];
            timestamp= av_rescale_q(timestamp, AV_TIME_BASE_Q, st->time_base);
        }
        //FIXME fully test the new seek API
        if(i&1) ret = avformat_seek_file(ic, stream_id, INT64_MIN, timestamp, timestamp, 0);
        else    ret = avformat_seek_file(ic, stream_id, timestamp, timestamp, INT64_MAX, 0);
        ts_str(ts_buf, timestamp, stream_id < 0 ? AV_TIME_BASE_Q : st->time_base);
        printf("ret:%-10s st:%2d flags:%d  ts:%s\n", ret_str(ret), stream_id, i & 1, ts_buf);
    }

    av_close_input_file(ic);

    return 0;
}
