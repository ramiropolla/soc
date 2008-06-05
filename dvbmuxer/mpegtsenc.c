/*
 * MPEG2 transport stream (aka DVB) muxer
 * Copyright (c) 2003 Fabrice Bellard.
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

#include "libavutil/crc.h"
#include "avformat.h"
#include "mpegts.h"
#include "libavcodec/bytestream.h"
#include "mpegpes.h"

/* write DVB SI sections */

/*********************************************/
/* mpegts section writer */

typedef struct MpegTSSection {
    int pid;
    int cc;
    void (*write_packet)(struct MpegTSSection *s, const uint8_t *packet);
    void *opaque;
} MpegTSSection;

/* NOTE: 4 bytes must be left at the end for the crc32 */
static void mpegts_write_section(MpegTSSection *s, uint8_t *buf, int len)
{
    unsigned int crc;
    unsigned char packet[TS_PACKET_SIZE];
    const unsigned char *buf_ptr;
    unsigned char *q;
    int first, b, len1, left;

    crc = bswap_32(av_crc(av_crc_get_table(AV_CRC_32_IEEE), -1, buf, len - 4));
    buf[len - 4] = (crc >> 24) & 0xff;
    buf[len - 3] = (crc >> 16) & 0xff;
    buf[len - 2] = (crc >> 8) & 0xff;
    buf[len - 1] = (crc) & 0xff;

    /* send each packet */
    buf_ptr = buf;
    while (len > 0) {
        first = (buf == buf_ptr);
        q = packet;
        *q++ = 0x47;
        b = (s->pid >> 8);
        if (first)
            b |= 0x40;
        *q++ = b;
        *q++ = s->pid;
        *q++ = 0x10 | s->cc;
        s->cc = (s->cc + 1) & 0xf;
        if (first)
            *q++ = 0; /* 0 offset */
        len1 = TS_PACKET_SIZE - (q - packet);
        if (len1 > len)
            len1 = len;
        memcpy(q, buf_ptr, len1);
        q += len1;
        /* add known padding data */
        left = TS_PACKET_SIZE - (q - packet);
        if (left > 0)
            memset(q, 0xff, left);

        s->write_packet(s, packet);

        buf_ptr += len1;
        len -= len1;
    }
}

static inline void put16(uint8_t **q_ptr, int val)
{
    uint8_t *q;
    q = *q_ptr;
    *q++ = val >> 8;
    *q++ = val;
    *q_ptr = q;
}

static int mpegts_write_section1(MpegTSSection *s, int tid, int id,
                          int version, int sec_num, int last_sec_num,
                          uint8_t *buf, int len)
{
    uint8_t section[1024], *q;
    unsigned int tot_len;

    tot_len = 3 + 5 + len + 4;
    /* check if not too big */
    if (tot_len > 1024)
        return -1;

    q = section;
    *q++ = tid;
    put16(&q, 0xb000 | (len + 5 + 4)); /* 5 byte header + 4 byte CRC */
    put16(&q, id);
    *q++ = 0xc1 | (version << 1); /* current_next_indicator = 1 */
    *q++ = sec_num;
    *q++ = last_sec_num;
    memcpy(q, buf, len);

    mpegts_write_section(s, section, tot_len);
    return 0;
}

/*********************************************/
/* mpegts writer */

#define DEFAULT_PMT_START_PID   0x1000
#define DEFAULT_START_PID       0x0100
#define DEFAULT_PROVIDER_NAME   "FFmpeg"
#define DEFAULT_SERVICE_NAME    "Service01"

/* default network id, transport stream and service identifiers */
#define DEFAULT_ONID            0x0001
#define DEFAULT_TSID            0x0001
#define DEFAULT_SID             0x0001

/* a PES packet header is generated every DEFAULT_PES_HEADER_FREQ packets */
#define DEFAULT_PES_HEADER_FREQ 16
#define DEFAULT_PES_PAYLOAD_SIZE ((DEFAULT_PES_HEADER_FREQ - 1) * 184 + 170)

/* we retransmit the SI info at this rate */
#define SDT_RETRANS_TIME 500
#define PAT_RETRANS_TIME 100
#define MAX_DELTA_PCR 9000 /**< 0.1s according to ISO 13818-1 */

typedef struct MpegTSWriteStream {
    StreamInfo pes;
    struct MpegTSService *service;
    int pid; /* stream associated pid */
    int cc;
    int payload_index;
    int64_t payload_pts;
    int64_t payload_dts;
    uint8_t payload[DEFAULT_PES_PAYLOAD_SIZE];
} MpegTSWriteStream;

typedef struct MpegTSService {
    MpegTSSection pmt; /* MPEG2 pmt table context */
    int sid;           /* service ID */
    char *name;
    char *provider_name;
    int pcr_pid;
} MpegTSService;

typedef struct MpegTSWrite {
    MpegTSSection pat; /* MPEG2 pat table */
    MpegTSSection sdt; /* MPEG2 sdt table context */
    MpegTSService **services;
    int sdt_packet_count;
    int sdt_packet_freq;
    int pat_packet_count;
    int pat_packet_freq;
    int nb_services;
    int onid;
    int tsid;
    int64_t last_pcr; ///< last program clock reference */
    int64_t cur_pcr;  ///< last program clock reference */
    int mux_rate;
} MpegTSWrite;

static void mpegts_write_pat(AVFormatContext *s)
{
    MpegTSWrite *ts = s->priv_data;
    MpegTSService *service;
    uint8_t data[1012], *q;
    int i;

    q = data;
    for(i = 0; i < ts->nb_services; i++) {
        service = ts->services[i];
        put16(&q, service->sid);
        put16(&q, 0xe000 | service->pmt.pid);
    }
    mpegts_write_section1(&ts->pat, PAT_TID, ts->tsid, 0, 0, 0,
                          data, q - data);
}

static void mpegts_write_pmt(AVFormatContext *s, MpegTSService *service)
{
    //    MpegTSWrite *ts = s->priv_data;
    uint8_t data[1012], *q, *desc_length_ptr, *program_info_length_ptr;
    int val, stream_type, i;

    q = data;
    put16(&q, 0xe000 | service->pcr_pid);

    program_info_length_ptr = q;
    q += 2; /* patched after */

    /* put program info here */

    val = 0xf000 | (q - program_info_length_ptr - 2);
    program_info_length_ptr[0] = val >> 8;
    program_info_length_ptr[1] = val;

    for(i = 0; i < s->nb_streams; i++) {
        AVStream *st = s->streams[i];
        MpegTSWriteStream *ts_st = st->priv_data;
        switch(st->codec->codec_id) {
        case CODEC_ID_MPEG1VIDEO:
        case CODEC_ID_MPEG2VIDEO:
            stream_type = STREAM_TYPE_VIDEO_MPEG2;
            break;
        case CODEC_ID_MPEG4:
            stream_type = STREAM_TYPE_VIDEO_MPEG4;
            break;
        case CODEC_ID_H264:
            stream_type = STREAM_TYPE_VIDEO_H264;
            break;
        case CODEC_ID_MP2:
        case CODEC_ID_MP3:
            stream_type = STREAM_TYPE_AUDIO_MPEG1;
            break;
        case CODEC_ID_AAC:
            stream_type = STREAM_TYPE_AUDIO_AAC;
            break;
        case CODEC_ID_AC3:
            stream_type = STREAM_TYPE_AUDIO_AC3;
            break;
        default:
            stream_type = STREAM_TYPE_PRIVATE_DATA;
            break;
        }
        *q++ = stream_type;
        put16(&q, 0xe000 | ts_st->pid);
        desc_length_ptr = q;
        q += 2; /* patched after */

        /* write optional descriptors here */
        switch(st->codec->codec_type) {
        case CODEC_TYPE_AUDIO:
            if (strlen(st->language) == 3) {
                *q++ = 0x0a; /* ISO 639 language descriptor */
                *q++ = 4;
                *q++ = st->language[0];
                *q++ = st->language[1];
                *q++ = st->language[2];
                *q++ = 0; /* undefined type */
            }
            break;
        case CODEC_TYPE_SUBTITLE:
            {
                const char *language;
                language = st->language;
                if (strlen(language) != 3)
                    language = "eng";
                *q++ = 0x59;
                *q++ = 8;
                *q++ = language[0];
                *q++ = language[1];
                *q++ = language[2];
                *q++ = 0x10; /* normal subtitles (0x20 = if hearing pb) */
                put16(&q, 1); /* page id */
                put16(&q, 1); /* ancillary page id */
            }
            break;
        }

        val = 0xf000 | (q - desc_length_ptr - 2);
        desc_length_ptr[0] = val >> 8;
        desc_length_ptr[1] = val;
    }
    mpegts_write_section1(&service->pmt, PMT_TID, service->sid, 0, 0, 0,
                          data, q - data);
}

/* NOTE: str == NULL is accepted for an empty string */
static void putstr8(uint8_t **q_ptr, const char *str)
{
    uint8_t *q;
    int len;

    q = *q_ptr;
    if (!str)
        len = 0;
    else
        len = strlen(str);
    *q++ = len;
    memcpy(q, str, len);
    q += len;
    *q_ptr = q;
}

static void mpegts_write_sdt(AVFormatContext *s)
{
    MpegTSWrite *ts = s->priv_data;
    MpegTSService *service;
    uint8_t data[1012], *q, *desc_list_len_ptr, *desc_len_ptr;
    int i, running_status, free_ca_mode, val;

    q = data;
    put16(&q, ts->onid);
    *q++ = 0xff;
    for(i = 0; i < ts->nb_services; i++) {
        service = ts->services[i];
        put16(&q, service->sid);
        *q++ = 0xfc | 0x00; /* currently no EIT info */
        desc_list_len_ptr = q;
        q += 2;
        running_status = 4; /* running */
        free_ca_mode = 0;

        /* write only one descriptor for the service name and provider */
        *q++ = 0x48;
        desc_len_ptr = q;
        q++;
        *q++ = 0x01; /* digital television service */
        putstr8(&q, service->provider_name);
        putstr8(&q, service->name);
        desc_len_ptr[0] = q - desc_len_ptr - 1;

        /* fill descriptor length */
        val = (running_status << 13) | (free_ca_mode << 12) |
            (q - desc_list_len_ptr - 2);
        desc_list_len_ptr[0] = val >> 8;
        desc_list_len_ptr[1] = val;
    }
    mpegts_write_section1(&ts->sdt, SDT_TID, ts->tsid, 0, 0, 0,
                          data, q - data);
}

static MpegTSService *mpegts_add_service(MpegTSWrite *ts,
                                         int sid,
                                         const char *provider_name,
                                         const char *name)
{
    MpegTSService *service;

    service = av_mallocz(sizeof(MpegTSService));
    if (!service)
        return NULL;
    service->pmt.pid = DEFAULT_PMT_START_PID + ts->nb_services - 1;
    service->sid = sid;
    service->provider_name = av_strdup(provider_name);
    service->name = av_strdup(name);
    service->pcr_pid = 0x1fff;
    dynarray_add(&ts->services, &ts->nb_services, service);
    return service;
}

static void section_write_packet(MpegTSSection *s, const uint8_t *packet)
{
    AVFormatContext *ctx = s->opaque;
    put_buffer(ctx->pb, packet, TS_PACKET_SIZE);
}

static int mpegts_write_header(AVFormatContext *s)
{
    MpegTSWrite *ts = s->priv_data;
    MpegTSWriteStream *ts_st;
    MpegTSService *service;
    AVStream *st;
    int i, total_bit_rate;
    const char *service_name;

    ts->tsid = DEFAULT_TSID;
    ts->onid = DEFAULT_ONID;
    /* allocate a single DVB service */
    service_name = s->title;
    if (service_name[0] == '\0')
        service_name = DEFAULT_SERVICE_NAME;
    service = mpegts_add_service(ts, DEFAULT_SID,
                                 DEFAULT_PROVIDER_NAME, service_name);
    service->pmt.write_packet = section_write_packet;
    service->pmt.opaque = s;

    ts->pat.pid = PAT_PID;
    ts->pat.cc = 0;
    ts->pat.write_packet = section_write_packet;
    ts->pat.opaque = s;

    ts->sdt.pid = SDT_PID;
    ts->sdt.cc = 0;
    ts->sdt.write_packet = section_write_packet;
    ts->sdt.opaque = s;

    /* assign pids to each stream */
    total_bit_rate = 0;
    for(i = 0;i < s->nb_streams; i++) {
        int codec_rate;
        st = s->streams[i];
        ts_st = av_mallocz(sizeof(MpegTSWriteStream));
        if (!ts_st)
            goto fail;
        st->priv_data = ts_st;
        ts_st->service = service;
        ts_st->pid = DEFAULT_START_PID + i;
        ts_st->payload_pts = AV_NOPTS_VALUE;
        ts_st->payload_dts = AV_NOPTS_VALUE;
        /* update PCR pid by using the first video stream */
        if (st->codec->codec_type == CODEC_TYPE_VIDEO &&
            service->pcr_pid == 0x1fff)
            service->pcr_pid = ts_st->pid;

        if (st->codec->codec_type == CODEC_TYPE_VIDEO) {
            ts_st->pes.id = 0xe0;
        } else if (st->codec->codec_type == CODEC_TYPE_AUDIO &&
                   (st->codec->codec_id == CODEC_ID_MP2 ||
                    st->codec->codec_id == CODEC_ID_MP3)) {
            ts_st->pes.id = 0xc0;
        } else {
            ts_st->pes.id = 0xbd;
            if (st->codec->codec_type == CODEC_TYPE_SUBTITLE) {
                //private_code = 0x20;
            }
        }

        ts_st->pes.format = PES_FMT_TS;

        if(st->codec->rc_max_rate)
            codec_rate= st->codec->rc_max_rate;
        else
            codec_rate= st->codec->bit_rate;
        if(!codec_rate)
            codec_rate= (1<<21)*8/s->nb_streams;
        total_bit_rate += codec_rate;
    }

    /* if no video stream, use the first stream as PCR */
    if (service->pcr_pid == 0x1fff && s->nb_streams > 0) {
        ts_st = s->streams[0]->priv_data;
        service->pcr_pid = ts_st->pid;
    }

    if(ff_pes_muxer_init(s) != 0)
        goto fail;

    if(s->mux_rate) {
        ts->mux_rate= s->mux_rate;
    } else {
        total_bit_rate +=
            total_bit_rate * 25 / (8 * DEFAULT_PES_PAYLOAD_SIZE) + /* PES header size */
            total_bit_rate * 4 / (8 * TS_PACKET_SIZE) +            /* TS  header size */
            SDT_RETRANS_TIME * 12 +                                /* SDT size */
            PAT_RETRANS_TIME * 16;                                 /* PAT size */
        ts->mux_rate = total_bit_rate;
    }

    ts->sdt_packet_freq = (total_bit_rate * SDT_RETRANS_TIME) /
        (TS_PACKET_SIZE * 8 * 1000);
    ts->pat_packet_freq = (total_bit_rate * PAT_RETRANS_TIME) /
        (TS_PACKET_SIZE * 8 * 1000);
#if 0
    printf("%d %d %d\n",
           total_bit_rate, ts->sdt_packet_freq, ts->pat_packet_freq);
#endif

    /* write info at the start of the file, so that it will be fast to
       find them */
    mpegts_write_sdt(s);
    mpegts_write_pat(s);
    for(i = 0; i < ts->nb_services; i++) {
        mpegts_write_pmt(s, ts->services[i]);
    }
    put_flush_packet(s->pb);

    return 0;

 fail:
    for(i = 0;i < s->nb_streams; i++) {
        st = s->streams[i];
        av_free(st->priv_data);
    }
    return -1;
}

/* send SDT, PAT and PMT tables regulary */
static void retransmit_si_info(AVFormatContext *s)
{
    MpegTSWrite *ts = s->priv_data;
    int i;

    if (++ts->sdt_packet_count == ts->sdt_packet_freq) {
        ts->sdt_packet_count = 0;
        mpegts_write_sdt(s);
    }
    if (++ts->pat_packet_count == ts->pat_packet_freq) {
        ts->pat_packet_count = 0;
        mpegts_write_pat(s);
        for(i = 0; i < ts->nb_services; i++) {
            mpegts_write_pmt(s, ts->services[i]);
        }
    }
}

/* NOTE: pes_data contains all the PES packet */
static void mpegts_write_pes(AVFormatContext *s, AVStream *st, int64_t pcr,
                             const uint8_t *payload, int payload_size)
{
    MpegTSWriteStream *ts_st = st->priv_data;
    MpegTSWrite *ts = s->priv_data;
    uint8_t buf[TS_PACKET_SIZE];
    uint8_t *q;
    int val, is_start, len, header_len, write_pcr;
    int afc_len, stuffing_len;

    is_start = 1;
    ts->cur_pcr = pcr;
    while (payload_size > 0) {
        retransmit_si_info(s);

        write_pcr = !ts->cur_pcr;
        if (ts_st->pid == ts_st->service->pcr_pid) {
            pcr = ts->cur_pcr + (TS_PACKET_SIZE+4+7)*8*90000LL / ts->mux_rate;
            if (pcr - ts->last_pcr > MAX_DELTA_PCR)
                write_pcr = 1;
        }

        /* prepare packet header */
        q = buf;
        *q++ = 0x47;
        val = (ts_st->pid >> 8);
        if (is_start) {
            val |= 0x40;
            is_start = 0;
        }
        *q++ = val;
        *q++ = ts_st->pid;
        *q++ = 0x10 | ts_st->cc | (write_pcr ? 0x20 : 0);
        ts_st->cc = (ts_st->cc + 1) & 0xf;
        if (write_pcr) {
            /* add header and pcr bytes to pcr according to specs */
            pcr = ts->cur_pcr + (4+7)*8*90000LL / ts->mux_rate;
            *q++ = 7; /* AFC length */
            *q++ = 0x10; /* flags: PCR present */
            *q++ = pcr >> 25;
            *q++ = pcr >> 17;
            *q++ = pcr >> 9;
            *q++ = pcr >> 1;
            *q++ = (pcr & 1) << 7;
            *q++ = 0;
            ts->last_pcr = pcr;
        }
        /* header size */
        header_len = q - buf;
        /* data len */
        len = TS_PACKET_SIZE - header_len;
        if (len > payload_size)
            len = payload_size;
        stuffing_len = TS_PACKET_SIZE - header_len - len;
        if (stuffing_len > 0) {
            /* add stuffing with AFC */
            if (buf[3] & 0x20) {
                /* stuffing already present: increase its size */
                afc_len = buf[4] + 1;
                memmove(buf + 4 + afc_len + stuffing_len,
                        buf + 4 + afc_len,
                        header_len - (4 + afc_len));
                buf[4] += stuffing_len;
                memset(buf + 4 + afc_len, 0xff, stuffing_len);
            } else {
                /* add stuffing */
                memmove(buf + 4 + stuffing_len, buf + 4, header_len - 4);
                buf[3] |= 0x20;
                buf[4] = stuffing_len - 1;
                if (stuffing_len >= 2) {
                    buf[5] = 0x00;
                    memset(buf + 6, 0xff, stuffing_len - 2);
                }
            }
        }
        memcpy(buf + TS_PACKET_SIZE - len, payload, len);
        payload += len;
        payload_size -= len;
        put_buffer(s->pb, buf, TS_PACKET_SIZE);
        ts->cur_pcr += (TS_PACKET_SIZE+write_pcr)*8*90000LL / ts->mux_rate;
    }
    put_flush_packet(s->pb);
}

static void put_padding_packet(uint8_t **pes_payload, int packet_bytes)
{
    int i;

    bytestream_put_be32(pes_payload, PADDING_STREAM);
    bytestream_put_be16(pes_payload, packet_bytes - 6);
    packet_bytes -= 6;

    for(i=0;i<packet_bytes;i++)
        bytestream_put_byte(pes_payload, 0xff);
}

/* flush the packet on stream stream_index */
static int flush_packet(AVFormatContext *ctx, int stream_index,
                        int64_t pts, int64_t dts, int64_t pcr, int trailer_size)
{
    MpegTSWriteStream *stream = ctx->streams[stream_index]->priv_data;
    int payload_size, stuffing_size, i;
    int packet_size;
    int zero_trail_bytes = 0;
    int pad_packet_bytes = 0;
    uint8_t *q = stream->payload;

    packet_size = DEFAULT_PES_PAYLOAD_SIZE;

    if (packet_size > 0) {
        int pes_size = ff_pes_write_buf(ctx, stream_index, stream->payload,
                                        &pts, &dts, trailer_size,
                                        &packet_size, &pad_packet_bytes,
                                        &payload_size, &stuffing_size);
        if(pes_size < 0)
            return -1;
        q += pes_size;
    }else{
        payload_size=
        stuffing_size= 0;
    }

    if (pad_packet_bytes > 0)
        put_padding_packet(&q, pad_packet_bytes);

    for(i=0;i<zero_trail_bytes;i++)
        bytestream_put_byte(&q, 0x00);

    mpegts_write_pes(ctx, ctx->streams[stream_index], pcr,
                     stream->payload, q - stream->payload);

    return payload_size - stuffing_size;
}

static int output_packet(AVFormatContext *ctx, int flush){
    MpegTSWrite *s = ctx->priv_data;
    StreamInfo *stream;
    int es_size;
    int best_i= -1;
    int64_t pcr= s->last_pcr;

    if ((es_size = ff_pes_output_packet(ctx, DEFAULT_PES_PAYLOAD_SIZE,
                                        &pcr, &best_i, flush, flush_packet)) <= 0)
        return es_size;

    stream= ctx->streams[best_i]->priv_data;
    stream->buffer_index += es_size;
    while(stream->premux_packet && stream->premux_packet->unwritten_size <= es_size){
        es_size -= stream->premux_packet->unwritten_size;
        stream->premux_packet= stream->premux_packet->next;
    }
    if(es_size)
        stream->premux_packet->unwritten_size -= es_size;

    if(ff_pes_remove_decoded_packets(ctx, s->last_pcr) < 0)
        return -1;

    return 1;
}


static int mpegts_write_packet(AVFormatContext *ctx, AVPacket *pkt)
{
    ff_pes_write_packet(ctx, pkt, 0);
    for(;;){
        int ret = output_packet(ctx, 0);
        if(ret<=0)
            return ret;
    }
}

static int mpegts_write_end(AVFormatContext *s)
{
    MpegTSWrite *ts = s->priv_data;
    MpegTSService *service;
    int i;

    for(;;){
        int ret= output_packet(s, 1);
        if(ret<0)
            return ret;
        else if(ret==0)
            break;
    }

    ff_pes_muxer_end(s);
    for(i = 0; i < ts->nb_services; i++) {
        service = ts->services[i];
        av_freep(&service->provider_name);
        av_freep(&service->name);
        av_free(service);
    }
    av_free(ts->services);

    return 0;
}

AVOutputFormat mpegts_muxer = {
    "mpegts",
    NULL_IF_CONFIG_SMALL("MPEG-2 transport stream format"),
    "video/x-mpegts",
    "ts,m2t",
    sizeof(MpegTSWrite),
    CODEC_ID_MP2,
    CODEC_ID_MPEG2VIDEO,
    mpegts_write_header,
    mpegts_write_packet,
    mpegts_write_end,
};
