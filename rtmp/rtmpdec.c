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
#include "libavutil/lfg.h"
#include "libavutil/sha.h"
#include "avformat.h"

#include <unistd.h>
#include <sys/time.h>
#include "network.h"

#include "flv.h"
#include "rtmp.h"
#include "rtmppkt.h"

typedef enum {
    STATE_START,
    STATE_HANDSHAKED,
    STATE_CONNECTING,
    STATE_PLAYING,
} ClientState;

typedef struct RTMPState {
    URLContext*   rtmp_hd;
    RTMPPacket    prev_pkt[2][RTMP_CHANNELS];
    int           chunk_size;
    char          playpath[256];
    ClientState   state;
    int           main_stream_id;
    uint8_t*      flv_data;
    int           flv_size;
    ByteIOContext pb;
    int           wrong_dts;
} RTMPState;

#define PLAYER_KEY_OPEN_PART_LEN 30
static const uint8_t rtmp_player_key[] =
{
    0x47, 0x65, 0x6E, 0x75, 0x69, 0x6E, 0x65, 0x20, 0x41, 0x64, 0x6F, 0x62, 0x65,
    0x20, 0x46, 0x6C, 0x61, 0x73, 0x68, 0x20, 0x50, 0x6C, 0x61, 0x79, 0x65, 0x72,
    0x20, 0x30, 0x30, 0x31, 0xF0, 0xEE, 0xC2, 0x4A, 0x80, 0x68, 0xBE, 0xE8, 0x2E,
    0x00, 0xD0, 0xD1, 0x02, 0x9E, 0x7E, 0x57, 0x6E, 0xEC, 0x5D, 0x2D, 0x29, 0x80,
    0x6F, 0xAB, 0x93, 0xB8, 0xE6, 0x36, 0xCF, 0xEB, 0x31, 0xAE
};

#define SERVER_KEY_OPEN_PART_LEN 36
static const uint8_t rtmp_server_key[] =
{
    0x47, 0x65, 0x6E, 0x75, 0x69, 0x6E, 0x65, 0x20, 0x41, 0x64, 0x6F, 0x62, 0x65,
    0x20, 0x46, 0x6C, 0x61, 0x73, 0x68, 0x20, 0x4D, 0x65, 0x64, 0x69, 0x61, 0x20,
    0x53, 0x65, 0x72, 0x76, 0x65, 0x72, 0x20, 0x30, 0x30, 0x31, // Genuine Adobe Flash Media Server 001

    0xF0, 0xEE, 0xC2, 0x4A, 0x80, 0x68, 0xBE, 0xE8, 0x2E, 0x00, 0xD0, 0xD1, 0x02,
    0x9E, 0x7E, 0x57, 0x6E, 0xEC, 0x5D, 0x2D, 0x29, 0x80, 0x6F, 0xAB, 0x93, 0xB8,
    0xE6, 0x36, 0xCF, 0xEB, 0x31, 0xAE
};

static void gen_connect(AVFormatContext *s, RTMPState *rt, const char *proto,
                        const char *host, int port, const char *app)
{
    RTMPPacket pkt;
    uint8_t ver[32], *p;
    char tcurl[512];
    double num = 1.0;
    uint8_t bool;

    rtmp_packet_create(&pkt, RTMP_VIDEO_CHANNEL, RTMP_PT_INVOKE, 0, 4096);
    p = pkt.data;

    snprintf(tcurl, sizeof(tcurl), "%s://%s:%d/%s", proto, host, port, app);
    rtmp_amf_write_tag(&p, AMF_STRING, "connect");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);
    rtmp_amf_write_tag(&p, AMF_OBJECT, NULL);
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "app");
    rtmp_amf_write_tag(&p, AMF_STRING, app);

    snprintf(ver, sizeof(ver), "%s %d,%d,%d,%d", RTMP_CLIENT_PLATFORM, RTMP_CLIENT_VER1,
             RTMP_CLIENT_VER2, RTMP_CLIENT_VER3, RTMP_CLIENT_VER4);
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "flashVer");
    rtmp_amf_write_tag(&p, AMF_STRING, ver);
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "tcUrl");
    rtmp_amf_write_tag(&p, AMF_STRING, tcurl);
    bool = 0;
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "fpad");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &bool);
    num = 15.0;
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "capabilities");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);
    num = 1639.0;
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "audioCodecs");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);
    num = 252.0;
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "videoCodecs");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);
    num = 1.0;
    rtmp_amf_write_tag(&p, AMF_STRING_IN_OBJECT, "videoFunction");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);
    rtmp_amf_write_tag(&p, AMF_OBJECT_END, NULL);

    pkt.data_size = p - pkt.data;

    rtmp_packet_write(s, rt->rtmp_hd, &pkt, rt->chunk_size, rt->prev_pkt[1]);
}

static void gen_create_stream(AVFormatContext *s, RTMPState *rt)
{
    RTMPPacket pkt;
    uint8_t *p;
    double num;

    av_log(s, AV_LOG_DEBUG, "Creating stream...\n");
    rtmp_packet_create(&pkt, RTMP_VIDEO_CHANNEL, RTMP_PT_INVOKE, 0, 25);

    num = 3.0;
    p = pkt.data;
    rtmp_amf_write_tag(&p, AMF_STRING, "createStream");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);
    rtmp_amf_write_tag(&p, AMF_NULL, NULL);

    rtmp_packet_write(s, rt->rtmp_hd, &pkt, rt->chunk_size, rt->prev_pkt[1]);
    rtmp_packet_destroy(&pkt);
}

static void gen_play(AVFormatContext *s, RTMPState *rt)
{
    RTMPPacket pkt;
    uint8_t *p;
    double num;

    av_log(s, AV_LOG_DEBUG, "Sending play command for '%s'\n", rt->playpath);
    rtmp_packet_create(&pkt, RTMP_VIDEO_CHANNEL, RTMP_PT_INVOKE, 0,
                       29 + strlen(rt->playpath));
    pkt.extra = rt->main_stream_id;

    num = 0.0;
    p = pkt.data;
    rtmp_amf_write_tag(&p, AMF_STRING, "play");
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);
    rtmp_amf_write_tag(&p, AMF_NULL, NULL);
    rtmp_amf_write_tag(&p, AMF_STRING, rt->playpath);
    num = 0.0;
    rtmp_amf_write_tag(&p, AMF_NUMBER, &num);

    rtmp_packet_write(s, rt->rtmp_hd, &pkt, rt->chunk_size, rt->prev_pkt[1]);
    rtmp_packet_destroy(&pkt);

    // set client buffer time disguised in ping packet
    rtmp_packet_create(&pkt, RTMP_NETWORK_CHANNEL, RTMP_PT_PING, 1, 10);

    p = pkt.data;
    bytestream_put_be16(&p, 3);
    bytestream_put_be32(&p, 1);
    bytestream_put_be32(&p, 256); //TODO: what is a good value here?

    rtmp_packet_write(s, rt->rtmp_hd, &pkt, rt->chunk_size, rt->prev_pkt[1]);
    rtmp_packet_destroy(&pkt);
}

static void gen_pong(AVFormatContext *s, RTMPState *rt, RTMPPacket *ppkt)
{
    RTMPPacket pkt;
    uint8_t *p;

    rtmp_packet_create(&pkt, RTMP_NETWORK_CHANNEL, RTMP_PT_PING, ppkt->timestamp + 1, 6);
    p = pkt.data;
    bytestream_put_be16(&p, 7);
    bytestream_put_be32(&p, AV_RB32(ppkt->data+2) + 1);
    rtmp_packet_write(s, rt->rtmp_hd, &pkt, rt->chunk_size, rt->prev_pkt[1]);
    rtmp_packet_destroy(&pkt);
}

//TODO: Move HMAC code somewhere. Eventually.
#define HMAC_IPAD_VAL 0x36
#define HMAC_OPAD_VAL 0x5C

static void rtmp_calc_digest(const uint8_t *src, int len, int gap,
                             const uint8_t *key, int keylen, uint8_t *dst)
{
    struct AVSHA *sha;
    uint8_t hmac_buf[64+32];
    int i;

    sha = av_mallocz(av_sha_size);

    memset(hmac_buf, 0, 64);
    if (keylen < 64)
        memcpy(hmac_buf, key, keylen);
    else {
        av_sha_init(sha, 256);
        av_sha_update(sha,key, keylen);
        av_sha_final(sha, hmac_buf);
    }
    for (i = 0; i < 64; i++)
        hmac_buf[i] ^= HMAC_IPAD_VAL;

    av_sha_init(sha, 256);
    av_sha_update(sha, hmac_buf, 64);
    if (gap <= 0)
        av_sha_update(sha, src, len);
    else { //skip 32 bytes used for storing digest
        av_sha_update(sha, src, gap);
        av_sha_update(sha, src + gap + 32, len - gap - 32);
    }
    av_sha_final(sha, hmac_buf + 64);

    for (i = 0; i < 64; i++)
        hmac_buf[i] ^= HMAC_IPAD_VAL ^ HMAC_OPAD_VAL; //reuse XORed key for opad
    av_sha_init(sha, 256);
    av_sha_update(sha, hmac_buf, 64+32);
    av_sha_final(sha, dst);

    av_free(sha);
}

static int rtmp_handshake_imprint_with_digest(uint8_t *buf)
{
    int i, digest_pos = 0;

    for (i = 8; i < 12; i++)
        digest_pos += buf[i];
    digest_pos = (digest_pos % 728) + 12;

    rtmp_calc_digest(buf, RTMP_HANDSHAKE_PACKET_SIZE, digest_pos,
                     rtmp_player_key, PLAYER_KEY_OPEN_PART_LEN,
                     buf + digest_pos);
    return digest_pos;
}

static int rtmp_validate_digest(uint8_t *buf, int off)
{
    int i, digest_pos = 0;
    uint8_t digest[32];

    for (i = 0; i < 4; i++)
        digest_pos += buf[i + off];
    digest_pos = (digest_pos % 728) + off + 4;

    rtmp_calc_digest(buf, RTMP_HANDSHAKE_PACKET_SIZE, digest_pos,
                     rtmp_server_key, SERVER_KEY_OPEN_PART_LEN,
                     digest);
    if (!memcmp(digest, buf + digest_pos, 32))
        return digest_pos;
    return 0;
}

static int rtmp_handshake(AVFormatContext *s, RTMPState *rt)
{
    AVLFG rnd;
    uint8_t tosend    [RTMP_HANDSHAKE_PACKET_SIZE+1];
    uint8_t clientdata[RTMP_HANDSHAKE_PACKET_SIZE];
    uint8_t serverdata[RTMP_HANDSHAKE_PACKET_SIZE+1];
    int i;
    int server_pos, client_pos;
    uint8_t digest[32];

    av_log(s, AV_LOG_DEBUG, "Handshaking...\n");

    av_lfg_init(&rnd, 0xDEADC0DE);
    // generate handshake packet - 1536 bytes of pseudorandom data
    tosend[0] = 3; //unencrypted data
    memset(tosend+1, 0, 4);
    //write client "version"
    tosend[5] = RTMP_CLIENT_VER1;
    tosend[6] = RTMP_CLIENT_VER2;
    tosend[7] = RTMP_CLIENT_VER3;
    tosend[8] = RTMP_CLIENT_VER4;
    for (i = 9; i <= RTMP_HANDSHAKE_PACKET_SIZE; i++)
        tosend[i] = av_lfg_get(&rnd) >> 24;
    client_pos = rtmp_handshake_imprint_with_digest(tosend + 1);

    url_write(rt->rtmp_hd, tosend, RTMP_HANDSHAKE_PACKET_SIZE + 1);
    i = url_read_complete(rt->rtmp_hd, serverdata, RTMP_HANDSHAKE_PACKET_SIZE + 1);
    if (i != RTMP_HANDSHAKE_PACKET_SIZE + 1) {
        av_log(s, AV_LOG_ERROR, "Cannot read RTMP handshake response\n");
        return -1;
    }
    i = url_read_complete(rt->rtmp_hd, clientdata, RTMP_HANDSHAKE_PACKET_SIZE);
    if (i != RTMP_HANDSHAKE_PACKET_SIZE) {
        av_log(s, AV_LOG_ERROR, "Cannot read RTMP handshake response\n");
        return -1;
    }

    av_log(s, AV_LOG_DEBUG, "Server version %d.%d.%d.%d\n",
           serverdata[5], serverdata[6], serverdata[7], serverdata[8]);

    server_pos = rtmp_validate_digest(serverdata + 1, 772);
    if (!server_pos) {
        server_pos = rtmp_validate_digest(serverdata + 1, 8);
        if (!server_pos) {
            av_log(s, AV_LOG_ERROR, "Server response validating failed\n");
            return -1;
        }
    }

    rtmp_calc_digest(tosend + 1 + client_pos, 32, 0,
                     rtmp_server_key, sizeof(rtmp_server_key),
                     digest);
    rtmp_calc_digest(clientdata, RTMP_HANDSHAKE_PACKET_SIZE-32, 0,
                     digest, 32,
                     digest);
    if (memcmp(digest, clientdata + RTMP_HANDSHAKE_PACKET_SIZE - 32, 32)) {
        av_log(s, AV_LOG_ERROR, "Signature mismatch\n");
        return -1;
    }

    for (i = 0; i < RTMP_HANDSHAKE_PACKET_SIZE; i++)
        tosend[i] = av_lfg_get(&rnd) >> 24;
    rtmp_calc_digest(serverdata + 1 + server_pos, 32, 0,
                     rtmp_player_key, sizeof(rtmp_player_key),
                     digest);
    rtmp_calc_digest(tosend,  RTMP_HANDSHAKE_PACKET_SIZE - 32, 0,
                     digest, 32,
                     tosend + RTMP_HANDSHAKE_PACKET_SIZE - 32);

    // write reply back to server
    url_write(rt->rtmp_hd, tosend, RTMP_HANDSHAKE_PACKET_SIZE);
    return 0;
}

static int rtmp_probe(AVProbeData *p)
{
    if (av_strstart(p->filename, "rtmp:", NULL))
        return AVPROBE_SCORE_MAX;
    return 0;
}

static int rtmp_read_header(AVFormatContext *s,
                            AVFormatParameters *ap)
{
    RTMPState *rt = s->priv_data;
    char proto[8], hostname[256], path[512], app[128], *fname;
    int port;
    uint8_t buf[2048];
    AVStream *st;

    url_split(proto, sizeof(proto), NULL, 0, hostname, sizeof(hostname), &port,
              path, sizeof(path), s->filename);

    if(port == -1)
        port = RTMP_DEFAULT_PORT;
    snprintf(buf, sizeof(buf), "tcp://%s:%d", hostname, port);
    url_open(&rt->rtmp_hd, buf, URL_RDWR);
    rt->state = STATE_START;
    if (rtmp_handshake(s, rt))
        return -1;

    rt->chunk_size = 128;
    rt->state = STATE_HANDSHAKED;
    //extract "app" part from path
    if (!strncmp(path, "/ondemand/", 10)) {
        fname = path + 10;
        memcpy(app, "ondemand", 9);
    } else {
        char *p = strchr(path + 1, '/');
        if (!p) {
            fname = path + 1;
            app[0] = '\0';
        } else {
            fname = strchr(p + 1, '/');
            if (!fname) {
                fname = p + 1;
                strncpy(app, path + 1, p - path - 1);
            } else {
                fname++;
                strncpy(app, path + 1, fname - path - 2);
            }
        }
    }
    if (!strcmp(fname + strlen(fname) - 4, ".f4v") ||
        !strcmp(fname + strlen(fname) - 4, ".mp4")) {
        memcpy(rt->playpath, "mp4:", 5);
    } else {
        rt->playpath[0] = ':';
        rt->playpath[0] = 0;
    }
    strncat(rt->playpath, fname, sizeof(rt->playpath) - 5);

    av_log(s, AV_LOG_DEBUG, "Proto = %s, path = %s, app = %s, fname = %s\n",
           proto, path, app, rt->playpath);
    gen_connect(s, rt, proto, hostname, port, app);

    st = av_new_stream(s, 0);
    if (!st)
        return -1;
    st->codec->codec_type = CODEC_TYPE_VIDEO;
    av_set_pts_info(st, 32, 1, 1000); /* 32 bit pts in ms */

    st = av_new_stream(s, 1);
    if (!st)
        return -1;
    st->codec->codec_type = CODEC_TYPE_AUDIO;
    av_set_pts_info(st, 32, 1, 1000); /* 32 bit pts in ms */

    init_put_byte(&rt->pb, NULL, 0, 0, NULL, NULL, NULL, NULL);
    return 0;
}

static int rtmp_parse_result(AVFormatContext *s, RTMPState *rt, RTMPPacket *pkt)
{
    int i, t;

    switch (pkt->type) {
    case RTMP_PT_CHUNK_SIZE:
        if (pkt->data_size != 4) {
            av_log(s, AV_LOG_ERROR, "Chunk size change packet is not 4 (%d)\n",
                   pkt->data_size);
            return -1;
        }
        rt->chunk_size = AV_RB32(pkt->data);
        av_log(s, AV_LOG_DEBUG, "New chunk size = %d\n", rt->chunk_size);
        break;
    case RTMP_PT_PING:
        t = AV_RB16(pkt->data);
        if (t == 6)
            gen_pong(s, rt, pkt);
        break;
    case RTMP_PT_INVOKE:
        if (!memcmp(pkt->data, "\002\000\006_error", 9)) {//TODO: search data for error description
            return -1;
        }
        if (!memcmp(pkt->data, "\002\000\007_result", 10)) {
            switch (rt->state) {
            case STATE_HANDSHAKED:
                gen_create_stream(s, rt);
                rt->state = STATE_CONNECTING;
                break;
            case STATE_CONNECTING:
                //extract a number from result
                if (pkt->data[10] || pkt->data[19] != 5 || pkt->data[20])
                    av_log(s, AV_LOG_WARNING, "Unexpected reply on connect()\n");
                else
                    rt->main_stream_id = (int) av_int2dbl(AV_RB64(pkt->data + 21));
                gen_play(s, rt);
                rt->state = STATE_PLAYING;
                break;
            }
        }
        if (!memcmp(pkt->data, "\002\000\008onStatus", 11)) {
            //TODO: catch stream close event
        }
        break;
    }
    return 0;
}

static int rtmp_read_packet(AVFormatContext *s, AVPacket *pkt)
{
    RTMPState *rt = s->priv_data;
    struct timespec ts;
    int ret;

    ts.tv_sec = 0;
    ts.tv_nsec = 500000000;

    while (url_ftell(&rt->pb) == rt->flv_size) {
        RTMPPacket rpkt;
        int has_data = 0;
        if ((ret = rtmp_packet_read(s, rt->rtmp_hd, &rpkt,
                                    rt->chunk_size, rt->prev_pkt[0])) != 0) {
            if (ret > 0) {
                nanosleep(&ts, NULL);
                continue;
            } else {
                return AVERROR(EIO);
            }
        }

        ret = rtmp_parse_result(s, rt, &rpkt);
        if (ret < 0) {//serious error in packet
            rtmp_packet_destroy(&rpkt);
            return -1;
        }
        if (!rpkt.data_size) {
            rtmp_packet_destroy(&rpkt);
            continue;
        }
        if (rpkt.type == RTMP_PT_VIDEO || rpkt.type == RTMP_PT_AUDIO
         || rpkt.type == RTMP_PT_NOTIFY) {
            uint8_t *p;

            rt->flv_size = rpkt.data_size + 15;
            rt->flv_data = p = av_realloc(rt->flv_data, rt->flv_size);
            bytestream_put_be32(&p, 0);
            bytestream_put_byte(&p, rpkt.type);
            bytestream_put_be24(&p, rpkt.data_size);
            bytestream_put_be24(&p, rpkt.extra);
            bytestream_put_byte(&p, rpkt.extra >> 24);
            bytestream_put_be24(&p, 0);
            bytestream_put_buffer(&p, rpkt.data, rpkt.data_size);
            has_data = 1;
        } else if (rpkt.type == RTMP_PT_METADATA) {
            rt->flv_size = rpkt.data_size + 4;
            rt->flv_data = av_realloc(rt->flv_data, rt->flv_size);
            memset(rt->flv_data, 0, 4); //previous packet size
            memcpy(rt->flv_data + 4, rpkt.data, rpkt.data_size);
            has_data = 1;
        }
        rtmp_packet_destroy(&rpkt);
        if (has_data) {
            init_put_byte(&rt->pb, rt->flv_data, rt->flv_size, 0, NULL, NULL, NULL, NULL);
            break;
        }
    }
    ret = ff_flv_read_packet(s, &rt->pb, pkt, &rt->wrong_dts);
    if (ret == AVERROR_EOF)
        return AVERROR(EAGAIN);
    return ret;
}

static int rtmp_read_close(AVFormatContext *s)
{
    RTMPState *rt = s->priv_data;

    av_freep(&rt->flv_data);
    url_close(rt->rtmp_hd);
    return 0;
}

AVInputFormat rtmp_demuxer = {
    "rtmp",
    NULL_IF_CONFIG_SMALL("RTMP"),
    sizeof(RTMPState),
    rtmp_probe,
    rtmp_read_header,
    rtmp_read_packet,
    rtmp_read_close,
};
