/*
 * MMS protocol over HTTP
 * Copyright (c) 2010 Zhentan Feng <spyfeng at gmail dot com>
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
#include <string.h>

#define CHUNK_TYPE_DATA        0x4424
#define CHUNK_TYPE_ASF_HEADER  0x4824

#define USERAGENT "User-Agent: NSPlayer/4.1.0.3856\r\n"
#define CLIENTGUID "Pragma: xClientGUID={c77e7400-738a-11d2-9add-0020af0a3278}\r\n"

static const char* mmsh_first_request =
    "GET %s HTTP/1.0\r\n"
    "Accept: */*\r\n"
    USERAGENT
    "Host: %s:%d\r\n"
    "Pragma: no-cache,rate=1.000000,stream-time=0,stream-offset=0:0,request-context=%u,max-duration=0\r\n"
    CLIENTGUID
    "Connection: Close\r\n\r\n";

static const char* mmsh_seekable_request =
    "GET %s HTTP/1.0\r\n"
    "Accept: */*\r\n"
    USERAGENT
    "Host: %s:%d\r\n"
    "Pragma: no-cache,rate=1.000000,stream-time=%u,stream-offset=%u:%u,request-context=%u,max-duration=%u\r\n"
    CLIENTGUID
    "Pragma: xPlayStrm=1\r\n"
    "Pragma: stream-switch-count=%d\r\n"
    "Pragma: stream-switch-entry=%s\r\n" /*  ffff:1:0 ffff:2:0 */
    "Connection: Close\r\n\r\n";

static const char* mmsh_live_request =
    "GET %s HTTP/1.0\r\n"
    "Accept: */*\r\n"
    USERAGENT
    "Host: %s:%d\r\n"
    "Pragma: no-cache,rate=1.000000,request-context=%u\r\n"
    "Pragma: xPlayStrm=1\r\n"
    CLIENTGUID
    "Pragma: stream-switch-count=%d\r\n"
    "Pragma: stream-switch-entry=%s\r\n"
    "Connection: Close\r\n\r\n";

typedef struct
{
    URLContext *mms_hd;
    uint8_t out_buffer[1024];             ///< Buffer for outgoing packet.
    uint8_t in_buffer[1024]; //TODO, maybe reused by out_buffer.
    uint8_t *asf_header_pos;
    uint8_t *http_header_data;
    int content_length;
    int asf_header_len;
    int asf_header_remaining_len;
    int asf_data_remaining_len;

    char path[256];
    char host[128];
    int seekable;
    int stream_num;

    char stream_selection[10 * MAX_STREAMS];
}MMSHContext;

static int mmsh_close(URLContext *h)
{
    MMSHContext *mms = (MMSHContext *)h->priv_data;
    if(mms->mms_hd)
        url_close(mms->mms_hd);
    av_freep(&h->priv_data);
    av_freep(&mms->http_header_data);
    //TODO free other alloced mem.
    return 0;
}

static int send_pack(MMSHContext *mms)
{
    int len, result;
    len = strlen(mms->out_buffer);
    result = url_write(mms->mms_hd, mms->out_buffer, len)
    if(result != len) {
        dprintf(NULL,"send pack failed!return len %d != %d\n", result, len);
        return AVERROR_IO;
    }
    return 0;
}

static char* find_str(char * dst, const char *str, const int len)
{
    char *p = NULL;
    if(strncasecmp(dst, str, len) == 0) {
        p=dst[len];
        while(isspace(*p))
            p++;
    }
    return p;
}

static int get_and_parse_http_header(MMSHContext *mms)
{
    int len = 0, line_num = 0;
    int http_code;
    char content_type[128]={'\0'};
    char *p, *pos;
    for(;;) {
        if(url_read(mms->mms_hd, mms->in_buffer[len], 1) != 1) {
            dprintf(NULL, "recv http header failed!\n");
            return AVERROR(EIO);
        }

        if(mms->in_buffer[len] != 0x0A) {
            len++;
            if(len >= sizeof(mms->in_buffer)) {
                dprintf(NULL, "recv http header overwrite the buffer!\n");
                return -1;
            }
        } else {
            mms->in_buffer[len--] = '\0';
            if (len >= 0 && mms->in_buffer[len] == 0x0D) {
                line_num++;
                mms->in_buffer[len--] = '\0';
                if(len < 0) {
                    return 0; // \r\n\r\n is the end of http header
                } else {
                    len = 0;// begin to read next http header line
                }
            }
            if (line_num == 1) {
                p = mms->in_buffer;
                while(!isspace(*p) && *p != '\0')
                    p++;
                while(isspace(*p))
                    p++;
                http_code = strtol(p, NULL, 10);
                dprintf(NULL, "mmsh protocol http_code=%d\n", http_code);
                if(http_code != 200) {
                    return -1;
                }
            } else {
                if ((p = find_str(mms->in_buffer, "Content-Type:", 13)) != NULL) {
                    strncpy(content_type, p, sizeof(content_type));
                    dprintf(NULL, "Content-Type:%s\n", content_type);
                    if(strcmp(content_type, "application/x-mms-framed") != 0
                        && strcmp(content_type, "application/octet-stream") != 0
                        && strcmp(content_type, "application/vnd.ms.wms-hdr.asfv1") != 0) {
                        return -1;
                    }
                } else if((p = find_str(mms->in_buffer, "Content-Length:", 15)) != NULL) {
                    mms->content_length = atoi(*p);
                } else if(p = find_str(mms->in_buffer, "Pragma:", 7) != NULL) {
                    pos = strstr(p, "features=");
                    if (pos){
                        if(strstr(pos, "seekable")) {
                            mms->seekable = 1;
                        } else if (strstr(pos, "broadcast")) {
                            mms->seekable = 0;
                        }
                    } else {
                        dprintf(NULL, "Can't find features!\n");
                    }
                }
            }
        }
    }
}

static uint16_t http_header_data_parser(MMSHContext *mms)
{
    uint16_t chunk_type;
    int chunk_len;
    int data_len = mms->content_length;
    char *pos = mms->http_header_data;

    while(data_len) {
        chunk_type = AV_RL16(pos);
        chunk_len = AV_RL16(pos + 2);
        if(chunk_type == CHUNK_TYPE_ASF_HEADER) {
            mms->asf_header_pos = pos + 4; // start from $H, ox4824
            mms->asf_header_len = chunk_len;
        }
        data_len -= chunk_len + 4;
        pos += chunk_len + 4;
        if (data_len <= 0) {
            dprintf(NULL, "http header data len is %d for type %x\n", chunk_len, chunk_type);
            return -1;
        }
        if (chunk_type == CHUNK_TYPE_ASF_HEADER) {
            //TODO parse asf header
        }
    }
    return 0;
}

static int get_http_header_data(MMSHContext *mms, const int flag)
{
    int res, len;
    if(mms->content_length && flag == 1) {
        mms->http_header_data = av_mallocz(mms->content_length);
        if (!mms->http_header_data)
            return AVERROR(ENOMEM);
        // read the http header data, may contain $H, $M, $P packet.
        // In this situation, it should only has $H packet, ie asf header data.
        res = url_read_complete(mms->mms_hd, mms->http_header_data, mms->content_length);
        if (res != mms->content_length) {
            dprintf(NULL, "recv header data len %d != %d", res, mms->content_length);
            return -1;
        }
    } else  if(flag == 2){
        // skip asf header
        uint16_t type;
        char *tmp = av_mallocz(mms->asf_header_len);
        if (!tmp)
            return AVERROR(ENOMEM);
        res = url_read_complete(mms->mms_hd, tmp, mms->asf_header_len);
        if (res != mms->asf_header_len) {
           dprintf(NULL, "read skipped asf header failed!\n");
           av_free(tmp);
           return -1;
        }
        type = AV_RL16(tmp);
        if (type != CHUNK_TYPE_ASF_HEADER) {
            dprintf(NULL, "cann't skip asf header because we didn't recv it!\n");
            av_free(tmp);
            return -1;
        }
        av_free(tmp);
    } else {
        dprintf(NULL, "http response has no data!\n");
        return -1;
    }
    return 0;
}

static int get_http_answer(MMSHContext *mms, const int flag)
{
    int result;
    result = get_and_parse_http_header(mms);
    if (result) {
        dprintf(NULL, "http header parser failed!\n");
        return result;
    }

    result = get_http_header_data(mms, flag);
    if (result) {
        dprintf(NULL, "get http header data fialed!\n");
        return result;
    }
    if (flag == 1) {
        result = http_header_data_parser(mms);
        if (result) {
            dprintf(NULL, "http header data parser failed!\n");
            return result;
        }
    }
    return 0;
}

static int mmsh_open(URLContext *h, const char *uri, int flags)
{
    MMSHContext *mms;
    int port, err;
    char tcpname[256];

    h->is_streamed = 1;
    mms = h->priv_data = av_mallocz(sizeof(MMSHContext));
    if (!h->priv_data)
        return AVERROR(ENOMEM);

    ff_url_split(NULL, 0, NULL, 0,
            mms->host, sizeof(mms->host), &port, mms->path,
            sizeof(mms->path), uri);
    if(port<0)
        port = 80; // defaut mmsh protocol port

    ff_url_join(tcpname, sizeof(tcpname), "tcp", NULL, mms->host, port, NULL);
    err = url_open(&mms->mms_hd, tcpname, URL_RDWR);
    if (err)
        goto fail;
    // send describe request
    snprintf (mms->out_buffer, sizeof(mms->out_buffer), mmsh_first_request, mms->path,
            mms->host, port, this->http_request_number++);
    err = send_pack(mms);
    if (err)
        goto fail;
    err = get_http_answer(mms, 1); // TODO match with the first request
    if(err)
        goto fail;

    // send paly request
    if (mms->seekable) {
        snprintf(mms->out_buffer, sizeof(mms->out_buffer), mmsh_seekable_request, mms->path,
            mms->host, port, 0, 0, 0, this->http_request_number++, 0, mms->stream_num, mms->stream_selection);
    } else {
        snprintf(mms->out_buffer, sizeof(mms->out_buffer),, mmsh_live_request, mms->path,
            mms->host, port, this->http_request_number++, mms->stream_num, mms->stream_selection);
    }
    err = send_pack(mms);
    if (err)
        goto fail;
    err = get_http_answer(mms, 2);// TODO mathc with the second request
    if(err)
        goto fail;
    dprintf(NULL, "Leaving mmsh open success.\n");
fail:
    mmsh_close(h);
    dprintf(NULL, "Leaving mmsh open (failure: %d)\n", err);
    return err;
}

static int read_data(MMSHContext *mms, char *buf, int size)
{
}

static int read_data_packet(MMSHContext *mms)
{
}

static int mmsh_read(URLContext *h, uint8_t *buf, int size)
{
    int size_to_copy;
    MMSHContext *mms = h->priv_data;

    if (mms->asf_header_remaining_len != 0) {
        // copy asf header into buffer
        char *pos;
        size_to_copy = FFMIN(size, mms->asf_header_remaining_len);
        pos = mms->asf_header_pos + (mms->asf_header_len - mms->asf_header_remaining_len);
        memcpy(buf, pos, size_to_copy);
        mms->asf_header_remaining_len -= size_to_copy;
        if (mms->asf_header_remaining_len == 0) {
            av_freep(&mms->http_header_data); // which contains asf header
        }
    } else if (mms->asf_data_remaining_len){
        read_data(mms, buf, size);
    } else {
        // read data packet from network
        read_data_packet(mms);
    }
}

URLProtocol mmsh_protocol = {
    "mmsh, mms",
    mmsh_open,
    mmsh_read,
    NULL, // write
    NULL, // seek
    mmsh_close,
};

