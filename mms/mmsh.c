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
#include "mms.h"
#include "internal.h"
#include "libavutil/intreadwrite.h"
#include <string.h>
#include "libavutil/avstring.h"
#include "asf.h"
#include "http.h"

#define DEBUG
#define CHUNK_TYPE_DATA        0x4424
#define CHUNK_TYPE_ASF_HEADER  0x4824
#define CHUNK_TYPE_END         0x4524
#define CHUNK_TYPE_STREAM_CHANGE 0x4324

#define CHUNK_HEADER_LENGTH 4
#define EXT_HEADER_LENGTH   8

#define USERAGENT "User-Agent: NSPlayer/4.1.0.3856\r\n"
#define CLIENTGUID "Pragma: xClientGUID={c77e7400-738a-11d2-9add-0020af0a3278}\r\n"

static const char* mmsh_first_request =
    "Accept: */*\r\n"
    USERAGENT
    "Host: %s:%d\r\n"
    "Pragma: no-cache,rate=1.000000,stream-time=0,stream-offset=0:0,request-context=%u,max-duration=0\r\n"
    CLIENTGUID
    "Connection: Close\r\n\r\n";

static const char* mmsh_seekable_request =
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
    MMSContext *ff_ctx;
    char location[1024];
    int seekable;
    int stream_num;
    int request_seq;
    int chunk_seq;
}MMSHContext;

static int mmsh_close(URLContext *h)
{
    MMSHContext *mmsh_ctx = (MMSHContext *)h->priv_data;
    MMSContext *mms = mmsh_ctx->ff_ctx;
    if(mms->mms_hd)
        url_close(mms->mms_hd);
    av_freep(&mms->asf_header);
    av_freep(&mms);
    av_freep(&h->priv_data);
    return 0;
}

static int get_chunk_header(MMSHContext *mmsh_ctx, int *len)
{
    MMSContext *mms = mmsh_ctx->ff_ctx;
    uint8_t chunk_header[CHUNK_HEADER_LENGTH];
    uint8_t ext_header[EXT_HEADER_LENGTH];
    int chunk_type;
    int chunk_len, res, ext_header_len = 0;

    res = url_read(mms->mms_hd, chunk_header, CHUNK_HEADER_LENGTH);
    if (res != CHUNK_HEADER_LENGTH) { // TODO extact common log code as macro define
        dprintf(NULL, "read data packet  header failed!\n");
        return AVERROR(EIO);
    }
    chunk_type = AV_RL16(chunk_header);
    chunk_len = AV_RL16(chunk_header + 2);
    if (chunk_type == CHUNK_TYPE_END ||chunk_type == CHUNK_TYPE_STREAM_CHANGE) {
        ext_header_len = 4;
    } else if (chunk_type == CHUNK_TYPE_ASF_HEADER || chunk_type == CHUNK_TYPE_DATA) {
        ext_header_len = 8;
    }
    if (ext_header_len) {
        res = url_read(mms->mms_hd, ext_header, ext_header_len);
        if (res != ext_header_len) {
            dprintf(NULL, "read ext header failed!\n");
            return AVERROR(EIO);
        }
    } else {
        dprintf(NULL, "strange chunk type %d\n", chunk_type);
        return -1;
    }
    *len = chunk_len - ext_header_len;
    if (chunk_type == CHUNK_TYPE_END || chunk_type == CHUNK_TYPE_DATA)
        mmsh_ctx->chunk_seq = AV_RL32(ext_header);
    return chunk_type;
}

static int read_data_packet(MMSHContext *mmsh_ctx, const int len)
{
    MMSContext *mms = mmsh_ctx->ff_ctx;
    int res, pad_size = 0;
    res = url_read_complete(mms->mms_hd, mms->in_buffer, len);
    dprintf(NULL, "data packet len = %d\n", len);
    if (res != len) {
        dprintf(NULL, "read data packet failed!\n");
        return AVERROR(EIO);
    }
    if (len > mms->asf_packet_len) {
        dprintf(NULL, "chunk length %d exceed packet length %d\n", len, mms->asf_packet_len);
        return -1;
    } else {
        pad_size = mms->asf_packet_len - len;
        memset(mms->in_buffer + len, 0, pad_size);
    }
    mms->read_in_ptr = mms->in_buffer;
    mms->remaining_in_len = mms->asf_packet_len;
    return 0;
}

static int get_http_header_data(MMSHContext *mmsh_ctx)
{
    MMSContext *mms = mmsh_ctx->ff_ctx;
    int res, len;
    int chunk_type;

    for (;;) {
        len = 0;
        chunk_type = get_chunk_header(mmsh_ctx, &len);
        if (chunk_type < 0) {
            return chunk_type;
        } else if (chunk_type == CHUNK_TYPE_ASF_HEADER){
            // get asf header and stored it
            if (!mms->header_parsed) {
                if (mms->asf_header) {
                    if (len != mms->asf_header_size) {
                        mms->asf_header_size = len;
                        dprintf(NULL, "header len changed form %d to %d\n",
                                mms->asf_header_size, len);
                        av_freep(&mms->asf_header);
                    }
                }
                mms->asf_header = av_mallocz(len);
                if (!mms->asf_header) {
                    return AVERROR(ENOMEM);
                }
            }
            res = url_read_complete(mms->mms_hd, mms->asf_header, len);
            if (res != len) {
                dprintf(NULL, "recv asf header data len %d != %d", res, len);
                return -1;
            }
            mms->asf_header_size = len;
            if (!mms->header_parsed) {
                res = ff_mms_asf_header_parser(mms);
                mms->header_parsed = 1;
                return res;
            }
        } else if (chunk_type == CHUNK_TYPE_DATA) {
            // read data packet and do padding
            return read_data_packet(mmsh_ctx, len);
            break;
        } else {
            if (len) {
                res = url_read_complete(mms->mms_hd, mms->in_buffer, len);
                if (res != len) {
                    dprintf(NULL, "read other chunk type data failed!\n");
                    return AVERROR(EIO);
                } else {
                    dprintf(NULL, "skip chunk type %d \n", chunk_type);
                    continue;
                }
            }
        }
    }
    return 0;
}

static int mmsh_open_cnx(MMSHContext *mmsh_ctx)
{
    MMSContext *mms = mmsh_ctx->ff_ctx;
    int i, port, err, offset = 0;
    char httpname[256], path[256], host[128];
    char stream_selection[10 * MAX_STREAMS];
    char headers[1024];

    if (mms->mms_hd) {
        url_close(mms->mms_hd);
    }

    ff_url_split(NULL, 0, NULL, 0,
            host, sizeof(host), &port, path, sizeof(path), mmsh_ctx->location);
    if(port<0)
        port = 80; // defaut mmsh protocol port
    ff_url_join(httpname, sizeof(httpname), "http", NULL, host, port, path);

    if (url_alloc(&mms->mms_hd, httpname, URL_RDONLY) < 0) {
        return AVERROR(EIO);
    }

    snprintf (headers, sizeof(headers), mmsh_first_request,
        host, port, mmsh_ctx->request_seq++);
    ff_http_set_headers(mms->mms_hd, headers);

    if (url_connect(mms->mms_hd)) {
          return AVERROR(EIO);
    }
    err = get_http_header_data(mmsh_ctx);
    if (err) {
        dprintf(NULL, "get http header data fialed!\n");
        return (err);
    }

    // close the socket and then reopen it for sending the second play request.
    url_close(mms->mms_hd);
    memset(headers, 0, sizeof(headers));
    if (url_alloc(&mms->mms_hd, httpname, URL_RDONLY) < 0) {
        return AVERROR(EIO);
    }
    for (i = 0; i < mms->stream_num; i++) {
        err = snprintf(stream_selection + offset, sizeof(stream_selection) - offset,
                          "ffff:%d:0 ", mms->streams[i].id);
        if (err < 0)
            return err;
        offset += err;
    }
    // send paly request
    if (mmsh_ctx->seekable) {
        err = snprintf(headers, sizeof(headers), mmsh_seekable_request,
            host, port, 0, 0, 0, mmsh_ctx->request_seq++, 0, mms->stream_num, stream_selection);
    } else {
        err = snprintf(headers, sizeof(headers), mmsh_live_request,
            host, port, mmsh_ctx->request_seq++, mms->stream_num, stream_selection);
    }
    if (err < 0) {
        dprintf(NULL, "build play request failed!\n");
        return err;
    }
    dprintf(NULL, "out_buffer is %s", headers);
    ff_http_set_headers(mms->mms_hd, headers);

    if (url_connect(mms->mms_hd)) {
          return AVERROR(EIO);
    }

   err = get_http_header_data(mmsh_ctx);
    if (err) {
        dprintf(NULL, "get http header data fialed!\n");
        return (err);
    }

    return 0;
}

static int mmsh_open(URLContext *h, const char *uri, int flags)
{
    MMSHContext *mmsh_ctx;
    MMSContext *mms;
    int err;
    mmsh_ctx = h->priv_data = av_mallocz(sizeof(MMSHContext));
    if (!h->priv_data)
        return AVERROR(ENOMEM);
    mms = mmsh_ctx->ff_ctx = av_mallocz(sizeof(MMSContext));
    if (!mmsh_ctx->ff_ctx)
        return AVERROR(ENOMEM);
    mmsh_ctx->request_seq = h->is_streamed = 1;
    av_strlcpy(mmsh_ctx->location, uri, sizeof(mmsh_ctx->location));

    err =mmsh_open_cnx(mmsh_ctx);
    if (err) {
        dprintf(NULL, "Leaving mmsh open (failure: %d)\n", err);
        mmsh_close(h);
        return err;
    }
    dprintf(NULL, "Leaving mmsh open success.\n");
    return 0;
}

static int handle_chunk_type(MMSHContext *mmsh_ctx)
{
    MMSContext *mms = mmsh_ctx->ff_ctx;
    int res, len = 0;
    int chunk_type;
    chunk_type = get_chunk_header(mmsh_ctx, &len);

    if(chunk_type == CHUNK_TYPE_END) {
        if (mmsh_ctx->chunk_seq == 0) {
            dprintf(NULL, "The stream is end.\n");
            return -1;
        }
        // reconnect
        mmsh_ctx->request_seq = 1;
        if ((res = mmsh_open_cnx(mmsh_ctx)) !=0) {
            dprintf(NULL, "Reconnect failed!\n");
            return res;
        }
    } else if (chunk_type == CHUNK_TYPE_STREAM_CHANGE) {
        mms->header_parsed = 0;
        if ((res = get_http_header_data(mmsh_ctx)) !=0) {
            dprintf(NULL,"stream changed! get new header failed!\n");
            return res;
        }
    } else if (chunk_type == CHUNK_TYPE_DATA) {
        return read_data_packet(mmsh_ctx, len);
    } else {
        dprintf(NULL, "recv other type packet %d\n", chunk_type);
        return -1;
    }
    return 0;
}

static int mmsh_read(URLContext *h, uint8_t *buf, int size)
{
    int res = 0;
    MMSHContext *mmsh_ctx = h->priv_data;
    MMSContext *mms = mmsh_ctx->ff_ctx;
    do{
        if (mms->asf_header_read_size < mms->asf_header_size) {
            // copy asf header into buffer
            res = ff_mms_read_header(mms, buf, size);
        } else if (mms->remaining_in_len){
            res = ff_mms_read_data(mms, buf, size);
        } else {
             // read data packet from network
            res = handle_chunk_type(mmsh_ctx);
            if (res == 0) {
                res = ff_mms_read_data(mms, buf, size);
            } else {
                dprintf(NULL, "other situation!\n");
            }
        }
    }while(!res);
    return res;
}

URLProtocol mmsh_protocol = {
    "mmsh",
    mmsh_open,
    mmsh_read,
    NULL, // write
    NULL, // seek
    mmsh_close,
};

