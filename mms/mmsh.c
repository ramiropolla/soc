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
#include "internal.h"
#include "libavutil/intreadwrite.h"
#include <string.h>
#include "libavutil/avstring.h"
#include "asf.h"

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

typedef struct {
    int id;
}MMSStream;

typedef struct
{
    URLContext *mms_hd;
    uint8_t out_buffer[8192];             ///< Buffer for outgoing packet.
    uint8_t in_buffer[8192]; //TODO, maybe reused by out_buffer.
    uint8_t *read_in_ptr;
    MMSStream streams[MAX_STREAMS];

    uint8_t *asf_header;
    int asf_header_size;
    int asf_header_read_size;
    int asf_data_remaining_len;
    int asf_packet_len;


    char location[1024];
    int seekable;
    int stream_num;
    int request_seq;
    int chunk_seq;
    int is_header_parsed;
}MMSHContext;

static int mmsh_close(URLContext *h)
{
    MMSHContext *mms = (MMSHContext *)h->priv_data;
    if(mms->mms_hd)
        url_close(mms->mms_hd);
    av_freep(&h->priv_data);
    av_freep(&mms->asf_header);
    //TODO free other alloced mem.
    return 0;
}

static int send_pack(MMSHContext *mms)
{
    int len, result;
    len = strlen(mms->out_buffer);
    result = url_write(mms->mms_hd, mms->out_buffer, len);
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
        p = &dst[len];
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
    //uint8_t tmp_buf[8192];
    //url_read_complete(mms->mms_hd, tmp_buf, sizeof(tmp_buf));
    for(;;) {
        if(url_read(mms->mms_hd, &mms->in_buffer[len], 1) != 1) {
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
                } else if((p = find_str(mms->in_buffer, "Pragma:", 7)) != NULL) {
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

static int asf_header_parser(MMSHContext * mms)
{
    uint8_t *p = mms->asf_header;
    uint8_t *end;
    int flags, stream_id, real_header_size;
    mms->stream_num = 0;

    if (mms->asf_header_size < sizeof(ff_asf_guid) * 2 + 22 ||
        memcmp(p, ff_asf_header, sizeof(ff_asf_guid)))
        return -1;

    real_header_size = AV_RL64(p + sizeof(ff_asf_guid));
    end = mms->asf_header + real_header_size;

    p += sizeof(ff_asf_guid) + 14;
    while(end - p >= sizeof(ff_asf_guid) + 8) {
        uint64_t chunksize = AV_RL64(p + sizeof(ff_asf_guid));
        if (!chunksize || chunksize > end - p) {
            dprintf(NULL, "chunksize is exceptional value:%d!\n", chunksize);
            return -1;
        }
        if (!memcmp(p, ff_asf_file_header, sizeof(ff_asf_guid))) {
            /* read packet size */
            if (end - p > sizeof(ff_asf_guid) * 2 + 68) {
                mms->asf_packet_len = AV_RL32(p + sizeof(ff_asf_guid) * 2 + 64);
                if (mms->asf_packet_len <= 0 || mms->asf_packet_len > sizeof(mms->in_buffer)) {
                    dprintf(NULL,"Too large packet len:%d"
                        " may overwrite in_buffer when padding", mms->asf_packet_len);
                    return -1;
                }
            }
        } else if (!memcmp(p, ff_asf_stream_header, sizeof(ff_asf_guid))) {
            flags     = AV_RL16(p + sizeof(ff_asf_guid)*3 + 24);
            stream_id = flags & 0x7F;
            if (mms->stream_num < MAX_STREAMS ) {
                mms->streams[mms->stream_num].id = stream_id;
                mms->stream_num++;
            } else {
                dprintf(NULL, "Too many streams.\n");
                return -1;
            }
        }
        p += chunksize;
    }
    return 0;
}

static int get_chunk_header(MMSHContext *mms, int *len)
{
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
        mms->chunk_seq = AV_RL32(ext_header);
    return chunk_type;
}

static int read_data_packet(MMSHContext *mms, const int len)
{
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
    mms->asf_data_remaining_len = mms->asf_packet_len;
    return 0;
}

static int get_http_header_data(MMSHContext *mms)
{
    int res, len;
    int chunk_type;

    for (;;) {
        len = 0;
        chunk_type = get_chunk_header(mms, &len);
        if (chunk_type < 0) {
            return chunk_type;
        } else if (chunk_type == CHUNK_TYPE_ASF_HEADER){
            // get asf header and stored it
            if (!mms->is_header_parsed) {
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
            if (!mms->is_header_parsed) {
                res = asf_header_parser(mms);
                mms->is_header_parsed = 1;
                return res;
            }
        } else if (chunk_type == CHUNK_TYPE_DATA) {
            // read data packet and do padding
            return read_data_packet(mms, len);
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

static int get_http_answer(MMSHContext *mms)
{
    int result;
    result = get_and_parse_http_header(mms);
    if (result) {
        dprintf(NULL, "http header parser failed!\n");
        return result;
    }

    result = get_http_header_data(mms);
    if (result) {
        dprintf(NULL, "get http header data fialed!\n");
        return result;
    }
    return 0;
}

static int mmsh_open_cnx(MMSHContext *mms)
{
    int i, port, err, offset = 0;
    char tcpname[256], path[256], host[128];
    char stream_selection[10 * MAX_STREAMS];

    if (mms->mms_hd) {
        url_close(mms->mms_hd);
    }
    ff_url_split(NULL, 0, NULL, 0,
            host, sizeof(host), &port, path, sizeof(path), mms->location);
    if(port<0)
        port = 80; // defaut mmsh protocol port

    ff_url_join(tcpname, sizeof(tcpname), "tcp", NULL, host, port, NULL);
    err = url_open(&mms->mms_hd, tcpname, URL_RDWR);
    if (err)
        return err;
    // send describe request
    snprintf (mms->out_buffer, sizeof(mms->out_buffer), mmsh_first_request, path,
             host, port, mms->request_seq++);
    err = send_pack(mms);
    if (err)
        return err;
    err = get_http_answer(mms); // TODO match with the first request
    if(err)
        return err;
    // close the socket and then reopen it for sending the second play request.
    url_close(mms->mms_hd);

    for (i = 0; i < mms->stream_num; i++) {
        err = snprintf(stream_selection + offset, sizeof(stream_selection) - offset,
                          "ffff:%d:0 ", mms->streams[i].id);
        if (err < 0)
            return err;
        offset += err;
    }

    // send paly request
    if (mms->seekable) {
        err = snprintf(mms->out_buffer, sizeof(mms->out_buffer), mmsh_seekable_request, path,
            host, port, 0, 0, 0, mms->request_seq++, 0, mms->stream_num, stream_selection);
    } else {
        err = snprintf(mms->out_buffer, sizeof(mms->out_buffer), mmsh_live_request, path,
            host, port, mms->request_seq++, mms->stream_num, stream_selection);
    }
    if (err < 0) {
        dprintf(NULL, "build play request failed!\n");
        return err;
    }
    dprintf(NULL, "out_buffer is %s", mms->out_buffer);

    //reopen the connection.
    err = url_open(&mms->mms_hd, tcpname, URL_RDWR);
    if (err)
        return err;
    err = send_pack(mms);
    if (err)
        return err;
    err = get_http_answer(mms);// TODO mathc with the second request
    if(err)
        return err;
    return 0;
}

static int mmsh_open(URLContext *h, const char *uri, int flags)
{
    MMSHContext *mms;
    int err;
    mms = h->priv_data = av_mallocz(sizeof(MMSHContext));
    if (!h->priv_data)
        return AVERROR(ENOMEM);
    mms->request_seq = h->is_streamed = 1;
    av_strlcpy(mms->location, uri, sizeof(mms->location));

    err =mmsh_open_cnx(mms);
    if (err) {
        dprintf(NULL, "Leaving mmsh open (failure: %d)\n", err);
        mmsh_close(h);
        return err;
    }
    dprintf(NULL, "Leaving mmsh open success.\n");
    return 0;
}

static int read_data(MMSHContext *mms, char *buf, int size)
{
    int read_size;
    read_size = FFMIN(size, mms->asf_data_remaining_len);
    memcpy(buf, mms->read_in_ptr, read_size);
    mms->asf_data_remaining_len -= read_size;
    mms->read_in_ptr      += read_size;
    return read_size;
}

static int handle_chunk_type(MMSHContext *mms)
{
    int res, len = 0;
    int chunk_type;
    chunk_type = get_chunk_header(mms, &len);

    if(chunk_type == CHUNK_TYPE_END) {
        if (mms->chunk_seq == 0) {
            dprintf(NULL, "The stream is end.\n");
            return -1;
        }
        // reconnect
        mms->request_seq = 1;
        if ((res = mmsh_open_cnx(mms)) !=0) {
            dprintf(NULL, "Reconnect failed!\n");
            return res;
        }
    } else if (chunk_type == CHUNK_TYPE_STREAM_CHANGE) {
        mms->is_header_parsed = 0;
        if ((res = get_http_header_data(mms)) !=0) {
            dprintf(NULL,"stream changed! get new header failed!\n");
            return res;
        }
    } else if (chunk_type == CHUNK_TYPE_DATA) {
        return read_data_packet(mms, len);
    } else {
        dprintf(NULL, "recv other type packet %d\n", chunk_type);
        return -1;
    }
    return 0;
}

static int mmsh_read(URLContext *h, uint8_t *buf, int size)
{
    int res = 0;
    MMSHContext *mms = h->priv_data;

    do{
        if (mms->asf_header_read_size < mms->asf_header_size) {
            // copy asf header into buffer
            char *pos;
            int size_to_copy;
            int remaining_size = mms->asf_header_size - mms->asf_header_read_size;
            size_to_copy = FFMIN(size, remaining_size);
            pos = mms->asf_header + mms->asf_header_read_size;
            memcpy(buf, pos, size_to_copy);
            mms->asf_header_read_size += size_to_copy;
            res = size_to_copy;
            if (mms->asf_header_read_size == mms->asf_header_size) {
                av_freep(&mms->asf_header); // which contains asf header
            }
        } else if (mms->asf_data_remaining_len){
            res =read_data(mms, buf, size);
        } else {
             // read data packet from network
            res = handle_chunk_type(mms);
            if (res == 0) {
                res = read_data(mms, buf, size);
            } else {
                dprintf(NULL, "other situation!\n");
            }
        }
    }while(!res);
    return res;
}

URLProtocol mms_protocol = {
    "mms",
    mmsh_open,
    mmsh_read,
    NULL, // write
    NULL, // seek
    mmsh_close,
};

URLProtocol mmsh_protocol = {
    "mmsh",
    mmsh_open,
    mmsh_read,
    NULL, // write
    NULL, // seek
    mmsh_close,
};

