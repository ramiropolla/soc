/*
 * RTMP network protocol
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

/**
 * @file libavformat/rtmpproto.c
 * RTMP protocol
 */

#include "libavutil/avstring.h"
#include "avformat.h"

#include <unistd.h>
#include <stdarg.h>
#include "network.h"
#include "os_support.h"
#include <fcntl.h>
#if HAVE_SYS_SELECT_H
#include <sys/select.h>
#endif

#define RTMP_DEFAULT_PORT 1935
#define RTMP_HANDSHAKE_PACKET_SIZE 1536

typedef struct RTMPContext {
    URLContext *rtmp_hd;
    int wrote;
} RTMPContext;


/**
 * url syntax: rtp://host:port[?option=val...]
 * option: 'ttl=n'       : set the ttl value (for multicast only)
 *         'localport=n' : set the local port to n
 *
 */

static int rtmp_open(URLContext *h, const char *uri, int flags)
{
    RTMPContext *s;
    int port, is_input;
    char hostname[256];
    char buf[1024];
    char path[1024];

    is_input = !(flags & URL_WRONLY);

    s = av_mallocz(sizeof(RTMPContext));
    if (!s)
        return AVERROR(ENOMEM);
    h->priv_data = s;

    url_split(NULL, 0, NULL, 0, hostname, sizeof(hostname), &port,
              path, sizeof(path), uri);

    if (port == -1)
        port = RTMP_DEFAULT_PORT;
    snprintf(buf, sizeof(buf), "tcp://%s:%d", hostname, port);

    if (url_open(&s->rtmp_hd, buf, URL_RDWR) < 0)
        goto fail;

    h->max_packet_size = url_get_max_packet_size(s->rtmp_hd);
    h->is_streamed = 1;
    url_close(s->rtmp_hd);
    return 0;

 fail:
    if (s->rtmp_hd)
        url_close(s->rtmp_hd);
    av_free(s);
    return AVERROR(EIO);
}

static int rtmp_read(URLContext *h, uint8_t *buf, int size)
{
    return 0;
}

static int rtmp_write(URLContext *h, uint8_t *buf, int size)
{
    return 0;
}

static int rtmp_close(URLContext *h)
{
    RTMPContext *s = h->priv_data;

    av_free(s);
    return 0;
}

URLProtocol rtmp_protocol = {
    "rtmp",
    rtmp_open,
    rtmp_read,
    rtmp_write,
    NULL, /* seek */
    rtmp_close,
};
