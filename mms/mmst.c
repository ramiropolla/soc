/*
 * MMS protocol over TCP
 * Copyright (c) 2006,2007 Ryan Martell
 * Copyright (c) 2007 Bj�rn Axelsson
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

/* References
 * MMS protocol specification:
 *  [1]http://msdn.microsoft.com/en-us/library/cc234711(PROT.10).aspx
 * ASF specification. Revision 01.20.03.
 *  [2]http://msdn.microsoft.com/en-us/library/bb643323.aspx
 */

#include "mms.h"
#include "internal.h"
#include "libavutil/intreadwrite.h"
#include "libavcodec/bytestream.h"
#include "network.h"
#include "asf.h"

#define LOCAL_ADDRESS 0xc0a80081    // FIXME get and use correct local ip address.
#define LOCAL_PORT    1037          // as above.

typedef struct {
    MMSContext  *ff_ctx;
    int outgoing_packet_seq;             ///< Outgoing packet sequence number.
    char path[256];                      ///< Path of the resource being asked for.
    char host[128];                      ///< Host of the resources.
    int incoming_packet_seq;             ///< Incoming packet sequence number.
    int incoming_flags;                  ///< Incoming packet flags.
    int packet_id;                       ///< Identifier for packets in the current stream.
    unsigned int header_packet_id;       ///< default is 2.
} MMSTContext;


/** Client to server packet types. */
typedef enum {
    CS_PKT_INITIAL                  = 0x01,
    CS_PKT_PROTOCOL_SELECT          = 0x02,
    CS_PKT_MEDIA_FILE_REQUEST       = 0x05,
    CS_PKT_START_FROM_PKT_ID        = 0x07,
    CS_PKT_STREAM_PAUSE             = 0x09,
    CS_PKT_STREAM_CLOSE             = 0x0d,
    CS_PKT_MEDIA_HEADER_REQUEST     = 0x15,
    CS_PKT_TIMING_DATA_REQUEST      = 0x18,
    CS_PKT_USER_PASSWORD            = 0x1a,
    CS_PKT_KEEPALIVE                = 0x1b,
    CS_PKT_STREAM_ID_REQUEST        = 0x33,
} MMSCSPacketType;

/** Server to client packet types. */
typedef enum {
    /** Control packets. */
    /*@{*/
    SC_PKT_CLIENT_ACCEPTED          = 0x01,
    SC_PKT_PROTOCOL_ACCEPTED        = 0x02,
    SC_PKT_PROTOCOL_FAILED          = 0x03,
    SC_PKT_MEDIA_PKT_FOLLOWS        = 0x05,
    SC_PKT_MEDIA_FILE_DETAILS       = 0x06,
    SC_PKT_HEADER_REQUEST_ACCEPTED  = 0x11,
    SC_PKT_TIMING_TEST_REPLY        = 0x15,
    SC_PKT_PASSWORD_REQUIRED        = 0x1a,
    SC_PKT_KEEPALIVE                = 0x1b,
    SC_PKT_STREAM_STOPPED           = 0x1e,
    SC_PKT_STREAM_CHANGING          = 0x20,
    SC_PKT_STREAM_ID_ACCEPTED       = 0x21,
    /*@}*/

    /** Pseudo packets. */
    /*@{*/
    SC_PKT_CANCEL                   = -1,
    SC_PKT_NO_DATA                  = -2,
    /*@}*/

    /** Data packets. */
    /*@{*/
    SC_PKT_ASF_HEADER               = 0x010000,// make it bigger than 0xFF in case of
    SC_PKT_ASF_MEDIA                = 0x010001,// receiving false data packets.
    /*@}*/
} MMSSCPacketType;

/** Create MMST command packet header */
static void start_command_packet(MMSTContext *mmst_ctx, MMSCSPacketType packet_type)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    mms->write_out_ptr = mms->out_buffer;

    bytestream_put_le32(&mms->write_out_ptr, 1); // start sequence
    bytestream_put_le32(&mms->write_out_ptr, 0xb00bface);
    bytestream_put_le32(&mms->write_out_ptr, 0); // Length starts from after the protocol type bytes
    bytestream_put_le32(&mms->write_out_ptr, MKTAG('M','M','S',' '));
    bytestream_put_le32(&mms->write_out_ptr, 0);
    bytestream_put_le32(&mms->write_out_ptr, mmst_ctx->outgoing_packet_seq++);
    bytestream_put_le64(&mms->write_out_ptr, 0); // timestamp
    bytestream_put_le32(&mms->write_out_ptr, 0);
    bytestream_put_le16(&mms->write_out_ptr, packet_type);
    bytestream_put_le16(&mms->write_out_ptr, 3); // direction to server
}

/** Add prefixes to MMST command packet. */
static void insert_command_prefixes(MMSTContext *mmst_ctx,
        uint32_t prefix1, uint32_t prefix2)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    bytestream_put_le32(&mms->write_out_ptr, prefix1); // first prefix
    bytestream_put_le32(&mms->write_out_ptr, prefix2); // second prefix
}

/** Send a prepared MMST command packet. */
static int send_command_packet(MMSTContext *mmst_ctx)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    int len= mms->write_out_ptr - mms->out_buffer;
    int exact_length = (len + 7) & ~7;
    int first_length= exact_length - 16;
    int len8= first_length/8;
    int write_result;

    // update packet length fields.
    AV_WL32(mms->out_buffer + 8, first_length);
    AV_WL32(mms->out_buffer + 16, len8);
    AV_WL32(mms->out_buffer + 32, len8-2);
    memset(mms->write_out_ptr, 0, exact_length - len);

    // write it out.
    write_result= url_write(mms->mms_hd, mms->out_buffer, exact_length);
    if(write_result != exact_length) {
        av_log(NULL, AV_LOG_ERROR,
               "Failed to write data of length %d: %d (%s)\n",
               exact_length, write_result,
               write_result < 0 ? strerror(write_result) :
                   "The server closed the connection");
        return AVERROR_IO;
    }

    return 0;
}

static void mms_put_utf16(MMSTContext *mmst_ctx, uint8_t *src)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    ByteIOContext bic;
    int size = mms->write_out_ptr - mms->out_buffer;
    int len;
    init_put_byte(&bic, mms->write_out_ptr,
            sizeof(mms->out_buffer) - size, 1, NULL, NULL, NULL, NULL);

    len = ff_put_str16_nolen(&bic, src);
    mms->write_out_ptr += len;
}

static int send_time_test_data(MMSTContext *mmst_ctx)
{
    start_command_packet(mmst_ctx, CS_PKT_TIMING_DATA_REQUEST);
    insert_command_prefixes(mmst_ctx, 0xf0f0f0f1, 0x0004000b);
    return send_command_packet(mmst_ctx);
}

static int send_protocol_select(MMSTContext *mmst_ctx)
{
    char data_string[256];
    MMSContext *mms = mmst_ctx->ff_ctx;

    start_command_packet(mmst_ctx, CS_PKT_PROTOCOL_SELECT);
    insert_command_prefixes(mmst_ctx, 0, 0xffffffff);
    bytestream_put_le32(&mms->write_out_ptr, 0);          // maxFunnelBytes
    bytestream_put_le32(&mms->write_out_ptr, 0x00989680); // maxbitRate
    bytestream_put_le32(&mms->write_out_ptr, 2);          // funnelMode
    snprintf(data_string, sizeof(data_string), "\\\\%d.%d.%d.%d\\%s\\%d",
            (LOCAL_ADDRESS>>24)&0xff,
            (LOCAL_ADDRESS>>16)&0xff,
            (LOCAL_ADDRESS>>8)&0xff,
            LOCAL_ADDRESS&0xff,
            "TCP",                                        // or UDP
            LOCAL_PORT);

    mms_put_utf16(mmst_ctx, data_string);
    return send_command_packet(mmst_ctx);
}

static int send_media_file_request(MMSTContext *mmst_ctx)
{
    MMSContext *mms = mmst_ctx->ff_ctx;

    start_command_packet(mmst_ctx, CS_PKT_MEDIA_FILE_REQUEST);
    insert_command_prefixes(mmst_ctx, 1, 0xffffffff);
    bytestream_put_le32(&mms->write_out_ptr, 0);
    bytestream_put_le32(&mms->write_out_ptr, 0);
    mms_put_utf16(mmst_ctx, mmst_ctx->path + 1); // +1 for skip "/"

    return send_command_packet(mmst_ctx);
}

static void handle_packet_stream_changing_type(MMSTContext *mmst_ctx)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    dprintf(NULL, "Stream changing!\n");

    // 40 is the packet header size, 7 is the prefix size.
    mmst_ctx->header_packet_id= AV_RL32(mms->in_buffer + 40 + 7);
    dprintf(NULL, "Changed header prefix to 0x%x", mmst_ctx->header_packet_id);
}

static int send_keepalive_packet(MMSTContext *mmst_ctx)
{
    // respond to a keepalive with a keepalive...
    start_command_packet(mmst_ctx, CS_PKT_KEEPALIVE);
    insert_command_prefixes(mmst_ctx, 1, 0x100FFFF);
    return send_command_packet(mmst_ctx);
}

/** Pad media packets smaller than max_packet_size and/or adjust read position
  * after a seek. */
static void pad_media_packet(MMSTContext *mmst_ctx)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    if(mms->remaining_in_len<mms->asf_packet_len) {
        int padding_size = mms->asf_packet_len - mms->remaining_in_len;
        memset(mms->in_buffer + mms->remaining_in_len, 0, padding_size);
        mms->remaining_in_len += padding_size;
    }
}

/** Read incoming MMST media, header or command packet. */
static MMSSCPacketType get_tcp_server_response(MMSTContext *mmst_ctx)
{
    int read_result;
    MMSSCPacketType packet_type= -1;
    MMSContext *mms = mmst_ctx->ff_ctx;
    for(;;) {
        read_result = url_read_complete(mms->mms_hd, mms->in_buffer, 8);
        if (read_result != 8) {
            if(read_result < 0) {
                av_log(NULL, AV_LOG_ERROR,
                       "Error reading packet header: %d (%s)\n",
                       read_result, strerror(read_result));
                packet_type = SC_PKT_CANCEL;
            } else {
                av_log(NULL, AV_LOG_ERROR,
                       "The server closed the connection\n");
                packet_type = SC_PKT_NO_DATA;
            }
            return packet_type;
        }

        // handle command packet.
        if(AV_RL32(mms->in_buffer + 4)==0xb00bface) {
            int length_remaining, hr;

            mmst_ctx->incoming_flags= mms->in_buffer[3];
            read_result= url_read_complete(mms->mms_hd, mms->in_buffer+8, 4);
            if(read_result != 4) {
                av_log(NULL, AV_LOG_ERROR,
                       "Reading command packet length failed: %d (%s)\n",
                       read_result,
                       read_result < 0 ? strerror(read_result) :
                           "The server closed the connection");
                return read_result < 0 ? read_result : AVERROR_IO;
            }

            length_remaining= AV_RL32(mms->in_buffer+8) + 4;
            dprintf(NULL, "Length remaining is %d\n", length_remaining);
            // read the rest of the packet.
            if (length_remaining < 0
                || length_remaining > sizeof(mms->in_buffer) - 12) {
                av_log(NULL, AV_LOG_ERROR,
                       "Incoming packet length %d exceeds bufsize %zu\n",
                       length_remaining, sizeof(mms->in_buffer) - 12);
                return AVERROR_INVALIDDATA;
            }
            read_result = url_read_complete(mms->mms_hd, mms->in_buffer + 12,
                                            length_remaining) ;
            if (read_result != length_remaining) {
                av_log(NULL, AV_LOG_ERROR,
                       "Reading pkt data (length=%d) failed: %d (%s)\n",
                       length_remaining, read_result,
                       read_result < 0 ? strerror(read_result) :
                           "The server closed the connection");
                return read_result < 0 ? read_result : AVERROR_IO;
            }
            packet_type= AV_RL16(mms->in_buffer+36);
            hr = AV_RL32(mms->in_buffer + 40);
            if (hr) {
                av_log(NULL, AV_LOG_ERROR,
                       "Server sent an error status code: 0x%08x\n", hr);
                return AVERROR_UNKNOWN;
            }
        } else {
            int length_remaining;
            int packet_id_type;
            int tmp;

            // note we cache the first 8 bytes,
            // then fill up the buffer with the others
            tmp                            = AV_RL16(mms->in_buffer + 6);
            length_remaining               = (tmp - 8) & 0xffff;
            mmst_ctx->incoming_packet_seq  = AV_RL32(mms->in_buffer);
            packet_id_type                 = mms->in_buffer[4];
            mmst_ctx->incoming_flags       = mms->in_buffer[5];

            if (length_remaining < 0
                || length_remaining > sizeof(mms->in_buffer) - 8) {
                av_log(NULL, AV_LOG_ERROR,
                       "Data length %d is invalid or too large (max=%zu)\n",
                       length_remaining, sizeof(mms->in_buffer));
                return AVERROR_INVALIDDATA;
            }
            mms->remaining_in_len    = length_remaining;
            mms->read_in_ptr         = mms->in_buffer;
            read_result= url_read_complete(mms->mms_hd, mms->in_buffer, length_remaining);
            if(read_result != length_remaining) {
                av_log(NULL, AV_LOG_ERROR,
                       "Failed to read packet data of size %d: %d (%s)\n",
                       length_remaining, read_result,
                       read_result < 0 ? strerror(read_result) :
                           "The server closed the connection");
                return read_result < 0 ? read_result : AVERROR_IO;
            }

            // if we successfully read everything.
            if(packet_id_type == mmst_ctx->header_packet_id) {
                packet_type = SC_PKT_ASF_HEADER;
                // Store the asf header
                if(!mms->header_parsed) {
                    void *p = av_realloc(mms->asf_header,
                                  mms->asf_header_size + mms->remaining_in_len);
                    if (!p) {
                        av_freep(&mms->asf_header);
                        return AVERROR(ENOMEM);
                    }
                    mms->asf_header = p;
                    memcpy(mms->asf_header + mms->asf_header_size,
                           mms->read_in_ptr, mms->remaining_in_len);
                    mms->asf_header_size += mms->remaining_in_len;
                }
                // 0x04 means asf header is sent in multiple packets.
                if (mmst_ctx->incoming_flags == 0x04)
                    continue;
            } else if(packet_id_type == mmst_ctx->packet_id) {
                packet_type = SC_PKT_ASF_MEDIA;
            } else {
                dprintf(NULL, "packet id type %d is old.", packet_id_type);
                continue;
            }
        }

        // preprocess some packet type
        if(packet_type == SC_PKT_KEEPALIVE) {
            send_keepalive_packet(mmst_ctx);
            continue;
        } else if(packet_type == SC_PKT_STREAM_CHANGING) {
            handle_packet_stream_changing_type(mmst_ctx);
        } else if(packet_type == SC_PKT_ASF_MEDIA) {
            pad_media_packet(mmst_ctx);
        }
        return packet_type;
    }
}

static int mms_safe_send_recv(MMSTContext *mmst_ctx,
                              int (*send_fun)(MMSTContext *mmst_ctx),
                              const MMSSCPacketType expect_type)
{
    MMSSCPacketType type;
    if(send_fun) {
        int ret = send_fun(mmst_ctx);
        if (ret < 0) {
            dprintf(NULL, "Send Packet error before expecting recv packet %d\n", expect_type);
            return ret;
        }
    }

    if ((type = get_tcp_server_response(mmst_ctx)) != expect_type) {
        av_log(NULL, AV_LOG_ERROR,
               "Corrupt stream (unexpected packet type 0x%x, expected 0x%x)\n",
               type, expect_type);
        return AVERROR_INVALIDDATA;
    } else {
        return 0;
    }
}

static int send_media_header_request(MMSTContext *mmst_ctx)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    start_command_packet(mmst_ctx, CS_PKT_MEDIA_HEADER_REQUEST);
    insert_command_prefixes(mmst_ctx, 1, 0);
    bytestream_put_le32(&mms->write_out_ptr, 0);
    bytestream_put_le32(&mms->write_out_ptr, 0x00800000);
    bytestream_put_le32(&mms->write_out_ptr, 0xffffffff);
    bytestream_put_le32(&mms->write_out_ptr, 0);
    bytestream_put_le32(&mms->write_out_ptr, 0);
    bytestream_put_le32(&mms->write_out_ptr, 0);

    // the media preroll value in milliseconds?
    bytestream_put_le32(&mms->write_out_ptr, 0);
    bytestream_put_le32(&mms->write_out_ptr, 0x40AC2000);
    bytestream_put_le32(&mms->write_out_ptr, 2);
    bytestream_put_le32(&mms->write_out_ptr, 0);

    return send_command_packet(mmst_ctx);
}

/** Send the initial handshake. */
static int send_startup_packet(MMSTContext *mmst_ctx)
{
    char data_string[256];
    // SubscriberName is defined in MS specification linked below.
    // The guid value can be any valid value.
    // http://download.microsoft.com/
    // download/9/5/E/95EF66AF-9026-4BB0-A41D-A4F81802D92C/%5BMS-WMSP%5D.pdf
    snprintf(data_string, sizeof(data_string),
            "NSPlayer/7.0.0.1956; {%s}; Host: %s",
            "7E667F5D-A661-495E-A512-F55686DDA178", mmst_ctx->host);

    start_command_packet(mmst_ctx, CS_PKT_INITIAL);
    insert_command_prefixes(mmst_ctx, 0, 0x0004000b);
    bytestream_put_le32(&(mmst_ctx->ff_ctx->write_out_ptr), 0x0003001c);
    mms_put_utf16(mmst_ctx, data_string);
    return send_command_packet(mmst_ctx);
}

/** Send MMST stream selection command based on the AVStream->discard values. */
static int send_stream_selection_request(MMSTContext *mmst_ctx)
{
    int i;
    MMSContext *mms = mmst_ctx->ff_ctx;
    //  send the streams we want back...
    start_command_packet(mmst_ctx, CS_PKT_STREAM_ID_REQUEST);
    bytestream_put_le32(&mms->write_out_ptr, mms->stream_num);         // stream nums
    for(i= 0; i<mms->stream_num; i++) {
        bytestream_put_le16(&mms->write_out_ptr, 0xffff);              // flags
        bytestream_put_le16(&mms->write_out_ptr, mms->streams[i].id);  // stream id
        bytestream_put_le16(&mms->write_out_ptr, 0);                   // selection
    }
    return send_command_packet(mmst_ctx);
}

static int send_close_packet(MMSTContext *mmst_ctx)
{
    start_command_packet(mmst_ctx, CS_PKT_STREAM_CLOSE);
    insert_command_prefixes(mmst_ctx, 1, 1);

    return send_command_packet(mmst_ctx);
}

/** Close the MMSH/MMST connection */
static int mms_close(URLContext *h)
{
    MMSTContext *mmst_ctx = (MMSTContext *)h->priv_data;
    MMSContext *mms = mmst_ctx->ff_ctx;
    if(mms->mms_hd) {
        send_close_packet(mmst_ctx);
        url_close(mms->mms_hd);
    }

    /* free all separately allocated pointers in mms */
    av_free(mms->asf_header);
    av_free(mmst_ctx->ff_ctx);
    av_freep(&h->priv_data);

    return 0;
}

static int send_media_packet_request(MMSTContext *mmst_ctx)
{
    MMSContext *mms = mmst_ctx->ff_ctx;
    start_command_packet(mmst_ctx, CS_PKT_START_FROM_PKT_ID);
    insert_command_prefixes(mmst_ctx, 1, 0x0001FFFF);
    bytestream_put_le64(&mms->write_out_ptr, 0);          // seek timestamp
    bytestream_put_le32(&mms->write_out_ptr, 0xffffffff); // unknown
    bytestream_put_le32(&mms->write_out_ptr, 0xffffffff); // packet offset
    bytestream_put_byte(&mms->write_out_ptr, 0xff);       // max stream time limit
    bytestream_put_byte(&mms->write_out_ptr, 0xff);       // max stream time limit
    bytestream_put_byte(&mms->write_out_ptr, 0xff);       // max stream time limit
    bytestream_put_byte(&mms->write_out_ptr, 0x00);       // stream time limit flag

    mmst_ctx->packet_id++;                                     // new packet_id
    bytestream_put_le32(&mms->write_out_ptr, mmst_ctx->packet_id);
    return send_command_packet(mmst_ctx);
}

static void clear_stream_buffers(MMSTContext *mmst_ctx)
{
    MMSContext *mms       = mmst_ctx->ff_ctx;
    mms->remaining_in_len = 0;
    mms->read_in_ptr      = mms->in_buffer;
}

static int mms_open(URLContext *h, const char *uri, int flags)
{
    MMSTContext *mmst_ctx;
    MMSContext *mms;
    int port, err;
    char tcpname[256];

    h->is_streamed = 1;
    mmst_ctx = h->priv_data = av_mallocz(sizeof(MMSTContext));
    if (!h->priv_data)
        return AVERROR(ENOMEM);
    mms = mmst_ctx->ff_ctx = av_mallocz(sizeof(MMSContext));
    if (!mmst_ctx->ff_ctx )
        return AVERROR(ENOMEM);

    // only for MMS over TCP, so set proto = NULL
    av_url_split(NULL, 0, NULL, 0,
            mmst_ctx->host, sizeof(mmst_ctx->host), &port, mmst_ctx->path,
            sizeof(mmst_ctx->path), uri);

    if(port<0)
        port = 1755; // defaut mms protocol port

    // establish tcp connection.
    ff_url_join(tcpname, sizeof(tcpname), "tcp", NULL, mmst_ctx->host, port, NULL);
    err = url_open(&mms->mms_hd, tcpname, URL_RDWR);
    if (err)
        goto fail;

    mmst_ctx->packet_id        = 3;          // default, initial value.
    mmst_ctx->header_packet_id = 2;          // default, initial value.
    err = mms_safe_send_recv(mmst_ctx, send_startup_packet, SC_PKT_CLIENT_ACCEPTED);
    if (err)
        goto fail;
    err = mms_safe_send_recv(mmst_ctx, send_time_test_data, SC_PKT_TIMING_TEST_REPLY);
    if (err)
        goto fail;
    err = mms_safe_send_recv(mmst_ctx, send_protocol_select, SC_PKT_PROTOCOL_ACCEPTED);
    if (err)
        goto fail;
    err = mms_safe_send_recv(mmst_ctx, send_media_file_request, SC_PKT_MEDIA_FILE_DETAILS);
    if (err)
        goto fail;
    err = mms_safe_send_recv(mmst_ctx, send_media_header_request, SC_PKT_HEADER_REQUEST_ACCEPTED);
    if (err)
        goto fail;
    err = mms_safe_send_recv(mmst_ctx, NULL, SC_PKT_ASF_HEADER);
    if (err)
        goto fail;
    if((mmst_ctx->incoming_flags != 0X08) && (mmst_ctx->incoming_flags != 0X0C))
        goto fail;
    err = ff_mms_asf_header_parser(mms);
    if (err) {
        dprintf(NULL, "asf header parsed failed!\n");
        goto fail;
    }
    mms->header_parsed = 1;

    if (!mms->asf_packet_len || !mms->stream_num)
        goto fail;

    /* Since we read the header at open(), this shouldn't be possible */
    assert(mmst_ctx->ff_ctx->header_parsed);

    dprintf(NULL, "mms_read() before play().\n");
    clear_stream_buffers(mmst_ctx);
    err = mms_safe_send_recv(mmst_ctx, send_stream_selection_request, SC_PKT_STREAM_ID_ACCEPTED);
    if (err)
        goto fail;
    // send media packet request
    err = mms_safe_send_recv(mmst_ctx, send_media_packet_request, SC_PKT_MEDIA_PKT_FOLLOWS);
    if (err) {
        goto fail;
    }

    dprintf(NULL, "Leaving open (success)\n");
    return 0;
fail:
    mms_close(h);
    dprintf(NULL, "Leaving open (failure: %d)\n", err);
    return err;
}

/** Read ASF data through the protocol. */
static int mms_read(URLContext *h, uint8_t *buf, int size)
{
    /* TODO: see tcp.c:tcp_read() about a possible timeout scheme */
    MMSTContext *mmst_ctx = h->priv_data;
    int result            = 0;
    MMSContext *mms       = mmst_ctx->ff_ctx;
    do {
        if(mms->asf_header_read_size < mms->asf_header_size) {
           result =  ff_mms_read_header(mms, buf, size);
        } else if(mms->remaining_in_len) {
            /* Read remaining packet data to buffer.
             * the result can not be zero because remaining_in_len is positive.*/
            result = ff_mms_read_data(mms, buf, size);
        } else {
            /* Read from network */
            int err = mms_safe_send_recv(mmst_ctx, NULL, SC_PKT_ASF_MEDIA);
            if (err == 0) {
                if(mms->remaining_in_len>mms->asf_packet_len) {
                    av_log(NULL, AV_LOG_ERROR,
                           "Incoming pktlen %d is larger than ASF pktsize %d\n",
                           mms->remaining_in_len, mms->asf_packet_len);
                    result= AVERROR_IO;
                } else {
                    // copy the data to the packet buffer.
                    result = ff_mms_read_data(mms, buf, size);
                    if (result == 0) {
                        dprintf(NULL, "read asf media paket size is zero!\n");
                        break;
                    }
                }
            } else {
                dprintf(NULL, "read packet error!\n");
                break;
            }
        }
    } while(!result); // only return one packet.
    return result;
}

URLProtocol mmst_protocol = {
    "mmst",
    mms_open,
    mms_read,
    NULL, // write
    NULL, // seek
    mms_close,
};
