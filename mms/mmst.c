/*
 * MMS protocol over TCP
 * Copyright (c) 2006,2007 Ryan Martell
 * Copyright (c) 2007 Björn Axelsson
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
#include "libavutil/intreadwrite.h"
#include "libavutil/avstring.h"
#include "network.h"
#include "asf.h"

/** Client to server packet types. */
typedef enum {
    CS_PKT_INITIAL= 0x01,
    CS_PKT_PROTOCOL_SELECT= 0x02,
    CS_PKT_MEDIA_FILE_REQUEST= 0x05,
    CS_PKT_START_FROM_PKT_ID= 0x07,
    CS_PKT_STREAM_PAUSE= 0x09,
    CS_PKT_STREAM_CLOSE= 0x0d,
    CS_PKT_MEDIA_HEADER_REQUEST= 0x15,
    CS_PKT_TIMING_DATA_REQUEST= 0x18,
    CS_PKT_USER_PASSWORD= 0x1a,
    CS_PKT_KEEPALIVE= 0x1b,
    CS_PKT_STREAM_ID_REQUEST= 0x33,
} MMSCSPacketType;

/** Server to client packet types. */
typedef enum {
    /** Control packets. */
    /*@{*/
    SC_PKT_CLIENT_ACCEPTED= 0x01,
    SC_PKT_PROTOCOL_ACCEPTED= 0x02,
    SC_PKT_PROTOCOL_FAILED= 0x03,
    SC_PKT_MEDIA_PKT_FOLLOWS= 0x05,
    SC_PKT_MEDIA_FILE_DETAILS= 0x06,
    SC_PKT_HEADER_REQUEST_ACCEPTED= 0x11,
    SC_PKT_TIMING_TEST_REPLY= 0x15,
    SC_PKT_PASSWORD_REQUIRED= 0x1a,
    SC_PKT_KEEPALIVE= 0x1b,
    SC_PKT_STREAM_STOPPED= 0x1e,
    SC_PKT_STREAM_CHANGING= 0x20,
    SC_PKT_STREAM_ID_ACCEPTED= 0x21,
    /*@}*/

    /** Pseudo packets. */
    /*@{*/
    SC_PKT_CANCEL = -1,
    SC_PKT_NO_DATA = -2,
    SC_PKT_HTTP_CONTROL_ACKNOWLEDGE = -3,
    /*@}*/

    /** Data packets. */
    /*@{*/
    SC_PKT_ASF_HEADER= 0x81,
    SC_PKT_ASF_MEDIA= 0x82,
    /*@}*/
} MMSSCPacketType;

typedef struct {
    uint32_t local_ip_address;
    int local_port;                      ///< My local port (sent but not correct).
    int sequence_number;                 ///< Outgoing packet sequence number.
    char path[256];                      ///< Path of the resource being asked for.
    char host[128];                      ///< Host of the resources.
    int port;                            ///< Port of the resource.

    URLContext *mms_hd;                  ///< TCP connection handle

    /** Buffer for outgoing packets. */
    /*@{*/
    ByteIOContext outgoing_packet_data;  ///< Outgoing packet stream
    uint8_t outgoing_packet_buffer[512]; ///< Outgoing packet data
    /*@}*/

    /** Buffer for incoming control packets. */
    /*@{*/
    uint8_t incoming_buffer[8192];       ///< Incoming buffer location.
    int incoming_buffer_length;          ///< Incoming buffer length.
    /*@}*/

    /** Buffer for incoming media/header packets. */
    /*@{*/
    uint8_t pkt_buf[8192];               ///< header or media packet.
    uint8_t *pkt_read_ptr;               ///< Pointer for partial reads.
    int pkt_buf_len;                     ///< Buffer length.
    int pkt_offset;                      ///< offset in packet.
    /*@}*/

    int incoming_packet_seq;             ///< Incoming packet sequence number.
    int incoming_flags;                  ///< Incoming packet flags.

    int packet_id;                       ///< Identifier for packets in the current stream.
    unsigned int header_packet_id;       ///< default is 2.

    /** Internal handling of the ASF header */
    /*@{*/
    uint8_t *asf_header;                 ///< Stored ASF header.
    int asf_header_size;                 ///< Size of stored ASF header.
    int asf_header_read_pos;             ///< Current read position in header.
    int header_parsed;                   ///< The header has been received and parsed.
    int asf_packet_len;
    /*@}*/

    char location[4096];
    int stream_num;
    int streaming_flag;
} MMSContext;

/** Close the remote connection. */
static void close_connection(MMSContext *mms)
{
    url_close(mms->mms_hd);
}

/** Create MMST command packet header */
static void start_command_packet(MMSContext *mms, MMSCSPacketType packet_type)
{
    ByteIOContext *context= &mms->outgoing_packet_data;

    url_fseek(context, 0, SEEK_SET);
    put_le32(context, 1); // start sequence
    put_le32(context, 0xb00bface);
    put_le32(context, 0); // Length starts from after the protocol type bytes
    put_le32(context, MKTAG('M','M','S',' '));
    put_le32(context, 0);
    put_le32(context, mms->sequence_number++);
    put_le64(context, 0); // timestmamp
    put_le32(context, 0);
    put_le16(context, packet_type);
    put_le16(context, 3); // direction to server
}

/** Add prefixes to MMST command packet. */
static void insert_command_prefixes(MMSContext *mms,
        uint32_t prefix1, uint32_t prefix2)
{
    ByteIOContext *context= &mms->outgoing_packet_data;

    put_le32(context, prefix1); // first prefix
    put_le32(context, prefix2); // second prefix
}

static void put_le_utf16(ByteIOContext *pb, char *utf8)
{
    int val;

    while(*utf8) {
        GET_UTF8(val, *utf8++, break;);
        put_le16(pb, val);
    }

    put_le16(pb, 0x00);

    return;
}

/** Send a prepared MMST command packet. */
static int send_command_packet(MMSContext *mms)
{
    ByteIOContext *context= &mms->outgoing_packet_data;
    int exact_length= url_ftell(context);
    int first_length= exact_length - 16;
    int len8= first_length/8;
    int write_result;

    // update packet length fields.
    url_fseek(context, 8, SEEK_SET);
    put_le32(context, first_length);
    url_fseek(context, 16, SEEK_SET);
    put_le32(context, len8);
    url_fseek(context, 32, SEEK_SET);
    put_le32(context, len8-2);

    // seek back to the end.
    url_fseek(context, exact_length, SEEK_SET);

    // write it out.
    write_result= url_write(mms->mms_hd, context->buffer, exact_length);
    if(write_result != exact_length) {
        dprintf(NULL, "url_write returned: %d != %d\n",
                write_result, exact_length);
        return AVERROR_IO;
    }

    return 0;
}

static int send_protocol_select(MMSContext *mms)
{
    char data_string[256];

    start_command_packet(mms, CS_PKT_PROTOCOL_SELECT);
    insert_command_prefixes(mms, 0, 0xffffffff);
    put_le32(&mms->outgoing_packet_data, 0);          // maxFunnelBytes
    put_le32(&mms->outgoing_packet_data, 0x00989680); // maxbitRate
    put_le32(&mms->outgoing_packet_data, 2);          // funnelMode
    snprintf(data_string, sizeof(data_string), "\\\\%d.%d.%d.%d\\%s\\%d",
            (mms->local_ip_address>>24)&0xff,
            (mms->local_ip_address>>16)&0xff,
            (mms->local_ip_address>>8)&0xff,
            mms->local_ip_address&0xff,
            "TCP",                                    // or UDP
            mms->local_port);
    put_le_utf16(&mms->outgoing_packet_data, data_string);

    return send_command_packet(mms);
}

static int send_media_file_request(MMSContext *mms)
{
    start_command_packet(mms, CS_PKT_MEDIA_FILE_REQUEST);
    insert_command_prefixes(mms, 1, 0xffffffff);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le_utf16(&mms->outgoing_packet_data, mms->path+1); // +1 for skip "/".
    put_le32(&mms->outgoing_packet_data, 0);

    return send_command_packet(mms);
}

static int read_bytes(MMSContext *mms, uint8_t *buffer, int length_to_read)
{
    int len= 0;

    while(len<length_to_read) {
        int read_result= url_read(mms->mms_hd, buffer+len, length_to_read-len);
        if(read_result < 0)
            return read_result;
        if(read_result) {
            len+= read_result;
        } else
            return read_result;
    }

    return len;
}

static void handle_packet_stream_changing_type(MMSContext *mms)
{
    ByteIOContext pkt;
    dprintf(NULL, "Stream changing!\n");

    // 40 is the packet header size, without the prefixea.s
    init_put_byte(&pkt, mms->incoming_buffer+40,
            mms->incoming_buffer_length-40, 0, NULL, NULL, NULL, NULL);
    get_le32(&pkt);                                 // prefix 1
    mms->header_packet_id= (get_le32(&pkt) & 0xff); // prefix 2
    dprintf(NULL, "Changed header prefix to 0x%x", mms->header_packet_id);
}

static int send_keepalive_packet(MMSContext *mms)
{
    // respond to a keepalive with a keepalive...
    start_command_packet(mms, CS_PKT_KEEPALIVE);
    insert_command_prefixes(mms, 1, 0x100FFFF);
    return send_command_packet(mms);
}

/** Pad media packets smaller than max_packet_size and/or adjust read position
  * after a seek. */
static void pad_media_packet(MMSContext *mms)
{
    if(mms->pkt_buf_len<mms->asf_packet_len) {
        int padding_size = mms->asf_packet_len - mms->pkt_buf_len;
        memset(mms->pkt_buf + mms->pkt_buf_len, 0, padding_size);
        mms->pkt_buf_len += padding_size;
    }
    if(mms->pkt_offset) {
        mms->pkt_buf_len -= mms->pkt_offset;
        mms->pkt_read_ptr += mms->pkt_offset;
        mms->pkt_offset = 0;
    }
}

/** Read incoming MMST media, header or command packet. */
static MMSSCPacketType get_tcp_server_response(MMSContext *mms)
{
    int read_result;
    MMSSCPacketType packet_type= -1;
    int done;

    do {
        done= 1;
        if((read_result= read_bytes(mms, mms->incoming_buffer, 8))==8) {
            // handle command packet.
            if(AV_RL32(mms->incoming_buffer + 4)==0xb00bface) {
                mms->incoming_flags= mms->incoming_buffer[3];
                read_result= read_bytes(mms, mms->incoming_buffer+8, 4);
                if(read_result == 4) {
                    int length_remaining= AV_RL32(mms->incoming_buffer+8) + 4;

                    dprintf(NULL, "Length remaining is %d\n", length_remaining);
                    // read the rest of the packet.
                    read_result = read_bytes(mms, mms->incoming_buffer + 12,
                                                  length_remaining) ;
                    if (read_result == length_remaining) {
                        mms->incoming_buffer_length= length_remaining+12;
                        packet_type= AV_RL16(mms->incoming_buffer+36);

                    } else {
                        dprintf(NULL, "3 read returned %d!\n", read_result);
                    }
                } else {
                    dprintf(NULL, "2 read returned %d!\n", read_result);
                }
            } else {
                int length_remaining;
                int packet_id_type;
                int tmp;

                assert(mms->pkt_buf_len==0);

                //** VERIFY LENGTH REMAINING HAS SPACE
                // note we cache the first 8 bytes,
                // then fill up the buffer with the others
                tmp                       = AV_RL16(mms->incoming_buffer + 6);
                length_remaining          = (tmp - 8) & 0xffff;
                mms->incoming_packet_seq  = AV_RL32(mms->incoming_buffer);
                packet_id_type            = mms->incoming_buffer[4];
                mms->incoming_flags       = mms->incoming_buffer[5];
                mms->pkt_buf_len          = length_remaining;
                mms->pkt_read_ptr         = mms->pkt_buf;

                if(mms->pkt_buf_len >= sizeof(mms->pkt_buf)) {
                    dprintf(NULL, "Incoming Buffer Length overflow: %d>%d\n",
                    mms ->pkt_buf_len, (int) sizeof(mms->pkt_buf));
                }
                read_result= read_bytes(mms, mms->pkt_buf, length_remaining);
                if(read_result != length_remaining) {
                    dprintf(NULL, "read_bytes result: %d asking for %d\n",
                            read_result, length_remaining);
                    break;
                } else {
                    // if we successfully read everything.
                    if(packet_id_type == mms->header_packet_id) {
                        packet_type = SC_PKT_ASF_HEADER;
                        // Store the asf header
                        if(!mms->header_parsed) {
                            mms->asf_header = av_realloc(mms->asf_header,
                                              mms->asf_header_size
                                              + mms->pkt_buf_len);
                            memcpy(mms->asf_header + mms->asf_header_size,
                                                 mms->pkt_read_ptr,
                                                 mms->pkt_buf_len);
                            mms->asf_header_size += mms->pkt_buf_len;
                        }
                    } else if(packet_id_type == mms->packet_id) {
                        packet_type = SC_PKT_ASF_MEDIA;
                    } else {
                        dprintf(NULL, "packet id type %d is old.", packet_id_type);
                        done= 0;
                    }
                }
            }
        } else {
            if(read_result<0) {
                dprintf(NULL, "Read error (or cancelled) returned %d!\n", read_result);
                packet_type = SC_PKT_CANCEL;
            } else {
                dprintf(NULL, "Read result of zero?!\n");
                packet_type = SC_PKT_NO_DATA;
            }
            done = 1;
        }
    } while(!done);

    if (packet_type == SC_PKT_KEEPALIVE) {
        send_keepalive_packet(mms);
    } else if (packet_type == SC_PKT_STREAM_CHANGING) {
        handle_packet_stream_changing_type(mms);
        //TODO: Handle new header when change the stream type.
    } else if (packet_type == SC_PKT_ASF_MEDIA) {
        pad_media_packet(mms);
    }
    return packet_type;
}

static int send_media_header_request(MMSContext *mms)
{
    start_command_packet(mms, CS_PKT_MEDIA_HEADER_REQUEST);
    insert_command_prefixes(mms, 1, 0);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le32(&mms->outgoing_packet_data, 0x00800000);
    put_le32(&mms->outgoing_packet_data, 0xffffffff);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le32(&mms->outgoing_packet_data, 0);

    // the media preroll value in milliseconds?
    put_le32(&mms->outgoing_packet_data, 0);
    put_le32(&mms->outgoing_packet_data, 0x40AC2000);
    put_le32(&mms->outgoing_packet_data, 2);
    put_le32(&mms->outgoing_packet_data, 0);

    return send_command_packet(mms);
}

/** Send the initial handshake. */
static int send_startup_packet(MMSContext *mms)
{
    char data_string[256];
    // SubscriberName is defined in MS specification linked below.
    // The guid value can be any valid value.
    // http://download.microsoft.com/
    // download/9/5/E/95EF66AF-9026-4BB0-A41D-A4F81802D92C/%5BMS-WMSP%5D.pdf
    snprintf(data_string, sizeof(data_string),
            "NSPlayer/7.0.0.1956; {%s}; Host: %s",
            "7E667F5D-A661-495E-A512-F55686DDA178", mms->host);

    start_command_packet(mms, CS_PKT_INITIAL);
    insert_command_prefixes(mms, 0, 0x0004000b);
    put_le32(&mms->outgoing_packet_data, 0x0003001c);
    put_le_utf16(&mms->outgoing_packet_data, data_string);
    put_le16(&mms->outgoing_packet_data, 0); // double unicode ended string.

    return send_command_packet(mms);
}

static int asf_header_parser(MMSContext *mms)
{
    uint8_t *p = mms->asf_header, *end = mms->asf_header + mms->asf_header_size;
    mms->stream_num = 0;

    if (mms->asf_header_size < sizeof(ff_asf_guid) * 2 + 22 ||
        memcmp(p, ff_asf_header, sizeof(ff_asf_guid)))
        return -1;

    p += sizeof(ff_asf_guid) + 14;
    do {
        uint64_t chunksize = AV_RL64(p + sizeof(ff_asf_guid));
        if (!memcmp(p, ff_asf_file_header, sizeof(ff_asf_guid))) {
            /* read packet size */
            if (end - p > sizeof(ff_asf_guid) * 2 + 68) {
                mms->asf_packet_len = AV_RL32(p + sizeof(ff_asf_guid) * 2 + 64);
            }
        } else if (!memcmp(p, ff_asf_stream_header, sizeof(ff_asf_guid))) {
            mms->stream_num++;
        }
        if (chunksize > end - p)
            return -1;
        p += chunksize;
    } while (end - p >= sizeof(ff_asf_guid) + 8);

    return 0;
}

/** Send MMST stream selection command based on the AVStream->discard values. */
static int send_stream_selection_request(MMSContext *mms)
{
    int ii;

    //  send the streams we want back...
    start_command_packet(mms, CS_PKT_STREAM_ID_REQUEST);
    put_le32(&mms->outgoing_packet_data, mms->stream_num); // stream nums.
    for(ii= 0; ii<mms->stream_num; ii++) {
        put_le16(&mms->outgoing_packet_data, 0xffff);      // flags
        put_le16(&mms->outgoing_packet_data, ii +1);       // stream id
        put_le16(&mms->outgoing_packet_data, 0);           // selection
    }

    put_le16(&mms->outgoing_packet_data, 0);

    return send_command_packet(mms);
}

static void read_data(MMSContext *mms, uint8_t *buf, const int buf_size, int* result)
{
    int read_size;
    read_size = FFMIN(buf_size, mms->pkt_buf_len);
    memcpy(buf, mms->pkt_read_ptr, read_size);
    mms->pkt_buf_len -= read_size;
    mms->pkt_read_ptr+= read_size;
    *result += read_size;
}

/** Read at most one media packet (or a whole header). */
static int read_mms_packet(MMSContext *mms, uint8_t *buf, int buf_size)
{
    int result = 0;
    MMSSCPacketType packet_type;
    int size_to_copy;

    do {
        if(mms->asf_header_read_pos < mms->asf_header_size) {
            /* Read from ASF header buffer */
            size_to_copy= FFMIN(buf_size,
                                mms->asf_header_size - mms->asf_header_read_pos);
            memcpy(buf, mms->asf_header + mms->asf_header_read_pos, size_to_copy);
            mms->asf_header_read_pos += size_to_copy;
            result += size_to_copy;
            dprintf(NULL, "Copied %d bytes from stored header. left: %d\n",
                   size_to_copy, mms->asf_header_size - mms->asf_header_read_pos);
        } else if(mms->pkt_buf_len) {
            /* Read from media packet buffer */
            read_data(mms, buf, buf_size, &result);
        } else {
            /* Read from network */
            packet_type= get_tcp_server_response(mms);
            switch (packet_type) {
            case SC_PKT_ASF_MEDIA:
               if(mms->pkt_buf_len>mms->asf_packet_len) {
                    dprintf(NULL, "Incoming packet
                            larger than the asf packet size stated (%d>%d)\n",
                            mms->pkt_buf_len, mms->asf_packet_len);
                    result= AVERROR_IO;
                    break;
                }

                // copy the data to the packet buffer.
                read_data(mms, buf, buf_size, &result);
                break;
            case SC_PKT_ASF_HEADER:
                // copy the data to the packet buffer.
                read_data(mms, buf, buf_size, &result);
                break;
            default:
                dprintf(NULL, "Got a unkown Packet Type: 0x%x\n", packet_type);
                break;
            }
        }
    } while(!result); // only return one packet.
    return result;
}

static int send_close_packet(MMSContext *mms)
{
    start_command_packet(mms, CS_PKT_STREAM_CLOSE);
    insert_command_prefixes(mms, 1, 1);

    return send_command_packet(mms);
}

/** Close the MMSH/MMST connection */
static int mms_close(URLContext *h)
{
    MMSContext *mms = (MMSContext *)h->priv_data;

    if(mms->mms_hd) {
        send_close_packet(mms);
        close_connection(mms);
    }

    /* free all separately allocated pointers in mms */
    av_free(mms->asf_header);
    av_freep(&h->priv_data);

    return 0;
}

static int handle_mms_msg_pkt(MMSContext *mms, const MMSSCPacketType packet_type)
{
    int ret = -1;

    switch(packet_type) {
    case SC_PKT_CLIENT_ACCEPTED:
        ret = send_protocol_select(mms);
        break;
    case SC_PKT_PROTOCOL_ACCEPTED:
        ret = send_media_file_request(mms);
        break;
    case SC_PKT_MEDIA_FILE_DETAILS:
        ret = send_media_header_request(mms);
        break;
    case SC_PKT_HEADER_REQUEST_ACCEPTED:
        ret = 0;
        break;
    case SC_PKT_ASF_HEADER:
        if((mms->incoming_flags == 0X08) || (mms->incoming_flags == 0X0C)) {
            ret = asf_header_parser(mms);
            mms->header_parsed = 1;
        }
        break;
    default:
        dprintf(NULL, "Unhandled packet type %d\n", packet_type);
        break;
    }
    return ret;
}

static int mms_open_cnx(URLContext *h)
{
    MMSContext *mms = h->priv_data;
    MMSSCPacketType packet_type;

    char tcpname[256];
    int err = AVERROR(EIO);
    int ret;

    // only for MMS over TCP, so set proto = NULL
    url_split(NULL, 0, NULL, 0,
            mms->host, sizeof(mms->host), &mms->port, mms->path,
            sizeof(mms->path), mms->location);

    if(mms->port<0)
        mms->port = 1755; // defaut mms protocol port

    /* the outgoing packet buffer */
    init_put_byte(&mms->outgoing_packet_data, mms->outgoing_packet_buffer,
                  sizeof(mms->outgoing_packet_buffer), 1, NULL,
                  NULL, NULL, NULL);
    // establish tcp connection.
    close_connection(mms);
    snprintf(tcpname, sizeof(tcpname), "tcp://%s:%d", mms->host, mms->port);
    err = url_open(&mms->mms_hd, tcpname, URL_RDWR);
    if (err)
        goto fail;

    mms->local_ip_address = 0xc0a80081; // local IP address,server don't care.
    mms->local_port       = 1037;       // as above,could be arbitrary value.
    mms->packet_id        = 3;          // default, initial value.
    mms->header_packet_id = 2;          // default, initial value.

    send_startup_packet(mms);
    while (!mms->header_parsed) {
        packet_type = get_tcp_server_response(mms);
        ret = handle_mms_msg_pkt(mms, packet_type);
        if (ret < 0)
            break;
    }

    if (!mms->asf_packet_len || !mms->stream_num)
        goto fail;

    dprintf(NULL, "Leaving open (success)\n");
    return 0;
fail:
    mms_close(h);
    dprintf(NULL, "Leaving open (failure: %d)\n", err);
    return err;
}

static int mms_open(URLContext *h, const char *uri, int flags)
{
    MMSContext *mms;

    h->is_streamed = 1;
    mms = av_malloc(sizeof(MMSContext));
    if (!mms)
        return AVERROR(ENOMEM);
    memset(mms, 0, sizeof(MMSContext));
    h->priv_data = mms;
    av_strlcpy(mms->location, uri, sizeof(mms->location));

    return mms_open_cnx(h);
}

static int send_media_packet_request(MMSContext *mms)
{
    start_command_packet(mms, CS_PKT_START_FROM_PKT_ID);
    insert_command_prefixes(mms, 1, 0x0001FFFF);
    put_le64(&mms->outgoing_packet_data, 0);          // seek timestamp
    put_le32(&mms->outgoing_packet_data, 0xffffffff); // unknown
    put_le32(&mms->outgoing_packet_data, 0xffffffff); // packet offset
    put_byte(&mms->outgoing_packet_data, 0xff);       // max stream time limit
    put_byte(&mms->outgoing_packet_data, 0xff);       // max stream time limit
    put_byte(&mms->outgoing_packet_data, 0xff);       // max stream time limit
    put_byte(&mms->outgoing_packet_data, 0x00);       // stream time limit flag

    mms->packet_id++;                                 // new packet_id
    put_le32(&mms->outgoing_packet_data, mms->packet_id);
    return send_command_packet(mms);
}


static void clear_stream_buffers(MMSContext *mms)
{
    mms->pkt_buf_len = 0;
    mms->pkt_read_ptr = mms->pkt_buf;
}

/** Read ASF data through the protocol. */
static int mms_read(URLContext *h, uint8_t *buf, int size)
{
    /* TODO: see tcp.c:tcp_read() about a possible timeout scheme */
    MMSContext *mms = h->priv_data;
    int result = 0;

    /* Since we read the header at open(), this shouldn't be possible */
    assert(mms->header_parsed);

    if(mms->header_parsed) {
        if (mms->asf_header_read_pos >= mms->asf_header_size
            && !mms->streaming_flag) {
            dprintf(NULL, "mms_read() before play().\n");
            clear_stream_buffers(mms);
            result = send_stream_selection_request(mms);
            if(result < 0)
                return result;
            if (get_tcp_server_response(mms) != SC_PKT_STREAM_ID_ACCEPTED) {
                dprintf(NULL, "Can't get stream id accepted packet.\n");
                return 0;
            }

            // send media packet request
            send_media_packet_request(mms);
            if (get_tcp_server_response(mms) == SC_PKT_MEDIA_PKT_FOLLOWS) {
               mms->streaming_flag = 1;
            } else {
                dprintf(NULL, "Canot get media follows packet from server.\n");
                return 0;
            }
        }
        result = read_mms_packet(mms, buf, size);
    }

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
