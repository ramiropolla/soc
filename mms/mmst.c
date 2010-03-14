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
    CS_PACKET_INITIAL_TYPE= 0x01,
    CS_PACKET_PROTOCOL_SELECT_TYPE= 0x02,
    CS_PACKET_MEDIA_FILE_REQUEST_TYPE= 0x05,
    CS_PACKET_START_FROM_PACKET_ID_TYPE= 0x07,
    CS_PACKET_STREAM_PAUSE_TYPE= 0x09, // tcp left open, but data stopped.
    CS_PACKET_STREAM_CLOSE_TYPE= 0x0d,
    CS_PACKET_MEDIA_HEADER_REQUEST_TYPE= 0x15,
    CS_PACKET_TIMING_DATA_REQUEST_TYPE= 0x18,
    CS_PACKET_USER_PASSWORD_TYPE= 0x1a,
    CS_PACKET_KEEPALIVE_TYPE= 0x1b,
    CS_PACKET_STREAM_ID_REQUEST_TYPE= 0x33,
} MMSCSPacketType;

/** Server to client packet types. */
typedef enum {
    /** Control packets. */
    /*@{*/
    SC_PACKET_CLIENT_ACCEPTED= 0x01,
    SC_PACKET_PROTOCOL_ACCEPTED_TYPE= 0x02,
    SC_PACKET_PROTOCOL_FAILED_TYPE= 0x03,
    SC_PACKET_MEDIA_PACKET_FOLLOWS_TYPE= 0x05,
    SC_PACKET_MEDIA_FILE_DETAILS_TYPE= 0x06,
    SC_PACKET_HEADER_REQUEST_ACCEPTED_TYPE= 0x11,
    SC_PACKET_TIMING_TEST_REPLY_TYPE= 0x15,
    SC_PACKET_PASSWORD_REQUIRED_TYPE= 0x1a,
    SC_PACKET_KEEPALIVE_TYPE= 0x1b,
    SC_PACKET_STREAM_STOPPED_TYPE= 0x1e, // mmst, mmsh
    SC_PACKET_STREAM_CHANGING_TYPE= 0x20,
    SC_PACKET_STREAM_ID_ACCEPTED_TYPE= 0x21,
    /*@}*/

    /** Pseudo packets. */
    /*@{*/
    SC_PACKET_TYPE_CANCEL = -1, // mmst
    SC_PACKET_TYPE_NO_DATA = -2, // mmst
    SC_PACKET_HTTP_CONTROL_ACKNOWLEDGE = -3, // mmsh
    /*@}*/

    /** Data packets. */
    /*@{*/
    SC_PACKET_ASF_HEADER_TYPE= 0x81, // mmst, mmsh
    SC_PACKET_ASF_MEDIA_TYPE= 0x82, // mmst, mmsh
    /*@}*/
} MMSSCPacketType;

typedef struct {
    uint32_t local_ip_address; ///< Not ipv6 compatible, but neither is the protocol (sent, but not correct).
    int local_port; ///< My local port (sent but not correct).
    int sequence_number; ///< Outgoing packet sequence number.
    char path[256]; ///< Path of the resource being asked for.
    char host[128]; ///< Host of the resources.
    int port; ///< Port of the resource.

    URLContext *mms_hd; ///< TCP connection handle

    /** Buffer for outgoing packets. */
    /*@{*/
    ByteIOContext outgoing_packet_data; ///< Outgoing packet stream
    uint8_t outgoing_packet_buffer[512]; ///< Outgoing packet data
    /*@}*/

    /** Buffer for incoming control packets. */
    /*@{*/
    uint8_t incoming_buffer[8192]; ///< Incoming buffer location.
    int incoming_buffer_length; ///< Incoming buffer length.
    /*@}*/

    /** Buffer for incoming media/header packets. */
    /*@{*/
    uint8_t media_packet_incoming_buffer[8192]; ///< Either a header or media packet.
    uint8_t *media_packet_read_ptr; ///< Pointer for partial reads.
    int media_packet_buffer_length; ///< Buffer length.
    int media_packet_seek_offset;   ///< Additional offset into packet from seek.
    /*@}*/

    int incoming_packet_seq; ///< Incoming packet sequence number.
    int incoming_flags; ///< Incoming packet flags.

    int packet_id; ///< Identifier for packets in the current stream, incremented on stops, etc.
    unsigned int header_packet_id; ///< The header packet id (default is 2, can be reset)

    int seekable; ///< This tells you if the stream is seekable.

    /** Internal handling of the ASF header */
    /*@{*/
    uint8_t *asf_header; ///< Stored ASF header, for seeking into it and internal parsing.
    int asf_header_size; ///< Size of stored ASF header.
    int asf_header_read_pos; ///< Current read position in header. See read_packet().
    int header_parsed; ///< The header has been received and parsed.
    int asf_packet_len;
    /*@}*/

    int pause_resume_seq; ///< Last packet returned by mms_read. Useful for resuming pause.
    char location[4096];
    int stream_num;
} MMSContext;

/** Close the remote connection. */
static void close_connection(MMSContext *mms)
{
    url_close(mms->mms_hd);
}

/** Open or reopen (http) the remote connection. */
static int ff_mms_open_connection(MMSContext *mms)
{
    char tcpname[256];
    int err;

    close_connection(mms);

    snprintf(tcpname, sizeof(tcpname), "tcp://%s:%d", mms->host, mms->port);
    err = url_open(&mms->mms_hd, tcpname, URL_RDWR);
    return err;
}

/** Create MMST command packet header */
static void start_command_packet(MMSContext *mms, MMSCSPacketType packet_type)
{
    ByteIOContext *context= &mms->outgoing_packet_data;

    url_fseek(context, 0, SEEK_SET); // start at the beginning...
    put_le32(context, 1); // start sequence?
    put_le32(context, 0xb00bface);
    put_le32(context, 0); // Length of command until the end of all data  Value is in bytes and starts from after the protocol type bytes
    put_byte(context, 'M'); put_byte(context, 'M'); put_byte(context, 'S'); put_byte(context, ' ');
    put_le32(context, 0);
    put_le32(context, mms->sequence_number++);
    put_le64(context, 0); // timestmamp
    put_le32(context, 0);
    put_le16(context, packet_type);
    put_le16(context, 3); // direction- to server
}

/** Add prefixes to MMST command packet. */
static void insert_command_prefixes(MMSContext *mms,
        uint32_t prefix1, uint32_t prefix2)
{
    ByteIOContext *context= &mms->outgoing_packet_data;

    put_le32(context, prefix1); // first prefix
    put_le32(context, prefix2); // second prefix
}

/** Write a utf-16 string in little-endian order.
 * @note This is NOT the same as static int ascii_to_wc (ByteIOContext *pb, uint8_t *b), because ascii_to_wc is big endian, and this is little endian.
 */
static void put_le_utf16(ByteIOContext *pb, char *utf8)
{
    int val;

    while(*utf8) {
        GET_UTF8(val, *utf8++, break;); // goto's suck, but i want to make sure it's terminated.
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

    // first adjust the header fields (the lengths)...
    url_fseek(context, 8, SEEK_SET);
    put_le32(context, first_length);
    url_fseek(context, 16, SEEK_SET);
    put_le32(context, len8);
    url_fseek(context, 32, SEEK_SET);
    put_le32(context, len8-2);

    // seek back to the end (may not be necessary...)
    url_fseek(context, exact_length, SEEK_SET);

    // write it out...
    write_result= url_write(mms->mms_hd, context->buffer, exact_length);
    if(write_result != exact_length) {
        dprintf(NULL, "url_write returned: %d != %d\n", write_result, exact_length);
        return AVERROR_IO;
    }

    return 0;
}

static int send_protocol_select(MMSContext *mms)
{
    char data_string[256];

    // send the timing request packet...
    start_command_packet(mms, CS_PACKET_PROTOCOL_SELECT_TYPE);
    insert_command_prefixes(mms, 0, 0);
    put_le32(&mms->outgoing_packet_data, 0);  // timestamp?
    put_le32(&mms->outgoing_packet_data, 0);  // timestamp?
    put_le32(&mms->outgoing_packet_data, 2);
    snprintf(data_string, sizeof(data_string), "\\\\%d.%d.%d.%d\\%s\\%d",
            (mms->local_ip_address>>24)&0xff,
            (mms->local_ip_address>>16)&0xff,
            (mms->local_ip_address>>8)&0xff,
            mms->local_ip_address&0xff,
            "TCP", // or UDP
            mms->local_port);
    put_le_utf16(&mms->outgoing_packet_data, data_string);
    put_le16(&mms->outgoing_packet_data, 0x30);

    return send_command_packet(mms);
}

static int send_media_file_request(MMSContext *mms)
{
    start_command_packet(mms, CS_PACKET_MEDIA_FILE_REQUEST_TYPE);
    insert_command_prefixes(mms, 1, 0xffffffff);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le_utf16(&mms->outgoing_packet_data, mms->path+1); // +1 because we skip the leading /
    put_le32(&mms->outgoing_packet_data, 0); /* More zeroes */

    return send_command_packet(mms);
}

static int read_bytes(MMSContext *mms, uint8_t *buffer, int length_to_read)
{
    int len= 0;

    while(len<length_to_read)
    {
        int read_result= url_read(mms->mms_hd, buffer+len, length_to_read-len);
        if(read_result)
        {
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

    // read these from the incoming buffer.. (40 is the packet header size, without the prefixes)
    init_put_byte(&pkt, mms->incoming_buffer+40, mms->incoming_buffer_length-40, 0, NULL, NULL, NULL, NULL);
    get_le32(&pkt); // prefix 1
    mms->header_packet_id= (get_le32(&pkt) & 0xff); // prefix 2
    dprintf(NULL, "Changed header prefix to 0x%x", mms->header_packet_id);
}

static int send_keepalive_packet(MMSContext *mms)
{
    // respond to a keepalive with a keepalive...
    start_command_packet(mms, CS_PACKET_KEEPALIVE_TYPE);
    insert_command_prefixes(mms, 1, 0x100FFFF);
    return send_command_packet(mms);
}

/** Pad media packets smaller than max_packet_size and/or adjust read position
  * after a seek. */
static void pad_media_packet(MMSContext *mms)
{
    if(mms->media_packet_buffer_length<mms->asf_packet_len) {
        int padding_size = mms->asf_packet_len - mms->media_packet_buffer_length;
        memset(mms->media_packet_incoming_buffer+mms->media_packet_buffer_length, 0, padding_size);
        mms->media_packet_buffer_length += padding_size;
    }
    if(mms->media_packet_seek_offset) {
        mms->media_packet_buffer_length -= mms->media_packet_seek_offset;
        mms->media_packet_read_ptr += mms->media_packet_seek_offset;
        mms->media_packet_seek_offset = 0;
    }
}

/** Read incoming MMST media, header or command packet. */
static MMSSCPacketType get_tcp_server_response(MMSContext *mms)
{
    int read_result;
    MMSSCPacketType packet_type= -1;
    int done;

    do {
        done= 1; // assume we're going to get a valid packet.
        if((read_result= read_bytes(mms, mms->incoming_buffer, 8))==8) {
            // check if we are a command packet...
            if(AV_RL32(mms->incoming_buffer + 4)==0xb00bface) {
                mms->incoming_flags= mms->incoming_buffer[3];
                if((read_result= read_bytes(mms, mms->incoming_buffer+8, 4)) == 4) {
                    int length_remaining= AV_RL32(mms->incoming_buffer+8) + 4;

                    dprintf(NULL, "Length remaining is %d\n", length_remaining);
                    // FIXME? ** VERIFY LENGTH REMAINING HAS SPACE
                    // read the rest of the packet....
                    read_result = read_bytes(mms, mms->incoming_buffer + 12, length_remaining) ;
                    if (read_result == length_remaining) {
                        // we have it all; get the stuff out of it.
                        mms->incoming_buffer_length= length_remaining+12;

                        // get the packet type...
                        packet_type= AV_RL16(mms->incoming_buffer+36);

                    } else {
                        // read error...
                        dprintf(NULL, "3 read returned %d!\n", read_result);
                    }
                } else {
                    // read error...
                    dprintf(NULL, "2 read returned %d!\n", read_result);
                }
            } else {
                int length_remaining= (AV_RL16(mms->incoming_buffer + 6) - 8) & 0xffff;
                uint8_t *dst= mms->media_packet_incoming_buffer;
                int packet_id_type;

                assert(mms->media_packet_buffer_length==0); // assert all has been consumed.

                //** VERIFY LENGTH REMAINING HAS SPACE
                // note we cache the first 8 bytes, then fill up the buffer with the others
                mms->incoming_packet_seq        = AV_RL32(mms->incoming_buffer);
                packet_id_type                  = mms->incoming_buffer[4]; // NOTE: THIS IS THE ONE I CAN CHANGE
                mms->incoming_flags             = mms->incoming_buffer[5];
                mms->media_packet_buffer_length = length_remaining;
                mms->media_packet_read_ptr      = mms->media_packet_incoming_buffer;

                if(mms->media_packet_buffer_length>=sizeof(mms->media_packet_incoming_buffer)) {
                    dprintf(NULL, "Incoming Buffer Length exceeds buffer: %d>%d\n", mms->media_packet_buffer_length, (int) sizeof(mms->media_packet_incoming_buffer));
                }
                assert(mms->media_packet_buffer_length<sizeof(mms->media_packet_incoming_buffer));
                read_result= read_bytes(mms, dst, length_remaining);
                if(read_result != length_remaining) {
                    dprintf(NULL, "read_bytes result: %d asking for %d\n", read_result, length_remaining);
                    break;
                } else {
                    // if we successfully read everything....
                    if(packet_id_type == mms->header_packet_id) {
                        // asf header
                        packet_type = SC_PACKET_ASF_HEADER_TYPE;
                        // Store the asf header
                        if(!mms->header_parsed) {
                            mms->asf_header = av_realloc(mms->asf_header, mms->asf_header_size + mms->media_packet_buffer_length);
                            memcpy(mms->asf_header + mms->asf_header_size, mms->media_packet_read_ptr, mms->media_packet_buffer_length);
                            mms->asf_header_size += mms->media_packet_buffer_length;
                        }
                    } else if(packet_id_type == mms->packet_id) {
                        packet_type = SC_PACKET_ASF_MEDIA_TYPE;
                    } else {
                        dprintf(NULL, "packet id type %d which must be old, getting another one.", packet_id_type);
                        done= 0;
                    }
                }
            }
        } else {
            // read error...
            if(read_result<0) {
                dprintf(NULL, "Read error (or cancelled) returned %d!\n", read_result);
                packet_type = SC_PACKET_TYPE_CANCEL;
            } else {// 0 is okay, no data received.
                dprintf(NULL, "Read result of zero?!\n");
                packet_type = SC_PACKET_TYPE_NO_DATA;
            }
            done = 1;
        }
    } while(!done);

    if (packet_type == SC_PACKET_KEEPALIVE_TYPE) {
        send_keepalive_packet(mms);
    }
    if (packet_type == SC_PACKET_STREAM_CHANGING_TYPE) {
        handle_packet_stream_changing_type(mms);
        //TODO: Handle new header when change the stream type.
    }
    if (packet_type == SC_PACKET_ASF_MEDIA_TYPE) {
        pad_media_packet(mms);
    }
    return packet_type;
}

static void handle_packet_media_file_details(MMSContext *mms)
{
    ByteIOContext pkt;
    uint16_t broadcast_flags;
    int64_t total_file_length_in_seconds;
    uint32_t total_length_in_seconds;
    uint32_t packet_length;
    uint32_t total_packet_count;
    uint32_t highest_bit_rate;
    uint32_t header_size;
    uint32_t flags;
    double duration;

    // read these from the incoming buffer.. (48 is the packet header size)
    init_put_byte(&pkt, mms->incoming_buffer+48, mms->incoming_buffer_length-48, 0, NULL, NULL, NULL, NULL);
    flags= get_le32(&pkt); // flags?
    if(flags==0xffffffff) {
        // this is a permission denied event.
        dprintf(NULL, "Permission denied!\n");
    } else {
        get_le32(&pkt);
        get_le32(&pkt);
        get_le16(&pkt);
        broadcast_flags= get_le16(&pkt);

        total_file_length_in_seconds= get_le64(&pkt);
        duration= av_int2dbl(total_file_length_in_seconds);
        total_length_in_seconds= get_le32(&pkt);
        get_le32(&pkt);
        get_le32(&pkt);
        get_le32(&pkt);
        get_le32(&pkt);
        packet_length= get_le32(&pkt);
        total_packet_count= get_le32(&pkt);
        get_le32(&pkt);
        highest_bit_rate= get_le32(&pkt);
        header_size= get_le32(&pkt);
        dprintf(NULL, "Broadcast flags: 0x%x\n", broadcast_flags); // 8000= allow index, 01= prerecorded, 02= live 42= presentation with script commands
        dprintf(NULL, "File Time Point?: %lld double size: %d double value: %lf\n", total_file_length_in_seconds, (int) sizeof(double), duration);
        dprintf(NULL, "Total in Seconds: %d\n", total_length_in_seconds);
        dprintf(NULL, "Packet length: %d\n", packet_length);
        dprintf(NULL, "Total Packet Count: %d\n", total_packet_count);
        dprintf(NULL, "Highest Bit Rate: %d\n", highest_bit_rate);
        dprintf(NULL, "Header Size: %d\n", header_size);
        dprintf(NULL, "---- Done ----\n");

        /* Disable seeking in live broadcasts for now */
        if(! (broadcast_flags & 0x0200))
            mms->seekable = 1;
    }
}

static int send_media_header_request(MMSContext *mms)
{
    start_command_packet(mms, CS_PACKET_MEDIA_HEADER_REQUEST_TYPE);
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

    snprintf(data_string, sizeof(data_string), "NSPlayer/7.0.0.1956; {%s}; Host: %s",
            "7E667F5D-A661-495E-A512-F55686DDA178", mms->host);

    start_command_packet(mms, CS_PACKET_INITIAL_TYPE);
    insert_command_prefixes(mms, 0, 0x0004000b);
    put_le32(&mms->outgoing_packet_data, 0x0003001c);
    put_le_utf16(&mms->outgoing_packet_data, data_string);
    put_le16(&mms->outgoing_packet_data, 0); // double unicode ended string...

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
    start_command_packet(mms, CS_PACKET_STREAM_ID_REQUEST_TYPE);
    put_le32(&mms->outgoing_packet_data, mms->stream_num); // stream nums.
    for(ii= 0; ii<mms->stream_num; ii++) {
        put_le16(&mms->outgoing_packet_data, 0xffff); // flags
        put_le16(&mms->outgoing_packet_data, ii +1); // stream id
       put_le16(&mms->outgoing_packet_data, 0); // selection
    }

    put_le16(&mms->outgoing_packet_data, 0); /* Extra zeroes */

    return send_command_packet(mms);
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
            size_to_copy= FFMIN(buf_size, mms->asf_header_size - mms->asf_header_read_pos);
            memcpy(buf, mms->asf_header + mms->asf_header_read_pos, size_to_copy);
            mms->asf_header_read_pos += size_to_copy;
            result += size_to_copy;
            dprintf(NULL, "Copied %d bytes from stored header. left: %d\n", size_to_copy, mms->asf_header_size - mms->asf_header_read_pos);
        } else if(mms->media_packet_buffer_length) {
            /* Read from media packet buffer */
            size_to_copy = FFMIN(buf_size, mms->media_packet_buffer_length);
            memcpy(buf, mms->media_packet_read_ptr, size_to_copy);
            mms->media_packet_buffer_length -= size_to_copy;
            mms->media_packet_read_ptr+= size_to_copy;
            result += size_to_copy;
        } else {
            /* Read from network */
            packet_type= get_tcp_server_response(mms);
            switch (packet_type) {
            case SC_PACKET_ASF_MEDIA_TYPE:
               if(mms->media_packet_buffer_length>mms->asf_packet_len) {
                    dprintf(NULL, "Incoming packet larger than the asf packet size stated (%d>%d)\n", mms->media_packet_buffer_length, mms->asf_packet_len);
                    result= AVERROR_IO;
                    break;
                }

                // copy the data to the packet buffer...
                size_to_copy= FFMIN(buf_size, mms->media_packet_buffer_length);
                memcpy(buf, mms->media_packet_read_ptr, size_to_copy);
                mms->media_packet_buffer_length -= size_to_copy;
                mms->media_packet_read_ptr += size_to_copy;
                result += size_to_copy;
                break;
            case SC_PACKET_ASF_HEADER_TYPE:
                // copy the data to the packet buffer...
                size_to_copy= FFMIN(buf_size, mms->media_packet_buffer_length);
                memcpy(buf, mms->media_packet_read_ptr, size_to_copy);
                mms->media_packet_buffer_length -= size_to_copy;
                mms->media_packet_read_ptr+= size_to_copy;
                result+= size_to_copy;
                break;
            default:
                dprintf(NULL, "Got a unkown Packet Type: 0x%x\n", packet_type);
                break;
            }
        }
    } while(!result); // only return one packet...
    return result;
}

static int send_close_packet(MMSContext *mms)
{
    start_command_packet(mms, CS_PACKET_STREAM_CLOSE_TYPE);
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

    /* TODO: free all separately allocated pointers in mms */
    av_free(mms->asf_header);
    av_freep(&h->priv_data);

    return 0;
}

static int mms_open_cnx(URLContext *h)
{
    MMSContext *mms = h->priv_data;

    char authorization[64];
    int err = AVERROR(EIO);

    // only for MMS over TCP, so set proto = NULL
    url_split(NULL, 0, authorization, sizeof(authorization), mms->host, sizeof(mms->host),
              &mms->port, mms->path, sizeof(mms->path), mms->location);

    if(mms->port<0)
        mms->port = 1755; // defaut mms protocol port

    /* the outgoing packet buffer */
    init_put_byte(&mms->outgoing_packet_data, mms->outgoing_packet_buffer, sizeof(mms->outgoing_packet_buffer), 1, NULL, NULL, NULL, NULL);

    /* open the tcp connexion */
    if((err = ff_mms_open_connection(mms)))
        goto fail;

    // Fill in some parameters...
    mms->local_ip_address = 0xc0a80081; // This should be the local IP address; how do I get this from the url_ stuff?  (nothing is apparent)
    mms->local_port = 1037; // as above, this should be the port I am connected to; how do I get this frmo the url stuff? (Server doesn't really seem to care too much)
    mms->packet_id = 3; // default, initial value. (3 will be incremented to 4 before first use)
    mms->header_packet_id = 2; // default, initial value.

    send_startup_packet(mms);
    if (get_tcp_server_response(mms) == SC_PACKET_CLIENT_ACCEPTED) {
//        start_command_packet(mms, CS_PACKET_TIMING_DATA_REQUEST_TYPE);
//        insert_command_prefixes(mms, 0xf0f0f0f1, 0x0004000b);
//        send_command_packet(mms);
            send_protocol_select(mms);
    } else
        goto fail;

//    if (get_tcp_server_response(mms) == SC_PACKET_TIMING_TEST_REPLY_TYPE) {
//        send_protocol_select(mms);
//    } else
//        goto fail;

    if (get_tcp_server_response(mms) == SC_PACKET_PROTOCOL_ACCEPTED_TYPE) {
        send_media_file_request(mms);
    } else
        goto fail;

    if (get_tcp_server_response(mms) == SC_PACKET_MEDIA_FILE_DETAILS_TYPE) {
        handle_packet_media_file_details(mms);
        send_media_header_request(mms);
    } else
        goto fail;

    if (get_tcp_server_response(mms) == SC_PACKET_HEADER_REQUEST_ACCEPTED_TYPE) {
        // recv asf header data
        if (get_tcp_server_response(mms) == SC_PACKET_ASF_HEADER_TYPE) {
            if((mms->incoming_flags == 0X08) || (mms->incoming_flags == 0X0C)) {
                asf_header_parser(mms);
                mms->header_parsed = 1;
            }  else
                goto fail;
        } else
            goto fail;
    } else
        goto fail;

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
    start_command_packet(mms, CS_PACKET_START_FROM_PACKET_ID_TYPE);
    insert_command_prefixes(mms, 1, 0x0001FFFF);
    put_le64(&mms->outgoing_packet_data, 0); // seek timestamp
    put_le32(&mms->outgoing_packet_data, 0xffffffff);  // unknown
    put_le32(&mms->outgoing_packet_data, 0xffffffff);    // packet offset
    put_byte(&mms->outgoing_packet_data, 0xff); // max stream time limit
    put_byte(&mms->outgoing_packet_data, 0xff); // max stream time limit
    put_byte(&mms->outgoing_packet_data, 0xff); // max stream time limit
    put_byte(&mms->outgoing_packet_data, 0x00); // stream time limit flag

    mms->packet_id++; // new packet_id so we can separate new data from old data
    put_le32(&mms->outgoing_packet_data, mms->packet_id);
    return send_command_packet(mms);
}

/** Read ASF data through the protocol. */
static int mms_read(URLContext *h, uint8_t *buf, int size)
{
    /* TODO: see tcp.c:tcp_read() about a possible timeout scheme */
    MMSContext *mms = h->priv_data;
    int result = 0;

    /* Since we read the header at open(), this shouldn't be possible */
    assert(mms->header_parsed);

    /* Automatically start playing if the app wants to read before it has called play()
     * (helps with non-streaming aware apps) */
    if(mms->header_parsed) {
        if (mms->asf_header_read_pos >= mms->asf_header_size) {
            dprintf(NULL, "mms_read() before play(). Playing automatically.\n");
            result = send_stream_selection_request(mms);
            if(result < 0)
                return result;
            if (get_tcp_server_response(mms) != SC_PACKET_STREAM_ID_ACCEPTED_TYPE) {
                dprintf(NULL, "Canot get stream id accepted packet from server.\n");
                return 0;
            }

            // send media packet request
            send_media_packet_request(mms);
            if (get_tcp_server_response(mms) != SC_PACKET_MEDIA_PACKET_FOLLOWS_TYPE) {
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
