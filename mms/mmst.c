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
#include "network.h"
#include "asf.h"

#define MMS_DEBUG_LEVEL 2
#define MMS_MAXIMUM_PACKET_LENGTH 512
#define MMS_KILO                  1024
#define MMS_URL_SIZE          4096
#define DEFAULT_MMS_PORT      1755

/** State machine states. */
typedef enum {
    AWAITING_SC_PACKET_CLIENT_ACCEPTED= 0,
    AWAITING_SC_PACKET_TIMING_TEST_REPLY_TYPE,
    AWAITING_CS_PACKET_PROTOCOL_ACCEPTANCE,
    AWAITING_PASSWORD_QUERY_OR_MEDIA_FILE,
    AWAITING_PACKET_HEADER_REQUEST_ACCEPTED_TYPE,
    AWAITING_STREAM_ID_ACCEPTANCE,
    AWAITING_STREAM_START_PACKET,
    AWAITING_ASF_HEADER,
    ASF_HEADER_DONE,
    AWAITING_PAUSE_ACKNOWLEDGE,
    AWAITING_HTTP_PAUSE_CONTROL_ACKNOWLEDGE,
    STREAMING,
    STREAM_DONE,
    STATE_ERROR,
    STREAM_PAUSED,
    USER_CANCELLED
} MMSState;

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

#if (MMS_DEBUG_LEVEL>0)
static const char *state_names[]= {
    "AWAITING_SC_PACKET_CLIENT_ACCEPTED",
    "AWAITING_SC_PACKET_TIMING_TEST_REPLY_TYPE",
    "AWAITING_CS_PACKET_PROTOCOL_ACCEPTANCE",
    "AWAITING_PASSWORD_QUERY_OR_MEDIA_FILE",
    "AWAITING_PACKET_HEADER_REQUEST_ACCEPTED_TYPE",
    "AWAITING_STREAM_ID_ACCEPTANCE",
    "AWAITING_STREAM_START_PACKET",
    "AWAITING_ASF_HEADER",
    "ASF_HEADER_DONE",
    "AWAITING_PAUSE_ACKNOWLEDGE",
    "AWAITING_HTTP_PAUSE_CONTROL_ACKNOWLEDGE",
    "STREAMING",
    "STREAM_DONE",
    "STATE_ERROR",
    "STREAM_PAUSED",
    "USER_CANCELLED"
};
#endif

typedef struct {
    char local_guid[37]; ///< My randomly generated GUID.
    uint32_t local_ip_address; ///< Not ipv6 compatible, but neither is the protocol (sent, but not correct).
    int local_port; ///< My local port (sent but not correct).
    int sequence_number; ///< Outgoing packet sequence number.
    MMSState state; ///< Packet state machine current state.
    char path[256]; ///< Path of the resource being asked for.
    char host[128]; ///< Host of the resources.
    int port; ///< Port of the resource.

    URLContext *mms_hd; ///< TCP connection handle
    ByteIOContext incoming_io_buffer; ///< Incoming data on the socket

    /** Buffer for outgoing packets. */
    /*@{*/
    ByteIOContext outgoing_packet_data; ///< Outgoing packet stream
    uint8_t outgoing_packet_buffer[MMS_MAXIMUM_PACKET_LENGTH]; ///< Outgoing packet data
    /*@}*/

    /** Buffer for incoming control packets. */
    /*@{*/
    uint8_t incoming_buffer[8*MMS_KILO]; ///< Incoming buffer location.
    int incoming_buffer_length; ///< Incoming buffer length.
    /*@}*/

    /** Buffer for incoming media/header packets. */
    /*@{*/
    uint8_t media_packet_incoming_buffer[8*MMS_KILO]; ///< Either a header or media packet.
    uint8_t *media_packet_read_ptr; ///< Pointer for partial reads.
    int media_packet_buffer_length; ///< Buffer length.
    int media_packet_seek_offset;   ///< Additional offset into packet from seek.
    /*@}*/

    int incoming_packet_seq; ///< Incoming packet sequence number.
    int incoming_flags; ///< Incoming packet flags.

    int packet_id; ///< Identifier for packets in the current stream, incremented on stops, etc.
    unsigned int header_packet_id; ///< The header packet id (default is 2, can be reset)

    int seekable; ///< This tells you if the stream is seekable.

    int http_client_id; ///< HTTP's client id.
    int http_play_rate; ///< Rate of playback (1 for normal, 5 or -5 for ffwd/rewind)

    /** Internal handling of the ASF header */
    /*@{*/
    uint8_t *asf_header; ///< Stored ASF header, for seeking into it and internal parsing.
    int asf_header_size; ///< Size of stored ASF header.
    int asf_header_read_pos; ///< Current read position in header. See read_packet().
    int header_parsed; ///< The header has been received and parsed.
    AVFormatContext private_av_format_ctx; ///< Private parsed header data (generic).
    ASFContext      asf_context;           ///< Private parsed header data (ASF-specific).
    AVFormatContext *av_format_ctx; ///< Optional external format context (for stream selection).
    /*@}*/

    int pause_resume_seq; ///< Last packet returned by mms_read. Useful for resuming pause.
    // new added on 2010.2.21
    char location[MMS_URL_SIZE];
} MMSContext;

/** Perform state transition. */
static void ff_mms_set_state(MMSContext *mms, int new_state)
{
    /* Can't exit error state */
    if(mms->state==STATE_ERROR) {
#if (MMS_DEBUG_LEVEL>0)
        fprintf(stderr, "Trying to set state to %s from %s!\n", state_names[new_state], state_names[mms->state]);
#endif
        return;
    }

#if (MMS_DEBUG_LEVEL>0)
    fprintf(stderr, "Set state to %s (%d) from %s (%d)!\n", state_names[new_state], new_state, state_names[mms->state], mms->state);
#endif
    if(mms->state==new_state && new_state==USER_CANCELLED) {
        ff_mms_set_state(mms, STATE_ERROR);
        return;
    }

    mms->state= new_state;
}

/** Close the remote connection. */
static void close_connection(MMSContext *mms)
{
    av_freep(&mms->incoming_io_buffer.buffer);
    url_close(mms->mms_hd);
}

/** Open or reopen (http) the remote connection. */
static int ff_mms_open_connection(MMSContext *mms)
{
    char tcpname[256];
    int flags, err;

    close_connection(mms);

    snprintf(tcpname, sizeof(tcpname), "tcp://%s:%d", mms->host, mms->port);
    err = url_open(&mms->mms_hd, tcpname, URL_RDWR);
    if(err == 0) {
        /* open the incoming and outgoing connections; you can't open a single one with read/write, because it only has one buffer, not two. */
        /* you can't use url_fdopen if the flags of the mms_hd have a WR component, because it will screw up (returning data that is uninitialized) */
        flags = mms->mms_hd->flags;
        mms->mms_hd->flags = URL_RDONLY;
        err = url_fdopen(&mms->incoming_io_buffer, mms->mms_hd);
        mms->mms_hd->flags = flags;
        if(err) {
            // should have a url_fdclose()
            url_close(mms->mms_hd);
        }
    }

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
//    int len8= (first_length+7)/8;
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

    // print it out...
    //    print_command(context->buffer, exact_length);

    // write it out...
    write_result= url_write(mms->mms_hd, context->buffer, exact_length);
    if(write_result != exact_length) {
#if (MMS_DEBUG_LEVEL>0)
        fprintf(stderr, "url_write returned: %d != %d\n", write_result, exact_length);
#endif

        ff_mms_set_state(mms, STATE_ERROR);
        return AVERROR_IO;
    }

    return 0;
}


/** Log unexpected incoming packet */
void log_packet_in_wrong_state(MMSContext *mms, MMSSCPacketType packet_type)
{
#if (MMS_DEBUG_LEVEL>0)
    if(packet_type>=0) {
        fprintf(stderr, "Got a packet 0x%02x in the wrong state: %s (%d)!\n", packet_type, state_names[mms->state], mms->state);
    } else {
        fprintf(stderr, "Got a pseudo-packet %d in the wrong state: %s (%d)!\n", packet_type, state_names[mms->state], mms->state);
    }
#endif
}

static int send_protocol_select(MMSContext *mms)
{
    char data_string[256];
    int err= 0;

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

    err = send_command_packet(mms);
    ff_mms_set_state(mms, AWAITING_CS_PACKET_PROTOCOL_ACCEPTANCE);
    return err;
}

static int send_media_file_request(MMSContext *mms)
{
    int err= 0;
    start_command_packet(mms, CS_PACKET_MEDIA_FILE_REQUEST_TYPE);
    insert_command_prefixes(mms, 1, 0xffffffff);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le32(&mms->outgoing_packet_data, 0);
    put_le_utf16(&mms->outgoing_packet_data, mms->path+1); // +1 because we skip the leading /
    put_le32(&mms->outgoing_packet_data, 0); /* More zeroes */

    err = send_command_packet(mms);
    ff_mms_set_state(mms, AWAITING_PASSWORD_QUERY_OR_MEDIA_FILE);
    return err;
}

/** Read incoming MMST media, header or command packet. */
static MMSSCPacketType get_tcp_server_response(MMSContext *mms)
{
    // read the 8 byte header...
    int read_result;
    MMSSCPacketType packet_type= -1;
    int done;

    // use url_fdopen & url_fclose...
    do {
        done= 1; // assume we're going to get a valid packet.
        if((read_result= get_buffer(&mms->incoming_io_buffer, mms->incoming_buffer, 8))==8) {
            // check if we are a command packet...
            if(AV_RL32(mms->incoming_buffer + 4)==0xb00bface) {
                mms->incoming_flags= mms->incoming_buffer[3];
                if((read_result= get_buffer(&mms->incoming_io_buffer, mms->incoming_buffer+8, 4)) == 4) {
                    int length_remaining= AV_RL32(mms->incoming_buffer+8) + 4;

#if (MMS_DEBUG_LEVEL>0)
                    fprintf(stderr, "Length remaining is %d\n", length_remaining);
#endif
                    // FIXME? ** VERIFY LENGTH REMAINING HAS SPACE
                    // read the rest of the packet....
                    read_result = get_buffer(&mms->incoming_io_buffer, mms->incoming_buffer + 12, length_remaining) ;
                    if (read_result == length_remaining) {
                        // we have it all; get the stuff out of it.
                        mms->incoming_buffer_length= length_remaining+12;

                        // get the packet type...
                        packet_type= AV_RL16(mms->incoming_buffer+36);

                    } else {
#if (MMS_DEBUG_LEVEL>0)
                        // read error...
                        fprintf(stderr, "3 read returned %d!\n", read_result);
#endif
                    }
                } else {
#if (MMS_DEBUG_LEVEL>0)
                    // read error...
                    fprintf(stderr, "2 read returned %d!\n", read_result);
#endif
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
                    fprintf(stderr, "Incoming Buffer Length exceeds buffer: %d>%d\n", mms->media_packet_buffer_length, (int) sizeof(mms->media_packet_incoming_buffer));
                }
                assert(mms->media_packet_buffer_length<sizeof(mms->media_packet_incoming_buffer));
                read_result= get_buffer(&mms->incoming_io_buffer, dst, length_remaining);
                if(read_result != length_remaining) {
#if (MMS_DEBUG_LEVEL>0)
                    fprintf(stderr, "read_bytes result: %d asking for %d\n", read_result, length_remaining);
#endif
                    break;
                } else {
                    // if we successfully read everything....
                    if(packet_id_type == mms->header_packet_id) {
                        // asf header
                        // fprintf(stderr, "asf header: %d\n", mms->incoming_buffer_length);
                        packet_type = SC_PACKET_ASF_HEADER_TYPE;

                        // Store the asf header
                        if(!mms->header_parsed) {
                            mms->asf_header = av_realloc(mms->asf_header, mms->asf_header_size + mms->media_packet_buffer_length);
                            memcpy(mms->asf_header + mms->asf_header_size, mms->media_packet_read_ptr, mms->media_packet_buffer_length);
                            mms->asf_header_size += mms->media_packet_buffer_length;
                        }
                    } else if(packet_id_type == mms->packet_id) {
                        //                    fprintf(stderr, "asf packet: %d\n", mms->incoming_buffer_length);
                        packet_type = SC_PACKET_ASF_MEDIA_TYPE;
                    } else {
#if (MMS_DEBUG_LEVEL>0)
                        fprintf(stderr, "packet id type %d which must be old, getting another one.", packet_id_type);
#endif
                        done= 0;
                    }
                }
            }
        } else {
            // read error...
            if(read_result<0) {
#if (MMS_DEBUG_LEVEL>0)
                fprintf(stderr, "Read error (or cancelled) returned %d!\n", read_result);
#endif
                packet_type = SC_PACKET_TYPE_CANCEL;
            } else {// 0 is okay, no data received.
#if (MMS_DEBUG_LEVEL>0)
                fprintf(stderr, "Read result of zero?!\n");
#endif
                packet_type = SC_PACKET_TYPE_NO_DATA;
            }
            done = 1;
        }
    } while(!done);

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
#if (MMS_DEBUG_LEVEL>0)
        fprintf(stderr, "Permission denied!\n");
#endif
        ff_mms_set_state(mms, STATE_ERROR);
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
#if (MMS_DEBUG_LEVEL>0)
        fprintf(stderr, "Broadcast flags: 0x%x\n", broadcast_flags); // 8000= allow index, 01= prerecorded, 02= live 42= presentation with script commands
        fprintf(stderr, "File Time Point?: %lld double size: %d double value: %lf\n", total_file_length_in_seconds, (int) sizeof(double), duration);
        fprintf(stderr, "Total in Seconds: %d\n", total_length_in_seconds);
        fprintf(stderr, "Packet length: %d\n", packet_length);
        fprintf(stderr, "Total Packet Count: %d\n", total_packet_count);
        fprintf(stderr, "Highest Bit Rate: %d\n", highest_bit_rate);
        fprintf(stderr, "Header Size: %d\n", header_size);
        fprintf(stderr, "---- Done ----\n");
#endif

        /* Disable seeking in live broadcasts for now */
        if(! (broadcast_flags & 0x0200))
            mms->seekable = 1;
    }
}

static int send_media_header_request(MMSContext *mms)
{
    int err= 0;
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

    err = send_command_packet(mms);
    ff_mms_set_state(mms, AWAITING_PACKET_HEADER_REQUEST_ACCEPTED_TYPE);
    return err;
}

static void handle_packet_stream_changing_type(MMSContext *mms)
{
    ByteIOContext pkt;
#if (MMS_DEBUG_LEVEL>0)
    fprintf(stderr, "Stream changing!\n");
#endif

    // read these from the incoming buffer.. (40 is the packet header size, without the prefixes)
    init_put_byte(&pkt, mms->incoming_buffer+40, mms->incoming_buffer_length-40, 0, NULL, NULL, NULL, NULL);
    get_le32(&pkt); // prefix 1
    mms->header_packet_id= (get_le32(&pkt) & 0xff); // prefix 2

    fprintf(stderr, "Changed header prefix to 0x%x", mms->header_packet_id);
    // mms->asf_header_length= 0;

    ff_mms_set_state(mms, AWAITING_ASF_HEADER); // this is going to hork our avstreams.
}

static int send_keepalive_packet(MMSContext *mms)
{
    // respond to a keepalive with a keepalive...
    start_command_packet(mms, CS_PACKET_KEEPALIVE_TYPE);
    insert_command_prefixes(mms, 1, 0x100FFFF);
    return send_command_packet(mms);
}

/** Handling pf TCP-specific packets.
 * @param packet_type incoming packet type to process.
 * @return  0 if the packet_type wasn't handled by this function.
 *          1 if it expected and handled.
 *         -1 if it the packet was unexpected in the current state.
 */
static int tcp_packet_state_machine(MMSContext *mms, MMSSCPacketType packet_type)
{
    switch(packet_type) {
    case SC_PACKET_CLIENT_ACCEPTED:
#if 0
        if(mms->state==AWAITING_SC_PACKET_CLIENT_ACCEPTED) {
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "Transitioning from AWAITING_SC_PACKET_CLIENT_ACCEPTED to AWAITING_SC_PACKET_TIMING_TEST_REPLY_TYPE\n");
#endif
            // send the timing request packet...
            start_command_packet(mms, CS_PACKET_TIMING_DATA_REQUEST_TYPE);
            insert_command_prefixes(mms, 0xf0f0f0f1, 0x0004000b);
            send_command_packet(mms);

            ff_mms_set_state(mms, AWAITING_SC_PACKET_TIMING_TEST_REPLY_TYPE);
        } else {
            return -1;
        }
        break;
#endif

    case SC_PACKET_TIMING_TEST_REPLY_TYPE: // we may, or may not have timing tests.
        if(mms->state==AWAITING_SC_PACKET_TIMING_TEST_REPLY_TYPE || mms->state==AWAITING_SC_PACKET_CLIENT_ACCEPTED) {
            send_protocol_select(mms);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_PROTOCOL_ACCEPTED_TYPE:
        if(mms->state==AWAITING_CS_PACKET_PROTOCOL_ACCEPTANCE) {
            send_media_file_request(mms);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_PROTOCOL_FAILED_TYPE:
        if(mms->state==AWAITING_CS_PACKET_PROTOCOL_ACCEPTANCE) {
            // abort;
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "Protocol failed\n");
#endif
            ff_mms_set_state(mms, STATE_ERROR);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_PASSWORD_REQUIRED_TYPE:
        if(mms->state==AWAITING_PASSWORD_QUERY_OR_MEDIA_FILE) {
            // we don't support this right now.
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "Password required\n");
#endif
            ff_mms_set_state(mms, STATE_ERROR);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_MEDIA_FILE_DETAILS_TYPE:
        if(mms->state==AWAITING_PASSWORD_QUERY_OR_MEDIA_FILE) {
            handle_packet_media_file_details(mms);
            if(mms->state != STATE_ERROR)
                send_media_header_request(mms);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_HEADER_REQUEST_ACCEPTED_TYPE:
        if(mms->state==AWAITING_PACKET_HEADER_REQUEST_ACCEPTED_TYPE) {
            // reset (in case we are doing this more than once)
//                mms->asf_header_length= 0;

            // wait for the header (follows immediately)
            ff_mms_set_state(mms, AWAITING_ASF_HEADER);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_STREAM_CHANGING_TYPE:
        if(mms->state==STREAMING) {
            handle_packet_stream_changing_type(mms);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_STREAM_ID_ACCEPTED_TYPE:
        if(mms->state==AWAITING_STREAM_ID_ACCEPTANCE) {
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "Stream ID's accepted!\n");
#endif
            ff_mms_set_state(mms, STREAM_PAUSED); // only way to get out of this is to play...
        } else {
            return -1;
        }
        break;

    case SC_PACKET_MEDIA_PACKET_FOLLOWS_TYPE:
        if(mms->state==AWAITING_STREAM_START_PACKET) {
            // get the stream packets...
            ff_mms_set_state(mms, STREAMING);
        } else {
            return -1;
        }
        break;

    case SC_PACKET_KEEPALIVE_TYPE:
        if(mms->state==STREAMING || mms->state==STREAM_PAUSED) {
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "Got a Keepalive!\n");
#endif
            send_keepalive_packet(mms);
        } else {
            return -1;
        }
        break;

    default:
        return 0; // Not handled here
    }

    return 1; // Handled here
}

/** Pad media packets smaller than max_packet_size and/or adjust read position
 * after a seek. */
static void pad_media_packet(MMSContext *mms)
{
    if(mms->media_packet_buffer_length<mms->asf_context.packet_obj_size) {
        int padding_size = mms->asf_context.packet_obj_size - mms->media_packet_buffer_length;
        //  fprintf(stderr, "Incoming packet smaller than the asf packet size stated (%d<%d) Padding.\n", mms->media_packet_buffer_length, mms->asf_context.packet_size);
        memset(mms->media_packet_incoming_buffer+mms->media_packet_buffer_length, 0, padding_size);
        mms->media_packet_buffer_length += padding_size;
    }

    if(mms->media_packet_seek_offset) {
        mms->media_packet_buffer_length -= mms->media_packet_seek_offset;
        mms->media_packet_read_ptr += mms->media_packet_seek_offset;
        mms->media_packet_seek_offset = 0;
    }
}

/** Single-step the packet-pumping state machine.
 * @return The type of the last packet from the server.
 */
static MMSSCPacketType ff_mms_packet_state_machine(MMSContext *mms)
{
    MMSSCPacketType packet_type = get_tcp_server_response(mms);

    /* First, try protocol-specific packet handling */
    int ret = tcp_packet_state_machine(mms, packet_type);
    if(ret != 0) {
        if(ret == -1)
            log_packet_in_wrong_state(mms, packet_type);
        return packet_type;
    }

    /* Common packet handling */
    switch(packet_type) {
    case SC_PACKET_ASF_HEADER_TYPE:
        if(mms->state==AWAITING_ASF_HEADER) {
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "Got a SC_PACKET_ASF_HEADER: %d\n", mms->media_packet_buffer_length);
#endif
            if((mms->incoming_flags == 0X08) || (mms->incoming_flags == 0X0C))
            {
#if (MMS_DEBUG_LEVEL>0)
                fprintf(stderr, "Got the full header!\n");
#endif
                ff_mms_set_state(mms, ASF_HEADER_DONE);
            }
        } else {
            log_packet_in_wrong_state(mms, packet_type);
        }
        break;

    case SC_PACKET_ASF_MEDIA_TYPE:
        if(mms->state==STREAMING || mms->state==AWAITING_PAUSE_ACKNOWLEDGE || mms->state==AWAITING_HTTP_PAUSE_CONTROL_ACKNOWLEDGE) {
//                fprintf(stderr, "Got a stream packet of length %d!\n", mms->incoming_buffer_length);
            pad_media_packet(mms);
        } else {
            log_packet_in_wrong_state(mms, packet_type);
        }
        break;

    case SC_PACKET_STREAM_STOPPED_TYPE:
        if(mms->state==AWAITING_PAUSE_ACKNOWLEDGE) {
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "Server echoed stream pause\n");
#endif
            ff_mms_set_state(mms, STREAM_PAUSED);
        } else if(mms->state==STREAMING) {
            /*
            When echoing a start (from me):
             receive command 0x1e, 48 bytes
             start sequence 00000001
             command id     b00bface
             length               20
             protocol       20534d4d
             len8                  4
             sequence #     00000006
             len8  (II)            2
             dir | comm     0004001e
             prefix1        00000000
             prefix2        ffff0100

             When Ending on it's own:
             receive command 0x1e, 48 bytes
             start sequence 09000001
             command id     b00bface
             length               20
             protocol       20534d4d
             len8                  4
             sequence #     00000006
             len8  (II)            2
             dir | comm     0004001e
             prefix1        00000000
             prefix2        00000004
            */
#if (MMS_DEBUG_LEVEL>0)
            fprintf(stderr, "** Server hit end of stream (may be sending new header information)\n");
#endif
            // TODO: if this is a live stream, on the resumption of a pause, this happens, then it follows with a SC_PACKET_STREAM_CHANGING_TYPE
            // otherwise it means this stream is done.
            ff_mms_set_state(mms, STREAM_DONE);
        } else {
            log_packet_in_wrong_state(mms, packet_type);
        }
        break;

    case SC_PACKET_TYPE_CANCEL:
        fprintf(stderr, "Got a -1 packet type\n");
        // user cancelled; let us out so it gets closed down...
        ff_mms_set_state(mms, USER_CANCELLED);
        break;

    case SC_PACKET_TYPE_NO_DATA:
#if (MMS_DEBUG_LEVEL>0)
        fprintf(stderr, "Got no data (closed?)\n");
#endif
        ff_mms_set_state(mms, STREAM_DONE); //?
        break;

    case SC_PACKET_HTTP_CONTROL_ACKNOWLEDGE:
        ff_mms_set_state(mms, AWAITING_PAUSE_ACKNOWLEDGE);
        break;

    default:
        fprintf(stderr, "Unhandled packet type %d\n", packet_type);
        break;
    }

    return packet_type;
}

/** Send the initial handshake. */
static int send_startup_packet(MMSContext *mms)
{
    int err;
    char data_string[256];

    //snprintf(data_string, sizeof(data_string), "NSPlayer/7.0.0.1956; {%s}; Host: %s",
    //        my_guid, mms->host);
    snprintf(data_string, sizeof(data_string), "NSPlayer/7.0.0.1956; {%s}",
            "7E667F5D-A661-495E-A512-F55686DDA178");

    start_command_packet(mms, CS_PACKET_INITIAL_TYPE);
    insert_command_prefixes(mms, 0, 0x0004000b);
    put_le32(&mms->outgoing_packet_data, 0x0003001c);
    put_le_utf16(&mms->outgoing_packet_data, data_string);
//    put_le16(&mms->outgoing_packet_data, 0); // double unicode ended string...
    put_le32(&mms->outgoing_packet_data, 0); // double unicode ended string...

    err = send_command_packet(mms);
    ff_mms_set_state(mms, AWAITING_SC_PACKET_CLIENT_ACCEPTED);
    return err;
}

/** Read the whole mms header into a buffer of our own .*/
static int read_mms_header(MMSContext *mms)
{
    if(mms->state != AWAITING_ASF_HEADER) {
        fprintf(stderr, "cannot read header this state\n");
        ff_mms_set_state(mms, STATE_ERROR);
        return -1;
    }

    /* TODO: add timeout */
    /* This will run until the header is stored */
    while(mms->state != ASF_HEADER_DONE && mms->state != STATE_ERROR)
        ff_mms_packet_state_machine(mms);

    if(mms->state == STATE_ERROR)
        return -1;

    /* Parse the header */
    init_put_byte(&mms->private_av_format_ctx.pb, mms->asf_header, mms->asf_header_size, 0, NULL, NULL, NULL, NULL);
    mms->private_av_format_ctx.priv_data = &mms->asf_context;

    if(asf_demuxer.read_header(&mms->private_av_format_ctx, NULL) < 0) {
        fprintf(stderr, "read_header failed\n");
        return -1;
    }

    mms->av_format_ctx = &mms->private_av_format_ctx; // Default
    mms->header_parsed = 1;

    return 0;
}

/** Clear all buffers of partial and old packets after a seek or other discontinuity */
void clear_stream_buffers(MMSContext *mms)
{
    mms->incoming_io_buffer.buf_ptr = mms->incoming_io_buffer.buf_end;
    mms->media_packet_buffer_length = 0;
    mms->media_packet_read_ptr = mms->media_packet_incoming_buffer;
}

/** Convert from AVDISCARD_* values to MMS stream selection code for the stream. */
static int ff_mms_stream_selection_code(AVStream *st)
{
    switch(st->discard) {
    case AVDISCARD_NONE:
    case AVDISCARD_DEFAULT:
    default:
        return 0; // 00 = stream at full frame rate;

    case AVDISCARD_NONREF:
    case AVDISCARD_BIDIR:
    case AVDISCARD_NONKEY:
        return 1; // 01 = only stream key frames (doesn't work that well)

    case AVDISCARD_ALL:
        return 2; // 02 = no stream, switch it off.
    }
}


/** Send MMST stream selection command based on the AVStream->discard values. */
static int send_stream_selection_request(MMSContext *mms)
{
    int ii;
    int err;

    //  send the streams we want back...
    start_command_packet(mms, CS_PACKET_STREAM_ID_REQUEST_TYPE);
    put_le32(&mms->outgoing_packet_data, mms->av_format_ctx->nb_streams);

    for(ii= 0; ii<mms->av_format_ctx->nb_streams; ii++) {
        AVStream *st= mms->av_format_ctx->streams[ii];

        put_le16(&mms->outgoing_packet_data, 0xffff); // flags
        put_le16(&mms->outgoing_packet_data, st->id); // stream id
        put_le16(&mms->outgoing_packet_data, ff_mms_stream_selection_code(st)); // selection
    }

    put_le16(&mms->outgoing_packet_data, 0); /* Extra zeroes */

    err = send_command_packet(mms);
    ff_mms_set_state(mms, AWAITING_STREAM_ID_ACCEPTANCE);
    return err;
}

/** Request a MMST stream from a timestamp or packet offset.
 * @param rate Play rate. (Only 1 supported as for now)
 * @param byte_offset Byte position to seek to. Set to -1 when seeking by timestamp.
 *  The position is from the start of the media stream, i.e. not counting header size.
 * @param timestamp Time point in ms. Set to 0 when seeking from packet offsets.
 */
static int request_streaming_from(MMSContext *mms,
        int rate, int64_t byte_offset, int64_t timestamp)
{
    int32_t packet = -1;
    int result;

    if(byte_offset > 0)
        packet = byte_offset / mms->asf_context.packet_obj_size;

    /* Send a stream selection request if this is the first call to play */
    if(mms->state == ASF_HEADER_DONE) {
        result = send_stream_selection_request(mms);

        if(result==0) {
            while(mms->state != STREAM_PAUSED && mms->state != STATE_ERROR && mms->state != STREAM_DONE) {
                ff_mms_packet_state_machine(mms);
            }
        }
    }

    if(mms->state==STREAM_PAUSED || mms->state == ASF_HEADER_DONE) {
        timestamp= av_dbl2int((double)timestamp/1000.0); // is this needed?

        start_command_packet(mms, CS_PACKET_START_FROM_PACKET_ID_TYPE);
        insert_command_prefixes(mms, 1, mms->packet_id);
        put_le64(&mms->outgoing_packet_data, timestamp); // seek timestamp
        put_le32(&mms->outgoing_packet_data, 0xffffffff);  // unknown
        put_le32(&mms->outgoing_packet_data, packet);    // packet offset
        put_byte(&mms->outgoing_packet_data, 0xff); // max stream time limit
        put_byte(&mms->outgoing_packet_data, 0xff); // max stream time limit
        put_byte(&mms->outgoing_packet_data, 0xff); // max stream time limit
        put_byte(&mms->outgoing_packet_data, 0x00); // stream time limit flag

        mms->packet_id++; // new packet_id so we can separate new data from old data
        put_le32(&mms->outgoing_packet_data, mms->packet_id);
        send_command_packet(mms);

        ff_mms_set_state(mms, AWAITING_STREAM_START_PACKET);
        return 0;
    } else {
#if (MMS_DEBUG_LEVEL>0)
//        fprintf(stderr, "Tried a read_play when the state was not stream paused (%s)\n", state_names[mms->state]);
#endif
        return -1;
    }
}

/** Read at most one media packet (or a whole header). */
static int read_mms_packet(MMSContext *mms, uint8_t *buf, int buf_size)
{
    int result = 0;
    MMSSCPacketType packet_type;
    int size_to_copy;

//    fprintf(stderr, "mms_read_packet()\n");

//    fprintf(stderr, "*** read packet %p needs %d bytes at %lld...\n", buf, buf_size, url_ftell(&mms->av_format_ctx->pb));
    if(mms->state != STREAM_DONE && mms->state != STREAM_PAUSED && mms->state != STATE_ERROR) {
        do {
            if(mms->asf_header_read_pos < mms->asf_header_size) {
                /* Read from ASF header buffer */
                size_to_copy= FFMIN(buf_size, mms->asf_header_size - mms->asf_header_read_pos);
                memcpy(buf, mms->asf_header + mms->asf_header_read_pos, size_to_copy);
                mms->asf_header_read_pos += size_to_copy;
                result += size_to_copy;
#if (MMS_DEBUG_LEVEL > 0)
                fprintf(stderr, "Copied %d bytes from stored header. left: %d\n", size_to_copy, mms->asf_header_size - mms->asf_header_read_pos);
#endif
            } else if(mms->media_packet_buffer_length) {
                /* Read from media packet buffer */
                size_to_copy = FFMIN(buf_size, mms->media_packet_buffer_length);
                memcpy(buf, mms->media_packet_read_ptr, size_to_copy);
                mms->media_packet_buffer_length -= size_to_copy;
                mms->media_packet_read_ptr+= size_to_copy;
                result += size_to_copy;
                //fprintf(stderr, "Copied %d bytes from media_packet read pointer! (result: %d, left: %d)\n", size_to_copy, result, mms->media_packet_buffer_length);
            } else {
                /* Read from network */
//               fprintf(stderr, "Calling state machine...\n");
                packet_type= ff_mms_packet_state_machine(mms);
//                fprintf(stderr, "Type: 0x%x\n", packet_type);
                switch (packet_type) {
                case SC_PACKET_ASF_MEDIA_TYPE:
                    if(mms->media_packet_buffer_length>mms->asf_context.packet_obj_size) {
                        fprintf(stderr, "Incoming packet larger than the asf packet size stated (%d>%d)\n", mms->media_packet_buffer_length, mms->asf_context.packet_obj_size);
                        result= AVERROR_IO;
                        break;
                    }

                    // copy the data to the packet buffer...
                    size_to_copy= FFMIN(buf_size, mms->media_packet_buffer_length);
                    memcpy(buf, mms->media_packet_read_ptr, size_to_copy);
                    mms->media_packet_buffer_length -= size_to_copy;
                    mms->media_packet_read_ptr += size_to_copy;
                    result += size_to_copy;
//fprintf(stderr, "Copied %d bytes (Media) from media_packet read pointer! (result: %d)\n", size_to_copy, result);
                    break;
                case SC_PACKET_ASF_HEADER_TYPE:
                    // copy the data to the packet buffer...
                    size_to_copy= FFMIN(buf_size, mms->media_packet_buffer_length);
                    memcpy(buf, mms->media_packet_read_ptr, size_to_copy);
                    mms->media_packet_buffer_length -= size_to_copy;
                    mms->media_packet_read_ptr+= size_to_copy;
                    result+= size_to_copy;
//fprintf(stderr, "Copied %d bytes (header) from media_packet read pointer! (result: %d)\n", size_to_copy, result);
                    break;
                default:
                    if(mms->state==STREAM_PAUSED) {
                        result= 0;
                    } else if(mms->state==STREAM_DONE || mms->state==USER_CANCELLED) {
                        result=-1;
                    } else if(mms->state==AWAITING_ASF_HEADER) {
                        // we have reset the header; spin though the loop..
//                        fprintf(stderr, "****-- Spinning the loop!\n");
                        while(mms->state != STREAMING && mms->state != STATE_ERROR) {
                            ff_mms_packet_state_machine(mms);
                        }
//                        fprintf(stderr, "****-- Done Spinning the loop!\n");
                    } else {
#if (MMS_DEBUG_LEVEL>0)
                        fprintf(stderr, "Got a packet in odd state: %s Packet Type: 0x%x\n", state_names[mms->state], packet_type);
#endif
                    }
                    break;
                }
            }
        } while(result==0 && mms->state!=STREAM_PAUSED); // only return one packet...
    } else {
        if(mms->state==STREAM_PAUSED) {
            result= 0;
        } else {
            result= -1;
        }
    }
//    fprintf(stderr, "read packet %p needs %d bytes.  getting %d\n", buf, buf_size, result);

    return result;
}

static int send_close_packet(MMSContext *mms)
{
    int err;
    start_command_packet(mms, CS_PACKET_STREAM_CLOSE_TYPE);
    insert_command_prefixes(mms, 1, 1);

    err = send_command_packet(mms);
    ff_mms_set_state(mms, AWAITING_PAUSE_ACKNOWLEDGE);
    return err;
}

/** Close the MMSH/MMST connection */
static int mms_close(URLContext *h)
{
    MMSContext *mms = (MMSContext *)h->priv_data;

#if (MMS_DEBUG_LEVEL>0)
    fprintf(stderr, "mms_close\n");
#endif
    if(mms->mms_hd) {
        // send the close packet if we should...
        if(mms->state != STATE_ERROR) {
            send_close_packet(mms);
        }

        // need an url_fdclose()
        close_connection(mms);
    }
#if (MMS_DEBUG_LEVEL>0)
    fprintf(stderr, "done with mms_close\n");
#endif

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

//    mms->protocol = &mmst_mmsprotocol;

    if(mms->port<0)
        mms->port = DEFAULT_MMS_PORT;
//#if (MMS_DEBUG_LEVEL>0)
//    fprintf(stderr, "Opening %s\n Auth: %s\n Host: %s\n Port: %d\n Path: %s\n", mms->location, authorization, mms->host, mms->port, mms->path);
//#endif

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

    /* okay, now setup stuff.  working from unclear specifications is great! */
    //mms->protocol->send_startup_message(mms); //TOBEDONE
    send_startup_packet(mms);

    // TODO: add timeout here... (mms_watchdog_reset() ?)
    while(mms->state != AWAITING_ASF_HEADER &&
            mms->state != STATE_ERROR &&
            mms->state != STREAM_DONE) {
        ff_mms_packet_state_machine(mms); // TOBEDONE
    }

    /* We store the header internally */
    if(mms->state != STATE_ERROR)
        read_mms_header(mms); // TOBEDONE use asf_read_header?

    if(mms->state == STATE_ERROR) {
        err = AVERROR(EIO);
        goto fail;
    }

#if (MMS_DEBUG_LEVEL>0)
    fprintf(stderr, "Leaving open (success)\n");
#endif

    return 0;
fail:
    mms_close(h);

#if (MMS_DEBUG_LEVEL>0)
    fprintf(stderr, "Leaving open (failure: %d)\n", err);
#endif

    return err;
}

static int mms_open(URLContext *h, const char *uri, int flags)
{
    MMSContext *mms;
    int ret;

    h->is_streamed = 1;
    mms = av_malloc(sizeof(MMSContext));
    if (!mms)
        return AVERROR(ENOMEM);
    memset(mms, 0, sizeof(MMSContext));
    h->priv_data = mms;
    av_strlcpy(mms->location, uri, MMS_URL_SIZE);

    ret = mms_open_cnx(h);
    return ret;
}

/** Like AVInputFormat::read_play().
 * @see AVInputFormat::read_play()
 */
static int ff_mms_play(URLContext *h)
{
    MMSContext *mms = h->priv_data;
    int64_t stream_offset = 0;

    /* Early return if already playing. */
    if(mms->state == STREAMING || mms->state == AWAITING_STREAM_START_PACKET)
        return 0;

#if (MMS_DEBUG_LEVEL>0)
    fprintf(stderr, "MMS Play\n");
#endif

    /* If resuming from pause */
    if(mms->state==STREAM_PAUSED)
        stream_offset = (mms->pause_resume_seq+1) * mms->asf_context.packet_obj_size;

    clear_stream_buffers(mms);

    return request_streaming_from(mms, 1, stream_offset, 0);
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
    if(mms->state == ASF_HEADER_DONE && mms->asf_header_read_pos >= mms->asf_header_size) {
        fprintf(stderr, "mms_read() before play(). Playing automatically.\n");
        result = ff_mms_play(h);// TOBEDONE
        if(result < 0)
            return result;
    }

    /* We won't get any packets from the server if paused. Nothing else to do than
     * to return. FIXME: return AVERROR(EAGAIN)? */
    if(mms->state == STREAM_PAUSED) {
        fprintf(stderr, "mms_read in STREAM_PAUSED\n");
        return 0;
    } else if(mms->state==STREAMING || mms->state==AWAITING_STREAM_START_PACKET || mms->state == ASF_HEADER_DONE) {
        result = read_mms_packet(mms, buf, size);//TOBEDONE use asf_read_packet?

        /* Note which packet we last returned. FIXME: doesn't handle partially read packets */
        mms->pause_resume_seq = mms->incoming_packet_seq;
    } else {
#if (MMS_DEBUG_LEVEL>0)
        fprintf(stderr, "mms_read: wrong state %s, returning AVERROR_IO!\n", state_names[mms->state]);
#endif
        result = AVERROR(EIO);
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
