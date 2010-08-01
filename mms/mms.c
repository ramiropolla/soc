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
#include "mms.h"
#include "asf.h"
#include "libavutil/intreadwrite.h"

int ff_asf_header_parser(MMSContext *mms)
{
    uint8_t *p = mms->asf_header;
    uint8_t *end;
    int flags, stream_id;
    mms->stream_num = 0;

    if (mms->asf_header_size < sizeof(ff_asf_guid) * 2 + 22 ||
        memcmp(p, ff_asf_header, sizeof(ff_asf_guid))) {
        av_log(NULL, AV_LOG_ERROR,
               "Corrupt stream (invalid ASF header, size=%d)\n",
               mms->asf_header_size);
        return AVERROR_INVALIDDATA;
    }

    end = mms->asf_header + mms->asf_header_size;

    p += sizeof(ff_asf_guid) + 14;
    while(end - p >= sizeof(ff_asf_guid) + 8) {
        uint64_t chunksize = AV_RL64(p + sizeof(ff_asf_guid));
        if (!chunksize || chunksize > end - p) {
            av_log(NULL, AV_LOG_ERROR,
                   "Corrupt stream (header chunksize %"PRId64" is invalid)\n",
                   chunksize);
            return AVERROR_INVALIDDATA;
        }
        if (!memcmp(p, ff_asf_file_header, sizeof(ff_asf_guid))) {
            /* read packet size */
            if (end - p > sizeof(ff_asf_guid) * 2 + 68) {
                mms->asf_packet_len = AV_RL32(p + sizeof(ff_asf_guid) * 2 + 64);
                if (mms->asf_packet_len <= 0 || mms->asf_packet_len > sizeof(mms->in_buffer)) {
                    av_log(NULL, AV_LOG_ERROR,
                           "Corrupt stream (too large pkt_len %d)\n",
                           mms->asf_packet_len);
                    return AVERROR_INVALIDDATA;
                }
            }
        } else if (!memcmp(p, ff_asf_stream_header, sizeof(ff_asf_guid))) {
            flags     = AV_RL16(p + sizeof(ff_asf_guid)*3 + 24);
            stream_id = flags & 0x7F;
            //The second condition is for checking CS_PKT_STREAM_ID_REQUEST packet size,
            //we can calcuate the packet size by stream_num.
            //Please see function send_stream_selection_request().
            if (mms->stream_num < MAX_STREAMS &&
                    46 + mms->stream_num * 6 < sizeof(mms->out_buffer)) {
                mms->streams[mms->stream_num].id = stream_id;
                mms->stream_num++;
            } else {
                av_log(NULL, AV_LOG_ERROR,
                       "Corrupt stream (too many A/V streams)\n");
                return AVERROR_INVALIDDATA;
            }
        } else if (!memcmp(p, ff_asf_head1_guid, sizeof(ff_asf_guid))) {
            chunksize = 46; // see references [2] section 3.4. This should be set 46.
        }
        p += chunksize;
    }

    return 0;
}
