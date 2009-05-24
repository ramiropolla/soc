/*
 * ALS decoder
 * Copyright (c) 2009 Thilo Borgmann <thilo.borgmann _at_ googlemail.com>
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
 * @file libavcodec/alsdec.c
 * MPEG-4 ALS decoder
 * @author Thilo Borgmann <thilo.borgmann _at_ googlemail.com>
 */


#include "avcodec.h"


typedef struct {
} ALSDecContext;


/** Decodes an ALS frame.
 */
static av_cold int decode_frame(AVCodecContext *avctx,
                                void *data, int *data_size,
                                AVPacket *avpkt)
{
    return 0;
}


/** Initializes the ALS decoder.
 */
static av_cold int decode_init(AVCodecContext *avctx)
{
    return 0;
}


/** Uninitializes the ALS decoder.
 */
static av_cold int decode_end(AVCodecContext *avctx)
{
    return 0;
}


AVCodec als_decoder = {
    "als",
    CODEC_TYPE_AUDIO,
    CODEC_ID_MP4ALS,
    sizeof(ALSDecContext),
    decode_init,
    NULL,
    decode_end,
    decode_frame,
    .long_name = NULL_IF_CONFIG_SMALL("MPEG-4 Audio Lossless Coding"),
};

