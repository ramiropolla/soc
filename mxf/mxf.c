/*
 * MXF
 * Copyright (c) 2006 SmartJog S.A., Baptiste Coudurier <baptiste dot coudurier at smartjog dot com>.
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
 * SMPTE RP224 http://www.smpte-ra.org/mdd/index.html
 */

#include "mxf.h"

const MXFDataDefinitionUL ff_mxf_data_definition_uls[] = {
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x01,0x03,0x02,0x02,0x01,0x00,0x00,0x00 }, CODEC_TYPE_VIDEO },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x01,0x03,0x02,0x02,0x02,0x00,0x00,0x00 }, CODEC_TYPE_AUDIO },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x05,0x01,0x03,0x02,0x02,0x02,0x02,0x00,0x00 }, CODEC_TYPE_AUDIO },
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  CODEC_TYPE_DATA },
};

const MXFCodecUL ff_mxf_codec_uls[] = {
    /* PictureEssenceCoding */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x01,0x11,0x00 }, 14, CODEC_ID_MPEG2VIDEO }, /* MP@ML Long GoP */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x01,0x02,0x02,0x01,0x02,0x01,0x01 }, 14, CODEC_ID_MPEG2VIDEO }, /* D-10 50Mbps PAL */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x03,0x03,0x00 }, 14, CODEC_ID_MPEG2VIDEO }, /* MP@HL Long GoP */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x04,0x02,0x00 }, 14, CODEC_ID_MPEG2VIDEO }, /* 422P@HL I-Frame */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x20,0x02,0x03 }, 14,      CODEC_ID_MPEG4 }, /* XDCAM proxy_pal030926.mxf */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x01,0x02,0x02,0x02,0x01,0x02,0x00 }, 13,    CODEC_ID_DVVIDEO }, /* DV25 IEC PAL */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x07,0x04,0x01,0x02,0x02,0x03,0x01,0x01,0x00 }, 14,   CODEC_ID_JPEG2000 }, /* JPEG2000 Codestream */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x01,0x02,0x01,0x7F,0x00,0x00,0x00 }, 13,   CODEC_ID_RAWVIDEO }, /* Uncompressed */
    /* SoundEssenceCompression */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x01,0x00,0x00,0x00,0x00 }, 13,  CODEC_ID_PCM_S16LE }, /* Uncompressed */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x01,0x7F,0x00,0x00,0x00 }, 13,  CODEC_ID_PCM_S16LE },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x07,0x04,0x02,0x02,0x01,0x7E,0x00,0x00,0x00 }, 13,  CODEC_ID_PCM_S16BE }, /* From Omneon MXF file */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x04,0x04,0x02,0x02,0x02,0x03,0x01,0x01,0x00 }, 15,   CODEC_ID_PCM_ALAW }, /* XDCAM Proxy C0023S01.mxf */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x02,0x03,0x02,0x01,0x00 }, 15,        CODEC_ID_AC3 },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x02,0x03,0x02,0x05,0x00 }, 15,        CODEC_ID_MP2 }, /* MP2 or MP3 */
  //{ { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x02,0x03,0x02,0x1C,0x00 }, 15,    CODEC_ID_DOLBY_E }, /* Dolby-E */
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  0,       CODEC_ID_NONE },
};

const MXFCodecUL ff_mxf_essence_container_uls[] = {
    // video essence container uls
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x02,0x0D,0x01,0x03,0x01,0x02,0x04,0x60,0x01 }, 14, CODEC_ID_MPEG2VIDEO }, /* MPEG-ES Frame wrapped */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x03,0x01,0x02,0x02,0x41,0x01 }, 14,    CODEC_ID_DVVIDEO }, /* DV 625 25mbps */
    // sound essence container uls
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x03,0x01,0x02,0x06,0x01,0x00 }, 14, CODEC_ID_PCM_S16LE }, /* BWF Frame wrapped */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x02,0x0D,0x01,0x03,0x01,0x02,0x04,0x40,0x01 }, 14,       CODEC_ID_MP2 }, /* MPEG-ES Frame wrapped, 0x40 ??? stream id */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x03,0x01,0x02,0x01,0x01,0x01 }, 14, CODEC_ID_PCM_S16LE }, /* D-10 Mapping 50Mbps PAL Extended Template */
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  0,      CODEC_ID_NONE },
};

