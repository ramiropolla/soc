/*
 * AAC (LC) decoder
 * This code is developed as part of Google Summer of Code 2006 Program.
 *
 * Copyright (c) 2005 Oded Shimon( ods15 ods15 dyndns org )
 * Copyright (c) 2005-2006 Maxim Gavrilov ( maxim.gavrilov gmail com )
 * Copyright (c) 2007 Andreas Öman ( andreas lonelycoder com)
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef FFMPEG_AAC_H
#define FFMPEG_AAC_H

/**
 * Audio object types
 */
enum {
  AOT_NULL = 0x0,
  AOT_AAC_MAIN,
  AOT_AAC_LC,
  AOT_AAC_SSR,
  AOT_AAC_LTP,
  AOT_SBR,
  AOT_AAC_SCALABLE,
  AOT_TWINVQ,
  AOT_CELP,
  AOT_HVXC,
  AOT_TTSI = 12,
  AOT_MAINSYNTH,
  AOT_WAVESYNTH,
  AOT_MIDI,
  AOT_SAFX,
  AOT_ER_AAC_LC,
  AOT_ER_AAC_LTP = 19,
  AOT_ER_AAC_SCALABLE,
  AOT_ER_TWINVQ,
  AOT_ER_BSAC,
  AOT_ER_AAC_LD,
  AOT_ER_CELP,
  AOT_ER_HVXC,
  AOT_ER_HILN,
  AOT_ER_PARAM,
  AOT_SSC
};

/**
 * Syntactic elements
 * reference: Table 4.71
 */
enum {
    ID_SCE = 0x0,
    ID_CPE,
    ID_CCE,
    ID_LFE,
    ID_DSE,
    ID_PCE,
    ID_FIL,
    ID_END
};

/**
 * Window sequence types
 */
enum {
    ONLY_LONG_SEQUENCE = 0,
    LONG_START_SEQUENCE,
    EIGHT_SHORT_SEQUENCE,
    LONG_STOP_SEQUENCE
};

/**
 * Special codebooks
 */
#define ZERO_HCB 0
#define FIRST_PAIR_HCB 5
#define ESC_HCB 11
#define NOISE_HCB 13
#define INTENSITY_HCB2 14
#define INTENSITY_HCB 15
#define ESC_FLAG 16

#define TNS_MAX_ORDER 20

#endif /* FFMPEG_AACDEC_H */
