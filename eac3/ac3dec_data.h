/*
 * AC3 and E-AC3 decoder tables
 * Copyright (c) 2007 Bartlomiej Wolowiec <bartek.wolowiec@gmail.com>
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

#ifndef FFMPEG_AC3DEC_DATA_H
#define FFMPEG_AC3DEC_DATA_H

#include "common.h"

extern const uint8_t  ff_eac3_hebap_tab[64];
extern const uint8_t ff_eac3_bits_vs_hebap[20];
extern const int16_t ff_eac3_gaq_remap_1[12];
extern const int16_t ff_eac3_gaq_remap_2_4_a[9][2];
extern const int16_t ff_eac3_gaq_remap_2_4_b[9][2];

extern const int16_t (*ff_eac3_vq_hebap[8])[6];
extern const uint8_t ff_eac3_frm_expstr[32][6];
extern const uint8_t ff_eac3_default_cpl_band_struct[18];
extern const uint8_t ff_eac3_defspxbndstrc[17];
extern const uint8_t ff_eac3_defecplbndstrc[22];

extern const uint8_t ff_ac3_rematrix_band_tab[5];

extern const uint16_t ff_eac3_default_chmap[8];

/** Custom channel map locations bitmask
 *  Other channels described in documentation:
 *      Lc/Rc pair, Lrs/Rrs pair, Ts, Lsd/Rsd pair,
 *      Lw/Rw pair, Lvh/Rvh pair, Cvh, Reserved, LFE2
 */
enum CustomChannelMapLocation{
    AC3_CHMAP_L=        1<<(15-0),
    AC3_CHMAP_C=        1<<(15-1),
    AC3_CHMAP_R=        1<<(15-2),
    AC3_CHMAP_L_SUR=    1<<(15-3),
    AC3_CHMAP_R_SUR =   1<<(15-4),
    AC3_CHMAP_C_SUR=    1<<(15-7),
    AC3_CHMAP_LFE =     1<<(15-15)
};

#endif /* FFMPEG_AC3DEC_DATA_H */
