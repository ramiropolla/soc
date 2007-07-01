/*
 * RealVideo 4 decoder
 * copyright (c) 2007 Konstantin Shishkov
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
 * @file rv40data.h
 * Miscellaneous RV40 tables.
 */

#ifndef RV40DATA_H
#define RV40DATA_H

#include <stdint.h>

/**
 * Number of ones in nibble minus one
 */
static const uint8_t rv40_count_ones[16] = {
    0, 0, 0, 1, 0, 1, 1, 2, 0, 1, 1, 2, 1, 2, 2, 3
};

/**
 * Values used to reconstruct coded block pattern
 */
static const uint8_t rv40_cbp_code[16] = {
    0x00, 0x20, 0x10, 0x30, 0x02, 0x22, 0x12, 0x32,
    0x01, 0x21, 0x11, 0x31, 0x03, 0x23, 0x13, 0x33
};

/**
 * Standard widths and heights coded in RV40
 */
//@{
static const int rv40_standard_widths[]   = { 160, 172, 240, 320, 352, 640, 704, 0};
static const int rv40_standard_heights[]  = { 120, 132, 144, 240, 288, 480, 0, 0};
static const int rv40_standard_heights2[] = { 180, 360, 576, 0};
//@}

/**
 * Precalculated results of division by three and modulo three for values 0-107
 *
 * A lot of four-tuples in RV40 are represented as c0*27+c1*9+c2*3+c3
 * This table allows conversion from value back to vector
 */
static const uint8_t modulo_three_table[108][4] = {
 { 0, 0, 0, 0 }, { 0, 0, 0, 1 }, { 0, 0, 0, 2 }, { 0, 0, 1, 0 },
 { 0, 0, 1, 1 }, { 0, 0, 1, 2 }, { 0, 0, 2, 0 }, { 0, 0, 2, 1 },
 { 0, 0, 2, 2 }, { 0, 1, 0, 0 }, { 0, 1, 0, 1 }, { 0, 1, 0, 2 },
 { 0, 1, 1, 0 }, { 0, 1, 1, 1 }, { 0, 1, 1, 2 }, { 0, 1, 2, 0 },
 { 0, 1, 2, 1 }, { 0, 1, 2, 2 }, { 0, 2, 0, 0 }, { 0, 2, 0, 1 },
 { 0, 2, 0, 2 }, { 0, 2, 1, 0 }, { 0, 2, 1, 1 }, { 0, 2, 1, 2 },
 { 0, 2, 2, 0 }, { 0, 2, 2, 1 }, { 0, 2, 2, 2 }, { 1, 0, 0, 0 },
 { 1, 0, 0, 1 }, { 1, 0, 0, 2 }, { 1, 0, 1, 0 }, { 1, 0, 1, 1 },
 { 1, 0, 1, 2 }, { 1, 0, 2, 0 }, { 1, 0, 2, 1 }, { 1, 0, 2, 2 },
 { 1, 1, 0, 0 }, { 1, 1, 0, 1 }, { 1, 1, 0, 2 }, { 1, 1, 1, 0 },
 { 1, 1, 1, 1 }, { 1, 1, 1, 2 }, { 1, 1, 2, 0 }, { 1, 1, 2, 1 },
 { 1, 1, 2, 2 }, { 1, 2, 0, 0 }, { 1, 2, 0, 1 }, { 1, 2, 0, 2 },
 { 1, 2, 1, 0 }, { 1, 2, 1, 1 }, { 1, 2, 1, 2 }, { 1, 2, 2, 0 },
 { 1, 2, 2, 1 }, { 1, 2, 2, 2 }, { 2, 0, 0, 0 }, { 2, 0, 0, 1 },
 { 2, 0, 0, 2 }, { 2, 0, 1, 0 }, { 2, 0, 1, 1 }, { 2, 0, 1, 2 },
 { 2, 0, 2, 0 }, { 2, 0, 2, 1 }, { 2, 0, 2, 2 }, { 2, 1, 0, 0 },
 { 2, 1, 0, 1 }, { 2, 1, 0, 2 }, { 2, 1, 1, 0 }, { 2, 1, 1, 1 },
 { 2, 1, 1, 2 }, { 2, 1, 2, 0 }, { 2, 1, 2, 1 }, { 2, 1, 2, 2 },
 { 2, 2, 0, 0 }, { 2, 2, 0, 1 }, { 2, 2, 0, 2 }, { 2, 2, 1, 0 },
 { 2, 2, 1, 1 }, { 2, 2, 1, 2 }, { 2, 2, 2, 0 }, { 2, 2, 2, 1 },
 { 2, 2, 2, 2 }, { 3, 0, 0, 0 }, { 3, 0, 0, 1 }, { 3, 0, 0, 2 },
 { 3, 0, 1, 0 }, { 3, 0, 1, 1 }, { 3, 0, 1, 2 }, { 3, 0, 2, 0 },
 { 3, 0, 2, 1 }, { 3, 0, 2, 2 }, { 3, 1, 0, 0 }, { 3, 1, 0, 1 },
 { 3, 1, 0, 2 }, { 3, 1, 1, 0 }, { 3, 1, 1, 1 }, { 3, 1, 1, 2 },
 { 3, 1, 2, 0 }, { 3, 1, 2, 1 }, { 3, 1, 2, 2 }, { 3, 2, 0, 0 },
 { 3, 2, 0, 1 }, { 3, 2, 0, 2 }, { 3, 2, 1, 0 }, { 3, 2, 1, 1 },
 { 3, 2, 1, 2 }, { 3, 2, 2, 0 }, { 3, 2, 2, 1 }, { 3, 2, 2, 2 },
};

#endif /* RV40DATA_H */
