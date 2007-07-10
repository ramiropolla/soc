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

#define MODE2_PATTERNS_NUM 20
/**
 * Intra types table
 *
 * These values are actually coded 3-tuples
 * used for detecting standard block configurations
 */
static const uint16_t rv40_aic_table_index[MODE2_PATTERNS_NUM] = {
 0x000, 0x100, 0x200,
 0x011, 0x111, 0x211, 0x511, 0x611,
 0x022, 0x122, 0x222, 0x722,
 0x272, 0x227,
 0x822, 0x282, 0x228,
 0x112, 0x116, 0x221
};

/**
 * Luma quantizer values
 * Second table is used for inter blocks
 */
static uint8_t rv40_luma_quant[2][32] = {
 {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
   16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 22, 22, 22, 22 },
 {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
   16, 17, 18, 19, 20, 20, 21, 21, 22, 23, 23, 23, 24, 24, 24, 24 }
};

/**
 * Chroma quantizer values
 * Second table is used for DC-only blocks
 */
static uint8_t rv40_chroma_quant[2][32] = {
 {  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
   16, 17, 17, 18, 19, 20, 20, 21, 22, 22, 23, 23, 24, 24, 25, 25 },
 {  0,  0,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
   14, 15, 15, 16, 17, 18, 18, 19, 20, 20, 21, 21, 22, 22, 23, 23 }
};
#endif /* RV40DATA_H */
