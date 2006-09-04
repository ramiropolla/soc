/* AC3 Audio Decoder.
 *
 * Copyright (c) 2006 Kartikey Mahendra BHATT (bhattkm at gmail dot com)
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#ifndef _AC3_DECODER_H
#define _AC3_DECODER_H
static const int nfchans_tbl[8] = { 2, 1, 2, 3, 3, 4, 4, 5 };

/* table for exponent to scale_factor mapping
 * scale_factor[i] = 2 ^ -(i + 15)
 */
static float scale_factors[25];

static int16_t psdtab[25];

static int8_t exp_1[128];
static int8_t exp_2[128];
static int8_t exp_3[128];

static int16_t l3_quantizers_1[32];
static int16_t l3_quantizers_2[32];
static int16_t l3_quantizers_3[32];

static int16_t l5_quantizers_1[128];
static int16_t l5_quantizers_2[128];
static int16_t l5_quantizers_3[128];

static int16_t l7_quantizers[7];

static int16_t l11_quantizers_1[128];
static int16_t l11_quantizers_2[128];

static int16_t l15_quantizers[15];

static const uint8_t qntztab[16] = { 0, 5, 7, 3, 7, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 16 };

/* Adjustmens in dB gain */
#define LEVEL_MINUS_3DB         0.7071067811865476
#define LEVEL_MINUS_4POINT5DB   0.5946035575013605
#define LEVEL_MINUS_6DB         0.5000000000000000
#define LEVEL_PLUS_3DB          1.4142135623730951
#define LEVEL_PLUS_6DB          2.0000000000000000
#define LEVEL_ZERO              0.0000000000000000

static const float clevs[4] = { LEVEL_MINUS_3DB, LEVEL_MINUS_4POINT5DB, 
    LEVEL_MINUS_6DB, LEVEL_MINUS_4POINT5DB };

static const float slevs[4] = { LEVEL_MINUS_3DB, LEVEL_MINUS_6DB, LEVEL_ZERO, LEVEL_MINUS_6DB };
#endif
