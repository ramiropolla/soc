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
static const float scale_factors[25] = {
  0.000030517578125000000000000000000000000,
  0.000015258789062500000000000000000000000,
  0.000007629394531250000000000000000000000,
  0.000003814697265625000000000000000000000,
  0.000001907348632812500000000000000000000,
  0.000000953674316406250000000000000000000,
  0.000000476837158203125000000000000000000,
  0.000000238418579101562500000000000000000,
  0.000000119209289550781250000000000000000,
  0.000000059604644775390625000000000000000,
  0.000000029802322387695312500000000000000,
  0.000000014901161193847656250000000000000,
  0.000000007450580596923828125000000000000,
  0.000000003725290298461914062500000000000,
  0.000000001862645149230957031250000000000,
  0.000000000931322574615478515625000000000,
  0.000000000465661287307739257812500000000,
  0.000000000232830643653869628906250000000,
  0.000000000116415321826934814453125000000,
  0.000000000058207660913467407226562500000,
  0.000000000029103830456733703613281250000,
  0.000000000014551915228366851806640625000,
  0.000000000007275957614183425903320312500,
  0.000000000003637978807091712951660156250,
  0.000000000001818989403545856475830078125
};

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
