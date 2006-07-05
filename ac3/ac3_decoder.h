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

static const int8_t exp_1[128] = {
    -2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,-2,
    -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,
     0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
     2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
    25,25,25
};
static const int8_t exp_2[128] = {
    -2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    -2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    -2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    -2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    -2,-2,-2,-2,-2,-1,-1,-1,-1,-1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2,
    25,25,25
};
static const int8_t exp_3[128] = {
    -2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,
    -2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,
    -2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,
    -2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,
    -2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,-2,-1, 0, 1, 2,
    25,25,25
};

static const int16_t l3_q_tab[3] = { /* 3-level quantization table */
	(-2 << 15) / 3, 0, (2 << 15) / 3 
};

static const int16_t l5_q_tab[5] = { /* 5-level quantization table */
	(-4 << 15) / 5, (-2 << 15 ) / 5, 0, (2 << 15) / 5, (4 << 15) / 5
};

static const int16_t l7_q_tab[7] = { /* 7-level quantization table */
	(-6 << 15) / 7, (-4 << 15) / 7, (-2 << 15) / 7, 0,
	(2 << 15) / 7, (4 << 15) / 7, (6 << 15) / 7
};

static const int16_t l11_q_tab[11] = { /* 11-level quantization table */
	(-10 << 15) / 11, (-8 << 15) / 11, (-6 << 15) / 11, (-4 << 15) / 11, (-2 << 15) / 11, 0,
	(2 << 15) / 11, (4 << 15) / 11, (6 << 15) / 11, (8 << 15) / 11, (10 << 15) / 11
};

static const int16_t l15_q_tab[15] = { /* 15-level quantization table */
	(-14 << 15) / 15, (-12 << 15) / 15, (-10 << 15) / 15, (-8 << 15) / 15,
	(-6 << 15) / 15, (-4 << 15) / 15, (-2 << 15) / 15, 0,
	(2 << 15) / 15, (4 << 15) / 15, (6 << 15) / 15, (8 << 15) / 15,
	(10 << 15) / 15, (12 << 15) / 15, (14 << 15) / 15
};

static const uint8_t qntztab[16] = { 0, 5, 7, 3, 7, 4, 5, 6, 7, 8, 9, 10, 12, 12, 14, 16 };

static const float window[256] = {
        0.00014, 0.00024, 0.00037, 0.00051, 0.00067, 0.00086, 0.00107, 0.00130,
        0.00157, 0.00187, 0.00220, 0.00256, 0.00297, 0.00341, 0.00390, 0.00443,
        0.00501, 0.00564, 0.00632, 0.00706, 0.00785, 0.00871, 0.00962, 0.01061,
        0.01166, 0.01279, 0.01399, 0.01526, 0.01662, 0.01806, 0.01959, 0.02121,
        0.02292, 0.02472, 0.02662, 0.02863, 0.03073, 0.03294, 0.03527, 0.03770,
        0.04025, 0.04292, 0.04571, 0.04862, 0.05165, 0.05481, 0.05810, 0.06153,
        0.06508, 0.06878, 0.07261, 0.07658, 0.08069, 0.08495, 0.08935, 0.09389,
        0.09859, 0.10343, 0.10842, 0.11356, 0.11885, 0.12429, 0.12988, 0.13563,
        0.14152, 0.14757, 0.15376, 0.16011, 0.16661, 0.17325, 0.18005, 0.18699,
        0.19407, 0.20130, 0.20867, 0.21618, 0.22382, 0.23161, 0.23952, 0.24757,
        0.25574, 0.26404, 0.27246, 0.28100, 0.28965, 0.29841, 0.30729, 0.31626,
        0.32533, 0.33450, 0.34376, 0.35311, 0.36253, 0.37204, 0.38161, 0.39126,
        0.40096, 0.41072, 0.42054, 0.43040, 0.44030, 0.45023, 0.46020, 0.47019,
        0.48020, 0.49022, 0.50025, 0.51028, 0.52031, 0.53033, 0.54033, 0.55031,
        0.56026, 0.57019, 0.58007, 0.58991, 0.59970, 0.60944, 0.61912, 0.62873,
        0.63827, 0.64774, 0.65713, 0.66643, 0.67564, 0.68476, 0.69377, 0.70269,
        0.71150, 0.72019, 0.72877, 0.73723, 0.74557, 0.75378, 0.76186, 0.76981,
        0.77762, 0.78530, 0.79283, 0.80022, 0.80747, 0.81457, 0.82151, 0.82831,
        0.83496, 0.84145, 0.84779, 0.85398, 0.86001, 0.86588, 0.87160, 0.87716,
        0.88257, 0.88782, 0.89291, 0.89785, 0.90264, 0.90728, 0.91176, 0.91610,
        0.92028, 0.92432, 0.92822, 0.93197, 0.93558, 0.93906, 0.94240, 0.94560,
        0.94867, 0.95162, 0.95444, 0.95713, 0.95971, 0.96217, 0.96451, 0.96674,
        0.96887, 0.97089, 0.97281, 0.97463, 0.97635, 0.97799, 0.97953, 0.98099,
        0.98236, 0.98366, 0.98488, 0.98602, 0.98710, 0.98811, 0.98905, 0.98994,
        0.99076, 0.99153, 0.99225, 0.99291, 0.99353, 0.99411, 0.99464, 0.99513,
        0.99558, 0.99600, 0.99639, 0.99674, 0.99706, 0.99736, 0.99763, 0.99788,
        0.99811, 0.99831, 0.99850, 0.99867, 0.99882, 0.99895, 0.99908, 0.99919,
        0.99929, 0.99938, 0.99946, 0.99953, 0.99959, 0.99965, 0.99969, 0.99974,
        0.99978, 0.99981, 0.99984, 0.99986, 0.99988, 0.99990, 0.99992, 0.99993,
        0.99994, 0.99995, 0.99996, 0.99997, 0.99998, 0.99998, 0.99998, 0.99999,
        0.99999, 0.99999, 0.99999, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000,
        1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000 
};

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
