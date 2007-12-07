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

#ifndef FFMPEG_AACTAB_H
#define FFMPEG_AACTAB_H

#include "common.h"

typedef struct {
    const uint16_t (*a)[2];
    const unsigned int s;
} aac_codebook;

extern const int          aac_sample_rates[16];
extern const aac_codebook aac_codebooks[12];
extern const unsigned int aac_scalefactor_huffman_table[121][2];
extern const uint16_t    *aac_swb_offset_1024[12];
extern const uint16_t    *aac_swb_offset_128[12];
extern const uint8_t      aac_num_swb_1024[12];
extern const uint8_t      aac_num_swb_128[12];
extern const uint8_t      aac_tns_max_bands_1024[12];
extern const uint8_t      aac_tns_max_bands_128[12];
extern const float       *aac_tns_coeffs_table[4];

#endif /* FFMPEG_AACTAB_H */
