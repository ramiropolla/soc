/*
 * AAC data declarations
 * Copyright (c) 2005-2006 Oded Shimon ( ods15 ods15 dyndns org )
 * Copyright (c) 2006-2007 Maxim Gavrilov ( maxim.gavrilov gmail com )
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
 * @file aactab.h
 * AAC data declarations
 * @author Oded Shimon  ( ods15 ods15 dyndns org )
 * @author Maxim Gavrilov ( maxim.gavrilov gmail com )
 */

#ifndef FFMPEG_AACTAB_H
#define FFMPEG_AACTAB_H

#include "libavutil/mem.h"
#include "aac.h"

#include <stdint.h>

DECLARE_ALIGNED(16, extern float,  ff_aac_kbd_long_1024[1024]);
DECLARE_ALIGNED(16, extern float,  ff_aac_kbd_short_128[128]);
DECLARE_ALIGNED(16, extern float, ff_aac_sine_long_1024[1024]);
DECLARE_ALIGNED(16, extern float, ff_aac_sine_short_128[128]);

extern const uint8_t ff_aac_num_swb_1024[];
extern const uint8_t ff_aac_num_swb_128 [];

extern const uint32_t ff_aac_scalefactor_code[121];
extern const uint8_t  ff_aac_scalefactor_bits[121];

extern const uint16_t *ff_aac_spectral_codes[11];
extern const uint8_t  *ff_aac_spectral_bits [11];
extern const uint16_t  ff_aac_spectral_sizes[11];

extern const int8_t *ff_aac_codebook_vectors[];

#ifdef CONFIG_HARDCODED_TABLES
extern const float ff_aac_ivquant_tab[IVQUANT_SIZE];
extern const float  ff_aac_pow2sf_tab[316];
#endif /* CONFIG_HARDCODED_TABLES */

#endif /* FFMPEG_AACTAB_H */
