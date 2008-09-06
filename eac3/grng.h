/*
 * Gaussian PRNG using the Marsaglia polar method
 * Copyright (c) 2008 Justin Ruggles
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

#ifndef AVUTIL_GRNG_H
#define AVUTIL_GRNG_H

#include "lfg.h"

typedef struct {
    AVLFG state;
    int remaining;
    int y2;
} AVGRNG;

/**
 * Generates a random number, where the sequence has a standard normal
 * distribution (zero mean, unity standard deviation).
 * @return the next random number as a signed 4.28-bit fixed-point integer
 */
int av_grng_get(AVGRNG *c);

void av_grng_init(AVGRNG *c, unsigned int seed);

#endif /* AVUTIL_GRNG_H */
