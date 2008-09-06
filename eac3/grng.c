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

#include <math.h>
#include "grng.h"

int av_grng_get(AVGRNG *c)
{
    float w, x1, y1, x2, y2;

    /* return 2nd value from previous run */
    if (c->remaining) {
        c->remaining = 0;
        return c->y2;
    }

    /* generate a uniformly distributed random point which falls within the
       unit circle */
    do {
        x1 = (int)av_lfg_get(&c->state) / (float)0x80000000;
        y1 = (int)av_lfg_get(&c->state) / (float)0x80000000;
        w = x1 * x1 + y1 * y1;
    } while (w >= 1.0);

    /* scale the point using the inverse gaussian probability function */
    w = sqrtf((-2.0f * logf(w)) / w);
    x2 = x1 * w;
    y2 = y1 * w;

    /* convert to signed 4.28-bit fixed-point.
       save 'y' value in state. return 'x' value. */
    c->y2 = lrintf(y2 * 0x8000000);
    c->remaining = 1;
    return lrintf(x2 * 0x8000000);
}

void av_grng_init(AVGRNG *c, unsigned int seed)
{
    av_lfg_init(&c->state, seed);
    c->remaining = 0;
}

#ifdef TEST
#include "log.h"
#include "common.h"

#define TEST_COUNT 10000

int main(void)
{
    int x=0;
    int i, j;
    AVGRNG state;
    float mean, stddev;
    float data[TEST_COUNT]={0,};

    av_grng_init(&state, 0xDEADBEEF);

    for (j = 0; j < TEST_COUNT; j++) {
        START_TIMER
        for (i = 0; i < 624; i++) {
            x+=av_grng_get(&state);
        }
        STOP_TIMER("624 calls of av_pgrng_get");
    }
    av_log(NULL, AV_LOG_ERROR, "final value:%X\n", x);

    for (j = 0; j < 16; j++) {
        mean = 0;
        for(i=0; i<TEST_COUNT; i++) {
            data[i] = (float)av_grng_get(&state);
            mean += data[i];
        }
        mean /= TEST_COUNT;
        stddev = 0;
        for(i=0; i<TEST_COUNT; i++) {
            float diff = data[i] - mean;
            stddev += (diff * diff);
        }
        stddev /= TEST_COUNT;
        stddev = sqrtf(stddev);
        av_log(NULL, AV_LOG_ERROR, "mean:%+f  stddev:%f\n", mean / 0x8000000,
               stddev / 0x8000000);
    }

    return 0;
}
#endif
