/*
 * Normal distribution PRNG following Marsaglia and Tang's Ziggurat algorithm
 * Copyright (c) 2004 David Bateman
 *
 * The algorithm is described the article:
 * Marsaglia, George and Tsang, Wai Wan (2000), The Ziggurat Method for
 * Generating Random Variables, Journal of Statistical Software, Vol. 5,
 * Issue 8, Oct 2000.
 * http://www.jstatsoft.org/v05/i08/
 *
 * The original version written by David Bateman has been modified here to
 * only include the normal distribution, use floats instead of doubles, and
 * return signed 4.24-bit fixed-point integers.
 *
 * This file is part of FFmpeg.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. The names of its contributors may not be used to endorse or promote
 *      products derived from this software without specific prior written
 *      permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "common.h"

#include "znrng.h"

#define NMANTISSA 0x80000000

#define TAIL_CUTOFF 3.6541528853610088
#define SLICE_AREA  0.00492867323399

static int cutoff_area_ratio[256];
static float slice_width[256];
static float slice_y_min[256];

/** get a uniformly-distributed point between 0.0 and 1.0 */
static inline float random_float(AVRandomState *c)
{
    return ((av_random(c) + 0.5) / 4294967296.0);
}

/** Gaussian probabilty function */
static inline float gpf(float x)
{
    return expf(-0.5*x*x);
}

/** inverse Gaussian probabilty function */
static inline float igpf(float x)
{
    return sqrtf(-2.0*logf(x));
}

void av_znrng_init(AVRandomState *c, unsigned int seed)
{
    int i;
    float w, w0;

    av_init_random(seed, c);

    w0 = TAIL_CUTOFF;
    slice_width[255] = w0 / NMANTISSA;
    slice_y_min[255] = gpf(w0);
    cutoff_area_ratio[0] = w0 * slice_y_min[255] / SLICE_AREA * NMANTISSA;
    slice_width[0] = SLICE_AREA / slice_y_min[255] / NMANTISSA;
    slice_y_min[0] = 1.0;

    for (i = 254; i > 0; i--) {
        w = igpf(SLICE_AREA / w0 + slice_y_min[i+1]);
        cutoff_area_ratio[i+1] = (int)(w / w0 * NMANTISSA);
        slice_width[i] = w / NMANTISSA;
        slice_y_min[i] = gpf(w);
        w0 = w;
    }
    cutoff_area_ratio[1] = 0;
}

int av_znrng_get(AVRandomState *c)
{
    float x;
    do {
        for (;;) {
            int ux = av_random(c);
            int slice = ux & 0xFF;
            x = ux * slice_width[slice];
            if (abs(ux) < cutoff_area_ratio[slice]) {
                break;
            } else if (slice == 0) {
                /* tail */
                float y;
                do {
                    x = -logf(random_float(c)) / TAIL_CUTOFF;
                    y = -logf(random_float(c));
                } while (y+y <= x*x);
                x = (ux < 0) ? -TAIL_CUTOFF-x : TAIL_CUTOFF+x;
                break;
            } else if ((slice_y_min[slice-1]-slice_y_min[slice])*
                    random_float(c)+slice_y_min[slice] < gpf(x)) {
                break;
            }
            /* if the point is above the curve, try again */
        }
    } while (x < -15.0 || x > 15.0);
    return lrintf(x * 0x8000000);
}


#ifdef TEST
#include "log.h"

#define TEST_COUNT 10000

int main(void)
{
    AVRandomState c;
    int x=0;
    int i, j;
    float mean, stddev;
    float data[TEST_COUNT]={0,};

    av_znrng_init(&c, 0xDEADBEEF);

    for (j = 0; j < TEST_COUNT; j++) {
        START_TIMER
        for (i = 0; i < 624; i++) {
            x+=av_znrng_get(&c);
        }
        STOP_TIMER("624 calls of av_znrng_get");
    }
    av_log(NULL, AV_LOG_ERROR, "final value:%X\n", x);

    for (j = 0; j < 16; j++) {
        mean = 0;
        for(i=0; i<TEST_COUNT; i++) {
            data[i] = (float)av_znrng_get(&c);
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
