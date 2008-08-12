/*
 * Lowpass IIR filter
 * Copyright (c) 2008 Konstantin Shishkov
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
 * @file lowpass.h
 * lowpass filter interface
 */

#ifndef FFMPEG_LOWPASS_H
#define FFMPEG_LOWPASS_H

#include "avcodec.h"

/** filter order */
#define LOWPASS_FILTER_ORDER 4

/**
 * IIR filter global parameters
 */
typedef struct LPFilterCoeffs{
    float gain;
    float c[LOWPASS_FILTER_ORDER];
}LPFilterCoeffs;

/**
 * IIR filter state
 */
typedef struct LPFilterState{
    float x[LOWPASS_FILTER_ORDER + 1];
    float y[LOWPASS_FILTER_ORDER + 1];
}LPFilterState;

/**
 * Initialize filter coefficients.
 *
 * @param coeffs filter coefficients
 * @param freq   input frequency (sample rate/2)
 * @param cutoff cutoff frequency
 *
 * @return zero if filter creation succeeded, a negative value if filter could not be created
 */
int ff_lowpass_filter_init_coeffs(LPFilterCoeffs *coeffs, int freq, int cutoff);

/**
 * Filter input value.
 *
 * @param coeffs filter coefficients
 * @param s      filter state
 * @param in     input value
 *
 * @return filtered value
 */
static av_always_inline float ff_lowpass_filter(LPFilterCoeffs *coeffs, LPFilterState *s, float in)
{
    int i;
    for(i = 0; i < LOWPASS_FILTER_ORDER; i++){
        s->x[i] = s->x[i+1];
        s->y[i] = s->y[i+1];
    }
    s->x[LOWPASS_FILTER_ORDER] = in * coeffs->gain;
    //FIXME: made only for 4th order filter
    s->y[LOWPASS_FILTER_ORDER] = (s->x[0] + s->x[4])*1
                               + (s->x[1] + s->x[3])*4
                               +  s->x[2]           *6
                               + coeffs->c[0]*s->y[0] + coeffs->c[1]*s->y[1]
                               + coeffs->c[2]*s->y[2] + coeffs->c[3]*s->y[3];
    return s->y[LOWPASS_FILTER_ORDER];
}

#endif /* FFMPEG_LOWPASS_H */

