/*
 * Butterworth lowpass filter
 * This code is developed as part of Google Summer of Code 2008 Program.
 *
 * Copyright (c) 2008 Bartlomiej Wolowiec
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

#ifndef AVCODEC_LOWPASS_H
#define AVCODEC_LOWPASS_H

typedef struct LPFilterContext {
    float *filterCoeffs[2];
    float *buf[2];
    int N;
}LPFilterContext;

/**
 * Initialization of butterworth lowpass filter
 *
 * @param   s Lowpass filter context
 * @param   sample_rate Sample rate
 * @fpass   pass frequency
 * @fstop   stop frequency
 * @apass   distortion below pass frequency (dB)
 * @astop   stop frequency attenuation (dB)
 */
void ff_lowpass_init(LPFilterContext *s, float sample_rate, float fpass, float fstop, float apass, float astop);

void ff_lowpass_end(LPFilterContext *s);

void ff_lowpass_filter(LPFilterContext *s, int16_t *in, float *out, int n);

#endif /* AVCODEC_LOWPASS_H */
