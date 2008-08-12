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
 * @file lowpass.c
 * lowpass filter implementation
 */

#include "lowpass.h"

/**
 * filter data for 4th order IIR lowpass Butterworth filter
 *
 * data format:
 * normalized cutoff frequency | inverse filter gain | coefficients
 */
static const float lp_filter_data[][LOWPASS_FILTER_ORDER+2] = {
    { 0.5000000000, 9.398085e-01, -0.0176648009,  0.0000000000, -0.4860288221,  0.0000000000 },
    { 0.4535147392, 6.816645e-01, -0.4646665999, -2.2127207402, -3.9912017501, -3.2380429984 },
    { 0.4166666667, 4.998150e-01, -0.2498216698, -1.3392807613, -2.7693097862, -2.6386277439 },
    { 0.3628117914, 3.103469e-01, -0.0965076902, -0.5977763360, -1.4972580903, -1.7740085241 },
    { 0.3333333333, 2.346995e-01, -0.0557639007, -0.3623690447, -1.0304538354, -1.3066051440 },
    { 0.2916666667, 1.528432e-01, -0.0261686639, -0.1473794606, -0.6204721225, -0.6514716536 },
    { 0.2267573696, 6.917529e-02, -0.0202414073,  0.0780167640, -0.5277442247,  0.3631641670 },
    { 0.2187500000, 6.178391e-02, -0.0223681543,  0.1069446609, -0.5615167033,  0.4883976841 },
    { 0.2083333333, 5.298685e-02, -0.0261686639,  0.1473794606, -0.6204721225,  0.6514716536 },
    { 0.1587301587, 2.229030e-02, -0.0647354087,  0.4172275190, -1.1412129810,  1.4320761385 },
    { 0.1458333333, 1.693903e-02, -0.0823177861,  0.5192354923, -1.3444768251,  1.6365345642 },
    { 0.1133786848, 7.374053e-03, -0.1481421788,  0.8650973862, -1.9894244796,  2.1544844308 },
    { 0.1041666667, 5.541768e-03, -0.1742301048,  0.9921936565, -2.2090801108,  2.3024482658 },
};

int ff_lowpass_filter_init_coeffs(LPFilterCoeffs *coeffs, int freq, int cutoff)
{
    int i, j, size;
    float cutoff_ratio;

    //since I'm too lazy to calculate coefficients, I take more or less matching ones from the table
    //TODO: generic version
    size = sizeof(lp_filter_data) / sizeof(lp_filter_data[0]);
    cutoff_ratio = (float)cutoff / freq;
    if(cutoff_ratio > lp_filter_data[0][0])
        return -1;
    for(i = 0; i < size; i++){
        if(cutoff_ratio >= lp_filter_data[i][0])
            break;
    }
    if(i == size)
        i = size - 1;
    coeffs->gain     = lp_filter_data[i][1];
    for(j = 0; j < 4; j++)
        coeffs->c[j] = lp_filter_data[i][j+2];
    return 0;
}

