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
#include <complex.h>
#include <math.h>

/**
 * IIR filter global parameters
 */
typedef struct FFLPFilterCoeffs{
    int   order;
    float gain;
    int   *cx;
    float *cy;
}FFLPFilterCoeffs;

/**
 * IIR filter state
 */
typedef struct FFLPFilterState{
    float *x;
}FFLPFilterState;

struct FFLPFilterCoeffs* ff_lowpass_filter_init_coeffs(int order, float cutoff_ratio)
{
    int i, j, size;
    FFLPFilterCoeffs *c;
    double wa;
    complex *p, *zp;

    if(order <= 1 || (order & 1) || cutoff_ratio >= 1.0)
        return NULL;

    c = av_malloc(sizeof(FFLPFilterCoeffs));
    c->cx = av_mallocz(sizeof(c->cx[0]) * (order + 1));
    c->cy = av_malloc (sizeof(c->cy[0]) * order);
    c->order = order;

    p  = av_malloc(sizeof(p[0])  * (order + 1));
    zp = av_malloc(sizeof(zp[0]) * (order + 1));

    wa = 2 * tan(M_PI * cutoff_ratio / 2.0);

    for(i = 0; i <= order; i++){
        complex t;
        double th = (i + order/2 + 0.5) * M_PI / order;

        t = (cos(th) + I*sin(th)) * wa;
        zp[i] = (2.0 + t) / (2.0 - t);
    }

    c->cx[0] = 1;
    for(i = 0; i <= order; i++)
        for(j = i; j >= 1; j--)
            c->cx[j] += c->cx[j - 1];

    p[0] = 1.0;
    for(i = 1; i <= order; i++)
        p[i] = 0.0;
    for(i = 0; i < order; i++){
        for(j = order; j >= 1; j--)
            p[j] = -zp[i]*p[j] + p[j - 1];
        p[0] *= -zp[i];
    }
    c->gain = creal(p[order]);
    for(i = 0; i < order; i++){
        complex t = -p[i] / p[order];
        c->gain += creal(p[i]);
        c->cy[i] = creal(t);
    }
    c->gain /= 1 << order;

    av_free(p);
    av_free(zp);

    return c;
}

struct FFLPFilterState* ff_lowpass_filter_init_state(int order)
{
    FFLPFilterState *s = av_mallocz(sizeof(FFLPFilterState));
    s->x = av_mallocz(sizeof(s->x[0]) * order);
    return s;
}

#define FILTER(i0, i1, i2, i3)                    \
    in =   *src * c->gain                         \
         + c->cy[0]*s->x[i0] + c->cy[1]*s->x[i1]  \
         + c->cy[2]*s->x[i2] + c->cy[3]*s->x[i3]; \
    res =  (s->x[i0] + in      )*1                \
         + (s->x[i1] + s->x[i3])*4                \
         +  s->x[i2]            *6;               \
    *dst = av_clip_int16(lrintf(res));            \
    s->x[i0] = in;                                \
    src += sstep;                                 \
    dst += dstep;                                 \

void ff_lowpass_filter(const struct FFLPFilterCoeffs *c, struct FFLPFilterState *s, int size, int16_t *src, int sstep, int16_t *dst, int dstep)
{
    int i;

    if(c->order == 4){
        for(i = 0; i < size; i += 4){
            float in, res;

            FILTER(0, 1, 2, 3);
            FILTER(1, 2, 3, 0);
            FILTER(2, 3, 0, 1);
            FILTER(3, 0, 1, 2);
        }
    }else{
        for(i = 0; i < size; i++){
            int j;
            float in, res;
            in = *src * c->gain;
            for(j = 0; j < c->order; j++)
                in += c->cy[j] * s->x[j];
            res = s->x[0] + in + s->x[c->order/2] * c->cx[c->order/2];
            for(j = 1; j < (c->order / 2); j++)
                res += (s->x[j] + s->x[c->order - j]) * c->cx[j];
            for(j = 0; j < c->order - 1; j++)
                s->x[j] = s->x[j + 1];
            *dst = av_clip_int16(lrintf(res));
            s->x[c->order - 1] = in;
            src += sstep;
            dst += sstep;
        }
    }
}

void ff_lowpass_filter_free_state(struct FFLPFilterState *state)
{
    if(state)
        av_free(state->x);
    av_free(state);
}

void ff_lowpass_filter_free_coeffs(struct FFLPFilterCoeffs *coeffs)
{
    if(coeffs){
        av_free(coeffs->cx);
        av_free(coeffs->cy);
    }
    av_free(coeffs);
}

