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

#include "avcodec.h"
#include "lowpass2.h"
#include <complex.h>

void ff_lowpass_init(LPFilterContext *s, float sample_rate, float fpass, float fstop, float apass, float astop){
    int i, j;
    float fp = sample_rate*tanf(M_PI*fpass/sample_rate)/M_PI;
    float fs = sample_rate*tanf(M_PI*fstop/sample_rate)/M_PI;
    float ws = fs/fp;
    float vp = 2*M_PI*fp;
    int   N = (int)ceilf(log10f((pow(10, astop/10)-1) / (pow(10, apass/10)-1))/(2*log10f(ws)));
    float w0 = ws/pow(pow(10, astop/10)-1, 1.0/(2.0*N));
    float dfi0 = M_PI/N;
    complex *p, *pt;
    complex gain = 1.0;

    p = av_malloc(N*sizeof(*p));
    pt = av_malloc((N+1)*sizeof(*pt));
    for(i=0; i<2; i++){
        s->filterCoeffs[i] = av_malloc((N+1)*sizeof(*s->filterCoeffs[i]));
        s->buf[i] = av_malloc((N+1)*sizeof(*s->buf[i]));
    }
    for(i=0; i<N; i++){
        s->buf[0][i] = s->buf[1][i] = 0;
    }

    av_log(NULL, AV_LOG_DEBUG, "fp=%f fs=%f\n", fp, fs);
    av_log(NULL, AV_LOG_DEBUG, "vp=%f\n", vp);
    av_log(NULL, AV_LOG_DEBUG, "ws=%f\n", ws);
    av_log(NULL, AV_LOG_DEBUG, "N=%i w0=%f\n", N, w0);

    for(i=0; i<N; i++){
        p[i] = w0 * cexp(I*(M_PI/2.0 + (i+0.5)*dfi0));
        gain *= -p[i];
        p[i] *= vp;
        gain *= vp/(2.0*sample_rate-p[i]);
        p[i] = (2.0*sample_rate+p[i])/(2.0*sample_rate-p[i]);
        av_log(NULL, AV_LOG_DEBUG, "p[%i]=%f+%fI\n", i, creal(p[i]), cimag(p[i]));
    }

    av_log(NULL, AV_LOG_DEBUG, "gain=%f+%fI\n", creal(gain), cimag(gain));

    for(i=0; i<N; i++){
        pt[i] = 1;
        for(j=i; j>0; j--)
            pt[j] = -pt[j]*p[i] + pt[j-1];
        pt[0] *= -p[i];
    }
    for(i=0; i<=N; i++){
        av_log(NULL, AV_LOG_DEBUG, "a %i: %f\n", i, creal(pt[i]));
    }
    pt[N]=1;
    for(i=0; i<=N/2; i++){
        complex t;
        t=pt[i];
        pt[i] = pt[N-i];
        pt[N-i] = t;
    }
    for(i=0; i<=N; i++){
        av_log(NULL, AV_LOG_DEBUG, "%i: %f\n", i, creal(pt[i]));
    }

    for(i=0; i<N; i++)
        s->filterCoeffs[0][i] = creal(pt[i+1]);
    s->filterCoeffs[0][N] = 0;

    av_free(p);
    av_free(pt);

    for(i=0; i<N; i++){
        s->filterCoeffs[1][i] = gain;
        for(j=i; j>0; j--)
            s->filterCoeffs[1][j] = s->filterCoeffs[1][j] + s->filterCoeffs[1][j-1];
    }
    s->filterCoeffs[1][N] = gain;

    for(i=0; i<=N; i++){
        av_log(NULL, AV_LOG_DEBUG, "%i: ac=%f bc=%f\n", i, s->filterCoeffs[0][i], s->filterCoeffs[1][i]);
    }
    s->N = N;
}

void ff_lowpass_end(LPFilterContext *s){
    int i;
    for(i=0; i<2; i++){
        av_free(s->filterCoeffs[i]);
        av_free(s->buf[i]);
    }
}

void ff_lowpass_filter(LPFilterContext *s, int16_t *in, float *out, int n){
    int i, k;
    float tmp1, tmp2;

    //TODO cyclic buffer
    for(k=0; k<n; k++){
        for(i=s->N; i>0; i--)
            s->buf[0][i] = s->buf[0][i-1];

        s->buf[0][0] = in[k];
        tmp1 = tmp2 = 0;
        for(i=0; i<s->N; i++){
            tmp1 += (double)s->filterCoeffs[1][i]*(double)s->buf[0][i];
            tmp2 += (double)s->filterCoeffs[0][i]*(double)s->buf[1][i];
        }
        for(i=s->N; i>0; i--)
            s->buf[1][i] = s->buf[1][i-1];

        s->buf[1][0] = out[k] = tmp1-tmp2;
    }
}
