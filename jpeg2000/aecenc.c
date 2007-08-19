/*
 * Arithmetic entropy coder
 * Copyright (c) 2007 Kamil Nowosad
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
 * Arithmetic entropy coder
 * @file aecenc.c
 * @author Kamil Nowosad
 */

#include "aec.h"

static void byteout(AecState *aec)
{
retry:
    if (*aec->bp == 0xff){
        aec->bp++;
        *aec->bp = aec->c >> 20;
        aec->c &= 0xfffff;
        aec->ct = 7;
    } else if ((aec->c & 0x8000000)){
        (*aec->bp)++;
        aec->c &= 0x7ffffff;
        goto retry;
    } else{
        aec->bp++;
        *aec->bp = aec->c >> 19;
        aec->c &= 0x7ffff;
        aec->ct = 8;
    }
}

static void renorme(AecState *aec)
{
    do{
        aec->a += aec->a;
        aec->c += aec->c;
        if (!--aec->ct)
            byteout(aec);
    } while (!(aec->a & 0x8000));
}

static void setbits(AecState *aec)
{
    int tmp = aec->c + aec->a;
    aec->c |= 0xffff;
    if (aec->c >= tmp)
        aec->c -= 0x8000;
}

void ff_aec_initenc(AecState *aec, uint8_t *bp)
{
    ff_aec_init_contexts(aec);
    aec->a = 0x8000;
    aec->c = 0;
    aec->bp = bp-1;
    aec->bpstart = bp;
    aec->ct = 12 + (*aec->bp == 0xff);
}

void ff_aec_encode(AecState *aec, int cx, int d)
{
    int qe;

    aec->curcxstate = aec->cx_states + cx;
    qe = ff_aec_qe[*aec->curcxstate];
    aec->a -= qe;
    if (*aec->curcxstate & 1 == d){
        if (!(aec->a & 0x8000)){
            if (aec->a < qe)
                aec->a = qe;
            else
                aec->c += qe;
            *aec->curcxstate = ff_aec_nmps[*aec->curcxstate];
            renorme(aec);
        } else
            aec->c += qe;
    } else{
        if (aec->a < qe)
            aec->c += qe;
        else
            aec->a = qe;
        *aec->curcxstate = ff_aec_nlps[*aec->curcxstate];
        renorme(aec);
    }
}

int ff_aec_length(AecState *aec)
{
    return aec->bp - aec->bpstart;
}

int ff_aec_flush(AecState *aec)
{
    setbits(aec);
    aec->c = aec->c << aec->ct;
    byteout(aec);
    aec->c = aec->c << aec->ct;
    byteout(aec);
    if (*aec->bp != 0xff)
        aec->bp++;
    return aec->bp - aec->bpstart;
}
