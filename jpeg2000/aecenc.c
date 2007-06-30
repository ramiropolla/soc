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
 *
 */

/**
 * Arithmetic entropy coder
 * @file aecenc.c
 * @author Kamil Nowosad
 */

#include "aec.h"

static void byteout_l(AecState *aec)
{
    aec->bp++;
    *aec->bp = aec->c >> 19;
    aec->c &= 0x7ffff;
    aec->ct = 8;
}

static void byteout_r(AecState *aec)
{
    aec->bp++;
    *aec->bp = aec->c >> 20;
    aec->c &= 0xfffff;
    aec->ct = 7;
}

static void byteout(AecState *aec)
{
    if (*aec->bp == 0xff){
        byteout_r(aec);
    }
    else
    {
        if ((aec->c & 0x8000000) == 0){
            byteout_l(aec);
        }
        else{
            (*aec->bp)++;
            if (*aec->bp == 0xff){
                aec->c &= 0x7ffffff;
                byteout_r(aec);
            }
            else{
                byteout_l(aec);
            }
        }
    }
}

static void renorme(AecState *aec)
{
    do{
        aec->a = aec->a << 1;
        aec->c = aec->c << 1;
        aec->ct--;
        if (!aec->ct)
            byteout(aec);
    }
    while ((aec->a & 0x8000) == 0);
}

static void codelps(AecState *aec)
{
    int qe = cx_states[aec->curctx->state].qe;
    aec->a -= qe;
    if (aec->a < qe)
        aec->c += qe;
    else
        aec->a = qe;
    if (cx_states[aec->curctx->state].sw)
        aec->curctx->mps = 1 - aec->curctx->mps;
    aec->curctx->state = cx_states[aec->curctx->state].nlps;
    renorme(aec);
}

static void codemps(AecState *aec)
{
    int qe = cx_states[aec->curctx->state].qe;
    aec->a -= qe;
    if ((aec->a & 0x8000) == 0){
        if (aec->a < qe)
            aec->a = qe;
        else
            aec->c += qe;
        aec->curctx->state = cx_states[aec->curctx->state].nmps;
        renorme(aec);
    }
    else
        aec->c += qe;
}

static void setbits(AecState *aec)
{
    int tmp = aec->c + aec->a;
    aec->c |= 0xffff;
    if (aec->c >= tmp)
        aec->c -= 0x8000;
}

void aec_initenc(AecState *aec, uint8_t *bp)
{
    bzero(aec->contexts, 19*sizeof(AecContext));
    aec->contexts[AEC_CX_UNI].state = 46;
    aec->contexts[AEC_CX_RL].state = 3;
    aec->contexts[0].state = 4;
    aec->curctx = aec->contexts;

    aec->a = 0x8000;
    aec->c = 0;
    aec->bp = bp-1;
    aec->bpstart = bp;
    if (*aec->bp == 0xff)
        aec->ct = 13;
    else
        aec->ct = 12;
}

void aec_encode(AecState *aec, int cx, int d)
{
    aec->curctx = aec->contexts + cx;
    if (aec->curctx->mps == d){
        codemps(aec);
    }
    else{
        codelps(aec);
    }
}

int aec_flush(AecState *aec)
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
