/*
 * Arithmetic entropy decoder
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
 * Arithmetic entropy decoder
 * @file aecdec.c
 * @author Kamil Nowosad
 */

#include "aec.h"

static void bytein(AecState *aec)
{
    if (*aec->bp == 0xff){
        if (*(aec->bp+1) > 0x8f)
            aec->ct = 8;
        else{
            aec->bp++;
            aec->c += 0xfe00 - (*aec->bp << 9);
            aec->ct = 7;
        }
    }
    else{
        aec->bp++;
        aec->c += 0xff00 - (*aec->bp << 8);
        aec->ct = 8;
    }
}

static int exchange(AecState *aec, int lps)
{
    int d;
    if ((aec->a < cx_states[aec->curctx->state].qe) ^ (!lps)){
        if (lps)
            aec->a = cx_states[aec->curctx->state].qe;
        d = aec->curctx->mps;
        aec->curctx->state = cx_states[aec->curctx->state].nmps;
    }
    else{
        if (lps)
            aec->a = cx_states[aec->curctx->state].qe;
        d = 1 - aec->curctx->mps;
        if (cx_states[aec->curctx->state].sw)
            aec->curctx->mps ^= 1;
        aec->curctx->state = cx_states[aec->curctx->state].nlps;
    }
    // renormd:
    do{
        if (!aec->ct)
            bytein(aec);
        aec->a += aec->a;
        aec->c += aec->c;
        aec->ct--;
    } while (!(aec->a & 0x8000));
    return d;
}

void ff_aec_initdec(AecState *aec, uint8_t *bp)
{
    ff_aec_init_contexts(aec);
    aec->bp = bp;
    aec->c = (*aec->bp ^ 0xff) << 16;
    bytein(aec);
    aec->c = aec->c << 7;
    aec->ct -= 7;
    aec->a = 0x8000;
}

int ff_aec_decode(AecState *aec, int cx)
{
    aec->curctx = aec->contexts + cx;
    aec->a -= cx_states[aec->curctx->state].qe;
    if ((aec->c >> 16) < aec->a){
        if (aec->a & 0x8000)
            return aec->curctx->mps;
        else
            return exchange(aec, 0);
    } else {
        aec->c -= aec->a << 16;
        return exchange(aec, 1);
    }
}
