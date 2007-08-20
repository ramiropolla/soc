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
            aec->c++;
        else{
            aec->bp++;
            aec->c += 2 + 0xfe00 - (*aec->bp << 9);
        }
    }
    else{
        aec->bp++;
        aec->c += 1 + 0xff00 - (*aec->bp << 8);
    }
}

static int exchange(AecState *aec, uint8_t *cxstate, int lps)
{
    int d;
    if ((aec->a < ff_aec_qe[*cxstate]) ^ (!lps)){
        if (lps)
            aec->a = ff_aec_qe[*cxstate];
        d = *cxstate & 1;
        *cxstate = ff_aec_nmps[*cxstate];
    }
    else{
        if (lps)
            aec->a = ff_aec_qe[*cxstate];
        d = 1 - (*cxstate & 1);
        *cxstate = ff_aec_nlps[*cxstate];
    }
    // renormd:
    do{
        if (!(aec->c & 0xff)){
            aec->c -= 0x100;
            bytein(aec);
        }
        aec->a += aec->a;
        aec->c += aec->c;
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
    aec->a = 0x8000;
}

int ff_aec_decode(AecState *aec, uint8_t *cxstate)
{
    aec->a -= ff_aec_qe[*cxstate];
    if ((aec->c >> 16) < aec->a){
        if (aec->a & 0x8000)
            return *cxstate & 1;
        else
            return exchange(aec, cxstate, 0);
    } else {
        aec->c -= aec->a << 16;
        return exchange(aec, cxstate, 1);
    }
}
