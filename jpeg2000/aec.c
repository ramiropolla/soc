/*
 * Arithmetic entropy encoder and decoder common functions
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
 * Arithmetic entropy coder and decoder common functions
 * @file aec.c
 * @author Kamil Nowosad
 */

#include "aec.h"

const AecCxState cx_states[47] = {
    {0x5601,  1,  1, 1},
    {0x3401,  2,  6, 0},
    {0x1801,  3,  9, 0},
    {0x0AC1,  4, 12, 0},
    {0x0521,  5, 29, 0},
    {0x0221, 38, 33, 0},
    {0x5601,  7,  6, 1},
    {0x5401,  8, 14, 0},
    {0x4801,  9, 14, 0},
    {0x3801, 10, 14, 0},
    {0x3001, 11, 17, 0},
    {0x2401, 12, 18, 0},
    {0x1C01, 13, 20, 0},
    {0x1601, 29, 21, 0},
    {0x5601, 15, 14, 1},
    {0x5401, 16, 14, 0},
    {0x5101, 17, 15, 0},
    {0x4801, 18, 16, 0},
    {0x3801, 19, 17, 0},
    {0x3401, 20, 18, 0},
    {0x3001, 21, 19, 0},
    {0x2801, 22, 19, 0},
    {0x2401, 23, 20, 0},
    {0x2201, 24, 21, 0},
    {0x1C01, 25, 22, 0},
    {0x1801, 26, 23, 0},
    {0x1601, 27, 24, 0},
    {0x1401, 28, 25, 0},
    {0x1201, 29, 26, 0},
    {0x1101, 30, 27, 0},
    {0x0AC1, 31, 28, 0},
    {0x09C1, 32, 29, 0},
    {0x08A1, 33, 30, 0},
    {0x0521, 34, 31, 0},
    {0x0441, 35, 32, 0},
    {0x02A1, 36, 33, 0},
    {0x0221, 37, 34, 0},
    {0x0141, 38, 35, 0},
    {0x0111, 39, 36, 0},
    {0x0085, 40, 37, 0},
    {0x0049, 41, 38, 0},
    {0x0025, 42, 39, 0},
    {0x0015, 43, 40, 0},
    {0x0009, 44, 41, 0},
    {0x0005, 45, 42, 0},
    {0x0001, 45, 43, 0},
    {0x5601, 46, 46, 0}
};

void ff_aec_init_contexts(AecState *aec)
{
    memset(aec->contexts, 0, sizeof(aec->contexts));
    aec->contexts[AEC_CX_UNI].state = 46;
    aec->contexts[AEC_CX_RL].state = 3;
    aec->contexts[0].state = 4;
    aec->curctx = aec->contexts;
}
