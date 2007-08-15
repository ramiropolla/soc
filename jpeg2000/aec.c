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

void ff_aec_init_contexts(AecState *aec)
{
    memset(aec->contexts, 0, 19*sizeof(AecContext));
    aec->contexts[AEC_CX_UNI].state = 46;
    aec->contexts[AEC_CX_RL].state = 3;
    aec->contexts[0].state = 4;
    aec->curctx = aec->contexts;
}
