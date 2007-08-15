/*
 * JPEG2000 tables
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
 * JPEG2000 tables
 * @file j2k.h
 * @author Kamil Nowosad
 */

#ifndef J2K_H
#define J2K_H

#include "aec.h"

enum J2kMarkers{
    J2K_SOC = 0xff4f,
    J2K_SIZ = 0xff51,
    J2K_COD,
    J2K_COC,
    J2K_TLM = 0xff55,
    J2K_PLM = 0xff57,
    J2K_PLT,
    J2K_QCD = 0xff5c,
    J2K_QCC,
    J2K_RGN,
    J2K_POC,
    J2K_PPM,
    J2K_PPT,
    J2K_CRG = 0xff63,
    J2K_COM,
    J2K_SOT = 0xff90,
    J2K_SOP,
    J2K_EPH,
    J2K_SOD,
    J2K_EOC = 0xffd9,
};

enum J2kTransform{
    J2K_DWT97,
    J2K_DWT53
};

#define J2K_MAX_CBLKW 64
#define J2K_MAX_CBLKH 64

// T1 flags
// flags determining significance of neighbour coefficients
#define J2K_T1_SIG_N  0x0001
#define J2K_T1_SIG_E  0x0002
#define J2K_T1_SIG_W  0x0004
#define J2K_T1_SIG_S  0x0008
#define J2K_T1_SIG_NE 0x0010
#define J2K_T1_SIG_NW 0x0020
#define J2K_T1_SIG_SE 0x0040
#define J2K_T1_SIG_SW 0x0080
#define J2K_T1_SIG_NB (J2K_T1_SIG_N | J2K_T1_SIG_E | J2K_T1_SIG_S | J2K_T1_SIG_W \
                      |J2K_T1_SIG_NE | J2K_T1_SIG_NW | J2K_T1_SIG_SE | J2K_T1_SIG_SW)
// flags determining sign bit of neighbour coefficients
#define J2K_T1_SGN_N  0x0100
#define J2K_T1_SGN_S  0x0200
#define J2K_T1_SGN_W  0x0400
#define J2K_T1_SGN_E  0x0800

#define J2K_T1_VIS    0x1000
#define J2K_T1_SIG    0x2000
#define J2K_T1_REF    0x4000

typedef struct {
    int data[J2K_MAX_CBLKW][J2K_MAX_CBLKH];
    int flags[J2K_MAX_CBLKW+2][J2K_MAX_CBLKH+2];
    AecState aec;
} J2kT1Context;

typedef struct J2kTgtNode {
    uint8_t val;
    uint8_t vis;
    struct J2kTgtNode *parent;
} J2kTgtNode;

/** debug routines */
#if 0
#undef fprintf
#undef printf
void ff_j2k_printv(int *tab, int l);
void ff_j2k_printu(uint8_t *tab, int l);
#endif

/** misc tools */
static inline int ff_j2k_ceildivpow2(int a, int b)
{
    return (a + (1 << b) - 1)>> b;
}

static inline int ff_j2k_ceildiv(int a, int b)
{
    return (a + b - 1) / b;
}

/** tag tree routines */
J2kTgtNode *ff_j2k_tag_tree_init(int w, int h);

/** TIER-1 routines */
int ff_j2k_getnbctxno(int flag, int bandno);
int ff_j2k_getrefctxno(int flag);
int ff_j2k_getsgnctxno(int flag, int *xorbit);

void ff_j2k_set_significant(J2kT1Context *t1, int x, int y);

#endif
