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
 *
 */

/**
 * JPEG2000 tables
 * @file j2k.h
 * @author Kamil Nowosad
 */

#ifndef _J2K_H_
#define _J2K_H_

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
    J2K_CRG,
    J2K_COM,
    J2K_SOT = 0xff90,
    J2K_SOP,
    J2K_EPH,
    J2K_SOD,
    J2K_EOC = 0xffd9,
};


/* arithmetic entropy coder context */
//TODO: optimize [nice solution in openjpeg]
typedef struct {
        unsigned int qe;
        unsigned int nmps;
        unsigned int nlps;
        unsigned int sw;
} J2kAecState;

const static J2kAecState aec_cx_states[47] = {
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

typedef struct {
    unsigned int state;
    unsigned int mps;
} J2kAecContext;

typedef struct {
    uint8_t *bp, *bpstart;
    unsigned int a;
    unsigned int c;
    unsigned int ct;
    J2kAecContext contexts[19];
    J2kAecContext *curctx;
} J2kAec;

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

#define J2K_T1_CTX_RL  17
#define J2K_T1_CTX_UNI 18

typedef struct {
    int data[J2K_MAX_CBLKW][J2K_MAX_CBLKH];
    int flags[J2K_MAX_CBLKW+2][J2K_MAX_CBLKH+2];
    J2kAec aec;
} J2kT1Context;

#endif
