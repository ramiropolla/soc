/*
 * JPEG2000 encoder and decoder common functions
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
 * JPEG2000 image encoder and decoder common functions
 * @file j2k.c
 * @author Kamil Nowosad
 */

#include "j2k.h"
#include "avcodec.h"

#if 0
void ff_j2k_printv(int *tab, int l)
{
    int i;
    for (i = 0; i < l; i++)
        printf("%.3d ", tab[i]);
    printf("\n");
}

void ff_j2k_printu(uint8_t *tab, int l)
{
    int i;
    for (i = 0; i < l; i++)
        printf("%.3hd ", tab[i]);
    printf("\n");
}
#endif

/* tag tree routines */

/** allocate the memory for tag tree */
J2kTgtNode *ff_j2k_tag_tree_init(int w, int h)
{
    int size = 1, pw = w, ph = h;
    J2kTgtNode *res, *t, *t2;

    while (pw > 1 || ph > 1){
        size += pw*ph;
        pw = (pw+1) >> 1;
        ph = (ph+1) >> 1;
    }
    t = res = av_mallocz(size*sizeof(J2kTgtNode));

    if (res == NULL)
        return NULL;

    while (w > 1 || h > 1){
        int i, j;
        pw = w;
        ph = h;

        w = (w+1) >> 1;
        h = (h+1) >> 1;
        t2 = t + pw*ph;

        for (i = 0; i < ph; i++)
            for (j = 0; j < pw; j++){
                t[i*pw + j].parent = &t2[(i>>1)*w + (j>>1)];
            }
        t = t2;
    }
    t[0].parent = NULL;
    return res;
}

uint8_t ff_j2k_nbctxno_lut[256][4];

static int getnbctxno(int flag, int bandno)
{
    int h, v, d;

    h = ((flag & J2K_T1_SIG_E) ? 1:0)+
        ((flag & J2K_T1_SIG_W) ? 1:0);
    v = ((flag & J2K_T1_SIG_N) ? 1:0)+
        ((flag & J2K_T1_SIG_S) ? 1:0);
    d = ((flag & J2K_T1_SIG_NE) ? 1:0)+
        ((flag & J2K_T1_SIG_NW) ? 1:0)+
        ((flag & J2K_T1_SIG_SE) ? 1:0)+
        ((flag & J2K_T1_SIG_SW) ? 1:0);
    switch(bandno){
        case 0: // LL || LH
        case 2:
            if (h == 2) return 8;
            if (h == 1){
                if (v >= 1) return 7;
                if (d >= 1) return 6;
                return 5;
            }
            if (v == 2) return 4;
            if (v == 1) return 3;
            if (d >= 2) return 2;
            if (d == 1) return 1;
            return 0;
        case 1: // HL
            if (v == 2) return 8;
            if (v == 1){
                if (h >= 1) return 7;
                if (d >= 1) return 6;
                return 5;
            }
            if (h == 2) return 4;
            if (h == 1) return 3;
            if (d >= 2) return 2;
            if (d >= 1) return 1;
            return 0;
        case 3:
            if (d >= 3) return 8;
            if (d == 2){
                if (h+v >= 1) return 7;
                return 6;
            }
            if (d == 1){
                if (h+v >= 2) return 5;
                if (h+v == 1) return 4;
                return 3;
            }
            if (h+v >= 2) return 2;
            if (h+v == 1) return 1;
            return 0;
    }
    assert(0);
}

uint8_t ff_j2k_sgnctxno_lut[16][16], ff_j2k_xorbit_lut[16][16];

static int getsgnctxno(int flag, uint8_t *xorbit)
{
    int vcontrib, hcontrib;
    static const int contribtab[3][3] = {{0, -1, 1}, {-1, -1, 0}, {1, 0, 1}};
    static const int ctxlbltab[3][3] = {{13, 12, 11}, {10, 9, 10}, {11, 12, 13}};
    static const int xorbittab[3][3] = {{1, 1, 1,}, {1, 0, 0}, {0, 0, 0}};

    hcontrib = contribtab[flag & J2K_T1_SIG_E ? flag & J2K_T1_SGN_E ? 1:2:0]
                         [flag & J2K_T1_SIG_W ? flag & J2K_T1_SGN_W ? 1:2:0]+1;
    vcontrib = contribtab[flag & J2K_T1_SIG_S ? flag & J2K_T1_SGN_S ? 1:2:0]
                         [flag & J2K_T1_SIG_N ? flag & J2K_T1_SGN_N ? 1:2:0]+1;
    *xorbit = xorbittab[hcontrib][vcontrib];
    return ctxlbltab[hcontrib][vcontrib];
}

void ff_j2k_init_tier1_luts()
{
    int i, j;
    for (i = 0; i < 256; i++)
        for (j = 0; j < 4; j++)
            ff_j2k_nbctxno_lut[i][j] = getnbctxno(i, j);
    for (i = 0; i < 16; i++)
        for (j = 0; j < 16; j++)
            ff_j2k_sgnctxno_lut[i][j] = getsgnctxno(i + (j << 8), &ff_j2k_xorbit_lut[i][j]);
}

void ff_j2k_set_significant(J2kT1Context *t1, int x, int y)
{
    x++; y++;
    t1->flags[y][x] |= J2K_T1_SIG;
    if (t1->data[y-1][x-1] < 0){
        t1->flags[y][x+1] |= J2K_T1_SIG_W | J2K_T1_SGN_W;
        t1->flags[y][x-1] |= J2K_T1_SIG_E | J2K_T1_SGN_E;
        t1->flags[y+1][x] |= J2K_T1_SIG_N | J2K_T1_SGN_N;
        t1->flags[y-1][x] |= J2K_T1_SIG_S | J2K_T1_SGN_S;
    }
    else{
        t1->flags[y][x+1] |= J2K_T1_SIG_W;
        t1->flags[y][x-1] |= J2K_T1_SIG_E;
        t1->flags[y+1][x] |= J2K_T1_SIG_N;
        t1->flags[y-1][x] |= J2K_T1_SIG_S;
    }
    t1->flags[y+1][x+1] |= J2K_T1_SIG_NW;
    t1->flags[y+1][x-1] |= J2K_T1_SIG_NE;
    t1->flags[y-1][x+1] |= J2K_T1_SIG_SW;
    t1->flags[y-1][x-1] |= J2K_T1_SIG_SE;
}
