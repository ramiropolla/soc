/*
 * JPEG2000 encoder and decoder
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

int ff_j2k_ceildivpow2(int a, int b)
{
    return (a + (1 << b) - 1)>> b;
}

int ff_j2k_ceildiv(int a, int b)
{
    return (a + b - 1) / b;
}

/* tag tree routines */

/** allocate the memory for tag tree */
//TODO: optimize (too many mallocs)
static J2kTgtNode *tag_tree_alloc(int w, int h)
{
    int i;
   J2kTgtNode *t = av_malloc(w*h*sizeof(J2kTgtNode));
    if (t == NULL)
        return NULL;
    for (i = 0; i < w*h; i++){
        t[i].val = 0;
        t[i].vis = 0;
    }
    return t;
}

J2kTgtNode *ff_j2k_tag_tree_init(int w, int h)
{
    J2kTgtNode *res = tag_tree_alloc(w, h),
               *t = res;

    if (res == NULL)
        return NULL;

    while (w > 1 || h > 1){
        int pw = w, ph = h;
        int i, j;
        J2kTgtNode *t2;

        w = (w+1) >> 1;
        h = (h+1) >> 1;
        t2 = tag_tree_alloc(w, h);
        if (t2 == NULL)
            return NULL;

        for (i = 0; i < ph; i++)
            for (j = 0; j < pw; j++){
                t[i*pw + j].parent = &t2[(i>>1)*w + (j>>1)];
            }
        t = t2;
    }
    t[0].parent = NULL;
    return res;
}

void ff_j2k_tag_tree_destroy(J2kTgtNode *tree)
{
    while (tree != NULL){
        J2kTgtNode *parent = tree[0].parent;
        av_free(tree);
        tree = parent;
    }
}

int ff_j2k_getnbctxno(int flag, int bandno)
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

int ff_j2k_getrefctxno(int flag)
{
    if (!(flag & J2K_T1_REF)){
        if (flag & J2K_T1_SIG_NB)
            return 15;
        return 14;
    }
    return 16;
}

int ff_j2k_getsgnctxno(int flag, int *xorbit)
{
    int vcontrib, hcontrib;
    const int contribtab[3][3] = {{0, -1, 1}, {-1, -1, 0}, {1, 0, 1}};
    const int ctxlbltab[3][3] = {{13, 12, 11}, {10, 9, 10}, {11, 12, 13}};
    const int xorbittab[3][3] = {{1, 1, 1,}, {1, 0, 0}, {0, 0, 0}};

    hcontrib = contribtab[flag & J2K_T1_SIG_E ? flag & J2K_T1_SGN_E ? 1:2:0]
                         [flag & J2K_T1_SIG_W ? flag & J2K_T1_SGN_W ? 1:2:0]+1;
    vcontrib = contribtab[flag & J2K_T1_SIG_S ? flag & J2K_T1_SGN_S ? 1:2:0]
                         [flag & J2K_T1_SIG_N ? flag & J2K_T1_SGN_N ? 1:2:0]+1;
    *xorbit = xorbittab[hcontrib][vcontrib];
    return ctxlbltab[hcontrib][vcontrib];
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
