/*
 * JPEG2000 image encoder
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
 * JPEG2000 image encoder
 * @file j2kenc.c
 * @author Kamil Nowosad
 */

#include "avcodec.h"
#include "bytestream.h"
#include "j2k.h"
#include "common.h"

// TODO: doxygen-compatible comments

typedef struct {
    int length;
    int npassess;
    int zerobits;
    int zero;
    uint8_t data[8192];
} J2kCblk; // code block

typedef struct J2kTgtNode {
    uint8_t val;
    uint8_t vis;
    struct J2kTgtNode *parent;
} J2kTgtNode;

typedef struct {
    int xi0, xi1, yi0, yi1; /// indices of codeblocks ([xi0, xi1))
} J2kPrec; // precinct

typedef struct {
    int x0, x1, y0, y1;
    int cblkw, cblkh;
    int cblknx, cblkny;
    J2kPrec *prec;
    J2kCblk *cblk;
} J2kBand; // subband

typedef struct {
    int x0, x1, y0, y1;
    int nbands;
    int nprecw, nprech;
    J2kBand *band;
} J2kResLevel; // resolution level

typedef struct {
   J2kResLevel *reslevel;
   int *data;
   int x0, x1, y0, y1;
} J2kComponent;

typedef struct { // flatten with context
   J2kComponent *comp;
} J2kTile;

typedef struct {
    AVCodecContext *avctx;
    AVFrame *picture;

    int Xsiz, Ysiz; // image width and height
    unsigned int bpp;
    int ncomponents;
    int ppx, ppy; // exponent of the precinct size [global]
    int xcb, ycb; // exponent of the code block size
    int XTsiz, YTsiz; // tile size
    int numXtiles, numYtiles;

    int expn;
    int nguardbits;

//    int *samples[3];
    int nreslevels; // number of resolution levels
    uint8_t *buf_start;
    uint8_t *buf;
    uint8_t *buf_end;
    int bit_index;

    J2kTile *tile;
} J2kEncoderContext;


/* debug */
#if 0
#undef ifprintf
#undef printf

static void nspaces(FILE *fd, int n)
{
    while(n--) putc(' ', fd);
}

static void printv(int *tab, int l)
{
    int i;
    for (i = 0; i < l; i++)
        printf("%.3d ", tab[i]);
    printf("\n");
}

static void printu(uint8_t *tab, int l)
{
    int i;
    for (i = 0; i < l; i++)
        printf("%.3hd ", tab[i]);
    printf("\n");
}

static void printcomp(J2kComponent *comp)
{
    int i;
    for (i = 0; i < comp->y1 - comp->y0; i++)
        printv(comp->data + i * (comp->x1 - comp->x0), comp->x1 - comp->x0);
}

static void dump(J2kEncoderContext *s, FILE *fd)
{
    int tileno, compno, reslevelno, bandno, precno;
    fprintf(fd, "XSiz = %d, YSiz = %d, XTsiz = %d, YTsiz = %d\n"
                "numXtiles = %d, numYtiles = %d, ncomponents = %d\n"
                "tiles:\n",
            s->Xsiz, s->Ysiz, s->XTsiz, s->YTsiz,
            s->numXtiles, s->numYtiles, s->ncomponents);
    for (tileno = 0; tileno < s->numXtiles * s->numYtiles; tileno++){
        J2kTile *tile = s->tile + tileno;
        nspaces(fd, 2);
        fprintf(fd, "tile %d:\n", tileno);
        for(compno = 0; compno < s->ncomponents; compno++){
            J2kComponent *comp = tile->comp + compno;
            nspaces(fd, 4);
            fprintf(fd, "component %d:\n", compno);
            nspaces(fd, 4);
            fprintf(fd, "x0 = %d, x1 = %d, y0 = %d, y1 = %d\n",
                        comp->x0, comp->x1, comp->y0, comp->y1);
            for(reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
                J2kResLevel *reslevel = comp->reslevel + reslevelno;
                nspaces(fd, 6);
                fprintf(fd, "reslevel %d:\n", reslevelno);
                nspaces(fd, 6);
                fprintf(fd, "x0 = %d, x1 = %d, y0 = %d, y1 = %d, nbands = %d\n",
                        reslevel->x0, reslevel->x1, reslevel->y0,
                        reslevel->y1, reslevel->nbands);
                for(bandno = 0; bandno < reslevel->nbands; bandno++){
                    J2kBand *band = reslevel->band + bandno;
                    nspaces(fd, 8);
                    fprintf(fd, "band %d:\n", bandno);
                    nspaces(fd, 8);
                    fprintf(fd, "x0 = %d, x1 = %d, y0 = %d, y1 = %d,"
                                "cblkw = %d, cblkh = %d cblknx = %d cblkny = %d\n",
                                band->x0, band->x1,
                                band->y0, band->y1,
                                band->cblkw, band->cblkh,
                                band->cblknx, band->cblkny);
                    for (precno = 0; precno < reslevel->nprecw * reslevel->nprech; precno++){
                        J2kPrec *prec = band->prec + precno;
                        nspaces(fd, 10);
                        fprintf(fd, "prec %d:\n", precno);
                        nspaces(fd, 10);
                        fprintf(fd, "xi0 = %d, xi1 = %d, yi0 = %d, yi1 = %d\n",
                                     prec->xi0, prec->xi1, prec->yi0, prec->yi1);
                    }
                }
            }
        }
    }
}
#endif

/* misc tools */
static int ceildivpow2(int a, int b)
{
    return (a + (1 << b) - 1)>> b;
}

static int ceildiv(int a, int b)
{
    return (a + b - 1) / b;
}

/* bitstream routines */

/* put n times val bit */
static void put_bits(J2kEncoderContext *s, int val, int n) // TODO: optimize
{
    while (n-- > 0){
        if (s->bit_index == 8)
        {
            s->bit_index = *s->buf == 0xff ? 1:0;
            *(++s->buf) = 0;
        }
        *s->buf |= val << (7 - s->bit_index++);
    }
}

/* put n least significant bits of a number num */
static void put_num(J2kEncoderContext *s, int num, int n)
{
    while(--n >= 0)
        put_bits(s, (num & (1<<n)) ? 1:0, 1);
}

/* flush the bitstream */
static void j2k_flush(J2kEncoderContext *s)
{
    if (s->bit_index){
        s->bit_index = 0;
        s->buf++;
    }
}

/* tag tree routines */

/* allocate the memory for tag tree */
//TODO: optimize (too many mallocs)
static J2kTgtNode *tag_tree_alloc(int w, int h)
{
    int i;
   J2kTgtNode *t = av_malloc(w*h*sizeof(J2kTgtNode));
    if (t == NULL)
        return NULL;
    for (i = 0; i < w*h; i++){
        t[i].val = 0xff;
        t[i].vis = 0;
    }
    return t;
}

static J2kTgtNode *tag_tree_init(int w, int h)
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

/* code the value stored in node */
static void tag_tree_code(J2kEncoderContext *s, J2kTgtNode *node)
{
    J2kTgtNode *stack[30];
    int sp = 1, curval = 0;
    stack[0] = node;

    node = node->parent;
    while(node != NULL){
        if (node->vis){
            curval = node->val;
            break;
        }
        node->vis++;
        stack[sp++] = node;
        node = node->parent;
    }
    while(--sp >= 0){
        put_bits(s, 0, stack[sp]->val - curval);
        put_bits(s, 1, 1);
        curval = stack[sp]->val;
    }
}

/* update the value in node */
static void tag_tree_update(J2kTgtNode *node)
{
    int lev = 0;
    while (node->parent != NULL){
        if (node->parent->val <= node->val)
            break;
        node->parent->val = node->val;
        node = node->parent;
        lev++;
    }
}

static void tag_tree_destroy(J2kTgtNode *tree)
{
    while (tree != NULL){
        J2kTgtNode *parent = tree[0].parent;
        av_free(tree);
        tree = parent;
    }
}

/* marker segments */
static void put_marker(J2kEncoderContext *s, uint16_t marker)
{
    bytestream_put_be16(&s->buf, marker);
}

static void put_soc(J2kEncoderContext *s)
{
    put_marker(s, J2K_SOC);
}

static void put_siz(J2kEncoderContext *s)
{
    int i;

    put_marker(s, J2K_SIZ);
    bytestream_put_be16(&s->buf, 38 + 3 * s->ncomponents); // Lsiz
    bytestream_put_be16(&s->buf, 0); // Rsiz
    bytestream_put_be32(&s->buf, s->Xsiz); // Xsiz
    bytestream_put_be32(&s->buf, s->Ysiz); // Ysiz
    bytestream_put_be32(&s->buf, 0); // X0Siz
    bytestream_put_be32(&s->buf, 0); // Y0Siz

    bytestream_put_be32(&s->buf, s->XTsiz); // XTSiz
    bytestream_put_be32(&s->buf, s->YTsiz); // YTSiz
    bytestream_put_be32(&s->buf, 0); // XT0Siz
    bytestream_put_be32(&s->buf, 0); // YT0Siz
    bytestream_put_be16(&s->buf, s->ncomponents); // CSiz

    for (i = 0; i < 3; i++){ // Ssiz_i XRsiz_i, YRsiz_i
        bytestream_put_byte(&s->buf, 7);
        bytestream_put_byte(&s->buf, 1);
        bytestream_put_byte(&s->buf, 1);
    }
}

static void put_cod(J2kEncoderContext *s)
{
    put_marker(s, J2K_COD);
    bytestream_put_be16(&s->buf, 12); // Lcod
    bytestream_put_byte(&s->buf, 0);  // Scod
    // SGcod
    bytestream_put_byte(&s->buf, 0); // progression level
    bytestream_put_be16(&s->buf, 1); // num of layers
    bytestream_put_byte(&s->buf, 0); // multiple component transformation
    // SPcod
    bytestream_put_byte(&s->buf, s->nreslevels - 1); // num of decomp. levels
    bytestream_put_byte(&s->buf, s->xcb-2); // cblk width
    bytestream_put_byte(&s->buf, s->ycb-2); // cblk height
    bytestream_put_byte(&s->buf, 0); // cblk style
    bytestream_put_byte(&s->buf, 1); // transformation
}

static void put_qcd(J2kEncoderContext *s)
{
    int reslevelno;
    put_marker(s, J2K_QCD);
    bytestream_put_be16(&s->buf, 4+3*(s->nreslevels-1));  // LQcd
    bytestream_put_byte(&s->buf, s->nguardbits << 5);  // Sqcd
    for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
        int bandno, nbands = reslevelno == 0 ? 1:3;
        for (bandno = 0; bandno < nbands; bandno++)
            bytestream_put_byte(&s->buf, s->expn << 3);
    }
}

static uint8_t *put_sot(J2kEncoderContext *s, int tileno)
{
    uint8_t *psotptr;
    put_marker(s, J2K_SOT);
    bytestream_put_be16(&s->buf, 10); // Lsot
    bytestream_put_be16(&s->buf, tileno); // Isot

    psotptr = s->buf;
    bytestream_put_be32(&s->buf, 0); // Psot (filled in later)

    bytestream_put_byte(&s->buf, 0); // TPsot
    bytestream_put_byte(&s->buf, 1); // TNsot
    return psotptr;
}

/* compute the sizes of tiles, resolution levels, bands, etc.
 * allocate memory for them
 * divide the input image into tile-components
 */
static int init_tiles(J2kEncoderContext *s)
{
    // only one tile
    // only rgb24 supported now
    int y, x, tno, i;

    s->numXtiles = ceildiv(s->Xsiz, s->XTsiz);
    s->numYtiles = ceildiv(s->Ysiz, s->YTsiz);

    s->tile = av_malloc(s->numXtiles * s->numYtiles * sizeof(J2kTile));
    if (s->tile == NULL)
        return -1;
    for (tno = 0; tno < s->numXtiles * s->numYtiles; tno++){
        J2kTile *tile = s->tile + tno;
        int p = tno % s->numXtiles;
        int q = tno / s->numXtiles;
        int compno;

        tile->comp = av_malloc(s->ncomponents * sizeof(J2kComponent));
        if (tile->comp == NULL)
            return -1;
        for (compno = 0; compno < s->ncomponents; compno++){
            J2kComponent *comp = tile->comp + compno;
            int reslevelno;

            comp->x0 = p * s->XTsiz;
            comp->x1 = FFMIN((p+1)*s->XTsiz, s->Xsiz);
            comp->y0 = q * s->YTsiz;
            comp->y1 = FFMIN((q+1)*s->YTsiz, s->Ysiz);

            comp->data = av_malloc((comp->y1 - comp->y0) * (comp->x1 -comp->x0) * sizeof(int));
            if (comp->data == NULL)
                return -1;
            comp->reslevel = av_malloc(s->nreslevels * sizeof(J2kResLevel));
            if (comp->reslevel == NULL)
                return -1;
            for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
                int bandno;
                int n = s->nreslevels - reslevelno;
                J2kResLevel *reslevel = comp->reslevel + reslevelno;

                reslevel->x0 = ceildivpow2(comp->x0, s->nreslevels - reslevelno - 1);
                reslevel->x1 = ceildivpow2(comp->x1, s->nreslevels - reslevelno - 1);
                reslevel->y0 = ceildivpow2(comp->y0, s->nreslevels - reslevelno - 1);
                reslevel->y1 = ceildivpow2(comp->y1, s->nreslevels - reslevelno - 1);

                if (reslevelno == 0)
                    reslevel->nbands = 1;
                else
                    reslevel->nbands = 3;

                if (reslevel->x1 == reslevel->x0)
                    reslevel->nprecw = 0;
                else
                    reslevel->nprecw = ceildivpow2(reslevel->x1, s->ppx) - reslevel->x0 / (1<<s->ppx);

                if (reslevel->y1 == reslevel->y0)
                    reslevel->nprech = 0;
                else
                    reslevel->nprech = ceildivpow2(reslevel->y1, s->ppy) - reslevel->y0 / (1<<s->ppy);

                reslevel->band = av_malloc(reslevel->nbands * sizeof(J2kBand));
                if (reslevel->band == NULL)
                    return -1;
                for (bandno = 0; bandno < reslevel->nbands; bandno++){
                    J2kBand *band = reslevel->band + bandno;
                    int cblkno, precx, precy, precno;
                    int x0, y0, x1, y1;
                    int xi0, yi0, xi1, yi1;
                    int cblkperprecw, cblkperprech;

                    if (reslevelno == 0){  // the same everywhere
                        band->cblkw = 1 << FFMIN(s->xcb, s->ppx-1);
                        band->cblkh = 1 << FFMIN(s->ycb, s->ppy-1);

                        band->x0 = ceildivpow2(comp->x0, n-1);
                        band->x1 = ceildivpow2(comp->x1, n-1);
                        band->y0 = ceildivpow2(comp->y0, n-1);
                        band->y1 = ceildivpow2(comp->y1, n-1);
                    }
                    else{
                        band->cblkw = 1 << FFMIN(s->xcb, s->ppx);
                        band->cblkh = 1 << FFMIN(s->ycb, s->ppy);

                        band->x0 = ceildivpow2(comp->x0 - (1 << (n-1)) * ((bandno+1)&1), n);
                        band->x1 = ceildivpow2(comp->x1 - (1 << (n-1)) * ((bandno+1)&1), n);
                        band->y0 = ceildivpow2(comp->y0 - (1 << (n-1)) * (((bandno+1)&2)>>1), n);
                        band->y1 = ceildivpow2(comp->y1 - (1 << (n-1)) * (((bandno+1)&2)>>1), n);
                    }

                    band->cblknx = ceildiv(band->x1 - band->x0, band->cblkw);
                    band->cblkny = ceildiv(band->y1 - band->y0, band->cblkh);

                    band->cblk = av_malloc(band->cblknx * band->cblkny * sizeof(J2kCblk));
                    if (band->cblk == NULL)
                        return -1;
                    band->prec = av_malloc(reslevel->nprecw * reslevel->nprech * sizeof(J2kPrec));
                    if (band->prec == NULL)
                        return -1;

                    for (cblkno = 0; cblkno < band->cblknx * band->cblkny; cblkno++){
                        band->cblk[cblkno].zero = 0;
                    }

                    y0 = band->y0;
                    y1 = (band->y0 + (1<<s->ppy))/(1<<s->ppy)*(1<<s->ppy) - band->y0;
                    yi0 = 0;
                    yi1 = ceildiv(y1 - y0, 1<<s->ycb) * (1<<s->ycb);
                    yi1 = FFMIN(yi1, band->cblkny);
                    cblkperprech = 1<<(s->ppy - s->ycb);
                    for (precy = 0, precno = 0; precy < reslevel->nprech; precy++){
                        for (precx = 0; precx < reslevel->nprecw; precx++, precno++){
                            band->prec[precno].yi0 = yi0;
                            band->prec[precno].yi1 = yi1;
                        }
                        yi1 += cblkperprech;
                        yi0 = yi1 - cblkperprech;
                        yi1 = FFMIN(yi1, band->cblkny);
                    }
                    x0 = band->x0;
                    x1 = (band->x0 + (1<<s->ppx))/(1<<s->ppx)*(1<<s->ppx) - band->x0;
                    xi0 = 0;
                    xi1 = ceildiv(x1 - x0, 1<<s->xcb) * (1<<s->xcb);
                    xi1 = FFMIN(xi1, band->cblknx);
                    cblkperprecw = 1<<(s->ppx - s->xcb);
                    for (precx = 0, precno = 0; precx < reslevel->nprecw; precx++){
                        for (precy = 0; precy < reslevel->nprech; precy++, precno = 0){
                            band->prec[precno].xi0 = xi0;
                            band->prec[precno].xi1 = xi1;
                        }
                        xi1 += cblkperprecw;
                        xi0 = xi1 - cblkperprecw;
                        xi1 = FFMIN(xi1, band->cblknx);
                    }
                }
            }
        }
    }
    for (tno = 0; tno < s->numXtiles * s->numYtiles; tno++){
        J2kTile *tile = s->tile + tno;
        uint8_t *line = s->picture->data[0] + tile->comp[0].y0 * s->picture->linesize[0] + tile->comp[0].x0 * s->ncomponents;

        i = 0;
        for (y = tile->comp[0].y0; y < tile->comp[0].y1; y++){
            uint8_t *ptr = line;
            for (x = tile->comp[0].x0; x < tile->comp[0].x1; x++, i++){
                int compno;
                for (compno = 0; compno < s->ncomponents; compno++){
                    tile->comp[compno].data[i] = *ptr++  - (1 << 7);
                }
            }
            line += s->picture->linesize[0];
        }
    }
    return 0;
}

/* discrete wavelet transform routines */
static void sd_1d(int *p, int i0, int i1, int ileft, int iright)
{
#define PSE (i0 + FFMIN((i-i0+2*(i1-i0-1))%(2*(i1-i0-1)), 2*(i1-i0-1)-(i-i0+2*(i1-i0-1))%(2*(i1-i0-1))))
    int i;

    for (i = i0 - ileft; i < i0; i++){
        p[i] = p[PSE];
    }
    for (i = i1; i < i1+iright; i++){
        p[i] = p[PSE];
    }

    for (i = (i0+1)/2 - 1; i < (i1+1)/2; i++){
        p[2*i+1] -= (p[2*i] + p[2*i+2]) >> 1;
    }
    for (i = (i0+1)/2; i < (i1+1)/2; i++){
        p[2*i] += (p[2*i-1] + p[2*i+1] + 2) >> 2;
    }
#undef PSE
}

static void dwt_encode53(J2kEncoderContext *s, J2kComponent *comp)
{
    int lev = s->nreslevels,
        *t = comp->data, w = comp->x1 - comp->x0;
    int *ppv = av_malloc((comp->reslevel[lev-1].y1 + 4)*sizeof(int)), *pv = ppv+2;
    int *ppu = av_malloc((comp->reslevel[lev-1].x1 + 4)*sizeof(int)), *pu = ppu+2;

    while (--lev){
        int u0 = comp->reslevel[lev].x0,
            u1 = comp->reslevel[lev].x1,
            v0 = comp->reslevel[lev].y0,
            v1 = comp->reslevel[lev].y1,
            u = u0, v = v0;
        const static int tileft[2] = {2, 1}, tiright[2] = {1, 2};

        //VER_SD
        while (u < u1){
            int i, j;
            for (i = v0; i < v1; i++)
                pv[i] = t[w*(i-v0) + u-u0];
            sd_1d(pv, v0, v1, tileft[v0&1], tiright[v1&1]);

            // copy back and deinterleave
            for (i = (v0+1)/2, j = 0; i < (v1+1)/2; i++, j++){
                t[w*j + u-u0] = pv[2*i];
            }
            for (i = v0/2; i < v1/2; i++, j++){
                t[w*j + u-u0] = pv[2*i+1];
            }
            u++;
        }

        //HOR_SD
        while (v < v1){
            int i, j;
            for (i = u0; i < u1; i++)
                pu[i] = t[w*(v-v0) + i-u0];
            sd_1d(pu, u0, u1, tileft[u0&1], tiright[u1&1]);

            // copy back and deinterleave
            for (i = (u0+1)/2, j = 0; i < (u1+1)/2; i++, j++){
                t[w*(v-v0) + j] = pu[2*i];
            }
            for (i = u0/2; i < u1/2; i++, j++){
                t[w*(v-v0) + j] = pu[2*i + 1];
            }
            v++;
        }
    }
    av_free(ppv);
    av_free(ppu);
}

/* tier-1 routines */
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

static int getrefctxno(int flag)
{
    if (!(flag & J2K_T1_REF)){
        if (flag & J2K_T1_SIG_NB)
            return 15;
        return 14;
    }
    return 16;
}

static int getsgnctxno(int flag, int *xorbit)
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

static void set_significant(J2kT1Context *t1, int x, int y)
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


static void encode_sigpass(J2kT1Context *t1, int width, int height, int mask, int bandno)
{
    int i, j, k;
    for (i = 0; i < height; i += 4)
        for (j = 0; j < width; j++)
            for (k = i; k < height && k < i+4; k++){
                if (!(t1->flags[k+1][j+1] & J2K_T1_SIG) && (t1->flags[k+1][j+1] & J2K_T1_SIG_NB)){
                    int ctxno = getnbctxno(t1->flags[k+1][j+1], bandno),
                        bit = abs(t1->data[k][j]) & mask ? 1 : 0;
                    aec_encode(&t1->aec, ctxno, bit);
                    if (bit){
                        int xorbit;
                        int ctxno = getsgnctxno(t1->flags[k+1][j+1], &xorbit);
                        aec_encode(&t1->aec, ctxno, (t1->data[k][j] < 0 ? 1:0) ^ xorbit);
                        set_significant(t1, j, k);
                    }
                    t1->flags[k+1][j+1] |= J2K_T1_VIS;
                }
            }
}

static void encode_refpass(J2kT1Context *t1, int width, int height, int mask)
{
    int i, j, k;
    for (i = 0; i < height; i += 4)
        for (j = 0; j < width; j++)
            for (k = i; k < height && k < i+4; k++)
                if ((t1->flags[k+1][j+1] & (J2K_T1_SIG | J2K_T1_VIS)) == J2K_T1_SIG){
                    int ctxno = getrefctxno(t1->flags[k+1][j+1]);
                    aec_encode(&t1->aec, ctxno, abs(t1->data[k][j]) & mask ? 1:0);
                    t1->flags[k+1][j+1] |= J2K_T1_REF;
                }
}

static void encode_clnpass(J2kT1Context *t1, int width, int height, int mask, int bandno)
{
    int i, j, k;
    for (i = 0; i < height; i += 4)
        for (j = 0; j < width; j++){
            if (i + 3 < height && !(
            (t1->flags[i+1][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[i+2][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[i+3][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[i+4][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG))))
            {
                // aggregation mode
                int rlen;
                for (rlen = 0; rlen < 4; rlen++)
                    if (abs(t1->data[i+rlen][j]) & mask)
                        break;
                aec_encode(&t1->aec, AEC_CX_RL, rlen != 4);
                if (rlen == 4)
                    continue;
                aec_encode(&t1->aec, AEC_CX_UNI, rlen >> 1);
                aec_encode(&t1->aec, AEC_CX_UNI, rlen & 1);
                for (k = i + rlen; k < i + 4; k++){
                    if (!(t1->flags[k+1][j+1] & (J2K_T1_SIG | J2K_T1_VIS))){
                        int ctxno = getnbctxno(t1->flags[k+1][j+1], bandno);
                        if (k > i + rlen)
                            aec_encode(&t1->aec, ctxno, abs(t1->data[k][j]) & mask ? 1:0);
                        if (abs(t1->data[k][j]) & mask){ // newly significant
                            int xorbit;
                            int ctxno = getsgnctxno(t1->flags[k+1][j+1], &xorbit);
                            aec_encode(&t1->aec, ctxno, (t1->data[k][j] < 0 ? 1:0) ^ xorbit);
                            set_significant(t1, j, k);
                        }
                    }
                    t1->flags[k+1][j+1] &= ~J2K_T1_VIS;
                }
            }
            else{
                for (k = i; k < i + 4 && k < height; k++){
                    if (!(t1->flags[k+1][j+1] & (J2K_T1_SIG | J2K_T1_VIS))){
                        int ctxno = getnbctxno(t1->flags[k+1][j+1], bandno);
                        aec_encode(&t1->aec, ctxno, abs(t1->data[k][j]) & mask ? 1:0);
                        if (abs(t1->data[k][j]) & mask){ // newly significant
                            int xorbit;
                            int ctxno = getsgnctxno(t1->flags[k+1][j+1], &xorbit);
                            aec_encode(&t1->aec, ctxno, (t1->data[k][j] < 0 ? 1:0) ^ xorbit);
                            set_significant(t1, j, k);
                        }
                    }
                    t1->flags[k+1][j+1] &= ~J2K_T1_VIS;
                }
            }
        }
}

static void encode_cblk(J2kEncoderContext *s, J2kT1Context *t1, J2kCblk *cblk, int width, int height, int bandno)
{
    int pass_t = 2, passno, i, j, mask, max=0, nonzerobits;

    for (i = 0; i < height+2; i++)
        bzero(t1->flags[i], (width+2)*sizeof(int));

    for (i = 0; i < height; i++){
        for (j = 0; j < width; j++)
            max = FFMAX(max, abs(t1->data[i][j]));
    }

    if (max == 0){
        // XXX: both should be 0, but something goes wrong, when set so
        // - to be corrected
        nonzerobits = 1;
        mask = 1;
    }
    else{
        nonzerobits = av_log2(max) + 1;
        mask = 1 << (nonzerobits - 1);
    }

    aec_initenc(&t1->aec, cblk->data);

    cblk->zerobits = s->expn + s->nguardbits - 1 - nonzerobits;

    for (passno = 0; mask != 0; passno++){
        switch(pass_t){
            case 0: encode_sigpass(t1, width, height, mask, bandno);
                    break;
            case 1: encode_refpass(t1, width, height, mask);
                    break;
            case 2: encode_clnpass(t1, width, height, mask, bandno);
                    break;
        }
        if (++pass_t == 3){
            pass_t = 0;
            mask = mask >> 1;
        }
    }
    cblk->npassess = passno;

    // TODO: optional flush on each pass
    cblk->length = aec_flush(&t1->aec);
}

/* tier-2 routines: */

static void putnumpassess(J2kEncoderContext *s, int n)
{
    if (n == 1)
        put_num(s, 0, 1);
    else if (n == 2)
        put_num(s, 2, 2);
    else if (n <= 5)
        put_num(s, 0xc | (n-3), 4);
    else if (n <= 36)
        put_num(s, 0x1e0 | (n-6), 9);
    else
        put_num(s, 0xff80 | (n-37), 16);
}


static void encode_packet(J2kEncoderContext *s, J2kResLevel *rlevel, int precno)
{
    int bandno;

    // init bitstream
    *s->buf = 0;
    s->bit_index = 0;

    // header

    // is the packet empty?
    put_bits(s, 1, 1); // 1 - there are not any empty packets
    for (bandno = 0; bandno < rlevel->nbands; bandno++){
        J2kBand *band = rlevel->band + bandno;
        J2kTgtNode *cblkincl, *zerobits;
        int cblknw, cblknh, yi, xi, pos;

        cblknh = band->prec[precno].yi1 - band->prec[precno].yi0;
        cblknw = band->prec[precno].xi1 - band->prec[precno].xi0;

        cblkincl = tag_tree_init(cblknw, cblknh);
        zerobits = tag_tree_init(cblknw, cblknh);

        for (pos=0, yi = band->prec[precno].yi0; yi < band->prec[precno].yi1; yi++){
            for (xi = band->prec[precno].xi0; xi < band->prec[precno].xi1; xi++, pos++){
                cblkincl[pos].val = 0;
                tag_tree_update(cblkincl + pos);
                zerobits[pos].val = band->cblk[yi * cblknw + xi].zerobits;
                tag_tree_update(zerobits + pos);
            }
        }

        for (pos=0, yi = band->prec[precno].yi0; yi < band->prec[precno].yi1; yi++){
            for (xi = band->prec[precno].xi0; xi < band->prec[precno].xi1; xi++, pos++){
                int pad = 0, llen;
                J2kCblk *cblk = band->cblk + yi * cblknw + xi;

                // inclusion information
                tag_tree_code(s, cblkincl + pos);
                // zerobits information
                tag_tree_code(s, zerobits + pos);
                // number of passess
                putnumpassess(s, cblk->npassess);

                llen = av_log2(cblk->length) - av_log2(cblk->npassess) - 2;
                if (llen < 0){
                    pad = -llen;
                    llen = 0;
                }
                // length of code block
                put_bits(s, 1, llen);
                put_bits(s, 0, 1);
                put_num(s, cblk->length, av_log2(cblk->length)+1+pad);
            }
        }

        tag_tree_destroy(cblkincl);
        tag_tree_destroy(zerobits);
    }
    j2k_flush(s);
    for (bandno = 0; bandno < rlevel->nbands; bandno++){
        J2kBand *band = rlevel->band + bandno;
        int cblknw, yi;
        cblknw = band->prec[precno].xi1 - band->prec[precno].xi0;
        for (yi = band->prec[precno].yi0; yi < band->prec[precno].yi1; yi++){
            int xi;
            for (xi = band->prec[precno].xi0; xi < band->prec[precno].xi1; xi++){
                J2kCblk *cblk = band->cblk + yi * cblknw + xi;
                bytestream_put_buffer(&s->buf, cblk->data, cblk->length);
            }
        }
    }
}

static void encode_tile(J2kEncoderContext *s, int tileno)
{
    int compno, reslevelno, bandno;
    J2kT1Context t1;
    for (compno = 0; compno < s->ncomponents; compno++){
        J2kComponent *comp = s->tile[tileno].comp + compno;

        av_log(s->avctx, AV_LOG_DEBUG,"dwt\n");
        dwt_encode53(s, &s->tile[tileno].comp[compno]);
        av_log(s->avctx, AV_LOG_DEBUG,"after dwt -> tier1\n");

        for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
            J2kResLevel *reslevel = comp->reslevel + reslevelno;

            for (bandno = 0; bandno < reslevel->nbands ; bandno++){
                J2kBand *band = reslevel->band + bandno;
                int cblkx, cblky, cblkno=0, xx0, x0, xx1, y0, yy0, yy1;
                yy0 = bandno == 0 ? 0 : comp->reslevel[reslevelno-1].y1 - comp->reslevel[reslevelno-1].y0;
                y0 = yy0;
                yy1 = FFMIN(ceildiv(band->y0 + 1, band->cblkh) * band->cblkh, band->y1) - band->y0 + yy0;

                for (cblky = 0; cblky < band->cblkny; cblky++){
                    if (reslevelno == 0 || bandno == 1)
                        xx0 = 0;
                    else
                        xx0 = comp->reslevel[reslevelno-1].x1 - comp->reslevel[reslevelno-1].x0;
                    x0 = xx0;
                    xx1 = FFMIN(ceildiv(band->x0 + 1, band->cblkw) * band->cblkw, band->x1) - band->x0 + xx0;

                    for (cblkx = 0; cblkx < band->cblknx; cblkx++, cblkno++){
                        int y, x;
                        for (y = yy0; y < yy1; y++){
                            int *ptr = t1.data[y-yy0];
                            for (x = xx0; x < xx1; x++)
                                *ptr++ = comp->data[(comp->x1 - comp->x0) * y + x];
                        }
                        encode_cblk(s, &t1, band->cblk + cblkno, xx1 - xx0, yy1 - yy0, bandno + (reslevelno > 0?1:0));
                        xx0 = xx1;
                        xx1 = FFMIN(xx1 + band->cblkw, band->x1 - band->x0 + x0);
                    }
                    yy0 = yy1;
                    yy1 = FFMIN(yy1 + band->cblkh, band->y1 - band->y0 + y0);
                }
            }
        }
        av_free(comp->data);
        av_log(s->avctx, AV_LOG_DEBUG, "after tier1\n");
    }
    av_log(s->avctx, AV_LOG_DEBUG, "tier2\n");
    // lay-rlevel-comp-pos progression
    for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
        for (compno = 0; compno < s->ncomponents; compno++){
            int precno;
            J2kResLevel *reslevel = s->tile[tileno].comp[compno].reslevel + reslevelno;
            for (precno = 0; precno < reslevel->nprecw * reslevel->nprech; precno++){
                encode_packet(s, reslevel, precno);
            }
        }
    }
    av_log(s->avctx, AV_LOG_DEBUG, "after tier2\n");
}

void cleanup(J2kEncoderContext *s)
{
    int tileno, compno, reslevelno, bandno;
    for (tileno = 0; tileno < s->numXtiles * s->numYtiles; tileno++){
        for (compno = 0; compno < s->ncomponents; compno++){
            J2kComponent *comp = s->tile[tileno].comp + compno;

            for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
                J2kResLevel *reslevel = comp->reslevel + reslevelno;

                for (bandno = 0; bandno < reslevel->nbands ; bandno++){
                    J2kBand *band = reslevel->band + bandno;
                        av_free(band->cblk);
                        av_free(band->prec);
                    }
                av_free(reslevel->band);
            }
            av_free(comp->reslevel);
        }
        av_free(s->tile[tileno].comp);
    }
    av_free(s->tile);
}

static int encode_frame(AVCodecContext *avctx,
                        uint8_t *buf, int buf_size,
                        void *data)
{
    int tileno;
    J2kEncoderContext *s = avctx->priv_data;

    s->avctx = avctx;
    av_log(s->avctx, AV_LOG_DEBUG, "start\n");
    s->picture = data;

    // defaults:
    // TODO: implement setting non-standard precinct size
    s->ppx = 15; s->ppy = 15;

    s->XTsiz = 256; s->YTsiz = 256;
    s->nreslevels = 7;
    s->xcb = s->ycb = 4;

    // init:
    s->buf = s->buf_start = buf;
    s->buf_end = buf + buf_size;
    s->Xsiz = avctx->width;
    s->Ysiz = avctx->height;


    // TODO: other pixel formats
    if (avctx->pix_fmt == PIX_FMT_RGB24){
        s->ncomponents = 3;
        s->bpp = 24;

        // XXX: to beverified (after adding quantization)
        // (now it just works, but i'm not sure why ;-) )
        s->nguardbits = 1;
        s->expn = 8;
    }
    else{
        av_log(avctx, AV_LOG_ERROR, "only rgb24 supported\n");
    }

    put_soc(s);
    put_siz(s);
    put_cod(s);
    put_qcd(s);

    av_log(s->avctx, AV_LOG_DEBUG, "init\n");
    init_tiles(s);
    av_log(s->avctx, AV_LOG_DEBUG, "after init\n");

    for (tileno = 0; tileno < s->numXtiles * s->numYtiles; tileno++){
        uint8_t *psotptr;
        psotptr = put_sot(s, tileno);
        put_marker(s, J2K_SOD);
        encode_tile(s, tileno);
        bytestream_put_be32(&psotptr, s->buf - psotptr + 6);
    }
    put_marker(s, J2K_EOC);

    cleanup(s);
    av_log(s->avctx, AV_LOG_DEBUG, "end\n");
    return s->buf - s->buf_start;
}

AVCodec jpeg2000_encoder = {
    "j2k",
    CODEC_TYPE_VIDEO,
    CODEC_ID_JPEG2000,
    sizeof(J2kEncoderContext),
    NULL,
    encode_frame,
    NULL,
    NULL,
    0,
    .pix_fmts =
        (enum PixelFormat[]) {PIX_FMT_RGB24, -1}
};
