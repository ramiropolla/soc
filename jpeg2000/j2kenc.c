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
 */

/**
 * JPEG2000 image encoder
 * @file j2kenc.c
 * @author Kamil Nowosad
 */

#include <float.h>
#include "avcodec.h"
#include "bytestream.h"
#include "j2k.h"
#include "common.h"

#define NMSEDEC_BITS 7
#define NMSEDEC_FRACBITS (NMSEDEC_BITS-1)
#define WMSEDEC_SHIFT 13 ///< must be >= 13
#define LAMBDA_SCALE (100000000LL << (WMSEDEC_SHIFT - 13))

static int lut_nmsedec_ref [1<<NMSEDEC_BITS],
           lut_nmsedec_ref0[1<<NMSEDEC_BITS],
           lut_nmsedec_sig [1<<NMSEDEC_BITS],
           lut_nmsedec_sig0[1<<NMSEDEC_BITS];

typedef struct {
    uint16_t rate;
    int64_t disto;
} J2kPass;

typedef struct {
    uint8_t npasses;
    uint8_t ninclpasses; ///< number coding of passes included in codestream
    uint8_t nonzerobits;
    uint8_t zero;
    uint8_t data[8192];
    J2kPass passes[30];
} J2kCblk; ///< code block

typedef struct {
    uint16_t xi0, xi1, yi0, yi1; ///< indices of codeblocks ([xi0, xi1))
} J2kPrec; ///< precinct

typedef struct {
    uint16_t x0, x1, y0, y1;
    uint16_t codeblock_width, codeblock_height;
    uint16_t cblknx, cblkny;
    J2kPrec *prec;
    J2kCblk *cblk;
} J2kBand; ///< subband

typedef struct {
    uint16_t x0, x1, y0, y1;
    uint8_t nbands;
    uint16_t num_precincts_x, num_precincts_y; ///< number of precincts in x/y direction
    J2kBand *band;
} J2kResLevel; ///< resolution level

typedef struct {
   J2kResLevel *reslevel;
   int *data;
   uint16_t x0, x1, y0, y1;
} J2kComponent;

typedef struct {
   J2kComponent *comp;
} J2kTile;

typedef struct {
    AVCodecContext *avctx;
    AVFrame *picture;

    int width, height; ///< image width and height
    uint8_t cbps[4]; ///< numbps in components
    uint8_t bbps[4][32][3]; ///< numbps in bands
    uint8_t expn[4][32][3]; ///< quantization exponents
    int ncomponents;
    int log2_prec_width, log2_prec_height; ///< exponent of the precinct size [global]
    int log2_cblk_width, log2_cblk_height; ///< exponent of the code block size
    int tile_width, tile_height; ///< tile size
    int numXtiles, numYtiles;

    int nguardbits;

    int nreslevels; ///< number of resolution levels
    uint8_t *buf_start;
    uint8_t *buf;
    uint8_t *buf_end;
    int bit_index;

    int64_t lambda;

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
    fprintf(fd, "XSiz = %d, YSiz = %d, tile_width = %d, tile_height = %d\n"
                "numXtiles = %d, numYtiles = %d, ncomponents = %d\n"
                "tiles:\n",
            s->width, s->height, s->tile_width, s->tile_height,
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
                                "codeblock_width = %d, codeblock_height = %d cblknx = %d cblkny = %d\n",
                                band->x0, band->x1,
                                band->y0, band->y1,
                                band->codeblock_width, band->codeblock_height,
                                band->cblknx, band->cblkny);
                    for (precno = 0; precno < reslevel->num_precincts_x * reslevel->num_precincts_y; precno++){
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

/* bitstream routines */

/** put n times val bit */
static void put_bits(J2kEncoderContext *s, int val, int n) // TODO: optimize
{
    while (n-- > 0){
        if (s->bit_index == 8)
        {
            s->bit_index = *s->buf == 0xff;
            *(++s->buf) = 0;
        }
        *s->buf |= val << (7 - s->bit_index++);
    }
}

/** put n least significant bits of a number num */
static void put_num(J2kEncoderContext *s, int num, int n)
{
    while(--n >= 0)
        put_bits(s, (num >> n) & 1, 1);
}

/** flush the bitstream */
static void j2k_flush(J2kEncoderContext *s)
{
    if (s->bit_index){
        s->bit_index = 0;
        s->buf++;
    }
}

/* tag tree routines */

/** code the value stored in node */
static void tag_tree_code(J2kEncoderContext *s, J2kTgtNode *node, int threshold)
{
    J2kTgtNode *stack[30];
    int sp = 1, curval = 0;
    stack[0] = node;

    node = node->parent;
    while(node){
        if (node->vis){
            curval = node->val;
            break;
        }
        node->vis++;
        stack[sp++] = node;
        node = node->parent;
    }
    while(--sp >= 0){
        if (stack[sp]->val >= threshold){
            put_bits(s, 0, threshold - curval);
            break;
        }
        put_bits(s, 0, stack[sp]->val - curval);
        put_bits(s, 1, 1);
        curval = stack[sp]->val;
    }
}

/** update the value in node */
static void tag_tree_update(J2kTgtNode *node)
{
    int lev = 0;
    while (node->parent){
        if (node->parent->val <= node->val)
            break;
        node->parent->val = node->val;
        node = node->parent;
        lev++;
    }
}

static void put_siz(J2kEncoderContext *s)
{
    int i;

    bytestream_put_be16(&s->buf, J2K_SIZ);
    bytestream_put_be16(&s->buf, 38 + 3 * s->ncomponents); // Lsiz
    bytestream_put_be16(&s->buf, 0); // Rsiz
    bytestream_put_be32(&s->buf, s->width); // width
    bytestream_put_be32(&s->buf, s->height); // height
    bytestream_put_be32(&s->buf, 0); // X0Siz
    bytestream_put_be32(&s->buf, 0); // Y0Siz

    bytestream_put_be32(&s->buf, s->tile_width); // XTSiz
    bytestream_put_be32(&s->buf, s->tile_height); // YTSiz
    bytestream_put_be32(&s->buf, 0); // XT0Siz
    bytestream_put_be32(&s->buf, 0); // YT0Siz
    bytestream_put_be16(&s->buf, s->ncomponents); // CSiz

    for (i = 0; i < s->ncomponents; i++){ // Ssiz_i XRsiz_i, YRsiz_i
        bytestream_put_byte(&s->buf, 7);
        bytestream_put_byte(&s->buf, 1);
        bytestream_put_byte(&s->buf, 1);
    }
}

static void put_cod(J2kEncoderContext *s)
{
    bytestream_put_be16(&s->buf, J2K_COD);
    bytestream_put_be16(&s->buf, 12); // Lcod
    bytestream_put_byte(&s->buf, 0);  // Scod
    // SGcod
    bytestream_put_byte(&s->buf, 0); // progression level
    bytestream_put_be16(&s->buf, 1); // num of layers
    bytestream_put_byte(&s->buf, 0); // multiple component transformation
    // SPcod
    bytestream_put_byte(&s->buf, s->nreslevels - 1); // num of decomp. levels
    bytestream_put_byte(&s->buf, s->log2_cblk_width-2); // cblk width
    bytestream_put_byte(&s->buf, s->log2_cblk_height-2); // cblk height
    bytestream_put_byte(&s->buf, 0); // cblk style
    bytestream_put_byte(&s->buf, 1); // transformation
}

static void put_qcd(J2kEncoderContext *s, int compno)
{
    int reslevelno;
    bytestream_put_be16(&s->buf, J2K_QCD);
    bytestream_put_be16(&s->buf, 4+3*(s->nreslevels-1));  // LQcd
    bytestream_put_byte(&s->buf, s->nguardbits << 5);  // Sqcd
    for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
        int bandno, nbands = reslevelno == 0 ? 1:3;
        for (bandno = 0; bandno < nbands; bandno++)
            bytestream_put_byte(&s->buf, s->expn[compno][reslevelno][bandno] << 3);
    }
}

static uint8_t *put_sot(J2kEncoderContext *s, int tileno)
{
    uint8_t *psotptr;
    bytestream_put_be16(&s->buf, J2K_SOT);
    bytestream_put_be16(&s->buf, 10); // Lsot
    bytestream_put_be16(&s->buf, tileno); // Isot

    psotptr = s->buf;
    bytestream_put_be32(&s->buf, 0); // Psot (filled in later)

    bytestream_put_byte(&s->buf, 0); // TPsot
    bytestream_put_byte(&s->buf, 1); // TNsot
    return psotptr;
}

/**
 * compute the sizes of tiles, resolution levels, bands, etc.
 * allocate memory for them
 * divide the input image into tile-components
 */
static int init_tiles(J2kEncoderContext *s)
{
    int y, x, tno, compno, reslevelno, bandno, i;

    s->numXtiles = ff_j2k_ceildiv(s->width, s->tile_width);
    s->numYtiles = ff_j2k_ceildiv(s->height, s->tile_height);

    s->tile = av_malloc(s->numXtiles * s->numYtiles * sizeof(J2kTile));
    if (!s->tile)
        return AVERROR(ENOMEM);
    for (tno = 0; tno < s->numXtiles * s->numYtiles; tno++){
        J2kTile *tile = s->tile + tno;
        int p = tno % s->numXtiles;
        int q = tno / s->numXtiles;

        tile->comp = av_malloc(s->ncomponents * sizeof(J2kComponent));
        if (!tile->comp)
            return AVERROR(ENOMEM);
        for (compno = 0; compno < s->ncomponents; compno++){
            J2kComponent *comp = tile->comp + compno;

            comp->x0 = p * s->tile_width;
            comp->x1 = FFMIN((p+1)*s->tile_width, s->width);
            comp->y0 = q * s->tile_height;
            comp->y1 = FFMIN((q+1)*s->tile_height, s->height);
            comp->data = av_malloc((comp->y1 - comp->y0) * (comp->x1 -comp->x0) * sizeof(int));
            if (!comp->data)
                return AVERROR(ENOMEM);
            comp->reslevel = av_malloc(s->nreslevels * sizeof(J2kResLevel));
            if (!comp->reslevel)
                return AVERROR(ENOMEM);
            for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
                int n = s->nreslevels - reslevelno;
                J2kResLevel *reslevel = comp->reslevel + reslevelno;

                reslevel->x0 = ff_j2k_ceildivpow2(comp->x0, s->nreslevels - reslevelno - 1);
                reslevel->x1 = ff_j2k_ceildivpow2(comp->x1, s->nreslevels - reslevelno - 1);
                reslevel->y0 = ff_j2k_ceildivpow2(comp->y0, s->nreslevels - reslevelno - 1);
                reslevel->y1 = ff_j2k_ceildivpow2(comp->y1, s->nreslevels - reslevelno - 1);

                if (reslevelno == 0)
                    reslevel->nbands = 1;
                else
                    reslevel->nbands = 3;

                if (reslevel->x1 == reslevel->x0)
                    reslevel->num_precincts_x = 0;
                else
                    reslevel->num_precincts_x = ff_j2k_ceildivpow2(reslevel->x1, s->log2_prec_width) - reslevel->x0 / (1<<s->log2_prec_width);

                if (reslevel->y1 == reslevel->y0)
                    reslevel->num_precincts_y = 0;
                else
                    reslevel->num_precincts_y = ff_j2k_ceildivpow2(reslevel->y1, s->log2_prec_height) - reslevel->y0 / (1<<s->log2_prec_height);
                reslevel->band = av_malloc(reslevel->nbands * sizeof(J2kBand));
                if (!reslevel->band)
                    return AVERROR(ENOMEM);
                for (bandno = 0; bandno < reslevel->nbands; bandno++){
                    J2kBand *band = reslevel->band + bandno;
                    int cblkno, precx, precy, precno;
                    int x0, y0, x1, y1;
                    int xi0, yi0, xi1, yi1;
                    int cblkperprecw, cblkperprech;

                    if (reslevelno == 0){  // the same everywhere
                        band->codeblock_width = 1 << FFMIN(s->log2_cblk_width, s->log2_prec_width-1);
                        band->codeblock_height = 1 << FFMIN(s->log2_cblk_height, s->log2_prec_height-1);

                        band->x0 = ff_j2k_ceildivpow2(comp->x0, n-1);
                        band->x1 = ff_j2k_ceildivpow2(comp->x1, n-1);
                        band->y0 = ff_j2k_ceildivpow2(comp->y0, n-1);
                        band->y1 = ff_j2k_ceildivpow2(comp->y1, n-1);
                    }
                    else{
                        band->codeblock_width = 1 << FFMIN(s->log2_cblk_width, s->log2_prec_width);
                        band->codeblock_height = 1 << FFMIN(s->log2_cblk_height, s->log2_prec_height);

                        band->x0 = ff_j2k_ceildivpow2(comp->x0 - (1 << (n-1)) * ((bandno+1)&1), n);
                        band->x1 = ff_j2k_ceildivpow2(comp->x1 - (1 << (n-1)) * ((bandno+1)&1), n);
                        band->y0 = ff_j2k_ceildivpow2(comp->y0 - (1 << (n-1)) * (((bandno+1)&2)>>1), n);
                        band->y1 = ff_j2k_ceildivpow2(comp->y1 - (1 << (n-1)) * (((bandno+1)&2)>>1), n);
                    }

                    band->cblknx = ff_j2k_ceildiv(band->x1, band->codeblock_width) - band->x0 / band->codeblock_width;
                    band->cblkny = ff_j2k_ceildiv(band->y1, band->codeblock_height) - band->y0 / band->codeblock_height;

                    band->cblk = av_malloc(band->cblknx * band->cblkny * sizeof(J2kCblk));
                    if (!band->cblk)
                        return AVERROR(ENOMEM);
                    band->prec = av_malloc(reslevel->num_precincts_x * reslevel->num_precincts_y * sizeof(J2kPrec));
                    if (!band->prec)
                        return AVERROR(ENOMEM);

                    for (cblkno = 0; cblkno < band->cblknx * band->cblkny; cblkno++){
                        band->cblk[cblkno].zero = 0;
                    }

                    y0 = band->y0;
                    y1 = (band->y0 + (1<<s->log2_prec_height))/(1<<s->log2_prec_height)*(1<<s->log2_prec_height) - band->y0;
                    yi0 = 0;
                    yi1 = ff_j2k_ceildiv(y1 - y0, 1<<s->log2_cblk_height) * (1<<s->log2_cblk_height);
                    yi1 = FFMIN(yi1, band->cblkny);
                    cblkperprech = 1<<(s->log2_prec_height - s->log2_cblk_height);
                    for (precy = 0, precno = 0; precy < reslevel->num_precincts_y; precy++){
                        for (precx = 0; precx < reslevel->num_precincts_x; precx++, precno++){
                            band->prec[precno].yi0 = yi0;
                            band->prec[precno].yi1 = yi1;
                        }
                        yi1 += cblkperprech;
                        yi0 = yi1 - cblkperprech;
                        yi1 = FFMIN(yi1, band->cblkny);
                    }
                    x0 = band->x0;
                    x1 = (band->x0 + (1<<s->log2_prec_width))/(1<<s->log2_prec_width)*(1<<s->log2_prec_width) - band->x0;
                    xi0 = 0;
                    xi1 = ff_j2k_ceildiv(x1 - x0, 1<<s->log2_cblk_width) * (1<<s->log2_cblk_width);
                    xi1 = FFMIN(xi1, band->cblknx);
                    cblkperprecw = 1<<(s->log2_prec_width - s->log2_cblk_width);
                    for (precx = 0, precno = 0; precx < reslevel->num_precincts_x; precx++){
                        for (precy = 0; precy < reslevel->num_precincts_y; precy++, precno = 0){
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
                for (compno = 0; compno < s->ncomponents; compno++){
                    tile->comp[compno].data[i] = *ptr++  - (1 << 7);
                }
            }
            line += s->picture->linesize[0];
        }
    }
    // calculate band bps and exponents
    for (compno = 0; compno < s->ncomponents; compno++){
        for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
            int nbands;
            nbands = reslevelno ? 3 : 1;
            for (bandno = 0; bandno < nbands; bandno++){
                int expn;

                expn = ((bandno&2)>>1) + (reslevelno>0) + s->cbps[compno];
                s->bbps[compno][reslevelno][bandno] = expn + s->nguardbits - 1;
                s->expn[compno][reslevelno][bandno] = expn;
            }
        }
    }
    return 0;
}

static void init_luts()
{
    int i;
    double u, v, t, pfr;

    pfr = pow(2, NMSEDEC_FRACBITS);
    for (i = 0; i < (1 << NMSEDEC_BITS); i++){
        t = i / pfr;
        u = t;
        v = t - 1.5;
        lut_nmsedec_sig[i]  = FFMAX((int) (floor((u*u - v*v) * pfr + 0.5) / pfr * 8192.0), 0);
        lut_nmsedec_sig0[i] = FFMAX((int) (floor((u*u) * pfr + 0.5) / pfr * 8192.0), 0);

        u = t - 1.0;
        v = t - ((i & (1<<(NMSEDEC_BITS-1))) ? 1.5 : 0.5);
        lut_nmsedec_ref[i]  = FFMAX((int) (floor((u*u - v*v) * pfr + 0.5) / pfr * 8192.0), 0);
        lut_nmsedec_ref0[i] = FFMAX((int) (floor((u*u) * pfr + 0.5) / pfr * 8192.0), 0);
    }
}

/* discrete wavelet transform routines */
static void sd_1d(int *p, int i0, int i1)
{
    int i;

    if (i1 == i0 + 1)
        return;

    p[i0 - 1] = p[i0 + 1];
    p[i1    ] = p[i1 - 2];
    p[i0 - 2] = p[i0 + 2];
    p[i1 + 1] = p[i1 - 3];

    for (i = (i0+1)/2 - 1; i < (i1+1)/2; i++){
        p[2*i+1] -= (p[2*i] + p[2*i+2]) >> 1;
    }
    for (i = (i0+1)/2; i < (i1+1)/2; i++){
        p[2*i] += (p[2*i-1] + p[2*i+1] + 2) >> 2;
    }
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

        //VER_SD
        while (u < u1){
            int i, j;
            for (i = v0; i < v1; i++)
                pv[i] = t[w*(i-v0) + u-u0];
            sd_1d(pv, v0, v1);

            // copy back and deinterleave
            for (i = v0+v0%2, j = 0; i < v1; i+=2, j++){
                t[w*j + u-u0] = pv[i];
            }
            for (i = v0+1-v0%2; i < v1; i+=2, j++){
                t[w*j + u-u0] = pv[i];
            }
            u++;
        }

        //HOR_SD
        while (v < v1){
            int i, j;
            for (i = u0; i < u1; i++)
                pu[i] = t[w*(v-v0) + i-u0];
            sd_1d(pu, u0, u1);

            // copy back and deinterleave
            for (i = u0+u0%2, j = 0; i < u1; i+=2, j++){
                t[w*(v-v0) + j] = pu[i];
            }
            for (i = u0+1-u0%2; i < u1; i+=2, j++){
                t[w*(v-v0) + j] = pu[i];
            }
            v++;
        }
    }
    av_free(ppv);
    av_free(ppu);
}

/* tier-1 routines */
static int getnmsedec_sig(int x, int bpno)
{
    if (bpno > NMSEDEC_FRACBITS)
        return lut_nmsedec_sig[(x >> (bpno - NMSEDEC_FRACBITS)) & ((1 << NMSEDEC_BITS) - 1)];
    return lut_nmsedec_sig0[x & ((1 << NMSEDEC_BITS) - 1)];
}

static int getnmsedec_ref(int x, int bpno)
{
    if (bpno > NMSEDEC_FRACBITS)
        return lut_nmsedec_ref[(x >> (bpno - NMSEDEC_FRACBITS)) & ((1 << NMSEDEC_BITS) - 1)];
    return lut_nmsedec_ref0[x & ((1 << NMSEDEC_BITS) - 1)];
}

static void encode_sigpass(J2kT1Context *t1, int width, int height, int bandno, int *nmsedec, int bpno)
{
    int y0, x, y, mask = 1 << (bpno + NMSEDEC_FRACBITS);
    for (y0 = 0; y0 < height; y0 += 4)
        for (x = 0; x < width; x++)
            for (y = y0; y < height && y < y0+4; y++){
                if (!(t1->flags[y+1][x+1] & J2K_T1_SIG) && (t1->flags[y+1][x+1] & J2K_T1_SIG_NB)){
                    int ctxno = ff_j2k_getnbctxno(t1->flags[y+1][x+1], bandno),
                        bit = abs(t1->data[y][x]) & mask ? 1 : 0;
                    ff_aec_encode(&t1->aec, t1->aec.cx_states + ctxno, bit);
                    if (bit){
                        int xorbit;
                        int ctxno = ff_j2k_getsgnctxno(t1->flags[y+1][x+1], &xorbit);
                        ff_aec_encode(&t1->aec, t1->aec.cx_states + ctxno, (t1->data[y][x] < 0) ^ xorbit);
                        *nmsedec += getnmsedec_sig(abs(t1->data[y][x]), bpno + NMSEDEC_FRACBITS);
                        ff_j2k_set_significant(t1, x, y);
                    }
                    t1->flags[y+1][x+1] |= J2K_T1_VIS;
                }
            }
}

static void encode_refpass(J2kT1Context *t1, int width, int height, int *nmsedec, int bpno)
{
    int y0, x, y, mask = 1 << (bpno + NMSEDEC_FRACBITS);
    for (y0 = 0; y0 < height; y0 += 4)
        for (x = 0; x < width; x++)
            for (y = y0; y < height && y < y0+4; y++)
                if ((t1->flags[y+1][x+1] & (J2K_T1_SIG | J2K_T1_VIS)) == J2K_T1_SIG){
                    int ctxno = ff_j2k_getrefctxno(t1->flags[y+1][x+1]);
                    *nmsedec += getnmsedec_ref(abs(t1->data[y][x]), bpno + NMSEDEC_FRACBITS);
                    ff_aec_encode(&t1->aec, t1->aec.cx_states + ctxno, abs(t1->data[y][x]) & mask ? 1:0);
                    t1->flags[y+1][x+1] |= J2K_T1_REF;
                }
}

static void encode_clnpass(J2kT1Context *t1, int width, int height, int bandno, int *nmsedec, int bpno)
{
    int y0, x, y, mask = 1 << (bpno + NMSEDEC_FRACBITS);
    for (y0 = 0; y0 < height; y0 += 4)
        for (x = 0; x < width; x++){
            if (y0 + 3 < height && !(
            (t1->flags[y0+1][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[y0+2][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[y0+3][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[y0+4][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG))))
            {
                // aggregation mode
                int rlen;
                for (rlen = 0; rlen < 4; rlen++)
                    if (abs(t1->data[y0+rlen][x]) & mask)
                        break;
                ff_aec_encode(&t1->aec, t1->aec.cx_states + AEC_CX_RL, rlen != 4);
                if (rlen == 4)
                    continue;
                ff_aec_encode(&t1->aec, t1->aec.cx_states + AEC_CX_UNI, rlen >> 1);
                ff_aec_encode(&t1->aec, t1->aec.cx_states + AEC_CX_UNI, rlen & 1);
                for (y = y0 + rlen; y < y0 + 4; y++){
                    if (!(t1->flags[y+1][x+1] & (J2K_T1_SIG | J2K_T1_VIS))){
                        int ctxno = ff_j2k_getnbctxno(t1->flags[y+1][x+1], bandno);
                        if (y > y0 + rlen)
                            ff_aec_encode(&t1->aec, t1->aec.cx_states + ctxno, abs(t1->data[y][x]) & mask ? 1:0);
                        if (abs(t1->data[y][x]) & mask){ // newly significant
                            int xorbit;
                            int ctxno = ff_j2k_getsgnctxno(t1->flags[y+1][x+1], &xorbit);
                            *nmsedec += getnmsedec_sig(abs(t1->data[y][x]), bpno + NMSEDEC_FRACBITS);
                            ff_aec_encode(&t1->aec, t1->aec.cx_states + ctxno, (t1->data[y][x] < 0) ^ xorbit);
                            ff_j2k_set_significant(t1, x, y);
                        }
                    }
                    t1->flags[y+1][x+1] &= ~J2K_T1_VIS;
                }
            }
            else{
                for (y = y0; y < y0 + 4 && y < height; y++){
                    if (!(t1->flags[y+1][x+1] & (J2K_T1_SIG | J2K_T1_VIS))){
                        int ctxno = ff_j2k_getnbctxno(t1->flags[y+1][x+1], bandno);
                        ff_aec_encode(&t1->aec, t1->aec.cx_states + ctxno, abs(t1->data[y][x]) & mask ? 1:0);
                        if (abs(t1->data[y][x]) & mask){ // newly significant
                            int xorbit;
                            int ctxno = ff_j2k_getsgnctxno(t1->flags[y+1][x+1], &xorbit);
                            *nmsedec += getnmsedec_sig(abs(t1->data[y][x]), bpno + NMSEDEC_FRACBITS);
                            ff_aec_encode(&t1->aec, t1->aec.cx_states + ctxno, (t1->data[y][x] < 0) ^ xorbit);
                            ff_j2k_set_significant(t1, x, y);
                        }
                    }
                    t1->flags[y+1][x+1] &= ~J2K_T1_VIS;
                }
            }
        }
}

static void encode_cblk(J2kEncoderContext *s, J2kT1Context *t1, J2kCblk *cblk, J2kTile *tile,
                        int width, int height, int bandpos, int lev)
{
    int pass_t = 2, passno, x, y, max=0, nmsedec, bpno;
    int64_t wmsedec = 0;

    for (y = 0; y < height+2; y++)
        memset(t1->flags[y], 0, (width+2)*sizeof(int));

    for (y = 0; y < height; y++){
        for (x = 0; x < width; x++)
            max = FFMAX(max, abs(t1->data[y][x]));
    }

    if (max == 0){
        cblk->nonzerobits = 0;
        bpno = 0;
    }
    else{
        cblk->nonzerobits = av_log2(max) + 1 - NMSEDEC_FRACBITS;
        bpno = cblk->nonzerobits - 1;
    }

    ff_aec_initenc(&t1->aec, cblk->data);

    for (passno = 0; bpno >= 0; passno++){
        nmsedec=0;

        switch(pass_t){
            case 0: encode_sigpass(t1, width, height, bandpos, &nmsedec, bpno);
                    break;
            case 1: encode_refpass(t1, width, height, &nmsedec, bpno);
                    break;
            case 2: encode_clnpass(t1, width, height, bandpos, &nmsedec, bpno);
                    break;
        }

        cblk->passes[passno].rate = 3 + ff_aec_length(&t1->aec);
        wmsedec += (int64_t)nmsedec << (2*bpno);
        cblk->passes[passno].disto = wmsedec;

        if (++pass_t == 3){
            pass_t = 0;
            bpno--;
        }
    }
    cblk->npasses = passno;
    cblk->ninclpasses = passno;

    // TODO: optional flush on each pass
    cblk->passes[passno-1].rate = ff_aec_flush(&t1->aec);
}

/* tier-2 routines: */

static void putnumpasses(J2kEncoderContext *s, int n)
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


static void encode_packet(J2kEncoderContext *s, J2kResLevel *rlevel, int precno, int compno, int rlevelno)
{
    int bandno, empty = 1;

    // init bitstream
    *s->buf = 0;
    s->bit_index = 0;

    // header

    // is the packet empty?
    for (bandno = 0; bandno < rlevel->nbands; bandno++){
        if (rlevel->band[bandno].x0 < rlevel->band[bandno].x1
        &&  rlevel->band[bandno].y0 < rlevel->band[bandno].y1){
            empty = 0;
            break;
        }
    }
    if (empty){
        put_bits(s, 0, 1);
        j2k_flush(s);
        return;
    }

    put_bits(s, 1, 1);
    for (bandno = 0; bandno < rlevel->nbands; bandno++){
        J2kBand *band = rlevel->band + bandno;
        J2kTgtNode *cblkincl, *zerobits;
        int cblknw, cblknh, yi, xi, pos;

        cblknh = band->prec[precno].yi1 - band->prec[precno].yi0;
        cblknw = band->prec[precno].xi1 - band->prec[precno].xi0;

        cblkincl = ff_j2k_tag_tree_init(cblknw, cblknh);
        zerobits = ff_j2k_tag_tree_init(cblknw, cblknh);

        for (pos=0, yi = band->prec[precno].yi0; yi < band->prec[precno].yi1; yi++){
            for (xi = band->prec[precno].xi0; xi < band->prec[precno].xi1; xi++, pos++){
                cblkincl[pos].val = band->cblk[yi * cblknw + xi].ninclpasses == 0;
                tag_tree_update(cblkincl + pos);
                zerobits[pos].val = s->bbps[compno][rlevelno][bandno] - band->cblk[yi * cblknw + xi].nonzerobits;
                tag_tree_update(zerobits + pos);
            }
        }

        for (pos=0, yi = band->prec[precno].yi0; yi < band->prec[precno].yi1; yi++){
            for (xi = band->prec[precno].xi0; xi < band->prec[precno].xi1; xi++, pos++){
                int pad = 0, llen, length;
                J2kCblk *cblk = band->cblk + yi * cblknw + xi;

                // inclusion information
                tag_tree_code(s, cblkincl + pos, 1);
                if (!cblk->ninclpasses)
                    continue;
                // zerobits information
                tag_tree_code(s, zerobits + pos, 100);
                // number of passes
                putnumpasses(s, cblk->ninclpasses);

                length = cblk->passes[cblk->ninclpasses-1].rate;
                llen = av_log2(length) - av_log2(cblk->ninclpasses) - 2;
                if (llen < 0){
                    pad = -llen;
                    llen = 0;
                }
                // length of code block
                put_bits(s, 1, llen);
                put_bits(s, 0, 1);
                put_num(s, length, av_log2(length)+1+pad);
            }
        }

        av_free(cblkincl);
        av_free(zerobits);
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
                if (cblk->ninclpasses)
                    bytestream_put_buffer(&s->buf, cblk->data, cblk->passes[cblk->ninclpasses-1].rate);
            }
        }
    }
}

static void encode_packets(J2kEncoderContext *s, J2kTile *tile, int tileno)
{
    int compno, reslevelno;

    av_log(s->avctx, AV_LOG_DEBUG, "tier2\n");
    // lay-rlevel-comp-pos progression
    for (reslevelno = 0; reslevelno < s->nreslevels; reslevelno++){
        for (compno = 0; compno < s->ncomponents; compno++){
            int precno;
            J2kResLevel *reslevel = s->tile[tileno].comp[compno].reslevel + reslevelno;
            for (precno = 0; precno < reslevel->num_precincts_x * reslevel->num_precincts_y; precno++){
                encode_packet(s, reslevel, precno, compno, reslevelno);
            }
        }
    }
    av_log(s->avctx, AV_LOG_DEBUG, "after tier2\n");
}

static int getcut(J2kCblk *cblk, int64_t lambda, int dwt_norm)
{
    int passno, res = 0;
    for (passno = 0; passno < cblk->npasses; passno++){
        int dr;
        int64_t dd;

        dr = cblk->passes[passno].rate
           - (res ? cblk->passes[res-1].rate:0);
        dd = cblk->passes[passno].disto
           - (res ? cblk->passes[res-1].disto:0);

        if (((dd * dwt_norm) >> WMSEDEC_SHIFT) * dwt_norm >= dr * lambda)
            res = passno+1;
    }
    return res;
}

static void truncpasses(J2kEncoderContext *s, J2kTile *tile)
{
    static const int dwt_norms[4][10] = { // multiplied by 10000
        {10000, 15000, 27500, 53750, 106800, 213400, 426700, 853300, 1707000, 3413000},
        {10380, 15920, 29190, 57030, 113300, 226400, 452500, 904800, 1809000},
        {10380, 15920, 29190, 57030, 113300, 226400, 452500, 904800, 1809000},
        { 7186,  9218, 15860, 30430,  60190, 120100, 240000, 479700,  959300}};
    int compno, reslevelno, bandno, cblkno, lev;
    for (compno = 0; compno < s->ncomponents; compno++){
        J2kComponent *comp = tile->comp + compno;

        for (reslevelno = 0, lev = s->nreslevels-1; reslevelno < s->nreslevels; reslevelno++, lev--){
            J2kResLevel *reslevel = comp->reslevel + reslevelno;

            for (bandno = 0; bandno < reslevel->nbands ; bandno++){
                int bandpos = bandno + (reslevelno > 0);
                J2kBand *band = reslevel->band + bandno;

                for (cblkno = 0; cblkno < band->cblknx * band->cblkny; cblkno++){
                    J2kCblk *cblk = band->cblk + cblkno;

                    cblk->ninclpasses = getcut(cblk, s->lambda, dwt_norms[bandpos][lev]);
                }
            }
        }
    }
}

static void encode_tile(J2kEncoderContext *s, J2kTile *tile, int tileno)
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
                int cblkx, cblky, cblkno=0, xx0, x0, xx1, y0, yy0, yy1, bandpos;
                yy0 = bandno == 0 ? 0 : comp->reslevel[reslevelno-1].y1 - comp->reslevel[reslevelno-1].y0;
                y0 = yy0;
                yy1 = FFMIN(ff_j2k_ceildiv(band->y0 + 1, band->codeblock_height) * band->codeblock_height, band->y1) - band->y0 + yy0;

                if (band->x0 == band->x1 || band->y0 == band->y1)
                    continue;

                bandpos = bandno + (reslevelno > 0);

                for (cblky = 0; cblky < band->cblkny; cblky++){
                    if (reslevelno == 0 || bandno == 1)
                        xx0 = 0;
                    else
                        xx0 = comp->reslevel[reslevelno-1].x1 - comp->reslevel[reslevelno-1].x0;
                    x0 = xx0;
                    xx1 = FFMIN(ff_j2k_ceildiv(band->x0 + 1, band->codeblock_width) * band->codeblock_width, band->x1) - band->x0 + xx0;

                    for (cblkx = 0; cblkx < band->cblknx; cblkx++, cblkno++){
                        int y, x;
                        for (y = yy0; y < yy1; y++){
                            int *ptr = t1.data[y-yy0];
                            for (x = xx0; x < xx1; x++)
                                *ptr++ = comp->data[(comp->x1 - comp->x0) * y + x] << NMSEDEC_FRACBITS;
                        }
                        encode_cblk(s, &t1, band->cblk + cblkno, tile, xx1 - xx0, yy1 - yy0,
                                    bandpos, s->nreslevels - reslevelno - 1);
                        xx0 = xx1;
                        xx1 = FFMIN(xx1 + band->codeblock_width, band->x1 - band->x0 + x0);
                    }
                    yy0 = yy1;
                    yy1 = FFMIN(yy1 + band->codeblock_height, band->y1 - band->y0 + y0);
                }
            }
        }
        av_free(comp->data);
        av_log(s->avctx, AV_LOG_DEBUG, "after tier1\n");
    }

    av_log(s->avctx, AV_LOG_DEBUG, "rate control\n");
    truncpasses(s, tile);
    encode_packets(s, tile, tileno);
    av_log(s->avctx, AV_LOG_DEBUG, "after rate control\n");
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
    int tileno, i, ret;
    J2kEncoderContext *s = avctx->priv_data;

    s->avctx = avctx;
    av_log(s->avctx, AV_LOG_DEBUG, "start\n");
    s->picture = data;

    // defaults:
    // TODO: implement setting non-standard precinct size
    s->log2_prec_width = 15; s->log2_prec_height = 15;

    s->tile_width = 256; s->tile_height = 256;
    s->nreslevels = 7;
    s->log2_cblk_width = s->log2_cblk_height = 4;

    // init:
    s->buf = s->buf_start = buf;
    s->buf_end = buf + buf_size;
    s->width = avctx->width;
    s->height = avctx->height;

    s->nguardbits = 1;
    s->lambda = s->picture->quality * LAMBDA_SCALE;

    ff_j2k_init_tier1_luts();

    // TODO: other pixel formats
    for (i = 0; i < 3; i++)
        s->cbps[i] = 8;

    if (avctx->pix_fmt == PIX_FMT_RGB24){
        s->ncomponents = 3;
    } else if (avctx->pix_fmt == PIX_FMT_GRAY8){
        s->ncomponents = 1;
    }
    else{
        av_log(avctx, AV_LOG_ERROR, "only rgb24 and gray8 supported\n");
        return -1;
    }

    av_log(s->avctx, AV_LOG_DEBUG, "init\n");
    if (ret=init_tiles(s))
        return ret;
    init_luts();
    av_log(s->avctx, AV_LOG_DEBUG, "after init\n");

    bytestream_put_be16(&s->buf, J2K_SOC);
    put_siz(s);
    put_cod(s);
    put_qcd(s, 0);

    for (tileno = 0; tileno < s->numXtiles * s->numYtiles; tileno++){
        uint8_t *psotptr;
        psotptr = put_sot(s, tileno);
        bytestream_put_be16(&s->buf, J2K_SOD);
        encode_tile(s, s->tile + tileno, tileno);
        bytestream_put_be32(&psotptr, s->buf - psotptr + 6);
    }
    bytestream_put_be16(&s->buf, J2K_EOC);

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
        (enum PixelFormat[]) {PIX_FMT_GRAY8, PIX_FMT_RGB24, -1}
};
