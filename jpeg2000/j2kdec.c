/*
 * JPEG2000 image decoder
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
 * JPEG2000 image decoder
 * @file j2kdec.c
 * @author Kamil Nowosad
 */

#include <float.h>
#include "avcodec.h"
#include "bytestream.h"
#include "j2k.h"
#include "common.h"

/** code block */
typedef struct {
    uint8_t npassess;
    uint8_t nonzerobits;
    uint16_t lengthinc;
    uint16_t length;
    uint8_t lblock;
    uint8_t zero;
    uint8_t data[8192];
    uint8_t included;
} J2kCblk;

/** precinct */
typedef struct {
    uint16_t xi0, xi1, yi0, yi1;
    J2kTgtNode *zerobits;
    J2kTgtNode *cblkincl;
} J2kPrec;

/** subband */
typedef struct {
    uint16_t x0, x1, y0, y1;
    uint16_t cblkw, cblkh;
    uint16_t cblknx, cblkny;
    J2kPrec *prec;
    J2kCblk *cblk;
} J2kBand;

/** resolution level */
typedef struct {
    uint16_t x0, x1, y0, y1;
    uint8_t nbands;
    uint16_t nprecw, nprech;
    uint8_t ppx, ppy;
    J2kBand *band;
} J2kResLevel;

typedef struct {
   J2kResLevel *reslevel;
   int *data;
   uint16_t x0, x1, y0, y1;

   uint8_t properties;

   /// COx fields
   int nreslevels;
   int xcb, ycb;
   int csty;
   int ppx, ppy;

   /// QCx fields
   uint8_t expn[32 * 3];
   uint8_t bbps[32 * 3];
   uint8_t nguardbits;
} J2kComponent;

#define HAD_COC 0x01
#define HAD_QCC 0x02

typedef struct {
   J2kComponent *comp;
   int properties;
   int nlayers;
   int progression;
} J2kTile;

typedef struct {
    AVCodecContext *avctx;
    AVFrame picture;

    int Xsiz, Ysiz; ///< image width and height
    int X0siz, Y0siz;
    int XT0siz, YT0siz;
    uint8_t cbps[4]; ///< numbps in components
    uint8_t sgnd[4]; ///< if a component is signed
    uint8_t bbps[3 * 32][4]; ///< numbps in bands
    uint8_t expn[3 * 32][4]; ///< quantization exponents
    int properties[4];

    int ncomponents;
    int XTsiz, YTsiz; ///< tile size
    int numXtiles, numYtiles;
    int maxtilelen;

    int nguardbits[4];

    int nreslevels[4]; ///< number of resolution levels
    int xcb[4], ycb[4]; ///< exponent of the code block size
    int ppx, ppy; ///< exponent of the precinct size

    uint8_t *buf_start;
    uint8_t *buf;
    uint8_t *buf_end;
    int bit_index;

    int16_t curtileno;
    int nlayers;

    int csty[4]; ///< coding styles for components

    J2kTile *tile;
} J2kDecoderContext;

static int get_bits(J2kDecoderContext *s, int n)
{
    int res = 0;
    while (--n >= 0){
        res <<= 1;
        if (s->bit_index == 0){
            s->bit_index = 7 + (*s->buf != 0xff);
            s->buf++;
        }
        s->bit_index--;
        res |= (*s->buf >> s->bit_index) & 1;
    }
    return res;
}

void j2k_flush(J2kDecoderContext *s)
{
    if (s->bit_index != 8){
        s->bit_index = 8;
        s->buf++;
    }
}
#if 0
void printcomp(J2kComponent *comp)
{
    int i;
    for (i = 0; i < comp->y1 - comp->y0; i++)
        ff_j2k_printv(comp->data + i * (comp->x1 - comp->x0), comp->x1 - comp->x0);
}

static void nspaces(FILE *fd, int n)
{
    while(n--) putc(' ', fd);
}

static void dump(J2kDecoderContext *s, FILE *fd)
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
            for(reslevelno = 0; reslevelno < comp->nreslevels; reslevelno++){
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

/** decode the value stored in node */
static int tag_tree_decode(J2kDecoderContext *s, J2kTgtNode *node, int threshold)
{
    J2kTgtNode *stack[30];
    int sp = -1, curval = 0;

    while(node && !node->vis){
        stack[++sp] = node;
        node = node->parent;
    }

    if (stack[sp]->val >= threshold)
        return threshold;

    if (node)
        curval = node->val;
    else
        curval = stack[sp]->val;

    while(sp >= 0){
        while (!get_bits(s, 1)){
            curval++;
            if (curval == threshold)
                break;
        }
        stack[sp]->val = curval;
        if (curval == threshold)
            break;
        stack[sp]->vis++;
        sp--;
    }
    return curval;
}

static void copy_defaults(J2kDecoderContext *s, J2kTile *tile)
{
    int compno;

    tile->nlayers = s->nlayers;

    for (compno = 0; compno < s->ncomponents; compno++){
        J2kComponent *comp = tile->comp + compno;
        int i;

        comp->nreslevels = s->nreslevels[compno];
        comp->xcb = s->xcb[compno];
        comp->ycb = s->ycb[compno];
        for (i = 0; i < 3*32; i++){
            comp->expn[i] = s->expn[i][compno];
            comp->bbps[i] = s->bbps[i][compno];
        }
        comp->nguardbits = s->nguardbits[compno];
    }
}

/** marker segments */
static int get_siz(J2kDecoderContext *s)
{
    int i;

    bytestream_get_be16(&s->buf); ///< Rsiz (skipped)
    s->Xsiz = bytestream_get_be32(&s->buf); ///< Xsiz
    s->Ysiz = bytestream_get_be32(&s->buf); ///< Ysiz
    s->X0siz = bytestream_get_be32(&s->buf); ///< X0Siz
    s->Y0siz = bytestream_get_be32(&s->buf); ///< Y0Siz

    s->XTsiz = bytestream_get_be32(&s->buf); ///< XTSiz
    s->YTsiz = bytestream_get_be32(&s->buf); ///< YTSiz
    s->XT0siz = bytestream_get_be32(&s->buf); ///< XT0Siz
    s->YT0siz = bytestream_get_be32(&s->buf); ///< YT0Siz
    s->ncomponents = bytestream_get_be16(&s->buf); ///< CSiz

    for (i = 0; i < s->ncomponents; i++){ ///< Ssiz_i XRsiz_i, YRsiz_i
        uint8_t x = bytestream_get_byte(&s->buf);
        s->cbps[i] = (x & 0x7f) + 1;
        s->sgnd[i] = (x & 0x80) == 1;
        if (bytestream_get_byte(&s->buf) != 1)
            return -1;
        if (bytestream_get_byte(&s->buf) != 1)
            return -1;
    }

    s->numXtiles = ff_j2k_ceildiv(s->Xsiz - s->XT0siz, s->XTsiz);
    s->numYtiles = ff_j2k_ceildiv(s->Ysiz - s->YT0siz, s->YTsiz);

    s->tile = av_mallocz(s->numXtiles * s->numYtiles * sizeof(J2kTile));
    for (i = 0; i < s->numXtiles * s->numYtiles; i++){
        J2kTile *tile = s->tile + i;

        tile->comp = av_mallocz(s->ncomponents * sizeof(J2kComponent));
    }

    s->avctx->width = s->Xsiz - s->X0siz;
    s->avctx->height = s->Ysiz - s->Y0siz;

    s->avctx->get_buffer(s->avctx, &s->picture);

    return 0;
}

#define SETFIELD(field, val)\
    if (s->curtileno == -1)\
        s->field = val;\
    else\
        s->tile[s->curtileno].field = val;

#define GETFIELD(field)\
    (s->curtileno == -1 ? s->field : s->tile[s->curtileno].field)

#define SETFIELDC(field, val)\
    if (s->curtileno == -1)\
        s->field[compno] = val;\
    else\
        s->tile[s->curtileno].comp[compno].field = val;

#define GETFIELDC(field)\
    (s->curtileno == -1 ? s->field[compno] : s->tile[s->curtileno].comp[compno].field)

static int get_cox(J2kDecoderContext *s, int compno)
{
    SETFIELDC(nreslevels, bytestream_get_byte(&s->buf) + 1); ///< num of resolution levels - 1
    SETFIELDC(xcb, bytestream_get_byte(&s->buf) + 2); ///< cblk width
    SETFIELDC(ycb, bytestream_get_byte(&s->buf) + 2); ///< cblk height
    if (bytestream_get_byte(&s->buf) != 0){ ///< cblk style
        av_log(s->avctx, AV_LOG_ERROR, "no extra cblk styles supported\n");
        return -1;
    }
    if (bytestream_get_byte(&s->buf) != 1){ ///< transformation
        av_log(s->avctx, AV_LOG_ERROR, "only DWT 5-3 supported\n");
        return -1;
    }
    return 0;
}

static int get_cod(J2kDecoderContext *s)
{
    uint8_t *pos;
    int compno, csty;

    csty = bytestream_get_byte(&s->buf);
    for (compno = 0; compno < s->ncomponents; compno++)
        if (!(GETFIELDC(properties) & HAD_COC)){
            SETFIELDC(csty, csty); ///< Scod
        }

    if (bytestream_get_byte(&s->buf)){ ///< progression level
        av_log(s->avctx, AV_LOG_ERROR, "only LRCP progression supported\n");
        return -1;
    }

    SETFIELD(nlayers, bytestream_get_be16(&s->buf));
    if (bytestream_get_byte(&s->buf)){ ///< multiple component transformation
        av_log(s->avctx, AV_LOG_ERROR, "MCT not supported\n");
    }
    pos = s->buf;
    for (compno = 0; compno < s->ncomponents; compno++)
        if (!(GETFIELDC(properties) & HAD_COC)){
            s->buf = pos;
            get_cox(s, compno);
        }
    return 0;
}

static int get_coc(J2kDecoderContext *s)
{
    int compno = bytestream_get_byte(&s->buf);
    SETFIELDC(csty, bytestream_get_byte(&s->buf));
    get_cox(s, compno);

    if (s->curtileno == -1)
        s->properties[compno] |= HAD_COC;
    else
        s->tile[s->curtileno].comp[compno].properties |= HAD_COC;
    return 0;
}

static int get_qcx(J2kDecoderContext *s, int n, int compno)
{
    int i, x;
    x = bytestream_get_byte(&s->buf); ///< Sqcd
    SETFIELDC(nguardbits, x >> 5);

    if (x & 0x1f){
        av_log(s->avctx, AV_LOG_ERROR, "no quantization supported\n");
        return -1;
    }
    n -= 3;

    for (i = 0; i < n; i++){
        x = bytestream_get_byte(&s->buf);
        SETFIELDC(expn[i], x >> 3);
        SETFIELDC(bbps[i], (x >> 3) + GETFIELDC(nguardbits) - 1);
    }
    return 0;
}

static int get_qcd(J2kDecoderContext *s, int n)
{
    uint8_t *pos = s->buf;
    int compno;
    for (compno = 0; compno < s->ncomponents; compno++)
        if (!(GETFIELDC(properties) & HAD_QCC)){
            s->buf = pos;
            if (get_qcx(s, n, compno))
                return -1;
        }
    return 0;
}

static int get_qcc(J2kDecoderContext *s, int n)
{
    int compno = bytestream_get_byte(&s->buf);
    if (s->curtileno == -1)
        s->properties[compno] |= HAD_QCC;
    else
        s->tile[s->curtileno].comp[compno].properties |= HAD_QCC;
    return get_qcx(s, n-1, compno);
}

static uint8_t get_sot(J2kDecoderContext *s)
{
    s->curtileno = bytestream_get_be16(&s->buf); ///< Isot

    s->buf += 4; ///< Psot (ignored)

    if (!bytestream_get_byte(&s->buf)) ///< TPsot
        copy_defaults(s, s->tile + s->curtileno);
    bytestream_get_byte(&s->buf); ///< TNsot

    return 0;
}

#undef SETFIELDC
#undef SETFIELD
#undef GETFIELDC
#undef GETFIELD

static int init_tile(J2kDecoderContext *s, int tileno)
{
    int compno, reslevelno, bandno, p, q;
    J2kTile *tile = s->tile + tileno;

    p = tileno % s->numXtiles;
    q = tileno / s->numXtiles;

    if (tile->comp == NULL)
        return -1;
    for (compno = 0; compno < s->ncomponents; compno++){
        J2kComponent *comp = tile->comp + compno;

        comp->x0 = FFMAX(p * s->XTsiz + s->XT0siz, s->X0siz);
        comp->x1 = FFMIN((p+1)*s->XTsiz + s->XT0siz, s->Xsiz);
        comp->y0 = FFMAX(q * s->YTsiz + s->YT0siz, s->Y0siz);
        comp->y1 = FFMIN((q+1)*s->YTsiz + s->YT0siz, s->Ysiz);

        comp->data = av_malloc((comp->y1 - comp->y0) * (comp->x1 -comp->x0) * sizeof(int));
        if (comp->data == NULL)
            return -1;
        comp->reslevel = av_malloc(comp->nreslevels * sizeof(J2kResLevel));
        if (comp->reslevel == NULL)
            return -1;
        for (reslevelno = 0; reslevelno < comp->nreslevels; reslevelno++){
            int n = comp->nreslevels - reslevelno;
            J2kResLevel *reslevel = comp->reslevel + reslevelno;

            reslevel->x0 = ff_j2k_ceildivpow2(comp->x0, comp->nreslevels - reslevelno - 1);
            reslevel->x1 = ff_j2k_ceildivpow2(comp->x1, comp->nreslevels - reslevelno - 1);
            reslevel->y0 = ff_j2k_ceildivpow2(comp->y0, comp->nreslevels - reslevelno - 1);
            reslevel->y1 = ff_j2k_ceildivpow2(comp->y1, comp->nreslevels - reslevelno - 1);

            if (reslevelno == 0)
                reslevel->nbands = 1;
            else
                reslevel->nbands = 3;

            if (reslevel->x1 == reslevel->x0)
                reslevel->nprecw = 0;
            else
                reslevel->nprecw = ff_j2k_ceildivpow2(reslevel->x1, s->ppx) - reslevel->x0 / (1<<s->ppx);

            if (reslevel->y1 == reslevel->y0)
                reslevel->nprech = 0;
            else
                reslevel->nprech = ff_j2k_ceildivpow2(reslevel->y1, s->ppy) - reslevel->y0 / (1<<s->ppy);

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
                    band->cblkw = 1 << FFMIN(comp->xcb, s->ppx-1);
                    band->cblkh = 1 << FFMIN(comp->ycb, s->ppy-1);

                    band->x0 = ff_j2k_ceildivpow2(comp->x0, n-1);
                    band->x1 = ff_j2k_ceildivpow2(comp->x1, n-1);
                    band->y0 = ff_j2k_ceildivpow2(comp->y0, n-1);
                    band->y1 = ff_j2k_ceildivpow2(comp->y1, n-1);
                }
                else{
                    band->cblkw = 1 << FFMIN(comp->xcb, s->ppx);
                    band->cblkh = 1 << FFMIN(comp->ycb, s->ppy);

                    band->x0 = ff_j2k_ceildivpow2(comp->x0 - (1 << (n-1)) * ((bandno+1)&1), n);
                    band->x1 = ff_j2k_ceildivpow2(comp->x1 - (1 << (n-1)) * ((bandno+1)&1), n);
                    band->y0 = ff_j2k_ceildivpow2(comp->y0 - (1 << (n-1)) * (((bandno+1)&2)>>1), n);
                    band->y1 = ff_j2k_ceildivpow2(comp->y1 - (1 << (n-1)) * (((bandno+1)&2)>>1), n);
                }
                band->cblknx = ff_j2k_ceildiv(band->x1, band->cblkw) - band->x0 / band->cblkw;
                band->cblkny = ff_j2k_ceildiv(band->y1, band->cblkh) - band->y0 / band->cblkh;

                band->cblk = av_malloc(band->cblknx * band->cblkny * sizeof(J2kCblk));
                if (band->cblk == NULL)
                    return -1;
                band->prec = av_malloc(reslevel->nprecw * reslevel->nprech * sizeof(J2kPrec));
                if (band->prec == NULL)
                    return -1;

                for (cblkno = 0; cblkno < band->cblknx * band->cblkny; cblkno++){
                    J2kCblk *cblk = band->cblk + cblkno;
                    cblk->zero = 0;
                    cblk->lblock = 3;
                    cblk->length = 0;
                    cblk->lengthinc = 0;
                    cblk->npassess = 0;
                    cblk->included = 0;
                }

                y0 = band->y0;
                y1 = (band->y0 + (1<<s->ppy))/(1<<s->ppy)*(1<<s->ppy) - band->y0;
                yi0 = 0;
                yi1 = ff_j2k_ceildiv(y1 - y0, 1<<comp->ycb) * (1<<comp->ycb);
                yi1 = FFMIN(yi1, band->cblkny);
                cblkperprech = 1<<(s->ppy - comp->ycb);
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
                xi1 = ff_j2k_ceildiv(x1 - x0, 1<<comp->xcb) * (1<<comp->xcb);
                xi1 = FFMIN(xi1, band->cblknx);

                cblkperprecw = 1<<(s->ppx - comp->xcb);
                for (precx = 0, precno = 0; precx < reslevel->nprecw; precx++){
                    for (precy = 0; precy < reslevel->nprech; precy++, precno = 0){
                        J2kPrec *prec = band->prec + precno;
                        prec->xi0 = xi0;
                        prec->xi1 = xi1;
                        prec->cblkincl = ff_j2k_tag_tree_init(prec->xi1 - prec->xi0,
                                                              prec->yi1 - prec->yi0);
                        prec->zerobits = ff_j2k_tag_tree_init(prec->xi1 - prec->xi0,
                                                              prec->yi1 - prec->yi0);

                    }
                    xi1 += cblkperprecw;
                    xi0 = xi1 - cblkperprecw;
                    xi1 = FFMIN(xi1, band->cblknx);
                }
            }
        }
    }
    return 0;
}

static int getnpassess(J2kDecoderContext *s)
{
    int num;
    if (!get_bits(s, 1))
        return 1;
    if (!get_bits(s, 1))
        return 2;
    if ((num = get_bits(s, 2)) != 3)
        return 3 + num;
    if ((num = get_bits(s, 5)) != 31)
        return 6 + num;
    return 37 + get_bits(s, 7);
}

static int getlblockinc(J2kDecoderContext *s)
{
    int res = 0;
    while (get_bits(s, 1))
        res++;
    return res;
}

static int decode_packet(J2kDecoderContext *s, J2kResLevel *rlevel, int precno, int layno, uint8_t *bandbps)
{
    int bandno, cblkny, cblknx, cblkno;

    if (!get_bits(s, 1)){
        j2k_flush(s);
        return 0;
    }
    for (bandno = 0; bandno < rlevel->nbands; bandno++){
        J2kBand *band = rlevel->band + bandno;
        J2kPrec *prec = band->prec + precno;
        int pos = 0;

        for (cblkny = prec->yi0; cblkny < prec->yi1; cblkny++)
            for(cblknx = prec->xi0, cblkno = cblkny * band->cblknx + cblknx; cblknx < prec->xi1; cblknx++, cblkno++, pos++){
                J2kCblk *cblk = band->cblk + cblkno;
                int incl, newpassess, llen;

                if (cblk->included)
                    incl = get_bits(s, 1);
                else{
                    incl = tag_tree_decode(s, prec->cblkincl + pos, layno+1) == layno;
                }
                if (!incl){
                    continue;
                }
                if (!cblk->included){
                    cblk->included++;
                    cblk->nonzerobits = bandbps[bandno] - tag_tree_decode(s, prec->zerobits + pos, 100);
                }
                newpassess = getnpassess(s);
                llen = getlblockinc(s);
                cblk->lblock += llen;
                cblk->lengthinc = get_bits(s, av_log2(newpassess) + cblk->lblock);
                cblk->npassess += newpassess;
            }
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
                bytestream_get_buffer(&s->buf, cblk->data, cblk->lengthinc);
                cblk->length += cblk->lengthinc;
                cblk->lengthinc = 0;
            }
        }
    }
    return 0;
}

static int decode_packets(J2kDecoderContext *s, J2kTile *tile)
{
    int layno, reslevelno, compno, precno, ok_reslevel = 1;
    s->bit_index = 8;
    for (layno = 0; layno < tile->nlayers; layno++){
        for (reslevelno = 0; ok_reslevel; reslevelno++){
            ok_reslevel = 0;
            for (compno = 0; compno < s->ncomponents; compno++){
                J2kComponent *comp = tile->comp + compno;
                if (reslevelno < comp->nreslevels){
                    J2kResLevel *rlevel = tile->comp[compno].reslevel + reslevelno;
                    ok_reslevel = 1;
                    for (precno = 0; precno < rlevel->nprecw * rlevel->nprech; precno++){
                        decode_packet(s, rlevel, precno, layno, comp->bbps + (reslevelno ? 3*(reslevelno-1)+1 : 0));
                    }
                }
            }
        }
    }
    return 0;
}

/** TIER-1 routines */
static void decode_sigpass(J2kT1Context *t1, int width, int height, int bpno, int bandno)
{
    int mask = 3 << (bpno - 1), i, j, k;

    for (i = 0; i < height; i += 4)
        for (j = 0; j < width; j++)
            for (k = i; k < height && k < i+4; k++){
                if ((t1->flags[k+1][j+1] & J2K_T1_SIG_NB)
                && !(t1->flags[k+1][j+1] & (J2K_T1_SIG | J2K_T1_VIS))){
                    if (ff_aec_decode(&t1->aec, ff_j2k_getnbctxno(t1->flags[k+1][j+1], bandno))){
                        int xorbit, ctxno = ff_j2k_getsgnctxno(t1->flags[k+1][j+1], &xorbit);

                        t1->data[k][j] = (ff_aec_decode(&t1->aec, ctxno) ^ xorbit) ? -mask : mask;

                        ff_j2k_set_significant(t1, j, k);
                    }
                    t1->flags[k+1][j+1] |= J2K_T1_VIS;
                }
            }
}

static void decode_refpass(J2kT1Context *t1, int width, int height, int bpno)
{
    int phalf, nhalf;
    int i, j, k;

    phalf = 1 << (bpno - 1);
    nhalf = -phalf;

    for (i = 0; i < height; i += 4)
        for (j = 0; j < width; j++)
            for (k = i; k < height && k < i+4; k++){
                if ((t1->flags[k+1][j+1] & (J2K_T1_SIG | J2K_T1_VIS)) == J2K_T1_SIG){
                    int ctxno = ff_j2k_getrefctxno(t1->flags[k+1][j+1]);
                    int r = ff_aec_decode(&t1->aec, ctxno) ? phalf : nhalf;
                    t1->data[k][j] += t1->data[k][j] < 0 ? -r : r;
                    t1->flags[k+1][j+1] |= J2K_T1_REF;
                }
            }
}

static void decode_clnpass(J2kT1Context *t1, int width, int height, int bpno, int bandno)
{
    int mask = 3 << (bpno - 1), i, j, k, runlen, dec;

    for (i = 0; i < height; i += 4)
        for (j = 0; j < width; j++){
            if (i + 3 < height && !(
            (t1->flags[i+1][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[i+2][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[i+3][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[i+4][j+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)))){
                if (!ff_aec_decode(&t1->aec, AEC_CX_RL))
                    continue;
                runlen = ff_aec_decode(&t1->aec, AEC_CX_UNI);
                runlen = (runlen << 1) | ff_aec_decode(&t1->aec, AEC_CX_UNI);
                dec = 1;
            }
            else{
                runlen = 0;
                dec = 0;
            }

            for (k = i + runlen; k < i + 4 && k < height; k++){
                if (!dec){
                    if (!(t1->flags[k+1][j+1] & (J2K_T1_SIG | J2K_T1_VIS)))
                        dec = ff_aec_decode(&t1->aec, ff_j2k_getnbctxno(t1->flags[k+1][j+1], bandno));
                }
                if (dec){
                    int xorbit, ctxno = ff_j2k_getsgnctxno(t1->flags[k+1][j+1], &xorbit);
                    t1->data[k][j] = (ff_aec_decode(&t1->aec, ctxno) ^ xorbit) ? -mask : mask;
                    ff_j2k_set_significant(t1, j, k);
                }
                dec = 0;
                t1->flags[k+1][j+1] &= ~J2K_T1_VIS;
            }
        }
}

static int decode_cblk(J2kDecoderContext *s, J2kT1Context *t1, J2kCblk *cblk, int width, int height, int bandpos)
{
    int passno = cblk->npassess, pass_t = 2, bpno = cblk->nonzerobits - 1, i;

    for (i = 0; i < height+2; i++)
        memset(t1->flags[i], 0, (width+2)*sizeof(int));

    for (i = 0; i < height; i++)
        memset(t1->data[i], 0, width*sizeof(int));

    ff_aec_initdec(&t1->aec, cblk->data);
    cblk->data[cblk->length] = 0xff;
    cblk->data[cblk->length+1] = 0xff;

    while(passno--){
        switch(pass_t){
            case 0: decode_sigpass(t1, width, height, bpno+1, bandpos);
                    break;
            case 1: decode_refpass(t1, width, height, bpno+1);
                    break;
            case 2: decode_clnpass(t1, width, height, bpno+1, bandpos);
                    break;
        }

        pass_t++;
        if (pass_t == 3){
            bpno--;
            pass_t = 0;
        }
    }
    return 0;
}

/** inverse discrete wavelet transform routines */
static void sr_1d(int *p, int i0, int i1, int ileft, int iright)
{
#define PSE (i0 + FFMIN((i-i0+2*(i1-i0-1))%(2*(i1-i0-1)), 2*(i1-i0-1)-(i-i0+2*(i1-i0-1))%(2*(i1-i0-1))))
    int i;

    if (i1 == i0 + 1)
        return;

    for (i = i0 - ileft; i < i0; i++)
        p[i] = p[PSE];
    for (i = i1; i < i1+iright; i++)
        p[i] = p[PSE];

    for (i = i0/2; i < i1/2 + 1; i++)
        p[2*i] -= (p[2*i-1] + p[2*i+1] + 2) >> 2;
    for (i = i0/2; i < i1/2; i++)
        p[2*i+1] += (p[2*i] + p[2*i+2]) >> 1;
#undef PSE
}

static int dwt_decode53(J2kDecoderContext *s, J2kComponent *comp, int nreslevels)
{
    int lev = nreslevels, i,
        *t = comp->data, w = comp->x1 - comp->x0;
    int *ppv = av_malloc((comp->reslevel[lev-1].y1 + 4)*sizeof(int)), *pv = ppv+2;
    int *ppu = av_malloc((comp->reslevel[lev-1].x1 + 4)*sizeof(int)), *pu = ppu+2;

    for (i = 1; i < lev; i++){
        int u0 = comp->reslevel[i].x0,
            u1 = comp->reslevel[i].x1,
            v0 = comp->reslevel[i].y0,
            v1 = comp->reslevel[i].y1,
            u = u0, v = v0;
        const static int tileft[2] = {1, 2}, tiright[2] = {2, 1};
        u = u0;
        v = v0;

        /// HOR_SD
        while (v < v1){
            int i, j;
            /// copy with interleaving
            for (i = u0 + (u0 & 1), j = 0; i < u1; i+=2, j++){
                pu[i] = t[w*(v-v0) + j];
            }
            for (i = u0 + 1 - (u0 % 2); i < u1; i+=2, j++){
                pu[i] = t[w*(v-v0) + j];
            }

            sr_1d(pu, u0, u1, tileft[u0&1], tiright[u1&1]);

            for (i = u0; i < u1; i++)
                t[w*(v-v0) + i-u0] = pu[i];

            v++;
        }
        /// VER_SD
        while (u < u1){
            int i, j;
            /// copy with interleaving
            for (i = v0 + (v0 & 1), j = 0; i < v1; i+=2, j++){
                pv[i] = t[w*j + u-u0];
            }
            for (i = v0 + 1 - (v0 % 2); i < v1; i+=2, j++){
                pv[i] = t[w*j + u-u0];
            }

            sr_1d(pv, v0, v1, tileft[v0&1], tiright[v1&1]);

            for (i = v0; i < v1; i++)
                t[w*(i-v0) + u-u0] = pv[i];

            u++;
        }
    }
    av_free(ppv);
    av_free(ppu);
    return 0;
}

static int decode_tile(J2kDecoderContext *s, J2kTile *tile)
{
    int compno, reslevelno, bandno;
    int x, y, *src[4];
    uint8_t *line;
    J2kT1Context t1;

    for (compno = 0; compno < s->ncomponents; compno++){
        J2kComponent *comp = tile->comp + compno;

        for (reslevelno = 0; reslevelno < comp->nreslevels; reslevelno++){
            J2kResLevel *rlevel = comp->reslevel + reslevelno;
            for (bandno = 0; bandno < rlevel->nbands; bandno++){
                J2kBand *band = rlevel->band + bandno;
                int cblkx, cblky, cblkno=0, xx0, x0, xx1, y0, yy0, yy1, bandpos;

                bandpos = bandno + (reslevelno > 0);

                yy0 = bandno == 0 ? 0 : comp->reslevel[reslevelno-1].y1 - comp->reslevel[reslevelno-1].y0;
                y0 = yy0;
                yy1 = FFMIN(ff_j2k_ceildiv(band->y0 + 1, band->cblkh) * band->cblkh, band->y1) - band->y0 + yy0;

                if (band->x0 == band->x1 || band->y0 == band->y1)
                    continue;

                for (cblky = 0; cblky < band->cblkny; cblky++){
                    if (reslevelno == 0 || bandno == 1)
                        xx0 = 0;
                    else
                        xx0 = comp->reslevel[reslevelno-1].x1 - comp->reslevel[reslevelno-1].x0;
                    x0 = xx0;
                    xx1 = FFMIN(ff_j2k_ceildiv(band->x0 + 1, band->cblkw) * band->cblkw, band->x1) - band->x0 + xx0;

                    for (cblkx = 0; cblkx < band->cblknx; cblkx++, cblkno++){
                        int y, x;
                        decode_cblk(s, &t1, band->cblk + cblkno, xx1 - xx0, yy1 - yy0, bandpos);
                        for (y = yy0; y < yy1; y++){
                            int *ptr = t1.data[y-yy0];
                            for (x = xx0; x < xx1; x++){
                                comp->data[(comp->x1 - comp->x0) * y + x] = *ptr++ >> 1;
                            }
                        }
                        xx0 = xx1;
                        xx1 = FFMIN(xx1 + band->cblkw, band->x1 - band->x0 + x0);
                    }
                    yy0 = yy1;
                    yy1 = FFMIN(yy1 + band->cblkh, band->y1 - band->y0 + y0);
                }
            }
        }
        dwt_decode53(s, comp, comp->nreslevels);
        src[compno] = comp->data;
    }

    y = tile->comp[0].y0 - s->Y0siz;

    line = s->picture.data[0] + y * s->picture.linesize[0];

    for (; y < tile->comp[0].y1 - s->Y0siz; y++){
        uint8_t *dst;

        x = tile->comp[0].x0 - s->X0siz;
        dst = line + x * s->ncomponents;

        for (; x < tile->comp[0].x1 - s->X0siz; x++)
            for (compno = 0; compno < s->ncomponents; compno++){
                *src[compno] += 1 << (s->cbps[compno]-1);
                if (*src[compno] < 0)
                    *src[compno] = 0;
                else if (*src[compno] >= (1 << s->cbps[compno]))
                    *src[compno] = (1 << s->cbps[compno]) - 1;
                *dst++ = *src[compno]++;
            }

        line += s->picture.linesize[0];
    }
    return 0;
}

static void cleanup(J2kDecoderContext *s)
{
    int tileno, compno, reslevelno, bandno, precno;
    for (tileno = 0; tileno < s->numXtiles * s->numYtiles; tileno++){
        for (compno = 0; compno < s->ncomponents; compno++){
            J2kComponent *comp = s->tile[tileno].comp + compno;

            for (reslevelno = 0; reslevelno < comp->nreslevels; reslevelno++){
                J2kResLevel *reslevel = comp->reslevel + reslevelno;

                for (bandno = 0; bandno < reslevel->nbands ; bandno++){
                    J2kBand *band = reslevel->band + bandno;

                    for (precno = 0; precno < reslevel->nprecw * reslevel->nprech; precno++){
                        J2kPrec *prec = band->prec + precno;
                        av_free(prec->zerobits);
                        av_free(prec->cblkincl);
                    }
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

static int decode_frame(AVCodecContext *avctx,
                        void *data, int *data_size,
                        uint8_t *buf, int buf_size)
{
    J2kDecoderContext *s = avctx->priv_data;
    AVFrame *picture = data;
    int tileno;

    s->avctx = avctx;
    av_log(s->avctx, AV_LOG_DEBUG, "start\n");

    /// init
    s->buf = s->buf_start = buf;
    s->buf_end = buf + buf_size;
    s->curtileno = -1;

    avcodec_get_frame_defaults((AVFrame*)&s->picture);
    avctx->coded_frame = (AVFrame*)&s->picture;

    s->ppx = s->ppy = 15;

    if (bytestream_get_be16(&s->buf) != J2K_SOC){
        av_log(avctx, AV_LOG_ERROR, "SOC marker not present\n");
        return -1;
    }

    for (;;){
        int marker = bytestream_get_be16(&s->buf), len;
        uint8_t *oldbuf = s->buf;

        if (marker == J2K_SOD){
            J2kTile *tile = s->tile + s->curtileno;
            init_tile(s, s->curtileno);
            decode_packets(s, tile);
            continue;
        }
        if (marker == J2K_EOC)
            break;

        len = bytestream_get_be16(&s->buf);
        switch(marker){
            case J2K_SIZ:
                get_siz(s); break;
            case J2K_COC:
                get_coc(s); break;
            case J2K_COD:
                get_cod(s); break;
            case J2K_QCC:
                get_qcc(s, len); break;
            case J2K_QCD:
                get_qcd(s, len); break;
            case J2K_SOT:
                get_sot(s); break;
            case J2K_COM:
                /// the comment is ignored
                s->buf += len - 2; break;
            default:
                av_log(avctx, AV_LOG_ERROR, "unsupported marker 0x%.4X at pos 0x%x\n", marker, s->buf - s->buf_start - 4);
                return -1;
        }
        if (s->buf - oldbuf != len){
            av_log(avctx, AV_LOG_ERROR, "error during processing marker segment %.4x\n", marker);
            return -1;
        }
    }
    for (tileno = 0; tileno < s->numXtiles * s->numYtiles; tileno++)
        decode_tile(s, s->tile + tileno);

    cleanup(s);
    av_log(s->avctx, AV_LOG_DEBUG, "end\n");

    *data_size = sizeof(AVPicture);
    *picture = s->picture;

    return s->buf - s->buf_start;
}

AVCodec jpeg2000_decoder = {
    "j2k",
    CODEC_TYPE_VIDEO,
    CODEC_ID_JPEG2000,
    sizeof(J2kDecoderContext),
    NULL,
    NULL,
    NULL,
    decode_frame,
    0,
    .pix_fmts =
        (enum PixelFormat[]) {PIX_FMT_GRAY8, PIX_FMT_RGB24, -1}
};
