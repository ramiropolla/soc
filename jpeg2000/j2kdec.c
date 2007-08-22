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

#define SHL(a, n) ((n)>=0 ? (a) << (n) : (a) >> -(n))

typedef struct {
    uint8_t nreslevels;       ///< number of resolution levels
    uint8_t log2_cblk_width,
            log2_cblk_height; ///< exponent of codeblock size
    uint8_t transform;        ///< DWT type
    uint8_t csty;             ///< coding style
    uint8_t log2_prec_width,
            log2_prec_height; ///< precinct size
    uint8_t nlayers;          ///< number of layers
    uint8_t mct;              ///< multiple component transformation
} J2kCodingStyle;

typedef struct {
    uint8_t  expn[32 * 3]; ///< quantization exponent
    uint16_t mant[32 * 3]; ///< quantization mantissa
    uint8_t  quantsty;     ///< quantization style
    uint8_t  nguardbits;   ///< number of guard bits
} J2kQuantStyle;

typedef struct {
   J2kResLevel *reslevel;
   int *data;
   uint16_t x0, x1, y0, y1;
} J2kComponent;

#define HAD_COC 0x01
#define HAD_QCC 0x02

typedef struct {
   J2kComponent *comp;
   uint8_t properties[4];
   J2kCodingStyle codsty[4];
   J2kQuantStyle  qntsty[4];
} J2kTile;

typedef struct {
    AVCodecContext *avctx;
    AVFrame picture;

    int width, height; ///< image width and height
    int image_offset_x, image_offset_y;
    int tile_offset_x, tile_offset_y;
    uint8_t cbps[4]; ///< numbps in components
    uint8_t sgnd[4]; ///< if a component is signed
    uint8_t properties[4];

    int ncomponents;
    int tile_width, tile_height; ///< tile size
    int numXtiles, numYtiles;
    int maxtilelen;

    J2kCodingStyle codsty[4];
    J2kQuantStyle  qntsty[4];

    uint8_t *buf_start;
    uint8_t *buf;
    uint8_t *buf_end;
    int bit_index;

    int16_t curtileno;

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
    if (*s->buf == 0xff)
        s->buf++;
    s->bit_index = 8;
    s->buf++;
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
            for(reslevelno = 0; reslevelno < codsty->nreslevels; reslevelno++){
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

/** decode the value stored in node */
static int tag_tree_decode(J2kDecoderContext *s, J2kTgtNode *node, int threshold)
{
    J2kTgtNode *stack[30];
    int sp = -1, curval = 0;

    while(node && !node->vis){
        stack[++sp] = node;
        node = node->parent;
    }

    if (node)
        curval = node->val;
    else
        curval = stack[sp]->val;

    while(curval < threshold && sp >= 0){
        if (curval < stack[sp]->val)
            curval = stack[sp]->val;
        while (curval < threshold){
            if (get_bits(s, 1)){
                stack[sp]->vis++;
                break;
            } else
                curval++;
        }
        stack[sp]->val = curval;
        sp--;
    }
    return curval;
}

/* marker segments */
/** get sizes and offsets of image, tiles; number of components */
static int get_siz(J2kDecoderContext *s)
{
    int i;

                        bytestream_get_be16(&s->buf); // Rsiz (skipped)
             s->width = bytestream_get_be32(&s->buf); // width
            s->height = bytestream_get_be32(&s->buf); // height
    s->image_offset_x = bytestream_get_be32(&s->buf); // X0Siz
    s->image_offset_y = bytestream_get_be32(&s->buf); // Y0Siz

        s->tile_width = bytestream_get_be32(&s->buf); // XTSiz
       s->tile_height = bytestream_get_be32(&s->buf); // YTSiz
     s->tile_offset_x = bytestream_get_be32(&s->buf); // XT0Siz
     s->tile_offset_y = bytestream_get_be32(&s->buf); // YT0Siz
       s->ncomponents = bytestream_get_be16(&s->buf); // CSiz

    for (i = 0; i < s->ncomponents; i++){ // Ssiz_i XRsiz_i, YRsiz_i
        uint8_t x = bytestream_get_byte(&s->buf);
        s->cbps[i] = (x & 0x7f) + 1;
        s->sgnd[i] = (x & 0x80) == 1;
        if (bytestream_get_byte(&s->buf) != 1)
            return -1;
        if (bytestream_get_byte(&s->buf) != 1)
            return -1;
    }

    s->numXtiles = ff_j2k_ceildiv(s->width - s->tile_offset_x, s->tile_width);
    s->numYtiles = ff_j2k_ceildiv(s->height - s->tile_offset_y, s->tile_height);

    s->tile = av_mallocz(s->numXtiles * s->numYtiles * sizeof(J2kTile));
    if (!s->tile)
        return AVERROR(ENOMEM);

    for (i = 0; i < s->numXtiles * s->numYtiles; i++){
        J2kTile *tile = s->tile + i;

        tile->comp = av_mallocz(s->ncomponents * sizeof(J2kComponent));
        if (!tile->comp)
            return AVERROR(ENOMEM);
    }

    s->avctx->width = s->width - s->image_offset_x;
    s->avctx->height = s->height - s->image_offset_y;

    switch(s->ncomponents){
        case 1: s->avctx->pix_fmt = PIX_FMT_GRAY8; break;
        case 3: s->avctx->pix_fmt = PIX_FMT_RGB24; break;
        case 4: s->avctx->pix_fmt = PIX_FMT_BGRA; break;
    }

    if (s->picture.data[0])
        s->avctx->release_buffer(s->avctx, &s->picture);

    s->avctx->get_buffer(s->avctx, &s->picture);

    s->picture.pict_type = FF_I_TYPE;
    s->picture.key_frame = 1;

    return 0;
}

/** get common part for COD and COC segments */
static int get_cox(J2kDecoderContext *s, J2kCodingStyle *c)
{
          c->nreslevels = bytestream_get_byte(&s->buf) + 1; // num of resolution levels - 1
     c->log2_cblk_width = bytestream_get_byte(&s->buf) + 2; // cblk width
    c->log2_cblk_height = bytestream_get_byte(&s->buf) + 2; // cblk height

    if (bytestream_get_byte(&s->buf) != 0){ // cblk style
        av_log(s->avctx, AV_LOG_ERROR, "no extra cblk styles supported\n");
        return -1;
    }
    c->transform = bytestream_get_byte(&s->buf); // transformation
    return 0;
}

/** get coding parameters for a particular tile or whole image*/
static int get_cod(J2kDecoderContext *s, J2kCodingStyle *c, uint8_t *properties)
{
    J2kCodingStyle tmp;
    int compno;

    tmp.log2_prec_width  =
    tmp.log2_prec_height = 15;

    tmp.csty = bytestream_get_byte(&s->buf);

    if (bytestream_get_byte(&s->buf)){ // progression level
        av_log(s->avctx, AV_LOG_ERROR, "only LRCP progression supported\n");
        return -1;
    }

    tmp.nlayers = bytestream_get_be16(&s->buf);
        tmp.mct = bytestream_get_byte(&s->buf); // multiple component transformation

    get_cox(s, &tmp);
    for (compno = 0; compno < s->ncomponents; compno++){
        if (!(properties[compno] & HAD_COC))
            memcpy(c + compno, &tmp, sizeof(J2kCodingStyle));
    }
    return 0;
}

/** get coding parameters for a component in the whole image on a particular tile */
static int get_coc(J2kDecoderContext *s, J2kCodingStyle *c, uint8_t *properties)
{
    int compno = bytestream_get_byte(&s->buf);

    c += compno;
    c->csty = bytestream_get_byte(&s->buf);
    get_cox(s, c);

    properties[compno] |= HAD_COC;
    return 0;
}

/** get common part for QCD and QCC segments */
static int get_qcx(J2kDecoderContext *s, int n, J2kQuantStyle *q)
{
    int i, x;
    x = bytestream_get_byte(&s->buf); // Sqcd

    q->nguardbits = x >> 5;
      q->quantsty = x & 0x1f;

    if (q->quantsty == J2K_QSTY_NONE){
        n -= 3;
        for (i = 0; i < n; i++)
            q->expn[i] = bytestream_get_byte(&s->buf) >> 3;
    }
    else if (q->quantsty == J2K_QSTY_SI){
        x = bytestream_get_be16(&s->buf);
        q->expn[0] = x >> 11;
        q->mant[0] = x & 0x7ff;
        for (i = 1; i < 32 * 3; i++){
            int curexpn = FFMAX(0, q->expn[0] - (i-1)/3);
            q->expn[i] = curexpn;
            q->mant[i] = q->mant[0];
        }
    }
    else{
        n = (n - 3) >> 1;
        for (i = 0; i < n; i++){
            x = bytestream_get_be16(&s->buf);
            q->expn[i] = x >> 11;
            q->mant[i] = x & 0x7ff;
        }
    }
    return 0;
}

/** get quantization parameters for a particular tile or a whole image */
static int get_qcd(J2kDecoderContext *s, int n, J2kQuantStyle *q, uint8_t *properties)
{
    J2kQuantStyle tmp;
    int compno;

    if (get_qcx(s, n, &tmp))
        return -1;
    for (compno = 0; compno < s->ncomponents; compno++)
        if (!(properties[compno] & HAD_QCC))
            memcpy(q + compno, &tmp, sizeof(J2kQuantStyle));
    return 0;
}

/** get quantization paramteres for a component in the whole image on in a particular tile */
static int get_qcc(J2kDecoderContext *s, int n, J2kQuantStyle *q, uint8_t *properties)
{
    int compno = bytestream_get_byte(&s->buf);
    properties[compno] |= HAD_QCC;
    return get_qcx(s, n-1, q+compno);
}

/** get start of tile segment */
static uint8_t get_sot(J2kDecoderContext *s)
{
    s->curtileno = bytestream_get_be16(&s->buf); ///< Isot

    s->buf += 4; ///< Psot (ignored)

    if (!bytestream_get_byte(&s->buf)){ ///< TPsot
        J2kTile *tile = s->tile + s->curtileno;

        /* copy defaults */
        memcpy(tile->codsty, s->codsty, s->ncomponents * sizeof(J2kCodingStyle));
        memcpy(tile->qntsty, s->qntsty, s->ncomponents * sizeof(J2kQuantStyle));
    }
    bytestream_get_byte(&s->buf); ///< TNsot

    return 0;
}

static int init_tile(J2kDecoderContext *s, int tileno)
{
    int compno, reslevelno, bandno, p, q;
    J2kTile *tile = s->tile + tileno;

    p = tileno % s->numXtiles;
    q = tileno / s->numXtiles;

    if (!tile->comp)
        return AVERROR(ENOMEM);
    for (compno = 0; compno < s->ncomponents; compno++){
        J2kComponent *comp = tile->comp + compno;
        J2kCodingStyle *codsty = tile->codsty + compno;
        J2kQuantStyle  *qntsty = tile->qntsty + compno;
        int gbandno = 0; // global bandno

        comp->x0 = FFMAX(p * s->tile_width + s->tile_offset_x, s->image_offset_x);
        comp->x1 = FFMIN((p+1)*s->tile_width + s->tile_offset_x, s->width);
        comp->y0 = FFMAX(q * s->tile_height + s->tile_offset_y, s->image_offset_y);
        comp->y1 = FFMIN((q+1)*s->tile_height + s->tile_offset_y, s->height);

        comp->data = av_malloc((comp->y1 - comp->y0) * (comp->x1 -comp->x0) * sizeof(int));
        if (!comp->data)
            return AVERROR(ENOMEM);
        comp->reslevel = av_malloc(codsty->nreslevels * sizeof(J2kResLevel));

        if (!comp->reslevel)
            return AVERROR(ENOMEM);
        for (reslevelno = 0; reslevelno < codsty->nreslevels; reslevelno++){
            int n = codsty->nreslevels - reslevelno;
            J2kResLevel *reslevel = comp->reslevel + reslevelno;

            reslevel->x0 = ff_j2k_ceildivpow2(comp->x0, codsty->nreslevels - reslevelno - 1);
            reslevel->x1 = ff_j2k_ceildivpow2(comp->x1, codsty->nreslevels - reslevelno - 1);
            reslevel->y0 = ff_j2k_ceildivpow2(comp->y0, codsty->nreslevels - reslevelno - 1);
            reslevel->y1 = ff_j2k_ceildivpow2(comp->y1, codsty->nreslevels - reslevelno - 1);

            if (reslevelno == 0)
                reslevel->nbands = 1;
            else
                reslevel->nbands = 3;

            if (reslevel->x1 == reslevel->x0)
                reslevel->num_precincts_x = 0;
            else
                reslevel->num_precincts_x = ff_j2k_ceildivpow2(reslevel->x1, codsty->log2_prec_width) - reslevel->x0 / (1<<codsty->log2_prec_width);

            if (reslevel->y1 == reslevel->y0)
                reslevel->num_precincts_y = 0;
            else
                reslevel->num_precincts_y = ff_j2k_ceildivpow2(reslevel->y1, codsty->log2_prec_height) - reslevel->y0 / (1<<codsty->log2_prec_height);

            reslevel->band = av_malloc(reslevel->nbands * sizeof(J2kBand));
            if (!reslevel->band)
                return AVERROR(ENOMEM);
            for (bandno = 0; bandno < reslevel->nbands; bandno++, gbandno++){
                J2kBand *band = reslevel->band + bandno;
                int cblkno, precx, precy, precno;
                int x0, y0, x1, y1;
                int xi0, yi0, xi1, yi1;
                int cblkperprecw, cblkperprech;

                if (qntsty->quantsty != J2K_QSTY_NONE){
                    const static uint8_t lut_gain[2][4] = {{0, 0, 0, 0}, {0, 1, 1, 2}};
                    int numbps;

                    numbps = s->cbps[compno] + lut_gain[codsty->transform][bandno + reslevelno>0];
                    band->stepsize = SHL(2048 + qntsty->mant[gbandno], 2 + numbps - qntsty->expn[gbandno]);
                }
                else
                    band->stepsize = 1 << 13;

                if (reslevelno == 0){  // the same everywhere
                    band->codeblock_width = 1 << FFMIN(codsty->log2_cblk_width, codsty->log2_prec_width-1);
                    band->codeblock_height = 1 << FFMIN(codsty->log2_cblk_height, codsty->log2_prec_height-1);

                    band->x0 = ff_j2k_ceildivpow2(comp->x0, n-1);
                    band->x1 = ff_j2k_ceildivpow2(comp->x1, n-1);
                    band->y0 = ff_j2k_ceildivpow2(comp->y0, n-1);
                    band->y1 = ff_j2k_ceildivpow2(comp->y1, n-1);
                }
                else{
                    band->codeblock_width = 1 << FFMIN(codsty->log2_cblk_width, codsty->log2_prec_width);
                    band->codeblock_height = 1 << FFMIN(codsty->log2_cblk_height, codsty->log2_prec_height);

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
                    J2kCblk *cblk = band->cblk + cblkno;
                    cblk->zero = 0;
                    cblk->lblock = 3;
                    cblk->length = 0;
                    cblk->lengthinc = 0;
                    cblk->npasses = 0;
                }

                y0 = band->y0;
                y1 = (band->y0 + (1<<codsty->log2_prec_height))/(1<<codsty->log2_prec_height)*(1<<codsty->log2_prec_height) - band->y0;
                yi0 = 0;
                yi1 = ff_j2k_ceildiv(y1 - y0, 1<<codsty->log2_cblk_height) * (1<<codsty->log2_cblk_height);
                yi1 = FFMIN(yi1, band->cblkny);
                cblkperprech = 1<<(codsty->log2_prec_height - codsty->log2_cblk_height);
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
                x1 = (band->x0 + (1<<codsty->log2_prec_width))/(1<<codsty->log2_prec_width)*(1<<codsty->log2_prec_width) - band->x0;
                xi0 = 0;
                xi1 = ff_j2k_ceildiv(x1 - x0, 1<<codsty->log2_cblk_width) * (1<<codsty->log2_cblk_width);
                xi1 = FFMIN(xi1, band->cblknx);

                cblkperprecw = 1<<(codsty->log2_prec_width - codsty->log2_cblk_width);
                for (precx = 0, precno = 0; precx < reslevel->num_precincts_x; precx++){
                    for (precy = 0; precy < reslevel->num_precincts_y; precy++, precno = 0){
                        J2kPrec *prec = band->prec + precno;
                        prec->xi0 = xi0;
                        prec->xi1 = xi1;
                        prec->cblkincl = ff_j2k_tag_tree_init(prec->xi1 - prec->xi0,
                                                              prec->yi1 - prec->yi0);
                        prec->zerobits = ff_j2k_tag_tree_init(prec->xi1 - prec->xi0,
                                                              prec->yi1 - prec->yi0);
                        if (!prec->cblkincl || !prec->zerobits)
                            return AVERROR(ENOMEM);

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

/** read the number of coding passes */
static int getnpasses(J2kDecoderContext *s)
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

static int decode_packet(J2kDecoderContext *s, J2kResLevel *rlevel, int precno, int layno, uint8_t *expn, int numgbits)
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
                int incl, newpasses, llen;

                if (cblk->npasses)
                    incl = get_bits(s, 1);
                else{
                    incl = tag_tree_decode(s, prec->cblkincl + pos, layno+1) == layno;
                }
                if (!incl){
                    continue;
                }
                if (!cblk->npasses)
                    cblk->nonzerobits = expn[bandno] + numgbits - 1 - tag_tree_decode(s, prec->zerobits + pos, 100);
                newpasses = getnpasses(s);
                llen = getlblockinc(s);
                cblk->lblock += llen;
                cblk->lengthinc = get_bits(s, av_log2(newpasses) + cblk->lblock);
                cblk->npasses += newpasses;
            }
    }
    j2k_flush(s);
    for (bandno = 0; bandno < rlevel->nbands; bandno++){
        J2kBand *band = rlevel->band + bandno;
        int yi, cblknw = band->prec[precno].xi1 - band->prec[precno].xi0;
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
    int layno, reslevelno, compno, precno, ok_reslevel;
    s->bit_index = 8;
    for (layno = 0; layno < tile->codsty[0].nlayers; layno++){
        ok_reslevel = 1;
        for (reslevelno = 0; ok_reslevel; reslevelno++){
            ok_reslevel = 0;
            for (compno = 0; compno < s->ncomponents; compno++){
                J2kCodingStyle *codsty = tile->codsty + compno;
                J2kQuantStyle  *qntsty = tile->qntsty + compno;
                if (reslevelno < codsty->nreslevels){
                    J2kResLevel *rlevel = tile->comp[compno].reslevel + reslevelno;
                    ok_reslevel = 1;
                    for (precno = 0; precno < rlevel->num_precincts_x * rlevel->num_precincts_y; precno++){
                        if (decode_packet(s, rlevel, precno, layno, qntsty->expn + (reslevelno ? 3*(reslevelno-1)+1 : 0),
                                          qntsty->nguardbits))
                            return -1;
                    }
                }
            }
        }
    }
    return 0;
}

/* TIER-1 routines */
static void decode_sigpass(J2kT1Context *t1, int width, int height, int bpno, int bandno)
{
    int mask = 3 << (bpno - 1), y0, x, y;

    for (y0 = 0; y0 < height; y0 += 4)
        for (x = 0; x < width; x++)
            for (y = y0; y < height && y < y0+4; y++){
                if ((t1->flags[y+1][x+1] & J2K_T1_SIG_NB)
                && !(t1->flags[y+1][x+1] & (J2K_T1_SIG | J2K_T1_VIS))){
                    if (ff_aec_decode(&t1->aec, t1->aec.cx_states + ff_j2k_getnbctxno(t1->flags[y+1][x+1], bandno))){
                        int xorbit, ctxno = ff_j2k_getsgnctxno(t1->flags[y+1][x+1], &xorbit);

                        t1->data[y][x] = (ff_aec_decode(&t1->aec, t1->aec.cx_states + ctxno) ^ xorbit) ? -mask : mask;

                        ff_j2k_set_significant(t1, x, y);
                    }
                    t1->flags[y+1][x+1] |= J2K_T1_VIS;
                }
            }
}

static void decode_refpass(J2kT1Context *t1, int width, int height, int bpno)
{
    int phalf, nhalf;
    int y0, x, y;

    phalf = 1 << (bpno - 1);
    nhalf = -phalf;

    for (y0 = 0; y0 < height; y0 += 4)
        for (x = 0; x < width; x++)
            for (y = y0; y < height && y < y0+4; y++){
                if ((t1->flags[y+1][x+1] & (J2K_T1_SIG | J2K_T1_VIS)) == J2K_T1_SIG){
                    int ctxno = ff_j2k_getrefctxno(t1->flags[y+1][x+1]);
                    int r = ff_aec_decode(&t1->aec, t1->aec.cx_states + ctxno) ? phalf : nhalf;
                    t1->data[y][x] += t1->data[y][x] < 0 ? -r : r;
                    t1->flags[y+1][x+1] |= J2K_T1_REF;
                }
            }
}

static void decode_clnpass(J2kT1Context *t1, int width, int height, int bpno, int bandno)
{
    int mask = 3 << (bpno - 1), y0, x, y, runlen, dec;

    for (y0 = 0; y0 < height; y0 += 4)
        for (x = 0; x < width; x++){
            if (y0 + 3 < height && !(
            (t1->flags[y0+1][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[y0+2][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[y0+3][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)) ||
            (t1->flags[y0+4][x+1] & (J2K_T1_SIG_NB | J2K_T1_VIS | J2K_T1_SIG)))){
                if (!ff_aec_decode(&t1->aec, t1->aec.cx_states + AEC_CX_RL))
                    continue;
                runlen = ff_aec_decode(&t1->aec, t1->aec.cx_states + AEC_CX_UNI);
                runlen = (runlen << 1) | ff_aec_decode(&t1->aec, t1->aec.cx_states + AEC_CX_UNI);
                dec = 1;
            }
            else{
                runlen = 0;
                dec = 0;
            }

            for (y = y0 + runlen; y < y0 + 4 && y < height; y++){
                if (!dec){
                    if (!(t1->flags[y+1][x+1] & (J2K_T1_SIG | J2K_T1_VIS)))
                        dec = ff_aec_decode(&t1->aec, t1->aec.cx_states + ff_j2k_getnbctxno(t1->flags[y+1][x+1], bandno));
                }
                if (dec){
                    int xorbit, ctxno = ff_j2k_getsgnctxno(t1->flags[y+1][x+1], &xorbit);
                    t1->data[y][x] = (ff_aec_decode(&t1->aec, t1->aec.cx_states + ctxno) ^ xorbit) ? -mask : mask;
                    ff_j2k_set_significant(t1, x, y);
                }
                dec = 0;
                t1->flags[y+1][x+1] &= ~J2K_T1_VIS;
            }
        }
}

static int decode_cblk(J2kDecoderContext *s, J2kT1Context *t1, J2kCblk *cblk, int width, int height, int bandpos)
{
    int passno = cblk->npasses, pass_t = 2, bpno = cblk->nonzerobits - 1, y;

    for (y = 0; y < height+2; y++)
        memset(t1->flags[y], 0, (width+2)*sizeof(int));

    for (y = 0; y < height; y++)
        memset(t1->data[y], 0, width*sizeof(int));

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

/* inverse discrete wavelet transform routines */
static void sr_1d53(int *p, int i0, int i1)
{
    int i;

    if (i1 == i0 + 1)
        return;

    p[i0 - 1] = p[i0 + 1];
    p[i1    ] = p[i1 - 2];
    p[i0 - 2] = p[i0 + 2];
    p[i1 + 1] = p[i1 - 3];

    for (i = i0/2; i < i1/2 + 1; i++)
        p[2*i] -= (p[2*i-1] + p[2*i+1] + 2) >> 2;
    for (i = i0/2; i < i1/2; i++)
        p[2*i+1] += (p[2*i] + p[2*i+2]) >> 1;
}

static void sr_1d97(float *p, int i0, int i1)
{
    int i;

    if (i1 == i0 + 1)
        return;

    for (i = 1; i <= 4; i++){
        p[i0 - i] = p[i0 + i];
        p[i1 + i - 1] = p[i1 - i - 1];
    }

    for (i = i0/2 - 1; i < i1/2 + 2; i++)
        p[2*i] *= 1.230174;
    for (i = i0/2 - 2; i < i1/2 + 2; i++)
        p[2*i+1] *= 1.625786;
    for (i = i0/2 - 1; i < i1/2 + 2; i++)
        p[2*i] -= 0.443506 * (p[2*i-1] + p[2*i+1]);
    for (i = i0/2 - 1; i < i1/2 + 1; i++)
        p[2*i+1] -= 0.882911 * (p[2*i] + p[2*i+2]);
    for (i = i0/2; i < i1/2 + 1; i++)
        p[2*i] += 0.052980 * (p[2*i-1] + p[2*i+1]);
    for (i = i0/2; i < i1/2; i++)
        p[2*i+1] += 1.586134 * (p[2*i] + p[2*i+2]);
}

static int dwt_decode53(J2kDecoderContext *s, J2kComponent *comp, int nreslevels)
{
    int lev = nreslevels, i,
        *t = comp->data, w = comp->x1 - comp->x0;
    int *ppv = av_malloc((comp->reslevel[lev-1].y1 + 4)*sizeof(int)), *pv = ppv+2;
    int *ppu = av_malloc((comp->reslevel[lev-1].x1 + 4)*sizeof(int)), *pu = ppu+2;

    if (!ppv || !ppu)
        return -1;

    for (i = 1; i < lev; i++){
        int u0 = comp->reslevel[i].x0,
            u1 = comp->reslevel[i].x1,
            v0 = comp->reslevel[i].y0,
            v1 = comp->reslevel[i].y1,
            u = u0, v = v0;

        // HOR_SD
        while (v < v1){
            int i, j;
            // copy with interleaving
            for (i = u0 + (u0 & 1), j = 0; i < u1; i+=2, j++){
                pu[i] = t[w*(v-v0) + j];
            }
            for (i = u0 + 1 - (u0 % 2); i < u1; i+=2, j++){
                pu[i] = t[w*(v-v0) + j];
            }

            sr_1d53(pu, u0, u1);

            for (i = u0; i < u1; i++)
                t[w*(v-v0) + i-u0] = pu[i];

            v++;
        }
        // VER_SD
        while (u < u1){
            int i, j;
            // copy with interleaving
            for (i = v0 + (v0 & 1), j = 0; i < v1; i+=2, j++){
                pv[i] = t[w*j + u-u0];
            }
            for (i = v0 + 1 - (v0 % 2); i < v1; i+=2, j++){
                pv[i] = t[w*j + u-u0];
            }

            sr_1d53(pv, v0, v1);

            for (i = v0; i < v1; i++)
                t[w*(i-v0) + u-u0] = pv[i];

            u++;
        }
    }
    av_free(ppv);
    av_free(ppu);
    return 0;
}

static int dwt_decode97(J2kDecoderContext *s, J2kComponent *comp, int nreslevels)
{
    int lev = nreslevels, i,
        *t = comp->data, w = comp->x1 - comp->x0;
    float *ppv = av_malloc((comp->reslevel[lev-1].y1 + 8)*sizeof(float)), *pv = ppv+4;
    float *ppu = av_malloc((comp->reslevel[lev-1].x1 + 8)*sizeof(float)), *pu = ppu+4;

    if (!ppv || !ppu)
        return -1;

    for (i = 1; i < lev; i++){
        int u0 = comp->reslevel[i].x0,
            u1 = comp->reslevel[i].x1,
            v0 = comp->reslevel[i].y0,
            v1 = comp->reslevel[i].y1,
            u = u0, v = v0;

        // HOR_SD
        while (v < v1){
            int i, j;
            // copy with interleaving
            for (i = u0 + (u0 & 1), j = 0; i < u1; i+=2, j++){
                pu[i] = t[w*(v-v0) + j];
            }
            for (i = u0 + 1 - (u0 % 2); i < u1; i+=2, j++){
                pu[i] = t[w*(v-v0) + j];
            }

            sr_1d97(pu, u0, u1);

            for (i = u0; i < u1; i++)
                t[w*(v-v0) + i-u0] = pu[i];

            v++;
        }
        // VER_SD
        while (u < u1){
            int i, j;
            // copy with interleaving
            for (i = v0 + (v0 & 1), j = 0; i < v1; i+=2, j++){
                pv[i] = t[w*j + u-u0];
            }
            for (i = v0 + 1 - (v0 % 2); i < v1; i+=2, j++){
                pv[i] = t[w*j + u-u0];
            }

            sr_1d97(pv, v0, v1);

            for (i = v0; i < v1; i++)
                t[w*(i-v0) + u-u0] = pv[i];

            u++;
        }
    }
    av_free(ppv);
    av_free(ppu);
    return 0;
}

static void mct_decode(J2kDecoderContext *s, J2kTile *tile)
{
    int i, *src[3], i0, i1, i2;

    for (i = 0; i < 3; i++)
        src[i] = tile->comp[i].data;

    if (tile->codsty[0].transform == J2K_DWT97){
        for (i = 0; i < (tile->comp[0].y1 - tile->comp[0].y0) * (tile->comp[0].x1 - tile->comp[0].x0); i++){
            i0 = *src[0] + (*src[2] * 46802 >> 16);
            i1 = *src[0] - (*src[1] * 22553 + *src[2] * 46802 >> 16);
            i2 = *src[0] + (116130 * *src[1] >> 16);
            *src[0]++ = i0;
            *src[1]++ = i1;
            *src[2]++ = i2;
        }
    } else{
        for (i = 0; i < (tile->comp[0].y1 - tile->comp[0].y0) * (tile->comp[0].x1 - tile->comp[0].x0); i++){
            i1 = *src[0] - (*src[2] + *src[1] >> 2);
            i0 = i1 + *src[2];
            i2 = i1 + *src[1];
            *src[0]++ = i0;
            *src[1]++ = i1;
            *src[2]++ = i2;
        }
    }
}

static int decode_tile(J2kDecoderContext *s, J2kTile *tile)
{
    int compno, reslevelno, bandno;
    int x, y, *src[4];
    uint8_t *line;
    J2kT1Context t1;

    for (compno = 0; compno < s->ncomponents; compno++){
        J2kComponent *comp = tile->comp + compno;
        J2kCodingStyle *codsty = tile->codsty + compno;

        for (reslevelno = 0; reslevelno < codsty->nreslevels; reslevelno++){
            J2kResLevel *rlevel = comp->reslevel + reslevelno;
            for (bandno = 0; bandno < rlevel->nbands; bandno++){
                J2kBand *band = rlevel->band + bandno;
                int cblkx, cblky, cblkno=0, xx0, x0, xx1, y0, yy0, yy1, bandpos;

                bandpos = bandno + (reslevelno > 0);

                yy0 = bandno == 0 ? 0 : comp->reslevel[reslevelno-1].y1 - comp->reslevel[reslevelno-1].y0;
                y0 = yy0;
                yy1 = FFMIN(ff_j2k_ceildiv(band->y0 + 1, band->codeblock_height) * band->codeblock_height, band->y1) - band->y0 + yy0;

                if (band->x0 == band->x1 || band->y0 == band->y1)
                    continue;

                for (cblky = 0; cblky < band->cblkny; cblky++){
                    if (reslevelno == 0 || bandno == 1)
                        xx0 = 0;
                    else
                        xx0 = comp->reslevel[reslevelno-1].x1 - comp->reslevel[reslevelno-1].x0;
                    x0 = xx0;
                    xx1 = FFMIN(ff_j2k_ceildiv(band->x0 + 1, band->codeblock_width) * band->codeblock_width, band->x1) - band->x0 + xx0;

                    for (cblkx = 0; cblkx < band->cblknx; cblkx++, cblkno++){
                        int y, x;
                        decode_cblk(s, &t1, band->cblk + cblkno, xx1 - xx0, yy1 - yy0, bandpos);
                        if (codsty->transform == J2K_DWT53){
                            for (y = yy0; y < yy1; y++){
                                int *ptr = t1.data[y-yy0];
                                for (x = xx0; x < xx1; x++){
                                    comp->data[(comp->x1 - comp->x0) * y + x] = *ptr++ >> 1;
                                }
                            }
                        } else{
                            for (y = yy0; y < yy1; y++){
                                int *ptr = t1.data[y-yy0];
                                for (x = xx0; x < xx1; x++){
                                    int tmp = ((int64_t)*ptr++) * ((int64_t)band->stepsize) >> 13, tmp2;
                                    tmp2 = FFABS(tmp>>1) + FFABS(tmp&1);
                                    comp->data[(comp->x1 - comp->x0) * y + x] = tmp < 0 ? -tmp2 : tmp2;
                                }
                            }
                        }
                        xx0 = xx1;
                        xx1 = FFMIN(xx1 + band->codeblock_width, band->x1 - band->x0 + x0);
                    }
                    yy0 = yy1;
                    yy1 = FFMIN(yy1 + band->codeblock_height, band->y1 - band->y0 + y0);
                }
            }
        }
        if (codsty->transform == J2K_DWT53){
            if (dwt_decode53(s, comp, codsty->nreslevels))
                return -1;
        } else{
            if (dwt_decode97(s, comp, codsty->nreslevels))
                return -1;
        }
        src[compno] = comp->data;
    }
    if (tile->codsty[0].mct)
        mct_decode(s, tile);

    y = tile->comp[0].y0 - s->image_offset_y;

    line = s->picture.data[0] + y * s->picture.linesize[0];
    if (s->avctx->pix_fmt == PIX_FMT_BGRA) // RGBA -> BGRA
        FFSWAP(int *, src[0], src[2]);

    for (; y < tile->comp[0].y1 - s->image_offset_y; y++){
        uint8_t *dst;

        x = tile->comp[0].x0 - s->image_offset_x;
        dst = line + x * s->ncomponents;

        for (; x < tile->comp[0].x1 - s->image_offset_x; x++)
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
            J2kCodingStyle *codsty = s->tile[tileno].codsty + compno;

            for (reslevelno = 0; reslevelno < codsty->nreslevels; reslevelno++){
                J2kResLevel *reslevel = comp->reslevel + reslevelno;

                for (bandno = 0; bandno < reslevel->nbands ; bandno++){
                    J2kBand *band = reslevel->band + bandno;

                    for (precno = 0; precno < reslevel->num_precincts_x * reslevel->num_precincts_y; precno++){
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

static int decode_codestream(J2kDecoderContext *s)
{
    J2kCodingStyle *codsty = s->codsty;
    J2kQuantStyle  *qntsty = s->qntsty;
    uint8_t *properties = s->properties;

    for (;;){
        int marker = bytestream_get_be16(&s->buf), len, ret = 0;
        uint8_t *oldbuf = s->buf;

        if (marker == J2K_SOD){
            J2kTile *tile = s->tile + s->curtileno;
            if (ret = init_tile(s, s->curtileno))
                return ret;
            if (ret = decode_packets(s, tile))
                return ret;
            continue;
        }
        if (marker == J2K_EOC)
            break;

        len = bytestream_get_be16(&s->buf);
        switch(marker){
            case J2K_SIZ:
                ret = get_siz(s); break;
            case J2K_COC:
                ret = get_coc(s, codsty, properties); break;
            case J2K_COD:
                ret = get_cod(s, codsty, properties); break;
            case J2K_QCC:
                ret = get_qcc(s, len, qntsty, properties); break;
            case J2K_QCD:
                ret = get_qcd(s, len, qntsty, properties); break;
            case J2K_SOT:
                if (!(ret = get_sot(s))){
                    codsty = s->tile[s->curtileno].codsty;
                    qntsty = s->tile[s->curtileno].qntsty;
                    properties = s->tile[s->curtileno].properties;
                }
                break;
            case J2K_COM:
                // the comment is ignored
                s->buf += len - 2; break;
            default:
                av_log(s->avctx, AV_LOG_ERROR, "unsupported marker 0x%.4X at pos 0x%x\n", marker, s->buf - s->buf_start - 4);
                return -1;
        }
        if (s->buf - oldbuf != len || ret){
            av_log(s->avctx, AV_LOG_ERROR, "error during processing marker segment %.4x\n", marker);
            return ret ? ret : -1;
        }
    }
    return 0;
}

static int decode_frame(AVCodecContext *avctx,
                        void *data, int *data_size,
                        uint8_t *buf, int buf_size)
{
    J2kDecoderContext *s = avctx->priv_data;
    AVFrame *picture = data;
    int tileno, ret;

    s->avctx = avctx;
    av_log(s->avctx, AV_LOG_DEBUG, "start\n");

    // init
    s->buf = s->buf_start = buf;
    s->buf_end = buf + buf_size;
    s->curtileno = -1;

    ff_j2k_init_tier1_luts();

    if (bytestream_get_be16(&s->buf) != J2K_SOC){
        av_log(avctx, AV_LOG_ERROR, "SOC marker not present\n");
        return -1;
    }
    if (ret = decode_codestream(s))
        return ret;

    for (tileno = 0; tileno < s->numXtiles * s->numYtiles; tileno++)
        if (ret = decode_tile(s, s->tile + tileno))
            return ret;

    cleanup(s);
    av_log(s->avctx, AV_LOG_DEBUG, "end\n");

    *data_size = sizeof(AVPicture);
    *picture = s->picture;

    return s->buf - s->buf_start;
}

static void j2kdec_init(AVCodecContext *avctx)
{
    J2kDecoderContext *s = avctx->priv_data;

    avcodec_get_frame_defaults((AVFrame*)&s->picture);
    avctx->coded_frame = (AVFrame*)&s->picture;
}

AVCodec jpeg2000_decoder = {
    "j2k",
    CODEC_TYPE_VIDEO,
    CODEC_ID_JPEG2000,
    sizeof(J2kDecoderContext),
    j2kdec_init,
    NULL,
    NULL,
    decode_frame,
    0,
    .pix_fmts =
        (enum PixelFormat[]) {PIX_FMT_GRAY8, PIX_FMT_RGB24, -1}
};
