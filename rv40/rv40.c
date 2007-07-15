/*
 * RV40 decoder
 * Copyright (c) 2007 Mike Melanson, Konstantin Shishkov
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
 * @file rv40.c
 * RV40 decoder.
 */

#include "avcodec.h"
#include "dsputil.h"
#include "mpegvideo.h"

#include "rv40vlc.h"
#include "rv40vlc2.h"
#include "rv40data.h"

#include "h264pred.h"

//#define DEBUG


/**
 * VLC tables used by decoder
 *
 * intra frame VLC sets do not contain some of those tables
 */
typedef struct RV40VLC{
    VLC cbppattern[2];     ///< VLCs used for pattern of coded block patterns decoding
    VLC cbp[2][4];         ///< VLCs used for coded block patterns decoding
    VLC first_pattern[4];  ///< VLCs used for decoding coefficients in the first subblock
    VLC second_pattern[2]; ///< VLCs used for decoding coefficients in the subblocks 2 and 3
    VLC third_pattern[2];  ///< VLCs used for decoding coefficients in the last subblock
    VLC coefficient;       ///< VLCs used for decoding big coefficients
}RV40VLC;

/** Decoder context */
typedef struct RV40DecContext{
    MpegEncContext s;
    int mb_bits;             ///< bits needed to read MB offet in slice header
    int *intra_types_hist;   ///< old block types, used for prediction
    int *intra_types;        ///< block types
    int intra_types_stride;  ///< stride for block types data
    int block_start;         ///< start of slice in blocks
    int ptype;               ///< picture type
    int quant;               ///< quantizer

    RV40VLC *cur_vlcs;       ///< VLC set used for current frame decoding
    int bits;                ///< slice size in bits
    H264PredContext h;       ///< functions for 4x4 and 16x16 intra block prediction
}RV40DecContext;

static RV40VLC intra_vlcs[NUM_INTRA_TABLES], inter_vlcs[NUM_INTER_TABLES];
static VLC aic_top_vlc;
static VLC aic_mode1_vlc[AIC_MODE1_NUM], aic_mode2_vlc[AIC_MODE2_NUM];

/**
 * @defgroup vlc RV40 VLC generating functions
 * @{
 */

/**
 * Generate VLC from codeword lengths
 */
static int rv40_gen_vlc(const uint8_t *bits2, int size, VLC *vlc)
{
    int i;
    int counts[17], codes[17];
    uint16_t *cw, *syms;
    uint8_t *bits;
    int maxbits = 0, realsize;
    int ret;

    memset(counts, 0, 16 * sizeof(int));

    cw = av_mallocz(size * 2);
    syms = av_malloc(size * 2);
    bits = av_malloc(size);

    realsize = 0;
    for(i = 0; i < size; i++){
        if(bits2[i]){
            bits[realsize] = bits2[i];
            syms[realsize] = i;
            realsize++;
            if(bits2[i] > maxbits)
                maxbits = bits2[i];
        }
    }

    size = realsize;
    for(i = 0; i < size; i++)
        counts[bits[i]]++;
    codes[0] = 0;
    for(i = 0; i < 16; i++)
        codes[i+1] = (codes[i] + counts[i]) << 1;
    for(i = 0; i < realsize; i++)
        cw[i] = codes[bits[i]]++;

    ret = init_vlc_sparse(vlc, FFMIN(maxbits, 9), size,
                          bits, 1, 1,
                          cw,   2, 2,
                          syms, 2, 2, INIT_VLC_USE_STATIC);
    av_free(cw);
    av_free(syms);
    av_free(bits);
    return ret;
}

/**
 * Initialize all tables
 */
static void rv40_init_tables()
{
    int i, j, k;

    for(i = 0; i < NUM_INTRA_TABLES; i++){
        for(j = 0; j < 2; j++)
            rv40_gen_vlc(rv40_intra_cbppatvlc_pointers[i][j], CBPPAT_VLC_SIZE, &intra_vlcs[i].cbppattern[j]);
        for(j = 0; j < 2; j++)
            for(k = 0; k < 4; k++)
                rv40_gen_vlc(rv40_intra_cbpvlc_pointers[i][j][k], CBP_VLC_SIZE, &intra_vlcs[i].cbp[j][k]);
        for(j = 0; j < 4; j++)
            rv40_gen_vlc(rv40_intra_firstpatvlc_pointers[i][j], FIRSTBLK_VLC_SIZE, &intra_vlcs[i].first_pattern[j]);
        for(j = 0; j < 2; j++)
            rv40_gen_vlc(rv40_intra_secondpatvlc_pointers[i][j], OTHERBLK_VLC_SIZE, &intra_vlcs[i].second_pattern[j]);
        for(j = 0; j < 2; j++)
            rv40_gen_vlc(rv40_intra_thirdpatvlc_pointers[i][j], OTHERBLK_VLC_SIZE, &intra_vlcs[i].third_pattern[j]);
        rv40_gen_vlc(rv40_intra_coeffvlc_pointers[i], COEFF_VLC_SIZE, &intra_vlcs[i].coefficient);
    }

    for(i = 0; i < NUM_INTER_TABLES; i++){
        rv40_gen_vlc(rv40_inter_cbppatvlc_pointers[i], CBPPAT_VLC_SIZE, &inter_vlcs[i].cbppattern[0]);
        for(j = 0; j < 4; j++)
            rv40_gen_vlc(rv40_inter_cbpvlc_pointers[i][j], CBP_VLC_SIZE, &inter_vlcs[i].cbp[0][j]);
        for(j = 0; j < 2; j++)
            rv40_gen_vlc(rv40_inter_firstpatvlc_pointers[i][j], FIRSTBLK_VLC_SIZE, &inter_vlcs[i].first_pattern[j]);
        for(j = 0; j < 2; j++)
            rv40_gen_vlc(rv40_inter_secondpatvlc_pointers[i][j], OTHERBLK_VLC_SIZE, &inter_vlcs[i].second_pattern[j]);
        for(j = 0; j < 2; j++)
            rv40_gen_vlc(rv40_inter_thirdpatvlc_pointers[i][j], OTHERBLK_VLC_SIZE, &inter_vlcs[i].third_pattern[j]);
        rv40_gen_vlc(rv40_inter_coeffvlc_pointers[i], COEFF_VLC_SIZE, &inter_vlcs[i].coefficient);
    }

    init_vlc(&aic_top_vlc, AIC_TOP_BITS, AIC_TOP_SIZE,
             rv40_aic_top_vlc_bits,  1, 1,
             rv40_aic_top_vlc_codes, 1, 1, INIT_VLC_USE_STATIC);
    for(i = 0; i < AIC_MODE1_NUM; i++){
        // For some reason every tenth VLC table is empty
        // So skip it for consistency
        // XXX: redo without this hack
        if((i % 10) == 9) continue;
        init_vlc(&aic_mode1_vlc[i], AIC_MODE1_BITS, AIC_MODE1_SIZE,
                 aic_mode1_vlc_bits[i],  1, 1,
                 aic_mode1_vlc_codes[i], 1, 1, INIT_VLC_USE_STATIC);
    }
    for(i = 0; i < AIC_MODE2_NUM; i++){
        init_vlc(&aic_mode2_vlc[i], AIC_MODE2_BITS, AIC_MODE2_SIZE,
                 aic_mode2_vlc_bits[i],  1, 1,
                 aic_mode2_vlc_codes[i], 2, 2, INIT_VLC_USE_STATIC);
    }
}

/** @} */ // vlc group


/**
 * @defgroup transform RV40 inverse transform functions
 * @{
 */

/**
 * Real Video 4.0 inverse transform
 * Code is almost the same as in SVQ3, only scaling is different
 */
static void rv40_intra_inv_transform(DCTELEM *block, const int offset){
    int temp[16];
    unsigned int i;

    for(i=0; i<4; i++){
        const int z0= 13*(block[offset+i+8*0] +    block[offset+i+8*2]);
        const int z1= 13*(block[offset+i+8*0] -    block[offset+i+8*2]);
        const int z2=  7* block[offset+i+8*1] - 17*block[offset+i+8*3];
        const int z3= 17* block[offset+i+8*1] +  7*block[offset+i+8*3];

        temp[4*i+0]= z0+z3;
        temp[4*i+1]= z1+z2;
        temp[4*i+2]= z1-z2;
        temp[4*i+3]= z0-z3;
    }

    for(i=0; i<4; i++){
        const int z0= 13*(temp[4*0+i] +    temp[4*2+i]);
        const int z1= 13*(temp[4*0+i] -    temp[4*2+i]);
        const int z2=  7* temp[4*1+i] - 17*temp[4*3+i];
        const int z3= 17* temp[4*1+i] +  7*temp[4*3+i];

        block[offset+i*8+0]= ((z0 + z3) + 0x200)>>10;
        block[offset+i*8+1]= ((z1 + z2) + 0x200)>>10;
        block[offset+i*8+2]= ((z1 - z2) + 0x200)>>10;
        block[offset+i*8+3]= ((z0 - z3) + 0x200)>>10;
    }

}

/**
 * RealVideo 4.0 inverse transform - special version
 *
 * Code is almost the same but final coefficients are multiplied by 1.5
 * and have no rounding
 */
static void rv40_intra_inv_transform_noround(DCTELEM *block, const int offset){
    int temp[16];
    unsigned int i;

    for(i=0; i<4; i++){
        const int z0= 13*(block[offset+i+8*0] +    block[offset+i+8*2]);
        const int z1= 13*(block[offset+i+8*0] -    block[offset+i+8*2]);
        const int z2=  7* block[offset+i+8*1] - 17*block[offset+i+8*3];
        const int z3= 17* block[offset+i+8*1] +  7*block[offset+i+8*3];

        temp[4*i+0]= z0+z3;
        temp[4*i+1]= z1+z2;
        temp[4*i+2]= z1-z2;
        temp[4*i+3]= z0-z3;
    }

    for(i=0; i<4; i++){
        const int z0= 13*(temp[4*0+i] +    temp[4*2+i]);
        const int z1= 13*(temp[4*0+i] -    temp[4*2+i]);
        const int z2=  7* temp[4*1+i] - 17*temp[4*3+i];
        const int z3= 17* temp[4*1+i] +  7*temp[4*3+i];

        block[offset+i*8+0]= ((z0 + z3)*3)>>11;
        block[offset+i*8+1]= ((z1 + z2)*3)>>11;
        block[offset+i*8+2]= ((z1 - z2)*3)>>11;
        block[offset+i*8+3]= ((z0 - z3)*3)>>11;
    }

}

/** @} */ // transform


/**
 * @defgroup block RV40 4x4 block decoding functions
 * @{
 */

/**
 * Decode coded block pattern
 */
static int rv40_decode_cbp(GetBitContext *gb, RV40VLC *vlc, int table)
{
    int pattern, code, cbp=0;
    int table2;
    static const int cbp_masks[4] = {0x000000, 0x100000, 0x010000, 0x110000};
    static const int shifts[4] = { 0, 2, 8, 10 };
    int *curshift = shifts;
    int i, t, mask;

    code = get_vlc2(gb, vlc->cbppattern[table].table, 9, 2);
    pattern = code & 0xF;
    code >>= 4;

    table2 = rv40_count_ones[pattern];

    for(mask = 8; mask; mask >>= 1, curshift++){
        if(!(pattern & mask)) continue;
        t = get_vlc2(gb, vlc->cbp[table][table2].table, vlc->cbp[table][table2].bits, 1);
        cbp |= rv40_cbp_code[t] << curshift[0];
    }

    for(i = 0; i < 4; i++){
        t = modulo_three_table[code][i];
        if(t == 1)
            cbp |= cbp_masks[1+get_bits1(gb)] << i;
        if(t == 2)
            cbp |= cbp_masks[3] << i;
    }
    return cbp;
}

/**
 * Get one coefficient value from bistream and store it
 */
static inline void decode_coeff(DCTELEM *dst, int coef, int esc, GetBitContext *gb, VLC* vlc)
{
    if(coef){
        if(coef == esc){
            coef = get_vlc2(gb, vlc->table, 9, 2);
            if(coef > 23){
                coef -= 23;
                coef = 22 + ((1 << coef) | get_bits(gb, coef));
            }
            coef += esc;
        }
        if(get_bits1(gb))
            coef = -coef;
        *dst = coef;
    }
}

static inline void decode_subblock(DCTELEM *dst, int coeffs[4], GetBitContext *gb, VLC *vlc)
{
    decode_coeff(dst  , coeffs[0], 3, gb, vlc);
    decode_coeff(dst+1, coeffs[1], 2, gb, vlc);
    decode_coeff(dst+8, coeffs[2], 2, gb, vlc);
    decode_coeff(dst+9, coeffs[3], 2, gb, vlc);
}

// slightly modified version for the third subblock
static inline void decode_subblock2(DCTELEM *dst, int coeffs[4], GetBitContext *gb, VLC *vlc)
{
    decode_coeff(dst  , coeffs[0], 3, gb, vlc);
    decode_coeff(dst+8, coeffs[1], 2, gb, vlc);
    decode_coeff(dst+1, coeffs[2], 2, gb, vlc);
    decode_coeff(dst+9, coeffs[3], 2, gb, vlc);
}

/**
 * Decode coefficients for 4x4 block
 *
 * This is done by filling 2x2 subblocks with decoded coefficients
 * in this order (the same for subblocks and subblock coefficients):
 *  o--o
 *    /
 *   /
 *  o--o
 */

static inline void rv40_decode_block(DCTELEM *dst, GetBitContext *gb, RV40VLC *rvlc, int fc, int sc)
{
    int code, pattern;
    int coeffs[4];

    code = get_vlc2(gb, rvlc->first_pattern[fc].table, 9, 2);

    pattern = code & 0x7;

    code >>= 3;
    coeffs[0] = modulo_three_table[code][0];
    coeffs[1] = modulo_three_table[code][1];
    coeffs[2] = modulo_three_table[code][2];
    coeffs[3] = modulo_three_table[code][3];
    decode_subblock(dst, coeffs, gb, &rvlc->coefficient);

    if(pattern & 4){
        code = get_vlc2(gb, rvlc->second_pattern[sc].table, 9, 2);
        coeffs[0] = modulo_three_table[code][0];
        coeffs[1] = modulo_three_table[code][1];
        coeffs[2] = modulo_three_table[code][2];
        coeffs[3] = modulo_three_table[code][3];
        decode_subblock(dst + 2, coeffs, gb, &rvlc->coefficient);
    }
    if(pattern & 2){ // Looks like coefficients 1 and 2 are swapped for this block
        code = get_vlc2(gb, rvlc->second_pattern[sc].table, 9, 2);
        coeffs[0] = modulo_three_table[code][0];
        coeffs[1] = modulo_three_table[code][1];
        coeffs[2] = modulo_three_table[code][2];
        coeffs[3] = modulo_three_table[code][3];
        decode_subblock2(dst + 8*2, coeffs, gb, &rvlc->coefficient);
    }
    if(pattern & 1){
        code = get_vlc2(gb, rvlc->third_pattern[sc].table, 9, 2);
        coeffs[0] = modulo_three_table[code][0];
        coeffs[1] = modulo_three_table[code][1];
        coeffs[2] = modulo_three_table[code][2];
        coeffs[3] = modulo_three_table[code][3];
        decode_subblock(dst + 8*2+2, coeffs, gb, &rvlc->coefficient);
    }
}

/**
 * Dequantize ordinary 4x4 block
 * @todo optimize
 */
static inline void rv40_dequant4x4(DCTELEM *block, int offset, int Qdc, int Q)
{
    int i, j;

    block += offset;
    block[0] = (block[0] * Qdc + 8) >> 4;
    for(i = 0; i < 4; i++)
        for(j = 0; j < 4; j++)
            if(i || j)
                block[j + i*8] = (block[j + i*8] * Q + 8) >> 4;
}

/**
 * Dequantize 4x4 block of DC values for 16x16 macroblock
 * @todo optimize
 */
static inline void rv40_dequant4x4_16x16(DCTELEM *block, int offset, int Qdc, int Q)
{
    int i;

    block += offset;
    for(i = 0; i < 3; i++)
         block[rv40_dezigzag[i]] = (block[rv40_dezigzag[i]] * Qdc + 8) >> 4;
    for(; i < 16; i++)
         block[rv40_dezigzag[i]] = (block[rv40_dezigzag[i]] * Q + 8) >> 4;
}
/** @} */ //block functions


/**
 * @defgroup bitstream RV40 bitstream parsing
 * @{
 */

static inline int decode210(GetBitContext *gb){
    if (get_bits1(gb))
        return 0;
    else
        return 2 - get_bits1(gb);
}

/**
 * Get stored dimension from bitstream
 *
 * If the width/height is the standard one then it's coded as 3-bit index.
 * Otherwise it is coded as escaped 8-bit portions.
 */
static inline int get_dimension(GetBitContext *gb, const int *dim1, const int *dim2)
{
    int val, t;

    t = get_bits(gb, 3);
    val = dim1[t];
    if(!val && dim2)
        val = dim2[(t | get_bits1(gb)) & 3];
    if(!val){
        do{
            t = get_bits(gb, 8);
            val += t << 2;
        }while(t == 0xFF);
    }
    return val;
}

/**
 * Get encoded picture size - usually this is called from rv40_parse_slice_header
 */
static void rv40_parse_picture_size(GetBitContext *gb, int *w, int *h)
{
    *w = get_dimension(gb, rv40_standard_widths, NULL);
    *h = get_dimension(gb, rv40_standard_heights, rv40_standard_heights2);
}

static int rv40_parse_slice_header(RV40DecContext *r, GetBitContext *gb)
{
    int t, mb_bits;
    int w, h;

    if(get_bits1(gb))
        return -1;
    r->ptype = get_bits(gb, 2);
    r->quant = get_bits(gb, 5);
    if(get_bits(gb, 2))
        return -1;
    get_bits(gb, 2); /// ???
    if(get_bits1(gb))
        return -1;
    t = get_bits(gb, 13); /// ???
    rv40_parse_picture_size(gb, &w, &h);
    r->s.avctx->coded_width  = w;
    r->s.avctx->coded_height = h;
    mb_bits = av_log2((w + 7) >> 3) + av_log2((h + 7) >> 3);
    r->block_start = get_bits(gb, mb_bits);

    return 0;
}

/**
 * Decode 4x4 intra types array
 */
static int rv40_decode_intra_types(RV40DecContext *r, GetBitContext *gb, int *dst)
{
    MpegEncContext *s = &r->s;
    int i, j, k, v;
    int A, B, C;
    int pattern;
    int *ptr;

    for(i = 0; i < 4; i++, dst += r->intra_types_stride){
        if(!i && s->first_slice_line){
            pattern = get_vlc2(gb, aic_top_vlc.table, AIC_TOP_BITS, 1);
            dst[0] = (pattern >> 2) & 2;
            dst[1] = (pattern >> 1) & 2;
            dst[2] =  pattern       & 2;
            dst[3] = (pattern << 1) & 2;
            continue;
        }
        ptr = dst;
        for(j = 0; j < 4; j++){
            /* Coefficients are read using VLC chosen by prediction pattern
             * First one (used for retrieving a pair of coefficients) is
             * constructed from top, top right and left coefficients
             * Second one (used for retrieving only one coefficient) is
             * top + 10 * left
             */
            A = ptr[-r->intra_types_stride + 1]; // it won't be used for the last coefficient in a row
            B = ptr[-r->intra_types_stride];
            C = ptr[-1];
            pattern = A + (B << 4) + (C << 8);
            for(k = 0; k < MODE2_PATTERNS_NUM; k++)
                if(pattern == rv40_aic_table_index[k])
                    break;
            if(j < 3 && k < MODE2_PATTERNS_NUM){ //pattern is found, decoding 2 coefficients
                v = get_vlc2(gb, aic_mode2_vlc[k].table, AIC_MODE2_BITS, 2);
                *ptr++ = v/9;
                *ptr++ = v%9;
                j++;
            }else{
                if(B != -1 && C != -1)
                    v = get_vlc2(gb, aic_mode1_vlc[B + C*10].table, AIC_MODE1_BITS, 1);
                else{ // tricky decoding
                    v = 0;
                    switch(C){
                    case -1: // code 0 -> 1, 1 -> 0
                        if(B == -1 || B == 0 || B == 1)
                            v = get_bits1(gb) ^ 1;
                        break;
                    case  0:
                    case  2: // code 0 -> 2, 1 -> 0
                        v = (get_bits1(gb) ^ 1) << 1;
                        break;
                    }
                }
                *ptr++ = v;
            }
        }
    }
    return 0;
}

/** Mapping of RV40 intra prediction types to standard H.264 types */
static const int ittrans[9] = {
 HOR_PRED, // almost ok
 VERT_PRED, // ok
 LEFT_DC_PRED, // ok?
 DIAG_DOWN_RIGHT_PRED, // ok
 DIAG_DOWN_LEFT_PRED,
 VERT_RIGHT_PRED, // ok
 VERT_LEFT_PRED, //ok
 HOR_DOWN_PRED,
 HOR_UP_PRED,
};

static void rv40_output_macroblock(RV40DecContext *r, int *intra_types, int cbp, int is16)
{
    MpegEncContext *s = &r->s;
    DSPContext *dsp = &s->dsp;
    int i, j, x, y;
    uint8_t *Y, *YY, *PY;
    int no_up, no_left, itype;

    no_up = s->first_slice_line;
    Y = s->dest[0];
    if(!is16){
        for(j = 0; j < 4; j++){
            no_left = !s->mb_x || (s->mb_x == s->resync_mb_x && s->first_slice_line);
            for(YY = Y, i = 0; i < 4; i++, cbp >>= 1, no_left = 0, YY += 4){
                itype = ittrans[intra_types[i]];
                if(no_up && no_left)
                    itype = DC_128_PRED;
                else if(no_up)
                    itype = HOR_PRED;
                else if(no_left)
                    itype = VERT_PRED;
                /* silly RV40 predictor assumes that blocks are aligned in one row
                   so first block may reference block from the previous row
                   or last block from the previous macroblock
                 */
                if(!i && !j)
                    PY = YY - 1 + s->linesize * 12;
                else
                    PY = YY - 1;
                r->h.pred4x4[itype](YY, PY, YY - s->linesize + 4, s->linesize);
                if(!(cbp & 1)) continue;
                /* add_pixels_clamped for 4x4 block */
                for(y = 0; y < 4; y++)
                    for(x = 0; x < 4; x++)
                        YY[x + y * s->linesize] = av_clip_uint8(YY[x + y * s->linesize] + s->block[(i>>1)+(j&2)][(i&1)*4+(j&1)*32+x+y*8]);
            }
            no_up = 0;
            Y += s->linesize * 4;
            intra_types += r->intra_types_stride;
        }
    }else{
        no_left = !s->mb_x || (s->mb_x == s->resync_mb_x && s->first_slice_line);
        if(no_up && no_left)
            r->h.pred16x16[DC_128_PRED8x8](Y, s->linesize);
        else if(no_up)
            r->h.pred16x16[HOR_PRED8x8](Y, s->linesize);
        else if(no_left)
            r->h.pred16x16[VERT_PRED8x8](Y, s->linesize);
        else
            r->h.pred16x16[intra_types[0]](Y, s->linesize);
        dsp->add_pixels_clamped(s->block[0], Y, s->current_picture.linesize[0]);
        dsp->add_pixels_clamped(s->block[1], Y + 8, s->current_picture.linesize[0]);
        Y += s->current_picture.linesize[0] * 8;
        dsp->add_pixels_clamped(s->block[2], Y, s->current_picture.linesize[0]);
        dsp->add_pixels_clamped(s->block[3], Y + 8, s->current_picture.linesize[0]);
    }
}

static int rv40_decode_macroblock(RV40DecContext *r, int *intra_types)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int q, cbp, cbp2;
    int i, blknum, blkoff;
    int luma_vlc, chroma_vlc;
    int is16 = 0;
    DCTELEM block16[64];

    q = decode210(gb);
    switch(q){
    case 0: // 16x16 block
        is16 = 1;
        break;
    case 1:
        break;
    case 2:
        av_log(NULL,0,"Need DQUANT\n");
        // q = decode_dquant(gb);
        break;
    }
    if(!is16){
        rv40_decode_intra_types(r, gb, intra_types);
        chroma_vlc = 0;
        luma_vlc   = 1;
    }else{
        q = get_bits(gb, 2);
        for(i = 0; i < 16; i++)
            intra_types[(i & 3) + (i>>2) * r->intra_types_stride] = q;
        chroma_vlc = 0;
        luma_vlc   = 2;
    }
    cbp = cbp2 = rv40_decode_cbp(gb, r->cur_vlcs, is16);

    if(is16){
        memset(block16, 0, sizeof(block16));
        rv40_decode_block(block16, gb, r->cur_vlcs, 3, 0);
        rv40_dequant4x4_16x16(block16, 0, rv40_qscale_tab[r->quant],rv40_qscale_tab[r->quant]);
        rv40_intra_inv_transform_noround(block16, 0);
    }

    for(i = 0; i < 16; i++, cbp >>= 1){
        if(!is16 && !(cbp & 1)) continue;
        blknum = ((i & 2) >> 1) + ((i & 8) >> 2);
        blkoff = ((i & 1) << 2) + ((i & 4) << 3);
        if(cbp & 1)
            rv40_decode_block(s->block[blknum] + blkoff, gb, r->cur_vlcs, luma_vlc, 0);
        if((cbp & 1) || is16){
            if(is16) //FIXME: optimize
                s->block[blknum][blkoff] = block16[i];
            rv40_dequant4x4(s->block[blknum], blkoff, rv40_qscale_tab[r->quant],rv40_qscale_tab[r->quant]);
            rv40_intra_inv_transform(s->block[blknum], blkoff);
        }
    }
    for(; i < 24; i++, cbp >>= 1){
        if(!(cbp & 1)) continue;
        blknum = ((i & 4) >> 2) + 4;
        blkoff = ((i & 1) << 2) + ((i & 2) << 4);
        rv40_decode_block(s->block[blknum] + blkoff, gb, r->cur_vlcs, chroma_vlc, 1);
        rv40_dequant4x4(s->block[blknum], blkoff, rv40_qscale_tab[r->quant],rv40_qscale_tab[r->quant]);
        rv40_intra_inv_transform(s->block[blknum], blkoff);
    }
    rv40_output_macroblock(r, intra_types, cbp2, is16);

    return 0;
}

static int rv40_decode_slice(RV40DecContext *r)
{
    MpegEncContext *s = &r->s;
    int mb_pos;

    mb_pos = s->mb_x + s->mb_y * s->mb_width;
    if(r->block_start != mb_pos){
        av_log(s->avctx, AV_LOG_ERROR, "Slice indicates MB offset %d, got %d\n", r->block_start, mb_pos);
        s->mb_x = r->block_start % s->mb_width;
        s->mb_y = r->block_start / s->mb_width;
    }
    memset(r->intra_types_hist, -1, r->intra_types_stride * 4 * 2 * sizeof(int));
    s->first_slice_line = 1;
    s->resync_mb_x= s->mb_x;
    s->resync_mb_y= s->mb_y;
    ff_init_block_index(s);
    //XXX: better bounds detection?
    while((get_bits_count(&s->gb) + 5 < r->bits) && (s->mb_y < s->mb_height)) {
        ff_update_block_index(s);
        s->dsp.clear_blocks(s->block[0]);

        rv40_decode_macroblock(r, r->intra_types + (s->mb_x + 1) * 4);
        ff_draw_horiz_band(s, s->mb_y * 16, 16);
        if (++s->mb_x == s->mb_width) {
            s->mb_x = 0;
            s->mb_y++;
            ff_init_block_index(s);

            memmove(r->intra_types_hist, r->intra_types, r->intra_types_stride * 4 * sizeof(int));
            memset(r->intra_types, -1, r->intra_types_stride * 4 * sizeof(int));
        }
        if(s->mb_x == s->resync_mb_x)
            s->first_slice_line=0;
    }
    ff_er_add_slice(s, s->resync_mb_x, s->resync_mb_y, s->mb_x-1, s->mb_y, AC_END|DC_END|MV_END);

    return 0;
}

/** @} */ //bitstream functions


/**
 * Initialize decoder
 * @todo Maybe redone in some other way
 */
static int rv40_decode_init(AVCodecContext *avctx)
{
    RV40DecContext *r = avctx->priv_data;
    MpegEncContext *s = &r->s;

    static int tables_done = 0;

    MPV_decode_defaults(s);
    s->avctx= avctx;
    s->out_format = FMT_H264;
    s->codec_id= avctx->codec_id;

    s->width = avctx->width;
    s->height = avctx->height;

    r->s.avctx = avctx;
    avctx->flags |= CODEC_FLAG_EMU_EDGE;
    r->s.flags |= CODEC_FLAG_EMU_EDGE;
    avctx->pix_fmt = PIX_FMT_YUV420P;
    s->low_delay = 1;

    if (MPV_common_init(s) < 0)
        return -1;

    ff_h264_pred_init(&r->h);

    r->intra_types_stride = (s->mb_width + 1) * 4;
    r->intra_types_hist = av_malloc(r->intra_types_stride * 4 * 2 * sizeof(int));
    r->intra_types = r->intra_types_hist + r->intra_types_stride * 4;

    if(!tables_done){
        rv40_init_tables();
        tables_done = 1;
    }
    return 0;
}


static int rv40_decode_frame(AVCodecContext *avctx,
                            void *data, int *data_size,
                            uint8_t *buf, int buf_size)
{
    RV40DecContext *r = avctx->priv_data;
    MpegEncContext *s = &r->s;
    AVFrame *pict = data;
    int res;

    /* no supplementary picture */
    if (buf_size == 0) {
        return 0;
    }

    init_get_bits(&s->gb, buf, buf_size*8);
    r->bits = buf_size * 8;
    rv40_parse_slice_header(r, &r->s.gb);
    s->pict_type = r->ptype+1;

if(s->pict_type != I_TYPE)return -1;
    if ((s->mb_x == 0 && s->mb_y == 0) || s->current_picture_ptr==NULL) {
        if(s->current_picture_ptr){ //FIXME write parser so we always have complete frames?
            ff_er_frame_end(s);
            MPV_frame_end(s);
            s->mb_x= s->mb_y = s->resync_mb_x = s->resync_mb_y= 0;
        }
        if(MPV_frame_start(s, avctx) < 0)
            return -1;
        ff_er_frame_start(s);
        r->cur_vlcs = &intra_vlcs[2];
    }
    rv40_decode_slice(r);

    if(s->current_picture_ptr != NULL && s->mb_y>=s->mb_height){
        ff_er_frame_end(s);
        MPV_frame_end(s);
        if (s->pict_type == B_TYPE || s->low_delay) {
            *pict= *(AVFrame*)s->current_picture_ptr;
        } else if (s->last_picture_ptr != NULL) {
            *pict= *(AVFrame*)s->last_picture_ptr;
        }

        if(s->last_picture_ptr || s->low_delay){
            *data_size = sizeof(AVFrame);
            ff_print_debug_info(s, pict);
        }
        s->current_picture_ptr= NULL; //so we can detect if frame_end wasnt called (find some nicer solution...)
    }

    return buf_size;
}

static int rv40_decode_end(AVCodecContext *avctx)
{
    RV40DecContext *r = avctx->priv_data;

    MPV_common_end(&r->s);

    av_freep(&r->intra_types_hist);
    r->intra_types = NULL;

    return 0;
}

AVCodec rv40_decoder = {
    "rv40",
    CODEC_TYPE_VIDEO,
    CODEC_ID_RV40,
    sizeof(RV40DecContext),
    rv40_decode_init,
    NULL,
    rv40_decode_end,
    rv40_decode_frame,
};
