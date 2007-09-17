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
 * RV30 and RV40 decoder.
 */

#include "avcodec.h"
#include "dsputil.h"
#include "mpegvideo.h"

#include "rv40vlc.h"
#include "rv40vlc2.h"
#include "rv40data.h"
#include "rv30data.h"

#include "h264pred.h"

//#define DEBUG

/** Translation of RV40 macroblock types to lavc ones */
static const int rv40_mb_type_to_lavc[12] = {
    MB_TYPE_INTRA, MB_TYPE_INTRA16x16, MB_TYPE_16x16,   MB_TYPE_8x8,
    MB_TYPE_16x16, MB_TYPE_16x16,      MB_TYPE_SKIP,    MB_TYPE_DIRECT2,
    MB_TYPE_16x8,  MB_TYPE_8x16,       MB_TYPE_DIRECT2, MB_TYPE_16x16
};

/**
 * RV40 Macroblock types
 * @{
 */
enum RV40BlockTypes{
    RV40_MB_TYPE_INTRA,      ///< Intra macroblock
    RV40_MB_TYPE_INTRA16x16, ///< Intra macroblock with DCs in a separate 4x4 block
    RV40_MB_P_16x16,         ///< P-frame macroblock, one motion frame
    RV40_MB_P_8x8,           ///< P-frame macroblock, 8x8 motion compensation partitions
    RV40_MB_B_FORWARD,       ///< B-frame macroblock, forward prediction
    RV40_MB_B_BACKWARD,      ///< B-frame macroblock, backward prediction
    RV40_MB_SKIP,            ///< Skipped block
    RV40_MB_B_INTERP,        ///< Bidirectionally predicted B-frame macroblock, no motion vectors
    RV40_MB_P_16x8,          ///< P-frame macroblock, 16x8 motion compensation partitions
    RV40_MB_P_8x16,          ///< P-frame macroblock, 8x16 motion compensation partitions
    RV40_MB_B_DIRECT,        ///< Bidirectionally predicted B-frame macroblock, two motion vectors
    RV40_MB_P_MIX16x16,      ///< P-frame macroblock with DCs in a separate 4x4 block, one motion vector
    RV40_MB_TYPES
};
/** @} */

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

/** Essential slice information */
typedef struct SliceInfo{
    int type;              ///< slice type (intra, inter)
    int size;              ///< size of the slice in bits
    int quant;             ///< quantizer used for this slice
    int vlc_set;           ///< VLCs used for this slice
    int start, end;        ///< start and end macroblocks of the slice
    int header_size;       ///< header size in bits
    int width;             ///< coded width
    int height;            ///< coded height
}SliceInfo;

/** Slice information saved for truncated slices */
typedef struct SavedSliceInfo{
    uint8_t *data;         ///< bitstream data
    int data_size;         ///< data size
    int bits_used;         ///< bits used up to last decoded block
    int mb_x, mb_y;        ///< coordinates of the last decoded block
}SavedSliceInfo;

/** Decoder context */
typedef struct RV40DecContext{
    MpegEncContext s;
    int mb_bits;             ///< bits needed to read MB offet in slice header
    int *intra_types_hist;   ///< old block types, used for prediction
    int *intra_types;        ///< block types
    int intra_types_stride;  ///< stride for block types data
    int block_start;         ///< start of slice in blocks

    int vlc_set;             ///< index of currently selected VLC set
    RV40VLC *cur_vlcs;       ///< VLC set used for current frame decoding
    int bits;                ///< slice size in bits
    H264PredContext h;       ///< functions for 4x4 and 16x16 intra block prediction
    SliceInfo si;            ///< current slice information
    SliceInfo prev_si;       ///< info for the saved slice
    uint8_t *slice_data;     ///< saved slice data

    int *mb_type;            ///< internal macroblock types
    int block_type;          ///< current block type
    int luma_vlc;            ///< which VLC set will be used for luma blocks decoding
    int chroma_vlc;          ///< which VLC set will be used for chroma blocks decoding
    int is16;                ///< current block has additional 16x16 specific features or not
    int dmv[4][2];           ///< differential motion vectors for the current macroblock

    int truncated;           ///< flag signalling that slice ended prematurely
    SavedSliceInfo ssi;      ///< data for truncated slice

    int rv30;                ///< indicates which RV variasnt is currently decoded
    int rpr;                 ///< one field size in RV30 slice header

    int avail[4];            ///< whether left, top, top rights and top left MBs are available
}RV40DecContext;

static RV40VLC intra_vlcs[NUM_INTRA_TABLES], inter_vlcs[NUM_INTER_TABLES];
static VLC aic_top_vlc;
static VLC aic_mode1_vlc[AIC_MODE1_NUM], aic_mode2_vlc[AIC_MODE2_NUM];
static VLC mbinfo_vlc, ptype_vlc[NUM_PTYPE_VLCS], btype_vlc[NUM_BTYPE_VLCS];

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
    int counts[17] = {0}, codes[17];
    uint16_t cw[size], syms[size];
    uint8_t bits[size];
    int maxbits = 0, realsize;
    int ret;

    realsize = 0;
    for(i = 0; i < size; i++){
        if(bits2[i]){
            bits[realsize] = bits2[i];
            syms[realsize] = i;
            realsize++;
            maxbits = FFMAX(maxbits, bits2[i]);
            counts[bits2[i]]++;
        }
    }

    size = realsize;
    codes[0] = 0;
    for(i = 0; i < 16; i++)
        codes[i+1] = (codes[i] + counts[i]) << 1;
    for(i = 0; i < realsize; i++)
        cw[i] = codes[bits[i]]++;

    ret = init_vlc_sparse(vlc, FFMIN(maxbits, 9), size,
                          bits, 1, 1,
                          cw,   2, 2,
                          syms, 2, 2, INIT_VLC_USE_STATIC);
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
    init_vlc_sparse(&mbinfo_vlc, MBINFO_BITS, NUM_MBINFO,
                    mbinfo_vlc_bits,  1, 1,
                    mbinfo_vlc_codes, 1, 1,
                    mbinfo_vlc_syms,  1, 1, INIT_VLC_USE_STATIC);
    for(i = 0; i < NUM_PTYPE_VLCS; i++)
         init_vlc_sparse(&ptype_vlc[i], PTYPE_VLC_BITS, PTYPE_VLC_SIZE,
                         ptype_vlc_bits[i],  1, 1,
                         ptype_vlc_codes[i], 1, 1,
                         ptype_vlc_syms,     1, 1, INIT_VLC_USE_STATIC);
    for(i = 0; i < NUM_BTYPE_VLCS; i++)
         init_vlc_sparse(&btype_vlc[i], BTYPE_VLC_BITS, BTYPE_VLC_SIZE,
                         btype_vlc_bits[i],  1, 1,
                         btype_vlc_codes[i], 1, 1,
                         btype_vlc_syms,     1, 1, INIT_VLC_USE_STATIC);
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
        const int z0= 13*(temp[4*0+i] +    temp[4*2+i]) + 0x200;
        const int z1= 13*(temp[4*0+i] -    temp[4*2+i]) + 0x200;
        const int z2=  7* temp[4*1+i] - 17*temp[4*3+i];
        const int z3= 17* temp[4*1+i] +  7*temp[4*3+i];

        block[offset+i*8+0]= (z0 + z3)>>10;
        block[offset+i*8+1]= (z1 + z2)>>10;
        block[offset+i*8+2]= (z1 - z2)>>10;
        block[offset+i*8+3]= (z0 - z3)>>10;
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
    static const int cbp_masks[3] = {0x100000, 0x010000, 0x110000};
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
            cbp |= cbp_masks[get_bits1(gb)] << i;
        if(t == 2)
            cbp |= cbp_masks[2] << i;
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

/**
 * Decode 2x2 subblock of coefficients
 */
static inline void decode_subblock(DCTELEM *dst, int code, const int is_block2, GetBitContext *gb, VLC *vlc)
{
    int coeffs[4];

    coeffs[0] = modulo_three_table[code][0];
    coeffs[1] = modulo_three_table[code][1];
    coeffs[2] = modulo_three_table[code][2];
    coeffs[3] = modulo_three_table[code][3];
    decode_coeff(dst  , coeffs[0], 3, gb, vlc);
    if(!is_block2){
        decode_coeff(dst+1, coeffs[1], 2, gb, vlc);
        decode_coeff(dst+8, coeffs[2], 2, gb, vlc);
    }else{
        decode_coeff(dst+8, coeffs[1], 2, gb, vlc);
        decode_coeff(dst+1, coeffs[2], 2, gb, vlc);
    }
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

    code = get_vlc2(gb, rvlc->first_pattern[fc].table, 9, 2);

    pattern = code & 0x7;

    code >>= 3;
    decode_subblock(dst, code, 0, gb, &rvlc->coefficient);

    if(pattern & 4){
        code = get_vlc2(gb, rvlc->second_pattern[sc].table, 9, 2);
        decode_subblock(dst + 2, code, 0, gb, &rvlc->coefficient);
    }
    if(pattern & 2){ // Looks like coefficients 1 and 2 are swapped for this block
        code = get_vlc2(gb, rvlc->second_pattern[sc].table, 9, 2);
        decode_subblock(dst + 8*2, code, 1, gb, &rvlc->coefficient);
    }
    if(pattern & 1){
        code = get_vlc2(gb, rvlc->third_pattern[sc].table, 9, 2);
        decode_subblock(dst + 8*2+2, code, 0, gb, &rvlc->coefficient);
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
        for(j = !i; j < 4; j++)
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
static int get_dimension(GetBitContext *gb, const int *dim1, const int *dim2)
{
    int val, t;

    t = get_bits(gb, 3);
    val = dim1[t];
    if(!val && dim2)
        val = dim2[(t*2 | get_bits1(gb)) & 3];
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

/**
 * Select VLC set for decoding from current quantizer, modifier and frame type
 */
static inline RV40VLC* choose_vlc_set(int quant, int mod, int type)
{
    if(mod == 2){
        if(quant < 19) quant += 10;
        else if(quant < 26) quant += 5;
    }
    if(mod == 1)
        if(quant < 26) quant += 5;
    return type ? &inter_vlcs[rv40_quant_to_vlc_set[1][av_clip(quant, 0, 30)]]
                : &intra_vlcs[rv40_quant_to_vlc_set[0][av_clip(quant, 0, 30)]];
}

static int rv30_parse_slice_header(RV40DecContext *r, GetBitContext *gb, SliceInfo *si)
{
    int t, mb_bits;
    int w = r->s.width, h = r->s.height;
    int i, mb_size;

    memset(si, 0, sizeof(SliceInfo));
    get_bits(gb, 3);
    si->type = get_bits(gb, 2);
    if(si->type == 1) si->type = 0;
    if(get_bits1(gb))
        return -1;
    si->quant = get_bits(gb, 5);
    get_bits1(gb);
    t = get_bits(gb, 13);
    skip_bits(gb, r->rpr);
    si->vlc_set = 0;
    si->width  = w;
    si->height = h;
    mb_size = ((w + 15) >> 4) * ((h + 15) >> 4);
    for(i = 0; i < 5; i++)
        if(rv40_mb_max_sizes[i] > mb_size)
            break;
    mb_bits = rv40_mb_bits_sizes[i];
    si->start = get_bits(gb, mb_bits);
    get_bits1(gb);
    si->header_size = get_bits_count(gb);
    return 0;
}

static int rv40_parse_slice_header(RV40DecContext *r, GetBitContext *gb, SliceInfo *si)
{
    int t, mb_bits;
    int w = r->s.width, h = r->s.height;
    int i, mb_size;

    memset(si, 0, sizeof(SliceInfo));
    si->type = -1;
    if(get_bits1(gb))
        return -1;
    si->type = get_bits(gb, 2);
    if(si->type == 1) si->type = 0;
    si->quant = get_bits(gb, 5);
    if(get_bits(gb, 2))
        return -1;
    si->vlc_set = get_bits(gb, 2);
    get_bits1(gb);
    t = get_bits(gb, 13); /// ???
    if(!si->type || !get_bits1(gb))
        rv40_parse_picture_size(gb, &w, &h);
    si->width  = w;
    si->height = h;
    mb_size = ((w + 15) >> 4) * ((h + 15) >> 4);
    for(i = 0; i < 5; i++)
        if(rv40_mb_max_sizes[i] > mb_size)
            break;
    mb_bits = rv40_mb_bits_sizes[i];
    si->start = get_bits(gb, mb_bits);
    si->header_size = get_bits_count(gb);

    return 0;
}

static inline int get_omega(GetBitContext *gb);

/**
 * Decode 4x4 intra types array
 */
static int rv30_decode_intra_types(RV40DecContext *r, GetBitContext *gb, int *dst)
{
    int i, j, k;
    int A, B;
    int *ptr;
    int code;

    for(i = 0; i < 4; i++, dst += r->intra_types_stride){
        ptr = dst;
        for(j = 0; j < 4; j+= 2){
            code = (get_omega(gb) - 1) << 1;
            if(code >= 81*2){
                av_log(r->s.avctx, AV_LOG_ERROR, "Incorrect intra prediction code\n");
                return -1;
            }
            for(k = 0; k < 2; k++){
                A = ptr[-r->intra_types_stride] + 1;
                B = ptr[-1] + 1;
                *ptr++ = rv30_itype_from_context[A * 90 + B * 9 + rv30_itype_code[code + k]];
                if(ptr[-1] == 9){
                    av_log(r->s.avctx, AV_LOG_ERROR, "Incorrect intra prediction mode\n");
                    return -1;
                }
            }
        }
    }
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

/**
 * Decode quantizer difference and return modified quantizer
 */
static inline int rv40_decode_dquant(GetBitContext *gb, int quant)
{
    if(get_bits1(gb))
        return av_clip(quant + rv40_dquant_tab[quant * 2 + get_bits1(gb)], 0, 31);
    else
        return get_bits(gb, 5);
}

/**
 * Decode variable-length code constructed from variable-length codes
 * similar to Even-Rodeh and Elias Omega codes
 *
 * Code is constructed from bit chunks of even length (odd length means end of code)
 * and chunks are coded with variable-length codes too
 */
static inline int get_omega(GetBitContext *gb)
{
    int code = 1, t, tb;

    for(;;){
        t = get_vlc2(gb, mbinfo_vlc.table, MBINFO_BITS, 1);
        tb = t >> 5;
        code = (code << tb) | (t & 0xF);
        if(t & 0x10) break;
    }
    return code;
}

/**
 * Decode signed integer variable-length code constructed from variable-length codes
 * similar to Even-Rodeh and Elias Omega codes
 *
 * Code is constructed from bit chunks of even length (odd length means end of code)
 * and chunks are coded with variable-length codes too
 */
static inline int get_omega_signed(GetBitContext *gb)
{
    int code;

    code = get_omega(gb);
    if(code & 1)
        return -(code >> 1);
    else
        return code >> 1;
}

/**
 * Decode macroblock information
 */
static int rv30_decode_mb_info(RV40DecContext *r)
{
    static const int rv30_p_types[6] = { RV40_MB_SKIP, RV40_MB_P_16x16, RV40_MB_P_8x8, -1, RV40_MB_TYPE_INTRA, RV40_MB_TYPE_INTRA16x16 };
    static const int rv30_b_types[6] = { RV40_MB_SKIP, RV40_MB_B_INTERP, RV40_MB_B_FORWARD, RV40_MB_B_BACKWARD, RV40_MB_TYPE_INTRA, RV40_MB_TYPE_INTRA16x16 };
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int code;

    code = get_omega(gb) - 1;
    if(code > 11){
        av_log(s->avctx, AV_LOG_ERROR, "Incorrect MB type code\n");
        return -1;
    }
    if(code > 5){
        av_log(NULL,0, "dquant needed\n");
        code -= 6;
    }
    if(s->pict_type != B_TYPE)
        return rv30_p_types[code];
    else
        return rv30_b_types[code];
}

/**
 * Decode macroblock information
 */
static int rv40_decode_mb_info(RV40DecContext *r)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int q, i;
    int prev_type = 0;
    int mb_pos = s->mb_x + s->mb_y * s->mb_stride;
    int blocks[RV40_MB_TYPES];
    int count = 0;

    if(!r->s.mb_skip_run)
        r->s.mb_skip_run = get_omega(gb);

    if(--r->s.mb_skip_run)
         return RV40_MB_SKIP;

    memset(blocks, 0, sizeof(blocks));
    if(r->avail[0])
        blocks[r->mb_type[mb_pos - 1]]++;
    if(r->avail[1])
        blocks[r->mb_type[mb_pos - s->mb_stride]]++;
    if(r->avail[1] && r->avail[2])
        blocks[r->mb_type[mb_pos - s->mb_stride + 1]]++;
    if(r->avail[1] && r->avail[3])
        blocks[r->mb_type[mb_pos - s->mb_stride - 1]]++;

    for(i = 0; i < RV40_MB_TYPES; i++){
        if(blocks[i] > count){
            count = blocks[i];
            prev_type = i;
        }
    }
    if(s->pict_type == P_TYPE){
        if(prev_type == RV40_MB_SKIP) prev_type = RV40_MB_P_16x16;
        prev_type = block_num_to_ptype_vlc_num[prev_type];
        q = get_vlc2(gb, ptype_vlc[prev_type].table, PTYPE_VLC_BITS, 1);
        if(q < PBTYPE_ESCAPE)
            return q;
        q = get_vlc2(gb, ptype_vlc[prev_type].table, PTYPE_VLC_BITS, 1);
        av_log(NULL,0,"Dquant for P-frame\n");
    }else{
        prev_type = block_num_to_btype_vlc_num[prev_type];
        q = get_vlc2(gb, btype_vlc[prev_type].table, BTYPE_VLC_BITS, 1);
        if(q < PBTYPE_ESCAPE)
            return q;
        q = get_vlc2(gb, btype_vlc[prev_type].table, BTYPE_VLC_BITS, 1);
        av_log(NULL,0,"Dquant for B-frame\n");
    }
    return 0;
}

/** @} */ //bitstream functions

/**
 * @defgroup mv motion vector related code (prediction, reconstruction, motion compensation)
 * @{
 */

/** Macroblock partition width in 8x8 blocks */
static const uint8_t part_sizes_w[RV40_MB_TYPES] = { 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2 };

/** Macroblock partition height in 8x8 blocks */
static const uint8_t part_sizes_h[RV40_MB_TYPES] = { 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2 };

/**
 * Motion vectors prediction
 *
 * Motion prediction performed for the block by using median prediction of
 * motion vector from the left, top and right top blocks but in corener cases
 * some other vectors may be used instead
 */
static void rv40_pred_mv(RV40DecContext *r, int block_type, int subblock_no)
{
    MpegEncContext *s = &r->s;
    int mv_pos = s->mb_x * 2 + s->mb_y * 2 * s->b8_stride;
    int A[2], B[2], C[2];
    int no_A = 1, no_B = 1, no_C = 1;
    int i, j;
    int mx, my;

    memset(A, 0, sizeof(A));
    memset(B, 0, sizeof(B));
    memset(C, 0, sizeof(C));
    no_A = !r->avail[0];
    no_B = !r->avail[1];
    no_C = !r->avail[2];
    switch(block_type){
    case RV40_MB_P_16x16:
    case RV40_MB_P_MIX16x16:
        if(!no_C){
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][1];
        }
        break;
    case RV40_MB_P_8x8:
        mv_pos += (subblock_no & 1) + (subblock_no >> 1)*s->b8_stride;
        if(subblock_no & 1) no_A = 0;
        if(subblock_no & 2) no_B = 0;
        no_C |= (subblock_no == 3);
        if(subblock_no == 2) no_C = 0;
        if(!subblock_no) no_C = no_B;
        if(!no_C){
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+1][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+1][1];
        }
        if(subblock_no == 3){
            no_C = 0;
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride-1][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride-1][1];
        }
        break;
    case RV40_MB_P_16x8:
        mv_pos += subblock_no*s->b8_stride;
        no_B &= ~subblock_no;
        no_C |= subblock_no;
        if(!no_C){
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][1];
        }
        break;
    case RV40_MB_P_8x16:
        mv_pos += subblock_no;
        no_A &= ~subblock_no;
        if(!subblock_no) no_C = no_B;
        if(!no_C){
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+1][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+1][1];
        }
        break;
    default:
        no_A = no_B = no_C = 1;
    }
    if(!no_A){
        A[0] = s->current_picture_ptr->motion_val[0][mv_pos-1][0];
        A[1] = s->current_picture_ptr->motion_val[0][mv_pos-1][1];
    }
    if(!no_B){
        B[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride][0];
        B[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride][1];
    }else{
        B[0] = A[0];
        B[1] = A[1];
    }
    if(no_C){
        if(no_B || (no_A && !r->rv30)){
            C[0] = A[0];
            C[1] = A[1];
        }else{
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride-1][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride-1][1];
        }
    }
    mx = mid_pred(A[0], B[0], C[0]);
    my = mid_pred(A[1], B[1], C[1]);
    mx += r->dmv[subblock_no][0];
    my += r->dmv[subblock_no][1];
    for(j = 0; j < part_sizes_h[block_type]; j++){
        for(i = 0; i < part_sizes_w[block_type]; i++){
            s->current_picture_ptr->motion_val[0][mv_pos + i + j*s->b8_stride][0] = mx;
            s->current_picture_ptr->motion_val[0][mv_pos + i + j*s->b8_stride][1] = my;
        }
    }
}

/**
 * Predict motion vector for B-frame macroblock.
 */
static inline void rv40_pred_b_vector(int A[2], int B[2], int C[2], int no_A, int no_B, int no_C, int *mx, int *my)
{
    if(no_A + no_B + no_C){
        *mx = A[0] + B[0] + C[0];
        *my = A[1] + B[1] + C[1];
        if(no_A + no_B + no_C == 1){
            *mx /= 2;
            *my /= 2;
        }
    }else{
        *mx = mid_pred(A[0], B[0], C[0]);
        *my = mid_pred(A[1], B[1], C[1]);
    }
}

/**
 * Motion vector prediction for B-frames.
 */
static void rv40_pred_mv_b(RV40DecContext *r, int block_type)
{
    MpegEncContext *s = &r->s;
    int mb_pos = s->mb_x + s->mb_y * s->mb_stride;
    int mv_pos = s->mb_x * 2 + s->mb_y * 2 * s->b8_stride;
    int A[2][2], B[2][2], C[2][2];
    int no_A[2], no_B[2], no_C[2];
    int c_mv_pos;
    int mx[2], my[2];
    int i, j;

    memset(A, 0, sizeof(A));
    memset(B, 0, sizeof(B));
    memset(C, 0, sizeof(C));
    memset(mx, 0, sizeof(mx));
    memset(my, 0, sizeof(my));
    if(!r->avail[0])
        no_A[0] = no_A[1] = 1;
    else{
        no_A[0] = no_A[1] = 0;
        if(r->mb_type[mb_pos - 1] != RV40_MB_B_FORWARD  && r->mb_type[mb_pos - 1] != RV40_MB_B_DIRECT)
            no_A[0] = 1;
        if(r->mb_type[mb_pos - 1] != RV40_MB_B_BACKWARD && r->mb_type[mb_pos - 1] != RV40_MB_B_DIRECT)
            no_A[1] = 1;
        if(!no_A[0]){
            A[0][0] = s->current_picture_ptr->motion_val[0][mv_pos - 1][0];
            A[0][1] = s->current_picture_ptr->motion_val[0][mv_pos - 1][1];
        }
        if(!no_A[1]){
            A[1][0] = s->current_picture_ptr->motion_val[1][mv_pos - 1][0];
            A[1][1] = s->current_picture_ptr->motion_val[1][mv_pos - 1][1];
        }
    }
    if(!r->avail[1]){
        no_B[0] = no_B[1] = 1;
    }else{
        no_B[0] = no_B[1] = 0;
        if(r->mb_type[mb_pos - s->mb_stride] != RV40_MB_B_FORWARD  && r->mb_type[mb_pos - s->mb_stride] != RV40_MB_B_DIRECT)
            no_B[0] = 1;
        if(r->mb_type[mb_pos - s->mb_stride] != RV40_MB_B_BACKWARD && r->mb_type[mb_pos - s->mb_stride] != RV40_MB_B_DIRECT)
            no_B[1] = 1;
        if(!no_B[0]){
            B[0][0] = s->current_picture_ptr->motion_val[0][mv_pos - s->b8_stride][0];
            B[0][1] = s->current_picture_ptr->motion_val[0][mv_pos - s->b8_stride][1];
        }
        if(!no_B[1]){
            B[1][0] = s->current_picture_ptr->motion_val[1][mv_pos - s->b8_stride][0];
            B[1][1] = s->current_picture_ptr->motion_val[1][mv_pos - s->b8_stride][1];
        }
    }
    if(r->avail[2]){
        no_C[0] = no_C[1] = 0;
        if(r->mb_type[mb_pos - s->mb_stride + 1] != RV40_MB_B_FORWARD  && r->mb_type[mb_pos - s->mb_stride + 1] != RV40_MB_B_DIRECT)
            no_C[0] = 1;
        if(r->mb_type[mb_pos - s->mb_stride + 1] != RV40_MB_B_BACKWARD && r->mb_type[mb_pos - s->mb_stride + 1] != RV40_MB_B_DIRECT)
            no_C[1] = 1;
        c_mv_pos = mv_pos - s->b8_stride + 2;
    }else if(r->avail[3]){
        no_C[0] = no_C[1] = 0;
        if(r->mb_type[mb_pos - s->mb_stride - 1] != RV40_MB_B_FORWARD  && r->mb_type[mb_pos - s->mb_stride - 1] != RV40_MB_B_DIRECT)
            no_C[0] = 1;
        if(r->mb_type[mb_pos - s->mb_stride - 1] != RV40_MB_B_BACKWARD && r->mb_type[mb_pos - s->mb_stride - 1] != RV40_MB_B_DIRECT)
            no_C[1] = 1;
        c_mv_pos = mv_pos - s->b8_stride - 1;
    }else{
        no_C[0] = no_C[1] = 1;
        c_mv_pos = 0;
    }
    if(!no_C[0]){
        C[0][0] = s->current_picture_ptr->motion_val[0][c_mv_pos][0];
        C[0][1] = s->current_picture_ptr->motion_val[0][c_mv_pos][1];
    }
    if(!no_C[1]){
        C[1][0] = s->current_picture_ptr->motion_val[1][c_mv_pos][0];
        C[1][1] = s->current_picture_ptr->motion_val[1][c_mv_pos][1];
    }
    switch(block_type){
    case RV40_MB_B_FORWARD:
        rv40_pred_b_vector(A[0], B[0], C[0], no_A[0], no_B[0], no_C[0], &mx[0], &my[0]);
        r->dmv[1][0] = 0;
        r->dmv[1][1] = 0;
        break;
    case RV40_MB_B_BACKWARD:
        r->dmv[1][0] = r->dmv[0][0];
        r->dmv[1][1] = r->dmv[0][1];
        r->dmv[0][0] = 0;
        r->dmv[0][1] = 0;
        rv40_pred_b_vector(A[1], B[1], C[1], no_A[1], no_B[1], no_C[1], &mx[1], &my[1]);
        break;
    case RV40_MB_B_DIRECT:
        rv40_pred_b_vector(A[0], B[0], C[0], no_A[0], no_B[0], no_C[0], &mx[0], &my[0]);
        rv40_pred_b_vector(A[1], B[1], C[1], no_A[1], no_B[1], no_C[1], &mx[1], &my[1]);
        break;
    default:
        no_A[0] = no_A[1] = no_B[0] = no_B[1] = no_C[0] = no_C[1] = 1;
    }

    mx[0] += r->dmv[0][0];
    my[0] += r->dmv[0][1];
    mx[1] += r->dmv[1][0];
    my[1] += r->dmv[1][1];
    for(j = 0; j < 2; j++){
        for(i = 0; i < 2; i++){
            s->current_picture_ptr->motion_val[0][mv_pos + i + j*s->b8_stride][0] = mx[0];
            s->current_picture_ptr->motion_val[0][mv_pos + i + j*s->b8_stride][1] = my[0];
            s->current_picture_ptr->motion_val[1][mv_pos + i + j*s->b8_stride][0] = mx[1];
            s->current_picture_ptr->motion_val[1][mv_pos + i + j*s->b8_stride][1] = my[1];
        }
    }
}

/**
 * Generic motion compensation function - hopefully compiler will optimize it for each case
 *
 * @param r decoder context
 * @param block_type type of the current block
 * @param xoff horizontal offset from the start of the current block
 * @param yoff vertical offset from the start of the current block
 * @param mv_off offset to the motion vector information
 * @param width width of the current partition in 8x8 blocks
 * @param height height of the current partition in 8x8 blocks
 */
static inline void rv40_mc(RV40DecContext *r, const int block_type,
                          const int xoff, const int yoff, int mv_off,
                          const int width, const int height)
{
    MpegEncContext *s = &r->s;
    uint8_t *Y, *U, *V, *srcY, *srcU, *srcV;
    int dxy, mx, my, uvmx, uvmy, src_x, src_y, uvsrc_x, uvsrc_y;
    int mv_pos = s->mb_x * 2 + s->mb_y * 2 * s->b8_stride + mv_off;

    mx = s->current_picture_ptr->motion_val[0][mv_pos][0];
    my = s->current_picture_ptr->motion_val[0][mv_pos][1];
    srcY = s->last_picture_ptr->data[0];
    srcU = s->last_picture_ptr->data[1];
    srcV = s->last_picture_ptr->data[2];
    src_x = s->mb_x * 16 + xoff + (mx >> 2);
    src_y = s->mb_y * 16 + yoff + (my >> 2);
    uvsrc_x = s->mb_x * 8 + (xoff >> 1) + (mx >> 3);
    uvsrc_y = s->mb_y * 8 + (yoff >> 1) + (my >> 3);
    srcY += src_y * s->linesize + src_x;
    srcU += uvsrc_y * s->uvlinesize + uvsrc_x;
    srcV += uvsrc_y * s->uvlinesize + uvsrc_x;
    if(   (unsigned)(src_x - !!(mx&3)*2) > s->h_edge_pos - !!(mx&3)*2 - (width <<3) - 3
       || (unsigned)(src_y - !!(my&3)*2) > s->v_edge_pos - !!(my&3)*2 - (height<<3) - 3){
        uint8_t *uvbuf= s->edge_emu_buffer + 20 * s->linesize;

        srcY -= 2 + 2*s->linesize;
        ff_emulated_edge_mc(s->edge_emu_buffer, srcY, s->linesize, (width<<3)+4, (height<<3)+4,
                            src_x - 2, src_y - 2, s->h_edge_pos, s->v_edge_pos);
        srcY = s->edge_emu_buffer + 2 + 2*s->linesize;
        ff_emulated_edge_mc(uvbuf     , srcU, s->uvlinesize, (width<<2)+1, (height<<2)+1,
                            uvsrc_x, uvsrc_y, s->h_edge_pos >> 1, s->v_edge_pos >> 1);
        ff_emulated_edge_mc(uvbuf + 16, srcV, s->uvlinesize, (width<<2)+1, (height<<2)+1,
                            uvsrc_x, uvsrc_y, s->h_edge_pos >> 1, s->v_edge_pos >> 1);
        srcU = uvbuf;
        srcV = uvbuf + 16;
    }
    dxy = ((my & 3) << 2) | (mx & 3);
    uvmx = mx & 6;
    uvmy = my & 6;
    Y = s->dest[0] + xoff + yoff*s->linesize;
    U = s->dest[1] + (xoff>>1) + (yoff>>1)*s->uvlinesize;
    V = s->dest[2] + (xoff>>1) + (yoff>>1)*s->uvlinesize;
    if(block_type == RV40_MB_P_16x8){
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        Y    += 8;
        srcY += 8;
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        s->dsp.put_h264_chroma_pixels_tab[0]   (U, srcU, s->uvlinesize, 4, uvmx, uvmy);
        s->dsp.put_h264_chroma_pixels_tab[0]   (V, srcV, s->uvlinesize, 4, uvmx, uvmy);
    }else if(block_type == RV40_MB_P_8x16){
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        Y    += 8 * s->linesize;
        srcY += 8 * s->linesize;
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        s->dsp.put_h264_chroma_pixels_tab[1]   (U, srcU, s->uvlinesize, 8, uvmx, uvmy);
        s->dsp.put_h264_chroma_pixels_tab[1]   (V, srcV, s->uvlinesize, 8, uvmx, uvmy);
    }else if(block_type == RV40_MB_P_8x8){
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        s->dsp.put_h264_chroma_pixels_tab[1]   (U, srcU, s->uvlinesize, 4, uvmx, uvmy);
        s->dsp.put_h264_chroma_pixels_tab[1]   (V, srcV, s->uvlinesize, 4, uvmx, uvmy);
    }else{
        s->dsp.put_h264_qpel_pixels_tab[0][dxy](Y, srcY, s->linesize);
        s->dsp.put_h264_chroma_pixels_tab[0]   (U, srcU, s->uvlinesize, 8, uvmx, uvmy);
        s->dsp.put_h264_chroma_pixels_tab[0]   (V, srcV, s->uvlinesize, 8, uvmx, uvmy);
    }
}

/**
 * B-frame specific motion compensation function
 *
 * @param r decoder context
 * @param block_type type of the current block
 */
static inline void rv40_mc_b(RV40DecContext *r, const int block_type)
{
    MpegEncContext *s = &r->s;
    uint8_t *srcY, *srcU, *srcV;
    int dxy, mx, my, uvmx, uvmy, src_x, src_y, uvsrc_x, uvsrc_y;
    int mv_pos = s->mb_x * 2 + s->mb_y * 2 * s->b8_stride;

    if(block_type != RV40_MB_B_BACKWARD){
        mx = s->current_picture_ptr->motion_val[0][mv_pos][0];
        my = s->current_picture_ptr->motion_val[0][mv_pos][1];
        srcY = s->last_picture_ptr->data[0];
        srcU = s->last_picture_ptr->data[1];
        srcV = s->last_picture_ptr->data[2];
    }else{
        mx = s->current_picture_ptr->motion_val[1][mv_pos][0];
        my = s->current_picture_ptr->motion_val[1][mv_pos][1];
        srcY = s->next_picture_ptr->data[0];
        srcU = s->next_picture_ptr->data[1];
        srcV = s->next_picture_ptr->data[2];
    }
    if(block_type == RV40_MB_B_INTERP){
        mx += (s->next_picture_ptr->motion_val[0][mv_pos][0] + 1) >> 1;
        my += (s->next_picture_ptr->motion_val[0][mv_pos][1] + 1) >> 1;
    }
    src_x = s->mb_x * 16 + (mx >> 2);
    src_y = s->mb_y * 16 + (my >> 2);
    uvsrc_x = s->mb_x * 8 + (mx >> 3);
    uvsrc_y = s->mb_y * 8 + (my >> 3);
    srcY += src_y * s->linesize + src_x;
    srcU += uvsrc_y * s->uvlinesize + uvsrc_x;
    srcV += uvsrc_y * s->uvlinesize + uvsrc_x;
    if(   (unsigned)(src_x - !!(mx&3)*2) > s->h_edge_pos - !!(mx&3)*2 - 16 - 3
       || (unsigned)(src_y - !!(my&3)*2) > s->v_edge_pos - !!(my&3)*2 - 16 - 3){
        uint8_t *uvbuf= s->edge_emu_buffer + 20 * s->linesize;

        srcY -= 2 + 2*s->linesize;
        ff_emulated_edge_mc(s->edge_emu_buffer, srcY, s->linesize, 16+4, 16+4,
                            src_x - 2, src_y - 2, s->h_edge_pos, s->v_edge_pos);
        srcY = s->edge_emu_buffer + 2 + 2*s->linesize;
        ff_emulated_edge_mc(uvbuf     , srcU, s->uvlinesize, 8+1, 8+1,
                            uvsrc_x, uvsrc_y, s->h_edge_pos >> 1, s->v_edge_pos >> 1);
        ff_emulated_edge_mc(uvbuf + 16, srcV, s->uvlinesize, 8+1, 8+1,
                            uvsrc_x, uvsrc_y, s->h_edge_pos >> 1, s->v_edge_pos >> 1);
        srcU = uvbuf;
        srcV = uvbuf + 16;
    }
    dxy = ((my & 3) << 2) | (mx & 3);
    uvmx = mx & 6;
    uvmy = my & 6;
    s->dsp.put_h264_qpel_pixels_tab[0][dxy](s->dest[0], srcY, s->linesize);
    s->dsp.put_h264_chroma_pixels_tab[0]   (s->dest[1], srcU, s->uvlinesize, 8, uvmx, uvmy);
    s->dsp.put_h264_chroma_pixels_tab[0]   (s->dest[2], srcV, s->uvlinesize, 8, uvmx, uvmy);
}

/**
 * B-frame specific motion compensation function - for direct/interpolated blocks
 *
 * @param r decoder context
 * @param block_type type of the current block
 */
static inline void rv40_mc_b_interp(RV40DecContext *r, const int block_type)
{
    MpegEncContext *s = &r->s;
    uint8_t *srcY, *srcU, *srcV;
    int dxy, mx, my, uvmx, uvmy, src_x, src_y, uvsrc_x, uvsrc_y;
    int mv_pos = s->mb_x * 2 + s->mb_y * 2 * s->b8_stride;

    mx = s->current_picture_ptr->motion_val[1][mv_pos][0];
    my = s->current_picture_ptr->motion_val[1][mv_pos][1];
    if(block_type == RV40_MB_B_INTERP){
        mx -= s->next_picture_ptr->motion_val[0][mv_pos][0] >> 1;
        my -= s->next_picture_ptr->motion_val[0][mv_pos][1] >> 1;
    }
    srcY = s->next_picture_ptr->data[0];
    srcU = s->next_picture_ptr->data[1];
    srcV = s->next_picture_ptr->data[2];

    src_x = s->mb_x * 16 + (mx >> 2);
    src_y = s->mb_y * 16 + (my >> 2);
    uvsrc_x = s->mb_x * 8 + (mx >> 3);
    uvsrc_y = s->mb_y * 8 + (my >> 3);
    srcY += src_y * s->linesize + src_x;
    srcU += uvsrc_y * s->uvlinesize + uvsrc_x;
    srcV += uvsrc_y * s->uvlinesize + uvsrc_x;
    if(   (unsigned)(src_x - !!(mx&3)*2) > s->h_edge_pos - !!(mx&3)*2 - 16 - 3
       || (unsigned)(src_y - !!(my&3)*2) > s->v_edge_pos - !!(my&3)*2 - 16 - 3){
        uint8_t *uvbuf= s->edge_emu_buffer + 20 * s->linesize;

        srcY -= 2 + 2*s->linesize;
        ff_emulated_edge_mc(s->edge_emu_buffer, srcY, s->linesize, 16+4, 16+4,
                            src_x - 2, src_y - 2, s->h_edge_pos, s->v_edge_pos);
        srcY = s->edge_emu_buffer + 2 + 2*s->linesize;
        ff_emulated_edge_mc(uvbuf     , srcU, s->uvlinesize, 8+1, 8+1,
                            uvsrc_x, uvsrc_y, s->h_edge_pos >> 1, s->v_edge_pos >> 1);
        ff_emulated_edge_mc(uvbuf + 16, srcV, s->uvlinesize, 8+1, 8+1,
                            uvsrc_x, uvsrc_y, s->h_edge_pos >> 1, s->v_edge_pos >> 1);
        srcU = uvbuf;
        srcV = uvbuf + 16;
    }
    dxy = ((my & 3) << 2) | (mx & 3);
    uvmx = mx & 6;
    uvmy = my & 6;
    s->dsp.avg_h264_qpel_pixels_tab[0][dxy](s->dest[0], srcY, s->linesize);
    s->dsp.avg_h264_chroma_pixels_tab[0]   (s->dest[1], srcU, s->uvlinesize, 8, uvmx, uvmy);
    s->dsp.avg_h264_chroma_pixels_tab[0]   (s->dest[2], srcV, s->uvlinesize, 8, uvmx, uvmy);
}

/**
 * Decode motion vector differences
 * and perform motion vector reconstruction and motion compensation.
 */
static int rv40_decode_mv(RV40DecContext *r, int block_type)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int i, j;

    switch(block_type){
    case RV40_MB_TYPE_INTRA:
    case RV40_MB_TYPE_INTRA16x16:
        for(j = 0; j < 2; j++){
            for(i = 0; i < 2; i++){
                s->current_picture_ptr->motion_val[0][s->mb_x * 2 + s->mb_y * 2 * s->b8_stride + i + j*s->b8_stride][0] = 0;
                s->current_picture_ptr->motion_val[0][s->mb_x * 2 + s->mb_y * 2 * s->b8_stride + i + j*s->b8_stride][1] = 0;
            }
        }
        return 0;
    case RV40_MB_SKIP:
        r->dmv[0][0] = 0;
        r->dmv[0][1] = 0;
        if(s->pict_type == P_TYPE){
            rv40_pred_mv(r, block_type, 0);
            rv40_mc(r, block_type, 0, 0, 0, 2, 2);
            break;
        }
    case RV40_MB_B_INTERP:
        r->dmv[0][0] = 0;
        r->dmv[0][1] = 0;
        r->dmv[1][0] = 0;
        r->dmv[1][1] = 0;
        rv40_pred_mv_b  (r, RV40_MB_B_INTERP);
        rv40_mc_b       (r, RV40_MB_B_INTERP);
        rv40_mc_b_interp(r, RV40_MB_B_INTERP);
        break;
    case RV40_MB_P_16x16:
    case RV40_MB_P_MIX16x16:
        r->dmv[0][0] = get_omega_signed(gb);
        r->dmv[0][1] = get_omega_signed(gb);
        rv40_pred_mv(r, block_type, 0);
        rv40_mc(r, block_type, 0, 0, 0, 2, 2);
        break;
    case RV40_MB_B_FORWARD:
    case RV40_MB_B_BACKWARD:
        r->dmv[0][0] = get_omega_signed(gb);
        r->dmv[0][1] = get_omega_signed(gb);
        rv40_pred_mv_b  (r, block_type);
        rv40_mc_b       (r, block_type);
        break;
    case RV40_MB_P_16x8:
    case RV40_MB_P_8x16:
    case RV40_MB_B_DIRECT:
        r->dmv[0][0] = get_omega_signed(gb);
        r->dmv[0][1] = get_omega_signed(gb);
        r->dmv[1][0] = get_omega_signed(gb);
        r->dmv[1][1] = get_omega_signed(gb);
        rv40_pred_mv(r, block_type, 0);
        rv40_pred_mv(r, block_type, 1);
        if(block_type == RV40_MB_P_16x8){
            rv40_mc(r, block_type, 0, 0, 0,            2, 1);
            rv40_mc(r, block_type, 0, 8, s->b8_stride, 2, 1);
        }
        if(block_type == RV40_MB_P_8x16){
            rv40_mc(r, block_type, 0, 0, 0, 1, 2);
            rv40_mc(r, block_type, 8, 0, 1, 1, 2);
        }
        if(block_type == RV40_MB_B_DIRECT){
            rv40_pred_mv_b  (r, block_type);
            rv40_mc_b       (r, block_type);
            rv40_mc_b_interp(r, block_type);
        }
        break;
    case RV40_MB_P_8x8:
        for(i=0;i< 4;i++){
            r->dmv[i][0] = get_omega_signed(gb);
            r->dmv[i][1] = get_omega_signed(gb);
            rv40_pred_mv(r, block_type, i);
            rv40_mc(r, block_type, (i&1)<<3, (i&2)<<2, (i&1)+(i>>1)*s->b8_stride, 1, 1);
        }
        break;
    }

    return 0;
}
/** @} */ // mv group

/**
 * @defgroup recons Macroblock reconstruction functions
 * @{
 */
/** Mapping of RV40 intra prediction types to standard H.264 types */
static const int ittrans[9] = {
 DC_PRED, VERT_PRED, HOR_PRED, DIAG_DOWN_RIGHT_PRED, DIAG_DOWN_LEFT_PRED,
 VERT_RIGHT_PRED, VERT_LEFT_PRED, HOR_UP_PRED, HOR_DOWN_PRED,
};

/** Mapping of RV40 intra 16x16 prediction types to standard H.264 types */
static const int ittrans16[4] = {
 DC_PRED8x8, VERT_PRED8x8, HOR_PRED8x8, PLANE_PRED8x8,
};

/**
 * Perform 4x4 intra prediction
 */
static void rv40_pred_4x4_block(RV40DecContext *r, uint8_t *dst, int stride, int itype, int no_up, int no_left, int no_down, int no_right)
{
    uint8_t *prev = dst - stride + 4;
    uint32_t topleft;

    if(no_up && no_left)
        itype = DC_128_PRED;
    else if(no_up){
        if(itype == VERT_PRED) itype = HOR_PRED;
        if(itype == DC_PRED)   itype = LEFT_DC_PRED;
    }else if(no_left){
        if(itype == HOR_PRED)  itype = VERT_PRED;
        if(itype == DC_PRED)   itype = TOP_DC_PRED;
        if(itype == DIAG_DOWN_LEFT_PRED) itype = DIAG_DOWN_LEFT_PRED_RV40_NODOWN;
    }
    if(no_down){
        if(itype == DIAG_DOWN_LEFT_PRED) itype = DIAG_DOWN_LEFT_PRED_RV40_NODOWN;
        if(itype == HOR_UP_PRED) itype = HOR_UP_PRED_RV40_NODOWN;
    }
    if(no_right && !no_up){
        topleft = dst[-stride + 3] * 0x01010101;
        prev = &topleft;
    }
    r->h.pred4x4[itype](dst, prev, stride);
}

/** add_pixels_clamped for 4x4 block */
static void rv40_add_4x4_block(uint8_t *dst, int stride, DCTELEM block[64], int off)
{
    int x, y;
    for(y = 0; y < 4; y++)
        for(x = 0; x < 4; x++)
            dst[x + y*stride] = av_clip_uint8(dst[x + y*stride] + block[off + x+y*8]);
}

static void rv40_output_macroblock(RV40DecContext *r, int *intra_types, int cbp, int is16)
{
    MpegEncContext *s = &r->s;
    DSPContext *dsp = &s->dsp;
    int i, j;
    uint8_t *Y, *YY, *U, *V;
    int no_up, no_left, no_topright, itype;

    no_up = !r->avail[1];
    Y = s->dest[0];
    U = s->dest[1];
    V = s->dest[2];
    if(!is16){
        for(j = 0; j < 4; j++){
            no_left = !r->avail[0];
            YY = Y;
            for(i = 0; i < 4; i++, cbp >>= 1, YY += 4){
                no_topright = no_up || (i==3 && j) || (i==3 && !j && (s->mb_x-1) == s->mb_width);
                rv40_pred_4x4_block(r, YY, s->linesize, ittrans[intra_types[i]], no_up, no_left, i || (j==3), no_topright);
                no_left = 0;
                if(!(cbp & 1)) continue;
                rv40_add_4x4_block(YY, s->linesize, s->block[(i>>1)+(j&2)], (i&1)*4+(j&1)*32);
            }
            no_up = 0;
            Y += s->linesize * 4;
            intra_types += r->intra_types_stride;
        }
        intra_types -= r->intra_types_stride * 4;
        no_up = !r->avail[1];
        for(j = 0; j < 2; j++){
            no_left = !r->avail[0];
            for(i = 0; i < 2; i++, cbp >>= 1, no_left = 0){
                no_topright = no_up || (i && j) || (i && !j && (s->mb_x-1) == s->mb_width);
                rv40_pred_4x4_block(r, U + i*4 + j*4*s->uvlinesize, s->uvlinesize, ittrans[intra_types[i*2+j*2*r->intra_types_stride]], no_up, no_left, i || j, no_topright);
                rv40_pred_4x4_block(r, V + i*4 + j*4*s->uvlinesize, s->uvlinesize, ittrans[intra_types[i*2+j*2*r->intra_types_stride]], no_up, no_left, i || j, no_topright);
                if(cbp & 0x01)
                    rv40_add_4x4_block(U + i*4 + j*4*s->uvlinesize, s->uvlinesize, s->block[4], i*4+j*32);
                if(cbp & 0x10)
                    rv40_add_4x4_block(V + i*4 + j*4*s->uvlinesize, s->uvlinesize, s->block[5], i*4+j*32);
            }
            no_up = 0;
        }
    }else{
        no_left = !r->avail[0];
        itype = ittrans16[intra_types[0]];
        if(no_up && no_left)
            itype = DC_128_PRED8x8;
        else if(no_up){
            if(itype == PLANE_PRED8x8)itype = HOR_PRED8x8;
            if(itype == VERT_PRED8x8) itype = HOR_PRED8x8;
            if(itype == DC_PRED8x8)   itype = LEFT_DC_PRED8x8;
        }else if(no_left){
            if(itype == PLANE_PRED8x8)itype = VERT_PRED8x8;
            if(itype == HOR_PRED8x8)  itype = VERT_PRED8x8;
            if(itype == DC_PRED8x8)   itype = TOP_DC_PRED8x8;
        }
        r->h.pred16x16[itype](Y, s->linesize);
        dsp->add_pixels_clamped(s->block[0], Y, s->current_picture.linesize[0]);
        dsp->add_pixels_clamped(s->block[1], Y + 8, s->current_picture.linesize[0]);
        Y += s->current_picture.linesize[0] * 8;
        dsp->add_pixels_clamped(s->block[2], Y, s->current_picture.linesize[0]);
        dsp->add_pixels_clamped(s->block[3], Y + 8, s->current_picture.linesize[0]);

        itype = ittrans16[intra_types[0]];
        if(itype == PLANE_PRED8x8) itype = DC_PRED8x8;
        if(no_up && no_left)
            itype = DC_128_PRED8x8;
        else if(no_up){
            if(itype == VERT_PRED8x8) itype = HOR_PRED8x8;
            if(itype == DC_PRED8x8)   itype = LEFT_DC_PRED8x8;
        }else if(no_left){
            if(itype == HOR_PRED8x8)  itype = VERT_PRED8x8;
            if(itype == DC_PRED8x8)   itype = TOP_DC_PRED8x8;
        }
        r->h.pred8x8[itype](U, s->uvlinesize);
        dsp->add_pixels_clamped(s->block[4], U, s->uvlinesize);
        r->h.pred8x8[itype](V, s->uvlinesize);
        dsp->add_pixels_clamped(s->block[5], V, s->uvlinesize);
    }
}

/** @} */ // recons group

/**
 * @addtogroup bitstream
 * Decode macroblock header and return CBP in case of success, -1 otherwise.
 */
static int rv40_decode_mb_header(RV40DecContext *r, int *intra_types)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int mb_pos = s->mb_x + s->mb_y * s->mb_stride;
    int i, t;

    if(!r->si.type && !r->rv30){
        r->is16 = 0;
        switch(decode210(gb)){
        case 0: // 16x16 block
            r->is16 = 1;
            break;
        case 1:
            break;
        case 2:
            av_log(NULL,0,"Need DQUANT\n");
            // q = decode_dquant(gb);
            break;
        }
        s->current_picture_ptr->mb_type[mb_pos] = r->is16 ? MB_TYPE_INTRA16x16 : MB_TYPE_INTRA;
        r->block_type = r->is16 ? RV40_MB_TYPE_INTRA16x16 : RV40_MB_TYPE_INTRA;
    }else if(!r->si.type && r->rv30){
        r->is16 = get_bits1(gb);
        s->current_picture_ptr->mb_type[mb_pos] = r->is16 ? MB_TYPE_INTRA16x16 : MB_TYPE_INTRA;
        r->block_type = r->is16 ? RV40_MB_TYPE_INTRA16x16 : RV40_MB_TYPE_INTRA;
    }else{
        r->block_type = r->rv30 ? rv30_decode_mb_info(r) : rv40_decode_mb_info(r);
        if(r->block_type == -1)
            return -1;
        s->current_picture_ptr->mb_type[mb_pos] = rv40_mb_type_to_lavc[r->block_type];
        r->mb_type[mb_pos] = r->block_type;
        if(s->pict_type == P_TYPE && r->block_type == RV40_MB_SKIP)
            r->mb_type[mb_pos] = RV40_MB_P_16x16;
        if(s->pict_type == B_TYPE && r->block_type == RV40_MB_SKIP)
            r->mb_type[mb_pos] = RV40_MB_B_INTERP;
        r->is16 = !!IS_INTRA16x16(s->current_picture_ptr->mb_type[mb_pos]);
        rv40_decode_mv(r, r->block_type);
        if(r->block_type == RV40_MB_SKIP){
            for(i = 0; i < 16; i++)
                intra_types[(i & 3) + (i>>2) * r->intra_types_stride] = 0;
            return 0;
        }
        r->chroma_vlc = 1;
        r->luma_vlc   = 0;
    }
    if(IS_INTRA(s->current_picture_ptr->mb_type[mb_pos])){
        if(!r->is16){
            if(r->rv30){
                if(rv30_decode_intra_types(r, gb, intra_types) < 0)
                    return -1;
            }else{
                if(rv40_decode_intra_types(r, gb, intra_types) < 0)
                    return -1;
            }
            r->chroma_vlc = 0;
            r->luma_vlc   = 1;
        }else{
            t = get_bits(gb, 2);
            for(i = 0; i < 16; i++)
                intra_types[(i & 3) + (i>>2) * r->intra_types_stride] = t;
            r->chroma_vlc = 0;
            r->luma_vlc   = 2;
        }
        r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 0);
    }else{
        for(i = 0; i < 16; i++)
            intra_types[(i & 3) + (i>>2) * r->intra_types_stride] = 0;
        r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 1);
        if(r->mb_type[mb_pos] == RV40_MB_P_MIX16x16){
            r->is16 = 1;
            r->chroma_vlc = 1;
            r->luma_vlc   = 2;
            r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 0);
        }
    }
    return rv40_decode_cbp(gb, r->cur_vlcs, r->is16);
}

/**
 * @addtogroup recons
 * @{
 */
/** Mask for retrieving all bits in coded block pattern
 * corresponding to one 8x8 block.
 */
#define LUMA_CBP_BLOCK_MASK 0x303

#define U_CBP_MASK 0x0F0000
#define V_CBP_MASK 0xF00000


static void rv40_apply_differences(RV40DecContext *r, int cbp)
{
    static const int shifts[4] = { 0, 2, 8, 10 };
    MpegEncContext *s = &r->s;
    int i;

    for(i = 0; i < 4; i++)
        if(cbp & (LUMA_CBP_BLOCK_MASK << shifts[i]))
            s->dsp.add_pixels_clamped(s->block[i], s->dest[0] + (i & 1)*8 + (i&2)*4*s->linesize, s->linesize);
    if(cbp & U_CBP_MASK)
        s->dsp.add_pixels_clamped(s->block[4], s->dest[1], s->uvlinesize);
    if(cbp & V_CBP_MASK)
        s->dsp.add_pixels_clamped(s->block[5], s->dest[2], s->uvlinesize);
}

static int rv40_decode_macroblock(RV40DecContext *r, int *intra_types)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int cbp, cbp2;
    int i, blknum, blkoff;
    DCTELEM block16[64];
    int luma_dc_quant;

    // calculate which neighbours are available
    memset(r->avail, 0, sizeof(r->avail));
    if(s->mb_x && !(s->first_slice_line && s->mb_x == s->resync_mb_x))
        r->avail[0] = 1;
    if(!s->first_slice_line)
        r->avail[1] = 1;
    if((s->mb_x+1) < s->mb_width && (!s->first_slice_line || (s->first_slice_line && (s->mb_x+1) == s->resync_mb_x)))
        r->avail[2] = 1;
    if(s->mb_x && !s->first_slice_line && !((s->mb_y-1)==s->resync_mb_y && s->mb_x == s->resync_mb_x))
        r->avail[3] = 1;

    s->qscale = r->si.quant;
    cbp = cbp2 = rv40_decode_mb_header(r, intra_types);

    if(cbp == -1)
        return -1;

    luma_dc_quant = r->rv30 ? rv30_luma_dc_quant[s->qscale] : rv40_luma_quant[r->si.type>>1][s->qscale];
    if(r->is16){
        memset(block16, 0, sizeof(block16));
        rv40_decode_block(block16, gb, r->cur_vlcs, 3, 0);
        rv40_dequant4x4_16x16(block16, 0, rv40_qscale_tab[luma_dc_quant],rv40_qscale_tab[s->qscale]);
        rv40_intra_inv_transform_noround(block16, 0);
    }

    for(i = 0; i < 16; i++, cbp >>= 1){
        if(!r->is16 && !(cbp & 1)) continue;
        blknum = ((i & 2) >> 1) + ((i & 8) >> 2);
        blkoff = ((i & 1) << 2) + ((i & 4) << 3);
        if(cbp & 1)
            rv40_decode_block(s->block[blknum] + blkoff, gb, r->cur_vlcs, r->luma_vlc, 0);
        if((cbp & 1) || r->is16){
            rv40_dequant4x4(s->block[blknum], blkoff, rv40_qscale_tab[luma_dc_quant],rv40_qscale_tab[s->qscale]);
            if(r->is16) //FIXME: optimize
                s->block[blknum][blkoff] = block16[(i & 3) | ((i & 0xC) << 1)];
            rv40_intra_inv_transform(s->block[blknum], blkoff);
        }
    }
    if(r->block_type == RV40_MB_P_MIX16x16)
        r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 1);
    for(; i < 24; i++, cbp >>= 1){
        if(!(cbp & 1)) continue;
        blknum = ((i & 4) >> 2) + 4;
        blkoff = ((i & 1) << 2) + ((i & 2) << 4);
        rv40_decode_block(s->block[blknum] + blkoff, gb, r->cur_vlcs, r->chroma_vlc, 1);
        rv40_dequant4x4(s->block[blknum], blkoff, rv40_qscale_tab[rv40_chroma_quant[1][s->qscale]],rv40_qscale_tab[rv40_chroma_quant[0][s->qscale]]);
        rv40_intra_inv_transform(s->block[blknum], blkoff);
    }
    if(IS_INTRA(s->current_picture_ptr->mb_type[s->mb_x + s->mb_y*s->mb_stride]))
        rv40_output_macroblock(r, intra_types, cbp2, r->is16);
    else
        rv40_apply_differences(r, cbp2);

    return 0;
}

static int check_slice_end(RV40DecContext *r, MpegEncContext *s)
{
    int bits;
    if(r->s.mb_skip_run > 1)
        return 0;
    if(s->mb_y >= s->mb_height)
        return 1;
    bits = r->bits - get_bits_count(&s->gb);
    if(bits < 0 || (bits < 8 && !show_bits(&s->gb, bits)))
        return 1;
    return 0;
}

static inline int slice_compare(SliceInfo *si1, SliceInfo *si2)
{
    return si1->type   != si2->type ||
           si1->start  >= si2->start ||
           si1->width  != si2->width ||
           si1->height != si2->height;
}

static int rv40_decode_slice(RV40DecContext *r, int size, int end, int *last)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int mb_pos;
    int res;
    int old_mb_x = r->ssi.mb_x, old_mb_y = r->ssi.mb_y;
    *last = 1;

    init_get_bits(&r->s.gb, r->slice_data, r->si.size);
    res = r->rv30 ? rv30_parse_slice_header(r, gb, &r->si) : rv40_parse_slice_header(r, gb, &r->si);
    if((res < 0 && !s->current_picture_ptr) || (r->prev_si.type == -1 && r->si.start)){
        av_log(s->avctx, AV_LOG_ERROR, "Incorrect or unknown slice header\n");
        *last = 0;
        return -1;
    }
    if(res < 0 || (r->prev_si.type != -1 && slice_compare(&r->prev_si, &r->si))){
        r->ssi.data = av_realloc(r->ssi.data, r->ssi.data_size + (size>>3));
        memcpy(r->ssi.data + r->ssi.data_size, r->slice_data, size >> 3);
        r->ssi.data_size += size >> 3; // XXX: overflow check?
        size = r->ssi.data_size * 8;
        init_get_bits(&r->s.gb, r->ssi.data, r->ssi.data_size * 8);
        r->si = r->prev_si;
        skip_bits(gb, r->ssi.bits_used);
        s->mb_x = r->ssi.mb_x;
        s->mb_y = r->ssi.mb_y;
        r->si.start = s->mb_x + s->mb_y * s->mb_width;
        r->truncated = 1;
    }

    if ((s->mb_x == 0 && s->mb_y == 0) || s->current_picture_ptr==NULL) {
        if(s->width != r->si.width || s->height != r->si.height /*&& avcodec_check_dimensions(s->avctx, r->si.width, r->si.height) >= 0 */){
            av_log(s->avctx, AV_LOG_DEBUG, "Changing dimensions to %dx%d\n", r->si.width,r->si.height);
            MPV_common_end(s);
            s->width  = r->si.width;
            s->height = r->si.height;
            if(MPV_common_init(s) < 0)
                return -1;
        }
        s->pict_type = r->si.type ? r->si.type : I_TYPE;
        if(MPV_frame_start(s, s->avctx) < 0)
            return -1;
        ff_er_frame_start(s);
        s->current_picture_ptr = &s->current_picture;
        s->mb_x = s->mb_y = 0;
        r->truncated = 0;
    }

    r->si.size = size;
    r->si.end = end;
    s->qscale = r->si.quant;
    r->bits = r->si.size;
    r->block_start = r->si.start;
    s->mb_num_left = r->si.end - r->si.start;
    r->s.mb_skip_run = 0;

    r->prev_si = r->si;

    mb_pos = s->mb_x + s->mb_y * s->mb_width;
    if(!r->truncated && r->block_start != mb_pos){
        av_log(s->avctx, AV_LOG_ERROR, "Slice indicates MB offset %d, got %d\n", r->block_start, mb_pos);
        s->mb_x = r->block_start % s->mb_width;
        s->mb_y = r->block_start / s->mb_width;
    }
    if(!r->truncated){
        memset(r->intra_types_hist, -1, r->intra_types_stride * 4 * 2 * sizeof(int));
        s->first_slice_line = 1;
        s->resync_mb_x= s->mb_x;
        s->resync_mb_y= s->mb_y;
    }
    ff_init_block_index(s);
    while(!check_slice_end(r, s) && s->mb_num_left-- && s->mb_y < s->mb_height) {
        ff_update_block_index(s);
        s->dsp.clear_blocks(s->block[0]);

        /* save information about decoded position in case of truncated slice */
        if(r->bits > get_bits_count(gb)){
            r->ssi.bits_used = get_bits_count(gb);
            r->ssi.mb_x = s->mb_x;
            r->ssi.mb_y = s->mb_y;
        }

        if(rv40_decode_macroblock(r, r->intra_types + (s->mb_x + 1) * 4) < 0)
            break;
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
    if(!r->truncated)
        ff_er_add_slice(s, s->resync_mb_x, s->resync_mb_y, s->mb_x-1, s->mb_y, AC_END|DC_END|MV_END);
    else // add only additionally decoded blocks
        ff_er_add_slice(s, old_mb_x+1, old_mb_y, s->mb_x-1, s->mb_y, AC_END|DC_END|MV_END);
    r->truncated = 0;
    *last = 0;
    if(s->mb_y >= s->mb_height)
        *last = 1;
    if(r->bits > get_bits_count(gb) && show_bits(gb, r->bits-get_bits_count(gb)))
        *last = 1;

    r->ssi.data = av_realloc(r->ssi.data, r->bits >> 3);
    r->ssi.data_size = r->bits >> 3;
    memcpy(r->ssi.data, r->slice_data, r->bits >> 3);
    return 0;
}

#define CLIP_SYMM(a, b) av_clip(a, -(b), b)
/**
 * Weaker deblocking
 */
static inline void rv40_weak_loop_filter(uint8_t *src, const int step,
                            const int flag0, const int flag1, const int mult,
                            const int lim0, const int lim1, const int lim2, const int thr1,
                            const int S0, const int S1, const int S2, const int S3)
{
    uint8_t *cm = ff_cropTbl + MAX_NEG_CROP;
    int t, diff;

    t = src[0*step] - src[-1*step];
    if(!t) return;
    t = (mult * FFABS(t)) >> 7;
    if(t > 3) return;
    if(flag0 && flag1 && t > 2) return;
    t = src[-1*step] - src[0*step];

    if(flag0 && flag1)
        diff = (src[-2*step] - src[1*step] + t*4 + 4) >> 3;
    else
        diff = (t + 1) >> 1;
    diff = CLIP_SYMM(diff, lim2);
    src[-1*step] = cm[src[-1*step] + diff];
    src[ 0*step] = cm[src[ 0*step] - diff];
    if(FFABS(S1) <= thr1 && flag0){
        t = (S0 + S1 - diff) >> 1;
        src[-2*step] = cm[src[-2*step] - CLIP_SYMM(t, lim1)];
    }
    if(FFABS(S3) <= thr1 && flag1){
        t = (S2 + S3 + diff) >> 1;
        src[ 1*step] = cm[src[ 1*step] - CLIP_SYMM(t, lim0)];
    }
}

/**
 * This macro is used for calculating 25*x0+26*x1+26*x2+26*x3+25*x4
 * or 25*x0+26*x1+51*x2+26*x3
 * parameter  sub - index of the value with coefficient = 25
 * parameter last - index of the value with coefficient 25 or 51
 */
#define RV40_STRONG_FILTER(src, step, start, last, sub) \
     26*(src[start*step] + src[(start+1)*step] + src[(start+2)*step] + src[(start+3)*step] + src[last*step]) - src[last*step] - src[sub*step]
/**
 * Deblocking filter, the alternated version from JVT-A003r1 H.26L draft.
 */
static inline void rv40_loop_filter(uint8_t *src, const int step, const int stride, const int dmode, const int lim0, const int lim1, const int mult, const int thr0, const int thr1, const int chroma, const int edge)
{
    int diffs[4][4];
    int s0 = 0, s1 = 0, s2 = 0, s3 = 0;
    uint8_t *ptr;
    int flag0 = 1, flag1 = 1;
    int llim0 = 3, llim1 = 3;
    int i, t, sflag;
    int p0, p1;
    int v88;

    for(i = 0, ptr = src; i < 4; i++, ptr += stride){
        diffs[i][0] = ptr[-2*step] - ptr[-1*step];
        diffs[i][1] = ptr[-2*step] - ptr[-3*step];
        diffs[i][2] = ptr[ 1*step] - ptr[ 0*step];
        diffs[i][3] = ptr[ 1*step] - ptr[ 2*step];
        s0 += diffs[i][0];
        s1 += diffs[i][1];
        s2 += diffs[i][2];
        s3 += diffs[i][3];
    }
    if(FFABS(s0) >= (thr0<<2)){
        llim0 = 1;
        flag0 = 0;
    }
    if(FFABS(s2) >= (thr0<<2)){
        llim1 = 1;
        flag1 = 0;
    }
    if(llim0 + llim1 == 2)
        return;

    if(!edge)
        flag0 = flag1 = 0;
    if(flag0 && FFABS(s1) >= thr1)
        flag0 = 0;
    if(flag1 && FFABS(s3) >= thr1)
        flag1 = 0;

    v88 = (lim0 + lim1 + llim0 + llim1) >> 1;
    if(flag0 + flag1 == 2){ /* strong filtering */
        for(i = 0; i < 4; i++, src += stride){
            t = src[0*step] - src[-1*step];
            if(!t) continue;
            sflag = (mult * FFABS(t)) >> 7;
            if(sflag > 1) continue;

            p0 = (RV40_STRONG_FILTER(src, step, -3, 1, -3) + rv40_dither_l[dmode + i]) >> 7;
            p1 = (RV40_STRONG_FILTER(src, step, -1, 3, -1) + rv40_dither_r[dmode + i]) >> 7;
            if(!sflag){
                src[-1*step] = p0;
                src[ 0*step] = p1;
            }else{
                if((src[-1*step] - p0) >= -v88 && (src[-1*step] - p0) <= v88)
                    src[-1*step] = p0;
                else
                    src[-1*step] = p1;
                if((src[ 0*step] - p1) >= -v88 && (src[ 0*step] - p1) <= v88)
                    src[ 0*step] = p1;
                else
                    src[ 0*step] = src[-1*step];
            }
            p0 = (RV40_STRONG_FILTER(src, step, -4, 0, -4) + rv40_dither_l[dmode + i]) >> 7;
            p1 = (RV40_STRONG_FILTER(src, step, -1, 3, -1) + rv40_dither_r[dmode + i]) >> 7;
            if(!sflag){
                src[-2*step] = p0;
                src[ 1*step] = p1;
            }else{
                if((src[-2*step] - p0) >= -v88 && (src[-2*step] - p0) <= v88)
                    src[-2*step] = p0;
                else
                    src[-2*step] += v88;
                if((src[ 1*step] - p1) >= -v88 && (src[ 1*step] - p1) <= v88)
                    src[ 1*step] = p1;
                else
                    src[ 1*step] += v88;
            }
            if(!chroma){
                src[-3*step] = (RV40_STRONG_FILTER(src, step, -4, -1, -3) + 64) >> 7;
                src[ 2*step] = (RV40_STRONG_FILTER(src, step,  0,  0,  2) + 64) >> 7;
            }
        }
    }else if(llim0 == 3 && llim1 == 3)
        for(i = 0; i < 4; i++, src += stride)
            rv40_weak_loop_filter(src, step, 1, 1, mult, lim0, lim1, v88, thr1,
                                  diffs[i][0], diffs[i][1], diffs[i][2], diffs[i][3]);
    else
        for(i = 0; i < 4; i++, src += stride)
            rv40_weak_loop_filter(src, step, llim0==3, llim1==3, mult, lim0>>1, lim1>>1, v88>>1, thr1,
                                  diffs[i][0], diffs[i][1], diffs[i][2], diffs[i][3]);
}

static void rv40_v_loop_filter(uint8_t *src, int stride, int dmode, int lim0, int lim1, int mult, int thr0, int thr1, int chroma, int edge){
    rv40_loop_filter(src, 1, stride, dmode, lim0, lim1, mult, thr0, thr1, chroma, edge);
}
static void rv40_h_loop_filter(uint8_t *src, int stride, int dmode, int lim0, int lim1, int mult, int thr0, int thr1, int chroma, int edge){
    rv40_loop_filter(src, stride, 1, dmode, lim0, lim1, mult, thr0, thr1, chroma, edge);
}

static void rv40_postprocess(RV40DecContext *r)
{
    MpegEncContext *s = &r->s;
    int mb_pos;
    int i, j;
    int no_up, no_left;
    uint8_t *Y, *U, *V;
    const int alpha = rv40_alpha_tab[s->qscale], beta = rv40_beta_tab[s->qscale];
    //XXX these are probably not correct
    const int thr = s->qscale, lim0 = rv40_filter_clip_tbl[1][s->qscale], lim1 = rv40_filter_clip_tbl[2][s->qscale];

    mb_pos = s->resync_mb_x + s->resync_mb_y * s->mb_stride;
    memset(r->intra_types_hist, -1, r->intra_types_stride * 4 * 2 * sizeof(int));
    s->first_slice_line = 1;
    s->mb_x= s->resync_mb_x;
    s->mb_y= s->resync_mb_y;
    ff_init_block_index(s);
    while(s->mb_num_left-- && s->mb_y < s->mb_height) {
        ff_update_block_index(s);
        if(IS_INTRA(s->current_picture_ptr->mb_type[mb_pos])){
            no_up = s->first_slice_line || !IS_INTRA(s->current_picture_ptr->mb_type[mb_pos - s->mb_stride]);
            no_left = !s->mb_x || (s->first_slice_line && s->mb_x == s->resync_mb_x) || !IS_INTRA(s->current_picture_ptr->mb_type[mb_pos - 1]);
            for(j = 0; j < 4; j++){
                for(i = 0; i < 4; i++){
                    Y = s->dest[0] + i*4 + j*4*s->linesize;
                    if(!j && !no_up)
                        rv40_h_loop_filter(Y, s->linesize, i*4+j, lim0, lim1, alpha, beta, thr, 0, 1);
                    if(j != 3)
                        rv40_h_loop_filter(Y + 4*s->linesize, s->linesize, i*4+j, lim0, lim1, alpha, beta, thr, 0, 0);
                    if(i || !no_left)
                        rv40_v_loop_filter(Y, s->linesize, i*4+j, lim0, lim1, alpha, beta, thr, 0, !i);
                }
            }
        }
        if (++s->mb_x == s->mb_width) {
            s->mb_x = 0;
            s->mb_y++;
            ff_init_block_index(s);
            mb_pos = s->mb_x + s->mb_y * s->mb_stride;
        }
        if(s->mb_x == s->resync_mb_x)
            s->first_slice_line=0;
    }

}
/** @} */ // recons group end

/**
 * Initialize decoder
 * @todo Maybe redone in some other way
 */
static int rv40_decode_init(AVCodecContext *avctx)
{
    RV40DecContext *r = avctx->priv_data;
    MpegEncContext *s = &r->s;

    static int tables_done = 0;

    r->rv30 = (avctx->codec_id == CODEC_ID_RV30);
    MPV_decode_defaults(s);
    s->avctx= avctx;
    s->out_format = FMT_H263;
    s->codec_id= avctx->codec_id;

    s->width = avctx->width;
    s->height = avctx->height;

    r->s.avctx = avctx;
    avctx->flags |= CODEC_FLAG_EMU_EDGE;
    r->s.flags |= CODEC_FLAG_EMU_EDGE;
    avctx->pix_fmt = PIX_FMT_YUV420P;
    avctx->max_b_frames = 8;
    avctx->has_b_frames = 1;
    s->low_delay = 0;

    if (MPV_common_init(s) < 0)
        return -1;

    ff_h264_pred_init(&r->h, CODEC_ID_RV40);

    r->intra_types_stride = (s->mb_width + 1) * 4;
    r->intra_types_hist = av_malloc(r->intra_types_stride * 4 * 2 * sizeof(int));
    r->intra_types = r->intra_types_hist + r->intra_types_stride * 4;

    r->mb_type = av_mallocz(r->s.mb_stride * r->s.mb_height * sizeof(int));

    if(!tables_done){
        rv40_init_tables();
        tables_done = 1;
    }
    r->prev_si.type = -1;
    if(r->rv30){
        if(avctx->extradata_size < 2){
            av_log(avctx, AV_LOG_ERROR, "Extradata is too small\n");
            return -1;
        }
        r->rpr = (avctx->extradata[1] & 7) >> 1;
        r->rpr = FFMIN(r->rpr + 1, 3);
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
    SliceInfo si;
    int i;
    int slice_count, *slice_offset;
    int zero_offset = {0};
    int last = 0;

    /* no supplementary picture */
    if (buf_size == 0) {
        /* special case for last picture */
        if (s->low_delay==0 && s->next_picture_ptr) {
            *pict= *(AVFrame*)s->next_picture_ptr;
            s->next_picture_ptr= NULL;

            *data_size = sizeof(AVFrame);
        }
    }

    if(avctx->slice_count){
        slice_count = avctx->slice_count;
        slice_offset = avctx->slice_offset;
    }else{
        slice_count = 1;
        slice_offset = &zero_offset;
    }

    for(i=0; i<slice_count; i++){
        int offset= slice_offset[i];
        int size;

        if(i+1 == slice_count)
            size= buf_size - offset;
        else
            size= slice_offset[i+1] - offset;

        r->si.size = size * 8;
        r->si.end = s->mb_width * s->mb_height;
        if(i+1 < slice_count){
            init_get_bits(&s->gb, buf+slice_offset[i+1], (buf_size-slice_offset[i+1])*8);
            if(!r->rv30 && rv40_parse_slice_header(r, &r->s.gb, &si) < 0){
                if(i+2 < slice_count)
                    size = slice_offset[i+2] - offset;
                else
                    size = buf_size - offset;
                r->si.size = size * 8;
            }else if(!r->rv30)
                r->si.end = si.start;
            if(r->rv30 && rv30_parse_slice_header(r, &r->s.gb, &si) >= 0)
                r->si.end = si.start;
        }
        r->slice_data = buf + offset;
        if(r->rv30)
            rv40_decode_slice(r, r->si.size, r->si.end, &last);
        else
            rv40_decode_slice(r, r->si.size, r->si.end, &last);
        if(last)
            break;
        s->mb_num_left = r->si.end - r->si.start;
        //rv40_postprocess(r);
    }

    if(last){
        r->prev_si.type = -1;
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
    av_freep(&r->mb_type);
    av_freep(&r->ssi.data);

    return 0;
}

AVCodec rv30_decoder = {
    "rv30",
    CODEC_TYPE_VIDEO,
    CODEC_ID_RV30,
    sizeof(RV40DecContext),
    rv40_decode_init,
    NULL,
    rv40_decode_end,
    rv40_decode_frame,
};

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
