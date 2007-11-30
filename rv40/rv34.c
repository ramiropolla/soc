/*
 * RV30/40 decoder common data
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
 * @file rv34.c
 * RV30/40 decoder common data.
 */

#include "avcodec.h"
#include "dsputil.h"
#include "mpegvideo.h"

#include "rv34vlc.h"
#include "rv34data.h"
#include "rv34.h"

//#define DEBUG

/** Translation of RV30/40 macroblock types to lavc ones */
static const int rv34_mb_type_to_lavc[12] = {
    MB_TYPE_INTRA, MB_TYPE_INTRA16x16, MB_TYPE_16x16,   MB_TYPE_8x8,
    MB_TYPE_16x16, MB_TYPE_16x16,      MB_TYPE_SKIP,    MB_TYPE_DIRECT2,
    MB_TYPE_16x8,  MB_TYPE_8x16,       MB_TYPE_DIRECT2, MB_TYPE_16x16
};


static RV34VLC intra_vlcs[NUM_INTRA_TABLES], inter_vlcs[NUM_INTER_TABLES];
static VLC omega_part_vlc;

/**
 * @defgroup vlc RV30/40 VLC generating functions
 * @{
 */

/**
 * Generate VLC from codeword lengths
 * @param bits   codeword lengths (zeroes are accepted)
 * @param size   length of input data
 * @param insyms symbols for input codes (NULL for default ones)
 */
static void rv34_gen_vlc(const uint8_t *bits, int size, VLC *vlc, const uint8_t *insyms)
{
    int i;
    int counts[17] = {0}, codes[17];
    uint16_t cw[size], syms[size];
    uint8_t bits2[size];
    int maxbits = 0, realsize = 0;

    for(i = 0; i < size; i++){
        if(bits[i]){
            bits2[realsize] = bits[i];
            syms[realsize] = insyms ? insyms[i] : i;
            realsize++;
            maxbits = FFMAX(maxbits, bits[i]);
            counts[bits[i]]++;
        }
    }

    codes[0] = 0;
    for(i = 0; i < 16; i++)
        codes[i+1] = (codes[i] + counts[i]) << 1;
    for(i = 0; i < realsize; i++)
        cw[i] = codes[bits2[i]]++;

    init_vlc_sparse(vlc, FFMIN(maxbits, 9), realsize,
                    bits2, 1, 1,
                    cw,    2, 2,
                    syms,  2, 2, INIT_VLC_USE_STATIC);
}

/**
 * Initialize all tables
 */
static void rv34_init_tables()
{
    int i, j, k;

    for(i = 0; i < NUM_INTRA_TABLES; i++){
        for(j = 0; j < 2; j++){
            rv34_gen_vlc(rv34_table_intra_cbppat      [i][j], CBPPAT_VLC_SIZE,   &intra_vlcs[i].cbppattern[j],     NULL);
            rv34_gen_vlc(rv34_table_intra_secondpat[i][j], OTHERBLK_VLC_SIZE, &intra_vlcs[i].second_pattern[j], NULL);
            rv34_gen_vlc(rv34_table_intra_thirdpat [i][j], OTHERBLK_VLC_SIZE, &intra_vlcs[i].third_pattern[j],  NULL);
            for(k = 0; k < 4; k++)
                rv34_gen_vlc(rv34_table_intra_cbp[i][j+k*2],  CBP_VLC_SIZE,      &intra_vlcs[i].cbp[j][k],         rv34_cbp_code);
        }
        for(j = 0; j < 4; j++)
            rv34_gen_vlc(rv34_table_intra_firstpat[i][j], FIRSTBLK_VLC_SIZE, &intra_vlcs[i].first_pattern[j], NULL);
        rv34_gen_vlc(rv34_intra_coeffvlc[i], COEFF_VLC_SIZE, &intra_vlcs[i].coefficient, NULL);
    }

    for(i = 0; i < NUM_INTER_TABLES; i++){
        rv34_gen_vlc(rv34_inter_cbppatvlc[i], CBPPAT_VLC_SIZE, &inter_vlcs[i].cbppattern[0], NULL);
        for(j = 0; j < 4; j++)
            rv34_gen_vlc(rv34_inter_cbpvlc[i][j], CBP_VLC_SIZE, &inter_vlcs[i].cbp[0][j], rv34_cbp_code);
        for(j = 0; j < 2; j++){
            rv34_gen_vlc(rv34_table_inter_firstpat [i][j], FIRSTBLK_VLC_SIZE, &inter_vlcs[i].first_pattern[j],  NULL);
            rv34_gen_vlc(rv34_table_inter_secondpat[i][j], OTHERBLK_VLC_SIZE, &inter_vlcs[i].second_pattern[j], NULL);
            rv34_gen_vlc(rv34_table_intra_thirdpat [i][j], OTHERBLK_VLC_SIZE, &inter_vlcs[i].third_pattern[j],  NULL);
        }
        rv34_gen_vlc(rv34_inter_coeffvlc[i], COEFF_VLC_SIZE, &inter_vlcs[i].coefficient, NULL);
    }

    init_vlc_sparse(&omega_part_vlc, OMEGA_BITS, NUM_OMEGA,
                    omega_part_vlc_bits,  1, 1,
                    omega_part_vlc_codes, 1, 1,
                    omega_part_vlc_syms,  1, 1, INIT_VLC_USE_STATIC);
}

/** @} */ // vlc group


/**
 * @defgroup transform RV30/40 inverse transform functions
 * @{
 */

/**
 * Real Video 4.0 inverse transform
 * Code is almost the same as in SVQ3, only scaling is different
 */
static void rv34_intra_inv_transform(DCTELEM *block, const int offset){
    int temp[16];
    int i;

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
 * RealVideo 3.0/4.0 inverse transform - special version
 *
 * Code is almost the same but final coefficients are multiplied by 1.5
 * and have no rounding
 */
static void rv34_intra_inv_transform_noround(DCTELEM *block, const int offset){
    int temp[16];
    int i;

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
 * @defgroup block RV30/40 4x4 block decoding functions
 * @{
 */

/**
 * Decode coded block pattern
 */
static int rv34_decode_cbp(GetBitContext *gb, RV34VLC *vlc, int table)
{
    int pattern, code, cbp=0;
    int ones;
    static const int cbp_masks[3] = {0x100000, 0x010000, 0x110000};
    static const int shifts[4] = { 0, 2, 8, 10 };
    int *curshift = shifts;
    int i, t, mask;

    code = get_vlc2(gb, vlc->cbppattern[table].table, 9, 2);
    pattern = code & 0xF;
    code >>= 4;

    ones = rv34_count_ones[pattern];

    for(mask = 8; mask; mask >>= 1, curshift++){
        if(!(pattern & mask)) continue;
        cbp |= get_vlc2(gb, vlc->cbp[table][ones].table, vlc->cbp[table][ones].bits, 1) << curshift[0];
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

static inline void rv34_decode_block(DCTELEM *dst, GetBitContext *gb, RV34VLC *rvlc, int fc, int sc)
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
static inline void rv34_dequant4x4(DCTELEM *block, int offset, int Qdc, int Q)
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
static inline void rv34_dequant4x4_16x16(DCTELEM *block, int offset, int Qdc, int Q)
{
    int i;

    block += offset;
    for(i = 0; i < 3; i++)
         block[rv34_dezigzag[i]] = (block[rv34_dezigzag[i]] * Qdc + 8) >> 4;
    for(; i < 16; i++)
         block[rv34_dezigzag[i]] = (block[rv34_dezigzag[i]] * Q + 8) >> 4;
}
/** @} */ //block functions


/**
 * @defgroup bitstream RV30/40 bitstream parsing
 * @{
 */

static inline int decode210(GetBitContext *gb){
    if (get_bits1(gb))
        return 0;
    else
        return 2 - get_bits1(gb);
}

/**
 * Decode staring slice position
 * @todo maybe replace with ff_h263_decode_mba() ?
 */
int ff_rv34_get_start_offset(GetBitContext *gb, int mb_size)
{
    int i;
    for(i = 0; i < 5; i++)
        if(rv34_mb_max_sizes[i] > mb_size)
            break;
    return rv34_mb_bits_sizes[i];
}

/**
 * Select VLC set for decoding from current quantizer, modifier and frame type
 */
static inline RV34VLC* choose_vlc_set(int quant, int mod, int type)
{
    if(mod == 2 && quant < 19) quant += 10;
    else if(mod && quant < 26) quant += 5;
    return type ? &inter_vlcs[rv34_quant_to_vlc_set[1][av_clip(quant, 0, 30)]]
                : &intra_vlcs[rv34_quant_to_vlc_set[0][av_clip(quant, 0, 30)]];
}

/**
 * Decode variable-length code constructed from variable-length codes
 * similar to Even-Rodeh and Elias Omega codes
 *
 * Code is constructed from bit chunks of even length (odd length means end of code)
 * and chunks are coded with variable-length codes too
 */
int ff_rv34_get_omega(GetBitContext *gb)
{
    int code = 1, t, tb;

    for(;;){
        t = get_vlc2(gb, omega_part_vlc.table, OMEGA_BITS, 1);
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
int ff_rv34_get_omega_signed(GetBitContext *gb)
{
    int code;

    code = ff_rv34_get_omega(gb);
    if(code & 1)
        return -(code >> 1);
    else
        return code >> 1;
}

/**
 * Decode quantizer difference and return modified quantizer
 */
static inline int rv34_decode_dquant(GetBitContext *gb, int quant)
{
    if(get_bits1(gb))
        return quant + rv34_dquant_tab[quant * 2 + get_bits1(gb)];
    else
        return get_bits(gb, 5);
}

/** @} */ //bitstream functions

/**
 * @defgroup mv motion vector related code (prediction, reconstruction, motion compensation)
 * @{
 */

/** Macroblock partition width in 8x8 blocks */
static const uint8_t part_sizes_w[RV34_MB_TYPES] = { 2, 2, 2, 1, 2, 2, 2, 2, 2, 1, 2, 2 };

/** Macroblock partition height in 8x8 blocks */
static const uint8_t part_sizes_h[RV34_MB_TYPES] = { 2, 2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2 };

/**
 * Motion vectors prediction
 *
 * Motion prediction performed for the block by using median prediction of
 * motion vector from the left, top and right top blocks but in corener cases
 * some other vectors may be used instead
 */
static void rv34_pred_mv(RV34DecContext *r, int block_type, int subblock_no)
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
    case RV34_MB_P_16x16:
    case RV34_MB_P_MIX16x16:
        if(!no_C){
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][1];
        }
        break;
    case RV34_MB_P_8x8:
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
    case RV34_MB_P_16x8:
        mv_pos += subblock_no*s->b8_stride;
        no_B &= ~subblock_no;
        no_C |= subblock_no;
        if(!no_C){
            C[0] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][0];
            C[1] = s->current_picture_ptr->motion_val[0][mv_pos-s->b8_stride+2][1];
        }
        break;
    case RV34_MB_P_8x16:
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
static inline void rv34_pred_b_vector(int A[2], int B[2], int C[2], int no_A, int no_B, int no_C, int *mx, int *my)
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
static void rv34_pred_mv_b(RV34DecContext *r, int block_type)
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
        if(r->mb_type[mb_pos - 1] != RV34_MB_B_FORWARD  && r->mb_type[mb_pos - 1] != RV34_MB_B_BIDIR)
            no_A[0] = 1;
        if(r->mb_type[mb_pos - 1] != RV34_MB_B_BACKWARD && r->mb_type[mb_pos - 1] != RV34_MB_B_BIDIR)
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
        if(r->mb_type[mb_pos - s->mb_stride] != RV34_MB_B_FORWARD  && r->mb_type[mb_pos - s->mb_stride] != RV34_MB_B_BIDIR)
            no_B[0] = 1;
        if(r->mb_type[mb_pos - s->mb_stride] != RV34_MB_B_BACKWARD && r->mb_type[mb_pos - s->mb_stride] != RV34_MB_B_BIDIR)
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
        if(r->mb_type[mb_pos - s->mb_stride + 1] != RV34_MB_B_FORWARD  && r->mb_type[mb_pos - s->mb_stride + 1] != RV34_MB_B_BIDIR)
            no_C[0] = 1;
        if(r->mb_type[mb_pos - s->mb_stride + 1] != RV34_MB_B_BACKWARD && r->mb_type[mb_pos - s->mb_stride + 1] != RV34_MB_B_BIDIR)
            no_C[1] = 1;
        c_mv_pos = mv_pos - s->b8_stride + 2;
    }else if(r->avail[3]){
        no_C[0] = no_C[1] = 0;
        if(r->mb_type[mb_pos - s->mb_stride - 1] != RV34_MB_B_FORWARD  && r->mb_type[mb_pos - s->mb_stride - 1] != RV34_MB_B_BIDIR)
            no_C[0] = 1;
        if(r->mb_type[mb_pos - s->mb_stride - 1] != RV34_MB_B_BACKWARD && r->mb_type[mb_pos - s->mb_stride - 1] != RV34_MB_B_BIDIR)
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
    case RV34_MB_B_FORWARD:
        rv34_pred_b_vector(A[0], B[0], C[0], no_A[0], no_B[0], no_C[0], &mx[0], &my[0]);
        r->dmv[1][0] = 0;
        r->dmv[1][1] = 0;
        break;
    case RV34_MB_B_BACKWARD:
        r->dmv[1][0] = r->dmv[0][0];
        r->dmv[1][1] = r->dmv[0][1];
        r->dmv[0][0] = 0;
        r->dmv[0][1] = 0;
        rv34_pred_b_vector(A[1], B[1], C[1], no_A[1], no_B[1], no_C[1], &mx[1], &my[1]);
        break;
    case RV34_MB_B_BIDIR:
        rv34_pred_b_vector(A[0], B[0], C[0], no_A[0], no_B[0], no_C[0], &mx[0], &my[0]);
        rv34_pred_b_vector(A[1], B[1], C[1], no_A[1], no_B[1], no_C[1], &mx[1], &my[1]);
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
static inline void rv34_mc(RV34DecContext *r, const int block_type,
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
    if(block_type == RV34_MB_P_16x8){
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        Y    += 8;
        srcY += 8;
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        s->dsp.put_h264_chroma_pixels_tab[0]   (U, srcU, s->uvlinesize, 4, uvmx, uvmy);
        s->dsp.put_h264_chroma_pixels_tab[0]   (V, srcV, s->uvlinesize, 4, uvmx, uvmy);
    }else if(block_type == RV34_MB_P_8x16){
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        Y    += 8 * s->linesize;
        srcY += 8 * s->linesize;
        s->dsp.put_h264_qpel_pixels_tab[1][dxy](Y, srcY, s->linesize);
        s->dsp.put_h264_chroma_pixels_tab[1]   (U, srcU, s->uvlinesize, 8, uvmx, uvmy);
        s->dsp.put_h264_chroma_pixels_tab[1]   (V, srcV, s->uvlinesize, 8, uvmx, uvmy);
    }else if(block_type == RV34_MB_P_8x8){
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
static inline void rv34_mc_b(RV34DecContext *r, const int block_type)
{
    MpegEncContext *s = &r->s;
    uint8_t *srcY, *srcU, *srcV;
    int dxy, mx, my, uvmx, uvmy, src_x, src_y, uvsrc_x, uvsrc_y;
    int mv_pos = s->mb_x * 2 + s->mb_y * 2 * s->b8_stride;

    if(block_type != RV34_MB_B_BACKWARD){
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
    if(block_type == RV34_MB_B_DIRECT){
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
static inline void rv34_mc_b_interp(RV34DecContext *r, const int block_type)
{
    MpegEncContext *s = &r->s;
    uint8_t *srcY, *srcU, *srcV;
    int dxy, mx, my, uvmx, uvmy, src_x, src_y, uvsrc_x, uvsrc_y;
    int mv_pos = s->mb_x * 2 + s->mb_y * 2 * s->b8_stride;

    mx = s->current_picture_ptr->motion_val[1][mv_pos][0];
    my = s->current_picture_ptr->motion_val[1][mv_pos][1];
    if(block_type == RV34_MB_B_DIRECT){
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

/** Number of motion vectors in each macroblock type */
static const int num_mvs[RV34_MB_TYPES] = { 0, 0, 1, 4, 1, 1, 0, 0, 2, 2, 2, 1 };

/**
 * Decode motion vector differences
 * and perform motion vector reconstruction and motion compensation.
 */
static int rv34_decode_mv(RV34DecContext *r, int block_type)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int i, j;

    memset(r->dmv, 0, sizeof(r->dmv));
    for(i = 0; i < num_mvs[block_type]; i++){
        r->dmv[i][0] = ff_rv34_get_omega_signed(gb);
        r->dmv[i][1] = ff_rv34_get_omega_signed(gb);
    }
    switch(block_type){
    case RV34_MB_TYPE_INTRA:
    case RV34_MB_TYPE_INTRA16x16:
        for(j = 0; j < 2; j++){
            for(i = 0; i < 2; i++){
                s->current_picture_ptr->motion_val[0][s->mb_x * 2 + s->mb_y * 2 * s->b8_stride + i + j*s->b8_stride][0] = 0;
                s->current_picture_ptr->motion_val[0][s->mb_x * 2 + s->mb_y * 2 * s->b8_stride + i + j*s->b8_stride][1] = 0;
            }
        }
        return 0;
    case RV34_MB_SKIP:
        if(s->pict_type == P_TYPE){
            rv34_pred_mv(r, block_type, 0);
            rv34_mc(r, block_type, 0, 0, 0, 2, 2);
            break;
        }
    case RV34_MB_B_DIRECT:
        rv34_pred_mv_b  (r, RV34_MB_B_DIRECT);
        rv34_mc_b       (r, RV34_MB_B_DIRECT);
        rv34_mc_b_interp(r, RV34_MB_B_DIRECT);
        break;
    case RV34_MB_P_16x16:
    case RV34_MB_P_MIX16x16:
        rv34_pred_mv(r, block_type, 0);
        rv34_mc(r, block_type, 0, 0, 0, 2, 2);
        break;
    case RV34_MB_B_FORWARD:
    case RV34_MB_B_BACKWARD:
        rv34_pred_mv_b  (r, block_type);
        rv34_mc_b       (r, block_type);
        break;
    case RV34_MB_P_16x8:
    case RV34_MB_P_8x16:
    case RV34_MB_B_BIDIR:
        rv34_pred_mv(r, block_type, 0);
        rv34_pred_mv(r, block_type, 1);
        if(block_type == RV34_MB_P_16x8){
            rv34_mc(r, block_type, 0, 0, 0,            2, 1);
            rv34_mc(r, block_type, 0, 8, s->b8_stride, 2, 1);
        }
        if(block_type == RV34_MB_P_8x16){
            rv34_mc(r, block_type, 0, 0, 0, 1, 2);
            rv34_mc(r, block_type, 8, 0, 1, 1, 2);
        }
        if(block_type == RV34_MB_B_BIDIR){
            rv34_pred_mv_b  (r, block_type);
            rv34_mc_b       (r, block_type);
            rv34_mc_b_interp(r, block_type);
        }
        break;
    case RV34_MB_P_8x8:
        for(i=0;i< 4;i++){
            rv34_pred_mv(r, block_type, i);
            rv34_mc(r, block_type, (i&1)<<3, (i&2)<<2, (i&1)+(i>>1)*s->b8_stride, 1, 1);
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
/** Mapping of RV30/40 intra prediction types to standard H.264 types */
static const int ittrans[9] = {
 DC_PRED, VERT_PRED, HOR_PRED, DIAG_DOWN_RIGHT_PRED, DIAG_DOWN_LEFT_PRED,
 VERT_RIGHT_PRED, VERT_LEFT_PRED, HOR_UP_PRED, HOR_DOWN_PRED,
};

/** Mapping of RV30/40 intra 16x16 prediction types to standard H.264 types */
static const int ittrans16[4] = {
 DC_PRED8x8, VERT_PRED8x8, HOR_PRED8x8, PLANE_PRED8x8,
};

/**
 * Perform 4x4 intra prediction
 */
static void rv34_pred_4x4_block(RV34DecContext *r, uint8_t *dst, int stride, int itype, int no_up, int no_left, int no_down, int no_right)
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
static void rv34_add_4x4_block(uint8_t *dst, int stride, DCTELEM block[64], int off)
{
    int x, y;
    for(y = 0; y < 4; y++)
        for(x = 0; x < 4; x++)
            dst[x + y*stride] = av_clip_uint8(dst[x + y*stride] + block[off + x+y*8]);
}

static inline int adjust_pred16(int itype, int no_up, int no_left)
{
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
    return itype;
}

static void rv34_output_macroblock(RV34DecContext *r, int *intra_types, int cbp, int is16)
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
                rv34_pred_4x4_block(r, YY, s->linesize, ittrans[intra_types[i]], no_up, no_left, i || (j==3), no_topright);
                no_left = 0;
                if(!(cbp & 1)) continue;
                rv34_add_4x4_block(YY, s->linesize, s->block[(i>>1)+(j&2)], (i&1)*4+(j&1)*32);
            }
            no_up = 0;
            Y += s->linesize * 4;
            intra_types += s->b4_stride;
        }
        intra_types -= s->b4_stride * 4;
        no_up = !r->avail[1];
        for(j = 0; j < 2; j++){
            no_left = !r->avail[0];
            for(i = 0; i < 2; i++, cbp >>= 1, no_left = 0){
                no_topright = no_up || (i && j) || (i && !j && (s->mb_x-1) == s->mb_width);
                rv34_pred_4x4_block(r, U + i*4 + j*4*s->uvlinesize, s->uvlinesize, ittrans[intra_types[i*2+j*2*s->b4_stride]], no_up, no_left, i || j, no_topright);
                rv34_pred_4x4_block(r, V + i*4 + j*4*s->uvlinesize, s->uvlinesize, ittrans[intra_types[i*2+j*2*s->b4_stride]], no_up, no_left, i || j, no_topright);
                if(cbp & 0x01)
                    rv34_add_4x4_block(U + i*4 + j*4*s->uvlinesize, s->uvlinesize, s->block[4], i*4+j*32);
                if(cbp & 0x10)
                    rv34_add_4x4_block(V + i*4 + j*4*s->uvlinesize, s->uvlinesize, s->block[5], i*4+j*32);
            }
            no_up = 0;
        }
    }else{
        no_left = !r->avail[0];
        itype = ittrans16[intra_types[0]];
        itype = adjust_pred16(itype, no_up, no_left);
        r->h.pred16x16[itype](Y, s->linesize);
        dsp->add_pixels_clamped(s->block[0], Y, s->current_picture.linesize[0]);
        dsp->add_pixels_clamped(s->block[1], Y + 8, s->current_picture.linesize[0]);
        Y += s->current_picture.linesize[0] * 8;
        dsp->add_pixels_clamped(s->block[2], Y, s->current_picture.linesize[0]);
        dsp->add_pixels_clamped(s->block[3], Y + 8, s->current_picture.linesize[0]);

        itype = ittrans16[intra_types[0]];
        if(itype == PLANE_PRED8x8) itype = DC_PRED8x8;
        itype = adjust_pred16(itype, no_up, no_left);
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
static int rv34_decode_mb_header(RV34DecContext *r, int *intra_types)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int mb_pos = s->mb_x + s->mb_y * s->mb_stride;
    int i, t;

    if(!r->si.type){
        r->is16 = get_bits1(gb);
        if(!r->is16 && !r->rv30){
            if(!get_bits1(gb))
                av_log(s->avctx, AV_LOG_ERROR, "Need DQUANT\n");
        }
        s->current_picture_ptr->mb_type[mb_pos] = r->is16 ? MB_TYPE_INTRA16x16 : MB_TYPE_INTRA;
        r->block_type = r->is16 ? RV34_MB_TYPE_INTRA16x16 : RV34_MB_TYPE_INTRA;
    }else{
        r->block_type = r->decode_mb_info(r);
        if(r->block_type == -1)
            return -1;
        s->current_picture_ptr->mb_type[mb_pos] = rv34_mb_type_to_lavc[r->block_type];
        r->mb_type[mb_pos] = r->block_type;
        if(r->block_type == RV34_MB_SKIP){
            if(s->pict_type == P_TYPE)
                r->mb_type[mb_pos] = RV34_MB_P_16x16;
            if(s->pict_type == B_TYPE)
                r->mb_type[mb_pos] = RV34_MB_B_DIRECT;
        }
        r->is16 = !!IS_INTRA16x16(s->current_picture_ptr->mb_type[mb_pos]);
        rv34_decode_mv(r, r->block_type);
        if(r->block_type == RV34_MB_SKIP){
            for(i = 0; i < 16; i++)
                intra_types[(i & 3) + (i>>2) * s->b4_stride] = 0;
            return 0;
        }
        r->chroma_vlc = 1;
        r->luma_vlc   = 0;
    }
    if(IS_INTRA(s->current_picture_ptr->mb_type[mb_pos])){
        if(!r->is16){
            if(r->decode_intra_types(r, gb, intra_types) < 0)
                return -1;
            r->luma_vlc   = 1;
        }else{
            t = get_bits(gb, 2);
            for(i = 0; i < 16; i++)
                intra_types[(i & 3) + (i>>2) * s->b4_stride] = t;
            r->luma_vlc   = 2;
        }
        r->chroma_vlc = 0;
        r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 0);
    }else{
        for(i = 0; i < 16; i++)
            intra_types[(i & 3) + (i>>2) * s->b4_stride] = 0;
        r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 1);
        if(r->mb_type[mb_pos] == RV34_MB_P_MIX16x16){
            r->is16 = 1;
            r->chroma_vlc = 1;
            r->luma_vlc   = 2;
            r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 0);
        }
    }
    return rv34_decode_cbp(gb, r->cur_vlcs, r->is16);
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


static void rv34_apply_differences(RV34DecContext *r, int cbp)
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

static int rv34_decode_macroblock(RV34DecContext *r, int *intra_types)
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
    cbp = cbp2 = rv34_decode_mb_header(r, intra_types);

    if(cbp == -1)
        return -1;

    luma_dc_quant = r->si.type ? r->luma_dc_quant_p[s->qscale] : r->luma_dc_quant_i[s->qscale];
    if(r->is16){
        memset(block16, 0, sizeof(block16));
        rv34_decode_block(block16, gb, r->cur_vlcs, 3, 0);
        rv34_dequant4x4_16x16(block16, 0, rv34_qscale_tab[luma_dc_quant],rv34_qscale_tab[s->qscale]);
        rv34_intra_inv_transform_noround(block16, 0);
    }

    for(i = 0; i < 16; i++, cbp >>= 1){
        if(!r->is16 && !(cbp & 1)) continue;
        blknum = ((i & 2) >> 1) + ((i & 8) >> 2);
        blkoff = ((i & 1) << 2) + ((i & 4) << 3);
        if(cbp & 1)
            rv34_decode_block(s->block[blknum] + blkoff, gb, r->cur_vlcs, r->luma_vlc, 0);
        if((cbp & 1) || r->is16){
            rv34_dequant4x4(s->block[blknum], blkoff, rv34_qscale_tab[luma_dc_quant],rv34_qscale_tab[s->qscale]);
            if(r->is16) //FIXME: optimize
                s->block[blknum][blkoff] = block16[(i & 3) | ((i & 0xC) << 1)];
            rv34_intra_inv_transform(s->block[blknum], blkoff);
        }
    }
    if(r->block_type == RV34_MB_P_MIX16x16)
        r->cur_vlcs = choose_vlc_set(r->si.quant, r->si.vlc_set, 1);
    for(; i < 24; i++, cbp >>= 1){
        if(!(cbp & 1)) continue;
        blknum = ((i & 4) >> 2) + 4;
        blkoff = ((i & 1) << 2) + ((i & 2) << 4);
        rv34_decode_block(s->block[blknum] + blkoff, gb, r->cur_vlcs, r->chroma_vlc, 1);
        rv34_dequant4x4(s->block[blknum], blkoff, rv34_qscale_tab[rv34_chroma_quant[1][s->qscale]],rv34_qscale_tab[rv34_chroma_quant[0][s->qscale]]);
        rv34_intra_inv_transform(s->block[blknum], blkoff);
    }
    if(IS_INTRA(s->current_picture_ptr->mb_type[s->mb_x + s->mb_y*s->mb_stride]))
        rv34_output_macroblock(r, intra_types, cbp2, r->is16);
    else
        rv34_apply_differences(r, cbp2);

    return 0;
}

static int check_slice_end(RV34DecContext *r, MpegEncContext *s)
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

static int rv34_decode_slice(RV34DecContext *r, int size, int end, int *last)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int mb_pos;
    int res;
    *last = 1;

    init_get_bits(&r->s.gb, r->slice_data, r->si.size);
    res = r->parse_slice_header(r, gb, &r->si);
    if(res < 0){
        av_log(s->avctx, AV_LOG_ERROR, "Incorrect or unknown slice header\n");
        *last = 0;
        return -1;
    }

    if ((s->mb_x == 0 && s->mb_y == 0) || s->current_picture_ptr==NULL) {
        if(s->width != r->si.width || s->height != r->si.height /*&& avcodec_check_dimensions(s->avctx, r->si.width, r->si.height) >= 0 */){
            av_log(s->avctx, AV_LOG_DEBUG, "Changing dimensions to %dx%d\n", r->si.width,r->si.height);
            MPV_common_end(s);
            s->width  = r->si.width;
            s->height = r->si.height;
            if(MPV_common_init(s) < 0)
                return -1;
            r->intra_types_hist = av_realloc(r->intra_types_hist, s->b4_stride * 4 * 2 * sizeof(int));
            r->intra_types = r->intra_types_hist + s->b4_stride * 4;
            r->mb_type = av_realloc(r->mb_type, r->s.mb_stride * r->s.mb_height * sizeof(*r->mb_type));
        }
        s->pict_type = r->si.type ? r->si.type : I_TYPE;
        if(MPV_frame_start(s, s->avctx) < 0)
            return -1;
        ff_er_frame_start(s);
        s->current_picture_ptr = &s->current_picture;
        s->mb_x = s->mb_y = 0;
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
    if(r->block_start != mb_pos){
        av_log(s->avctx, AV_LOG_ERROR, "Slice indicates MB offset %d, got %d\n", r->block_start, mb_pos);
        s->mb_x = r->block_start % s->mb_width;
        s->mb_y = r->block_start / s->mb_width;
    }
    memset(r->intra_types_hist, -1, s->b4_stride * 4 * 2 * sizeof(int));
    s->first_slice_line = 1;
    s->resync_mb_x= s->mb_x;
    s->resync_mb_y= s->mb_y;

    ff_init_block_index(s);
    while(!check_slice_end(r, s) && s->mb_num_left-- && s->mb_y < s->mb_height) {
        ff_update_block_index(s);
        s->dsp.clear_blocks(s->block[0]);

        if(rv34_decode_macroblock(r, r->intra_types + s->mb_x * 4 + 1) < 0)
            break;
        if (++s->mb_x == s->mb_width) {
            s->mb_x = 0;
            s->mb_y++;
            ff_init_block_index(s);

            memmove(r->intra_types_hist, r->intra_types, s->b4_stride * 4 * sizeof(int));
            memset(r->intra_types, -1, s->b4_stride * 4 * sizeof(int));
        }
        if(s->mb_x == s->resync_mb_x)
            s->first_slice_line=0;
    }
    ff_er_add_slice(s, s->resync_mb_x, s->resync_mb_y, s->mb_x-1, s->mb_y, AC_END|DC_END|MV_END);
    *last = 0;
    if(s->mb_y >= s->mb_height)
        *last = 1;
    if(r->bits > get_bits_count(gb) && show_bits(gb, r->bits-get_bits_count(gb)))
        *last = 1;

    return 0;
}

/** @} */ // recons group end

/**
 * Initialize decoder
 * @todo Maybe redone in some other way
 */
int ff_rv34_decode_init(AVCodecContext *avctx)
{
    RV34DecContext *r = avctx->priv_data;
    MpegEncContext *s = &r->s;

    static int tables_done = 0;

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
    avctx->has_b_frames = 1;
    s->low_delay = 0;

    if (MPV_common_init(s) < 0)
        return -1;

    ff_h264_pred_init(&r->h, CODEC_ID_RV40);

    r->intra_types_hist = av_malloc(s->b4_stride * 4 * 2 * sizeof(int));
    r->intra_types = r->intra_types_hist + s->b4_stride * 4;

    r->mb_type = av_mallocz(r->s.mb_stride * r->s.mb_height * sizeof(*r->mb_type));

    if(!tables_done){
        rv34_init_tables();
        tables_done = 1;
    }
    r->prev_si.type = -1;
    return 0;
}

static int get_slice_offset(AVCodecContext *avctx, uint8_t *buf, int n)
{
    if(avctx->slice_count) return avctx->slice_offset[n];
    else                   return AV_RL32(buf + n*8 - 4) == 1 ? AV_RL32(buf + n*8) :  AV_RB32(buf + n*8);
}

int ff_rv34_decode_frame(AVCodecContext *avctx,
                            void *data, int *data_size,
                            uint8_t *buf, int buf_size)
{
    RV34DecContext *r = avctx->priv_data;
    MpegEncContext *s = &r->s;
    AVFrame *pict = data;
    SliceInfo si;
    int i;
    int slice_count;
    uint8_t *slices_hdr = NULL;
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

    if(!avctx->slice_count){
        slice_count = (*buf++) + 1;
        slices_hdr = buf + 4;
        buf += 8 * slice_count;
    }else
        slice_count = avctx->slice_count;

    for(i=0; i<slice_count; i++){
        int offset= get_slice_offset(avctx, slices_hdr, i);
        int size;
        if(i+1 == slice_count)
            size= buf_size - offset;
        else
            size= get_slice_offset(avctx, slices_hdr, i+1) - offset;

        r->si.size = size * 8;
        r->si.end = s->mb_width * s->mb_height;
        if(i+1 < slice_count){
            init_get_bits(&s->gb, buf+get_slice_offset(avctx, slices_hdr, i+1), (buf_size-get_slice_offset(avctx, slices_hdr, i+1))*8);
            if(r->parse_slice_header(r, &r->s.gb, &si) < 0){
                if(i+2 < slice_count)
                    size = get_slice_offset(avctx, slices_hdr, i+2) - offset;
                else
                    size = buf_size - offset;
                r->si.size = size * 8;
            }else
                r->si.end = si.start;
        }
        r->slice_data = buf + offset;
        rv34_decode_slice(r, r->si.size, r->si.end, &last);
        s->mb_num_left = r->s.mb_x + r->s.mb_y*r->s.mb_width - r->si.start;
        if(last)
            break;
    }

    if(last){
        r->prev_si.type = -1;
        if(r->loop_filter)
            r->loop_filter(r);
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

int ff_rv34_decode_end(AVCodecContext *avctx)
{
    RV34DecContext *r = avctx->priv_data;

    MPV_common_end(&r->s);

    av_freep(&r->intra_types_hist);
    r->intra_types = NULL;
    av_freep(&r->mb_type);

    return 0;
}
