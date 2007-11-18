/*
 * RV40 decoder
 * Copyright (c) 2007 Konstantin Shishkov
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

#include "rv34.h"
#include "rv40vlc2.h"
#include "rv40data.h"

static VLC aic_top_vlc;
static VLC aic_mode1_vlc[AIC_MODE1_NUM], aic_mode2_vlc[AIC_MODE2_NUM];
static VLC ptype_vlc[NUM_PTYPE_VLCS], btype_vlc[NUM_BTYPE_VLCS];

/**
 * Initialize all tables
 */
static void rv40_init_tables()
{
    int i;

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

static int rv40_parse_slice_header(RV34DecContext *r, GetBitContext *gb, SliceInfo *si)
{
    int t, mb_bits;
    int w = r->s.width, h = r->s.height;
    int mb_size;

    memset(si, 0, sizeof(SliceInfo));
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
    mb_bits = ff_rv34_get_start_offset(gb, mb_size);
    si->start = get_bits(gb, mb_bits);
    si->header_size = get_bits_count(gb);

    return 0;
}

/**
 * Decode 4x4 intra types array
 */
static int rv40_decode_intra_types(RV34DecContext *r, GetBitContext *gb, int *dst)
{
    MpegEncContext *s = &r->s;
    int i, j, k, v;
    int A, B, C;
    int pattern;
    int *ptr;

    for(i = 0; i < 4; i++, dst += s->b4_stride){
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
            A = ptr[-s->b4_stride + 1]; // it won't be used for the last coefficient in a row
            B = ptr[-s->b4_stride];
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
 * Decode macroblock information
 */
static int rv40_decode_mb_info(RV34DecContext *r)
{
    MpegEncContext *s = &r->s;
    GetBitContext *gb = &s->gb;
    int q, i;
    int prev_type = 0;
    int mb_pos = s->mb_x + s->mb_y * s->mb_stride;
    int blocks[RV34_MB_TYPES];
    int count = 0;

    if(!r->s.mb_skip_run)
        r->s.mb_skip_run = ff_rv34_get_omega(gb);

    if(--r->s.mb_skip_run)
         return RV34_MB_SKIP;

    memset(blocks, 0, sizeof(blocks));
    if(r->avail[0])
        blocks[r->mb_type[mb_pos - 1]]++;
    if(r->avail[1])
        blocks[r->mb_type[mb_pos - s->mb_stride]]++;
    if(r->avail[1] && r->avail[2])
        blocks[r->mb_type[mb_pos - s->mb_stride + 1]]++;
    if(r->avail[1] && r->avail[3])
        blocks[r->mb_type[mb_pos - s->mb_stride - 1]]++;

    for(i = 0; i < RV34_MB_TYPES; i++){
        if(blocks[i] > count){
            count = blocks[i];
            prev_type = i;
        }
    }
    if(s->pict_type == P_TYPE){
        if(prev_type == RV34_MB_SKIP) prev_type = RV34_MB_P_16x16;
        prev_type = block_num_to_ptype_vlc_num[prev_type];
        q = get_vlc2(gb, ptype_vlc[prev_type].table, PTYPE_VLC_BITS, 1);
        if(q < PBTYPE_ESCAPE)
            return q;
        q = get_vlc2(gb, ptype_vlc[prev_type].table, PTYPE_VLC_BITS, 1);
        av_log(s->avctx, AV_LOG_ERROR, "Dquant for P-frame\n");
    }else{
        prev_type = block_num_to_btype_vlc_num[prev_type];
        q = get_vlc2(gb, btype_vlc[prev_type].table, BTYPE_VLC_BITS, 1);
        if(q < PBTYPE_ESCAPE)
            return q;
        q = get_vlc2(gb, btype_vlc[prev_type].table, BTYPE_VLC_BITS, 1);
        av_log(s->avctx, AV_LOG_ERROR, "Dquant for B-frame\n");
    }
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
#define RV34_STRONG_FILTER(src, step, start, last, sub) \
     26*(src[start*step] + src[(start+1)*step] + src[(start+2)*step] + src[(start+3)*step] + src[last*step]) - src[last*step] - src[sub*step]
/**
 * Deblocking filter, the alternated version from JVT-A003r1 H.26L draft.
 */
static inline void rv40_adaptive_loop_filter(uint8_t *src, const int step, const int stride, const int dmode, const int lim0, const int lim1, const int mult, const int thr0, const int thr1, const int chroma, const int edge)
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

            p0 = (RV34_STRONG_FILTER(src, step, -3, 1, -3) + rv40_dither_l[dmode + i]) >> 7;
            p1 = (RV34_STRONG_FILTER(src, step, -1, 3, -1) + rv40_dither_r[dmode + i]) >> 7;
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
            p0 = (RV34_STRONG_FILTER(src, step, -4, 0, -4) + rv40_dither_l[dmode + i]) >> 7;
            p1 = (RV34_STRONG_FILTER(src, step, -1, 3, -1) + rv40_dither_r[dmode + i]) >> 7;
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
                src[-3*step] = (RV34_STRONG_FILTER(src, step, -4, -1, -3) + 64) >> 7;
                src[ 2*step] = (RV34_STRONG_FILTER(src, step,  0,  0,  2) + 64) >> 7;
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
    rv40_adaptive_loop_filter(src, 1, stride, dmode, lim0, lim1, mult, thr0, thr1, chroma, edge);
}
static void rv40_h_loop_filter(uint8_t *src, int stride, int dmode, int lim0, int lim1, int mult, int thr0, int thr1, int chroma, int edge){
    rv40_adaptive_loop_filter(src, stride, 1, dmode, lim0, lim1, mult, thr0, thr1, chroma, edge);
}

static void rv40_loop_filter(RV34DecContext *r)
{
    MpegEncContext *s = &r->s;
    int mb_pos;
    int i, j;
    int no_up, no_left;
    uint8_t *Y, *U, *V;
    const int alpha = rv40_alpha_tab[s->qscale], beta = rv40_beta_tab[s->qscale];
    //XXX these are probably not correct
    const int thr = s->qscale, lim0 = rv40_filter_clip_tbl[1][s->qscale], lim1 = rv40_filter_clip_tbl[2][s->qscale];

    s->first_slice_line = 1;
    s->mb_x= 0;
    s->mb_y= 0;
    mb_pos = 0;
    ff_init_block_index(s);
    s->mb_num_left = s->mb_width * s->mb_height;
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

/**
 * Initialize decoder
 */
static int rv40_decode_init(AVCodecContext *avctx)
{
    RV34DecContext *r = avctx->priv_data;
    static int tables_done = 0;

    r->rv30 = 0;
    ff_rv34_decode_init(avctx);
    if(!tables_done){
        rv40_init_tables();
        tables_done = 1;
    }
    r->parse_slice_header = rv40_parse_slice_header;
    r->decode_intra_types = rv40_decode_intra_types;
    r->decode_mb_info     = rv40_decode_mb_info;
    r->loop_filter        = rv40_loop_filter;
    r->luma_dc_quant_i = rv40_luma_quant[0];
    r->luma_dc_quant_p = rv40_luma_quant[1];
    return 0;
}

AVCodec rv40_decoder = {
    "rv40",
    CODEC_TYPE_VIDEO,
    CODEC_ID_RV40,
    sizeof(RV34DecContext),
    rv40_decode_init,
    NULL,
    ff_rv34_decode_end,
    ff_rv34_decode_frame,
};
