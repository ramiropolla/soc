/*
 * RV40 decoder motion compensation functions
 * Copyright (c) 2008 Konstantin Shishkov
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
 * @file rv40dsp.c
 * RV40 decoder motion compensation functions
 */

#include "avcodec.h"
#include "dsputil.h"

#define RV40_LOWPASS(OPNAME, OP) \
static av_unused void OPNAME ## rv40_qpel8_h_lowpass(uint8_t *dst, uint8_t *src, int dstStride, int srcStride,\
                                                     const int h, const int C1, const int C2, const int C3){\
    uint8_t *cm = ff_cropTbl + MAX_NEG_CROP;\
    int i;\
    for(i=0; i<h; i++)\
    {\
        OP(dst[0], (src[-2] + src[ 3] - 5*(src[-1]+src[2]) + src[0]*C1 + src[1]*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[1], (src[-1] + src[ 4] - 5*(src[ 0]+src[3]) + src[1]*C1 + src[2]*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[2], (src[ 0] + src[ 5] - 5*(src[ 1]+src[4]) + src[2]*C1 + src[3]*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[3], (src[ 1] + src[ 6] - 5*(src[ 2]+src[5]) + src[3]*C1 + src[4]*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[4], (src[ 2] + src[ 7] - 5*(src[ 3]+src[6]) + src[4]*C1 + src[5]*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[5], (src[ 3] + src[ 8] - 5*(src[ 4]+src[7]) + src[5]*C1 + src[6]*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[6], (src[ 4] + src[ 9] - 5*(src[ 5]+src[8]) + src[6]*C1 + src[7]*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[7], (src[ 5] + src[10] - 5*(src[ 6]+src[9]) + src[7]*C1 + src[8]*C2 + (1<<(C3-1))) >> C3);\
        dst+=dstStride;\
        src+=srcStride;\
    }\
}\
\
static void OPNAME ## rv40_qpel8_v_lowpass(uint8_t *dst, uint8_t *src, int dstStride, int srcStride,\
                                           const int w, const int C1, const int C2, const int C3){\
    uint8_t *cm = ff_cropTbl + MAX_NEG_CROP;\
    int i;\
    for(i=0; i<w; i++)\
    {\
        const int srcB = src[-2*srcStride];\
        const int srcA = src[-1*srcStride];\
        const int src0 = src[0 *srcStride];\
        const int src1 = src[1 *srcStride];\
        const int src2 = src[2 *srcStride];\
        const int src3 = src[3 *srcStride];\
        const int src4 = src[4 *srcStride];\
        const int src5 = src[5 *srcStride];\
        const int src6 = src[6 *srcStride];\
        const int src7 = src[7 *srcStride];\
        const int src8 = src[8 *srcStride];\
        const int src9 = src[9 *srcStride];\
        const int src10= src[10*srcStride];\
        OP(dst[0*dstStride], (srcB + src3  - 5*(srcA+src2) + src0*C1 + src1*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[1*dstStride], (srcA + src4  - 5*(src0+src3) + src1*C1 + src2*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[2*dstStride], (src0 + src5  - 5*(src1+src4) + src2*C1 + src3*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[3*dstStride], (src1 + src6  - 5*(src2+src5) + src3*C1 + src4*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[4*dstStride], (src2 + src7  - 5*(src3+src6) + src4*C1 + src5*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[5*dstStride], (src3 + src8  - 5*(src4+src7) + src5*C1 + src6*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[6*dstStride], (src4 + src9  - 5*(src5+src8) + src6*C1 + src7*C2 + (1<<(C3-1))) >> C3);\
        OP(dst[7*dstStride], (src5 + src10 - 5*(src6+src9) + src7*C1 + src8*C2 + (1<<(C3-1))) >> C3);\
        dst++;\
        src++;\
    }\
}\
\
static void OPNAME ## rv40_qpel8_hv_lowpass(uint8_t *dst, uint8_t *src, int dstStride, int srcStride){\
    const int h=8;\
    uint8_t *cm = ff_cropTbl + MAX_NEG_CROP;\
    int i;\
    for(i=0; i<h; i++)\
    {\
        OP(dst[0], (src[0] + src[1] + src[0 + srcStride] + src[1 + srcStride] + 2) >> 2);\
        OP(dst[1], (src[1] + src[2] + src[1 + srcStride] + src[2 + srcStride] + 2) >> 2);\
        OP(dst[2], (src[2] + src[3] + src[2 + srcStride] + src[3 + srcStride] + 2) >> 2);\
        OP(dst[3], (src[3] + src[4] + src[3 + srcStride] + src[4 + srcStride] + 2) >> 2);\
        OP(dst[4], (src[4] + src[5] + src[4 + srcStride] + src[5 + srcStride] + 2) >> 2);\
        OP(dst[5], (src[5] + src[6] + src[5 + srcStride] + src[6 + srcStride] + 2) >> 2);\
        OP(dst[6], (src[6] + src[7] + src[6 + srcStride] + src[7 + srcStride] + 2) >> 2);\
        OP(dst[7], (src[7] + src[8] + src[7 + srcStride] + src[8 + srcStride] + 2) >> 2);\
        dst+=dstStride;\
        src+=srcStride;\
    }\
}\
\
static void OPNAME ## rv40_qpel16_v_lowpass(uint8_t *dst, uint8_t *src, int dstStride, int srcStride,\
                                            const int w, const int C1, const int C2, const int C3){\
    OPNAME ## rv40_qpel8_v_lowpass(dst  , src  , dstStride, srcStride, 8, C1, C2, C3);\
    OPNAME ## rv40_qpel8_v_lowpass(dst+8, src+8, dstStride, srcStride, 8, C1, C2, C3);\
    src += 8*srcStride;\
    dst += 8*dstStride;\
    OPNAME ## rv40_qpel8_v_lowpass(dst  , src  , dstStride, srcStride, w-8, C1, C2, C3);\
    OPNAME ## rv40_qpel8_v_lowpass(dst+8, src+8, dstStride, srcStride, w-8, C1, C2, C3);\
}\
\
static void OPNAME ## rv40_qpel16_h_lowpass(uint8_t *dst, uint8_t *src, int dstStride, int srcStride,\
                                            const int h, const int C1, const int C2, const int C3){\
    OPNAME ## rv40_qpel8_h_lowpass(dst  , src  , dstStride, srcStride, 8, C1, C2, C3);\
    OPNAME ## rv40_qpel8_h_lowpass(dst+8, src+8, dstStride, srcStride, 8, C1, C2, C3);\
    src += 8*srcStride;\
    dst += 8*dstStride;\
    OPNAME ## rv40_qpel8_h_lowpass(dst  , src  , dstStride, srcStride, h-8, C1, C2, C3);\
    OPNAME ## rv40_qpel8_h_lowpass(dst+8, src+8, dstStride, srcStride, h-8, C1, C2, C3);\
}\
\
static void OPNAME ## rv40_qpel16_hv_lowpass(uint8_t *dst, uint8_t *src, int dstStride, int srcStride){\
    OPNAME ## rv40_qpel8_hv_lowpass(dst  , src  , dstStride, srcStride);\
    OPNAME ## rv40_qpel8_hv_lowpass(dst+8, src+8, dstStride, srcStride);\
    src += 8*srcStride;\
    dst += 8*dstStride;\
    OPNAME ## rv40_qpel8_hv_lowpass(dst  , src  , dstStride, srcStride);\
    OPNAME ## rv40_qpel8_hv_lowpass(dst+8, src+8, dstStride, srcStride);\
}\
\

#define RV40_MC(OPNAME, SIZE) \
static void OPNAME ## rv40_qpel ## SIZE ## _mc10_c(uint8_t *dst, uint8_t *src, int stride){\
    OPNAME ## rv40_qpel ## SIZE ## _h_lowpass(dst, src, stride, stride, SIZE, 52, 20, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc20_c(uint8_t *dst, uint8_t *src, int stride){\
    OPNAME ## rv40_qpel ## SIZE ## _h_lowpass(dst, src, stride, stride, SIZE, 20, 20, 5);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc30_c(uint8_t *dst, uint8_t *src, int stride){\
    OPNAME ## rv40_qpel ## SIZE ## _h_lowpass(dst, src, stride, stride, SIZE, 20, 52, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc01_c(uint8_t *dst, uint8_t *src, int stride){\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, src, stride, stride, SIZE, 52, 20, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc11_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 52, 20, 6);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 52, 20, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc21_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 20, 20, 5);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 52, 20, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc31_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 20, 52, 6);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 52, 20, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc02_c(uint8_t *dst, uint8_t *src, int stride){\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, src, stride, stride, SIZE, 20, 20, 5);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc12_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 52, 20, 6);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 20, 20, 5);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc22_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 20, 20, 5);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 20, 20, 5);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc32_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 20, 52, 6);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 20, 20, 5);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc03_c(uint8_t *dst, uint8_t *src, int stride){\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, src, stride, stride, SIZE, 20, 52, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc13_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 52, 20, 6);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 20, 52, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc23_c(uint8_t *dst, uint8_t *src, int stride){\
    uint8_t full[SIZE*(SIZE+5)];\
    uint8_t * const full_mid= full + SIZE*2;\
    put_rv40_qpel ## SIZE ## _h_lowpass(full, src - 2*stride, SIZE, stride, SIZE+5, 20, 20, 5);\
    OPNAME ## rv40_qpel ## SIZE ## _v_lowpass(dst, full_mid, stride, SIZE, SIZE, 20, 52, 6);\
}\
\
static void OPNAME ## rv40_qpel ## SIZE ## _mc33_c(uint8_t *dst, uint8_t *src, int stride){\
    OPNAME ## rv40_qpel ## SIZE ## _hv_lowpass(dst, src, stride, stride);\
}\
\

#define op_avg(a, b)  a = (((a)+cm[b]+1)>>1)
#define op_put(a, b)  a = cm[b]

RV40_LOWPASS(put_       , op_put)
RV40_LOWPASS(avg_       , op_avg)
RV40_MC(put_, 8)
RV40_MC(put_, 16)
RV40_MC(avg_, 8)
RV40_MC(avg_, 16)

void ff_rv40dsp_init(DSPContext* c, AVCodecContext *avctx) {
    c->put_rv40_qpel_pixels_tab[0][ 0] = c->put_h264_qpel_pixels_tab[0][0];
    c->put_rv40_qpel_pixels_tab[0][ 1] = put_rv40_qpel16_mc10_c;
    c->put_rv40_qpel_pixels_tab[0][ 2] = put_rv40_qpel16_mc20_c;
    c->put_rv40_qpel_pixels_tab[0][ 3] = put_rv40_qpel16_mc30_c;
    c->put_rv40_qpel_pixels_tab[0][ 4] = put_rv40_qpel16_mc01_c;
    c->put_rv40_qpel_pixels_tab[0][ 5] = put_rv40_qpel16_mc11_c;
    c->put_rv40_qpel_pixels_tab[0][ 6] = put_rv40_qpel16_mc21_c;
    c->put_rv40_qpel_pixels_tab[0][ 7] = put_rv40_qpel16_mc31_c;
    c->put_rv40_qpel_pixels_tab[0][ 8] = put_rv40_qpel16_mc02_c;
    c->put_rv40_qpel_pixels_tab[0][ 9] = put_rv40_qpel16_mc12_c;
    c->put_rv40_qpel_pixels_tab[0][10] = put_rv40_qpel16_mc22_c;
    c->put_rv40_qpel_pixels_tab[0][11] = put_rv40_qpel16_mc32_c;
    c->put_rv40_qpel_pixels_tab[0][12] = put_rv40_qpel16_mc03_c;
    c->put_rv40_qpel_pixels_tab[0][13] = put_rv40_qpel16_mc13_c;
    c->put_rv40_qpel_pixels_tab[0][14] = put_rv40_qpel16_mc23_c;
    c->put_rv40_qpel_pixels_tab[0][15] = put_rv40_qpel16_mc33_c;
    c->avg_rv40_qpel_pixels_tab[0][ 0] = c->avg_h264_qpel_pixels_tab[0][0];
    c->avg_rv40_qpel_pixels_tab[0][ 1] = avg_rv40_qpel16_mc10_c;
    c->avg_rv40_qpel_pixels_tab[0][ 2] = avg_rv40_qpel16_mc20_c;
    c->avg_rv40_qpel_pixels_tab[0][ 3] = avg_rv40_qpel16_mc30_c;
    c->avg_rv40_qpel_pixels_tab[0][ 4] = avg_rv40_qpel16_mc01_c;
    c->avg_rv40_qpel_pixels_tab[0][ 5] = avg_rv40_qpel16_mc11_c;
    c->avg_rv40_qpel_pixels_tab[0][ 6] = avg_rv40_qpel16_mc21_c;
    c->avg_rv40_qpel_pixels_tab[0][ 7] = avg_rv40_qpel16_mc31_c;
    c->avg_rv40_qpel_pixels_tab[0][ 8] = avg_rv40_qpel16_mc02_c;
    c->avg_rv40_qpel_pixels_tab[0][ 9] = avg_rv40_qpel16_mc12_c;
    c->avg_rv40_qpel_pixels_tab[0][10] = avg_rv40_qpel16_mc22_c;
    c->avg_rv40_qpel_pixels_tab[0][11] = avg_rv40_qpel16_mc32_c;
    c->avg_rv40_qpel_pixels_tab[0][12] = avg_rv40_qpel16_mc03_c;
    c->avg_rv40_qpel_pixels_tab[0][13] = avg_rv40_qpel16_mc13_c;
    c->avg_rv40_qpel_pixels_tab[0][14] = avg_rv40_qpel16_mc23_c;
    c->avg_rv40_qpel_pixels_tab[0][15] = avg_rv40_qpel16_mc33_c;
    c->put_rv40_qpel_pixels_tab[1][ 0] = c->put_h264_qpel_pixels_tab[1][0];
    c->put_rv40_qpel_pixels_tab[1][ 1] = put_rv40_qpel8_mc10_c;
    c->put_rv40_qpel_pixels_tab[1][ 2] = put_rv40_qpel8_mc20_c;
    c->put_rv40_qpel_pixels_tab[1][ 3] = put_rv40_qpel8_mc30_c;
    c->put_rv40_qpel_pixels_tab[1][ 4] = put_rv40_qpel8_mc01_c;
    c->put_rv40_qpel_pixels_tab[1][ 5] = put_rv40_qpel8_mc11_c;
    c->put_rv40_qpel_pixels_tab[1][ 6] = put_rv40_qpel8_mc21_c;
    c->put_rv40_qpel_pixels_tab[1][ 7] = put_rv40_qpel8_mc31_c;
    c->put_rv40_qpel_pixels_tab[1][ 8] = put_rv40_qpel8_mc02_c;
    c->put_rv40_qpel_pixels_tab[1][ 9] = put_rv40_qpel8_mc12_c;
    c->put_rv40_qpel_pixels_tab[1][10] = put_rv40_qpel8_mc22_c;
    c->put_rv40_qpel_pixels_tab[1][11] = put_rv40_qpel8_mc32_c;
    c->put_rv40_qpel_pixels_tab[1][12] = put_rv40_qpel8_mc03_c;
    c->put_rv40_qpel_pixels_tab[1][13] = put_rv40_qpel8_mc13_c;
    c->put_rv40_qpel_pixels_tab[1][14] = put_rv40_qpel8_mc23_c;
    c->put_rv40_qpel_pixels_tab[1][15] = put_rv40_qpel8_mc33_c;
    c->avg_rv40_qpel_pixels_tab[1][ 0] = c->avg_h264_qpel_pixels_tab[1][0];
    c->avg_rv40_qpel_pixels_tab[1][ 1] = avg_rv40_qpel8_mc10_c;
    c->avg_rv40_qpel_pixels_tab[1][ 2] = avg_rv40_qpel8_mc20_c;
    c->avg_rv40_qpel_pixels_tab[1][ 3] = avg_rv40_qpel8_mc30_c;
    c->avg_rv40_qpel_pixels_tab[1][ 4] = avg_rv40_qpel8_mc01_c;
    c->avg_rv40_qpel_pixels_tab[1][ 5] = avg_rv40_qpel8_mc11_c;
    c->avg_rv40_qpel_pixels_tab[1][ 6] = avg_rv40_qpel8_mc21_c;
    c->avg_rv40_qpel_pixels_tab[1][ 7] = avg_rv40_qpel8_mc31_c;
    c->avg_rv40_qpel_pixels_tab[1][ 8] = avg_rv40_qpel8_mc02_c;
    c->avg_rv40_qpel_pixels_tab[1][ 9] = avg_rv40_qpel8_mc12_c;
    c->avg_rv40_qpel_pixels_tab[1][10] = avg_rv40_qpel8_mc22_c;
    c->avg_rv40_qpel_pixels_tab[1][11] = avg_rv40_qpel8_mc32_c;
    c->avg_rv40_qpel_pixels_tab[1][12] = avg_rv40_qpel8_mc03_c;
    c->avg_rv40_qpel_pixels_tab[1][13] = avg_rv40_qpel8_mc13_c;
    c->avg_rv40_qpel_pixels_tab[1][14] = avg_rv40_qpel8_mc23_c;
    c->avg_rv40_qpel_pixels_tab[1][15] = avg_rv40_qpel8_mc33_c;
}
