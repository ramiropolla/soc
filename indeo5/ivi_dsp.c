/*
 * DSP functions for Indeo Video Interactive codecs (Indeo4 and Indeo5)
 *
 * Copyright (c) 2009 Maxim Poliakovski
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
 * @file libavcodec/ivi_dsp.c
 * DSP functions (inverse transforms, motion compensation, wavelet recompostions)
 * for Indeo Video Interactive codecs.
 */

#include "avcodec.h"
#include "dsputil.h"
#include "ivi_common.h"
#include "ivi_dsp.h"

/**
 *  5/3 wavelet recomposition filter for Indeo5
 */
void ff_ivi_recompose53(const IVIPlaneDesc *plane, uint8_t *dst,
                        const int dst_pitch, const int num_bands)
{
    int             x, y, indx;
    int32_t         p0, p1, p2, p3, tmp0, tmp1, tmp2;
    int32_t         b0_1, b0_2, b1_1, b1_2, b1_3, b2_1, b2_2, b2_3, b2_4, b2_5, b2_6;
    int32_t         b3_1, b3_2, b3_3, b3_4, b3_5, b3_6, b3_7, b3_8, b3_9;
    uint32_t        pitch, back_pitch;
    const IDWTELEM *b0_ptr, *b1_ptr, *b2_ptr, *b3_ptr;

    /* all bands should have the same pitch */
    pitch = plane->bands[0].pitch;

    /* pixels at the position "y-1" will be set to pixels at the "y" for the 1st iteration */
    back_pitch = 0;

    /* get pointers to the wavelet bands */
    b0_ptr = plane->bands[0].buf;
    b1_ptr = plane->bands[1].buf;
    b2_ptr = plane->bands[2].buf;
    b3_ptr = plane->bands[3].buf;

    for (y = 0; y < plane->height; y += 2) {
        /* load storage variables with values */
        if (num_bands > 0) {
            b0_1 = b0_ptr[0];
            b0_2 = b0_ptr[pitch];
        }

        if (num_bands > 1) {
            b1_1 = b1_ptr[back_pitch];
            b1_2 = b1_ptr[0];
            b1_3 = b1_1 - b1_2*6 + b1_ptr[pitch];
        }

        if (num_bands > 2) {
            b2_2 = b2_ptr[0];     // b2[x,  y  ]
            b2_3 = b2_2;          // b2[x+1,y  ] = b2[x,y]
            b2_5 = b2_ptr[pitch]; // b2[x  ,y+1]
            b2_6 = b2_5;          // b2[x+1,y+1] = b2[x,y+1]
        }

        if (num_bands > 3) {
            b3_2 = b3_ptr[back_pitch]; // b3[x  ,y-1]
            b3_3 = b3_2;               // b3[x+1,y-1] = b3[x  ,y-1]
            b3_5 = b3_ptr[0];          // b3[x  ,y  ]
            b3_6 = b3_5;               // b3[x+1,y  ] = b3[x  ,y  ]
            b3_8 = b3_2 - b3_5*6 + b3_ptr[pitch];
            b3_9 = b3_8;
        }

        for (x = 0, indx = 0; x < plane->width; x+=2, indx++) {
            /* some values calculated in the previous iterations can */
            /* be reused in the next ones, so do appropriate copying */
            b2_1 = b2_2; // b2[x-1,y  ] = b2[x,  y  ]
            b2_2 = b2_3; // b2[x  ,y  ] = b2[x+1,y  ]
            b2_4 = b2_5; // b2[x-1,y+1] = b2[x  ,y+1]
            b2_5 = b2_6; // b2[x  ,y+1] = b2[x+1,y+1]
            b3_1 = b3_2; // b3[x-1,y-1] = b3[x  ,y-1]
            b3_2 = b3_3; // b3[x  ,y-1] = b3[x+1,y-1]
            b3_4 = b3_5; // b3[x-1,y  ] = b3[x  ,y  ]
            b3_5 = b3_6; // b3[x  ,y  ] = b3[x+1,y  ]
            b3_7 = b3_8; // vert_HPF(x-1)
            b3_8 = b3_9; // vert_HPF(x  )

            p0 = p1 = p2 = p3 = 0;

            /* process the LL-band by applying LPF both vertically and horizontally */
            if (num_bands > 0) {
                tmp0 = b0_1;
                tmp2 = b0_2;
                b0_1 = b0_ptr[indx+1];
                b0_2 = b0_ptr[pitch+indx+1];
                tmp1 = tmp0 + b0_1;

                p0 =  tmp0 << 4;
                p1 =  tmp1 << 3;
                p2 = (tmp0 + tmp2) << 3;
                p3 = (tmp1 + tmp2 + b0_2) << 2;
            }

            /* process the HL-band by applying HPF vertically and LPF horizontally */
            if (num_bands > 1) {
                tmp0 = b1_2;
                tmp1 = b1_1;
                b1_2 = b1_ptr[indx+1];
                b1_1 = b1_ptr[back_pitch+indx+1];

                tmp2 = tmp1 - tmp0*6 + b1_3;
                b1_3 = b1_1 - b1_2*6 + b1_ptr[pitch+indx+1];

                p0 += (tmp0 + tmp1) << 3;
                p1 += (tmp0 + tmp1 + b1_1 + b1_2) << 2;
                p2 +=  tmp2 << 2;
                p3 += (tmp2 + b1_3) << 1;
            }

            /* process the LH-band by applying LPF vertically and HPF horizontally */
            if (num_bands > 2) {
                b2_3 = b2_ptr[indx+1];
                b2_6 = b2_ptr[pitch+indx+1];

                tmp0 = b2_1 + b2_2;
                tmp1 = b2_1 - b2_2*6 + b2_3;

                p0 += tmp0 << 3;
                p1 += tmp1 << 2;
                p2 += (tmp0 + b2_4 + b2_5) << 2;
                p3 += (tmp1 + b2_4 - b2_5*6 + b2_6) << 1;
            }

            /* process the HH-band by applying HPF both vertically and horizontally */
            if (num_bands > 3) {
                b3_6 = b3_ptr[indx+1];            // b3[x+1,y  ]
                b3_3 = b3_ptr[back_pitch+indx+1]; // b3[x+1,y-1]

                tmp0 = b3_1 + b3_4;
                tmp1 = b3_2 + b3_5;
                tmp2 = b3_3 + b3_6;

                b3_9 = b3_3 - b3_6*6 + b3_ptr[pitch+indx+1];

                p0 += (tmp0 + tmp1) << 2;
                p1 += (tmp0 - tmp1*6 + tmp2) << 1;
                p2 += (b3_7 + b3_8) << 1;
                p3 +=  b3_7 - b3_8*6 + b3_9;
            }

            /* output four pixels */
            dst[x]             = av_clip_uint8((p0 >> 6) + 128);
            dst[x+1]           = av_clip_uint8((p1 >> 6) + 128);
            dst[dst_pitch+x]   = av_clip_uint8((p2 >> 6) + 128);
            dst[dst_pitch+x+1] = av_clip_uint8((p3 >> 6) + 128);
        }// for x

        dst += dst_pitch << 1;

        back_pitch = -pitch;

        b0_ptr += pitch;
        b1_ptr += pitch;
        b2_ptr += pitch;
        b3_ptr += pitch;
    }
}

/** butterfly operation for the inverse slant transform */
#define IVI_SLANT_BFLY(x, y) t1 = x-y; x += y; y = t1;

/** This is a reflection a,b = 1/2, 5/4 for the inverse slant transform */
#define IVI_IREFLECT(s1, s2) {\
    t1 = s1 + ((s1 + s2*2 + 2) >> 2);\
    s2 = ((s1*2 - s2 + 2) >> 2) - s2;\
    s1 = t1;}

/** This is a reflection a,b = 1/2, 7/8 for the inverse slant transform */
#define IVI_SLANT_PART4(s1, s2) {\
    t1 = s2 + ((s1*4 - s2 + 4)  >> 3);\
    s2 = ((-s1 - s2*4 + 4) >> 3) + s1;\
    s1 = t1;}

/** inverse slant8 transform */
#define IVI_INV_SLANT8(s1, s4, s8, s5, s2, s6, s3, s7, d1, d2, d3, d4, d5, d6, d7, d8) {\
    IVI_SLANT_PART4(s4, s5);\
    IVI_SLANT_BFLY(s1, s5); IVI_SLANT_BFLY(s2, s6); IVI_SLANT_BFLY(s7, s3); IVI_SLANT_BFLY(s4, s8);\
    IVI_SLANT_BFLY(s1, s2); IVI_IREFLECT  (s4, s3); IVI_SLANT_BFLY(s5, s6); IVI_IREFLECT  (s8, s7);\
    IVI_SLANT_BFLY(s1, s4); IVI_SLANT_BFLY(s2, s3); IVI_SLANT_BFLY(s5, s8); IVI_SLANT_BFLY(s6, s7);\
    d1 = COMPENSATE(s1);\
    d2 = COMPENSATE(s2);\
    d3 = COMPENSATE(s3);\
    d4 = COMPENSATE(s4);\
    d5 = COMPENSATE(s5);\
    d6 = COMPENSATE(s6);\
    d7 = COMPENSATE(s7);\
    d8 = COMPENSATE(s8);}

/** inverse slant4 transform */
#define IVI_INV_SLANT4(s1, s4, s2, s3, d1, d2, d3, d4) {\
    IVI_SLANT_BFLY(s1, s2); IVI_IREFLECT  (s4, s3);\
    IVI_SLANT_BFLY(s1, s4); IVI_SLANT_BFLY(s2, s3);\
    d1 = COMPENSATE(s1);\
    d2 = COMPENSATE(s2);\
    d3 = COMPENSATE(s3);\
    d4 = COMPENSATE(s4);}


/**
 *  two-dimensional inverse slant 8x8 transform
 */
void ff_ivi_inverse_slant_8x8(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags)
{
    int     i, t1;
    int32_t *src, *dst, tmp[64];

    /* apply the InvSlant8 to all columns */
#define COMPENSATE(x) (x)
    src = in;
    dst = tmp;
    for (i = 0; i < 8; i++) {
        if (flags[i]) {
            IVI_INV_SLANT8(src[0], src[8], src[16], src[24], src[32], src[40], src[48], src[56],
                           dst[0], dst[8], dst[16], dst[24], dst[32], dst[40], dst[48], dst[56]);
        } else
            dst[0] = dst[8] = dst[16] = dst[24] = dst[32] = dst[40] = dst[48] = dst[56] = 0;

            src++;
            dst++;
    }
#undef COMPENSATE

    /* apply the InvSlant8 to all rows and output the resulting coeffs in the band buffer */
    /* the (x + 1)/2 term is needed to compensate the normalization of the forward transform */
#define COMPENSATE(x) ((x + 1)>>1)
    src = tmp;
    for (i = 0; i < 8; i++) {
        if (!src[0] && !src[1] && !src[2] && !src[3] && !src[4] && !src[5] && !src[6] && !src[7]) {
            memset(out, 0, 8*sizeof(int16_t));
        } else {
            IVI_INV_SLANT8(src[0], src[1], src[2], src[3], src[4], src[5], src[6], src[7],
                           out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7]);
        }
        src += 8;
        out += pitch;
    }
#undef COMPENSATE
}


/**
 *  two-dimensional inverse slant 4x4 transform
 */
void ff_ivi_inverse_slant_4x4(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags)
{
    int     i, t1;
    int32_t *src, *dst, tmp[16];

    /* apply the InvSlant4 to all columns */
#define COMPENSATE(x) (x)
    src = in;
    dst = tmp;
    for (i = 0; i < 4; i++) {
        if (flags[i]) {
            IVI_INV_SLANT4(src[0], src[4], src[8], src[12], dst[0], dst[4], dst[8], dst[12]);
        } else
            dst[0] = dst[4] = dst[8] = dst[12] = 0;

            src++;
            dst++;
    }
#undef COMPENSATE

    /* apply the InvSlant4 to all rows */
#define COMPENSATE(x) ((x + 1)>>1)
    src = tmp;
    for (i = 0; i < 4; i++) {
        if (!src[0] && !src[1] && !src[2] && !src[3]) {
            out[0] = out[1] = out[2] = out[3] = 0;
        } else {
            IVI_INV_SLANT4(src[0], src[1], src[2], src[3], out[0], out[1], out[2], out[3]);
        }
        src += 4;
        out += pitch;
    }
#undef COMPENSATE
}


/**
 *  DC-only inverse 2D slant transforms
 */
void ff_ivi_dc_slant_2d(int32_t *in, int16_t *out, uint32_t pitch, int blk_size)
{
    int     x, y;
    int16_t dc_coeff;

    dc_coeff = (*in + 1) >> 1;

    for (y = 0; y < blk_size; out += pitch, y++) {
        for (x = 0; x < blk_size; x++)
            out[x] = dc_coeff;
    }
}


/**
 *  inverse 1D row slant transform
 */
void ff_ivi_row_slant8(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags)
{
    int     i, t1;

    /* apply the InvSlant8 to all rows and output the resulting coeffs in the band buffer */
    /* the (x + 1)/2 term is needed to compensate the normalization of the forward transform */
#define COMPENSATE(x) ((x + 1)>>1)
    for (i = 0; i < 8; i++) {
        if (!in[0] && !in[1] && !in[2] && !in[3] && !in[4] && !in[5] && !in[6] && !in[7]) {
            memset(out, 0, 8*sizeof(int16_t));
        } else {
            IVI_INV_SLANT8( in[0],  in[1],  in[2],  in[3],  in[4],  in[5],  in[6],  in[7],
                           out[0], out[1], out[2], out[3], out[4], out[5], out[6], out[7]);
        }
        in += 8;
        out += pitch;
    }
#undef COMPENSATE
}


/**
 *  DC-only inverse row slant transform
 */
void ff_ivi_dc_row_slant(int32_t *in, int16_t *out, uint32_t pitch, int blk_size)
{
    int     x, y;
    int16_t dc_coeff;

    dc_coeff = (*in + 1) >> 1;

    for (x = 0; x < blk_size; x++)
        out[x] = dc_coeff;

    out += pitch;

    for (y = 1; y < blk_size; out += pitch, y++) {
        for (x = 0; x < blk_size; x++)
            out[x] = 0;
    }
}


/**
 *  inverse 1D column slant transform
 */
void ff_ivi_col_slant8(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags)
{
    int     i, t1, row2, row4, row8;

    row2 = pitch << 1;
    row4 = pitch << 2;
    row8 = pitch << 3;

    /* apply the InvSlant8 to all columns and output the resulting coeffs in the band buffer */
    /* the (x + 1)/2 term is needed to compensate the normalization of the forward transform */
#define COMPENSATE(x) ((x + 1)>>1)
    for (i = 0; i < 8; i++) {
        if (flags[i]) {
            IVI_INV_SLANT8(in[0], in[8], in[16], in[24], in[32], in[40], in[48], in[56],
                           out[0], out[pitch], out[row2], out[row2 + pitch], out[row4],
                           out[row4 + pitch],  out[row4 + row2], out[row8 - pitch]);
        } else {
            out[0] = out[pitch] = out[row2] = out[row2 + pitch] = out[row4] =
                    out[row4 + pitch] =  out[row4 + row2] = out[row8 - pitch] = 0;
        }

        in++;
        out++;
    }
#undef COMPENSATE
}


/**
 *  DC-only inverse column slant transform
 */
void ff_ivi_dc_col_slant(int32_t *in, int16_t *out, uint32_t pitch, int blk_size)
{
    int     x, y;
    int16_t dc_coeff;

    dc_coeff = (*in + 1) >> 1;

    for (y = 0; y < blk_size; out += pitch, y++) {
        out[0] = dc_coeff;
        for (x = 1; x < blk_size; x++)
            out[x] = 0;
    }
}


/**
 *  8x8 block motion compensation with adding delta
 */
void ff_ivi_mc_8x8_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type)
{
    int     i, j;
    int16_t *wptr;

    switch (mc_type) {
    case 0: /* fullpel (no interpolation) */
        for (i = 0; i < 8; i++, buf += pitch, ref_buf += pitch) {
            buf[0] += ref_buf[0];
            buf[1] += ref_buf[1];
            buf[2] += ref_buf[2];
            buf[3] += ref_buf[3];
            buf[4] += ref_buf[4];
            buf[5] += ref_buf[5];
            buf[6] += ref_buf[6];
            buf[7] += ref_buf[7];
        }
        break;
    case 1: /* horizontal halfpel interpolation */
        for (i = 0; i < 8; i++, buf += pitch, ref_buf += pitch)
            for (j = 0; j < 8; j++)
                buf[j] += (ref_buf[j] + ref_buf[j+1]) >> 1;
        break;
    case 2: /* vertical halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 8; i++, buf += pitch, wptr += pitch, ref_buf += pitch)
            for (j = 0; j < 8; j++)
                buf[j] += (ref_buf[j] + wptr[j]) >> 1;
        break;
    case 3: /* vertical and horizontal halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 8; i++, buf += pitch, wptr += pitch, ref_buf += pitch)
            for (j = 0; j < 8; j++)
                buf[j] += (ref_buf[j] + ref_buf[j+1] + wptr[j] + wptr[j+1]) >> 2;
        break;
    }
}


/**
 *  4x4 block motion compensation with adding delta
 */
void ff_ivi_mc_4x4_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type)
{
    int     i;
    int16_t *wptr;

    switch (mc_type) {
    case 0: /* fullpel (no interpolation) */
        for (i = 0; i < 4; i++, buf += pitch, ref_buf += pitch) {
            buf[0] += ref_buf[0];
            buf[1] += ref_buf[1];
            buf[2] += ref_buf[2];
            buf[3] += ref_buf[3];
        }
        break;
    case 1: /* horizontal halfpel interpolation */
        for (i = 0; i < 4; i++, buf += pitch, ref_buf += pitch) {
            buf[0] += (ref_buf[0] + ref_buf[1]) >> 1;
            buf[1] += (ref_buf[1] + ref_buf[2]) >> 1;
            buf[2] += (ref_buf[2] + ref_buf[3]) >> 1;
            buf[3] += (ref_buf[3] + ref_buf[4]) >> 1;
        }
        break;
    case 2: /* vertical halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 4; i++, buf += pitch, wptr += pitch, ref_buf += pitch) {
            buf[0] += (ref_buf[0] + wptr[0]) >> 1;
            buf[1] += (ref_buf[1] + wptr[1]) >> 1;
            buf[2] += (ref_buf[2] + wptr[2]) >> 1;
            buf[3] += (ref_buf[3] + wptr[3]) >> 1;
        }
        break;
    case 3: /* vertical and horizontal halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 4; i++, buf += pitch, wptr += pitch, ref_buf += pitch) {
            buf[0] += (ref_buf[0] + ref_buf[1] + wptr[0] + wptr[1]) >> 2;
            buf[1] += (ref_buf[1] + ref_buf[2] + wptr[1] + wptr[2]) >> 2;
            buf[2] += (ref_buf[2] + ref_buf[3] + wptr[2] + wptr[3]) >> 2;
            buf[3] += (ref_buf[3] + ref_buf[4] + wptr[3] + wptr[4]) >> 2;
        }
        break;
    }
}


/**
 *  motion compensation without adding delta
 */
void ff_ivi_mc_8x8_no_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type)
{
    int     i, j;
    int16_t *wptr;

    switch (mc_type) {
    case 0: /* fullpel (no interpolation, just copy) */
        for (i = 0; i < 8; i++, buf += pitch, ref_buf += pitch)
            memcpy(buf, ref_buf, 8*sizeof(int16_t)); /* FIXME: speed critical? */
        break;
    case 1: /* horizontal halfpel interpolation */
        for (i = 0; i < 8; i++, buf += pitch, ref_buf += pitch)
            for (j = 0; j < 8; j++)
                buf[j] = (ref_buf[j] + ref_buf[j+1]) >> 1;
        break;
    case 2: /* vertical halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 8; i++, buf += pitch, wptr += pitch, ref_buf += pitch)
            for (j = 0; j < 8; j++)
                buf[j] = (ref_buf[j] + wptr[j]) >> 1;
        break;
    case 3: /* vertical and horizontal halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 8; i++, buf += pitch, wptr += pitch, ref_buf += pitch)
            for (j = 0; j < 8; j++)
                buf[j] = (ref_buf[j] + ref_buf[j+1] + wptr[j] + wptr[j+1]) >> 2;
        break;
    }
}


/**
 *  4x4 block motion compensation without adding delta
 */
void ff_ivi_mc_4x4_no_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type)
{
    int     i;
    int16_t *wptr;

    switch (mc_type) {
    case 0: /* fullpel (no interpolation, just copy) */
        for (i = 0; i < 4; i++, buf += pitch, ref_buf += pitch) {
            buf[0] = ref_buf[0];
            buf[1] = ref_buf[1];
            buf[2] = ref_buf[2];
            buf[3] = ref_buf[3];
        }
        break;
    case 1: /* horizontal halfpel interpolation */
        for (i = 0; i < 4; i++, buf += pitch, ref_buf += pitch) {
            buf[0] = (ref_buf[0] + ref_buf[1]) >> 1;
            buf[1] = (ref_buf[1] + ref_buf[2]) >> 1;
            buf[2] = (ref_buf[2] + ref_buf[3]) >> 1;
            buf[3] = (ref_buf[3] + ref_buf[4]) >> 1;
        }
        break;
    case 2: /* vertical halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 4; i++, buf += pitch, wptr += pitch, ref_buf += pitch) {
            buf[0] = (ref_buf[0] + wptr[0]) >> 1;
            buf[1] = (ref_buf[1] + wptr[1]) >> 1;
            buf[2] = (ref_buf[2] + wptr[2]) >> 1;
            buf[3] = (ref_buf[3] + wptr[3]) >> 1;
        }
        break;
    case 3: /* vertical and horizontal halfpel interpolation */
        wptr = ref_buf + pitch;
        for (i = 0; i < 4; i++, buf += pitch, wptr += pitch, ref_buf += pitch) {
            buf[0] = (ref_buf[0] + ref_buf[1] + wptr[0] + wptr[1]) >> 2;
            buf[1] = (ref_buf[1] + ref_buf[2] + wptr[1] + wptr[2]) >> 2;
            buf[2] = (ref_buf[2] + ref_buf[3] + wptr[2] + wptr[3]) >> 2;
            buf[3] = (ref_buf[3] + ref_buf[4] + wptr[3] + wptr[4]) >> 2;
        }
        break;
    }
}
