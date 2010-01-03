/*
 * DSP functions for Indeo Video Interactive codecs (indeo4 and indeo5)
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
 * @file libavcodec/ivi_dsp.h
 * DSP functions (inverse transforms, motion compensations, wavelet recompostion)
 * for indeo video interactive codecs.
 */

#ifndef AVCODEC_IVI_DSP_H
#define AVCODEC_IVI_DSP_H

#include "avcodec.h"
#include "ivi_common.h"

/**
 *  5/3 - Wavelet recomposition filter for indeo5.
 *
 *  @param plane        [in]  pointer to the descriptor of the plane being processed
 *  @param dst          [out] pointer to the destination buffer
 *  @param dst_pitch    [in]  pitch of the destination buffer
 *  @param num_bands    [in]  number of wavelet bands to be processed
 */
void ff_ivi_recompose53(const IVIPlaneDesc *plane, uint8_t *dst, const int dst_pitch, const int num_bands);

/**
 *  Two-dimensional inverse slant 8x8 transform.
 *
 *  @param  in      [in]  pointer to the vector of transform coefficients
 *  @param  out     [out] pointer to the output buffer (frame)
 *  @param  pitch   [in]  pitch to move to the next y line
 *  @param  flags   [in]  pointer to the array of column flags:
 *                        != 0 - non_empty column, 0 - empty one
 *                        (this array must be filled by caller)
 */
void ff_ivi_inverse_slant_8x8(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags);

/**
 *  Two-dimensional inverse slant 4x4 transform.
 *
 *  @param  in      [in]  pointer to the vector of transform coefficients
 *  @param  out     [out] pointer to the output buffer (frame)
 *  @param  pitch   [in]  pitch to move to the next y line
 *  @param  flags   [in]  pointer to the array of column flags:
 *                        != 0 - non_empty column, 0 - empty one
 *                        (this array must be filled by caller)
 */
void ff_ivi_inverse_slant_4x4(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags);

/**
 *  This is a speed-up version of the inverse 2D slant transforms
 *  for the case if there is a non-zero DC coeff and all AC coeffs are zero.
 *  Performing the inverse slant transform in this case is equivalent to
 *  spreading (DC_coeff + 1)/2 over the whole block.
 *  It works much faster than performing the slant transform on a vector of zeroes.
 *
 *  @param  in          [in]  pointer to the dc coefficient
 *  @param  out         [out] pointer to the output buffer (frame)
 *  @param  pitch       [in]  pitch to move to the next y line
 *  @param  blk_size    [in]  transform block size
 */
void ff_ivi_dc_slant_2d(int32_t *in, int16_t *out, uint32_t pitch, int blk_size);

/**
 *  inverse 1D row slant transform
 *
 *  @param  in      [in]  pointer to the vector of transform coefficients
 *  @param  out     [out] pointer to the output buffer (frame)
 *  @param  pitch   [in]  pitch to move to the next y line
 *  @param  flags   [in]  pointer to the array of column flags (unused here)
 */
void ff_ivi_row_slant8(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags);

/**
 *  inverse 1D column slant transform
 *
 *  @param  in      [in]  pointer to the vector of transform coefficients
 *  @param  out     [out] pointer to the output buffer (frame)
 *  @param  pitch   [in]  pitch to move to the next y line
 *  @param  flags   [in]  pointer to the array of column flags:
 *                        != 0 - non_empty column, 0 - empty one
 *                        (this array must be filled by caller)
 */
void ff_ivi_col_slant8(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags);

void ff_ivi_dc_row_slant(int32_t *in, int16_t *out, uint32_t pitch, int blk_size);

void ff_ivi_dc_col_slant(int32_t *in, int16_t *out, uint32_t pitch, int blk_size);

/**
 *  8x8 block motion compensation with adding delta.
 *
 *  @param  buf     [in,out] pointer to the block in the current frame buffer containing delta
 *  @param  ref_buf [in]     pointer to the corresponding block in the reference frame
 *  @param  pitch   [in]     pitch for moving to the next y line
 *  @param  mc_type [in]     interpolation type
 */
void ff_ivi_mc_8x8_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type);

/**
 *  4x4 block motion compensation with adding delta.
 *
 *  @param  buf     [in,out] pointer to the block in the current frame buffer containing delta
 *  @param  ref_buf [in]     pointer to the corresponding block in the reference frame
 *  @param  pitch   [in]     pitch for moving to the next y line
 *  @param  mc_type [in]     interpolation type
 */
void ff_ivi_mc_4x4_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type);

/**
 *  Motion compensation without adding delta.
 *
 *  @param  buf     [in,out] pointer to the block in the current frame receiving the result
 *  @param  ref_buf [in]     pointer to the corresponding block in the reference frame
 *  @param  pitch   [in]     pitch for moving to the next y line
 *  @param  mc_type [in]     interpolation type
 */
void ff_ivi_mc_8x8_no_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type);

/**
 *  4x4 block motion compensation without adding delta.
 *
 *  @param  buf     [in,out] pointer to the block in the current frame receiving the result
 *  @param  ref_buf [in]     pointer to the corresponding block in the reference frame
 *  @param  pitch   [in]     pitch for moving to the next y line
 *  @param  mc_type [in]     interpolation type
 */
void ff_ivi_mc_4x4_no_delta(int16_t *buf, int16_t *ref_buf, uint32_t pitch, int mc_type);

#endif /* AVCODEC_IVI_DSP_H */
