/*
 * AAC encoder psychoacoustic model
 * Copyright (C) 2008 Konstantin Shishkov
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

#ifndef FFMPEG_AACPSY_H
#define FFMPEG_AACPSY_H

#include "avcodec.h"
#include "dsputil.h"

enum AACPsyModelType{
    AAC_PSY_NULL,              // do nothing on frequencies

    AAC_NB_PSY_MODELS
};

// data structures borrowed from aac.c with some minor modifications

/**
 * Individual Channel Stream
 */
typedef struct {
    int intensity_present;
    int max_sfb;
    int window_sequence;
    int window_shape;             ///< If set, use Kaiser-Bessel window, otherwise use a sinus window
    int window_shape_prev;
    int num_window_groups;
    uint8_t grouping;
    uint8_t group_len[8];
    const uint8_t *swb_sizes;
    int num_swb;
    int num_windows;
    int tns_max_bands;
} ics_struct;

/**
 * M/S joint channel coding
 */
typedef struct {
    int present;
    uint8_t mask[8][64];
} ms_struct;

/**
 * Single Channel Element
 * Used for both SCE and LFE elements
 */
typedef struct {
    int gain;                                 /**< Channel gain (not used by AAC bitstream).
                                               *   Note that this is applied before joint stereo decoding.
                                               *   Thus, when used inside CPE elements, both channels must have equal gain.
                                               */
    ics_struct ics;
    int zeroes[8][64];
    int sf_idx[8][64];
    int cb[8][64];                            ///< Codebooks
    float sf[8][64];                          ///< Scalefactors
    DECLARE_ALIGNED_16(float, coeffs[1024]);  ///< Coefficients for IMDCT
    DECLARE_ALIGNED_16(float, saved[1024]);   ///< Overlap
    DECLARE_ALIGNED_16(float, ret[1024]);     ///< PCM output
    DECLARE_ALIGNED_16(int,   icoefs[1024]);  ///< integer coefficients for coding
} sce_struct;

/**
 * Channel Pair Element
 */
typedef struct {
    int common_window;     ///< Set if channels share a common 'ics_struct' in bitstream
    ms_struct ms;
    sce_struct ch[2];
} cpe_struct;

// borrowing temporarily ends here

/**
 * context used by psychoacoustic model
 */
typedef struct AACPsyContext {
    AVCodecContext *avctx;
    DSPContext dsp;

    int window_type[2];
    int window_shape[2];
    const uint8_t *bands;
    int num_bands;

    const struct AACPsyModel *model;
    void* model_priv_data;
}AACPsyContext;

typedef struct AACPsyModel {
    const char *name;
    int   (*init)   (AACPsyContext *apc);
    void  (*window) (AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe);
    void  (*process)(AACPsyContext *apc, int16_t *audio, int channel, cpe_struct *cpe);
    void  (*end)    (AACPsyContext *apc);
}AACPsyModel;

int ff_aac_psy_init(AACPsyContext *ctx, AVCodecContext *avctx, int model, int flags,
                    const uint8_t *bands, int num_bands);
void ff_aac_psy_suggest_window(AACPsyContext *ctx, int16_t *audio, int channel, cpe_struct *cpe);
void ff_aac_psy_analyze(AACPsyContext *ctx, int16_t *audio, int channel, cpe_struct *cpe);
void ff_aac_psy_end(AACPsyContext *ctx);
#endif /* FFMPEG_AACPSY_H */

