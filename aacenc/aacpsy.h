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
    AAC_PSY_NULL8,             // do nothing on frequencies but work with short windows
    AAC_PSY_3GPP,              // model following recommendations from 3GPP TS 26.403

    AAC_NB_PSY_MODELS
};

// data structures borrowed from aac.c with some minor modifications
/**
 * window sequences
 */
enum WindowSequence {
    ONLY_LONG_SEQUENCE,
    LONG_START_SEQUENCE,
    EIGHT_SHORT_SEQUENCE,
    LONG_STOP_SEQUENCE,
};

/**
 * special codebooks
 */
enum Codebook {
    ZERO_HCB       = 0,
    FIRST_PAIR_HCB = 5,
    ESC_HCB        = 11,
    NOISE_HCB      = 13,
    INTENSITY_HCB2 = 14,
    INTENSITY_HCB  = 15,
    ESC_FLAG       = 16,
};

/**
 * pulse tool
 */
typedef struct {
    int present;
    int num_pulse;
    int start;
    int offset[4];
    int amp[4];
} Pulse;

#define MAX_TAGID 16

/**
 * Program configuration - describes how channels are arranged. Either read from
 * stream (ID_PCE) or created based on a default fixed channel arrangement.
 */
typedef struct {
    int che_type[4][MAX_TAGID]; ///< channel element type with the first index as the first 4 raw_data_block IDs
    int mono_mixdown;           ///< The SCE tag to use if user requests mono   output, -1 if not available.
    int stereo_mixdown;         ///< The CPE tag to use if user requests stereo output, -1 if not available.
    int matrix_mixdown;         ///< The CPE tag to use if user requests matrixed stereo output, -1 if not available.
    int mixdown_coeff_index;    ///< 0-3
    int pseudo_surround;        ///< Mix surround channels out of phase.
} ProgramConfig;

/**
 * Individual Channel Stream
 */
typedef struct {
    int intensity_present;
    uint8_t max_sfb;            ///< number of scalefactor bands per group
    enum WindowSequence window_sequence;
    enum WindowSequence window_sequence_prev;
    uint8_t use_kb_window[2];   ///< If set, use Kaiser-Bessel window, otherwise use a sinus window.
    int num_window_groups;
    uint8_t grouping;
    uint8_t group_len[8];
    const uint8_t *swb_sizes;
    int num_swb;
    int num_windows;
    int tns_max_bands;
} IndividualChannelStream;

#define TNS_MAX_ORDER 20
/**
 * Temporal Noise Shaping
 */
typedef struct {
    int present;
    int n_filt[8];
    int length[8][4];
    int direction[8][4];
    int order[8][4];
    int coef_res[8];
    int coef_compress[8][4];
    int coef_len[8][4];
    const float *tmp2_map[8][4];
    int coef[8][4][TNS_MAX_ORDER];
} TemporalNoiseShaping;

/**
 * M/S joint channel coding
 */
typedef struct {
    int present;
    uint8_t mask[8][64];
} MidSideStereo;

/**
 * Single Channel Element
 * Used for both SCE and LFE elements
 */
typedef struct {
    int gain;                                 /**< Channel gain (not used by AAC bitstream).
                                               *   Note that this is applied before joint stereo decoding.
                                               *   Thus, when used inside CPE elements, both channels must have equal gain.
                                               */
    IndividualChannelStream ics;
    TemporalNoiseShaping tns;
    Pulse pulse;
    int zeroes[8][64];
    int sf_idx[8][64];
    enum Codebook cb[8][64];                  ///< codebooks
    int cb_run_end[8][64];                    ///< codebook run end points
    float sf[8][64];                          ///< scalefactors
    DECLARE_ALIGNED_16(float, coeffs[1024]);  ///< coefficients for IMDCT
    DECLARE_ALIGNED_16(float, saved[1024]);   ///< overlap
    DECLARE_ALIGNED_16(float, ret[1024]);     ///< PCM output
    DECLARE_ALIGNED_16(int,   icoefs[1024]);  ///< integer coefficients for coding
} SingleChannelElement;

/**
 * channel element - generic struct for SCE/CPE/CCE/LFE
 */
typedef struct {
    // CPE specific
    int common_window;     ///< Set if channels share a common 'IndividualChannelStream' in bitstream.
    MidSideStereo ms;
    // shared
    SingleChannelElement ch[2];
    // CCE specific
//    ChannelCoupling coup;
} ChannelElement;

// borrowing temporarily ends here

/**
 * context used by psychoacoustic model
 */
typedef struct AACPsyContext {
    AVCodecContext *avctx;
    DSPContext dsp;

    int window_type[2];
    int window_shape[2];
    const uint8_t *bands1024;
    int num_bands1024;
    const uint8_t *bands128;
    int num_bands128;

    const struct AACPsyModel *model;
    void* model_priv_data;
}AACPsyContext;

typedef struct AACPsyModel {
    const char *name;
    int   (*init)   (AACPsyContext *apc);
    void  (*window) (AACPsyContext *apc, int16_t *audio, int16_t *la, int channel, ChannelElement *cpe);
    void  (*process)(AACPsyContext *apc,int channel, ChannelElement *cpe);
    void  (*end)    (AACPsyContext *apc);
}AACPsyModel;

int ff_aac_psy_init(AACPsyContext *ctx, AVCodecContext *avctx, int model, int flags,
                    const uint8_t *bands1024, int num_bands1024,
                    const uint8_t *bands128,  int num_bands128);
void ff_aac_psy_suggest_window(AACPsyContext *ctx, int16_t *audio, int16_t *la, int channel, ChannelElement *cpe);
void ff_aac_psy_analyze(AACPsyContext *ctx, int channel, ChannelElement *cpe);
void ff_aac_psy_end(AACPsyContext *ctx);
#endif /* FFMPEG_AACPSY_H */

