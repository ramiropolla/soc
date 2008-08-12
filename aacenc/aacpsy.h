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
#include "aac.h"

enum AACPsyModelType{
    AAC_PSY_NULL,              ///< do nothing on frequencies
    AAC_PSY_NULL8,             ///< do nothing on frequencies but work with short windows
    AAC_PSY_3GPP,              ///< model following recommendations from 3GPP TS 26.403

    AAC_NB_PSY_MODELS          ///< total number of psychoacoustic models
};

enum AACPsyModelMode{
    PSY_MODE_CBR,              ///< follow bitrate as closely as possible
    PSY_MODE_ABR,              ///< try to achieve bitrate but actual bitrate may differ significantly
    PSY_MODE_QUALITY,          ///< try to achieve set quality instead of bitrate
};

#define PSY_MODEL_MODE_MASK  0x0000000F ///< bit fields for storing mode (CBR, ABR, VBR)
#define PSY_MODEL_NO_PULSE   0x00000010 ///< disable pulse searching
#define PSY_MODEL_NO_SWITCH  0x00000020 ///< disable window switching
#define PSY_MODEL_NO_ST_ATT  0x00000040 ///< disable stereo attenuation
#define PSY_MODEL_NO_LOWPASS 0x00000080 ///< disable low-pass filtering

#define PSY_MODEL_NO_PREPROC (PSY_MODEL_NO_ST_ATT | PSY_MODEL_NO_LOWPASS)

#define PSY_MODEL_MODE(a)  ((a) & PSY_MODEL_MODE_MASK)

/**
 * context used by psychoacoustic model
 */
typedef struct AACPsyContext {
    AVCodecContext *avctx;

    int flags;
    const uint8_t *bands1024;
    int num_bands1024;
    const uint8_t *bands128;
    int num_bands128;

    const struct AACPsyModel *model;
    void* model_priv_data;

    float stereo_att;
    int   cutoff;
    void* lp_state;
}AACPsyContext;

typedef struct AACPsyModel {
    const char *name;
    int   (*init)   (AACPsyContext *apc, int elements);
    void  (*window) (AACPsyContext *apc, int16_t *audio, int16_t *la, int tag, int type, ChannelElement *cpe);
    void  (*process)(AACPsyContext *apc, int tag, int type, ChannelElement *cpe);
    void  (*end)    (AACPsyContext *apc);
}AACPsyModel;

/**
 * Initialize psychoacoustic model.
 *
 * @param ctx           model context
 * @param avctx         codec context
 * @param model         model implementation that will be used
 * @param elements      number of channel elements (single channel or channel pair) to handle by  model
 * @param flags         model flags, may be ignored by model if unsupported
 * @param bands1024     scalefactor band lengths for long (1024 samples) frame
 * @param num_bands1024 number of scalefactor bands for long frame
 * @param bands128      scalefactor band lengths for short (128 samples) frame
 * @param num_bands128  number of scalefactor bands for short frame
 *
 * @return zero if successful, a negative value if not
 */
int ff_aac_psy_init(AACPsyContext *ctx, AVCodecContext *avctx,
                    enum AACPsyModelType model, int elements, int flags,
                    const uint8_t *bands1024, int num_bands1024,
                    const uint8_t *bands128,  int num_bands128);

/**
 * Preprocess audio frame in order to compress it better.
 *
 * @param ctx   model context
 * @param audio samples to preprocess
 * @param dest  place to put filtered samples
 * @param tag   number of channel element to analyze
 * @param type  channel element type (e.g. ID_SCE or ID_CPE)
 */
void ff_aac_psy_preprocess(AACPsyContext *ctx, int16_t *audio, int16_t *dest, int tag, int type);

/**
 * Set window sequence and related parameters for channel element.
 *
 * @param ctx   model context
 * @param audio samples for the current frame
 * @param la    lookahead samples (NULL when unavailable)
 * @param tag   number of channel element to analyze
 * @param type  channel element type (e.g. ID_SCE or ID_CPE)
 * @param cpe   pointer to the current channel element
 */
void ff_aac_psy_suggest_window(AACPsyContext *ctx, int16_t *audio, int16_t *la, int tag, int type, ChannelElement *cpe);

/**
 * Perform psychoacoustic analysis and output coefficients in integer form
 * along with scalefactors, M/S flags, etc.
 *
 * @param ctx   model context
 * @param tag   number of channel element to analyze
 * @param type  channel element type (e.g. ID_SCE or ID_CPE)
 * @param cpe   pointer to the current channel element
 */
void ff_aac_psy_analyze(AACPsyContext *ctx, int tag, int type, ChannelElement *cpe);

/**
 * Cleanup model context at the end.
 *
 * @param ctx model context
 */
void ff_aac_psy_end(AACPsyContext *ctx);
#endif /* FFMPEG_AACPSY_H */

