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
#include "aac.h"

enum AACPsyModelType{
    AAC_PSY_NULL,              // do nothing on frequencies
    AAC_PSY_NULL8,             // do nothing on frequencies but work with short windows
    AAC_PSY_3GPP,              // model following recommendations from 3GPP TS 26.403

    AAC_NB_PSY_MODELS
};

/**
 * context used by psychoacoustic model
 */
typedef struct AACPsyContext {
    AVCodecContext *avctx;
    DSPContext dsp;

    int flags;
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
    int   (*init)   (AACPsyContext *apc, int elements);
    void  (*window) (AACPsyContext *apc, int16_t *audio, int16_t *la, int tag, int type, ChannelElement *cpe);
    void  (*process)(AACPsyContext *apc, int tag, int type, ChannelElement *cpe);
    void  (*end)    (AACPsyContext *apc);
}AACPsyModel;

int ff_aac_psy_init(AACPsyContext *ctx, AVCodecContext *avctx,
                    enum AACPsyModelType model, int elements, int flags,
                    const uint8_t *bands1024, int num_bands1024,
                    const uint8_t *bands128,  int num_bands128);
void ff_aac_psy_suggest_window(AACPsyContext *ctx, int16_t *audio, int16_t *la, int tag, int type, ChannelElement *cpe);
void ff_aac_psy_analyze(AACPsyContext *ctx, int tag, int type, ChannelElement *cpe);
void ff_aac_psy_end(AACPsyContext *ctx);
#endif /* FFMPEG_AACPSY_H */

