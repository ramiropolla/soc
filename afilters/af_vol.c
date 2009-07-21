/*
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
 * @file libavfilter/af_vol.c
 * null filter. used as an example, or for development
 */


#include <stdio.h>
#include "avfilter.h"

static int vol_filter(AVFilterLink *link, AVFilterBufferRef *sample_ref);
static int query_af_vol_formats(AVFilterContext *ctx);

AVFilter avfilter_af_volume =
{
    .name      = "audio_volume",

    .priv_size = 0,

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO,
                                    .filter_buffer    = vol_filter },
                                  { .name = NULL}},

    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO, },
                                  { .name = NULL}},
    .query_formats = query_af_vol_formats
};


static int query_af_vol_formats(AVFilterContext *ctx)
{
    AVFilterFormats *formats;
    formats = avfilter_make_format_list(2, CODEC_ID_PCM_S16LE,
                                            CODEC_ID_PCM_S16BE);
    avfilter_set_common_formats(ctx,formats);
}

static int vol_filter(AVFilterLink *link, AVFilterBufferRef *sample_ref)
{
    av_log(0,0, "Volume filter\n");
    int i;
    int16_t *data;
    int num_samples = sample_ref->buffer->n_samples;

    data = (int16_t*) sample_ref->buffer->data;
    for (i=0; i < num_samples; i++)
    {
        data[i]  = data[i] * 2;
    }

    return 0;
}
