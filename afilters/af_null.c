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
 * @file libavfilter/vf_null.c
 * null filter. used as an example, or for development
 */


#include <stdio.h>
#include "avfilter.h"
static int filter(AVFilterLink *link, AVFilterBufferRef *sample_ref);

typedef struct
{
    int history[100]; /*just an example */

} af_null_priv_t;

AVFilter avfilter_af_null =
{
    .name      = "audio_null",

    .priv_size = sizeof(af_null_priv_t),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO,
                                    .filter_buffer    = filter },
                                  { .name = NULL}},

    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO, },
                                  { .name = NULL}},
};


static int filter(AVFilterLink *link, AVFilterBufferRef *sample_ref)
{
    av_log(0,0, "Filter buffer\n");
    int num_samples = sample_ref->buffer->n_samples;
    int i;

    int16_t *data;
    data = (int16_t*) sample_ref->buffer->data;
    for (i=0; i < num_samples; i++)
    {
        data[i]  = data[i] +1;
    }

    return 0;
}


