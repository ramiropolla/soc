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
 * null filter
 */

#include <stdio.h>
#include "avfilter.h"


typedef struct
{
    int history[100]; /*just an example */

} af_null_priv_t;


static int start_buf(AVFilterLink *link, AVFilterSamplesRef *sample_ref)
{
    av_log(0,0, "Starting buffer\n");
    return;
}

static int end_buf(AVFilterLink *link, AVFilterSamplesRef *sample_ref)
{
    av_log(0,0, "Ending buffer\n");
    return;
}

AVFilter avfilter_af_null =
{
    .name      = "audio_null",

    .priv_size = sizeof(af_null_priv_t),

    .inputs    = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO,
                                    .start_buffer    = start_buf,
                                    .end_buffer      = end_buf },
                                  { .name = NULL}},

    .outputs   = (AVFilterPad[]) {{ .name            = "default",
                                    .type            = CODEC_TYPE_AUDIO, },
                                  { .name = NULL}},
};
