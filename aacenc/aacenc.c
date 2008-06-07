/*
 * AAC encoder
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

/**
 * @file aacenc.c
 * AAC encoder
 */

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"
#include "mpeg4audio.h"

// XXX: borrowed from aac.c, move to some header eventually

#include "aactab.h"
/**
 * IDs for raw_data_block
 */
enum {
    ID_SCE = 0x0,
    ID_CPE,
    ID_CCE,
    ID_LFE,
    ID_DSE,
    ID_PCE,
    ID_FIL,
    ID_END
};

static const uint8_t swb_size_1024_96[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8,
    12, 12, 12, 12, 12, 16, 16, 24, 28, 36, 44,
    64, 64, 64, 64, 64, 64, 64, 64, 64, 64, 64
};

static const uint8_t swb_size_1024_64[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8,
    12, 12, 12, 16, 16, 16, 20, 24, 24, 28, 36,
    40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40, 40
};

static const uint8_t swb_size_1024_48[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8,
    12, 12, 12, 12, 16, 16, 20, 20, 24, 24, 28, 28,
    32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
    96
};

static const uint8_t swb_size_1024_32[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8,
    12, 12, 12, 12, 16, 16, 20, 20, 24, 24, 28, 28,
    32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32
};

static const uint8_t swb_size_1024_24[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    12, 12, 12, 12, 16, 16, 16, 20, 20, 24, 24, 28, 28,
    32, 36, 36, 40, 44, 48, 52, 52, 64, 64, 64, 64, 64
};

static const uint8_t swb_size_1024_16[] = {
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    12, 12, 12, 12, 12, 12, 12, 12, 12, 16, 16, 16, 16, 20, 20, 20, 24, 24, 28, 28,
    32, 36, 40, 40, 44, 48, 52, 56, 60, 64, 64, 64
};

static const uint8_t swb_size_1024_8[] = {
    12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
    16, 16, 16, 16, 16, 16, 16, 20, 20, 20, 20, 24, 24, 24, 28, 28,
    32, 36, 36, 40, 44, 48, 52, 56, 60, 64, 80
};

static const uint8_t *swb_size_1024[] = {
    swb_size_1024_96, swb_size_1024_96, swb_size_1024_64,
    swb_size_1024_48, swb_size_1024_48, swb_size_1024_32,
    swb_size_1024_24, swb_size_1024_24, swb_size_1024_16,
    swb_size_1024_16, swb_size_1024_16, swb_size_1024_8
};

#define MAX_SWB_SIZE  51

//borrowed data ends here

#define CB_UNSIGNED 0x01
#define CB_PAIRS    0x02
#define CB_ESCAPE   0x04

/** Codebook information */
static const struct {
    int16_t maxval;         ///< maximum possible value

    const uint8_t  *bits;   ///< codeword lengths
    const uint16_t *codes;  ///< codewords

    uint8_t flags;          ///< codebook features
} aac_cb_info[] = {
    {    0, NULL  , NULL  , CB_UNSIGNED }, // zero codebook
    {    1, bits1 , code1 , 0 },
    {    1, bits2 , code2 , 0 },
    {    2, bits3 , code3 , CB_UNSIGNED },
    {    2, bits4 , code4 , CB_UNSIGNED },
    {    4, bits5 , code5 , CB_PAIRS },
    {    4, bits6 , code6 , CB_PAIRS },
    {    7, bits7 , code7 , CB_PAIRS | CB_UNSIGNED },
    {    7, bits8 , code8 , CB_PAIRS | CB_UNSIGNED },
    {   12, bits9 , code9 , CB_PAIRS | CB_UNSIGNED },
    {   12, bits10, code10, CB_PAIRS | CB_UNSIGNED },
    { 8191, bits11, code11, CB_PAIRS | CB_UNSIGNED | CB_ESCAPE },
    {   -1, NULL  , NULL  , 0 }, // reserved
    {   -1, NULL  , NULL  , 0 }, // perceptual noise substitution
    {   -1, NULL  , NULL  , 0 }, // intensity out-of-phase
    {   -1, NULL  , NULL  , 0 }, // intensity in-phase
};

typedef struct {
    PutBitContext pb;
    MDCTContext mdct;
    DECLARE_ALIGNED_16(float, kbd_long_1024[1024]);
    DECLARE_ALIGNED_16(FFTSample, output[2048]);
    DECLARE_ALIGNED_16(FFTSample, frame_out[2][2048]);
    DECLARE_ALIGNED_16(FFTSample, coefs[2][1024]);
    DECLARE_ALIGNED_16(FFTSample, tmp[1024]);
    DECLARE_ALIGNED_16(int, icoefs[2][1024]);

    int samplerate_index;
    uint8_t *swb_sizes;
    int swb_num;
    int coded_swb_num;
} AACEncContext;

/**
 * Make AAC audio config object.
 * @see 1.6.2.1
 */
static void put_audio_specific_config(AVCodecContext *avctx)
{
    PutBitContext pb;
    AACEncContext *s = avctx->priv_data;

    init_put_bits(&pb, avctx->extradata, avctx->extradata_size*8);
    put_bits(&pb, 5, 2); //object type - AAC-LC
    put_bits(&pb, 4, s->samplerate_index); //sample rate index
    put_bits(&pb, 4, avctx->channels); //channel config - stereo
    //GASpecificConfig
    put_bits(&pb, 1, 0); //frame length - 1024 samples
    put_bits(&pb, 1, 0); //does not depend on core coder
    put_bits(&pb, 1, 0); //is not extension
    flush_put_bits(&pb);
}

static int aac_encode_init(AVCodecContext *avctx)
{
    AACEncContext *s = avctx->priv_data;
    int i;

    avctx->frame_size = 1024;

    for(i = 0; i < 16; i++)
        if(avctx->sample_rate == ff_mpeg4audio_sample_rates[i])
            break;
    if(i == 16){
        av_log(avctx, AV_LOG_ERROR, "Unsupported sample rate %d\n", avctx->sample_rate);
        return -1;
    }
    s->samplerate_index = i;
    s->swb_sizes = swb_size_1024[i];
    s->swb_num = num_swb_1024[i];

    ff_mdct_init(&s->mdct, 11, 1);
    // window init
    ff_kbd_window_init(s->kbd_long_1024, 4.0, 1024);

    avctx->extradata = av_malloc(2);
    avctx->extradata_size = 2;
    put_audio_specific_config(avctx);
    return 0;
}

/* BIG FAT TODO! */
/* for now it just converts spectra to integer form */
static void apply_psychoacoustics(AVCodecContext *avctx, int channel)
{
    AACEncContext *s = avctx->priv_data;
    int i;

    for(i = 0; i < 1024; i++)
        s->icoefs[channel][i] = (int)s->coefs[channel][i];
}

static void analyze(AVCodecContext *avctx, AACEncContext *s, short *audio, int channel)
{
    int i, j;

    // perform MDCT
    memcpy(s->output, s->frame_out[channel], sizeof(float)*1024);
    j = channel;
    for (i = 0; i < 1024; i++, j += avctx->channels){
        s->output[i+1024]        = audio[j] / 512 * s->kbd_long_1024[1024 - i - 1];
        s->frame_out[channel][i] = audio[j] / 512 * s->kbd_long_1024[i];
    }
    ff_mdct_calc(&s->mdct, s->coefs[channel], s->output, s->tmp);

    apply_psychoacoustics(avctx, channel);
}

static int aac_encode_frame(AVCodecContext *avctx,
                            uint8_t *frame, int buf_size, void *data)
{
    int i,k,channel;
    AACEncContext *s = avctx->priv_data;
    int16_t *samples = data;

    init_put_bits(&s->pb, frame, buf_size*8);
    put_bits(&s->pb, 8, 0xAB);

    flush_put_bits(&s->pb);
    return put_bits_count(&s->pb)>>3;
}

AVCodec aac_encoder = {
    "aac",
    CODEC_TYPE_AUDIO,
    CODEC_ID_AAC,
    sizeof(AACEncContext),
    aac_encode_init,
    aac_encode_frame,
};
