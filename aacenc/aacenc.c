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

#include "aacpsy.h"

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
    DECLARE_ALIGNED_16(FFTSample, tmp[1024]);

    int samplerate_index;
    uint8_t *swb_sizes;
    int swb_num;
    cpe_struct cpe;
    AACPsyContext psy;
} AACEncContext;

#define SCALE_ONE_POS   140
#define SCALE_MAX_POS   255
#define SCALE_MAX_DIFF   60
#define SCALE_DIFF_ZERO  60

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

    ff_mdct_init(&s->mdct, 11, 0);
    // window init
    ff_kbd_window_init(s->kbd_long_1024, 4.0, 1024);

    ff_aac_psy_init(&s->psy, avctx, AAC_PSY_NULL, 0, s->swb_sizes, s->swb_num);
    avctx->extradata = av_malloc(2);
    avctx->extradata_size = 2;
    put_audio_specific_config(avctx);

    return 0;
}

static void analyze(AVCodecContext *avctx, AACEncContext *s, cpe_struct *cpe, short *audio, int channel)
{
    int i, j;

    // perform MDCT
    memcpy(s->output, cpe->ch[channel].saved, sizeof(float)*1024);
    j = channel;
    for (i = 0; i < 1024; i++, j += avctx->channels){
        s->output[i+1024]         = audio[j] / 512 * s->kbd_long_1024[1024 - i - 1];
        cpe->ch[channel].saved[i] = audio[j] / 512 * s->kbd_long_1024[i];
    }
    ff_mdct_calc(&s->mdct, cpe->ch[channel].coeffs, s->output, s->tmp);
    //convert coefficients into form used by AAC
    for(i = 0; i < 1024; i++)
        cpe->ch[channel].coeffs[i] = -copysignf(pow(fabsf(cpe->ch[channel].coeffs[i]), 0.75f), cpe->ch[channel].coeffs[i]);
}

/**
 * Encode ics_info element.
 * @see Table 4.6
 */
static void put_ics_info(AVCodecContext *avctx, ics_struct *info)
{
    AACEncContext *s = avctx->priv_data;

    put_bits(&s->pb, 1, 0);                // ics_reserved bit
    put_bits(&s->pb, 2, info->window_sequence);
    put_bits(&s->pb, 1, info->window_shape);
    put_bits(&s->pb, 6, info->max_sfb);
    put_bits(&s->pb, 1, 0);                // no prediction
}

/**
 * Scan spectral band and determine optimal codebook for it.
 */
static int determine_section_info(AACEncContext *s, cpe_struct *cpe, int channel, int start, int size)
{
    int i;
    int maxval, sign;

    maxval = 0;
    sign = 0;
    for(i = start; i < start + size; i++){
        maxval = FFMAX(maxval, FFABS(cpe->ch[channel].icoefs[i]));
        if(cpe->ch[channel].icoefs[i] < 0) sign = 1;
    }

    ///TODO: better decision
    if(!maxval) return 0; //zero codebook
    if(maxval == 1) return 2;
    return 11; //escape codebook
}

static void encode_codebook(AACEncContext *s, cpe_struct *cpe, int channel, int start, int size, int cb)
{
    const uint8_t *bits = aac_cb_info[cb].bits;
    const uint16_t *codes = aac_cb_info[cb].codes;
    const int dim = (aac_cb_info[cb].flags & CB_PAIRS) ? 2 : 4;
    int i, j, idx;

    if(!bits || !codes) return;

    //TODO: factorize?
    if(aac_cb_info[cb].flags & CB_ESCAPE){
        for(i = start; i < start + size; i += dim){
            idx = 0;
            for(j = 0; j < dim; j++)
                idx = idx*17 + FFMIN(FFABS(cpe->ch[channel].icoefs[i+j]), 16);
            put_bits(&s->pb, bits[idx], codes[idx]);
            //output signs
            for(j = 0; j < dim; j++)
                if(cpe->ch[channel].icoefs[i+j])
                    put_bits(&s->pb, 1, cpe->ch[channel].icoefs[i+j] < 0);
            //output escape values
            for(j = 0; j < dim; j++)
                if(FFABS(cpe->ch[channel].icoefs[i+j]) > 15){
                    int l = av_log2(FFABS(cpe->ch[channel].icoefs[i+j]));

                    put_bits(&s->pb, l - 4 + 1, (1 << (l - 4 + 1)) - 2);
                    put_bits(&s->pb, l, FFABS(cpe->ch[channel].icoefs[i+j]) & ((1 << l) - 1));
                }
        }
    }else if(aac_cb_info[cb].flags & CB_UNSIGNED){
        for(i = start; i < start + size; i += dim){
            idx = 0;
            for(j = 0; j < dim; j++)
                idx = idx * aac_cb_info[cb].maxval + FFABS(cpe->ch[channel].icoefs[i+j]);
            put_bits(&s->pb, bits[idx], codes[idx]);
            //output signs
            for(j = 0; j < dim; j++)
                if(cpe->ch[channel].icoefs[i+j])
                    put_bits(&s->pb, 1, cpe->ch[channel].icoefs[i+j] < 0);
        }
    }else{
        for(i = start; i < start + size; i += dim){
            idx = 0;
            for(j = 0; j < dim; j++)
                idx = idx * (aac_cb_info[cb].maxval*2 + 1) + cpe->ch[channel].icoefs[i+j] + aac_cb_info[cb].maxval;
            put_bits(&s->pb, bits[idx], codes[idx]);
        }
    }
}

static void encode_section_data(AVCodecContext *avctx, AACEncContext *s, cpe_struct *cpe, int channel)
{
    int i;
    int bits = 5; //for long window
    int count = 0;

    for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
        if(!i || cpe->ch[channel].cb[0][i] != cpe->ch[channel].cb[0][i-1]){
            if(count){
                while(count >= (1 << bits) - 1){
                    put_bits(&s->pb, bits, (1 << bits) - 1);
                    count -= (1 << bits) - 1;
                }
                put_bits(&s->pb, bits, count);
            }
            put_bits(&s->pb, 4, cpe->ch[channel].cb[0][i]);
            count = 1;
        }else
            count++;
    }
    if(count){
        while(count >= (1 << bits) - 1){
            put_bits(&s->pb, bits, (1 << bits) - 1);
            count -= (1 << bits) - 1;
        }
        put_bits(&s->pb, bits, count);
    }
}

static void encode_scale_factor_data(AVCodecContext *avctx, AACEncContext *s, cpe_struct *cpe, int channel)
{
    int off = cpe->ch[channel].gain, diff;
    int i;

    for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
        if(!cpe->ch[channel].zeroes[i]){
            diff = cpe->ch[channel].sf_idx[i] - off + SCALE_DIFF_ZERO;
            off = cpe->ch[channel].sf_idx[i];
            put_bits(&s->pb, bits[diff], code[diff]);
        }
    }
}

static void encode_spectral_data(AVCodecContext *avctx, AACEncContext *s, cpe_struct *cpe, int channel)
{
    int start = 0, i;

    for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
        if(!cpe->ch[channel].zeroes[i])
            encode_codebook(s, cpe, channel, start, s->swb_sizes[i], cpe->ch[channel].cb[0][i]);
        start += s->swb_sizes[i];
    }
}

/**
 * Encode one channel of audio data.
 */
static int encode_individual_channel(AVCodecContext *avctx, cpe_struct *cpe, int channel)
{
    AACEncContext *s = avctx->priv_data;
    int i, j, g = 0;

    i = 0;
    while(i < 1024){
        if(!cpe->ch[channel].zeroes[g]){
            cpe->ch[channel].cb[0][g] = determine_section_info(s, cpe, channel, i, s->swb_sizes[g]);
            cpe->ch[channel].zeroes[g] = !cpe->ch[channel].cb[0][g];
        }else
            cpe->ch[channel].cb[0][g] = 0;
        i += s->swb_sizes[g];
        g++;
    }

    put_bits(&s->pb, 8, cpe->ch[channel].gain); //global gain
    if(!cpe->common_window) put_ics_info(avctx, &cpe->ch[channel].ics);
    encode_section_data(avctx, s, cpe, channel);
    encode_scale_factor_data(avctx, s, cpe,channel);
    put_bits(&s->pb, 1, 0); //pulse
    put_bits(&s->pb, 1, 0); //tns
    put_bits(&s->pb, 1, 0); //ssr
    encode_spectral_data(avctx, s, cpe, channel);
    return 0;
}

static int aac_encode_frame(AVCodecContext *avctx,
                            uint8_t *frame, int buf_size, void *data)
{
    AACEncContext *s = avctx->priv_data;
    int16_t *samples = data;

    ff_aac_psy_suggest_window(&s->psy, samples, 0, &s->cpe);

    analyze(avctx, s, &s->cpe, samples, 0);
    if(avctx->channels > 1)
        analyze(avctx, s, &s->cpe, samples, 1);

    ff_aac_psy_analyze(&s->psy, samples, 0, &s->cpe);
    if(avctx->channels > 1){
        s->cpe.common_window = s->cpe.ch[0].ics.window_shape == s->cpe.ch[1].ics.window_shape;
        if(s->cpe.common_window){
            s->cpe.ch[0].ics.max_sfb = FFMAX(s->cpe.ch[0].ics.max_sfb, s->cpe.ch[1].ics.max_sfb);
            s->cpe.ch[1].ics.max_sfb = s->cpe.ch[0].ics.max_sfb;
        }
    }

    init_put_bits(&s->pb, frame, buf_size*8);
    //output encoded
    switch(avctx->channels){
    case 1:
        put_bits(&s->pb, 3, ID_SCE);
        put_bits(&s->pb, 4, 0); //tag
        encode_individual_channel(avctx, &s->cpe, 0);
        break;
    case 2:
        put_bits(&s->pb, 3, ID_CPE);
        put_bits(&s->pb, 4, 0); //tag
        put_bits(&s->pb, 1, s->cpe.common_window);
        if(s->cpe.common_window){
            put_ics_info(avctx, &s->cpe.ch[0].ics);
            put_bits(&s->pb, 2, 0); //no MS mode for now
        }
        encode_individual_channel(avctx, &s->cpe, 0);
        encode_individual_channel(avctx, &s->cpe, 1);
        break;
    default:
        av_log(NULL,0,"?");
    }

    put_bits(&s->pb, 3, ID_END);
    flush_put_bits(&s->pb);
    return put_bits_count(&s->pb)>>3;
}

static int aac_encode_end(AVCodecContext *avctx)
{
    AACEncContext *s = avctx->priv_data;

    ff_mdct_end(&s->mdct);
    ff_aac_psy_end(&s->psy);
    return 0;
}

AVCodec aac_encoder = {
    "aac",
    CODEC_TYPE_AUDIO,
    CODEC_ID_AAC,
    sizeof(AACEncContext),
    aac_encode_init,
    aac_encode_frame,
    aac_encode_end,
};
