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
#include "aac.h"
#include "aactab.h"

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

static const uint8_t swb_size_128_96[] = {
    4, 4, 4, 4, 4, 4, 8, 8, 8, 16, 28, 36
};

static const uint8_t swb_size_128_48[] = {
    4, 4, 4, 4, 4, 8, 8, 8, 12, 12, 12, 16, 16, 16
};

static const uint8_t swb_size_128_24[] = {
    4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 12, 12, 16, 16, 20
};

static const uint8_t swb_size_128_16[] = {
    4, 4, 4, 4, 4, 4, 4, 4, 8, 8, 12, 12, 16, 20, 20
};

static const uint8_t swb_size_128_8[] = {
    4, 4, 4, 4, 4, 4, 4, 8, 8, 8, 8, 12, 16, 20, 20
};

static const uint8_t *swb_size_128[] = {
    /* the last entry on the following row is swb_size_128_64 but is a
       duplicate of swb_size_128_96 */
    swb_size_128_96, swb_size_128_96, swb_size_128_96,
    swb_size_128_48, swb_size_128_48, swb_size_128_48,
    swb_size_128_24, swb_size_128_24, swb_size_128_16,
    swb_size_128_16, swb_size_128_16, swb_size_128_8
};

#define CB_UNSIGNED 0x01    ///< coefficients are coded as absolute values
#define CB_PAIRS    0x02    ///< coefficients are grouped into pairs before coding (quads by default)
#define CB_ESCAPE   0x04    ///< codebook allows escapes

/** spectral coefficients codebook information */
static const struct {
    int16_t maxval;         ///< maximum possible value
     int8_t cb_num;         ///< codebook number
    uint8_t flags;          ///< codebook features
} aac_cb_info[] = {
    {    0, -1, CB_UNSIGNED }, // zero codebook
    {    1,  0, 0 },
    {    1,  1, 0 },
    {    2,  2, CB_UNSIGNED },
    {    2,  3, CB_UNSIGNED },
    {    4,  4, CB_PAIRS },
    {    4,  5, CB_PAIRS },
    {    7,  6, CB_PAIRS | CB_UNSIGNED },
    {    7,  7, CB_PAIRS | CB_UNSIGNED },
    {   12,  8, CB_PAIRS | CB_UNSIGNED },
    {   12,  9, CB_PAIRS | CB_UNSIGNED },
    { 8191, 10, CB_PAIRS | CB_UNSIGNED | CB_ESCAPE },
    {   -1, -1, 0 }, // reserved
    {   -1, -1, 0 }, // perceptual noise substitution
    {   -1, -1, 0 }, // intensity out-of-phase
    {   -1, -1, 0 }, // intensity in-phase
};

/** default channel configurations */
static const uint8_t aac_chan_configs[6][5] = {
 {1, ID_SCE},                         // 1 channel  - single channel element
 {1, ID_CPE},                         // 2 channels - channel pair
 {2, ID_SCE, ID_CPE},                 // 3 channels - center + stereo
 {3, ID_SCE, ID_CPE, ID_SCE},         // 4 channels - front center + stereo + back center
 {3, ID_SCE, ID_CPE, ID_CPE},         // 5 channels - front center + stereo + back stereo
 {4, ID_SCE, ID_CPE, ID_CPE, ID_LFE}, // 6 channels - front center + stereo + back stereo + LFE
};

/**
 * AAC encoder context
 */
typedef struct {
    PutBitContext pb;
    MDCTContext mdct1024;                        ///< long (1024 samples) frame transform context
    MDCTContext mdct128;                         ///< short (128 samples) frame transform context
    DSPContext  dsp;
    DECLARE_ALIGNED_16(FFTSample, output[2048]); ///< temporary buffer for MDCT input coefficients
    DECLARE_ALIGNED_16(FFTSample, tmp[1024]);    ///< temporary buffer used by MDCT
    int16_t* samples;                            ///< saved preprocessed input

    int samplerate_index;                        ///< MPEG-4 samplerate index
    const uint8_t *swb_sizes1024;                ///< scalefactor band sizes for long frame
    int swb_num1024;                             ///< number of scalefactor bands for long frame
    const uint8_t *swb_sizes128;                 ///< scalefactor band sizes for short frame
    int swb_num128;                              ///< number of scalefactor bands for short frame

    ChannelElement *cpe;                         ///< channel elements
    AACPsyContext psy;                           ///< psychoacoustic model context
} AACEncContext;

/**
 * Make AAC audio config object.
 * @see 1.6.2.1 "Syntax - AudioSpecificConfig"
 */
static void put_audio_specific_config(AVCodecContext *avctx)
{
    PutBitContext pb;
    AACEncContext *s = avctx->priv_data;

    init_put_bits(&pb, avctx->extradata, avctx->extradata_size*8);
    put_bits(&pb, 5, 2); //object type - AAC-LC
    put_bits(&pb, 4, s->samplerate_index); //sample rate index
    put_bits(&pb, 4, avctx->channels);
    //GASpecificConfig
    put_bits(&pb, 1, 0); //frame length - 1024 samples
    put_bits(&pb, 1, 0); //does not depend on core coder
    put_bits(&pb, 1, 0); //is not extension
    flush_put_bits(&pb);
}

static av_cold int aac_encode_init(AVCodecContext *avctx)
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
    if(avctx->channels > 6){
        av_log(avctx, AV_LOG_ERROR, "Unsupported number of channels: %d\n", avctx->channels);
        return -1;
    }
    s->samplerate_index = i;
    s->swb_sizes1024 = swb_size_1024[i];
    s->swb_num1024 = ff_aac_num_swb_1024[i];
    s->swb_sizes128 = swb_size_128[i];
    s->swb_num128 = ff_aac_num_swb_128[i];

    dsputil_init(&s->dsp, avctx);
    ff_mdct_init(&s->mdct1024, 11, 0);
    ff_mdct_init(&s->mdct128,   8, 0);
    // window init
    ff_kbd_window_init(ff_aac_kbd_long_1024, 4.0, 1024);
    ff_kbd_window_init(ff_aac_kbd_short_128, 6.0, 128);
    ff_sine_window_init(ff_aac_sine_long_1024, 1024);
    ff_sine_window_init(ff_aac_sine_short_128, 128);

    s->samples = av_malloc(2 * 1024 * avctx->channels * sizeof(s->samples[0]));
    s->cpe = av_mallocz(sizeof(ChannelElement) * aac_chan_configs[avctx->channels-1][0]);
    //TODO: psy model selection with some option
    if(ff_aac_psy_init(&s->psy, avctx, AAC_PSY_3GPP, aac_chan_configs[avctx->channels-1][0], 0, s->swb_sizes1024, s->swb_num1024, s->swb_sizes128, s->swb_num128) < 0){
        av_log(avctx, AV_LOG_ERROR, "Cannot initialize selected model.\n");
        return -1;
    }
    avctx->extradata = av_malloc(2);
    avctx->extradata_size = 2;
    put_audio_specific_config(avctx);
    return 0;
}

/**
 * Perform windowing and MDCT.
 */
static void analyze(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, short *audio, int channel)
{
    int i, j, k;
    const float * lwindow = cpe->ch[channel].ics.use_kb_window[0] ? ff_aac_kbd_long_1024 : ff_aac_sine_long_1024;
    const float * swindow = cpe->ch[channel].ics.use_kb_window[0] ? ff_aac_kbd_short_128 : ff_aac_sine_short_128;
    const float * pwindow = cpe->ch[channel].ics.use_kb_window[1] ? ff_aac_kbd_short_128 : ff_aac_sine_short_128;

    if (cpe->ch[channel].ics.window_sequence[0] != EIGHT_SHORT_SEQUENCE) {
        memcpy(s->output, cpe->ch[channel].saved, sizeof(float)*1024);
        if(cpe->ch[channel].ics.window_sequence[0] == LONG_STOP_SEQUENCE){
            memset(s->output, 0, sizeof(s->output[0]) * 448);
            for(i = 448; i < 576; i++)
                s->output[i] = cpe->ch[channel].saved[i] * pwindow[i - 448];
            for(i = 576; i < 704; i++)
                s->output[i] = cpe->ch[channel].saved[i];
        }
        if(cpe->ch[channel].ics.window_sequence[0] != LONG_START_SEQUENCE){
            j = channel;
            for (i = 0; i < 1024; i++, j += avctx->channels){
                s->output[i+1024]         = audio[j] / 512.0 * lwindow[1024 - i - 1];
                cpe->ch[channel].saved[i] = audio[j] / 512.0 * lwindow[i];
            }
        }else{
            j = channel;
            for(i = 0; i < 448; i++, j += avctx->channels)
                s->output[i+1024]         = audio[j] / 512.0;
            for(i = 448; i < 576; i++, j += avctx->channels)
                s->output[i+1024]         = audio[j] / 512.0 * swindow[576 - i - 1];
            memset(s->output+1024+576, 0, sizeof(s->output[0]) * 448);
            j = channel;
            for(i = 0; i < 1024; i++, j += avctx->channels)
                cpe->ch[channel].saved[i] = audio[j] / 512.0;
        }
        ff_mdct_calc(&s->mdct1024, cpe->ch[channel].coeffs, s->output, s->tmp);
    }else{
        j = channel;
        for (k = 0; k < 1024; k += 128) {
            for(i = 448 + k; i < 448 + k + 256; i++)
                s->output[i - 448 - k] = (i < 1024) ? cpe->ch[channel].saved[i] : audio[channel + (i-1024)*avctx->channels] / 512.0;
            s->dsp.vector_fmul        (s->output,     k ?  swindow : pwindow, 128);
            s->dsp.vector_fmul_reverse(s->output+128, s->output+128, swindow, 128);
            ff_mdct_calc(&s->mdct128, cpe->ch[channel].coeffs + k, s->output, s->tmp);
        }
        j = channel;
        for(i = 0; i < 1024; i++, j += avctx->channels)
            cpe->ch[channel].saved[i] = audio[j] / 512.0;
    }
}

/**
 * Encode ics_info element.
 * @see Table 4.6 (syntax of ics_info)
 */
static void put_ics_info(AVCodecContext *avctx, IndividualChannelStream *info)
{
    AACEncContext *s = avctx->priv_data;
    int i;

    put_bits(&s->pb, 1, 0);                // ics_reserved bit
    put_bits(&s->pb, 2, info->window_sequence[0]);
    put_bits(&s->pb, 1, info->use_kb_window[0]);
    if(info->window_sequence[0] != EIGHT_SHORT_SEQUENCE){
        put_bits(&s->pb, 6, info->max_sfb);
        put_bits(&s->pb, 1, 0);            // no prediction
    }else{
        put_bits(&s->pb, 4, info->max_sfb);
        for(i = 1; i < info->num_windows; i++)
            put_bits(&s->pb, 1, info->group_len[i]);
    }
}

/**
 * Encode MS data.
 * @see 4.6.8.1 "Joint Coding - M/S Stereo"
 */
static void encode_ms_info(PutBitContext *pb, ChannelElement *cpe)
{
    int i, w;

    put_bits(pb, 2, cpe->ms.present);
    if(cpe->ms.present == 1)
        for(w = 0; w < cpe->ch[0].ics.num_windows; w++){
            if(cpe->ch[0].ics.group_len[w]) continue;
            for(i = 0; i < cpe->ch[0].ics.max_sfb; i++)
                put_bits(pb, 1, cpe->ms.mask[w][i]);
        }
}

/**
 * Scan scalefactor band and determine optimal codebook for it.
 *
 * @param s       encoder context
 * @param cpe     channel element
 * @param channel channel number inside channel pair
 * @param win     window group start number
 * @param band    scalefactor band to analyze
 * @param start   scalefactor band position in spectral coefficients
 * @param size    scalefactor band size
 */
static int determine_section_info(AACEncContext *s, ChannelElement *cpe, int channel, int win, int band, int start, int size)
{
    int i, j, w;
    int maxval, sign;
    int score, best, cb, bestcb, dim, idx;

    maxval = 0;
    sign = 0;
    w = win;
    do{
        for(i = start + (w-win)*128; i < start + (w-win)*128 + size; i++){
            maxval = FFMAX(maxval, FFABS(cpe->ch[channel].icoefs[i]));
            if(cpe->ch[channel].icoefs[i] < 0) sign = 1;
        }
        w++;
    }while(w < cpe->ch[channel].ics.num_windows && cpe->ch[channel].ics.group_len[w]);

    if(maxval > 12) return 11;
    if(!maxval) return 0;

    for(cb = 0; cb < 12; cb++)
        if(aac_cb_info[cb].maxval >= maxval)
            break;
    best = 9999;
    bestcb = 11;
    for(; cb < 12; cb++){
        score = 0;
        dim = (aac_cb_info[cb].flags & CB_PAIRS) ? 2 : 4;
        if(!band || cpe->ch[channel].band_type[win][band - 1] != cb)
            score += 9; //that's for new codebook entry
        w = win;
        if(aac_cb_info[cb].flags & CB_UNSIGNED){
            do{
                for(i = start + (w-win)*128; i < start + (w-win)*128 + size; i += dim){
                    idx = 0;
                    for(j = 0; j < dim; j++)
                        idx = idx * aac_cb_info[cb].maxval + FFABS(cpe->ch[channel].icoefs[i+j]);
                    score += ff_aac_spectral_bits[aac_cb_info[cb].cb_num][idx];
                    for(j = 0; j < dim; j++)
                        if(cpe->ch[channel].icoefs[i+j])
                            score++;
                }
                w++;
            }while(w < cpe->ch[channel].ics.num_windows && cpe->ch[channel].ics.group_len[w]);
        }else{
            do{
                for(i = start + (w-win)*128; i < start + (w-win)*128 + size; i += dim){
                    idx = 0;
                    for(j = 0; j < dim; j++)
                        idx = idx * (aac_cb_info[cb].maxval*2 + 1) + cpe->ch[channel].icoefs[i+j] + aac_cb_info[cb].maxval;
                    score += ff_aac_spectral_bits[aac_cb_info[cb].cb_num][idx];
                }
                w++;
            }while(w < cpe->ch[channel].ics.num_windows && cpe->ch[channel].ics.group_len[w]);
        }
        if(score < best){
            best = score;
            bestcb = cb;
        }
    }
    return bestcb;
}

/**
 * Encode one scalefactor band with selected codebook.
 */
static void encode_codebook(AACEncContext *s, ChannelElement *cpe, int channel, int start, int size, int cb)
{
    const uint8_t *bits = ff_aac_spectral_bits[aac_cb_info[cb].cb_num];
    const uint16_t *codes = ff_aac_spectral_codes[aac_cb_info[cb].cb_num];
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
                idx = idx * (aac_cb_info[cb].maxval + 1) + FFABS(cpe->ch[channel].icoefs[i+j]);
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

/**
 * Encode scalefactor band coding type.
 */
static void encode_section_data(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel)
{
    int i, w;
    int bits = cpe->ch[channel].ics.num_windows == 1 ? 5 : 3;
    int esc = (1 << bits) - 1;
    int count;

    for(w = 0; w < cpe->ch[channel].ics.num_windows; w++){
        if(cpe->ch[channel].ics.group_len[w]) continue;
        count = 0;
        for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
            if(!i || cpe->ch[channel].band_type[w][i] != cpe->ch[channel].band_type[w][i-1]){
                if(count){
                    while(count >= esc){
                        put_bits(&s->pb, bits, esc);
                        count -= esc;
                    }
                    put_bits(&s->pb, bits, count);
                }
                put_bits(&s->pb, 4, cpe->ch[channel].band_type[w][i]);
                count = 1;
            }else
                count++;
        }
        if(count){
            while(count >= esc){
                put_bits(&s->pb, bits, esc);
                count -= esc;
            }
            put_bits(&s->pb, bits, count);
        }
    }
}

/**
 * Encode scalefactors.
 */
static void encode_scale_factor_data(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel)
{
    int off = cpe->ch[channel].mixing_gain, diff;
    int i, w;

    for(w = 0; w < cpe->ch[channel].ics.num_windows; w++){
        if(cpe->ch[channel].ics.group_len[w]) continue;
        for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
            if(!cpe->ch[channel].zeroes[w][i]){
                diff = cpe->ch[channel].sf_idx[w][i] - off + SCALE_DIFF_ZERO;
                if(diff < 0 || diff > 120) av_log(avctx, AV_LOG_ERROR, "Scalefactor difference is too big to be coded\n");
                off = cpe->ch[channel].sf_idx[w][i];
                put_bits(&s->pb, ff_aac_scalefactor_bits[diff], ff_aac_scalefactor_code[diff]);
            }
        }
    }
}

/**
 * Encode pulse data.
 */
static void encode_pulse_data(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel)
{
    int i;

    put_bits(&s->pb, 1, cpe->ch[channel].pulse.present);
    if(!cpe->ch[channel].pulse.present) return;

    put_bits(&s->pb, 2, cpe->ch[channel].pulse.num_pulse - 1);
    put_bits(&s->pb, 6, cpe->ch[channel].pulse.start);
    for(i = 0; i < cpe->ch[channel].pulse.num_pulse; i++){
        put_bits(&s->pb, 5, cpe->ch[channel].pulse.offset[i]);
        put_bits(&s->pb, 4, cpe->ch[channel].pulse.amp[i]);
    }
}

/**
 * Encode temporal noise shaping data.
 */
static void encode_tns_data(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel)
{
    int i, w;

    put_bits(&s->pb, 1, cpe->ch[channel].tns.present);
    if(!cpe->ch[channel].tns.present) return;
    if(cpe->ch[channel].ics.window_sequence[0] == EIGHT_SHORT_SEQUENCE){
        for(w = 0; w < cpe->ch[channel].ics.num_windows; w++){
            put_bits(&s->pb, 1, cpe->ch[channel].tns.n_filt[w]);
            if(!cpe->ch[channel].tns.n_filt[w]) continue;
            put_bits(&s->pb, 1, cpe->ch[channel].tns.coef_res[w] - 3);
            put_bits(&s->pb, 4, cpe->ch[channel].tns.length[w][0]);
            put_bits(&s->pb, 3, cpe->ch[channel].tns.order[w][0]);
            if(cpe->ch[channel].tns.order[w][0]){
                put_bits(&s->pb, 1, cpe->ch[channel].tns.direction[w][0]);
                put_bits(&s->pb, 1, cpe->ch[channel].tns.coef_compress[w][0]);
                for(i = 0; i < cpe->ch[channel].tns.order[w][0]; i++)
                     put_bits(&s->pb, cpe->ch[channel].tns.coef_len[w][0], cpe->ch[channel].tns.coef[w][0][i]);
            }
        }
    }else{
        put_bits(&s->pb, 1, cpe->ch[channel].tns.n_filt[0]);
        if(!cpe->ch[channel].tns.n_filt[0]) return;
        put_bits(&s->pb, 1, cpe->ch[channel].tns.coef_res[0] - 3);
        for(w = 0; w < cpe->ch[channel].tns.n_filt[0]; w++){
            put_bits(&s->pb, 6, cpe->ch[channel].tns.length[0][w]);
            put_bits(&s->pb, 5, cpe->ch[channel].tns.order[0][w]);
            if(cpe->ch[channel].tns.order[0][w]){
                put_bits(&s->pb, 1, cpe->ch[channel].tns.direction[0][w]);
                put_bits(&s->pb, 1, cpe->ch[channel].tns.coef_compress[0][w]);
                for(i = 0; i < cpe->ch[channel].tns.order[0][w]; i++)
                     put_bits(&s->pb, cpe->ch[channel].tns.coef_len[0][w], cpe->ch[channel].tns.coef[0][w][i]);
            }
        }
    }
}

/**
 * Encode spectral coefficients processed by psychoacoustic model.
 */
static void encode_spectral_data(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel)
{
    int start, i, w, w2;

    for(w = 0; w < cpe->ch[channel].ics.num_windows; w++){
        if(cpe->ch[channel].ics.group_len[w]) continue;
        start = 0;
        for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
            if(cpe->ch[channel].zeroes[w][i]){
                start += cpe->ch[channel].ics.swb_sizes[i];
                continue;
            }
            w2 = w;
            do{
                encode_codebook(s, cpe, channel, start + w2*128, cpe->ch[channel].ics.swb_sizes[i], cpe->ch[channel].band_type[w][i]);
                w2++;
            }while(w2 < cpe->ch[channel].ics.num_windows && cpe->ch[channel].ics.group_len[w2]);
            start += cpe->ch[channel].ics.swb_sizes[i];
        }
    }
}

/**
 * Encode one channel of audio data.
 */
static int encode_individual_channel(AVCodecContext *avctx, ChannelElement *cpe, int channel)
{
    AACEncContext *s = avctx->priv_data;
    int i, g, w;

    for(w = 0; w < cpe->ch[channel].ics.num_windows; w++){
        i = w << 7;
        if(cpe->ch[channel].ics.group_len[w]) continue;
        for(g = 0; g < cpe->ch[channel].ics.max_sfb; g++){
            if(!cpe->ch[channel].zeroes[w][g]){
                cpe->ch[channel].band_type[w][g] = determine_section_info(s, cpe, channel, w, g, i, cpe->ch[channel].ics.swb_sizes[g]);
                cpe->ch[channel].zeroes[w][g] = !cpe->ch[channel].band_type[w][g];
            }else
                cpe->ch[channel].band_type[w][g] = 0;
            i += cpe->ch[channel].ics.swb_sizes[g];
        }
    }

    put_bits(&s->pb, 8, cpe->ch[channel].mixing_gain); //global gain
    if(!cpe->common_window) put_ics_info(avctx, &cpe->ch[channel].ics);
    encode_section_data(avctx, s, cpe, channel);
    encode_scale_factor_data(avctx, s, cpe, channel);
    encode_pulse_data(avctx, s, cpe, channel);
    encode_tns_data(avctx, s, cpe, channel);
    put_bits(&s->pb, 1, 0); //ssr
    encode_spectral_data(avctx, s, cpe, channel);
    return 0;
}

/**
 * Write some auxiliary information about the created AAC file.
 */
static void put_bitstream_info(AVCodecContext *avctx, AACEncContext *s, const char *name)
{
    int i, namelen, padbits;

    namelen = strlen(name) + 2;
    put_bits(&s->pb, 3, ID_FIL);
    put_bits(&s->pb, 4, FFMIN(namelen, 15));
    if(namelen >= 15)
        put_bits(&s->pb, 8, namelen - 16);
    put_bits(&s->pb, 4, 0); //extension type - filler
    padbits = 8 - (put_bits_count(&s->pb) & 7);
    align_put_bits(&s->pb);
    for(i = 0; i < namelen - 2; i++)
        put_bits(&s->pb, 8, name[i]);
    put_bits(&s->pb, 12 - padbits, 0);
}

static int aac_encode_frame(AVCodecContext *avctx,
                            uint8_t *frame, int buf_size, void *data)
{
    AACEncContext *s = avctx->priv_data;
    int16_t *samples = s->samples, *samples2, *la;
    ChannelElement *cpe;
    int i, j, chans, tag, start_ch;
    const uint8_t *chan_map = aac_chan_configs[avctx->channels-1];
    int chan_el_counter[4];

    if(data){
        start_ch = 0;
        samples2 = s->samples + 1024 * avctx->channels;
        for(i = 0; i < chan_map[0]; i++){
            tag = chan_map[i+1];
            chans = tag == ID_CPE ? 2 : 1;
            ff_aac_psy_preprocess(&s->psy, (uint16_t*)data + start_ch, samples2 + start_ch, i, tag);
            start_ch += chans;
        }
    }
    if(!avctx->frame_number){
        memmove(s->samples, s->samples + 1024 * avctx->channels, 1024 * avctx->channels * sizeof(s->samples[0]));
        return 0;
    }

    init_put_bits(&s->pb, frame, buf_size*8);
    if(avctx->frame_number==1 && !(avctx->flags & CODEC_FLAG_BITEXACT)){
        put_bitstream_info(avctx, s, LIBAVCODEC_IDENT);
    }
    start_ch = 0;
    memset(chan_el_counter, 0, sizeof(chan_el_counter));
    for(i = 0; i < chan_map[0]; i++){
        tag = chan_map[i+1];
        chans = tag == ID_CPE ? 2 : 1;
        cpe = &s->cpe[i];
        samples2 = samples + start_ch;
        la = samples2 + 1024 * avctx->channels + start_ch;
        if(!data) la = NULL;
        ff_aac_psy_suggest_window(&s->psy, samples2, la, i, tag, cpe);
        for(j = 0; j < chans; j++){
            analyze(avctx, s, cpe, samples2, j);
        }
        ff_aac_psy_analyze(&s->psy, i, tag, cpe);
        put_bits(&s->pb, 3, tag);
        put_bits(&s->pb, 4, chan_el_counter[tag]++);
        if(chans == 2){
            put_bits(&s->pb, 1, cpe->common_window);
            if(cpe->common_window){
                put_ics_info(avctx, &cpe->ch[0].ics);
                encode_ms_info(&s->pb, cpe);
            }
        }
        for(j = 0; j < chans; j++){
            encode_individual_channel(avctx, cpe, j);
        }
        start_ch += chans;
    }

    put_bits(&s->pb, 3, ID_END);
    flush_put_bits(&s->pb);
    avctx->frame_bits = put_bits_count(&s->pb);

    memmove(s->samples, s->samples + 1024 * avctx->channels, 1024 * avctx->channels * sizeof(s->samples[0]));
    return put_bits_count(&s->pb)>>3;
}

static av_cold int aac_encode_end(AVCodecContext *avctx)
{
    AACEncContext *s = avctx->priv_data;

    ff_mdct_end(&s->mdct1024);
    ff_mdct_end(&s->mdct128);
    ff_aac_psy_end(&s->psy);
    av_freep(&s->samples);
    av_freep(&s->cpe);
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
    .capabilities = CODEC_CAP_SMALL_LAST_FRAME,
    .sample_fmts = (enum SampleFormat[]){SAMPLE_FMT_S16,SAMPLE_FMT_NONE},
    .long_name = NULL_IF_CONFIG_SMALL("Advanced Audio Coding"),
};
