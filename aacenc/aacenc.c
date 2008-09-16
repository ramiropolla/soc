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

/***********************************
 *              TODOs:
 * speedup quantizer selection
 * add sane pulse detection
 * add temporal noise shaping
 ***********************************/

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"
#include "mpeg4audio.h"

#include "aac.h"
#include "aactab.h"

#include "psymodel.h"

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

/** bits needed to code codebook run value for long windows */
static const uint8_t run_value_bits_long[64] = {
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,
     5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5,  5, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
    10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 15
};

/** bits needed to code codebook run value for short windows */
static const uint8_t run_value_bits_short[16] = {
    3, 3, 3, 3, 3, 3, 3, 6, 6, 6, 6, 6, 6, 6, 6, 9
};

static const uint8_t* run_value_bits[2] = {
    run_value_bits_long, run_value_bits_short
};

/** default channel configurations */
static const uint8_t aac_chan_configs[6][5] = {
 {1, TYPE_SCE},                               // 1 channel  - single channel element
 {1, TYPE_CPE},                               // 2 channels - channel pair
 {2, TYPE_SCE, TYPE_CPE},                     // 3 channels - center + stereo
 {3, TYPE_SCE, TYPE_CPE, TYPE_SCE},           // 4 channels - front center + stereo + back center
 {3, TYPE_SCE, TYPE_CPE, TYPE_CPE},           // 5 channels - front center + stereo + back stereo
 {4, TYPE_SCE, TYPE_CPE, TYPE_CPE, TYPE_LFE}, // 6 channels - front center + stereo + back stereo + LFE
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
    int16_t* samples;                            ///< saved preprocessed input

    int samplerate_index;                        ///< MPEG-4 samplerate index

    ChannelElement *cpe;                         ///< channel elements
    FFPsyContext psy;
    struct FFPsyPreprocessContext* psypp;
    int cur_channel;
    int last_frame;
    float lambda;
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
    const uint8_t *sizes[2];
    int lengths[2];

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

    dsputil_init(&s->dsp, avctx);
    ff_mdct_init(&s->mdct1024, 11, 0);
    ff_mdct_init(&s->mdct128,   8, 0);
    // window init
    ff_kbd_window_init(ff_aac_kbd_long_1024, 4.0, 1024);
    ff_kbd_window_init(ff_aac_kbd_short_128, 6.0, 128);
    ff_sine_window_init(ff_sine_1024, 1024);
    ff_sine_window_init(ff_sine_128, 128);

    s->samples = av_malloc(2 * 1024 * avctx->channels * sizeof(s->samples[0]));
    s->cpe = av_mallocz(sizeof(ChannelElement) * aac_chan_configs[avctx->channels-1][0]);
    avctx->extradata = av_malloc(2);
    avctx->extradata_size = 2;
    put_audio_specific_config(avctx);

    sizes[0] = swb_size_1024[i];
    sizes[1] = swb_size_128[i];
    lengths[0] = ff_aac_num_swb_1024[i];
    lengths[1] = ff_aac_num_swb_128[i];
    ff_psy_init(&s->psy, avctx, 2, sizes, lengths);
    s->psypp = ff_psy_preprocess_init(avctx);

#ifndef CONFIG_HARDCODED_TABLES
    for (i = 0; i < 316; i++)
        ff_aac_pow2sf_tab[i] = pow(2, (i - 200)/4.);
#endif /* CONFIG_HARDCODED_TABLES */

    return 0;
}

static void apply_window_and_mdct(AVCodecContext *avctx, AACEncContext *s,
                                  SingleChannelElement *sce, short *audio, int channel)
{
    int i, j, k;
    const float * lwindow = sce->ics.use_kb_window[0] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow = sce->ics.use_kb_window[0] ? ff_aac_kbd_short_128 : ff_sine_128;
    const float * pwindow = sce->ics.use_kb_window[1] ? ff_aac_kbd_short_128 : ff_sine_128;

    if (sce->ics.window_sequence[0] != EIGHT_SHORT_SEQUENCE) {
        memcpy(s->output, sce->saved, sizeof(float)*1024);
        if(sce->ics.window_sequence[0] == LONG_STOP_SEQUENCE){
            memset(s->output, 0, sizeof(s->output[0]) * 448);
            for(i = 448; i < 576; i++)
                s->output[i] = sce->saved[i] * pwindow[i - 448];
            for(i = 576; i < 704; i++)
                s->output[i] = sce->saved[i];
        }
        if(sce->ics.window_sequence[0] != LONG_START_SEQUENCE){
            j = channel;
            for (i = 0; i < 1024; i++, j += avctx->channels){
                s->output[i+1024]         = audio[j] * lwindow[1024 - i - 1];
                sce->saved[i] = audio[j] * lwindow[i];
            }
        }else{
            j = channel;
            for(i = 0; i < 448; i++, j += avctx->channels)
                s->output[i+1024]         = audio[j];
            for(i = 448; i < 576; i++, j += avctx->channels)
                s->output[i+1024]         = audio[j] * swindow[576 - i - 1];
            memset(s->output+1024+576, 0, sizeof(s->output[0]) * 448);
            j = channel;
            for(i = 0; i < 1024; i++, j += avctx->channels)
                sce->saved[i] = audio[j];
        }
        ff_mdct_calc(&s->mdct1024, sce->coeffs, s->output);
    }else{
        j = channel;
        for (k = 0; k < 1024; k += 128) {
            for(i = 448 + k; i < 448 + k + 256; i++)
                s->output[i - 448 - k] = (i < 1024)
                                         ? sce->saved[i]
                                         : audio[channel + (i-1024)*avctx->channels];
            s->dsp.vector_fmul        (s->output,     k ?  swindow : pwindow, 128);
            s->dsp.vector_fmul_reverse(s->output+128, s->output+128, swindow, 128);
            ff_mdct_calc(&s->mdct128, sce->coeffs + k, s->output);
        }
        j = channel;
        for(i = 0; i < 1024; i++, j += avctx->channels)
            sce->saved[i] = audio[j];
    }
}

/**
 * Encode ics_info element.
 * @see Table 4.6 (syntax of ics_info)
 */
static void put_ics_info(AACEncContext *s, IndividualChannelStream *info)
{
    int w;

    put_bits(&s->pb, 1, 0);                // ics_reserved bit
    put_bits(&s->pb, 2, info->window_sequence[0]);
    put_bits(&s->pb, 1, info->use_kb_window[0]);
    if(info->window_sequence[0] != EIGHT_SHORT_SEQUENCE){
        put_bits(&s->pb, 6, info->max_sfb);
        put_bits(&s->pb, 1, 0);            // no prediction
    }else{
        put_bits(&s->pb, 4, info->max_sfb);
        for(w = 1; w < 8; w++){
            put_bits(&s->pb, 1, !info->group_len[w]);
        }
    }
}

/**
 * Encode MS data.
 * @see 4.6.8.1 "Joint Coding - M/S Stereo"
 */
static void encode_ms_info(PutBitContext *pb, ChannelElement *cpe)
{
    int i, w;

    put_bits(pb, 2, cpe->ms_mode);
    if(cpe->ms_mode == 1){
        for(w = 0; w < cpe->ch[0].ics.num_windows; w += cpe->ch[0].ics.group_len[w]){
            for(i = 0; i < cpe->ch[0].ics.max_sfb; i++)
                put_bits(pb, 1, cpe->ms_mask[w*16 + i]);
        }
    }
}

/**
 * Quantize one coefficient.
 * @return absolute value of the quantized coefficient
 * @see 3GPP TS26.403 5.6.2 "Scalefactor determination"
 */
static av_always_inline int quant(float coef, const float Q)
{
    return av_clip((int)(pow(fabsf(coef) * Q, 0.75) + 0.4054), 0, 8191);
}

/**
 * Calculate rate distortion cost for quantizing with given codebook
 *
 * @return quantization distortion
 */
static float quantize_band_cost(const float *in, int size, int scale_idx, int cb,
                                 const float lambda, const float uplim)
{
    const float Q = ff_aac_pow2sf_tab[200 + scale_idx - SCALE_ONE_POS + SCALE_DIV_512];
    int i, j, k;
    float cost = 0;
    const int dim = cb < FIRST_PAIR_BT ? 4 : 2;

    if(!cb){
        for(i = 0; i < size; i++)
            cost += in[i]*in[i]*lambda;
        return cost;
    }
    for(i = 0; i < size; i += dim){
        float mincost = INFINITY;
        int minidx = 0;
        int minbits = 0;
        const float *vec = ff_aac_codebook_vectors[cb-1];
        for(j = 0; j < ff_aac_spectral_sizes[cb-1]; j++, vec += dim){
            float rd = 0.0f;
            int curbits = ff_aac_spectral_bits[cb-1][minidx];
            if(IS_CODEBOOK_UNSIGNED(cb)){
                for(k = 0; k < dim; k++){
                    float t = fabsf(in[i+k]);
                    float di;
                    //do not code with escape sequence small values
                    if(vec[k] == 64.0f && t < 39.0f*Q){
                        rd = INFINITY;
                        break;
                    }
                    if(vec[k] == 64.0f){//FIXME: slow
                        if(t >= 165140.0f*Q){ // clipped value
                            di = t - 165140.0f;
                            curbits += 21;
                        }else{
                            int c = quant(t, 1.0/Q);
                            di = t - c*cbrt(c)*Q;
                            curbits += av_log2(c)*2 - 4 + 1;
                        }
                    }else{
                        di = t - vec[k]*Q;
                    }
                    if(vec[k] != 0.0f)
                        curbits++;
                    rd += di*di*lambda;
                }
            }else{
                for(k = 0; k < dim; k++){
                    float di = in[i+k] - vec[k]*Q;
                    rd += di*di*lambda;
                }
            }
            rd += curbits;
            if(rd < mincost){
                mincost = rd;
                minidx = j;
                minbits = curbits;
            }
        }
        cost += mincost;
        if(cost >= uplim)
            return uplim;
    }

    return cost;
}

/**
 * Prepare coefficients for encoding.
 *
 * @return sum of coefficient absolute values
 */
static void quantize_and_encode_band(PutBitContext *pb, const float *in, int size,
                                      int scale_idx, int cb, const float lambda)
{
    const float Q  = ff_aac_pow2sf_tab[200 + scale_idx - SCALE_ONE_POS + SCALE_DIV_512];
    const float IQ = ff_aac_pow2sf_tab[200 - scale_idx + SCALE_ONE_POS - SCALE_DIV_512];
    int i, j, k;
    const int dim = cb < FIRST_PAIR_BT ? 4 : 2;
    if(!cb)
        return;

    for(i = 0; i < size; i += dim){
        float mincost = INFINITY;
        int minidx = 0;
        int minbits = 0;
        const float *vec = ff_aac_codebook_vectors[cb-1];
        for(j = 0; j < ff_aac_spectral_sizes[cb-1]; j++, vec += dim){
            float rd = 0.0f;
            int curbits = ff_aac_spectral_bits[cb-1][minidx];
            if(IS_CODEBOOK_UNSIGNED(cb)){
                for(k = 0; k < dim; k++){
                    float t = fabsf(in[i+k]);
                    float di;
                    if(vec[k] != 0.0f)
                        curbits++;
                    //do not code with escape sequence small values
                    if(vec[k] == 64.0f && t < 39.0f*Q){
                        rd = INFINITY;
                        break;
                    }
                    if(vec[k] == 64.0f){//FIXME: slow
                        if(t*IQ >= 165140.0f){ // clipped value
                            di = t - 165140.0f;
                            curbits += 21;
                        }else{
                            int c = quant(t, IQ);
                            di = t - c*cbrt(c)*Q;
                            curbits += av_log2(c)*2 - 4 + 1;
                        }
                    }else{
                        di = t - vec[k]*Q;
                    }
                    rd += di*di*lambda;
                }
            }else{
                for(k = 0; k < dim; k++){
                    float di = in[i+k] - vec[k]*Q;
                    rd += di*di*lambda;
                }
            }
            rd += curbits;
            if(rd < mincost){
                mincost = rd;
                minidx = j;
                minbits = curbits;
            }
        }
        put_bits(pb, ff_aac_spectral_bits[cb-1][minidx], ff_aac_spectral_codes[cb-1][minidx]);
        if(IS_CODEBOOK_UNSIGNED(cb))
            for(j = 0; j < dim; j++)
                if(ff_aac_codebook_vectors[cb-1][minidx*dim+j] != 0.0f)
                    put_bits(pb, 1, in[i+j] < 0.0f);
        if(cb == ESC_BT){
            for(j = 0; j < 2; j++){
                if(ff_aac_codebook_vectors[cb-1][minidx*2+j] == 64.0f){
                    int coef = quant(in[i+j], IQ);
                    int len = av_log2(coef);

                    put_bits(pb, len - 4 + 1, (1 << (len - 4 + 1)) - 2);
                    put_bits(pb, len, coef & ((1 << len) - 1));
                }
            }
        }
    }
}

/**
 * structure used in optimal codebook search
 */
typedef struct BandCodingPath {
    int prev_idx; ///< pointer to the previous path point
    int codebook; ///< codebook for coding band run
    float cost;   ///< path cost
    int run;
} BandCodingPath;

/**
 * Encode band info for single window group bands.
 */
static void encode_window_bands_info(AACEncContext *s, SingleChannelElement *sce,
                                     int win, int group_len)
{
    BandCodingPath path[120][12];
    int w, swb, cb, start, start2, size;
    int i, j;
    const int max_sfb = sce->ics.max_sfb;
    const int run_bits = sce->ics.num_windows == 1 ? 5 : 3;
    const int run_esc = (1 << run_bits) - 1;
    int idx, ppos, count;
    int stackrun[120], stackcb[120], stack_len;

    start = win*128;
    for(cb = 0; cb < 12; cb++){
        path[0][cb].cost = 0.0f;
        path[0][cb].prev_idx = -1;
        path[0][cb].run = 0;
    }
    for(swb = 0; swb < max_sfb; swb++){
        start2 = start;
        size = sce->ics.swb_sizes[swb];
        if(sce->zeroes[win*16 + swb]){
            for(cb = 0; cb < 12; cb++){
                path[swb+1][cb].prev_idx = cb;
                path[swb+1][cb].cost = path[swb][cb].cost;
                path[swb+1][cb].run = path[swb][cb].run + 1;
            }
        }else{
            float minrd = INFINITY;
            int mincb = 0;
            for(cb = 0; cb < 12; cb++){
                float rd = 0.0f;
                for(w = 0; w < group_len; w++){
                    FFPsyBand *band = &s->psy.psy_bands[s->cur_channel*PSY_MAX_BANDS+(win+w)*16+swb];
                    rd += quantize_band_cost(sce->coeffs + start + w*128, size,
                                             sce->sf_idx[(win+w)*16+swb], cb,
                                             s->lambda / band->threshold, INFINITY);
                }
                if(   run_value_bits[sce->ics.num_windows == 8][path[swb][cb].run]
                   != run_value_bits[sce->ics.num_windows == 8][path[swb][cb].run+1])
                    rd += run_bits;
                path[swb+1][cb].prev_idx = cb;
                path[swb+1][cb].cost = path[swb][cb].cost + rd;
                path[swb+1][cb].run = path[swb][cb].run + 1;
                if(rd < minrd){
                    minrd = rd;
                    mincb = cb;
                }
            }
            for(cb = 0; cb < 12; cb++){
                float cost = path[swb][cb].cost + minrd + run_bits + 4;
                if(cost < path[swb+1][cb].cost){
                    path[swb+1][cb].prev_idx = mincb;
                    path[swb+1][cb].cost = cost;
                    path[swb+1][cb].run = 1;
                }
            }
        }
        start += sce->ics.swb_sizes[swb];
    }

    //convert resulting path from backward-linked list
    stack_len = 0;
    idx = 0;
    for(cb = 1; cb < 12; cb++){
        if(path[max_sfb][cb].cost < path[max_sfb][idx].cost)
            idx = cb;
    }
    ppos = max_sfb;
    while(ppos > 0){
        cb = idx;
        stackrun[stack_len] = path[ppos][cb].run;
        stackcb [stack_len] = cb;
        idx = path[ppos][cb].prev_idx;
        ppos -= path[ppos][cb].run;
        stack_len++;
    }
    //perform actual band info encoding
    start = 0;
    for(i = stack_len - 1; i >= 0; i--){
        put_bits(&s->pb, 4, stackcb[i]);
        count = stackrun[i];
        memset(sce->zeroes + win*16 + start, !stackcb[i], count);
        //XXX: memset when band_type is also uint8_t
        for(j = 0; j < count; j++){
            sce->band_type[win*16 + start] =  stackcb[i];
            start++;
        }
        while(count >= run_esc){
            put_bits(&s->pb, run_bits, run_esc);
            count -= run_esc;
        }
        put_bits(&s->pb, run_bits, count);
    }
}

/**
 * Produce integer coefficients from scalefactors provided by the model.
 */
static void adjust_frame_information(AACEncContext *apc, ChannelElement *cpe, int chans)
{
    int i, w, w2, g, ch;
    int start, sum, maxsfb, cmaxsfb;

    for(ch = 0; ch < chans; ch++){
        IndividualChannelStream *ics = &cpe->ch[ch].ics;
        start = 0;
        maxsfb = 0;
        cpe->ch[ch].pulse.num_pulse = 0;
        for(w = 0; w < ics->num_windows*16; w += 16){
            for(g = 0; g < ics->num_swb; g++){
                sum = 0;
                //apply M/S
                if(!ch && cpe->ms_mask[w + g]){
                    for(i = 0; i < ics->swb_sizes[g]; i++){
                        cpe->ch[0].coeffs[start+i] = (cpe->ch[0].coeffs[start+i] + cpe->ch[1].coeffs[start+i]) / 2.0;
                        cpe->ch[1].coeffs[start+i] =  cpe->ch[0].coeffs[start+i] - cpe->ch[1].coeffs[start+i];
                    }
                }
                start += ics->swb_sizes[g];
            }
            for(cmaxsfb = ics->num_swb; cmaxsfb > 0 && cpe->ch[ch].zeroes[w+cmaxsfb-1]; cmaxsfb--);
            maxsfb = FFMAX(maxsfb, cmaxsfb);
        }
        ics->max_sfb = maxsfb;

        //adjust zero bands for window groups
        for(w = 0; w < ics->num_windows; w += ics->group_len[w]){
            for(g = 0; g < ics->max_sfb; g++){
                i = 1;
                for(w2 = w; w2 < w + ics->group_len[w]; w2++){
                    if(!cpe->ch[ch].zeroes[w2*16 + g]){
                        i = 0;
                        break;
                    }
                }
                cpe->ch[ch].zeroes[w*16 + g] = i;
            }
        }
    }

    if(chans > 1 && cpe->common_window){
        IndividualChannelStream *ics0 = &cpe->ch[0].ics;
        IndividualChannelStream *ics1 = &cpe->ch[1].ics;
        int msc = 0;
        ics0->max_sfb = FFMAX(ics0->max_sfb, ics1->max_sfb);
        ics1->max_sfb = ics0->max_sfb;
        for(w = 0; w < ics0->num_windows*16; w += 16)
            for(i = 0; i < ics0->max_sfb; i++)
                if(cpe->ms_mask[w+i]) msc++;
        if(msc == 0 || ics0->max_sfb == 0) cpe->ms_mode = 0;
        else cpe->ms_mode = msc < ics0->max_sfb ? 1 : 2;
    }
}

typedef struct TrellisPath {
    float cost;
    int prev;
    int min_val;
    int max_val;
} TrellisPath;

static void search_for_quantizers_anmr(AACEncContext *s, SingleChannelElement *sce, const float lambda)
{
    int q, w, w2, g, start = 0;
    int i;
    int idx;
    TrellisPath paths[256*121];
    int bandaddr[121];
    int minq;
    float mincost;

    for(i = 0; i < 256; i++){
        paths[i].cost = 0.0f;
        paths[i].prev = -1;
        paths[i].min_val = i;
        paths[i].max_val = i;
    }
    for(i = 256; i < 256*121; i++){
        paths[i].cost = INFINITY;
        paths[i].prev = -2;
        paths[i].min_val = INT_MAX;
        paths[i].max_val = 0;
    }
    idx = 256;
    for(w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]){
        start = w*128;
        for(g = 0; g < sce->ics.num_swb; g++){
            const float *coefs = sce->coeffs + start;
            float qmin, qmax;
            int nz = 0;

            bandaddr[idx >> 8] = w*16+g;
            qmin = INT_MAX;
            qmax = 0.0f;
            for(w2 = 0; w2 < sce->ics.group_len[w]; w2++){
                FFPsyBand *band = &s->psy.psy_bands[s->cur_channel*PSY_MAX_BANDS+(w+w2)*16+g];
                if(band->energy <= band->threshold || band->threshold == 0.0f){
                    sce->zeroes[(w+w2)*16+g] = 1;
                    continue;
                }
                sce->zeroes[(w+w2)*16+g] = 0;
                nz = 1;
                for(i = 0; i < sce->ics.swb_sizes[g]; i++){
                    float t = fabsf(coefs[w2*128+i]);
                    if(t > 0.0f) qmin = fminf(qmin, t);
                    qmax = fmaxf(qmax, t);
                }
            }
            if(nz){
                int minscale, maxscale;
                float minrd = INFINITY;
                //minimum scalefactor index is when minimum nonzero coefficient after quantizing is not clipped
                minscale = av_clip_uint8(log2(qmin)*4 - 69 + SCALE_ONE_POS - SCALE_DIV_512);
                //maximum scalefactor index is when maximum coefficient after quantizing is still not zero
                maxscale = av_clip_uint8(log2(qmax)*4 +  6 + SCALE_ONE_POS - SCALE_DIV_512);
                for(q = minscale; q < maxscale; q++){
                    float dists[12], dist;
                    memset(dists, 0, sizeof(dists));
                    for(w2 = 0; w2 < sce->ics.group_len[w]; w2++){
                        FFPsyBand *band = &s->psy.psy_bands[s->cur_channel*PSY_MAX_BANDS+(w+w2)*16+g];
                        int cb;
                        for(cb = 0; cb <= ESC_BT; cb++){
                            dists[cb] += quantize_band_cost(coefs + w2*128, sce->ics.swb_sizes[g],
                                                            q, cb, s->lambda / band->threshold, INFINITY);
                        }
                    }
                    dist = dists[0];
                    for(i = 1; i <= ESC_BT; i++)
                        dist = fminf(dist, dists[i]);
                    minrd = fminf(minrd, dist);

                    for(i = FFMAX(q - SCALE_MAX_DIFF, 0); i < FFMIN(q + SCALE_MAX_DIFF, 256); i++){
                        float cost;
                        int minv, maxv;
                        if(isinf(paths[idx - 256 + i].cost))
                            continue;
                        cost = paths[idx - 256 + i].cost + dist
                               + ff_aac_scalefactor_bits[q - i + SCALE_DIFF_ZERO];
                        minv = FFMIN(paths[idx - 256 + i].min_val, q);
                        maxv = FFMAX(paths[idx - 256 + i].max_val, q);
                        if(cost < paths[idx + q].cost && maxv-minv < SCALE_MAX_DIFF){
                            paths[idx + q].cost = cost;
                            paths[idx + q].prev = idx - 256 + i;
                            paths[idx + q].min_val = minv;
                            paths[idx + q].max_val = maxv;
                        }
                    }
                }
            }else{
                for(q = 0; q < 256; q++){
                    if(!isinf(paths[idx - 256 + q].cost)){
                        paths[idx + q].cost = paths[idx - 256 + q].cost + 1;
                        paths[idx + q].prev = idx - 256 + q;
                        paths[idx + q].min_val = FFMIN(paths[idx - 256 + q].min_val, q);
                        paths[idx + q].max_val = FFMAX(paths[idx - 256 + q].max_val, q);
                        continue;
                    }
                    for(i = FFMAX(q - SCALE_MAX_DIFF, 0); i < FFMIN(q + SCALE_MAX_DIFF, 256); i++){
                        float cost;
                        int minv, maxv;
                        if(isinf(paths[idx - 256 + i].cost))
                            continue;
                        cost = paths[idx - 256 + i].cost + ff_aac_scalefactor_bits[q - i + SCALE_DIFF_ZERO];
                        minv = FFMIN(paths[idx - 256 + i].min_val, q);
                        maxv = FFMAX(paths[idx - 256 + i].max_val, q);
                        if(cost < paths[idx + q].cost && maxv-minv < SCALE_MAX_DIFF){
                            paths[idx + q].cost = cost;
                            paths[idx + q].prev = idx - 256 + i;
                            paths[idx + q].min_val = minv;
                            paths[idx + q].max_val = maxv;
                        }
                    }
                }
            }
            sce->zeroes[w*16+g] = !nz;
            start += sce->ics.swb_sizes[g];
            idx += 256;
        }
    }
    idx -= 256;
    mincost = paths[idx].cost;
    minq = idx;
    for(i = 1; i < 256; i++){
        if(paths[idx + i].cost < mincost){
            mincost = paths[idx + i].cost;
            minq = idx + i;
        }
    }
    while(minq >= 256){
        sce->sf_idx[bandaddr[minq>>8]] = minq & 0xFF;
        minq = paths[minq].prev;
    }
    //set the same quantizers inside window groups
    for(w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w])
        for(g = 0;  g < sce->ics.num_swb; g++)
            for(w2 = 1; w2 < sce->ics.group_len[w]; w2++)
                sce->sf_idx[(w+w2)*16+g] = sce->sf_idx[w*16+g];
}

static void search_for_ms(AACEncContext *s, ChannelElement *cpe, const float lambda)
{
    int start = 0, i, w, w2, g;
    float M[128], S[128];
    SingleChannelElement *sce0 = &cpe->ch[0];
    SingleChannelElement *sce1 = &cpe->ch[1];
    if(!cpe->common_window)
        return;
    for(w = 0; w < sce0->ics.num_windows; w += sce0->ics.group_len[w]){
        for(g = 0;  g < sce0->ics.num_swb; g++){
            if(!cpe->ch[0].zeroes[w*16+g] && !cpe->ch[1].zeroes[w*16+g]){
                float dist1 = 0.0f, dist2 = 0.0f;
                for(w2 = 0; w2 < sce0->ics.group_len[w]; w2++){
                    FFPsyBand *band0 = &s->psy.psy_bands[(s->cur_channel+0)*PSY_MAX_BANDS+(w+w2)*16+g];
                    FFPsyBand *band1 = &s->psy.psy_bands[(s->cur_channel+1)*PSY_MAX_BANDS+(w+w2)*16+g];
                    float minthr = fminf(band0->threshold, band1->threshold);
                    float maxthr = fmaxf(band0->threshold, band1->threshold);
                    for(i = 0; i < sce0->ics.swb_sizes[g]; i++){
                        M[i] = (sce0->coeffs[start+w2*128+i]
                              + sce1->coeffs[start+w2*128+i])*0.5;
                        S[i] =  sce0->coeffs[start+w2*128+i]
                              - sce1->coeffs[start+w2*128+i];
                    }
                    dist1 += quantize_band_cost(sce0->coeffs + start + w2*128,
                                                sce0->ics.swb_sizes[g],
                                                sce0->sf_idx[(w+w2)*16+g],
                                                sce0->band_type[(w+w2)*16+g],
                                                lambda / band0->threshold, INFINITY);
                    dist1 += quantize_band_cost(sce1->coeffs + start + w2*128,
                                                sce1->ics.swb_sizes[g],
                                                sce1->sf_idx[(w+w2)*16+g],
                                                sce1->band_type[(w+w2)*16+g],
                                                lambda / band1->threshold, INFINITY);
                    dist2 += quantize_band_cost(M,
                                                sce0->ics.swb_sizes[g],
                                                sce0->sf_idx[(w+w2)*16+g],
                                                sce0->band_type[(w+w2)*16+g],
                                                lambda / maxthr, INFINITY);
                    dist2 += quantize_band_cost(S,
                                                sce1->ics.swb_sizes[g],
                                                sce1->sf_idx[(w+w2)*16+g],
                                                sce1->band_type[(w+w2)*16+g],
                                                lambda / minthr, INFINITY);
                }
                cpe->ms_mask[w*16+g] = dist2 < dist1;
            }
            start += sce0->ics.swb_sizes[g];
        }
    }
}

/**
 * Encode scalefactor band coding type.
 */
static void encode_band_info(AACEncContext *s, SingleChannelElement *sce)
{
    int w;

    for(w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]){
        encode_window_bands_info(s, sce, w, sce->ics.group_len[w]);
    }
}

/**
 * Encode scalefactors.
 */
static void encode_scale_factors(AVCodecContext *avctx, AACEncContext *s, SingleChannelElement *sce)
{
    int off = sce->sf_idx[0], diff;
    int i, w;

    for(w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]){
        for(i = 0; i < sce->ics.max_sfb; i++){
            if(!sce->zeroes[w*16 + i]){
                diff = sce->sf_idx[w*16 + i] - off + SCALE_DIFF_ZERO;
                if(diff < 0 || diff > 120) av_log(avctx, AV_LOG_ERROR, "Scalefactor difference is too big to be coded\n");
                off = sce->sf_idx[w*16 + i];
                put_bits(&s->pb, ff_aac_scalefactor_bits[diff], ff_aac_scalefactor_code[diff]);
            }
        }
    }
}

/**
 * Encode pulse data.
 */
static void encode_pulses(AACEncContext *s, Pulse *pulse)
{
    int i;

    put_bits(&s->pb, 1, !!pulse->num_pulse);
    if(!pulse->num_pulse) return;

    put_bits(&s->pb, 2, pulse->num_pulse - 1);
    put_bits(&s->pb, 6, pulse->start);
    for(i = 0; i < pulse->num_pulse; i++){
        put_bits(&s->pb, 5, pulse->pos[i]);
        put_bits(&s->pb, 4, pulse->amp[i]);
    }
}

/**
 * Encode spectral coefficients processed by psychoacoustic model.
 */
static void encode_spectral_coeffs(AACEncContext *s, SingleChannelElement *sce)
{
    int start, i, w, w2;

    for(w = 0; w < sce->ics.num_windows; w += sce->ics.group_len[w]){
        start = 0;
        for(i = 0; i < sce->ics.max_sfb; i++){
            if(sce->zeroes[w*16 + i]){
                start += sce->ics.swb_sizes[i];
                continue;
            }
            for(w2 = w; w2 < w + sce->ics.group_len[w]; w2++){
                quantize_and_encode_band(&s->pb, sce->coeffs + start + w2*128,
                                         sce->ics.swb_sizes[i],
                                         sce->sf_idx[w*16 + i],
                                         sce->band_type[w*16 + i],
                                         s->lambda);
            }
            start += sce->ics.swb_sizes[i];
        }
    }
}

/**
 * Encode one channel of audio data.
 */
static int encode_individual_channel(AVCodecContext *avctx, AACEncContext *s, SingleChannelElement *sce, int common_window)
{
    put_bits(&s->pb, 8, sce->sf_idx[0]);
    if(!common_window) put_ics_info(s, &sce->ics);
    encode_band_info(s, sce);
    encode_scale_factors(avctx, s, sce);
    encode_pulses(s, &sce->pulse);
    put_bits(&s->pb, 1, 0); //tns
    put_bits(&s->pb, 1, 0); //ssr
    encode_spectral_coeffs(s, sce);
    return 0;
}

/**
 * Write some auxiliary information about the created AAC file.
 */
static void put_bitstream_info(AVCodecContext *avctx, AACEncContext *s, const char *name)
{
    int i, namelen, padbits;

    namelen = strlen(name) + 2;
    put_bits(&s->pb, 3, TYPE_FIL);
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

    if(s->last_frame)
        return 0;
    if(data){
        if(!s->psypp){
            memcpy(s->samples + 1024 * avctx->channels, data, 1024 * avctx->channels * sizeof(s->samples[0]));
        }else{
            start_ch = 0;
            samples2 = s->samples + 1024 * avctx->channels;
            for(i = 0; i < chan_map[0]; i++){
                tag = chan_map[i+1];
                chans = tag == TYPE_CPE ? 2 : 1;
                ff_psy_preprocess(s->psypp, (uint16_t*)data + start_ch, samples2 + start_ch, start_ch + i, chans);
                start_ch += chans;
            }
        }
    }
    if(!avctx->frame_number){
        memcpy(s->samples, s->samples + 1024 * avctx->channels, 1024 * avctx->channels * sizeof(s->samples[0]));
        return 0;
    }

    init_put_bits(&s->pb, frame, buf_size*8);
    if((avctx->frame_number & 0xFF)==1 && !(avctx->flags & CODEC_FLAG_BITEXACT)){
        put_bitstream_info(avctx, s, LIBAVCODEC_IDENT);
    }
    start_ch = 0;
    memset(chan_el_counter, 0, sizeof(chan_el_counter));
    for(i = 0; i < chan_map[0]; i++){
        FFPsyWindowInfo wi[2];
        tag = chan_map[i+1];
        chans = tag == TYPE_CPE ? 2 : 1;
        cpe = &s->cpe[i];
        samples2 = samples + start_ch;
        la = samples2 + 1024 * avctx->channels + start_ch;
        if(!data) la = NULL;
        s->lambda = 5e-7f;
        for(j = 0; j < chans; j++){
            IndividualChannelStream *ics = &cpe->ch[j].ics;
            int k;
            wi[j] = ff_psy_suggest_window(&s->psy, samples2, la, start_ch + j, ics->window_sequence[0]);
            ics->window_sequence[1] = ics->window_sequence[0];
            ics->window_sequence[0] = wi[j].window_type[0];
            ics->use_kb_window[1]   = ics->use_kb_window[0];
            ics->use_kb_window[0]   = wi[j].window_shape;
            ics->num_windows        = wi[j].num_windows;
            ics->swb_sizes          = s->psy.bands    [ics->num_windows == 8];
            ics->num_swb            = s->psy.num_bands[ics->num_windows == 8];
            for(k = 0; k < ics->num_windows; k++)
                ics->group_len[k] = wi[j].grouping[k];

            s->cur_channel = start_ch + j;
            apply_window_and_mdct(avctx, s, &cpe->ch[j], samples2, j);
            search_for_quantizers_anmr(s, &cpe->ch[j], s->lambda);
        }
        cpe->common_window = 0;
        if(chans > 1
            && wi[0].window_type[0] == wi[1].window_type[0]
            && wi[0].window_shape   == wi[1].window_shape){

            cpe->common_window = 1;
            for(j = 0; j < wi[0].num_windows; j++){
                if(wi[0].grouping[j] != wi[1].grouping[j]){
                    cpe->common_window = 0;
                    break;
                }
            }
        }
//        search_for_ms(s, cpe, s->lambda);
        adjust_frame_information(s, cpe, chans);
        put_bits(&s->pb, 3, tag);
        put_bits(&s->pb, 4, chan_el_counter[tag]++);
        if(chans == 2){
            put_bits(&s->pb, 1, cpe->common_window);
            if(cpe->common_window){
                put_ics_info(s, &cpe->ch[0].ics);
                encode_ms_info(&s->pb, cpe);
            }
        }
        for(j = 0; j < chans; j++){
            s->cur_channel = start_ch + j;
            ff_psy_set_band_info(&s->psy, s->cur_channel, cpe->ch[j].coeffs, &wi[j]);
            encode_individual_channel(avctx, s, &cpe->ch[j], cpe->common_window);
        }
        start_ch += chans;
    }

    put_bits(&s->pb, 3, TYPE_END);
    flush_put_bits(&s->pb);
    avctx->frame_bits = put_bits_count(&s->pb);

    if(!data)
        s->last_frame = 1;
    memcpy(s->samples, s->samples + 1024 * avctx->channels, 1024 * avctx->channels * sizeof(s->samples[0]));
    return put_bits_count(&s->pb)>>3;
}

static av_cold int aac_encode_end(AVCodecContext *avctx)
{
    AACEncContext *s = avctx->priv_data;

    ff_mdct_end(&s->mdct1024);
    ff_mdct_end(&s->mdct128);
    ff_psy_end(&s->psy);
    ff_psy_preprocess_end(s->psypp);
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
    .capabilities = CODEC_CAP_SMALL_LAST_FRAME | CODEC_CAP_DELAY,
    .sample_fmts = (enum SampleFormat[]){SAMPLE_FMT_S16,SAMPLE_FMT_NONE},
    .long_name = NULL_IF_CONFIG_SMALL("Advanced Audio Coding"),
};
