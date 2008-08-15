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
 * psy model selection with some option
 * add sane pulse detection
 ***********************************/

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
 * structure used in optimal codebook search
 */
typedef struct BandCodingPath {
    int prev_idx; ///< pointer to the previous path point
    int codebook; ///< codebook for coding band run
    int bits;     ///< number of bit needed to code given number of bands
} BandCodingPath;

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
    const uint8_t *swb_sizes1024;                ///< scalefactor band sizes for long frame
    int swb_num1024;                             ///< number of scalefactor bands for long frame
    const uint8_t *swb_sizes128;                 ///< scalefactor band sizes for short frame
    int swb_num128;                              ///< number of scalefactor bands for short frame

    ChannelElement *cpe;                         ///< channel elements
    AACPsyContext psy;                           ///< psychoacoustic model context
    int last_frame;
    BandCodingPath path[64];                     ///< auxiliary data needed for optimal band info coding
    int band_bits[64][12];                       ///< bits needed to encode each band with each codebook
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
    s->swb_num1024   = ff_aac_num_swb_1024[i];
    s->swb_sizes128  = swb_size_128[i];
    s->swb_num128    = ff_aac_num_swb_128[i];

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
    if(ff_aac_psy_init(&s->psy, avctx, AAC_PSY_3GPP, aac_chan_configs[avctx->channels-1][0], 0, s->swb_sizes1024, s->swb_num1024, s->swb_sizes128, s->swb_num128) < 0){
        av_log(avctx, AV_LOG_ERROR, "Cannot initialize selected model.\n");
        return -1;
    }
    avctx->extradata = av_malloc(2);
    avctx->extradata_size = 2;
    put_audio_specific_config(avctx);
    return 0;
}

static void apply_window_and_mdct(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, short *audio, int channel)
{
    int i, j, k;
    const float * lwindow = cpe->ch[channel].ics.use_kb_window[0] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow = cpe->ch[channel].ics.use_kb_window[0] ? ff_aac_kbd_short_128 : ff_sine_128;
    const float * pwindow = cpe->ch[channel].ics.use_kb_window[1] ? ff_aac_kbd_short_128 : ff_sine_128;

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
                s->output[i+1024]         = audio[j] * lwindow[1024 - i - 1];
                cpe->ch[channel].saved[i] = audio[j] * lwindow[i];
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
                cpe->ch[channel].saved[i] = audio[j];
        }
        ff_mdct_calc(&s->mdct1024, cpe->ch[channel].coeffs, s->output);
    }else{
        j = channel;
        for (k = 0; k < 1024; k += 128) {
            for(i = 448 + k; i < 448 + k + 256; i++)
                s->output[i - 448 - k] = (i < 1024) ? cpe->ch[channel].saved[i] : audio[channel + (i-1024)*avctx->channels] / 512.0;
            s->dsp.vector_fmul        (s->output,     k ?  swindow : pwindow, 128);
            s->dsp.vector_fmul_reverse(s->output+128, s->output+128, swindow, 128);
            ff_mdct_calc(&s->mdct128, cpe->ch[channel].coeffs + k, s->output);
        }
        j = channel;
        for(i = 0; i < 1024; i++, j += avctx->channels)
            cpe->ch[channel].saved[i] = audio[j];
    }
}

/**
 * Encode ics_info element.
 * @see Table 4.6 (syntax of ics_info)
 */
static void put_ics_info(AVCodecContext *avctx, IndividualChannelStream *info)
{
    AACEncContext *s = avctx->priv_data;
    int wg;

    put_bits(&s->pb, 1, 0);                // ics_reserved bit
    put_bits(&s->pb, 2, info->window_sequence[0]);
    put_bits(&s->pb, 1, info->use_kb_window[0]);
    if(info->window_sequence[0] != EIGHT_SHORT_SEQUENCE){
        put_bits(&s->pb, 6, info->max_sfb);
        put_bits(&s->pb, 1, 0);            // no prediction
    }else{
        put_bits(&s->pb, 4, info->max_sfb);
        for(wg = 0; wg < info->num_window_groups; wg++){
            if(wg)
                put_bits(&s->pb, 1, 0);
            if(info->group_len[wg] > 1)
                put_sbits(&s->pb, info->group_len[wg] - 1, 0xFF);
        }
    }
}

/**
 * Encode MS data.
 * @see 4.6.8.1 "Joint Coding - M/S Stereo"
 */
static void encode_ms_info(PutBitContext *pb, ChannelElement *cpe)
{
    int i, w, wg;

    put_bits(pb, 2, cpe->ms.present);
    if(cpe->ms.present == 1){
        w = 0;
        for(wg = 0; wg < cpe->ch[0].ics.num_window_groups; wg++){
            for(i = 0; i < cpe->ch[0].ics.max_sfb; i++)
                put_bits(pb, 1, cpe->ms.mask[w][i]);
            w += cpe->ch[0].ics.group_len[wg];
        }
    }
}

/**
 * Return number of bits needed to write codebook run length value.
 *
 * @param run     run length
 * @param bits    number of bits used to code value (5 for long frames, 3 for short frames)
 */
static av_always_inline int calculate_run_bits(int run, const int bits)
{
    int esc = (1 << bits) - 1;
    return (1 + (run >= esc)) * bits;
}

/**
 * Calculate the number of bits needed to code given band with given codebook.
 *
 * @param s       encoder context
 * @param cpe     channel element
 * @param channel channel number inside channel pair
 * @param win     window group start number
 * @param start   scalefactor band position in spectral coefficients
 * @param size    scalefactor band size
 * @param cb      codebook number
 */
static int calculate_band_bits(AACEncContext *s, ChannelElement *cpe, int channel, int win, int group_len, int start, int size, int cb)
{
    int i, j, w;
    int score = 0, dim, idx, start2;
    int range;

    if(!cb) return 0;
    cb--;
    dim = (aac_cb_info[cb].flags & CB_PAIRS) ? 2 : 4;
    if(aac_cb_info[cb].flags & CB_UNSIGNED)
        range = aac_cb_info[cb].maxval + 1;
    else
        range = aac_cb_info[cb].maxval*2 + 1;

    start2 = start;
    if(aac_cb_info[cb].flags & CB_ESCAPE){
        int coef_abs[2];
        for(w = win; w < win + group_len; w++){
            for(i = start2; i < start2 + size; i += dim){
                idx = 0;
                for(j = 0; j < dim; j++)
                    coef_abs[j] = FFABS(cpe->ch[channel].icoefs[i+j]);
                for(j = 0; j < dim; j++)
                    idx = idx*17 + FFMIN(coef_abs[j], 16);
                score += ff_aac_spectral_bits[cb][idx];
                for(j = 0; j < dim; j++)
                    if(cpe->ch[channel].icoefs[i+j])
                        score++;
                for(j = 0; j < dim; j++)
                    if(coef_abs[j] > 15)
                        score += av_log2(coef_abs[j]) * 2 - 4 + 1;
            }
            start2 += 128;
       }
    }else if(aac_cb_info[cb].flags & CB_UNSIGNED){
        for(w = win; w < win + group_len; w++){
            for(i = start2; i < start2 + size; i += dim){
                idx = 0;
                for(j = 0; j < dim; j++)
                    idx = idx * range + FFABS(cpe->ch[channel].icoefs[i+j]);
                score += ff_aac_spectral_bits[cb][idx];
                for(j = 0; j < dim; j++)
                     if(cpe->ch[channel].icoefs[i+j])
                         score++;
            }
            start2 += 128;
        }
    }else{
        for(w = win; w < win + group_len; w++){
            for(i = start2; i < start2 + size; i += dim){
                idx = 0;
                for(j = 0; j < dim; j++)
                    idx = idx * range + cpe->ch[channel].icoefs[i+j] + aac_cb_info[cb].maxval;
                score += ff_aac_spectral_bits[cb][idx];
            }
            start2 += 128;
        }
    }
    return score;
}

/**
 * Encode band info for single window group bands.
 */
static void encode_window_bands_info(AACEncContext *s, ChannelElement *cpe, int channel, int win, int group_len){
    int maxval;
    int w, swb, cb, ccb, start, start2, size;
    int i, j, k;
    const int max_sfb = cpe->ch[channel].ics.max_sfb;
    const int run_bits = cpe->ch[channel].ics.num_windows == 1 ? 5 : 3;
    const int run_esc = (1 << run_bits) - 1;
    int bits, idx, count;
    int stack[64], stack_len;

    start = win*128;
    for(swb = 0; swb < max_sfb; swb++){
        maxval = 0;
        start2 = start;
        size = cpe->ch[channel].ics.swb_sizes[swb];
        if(cpe->ch[channel].zeroes[win][swb])
            maxval = 0;
        else{
            for(w = win; w < win + group_len; w++){
                for(i = start2; i < start2 + size; i++){
                    maxval = FFMAX(maxval, FFABS(cpe->ch[channel].icoefs[i]));
                }
                start2 += 128;
            }
        }
        for(cb = 0; cb < 12; cb++){
            if(aac_cb_info[cb].maxval < maxval)
                s->band_bits[swb][cb] = INT_MAX;
            else
                s->band_bits[swb][cb] = calculate_band_bits(s, cpe, channel, win, group_len, start, size, cb);
        }
        start += cpe->ch[channel].ics.swb_sizes[swb];
    }
    s->path[0].bits = 0;
    for(i = 1; i <= max_sfb; i++)
        s->path[i].bits = INT_MAX;
    for(i = 0; i < max_sfb; i++){
        for(j = 1; j <= max_sfb - i; j++){
            bits = INT_MAX;
            ccb = 0;
            for(cb = 0; cb < 12; cb++){
                int sum = 0;
                for(k = 0; k < j; k++){
                    if(s->band_bits[i + k][cb] == INT_MAX){
                        sum = INT_MAX;
                        break;
                    }
                    sum += s->band_bits[i + k][cb];
                }
                if(sum < bits){
                    bits = sum;
                    ccb  = cb;
                }
            }
            assert(bits != INT_MAX);
            bits += s->path[i].bits + calculate_run_bits(j, run_bits);
            if(bits < s->path[i+j].bits){
                s->path[i+j].bits     = bits;
                s->path[i+j].codebook = ccb;
                s->path[i+j].prev_idx = i;
            }
        }
    }

    //convert resulting path from backward-linked list
    stack_len = 0;
    idx = max_sfb;
    while(idx > 0){
        stack[stack_len++] = idx;
        idx = s->path[idx].prev_idx;
    }

    //perform actual band info encoding
    start = 0;
    for(i = stack_len - 1; i >= 0; i--){
        put_bits(&s->pb, 4, s->path[stack[i]].codebook);
        count = stack[i] - s->path[stack[i]].prev_idx;
        for(j = 0; j < count; j++){
            cpe->ch[channel].band_type[win][start] =  s->path[stack[i]].codebook;
            cpe->ch[channel].zeroes[win][start]    = !s->path[stack[i]].codebook;
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
 * Encode one scalefactor band with selected codebook.
 */
static void encode_band_coeffs(AACEncContext *s, ChannelElement *cpe, int channel, int start, int size, int cb)
{
    const uint8_t  *bits  = ff_aac_spectral_bits [aac_cb_info[cb].cb_num];
    const uint16_t *codes = ff_aac_spectral_codes[aac_cb_info[cb].cb_num];
    const int dim = (aac_cb_info[cb].flags & CB_PAIRS) ? 2 : 4;
    int i, j, idx, range;

    if(!bits) return;

    if(aac_cb_info[cb].flags & CB_UNSIGNED)
        range = aac_cb_info[cb].maxval + 1;
    else
        range = aac_cb_info[cb].maxval*2 + 1;

    if(aac_cb_info[cb].flags & CB_ESCAPE){
        int coef_abs[2];
        for(i = start; i < start + size; i += dim){
            idx = 0;
            for(j = 0; j < dim; j++)
                coef_abs[j] = FFABS(cpe->ch[channel].icoefs[i+j]);
            for(j = 0; j < dim; j++)
                idx = idx*17 + FFMIN(coef_abs[j], 16);
            put_bits(&s->pb, bits[idx], codes[idx]);
            //output signs
            for(j = 0; j < dim; j++)
                if(cpe->ch[channel].icoefs[i+j])
                    put_bits(&s->pb, 1, cpe->ch[channel].icoefs[i+j] < 0);
            //output escape values
            for(j = 0; j < dim; j++)
                if(coef_abs[j] > 15){
                    int len = av_log2(coef_abs[j]);

                    put_bits(&s->pb, len - 4 + 1, (1 << (len - 4 + 1)) - 2);
                    put_bits(&s->pb, len, coef_abs[j] & ((1 << len) - 1));
                }
        }
    }else if(aac_cb_info[cb].flags & CB_UNSIGNED){
        for(i = start; i < start + size; i += dim){
            idx = 0;
            for(j = 0; j < dim; j++)
                idx = idx * range + FFABS(cpe->ch[channel].icoefs[i+j]);
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
                idx = idx * range + cpe->ch[channel].icoefs[i+j] + aac_cb_info[cb].maxval;
            put_bits(&s->pb, bits[idx], codes[idx]);
        }
    }
}

/**
 * Encode scalefactor band coding type.
 */
static void encode_band_info(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel)
{
    int w, wg;

    w = 0;
    for(wg = 0; wg < cpe->ch[channel].ics.num_window_groups; wg++){
        encode_window_bands_info(s, cpe, channel, w, cpe->ch[channel].ics.group_len[wg]);
        w += cpe->ch[channel].ics.group_len[wg];
    }
}

/**
 * Encode scalefactors.
 */
static void encode_scale_factors(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel, int global_gain)
{
    int off = global_gain, diff;
    int i, w, wg;

    w = 0;
    for(wg = 0; wg < cpe->ch[channel].ics.num_window_groups; wg++){
        for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
            if(!cpe->ch[channel].zeroes[w][i]){
                if(cpe->ch[channel].sf_idx[w][i] == 256) cpe->ch[channel].sf_idx[w][i] = off;
                diff = cpe->ch[channel].sf_idx[w][i] - off + SCALE_DIFF_ZERO;
                if(diff < 0 || diff > 120) av_log(avctx, AV_LOG_ERROR, "Scalefactor difference is too big to be coded\n");
                off = cpe->ch[channel].sf_idx[w][i];
                put_bits(&s->pb, ff_aac_scalefactor_bits[diff], ff_aac_scalefactor_code[diff]);
            }
        }
        w += cpe->ch[channel].ics.group_len[wg];
    }
}

/**
 * Encode pulse data.
 */
static void encode_pulses(AVCodecContext *avctx, AACEncContext *s, Pulse *pulse, int channel)
{
    int i;

    put_bits(&s->pb, 1, !!pulse->num_pulse);
    if(!pulse->num_pulse) return;

    put_bits(&s->pb, 2, pulse->num_pulse - 1);
    put_bits(&s->pb, 6, pulse->start);
    for(i = 0; i < pulse->num_pulse; i++){
        put_bits(&s->pb, 5, pulse->offset[i]);
        put_bits(&s->pb, 4, pulse->amp[i]);
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
static void encode_spectral_coeffs(AVCodecContext *avctx, AACEncContext *s, ChannelElement *cpe, int channel)
{
    int start, i, w, w2, wg;

    w = 0;
    for(wg = 0; wg < cpe->ch[channel].ics.num_window_groups; wg++){
        start = 0;
        for(i = 0; i < cpe->ch[channel].ics.max_sfb; i++){
            if(cpe->ch[channel].zeroes[w][i]){
                start += cpe->ch[channel].ics.swb_sizes[i];
                continue;
            }
            for(w2 = w; w2 < w + cpe->ch[channel].ics.group_len[wg]; w2++){
                encode_band_coeffs(s, cpe, channel, start + w2*128, cpe->ch[channel].ics.swb_sizes[i], cpe->ch[channel].band_type[w][i]);
            }
            start += cpe->ch[channel].ics.swb_sizes[i];
        }
        w += cpe->ch[channel].ics.group_len[wg];
    }
}

/**
 * Encode one channel of audio data.
 */
static int encode_individual_channel(AVCodecContext *avctx, ChannelElement *cpe, int channel)
{
    AACEncContext *s = avctx->priv_data;
    int g, w, wg;
    int global_gain;

    //determine global gain as standard recommends - the first scalefactor value
    global_gain = 0;
    w = 0;
    for(wg = 0; wg < cpe->ch[channel].ics.num_window_groups; wg++){
        for(g = 0; g < cpe->ch[channel].ics.max_sfb; g++){
            if(!cpe->ch[channel].zeroes[w][g]){
                global_gain = cpe->ch[channel].sf_idx[w][g];
                break;
            }
        }
        if(global_gain) break;
        w += cpe->ch[channel].ics.group_len[wg];
    }

    put_bits(&s->pb, 8, global_gain);
    if(!cpe->common_window) put_ics_info(avctx, &cpe->ch[channel].ics);
    encode_band_info(avctx, s, cpe, channel);
    encode_scale_factors(avctx, s, cpe, channel, global_gain);
    encode_pulses(avctx, s, &cpe->ch[channel].pulse, channel);
    encode_tns_data(avctx, s, cpe, channel);
    put_bits(&s->pb, 1, 0); //ssr
    encode_spectral_coeffs(avctx, s, cpe, channel);
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

    if(s->last_frame)
        return 0;
    if(data){
        if((s->psy.flags & PSY_MODEL_NO_PREPROC) == PSY_MODEL_NO_PREPROC){
            memcpy(s->samples + 1024 * avctx->channels, data, 1024 * avctx->channels * sizeof(s->samples[0]));
        }else{
            start_ch = 0;
            samples2 = s->samples + 1024 * avctx->channels;
            for(i = 0; i < chan_map[0]; i++){
                tag = chan_map[i+1];
                chans = tag == ID_CPE ? 2 : 1;
                ff_aac_psy_preprocess(&s->psy, (uint16_t*)data + start_ch, samples2 + start_ch, i, tag);
                start_ch += chans;
            }
        }
    }
    if(!avctx->frame_number){
        memcpy(s->samples, s->samples + 1024 * avctx->channels, 1024 * avctx->channels * sizeof(s->samples[0]));
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
            apply_window_and_mdct(avctx, s, cpe, samples2, j);
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
    .capabilities = CODEC_CAP_SMALL_LAST_FRAME | CODEC_CAP_DELAY,
    .sample_fmts = (enum SampleFormat[]){SAMPLE_FMT_S16,SAMPLE_FMT_NONE},
    .long_name = NULL_IF_CONFIG_SMALL("Advanced Audio Coding"),
};
