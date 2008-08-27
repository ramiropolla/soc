/*
 * Nellymoser encoder
 * This code is developed as part of Google Summer of Code 2008 Program.
 *
 * Copyright (c) 2008 Bartlomiej Wolowiec
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
 * @file nellymoserenc.c
 * Nellymoser encoder
 * by Bartlomiej Wolowiec
 *
 * Generic codec information: libavcodec/nellymoserdec.c
 * Log search algorithm idea: http://www1.mplayerhq.hu/ASAO/ASAO.zip
 * (Copyright Joseph Artsimovich and UAB "DKD")
 *
 * for more information about nellymoser format, visit:
 * http://wiki.multimedia.cx/index.php?title=Nellymoser
 */

#include "nellymoser.h"
#include "avcodec.h"
#include "dsputil.h"

#define BITSTREAM_WRITER_LE
#include "bitstream.h"

#define MAX_POW_CACHED (1<<11)
#define LOWPASS

#ifdef LOWPASS
#include "lowpass2.h"
#endif

typedef struct NellyMoserEncodeContext {
    AVCodecContext *avctx;
    int last_frame;
    int bufsize;
    int bits[NELLY_BUF_LEN];
    float pows[NELLY_FILL_LEN];
    DSPContext dsp;
    MDCTContext mdct_ctx;
#ifdef LOWPASS
    LPFilterContext lp;
#endif
     DECLARE_ALIGNED_16(float, mdct_out[NELLY_SAMPLES]);
     DECLARE_ALIGNED_16(float, buf[2 * NELLY_SAMPLES]);
} NellyMoserEncodeContext;

static DECLARE_ALIGNED_16(float, sine_window[NELLY_SAMPLES]);
static float pow_table[MAX_POW_CACHED];

void apply_mdct(NellyMoserEncodeContext *s, float *in, float *coefs)
{
    DECLARE_ALIGNED_16(float, in_buff[NELLY_SAMPLES]);

    memcpy(&in_buff[0], &in[0], NELLY_SAMPLES * sizeof(float));
    s->dsp.vector_fmul(in_buff, sine_window, NELLY_SAMPLES);
    memset(coefs, 0, NELLY_BUF_LEN * sizeof(float));
    ff_mdct_calc(&s->mdct_ctx, coefs, in_buff);
}

static av_cold int encode_init(AVCodecContext *avctx)
{
    NellyMoserEncodeContext *s = avctx->priv_data;
    int i;

    if (avctx->channels != 1) {
        av_log(avctx, AV_LOG_ERROR, "Nellymoser supports only 1 channel\n");
        return -1;
    }

    switch (avctx->sample_rate) {
    case 8000:
#ifdef LOWPASS
        ff_lowpass_init(&s->lp, 8000, 2000, 3800, 1, 60);
#endif
        break;
    case 11025:
#ifdef LOWPASS
        ff_lowpass_init(&s->lp, 11025, 3000, 5000, 1, 60);
#endif
        break;
    case 22050:
#ifdef LOWPASS
        ff_lowpass_init(&s->lp, 22025, 6000, 10000, 1, 60);
#endif
        break;
    case 44100:
#ifdef LOWPASS
        ff_lowpass_init(&s->lp, 44100, 12000, 20000, 1, 60);
#endif
        break;
    default:
        av_log(avctx, AV_LOG_ERROR,
               "Nellymoser works only with 8000, 11025, 22050 and 44100 sample rate\n");
        return -1;
    }

    avctx->frame_size = NELLY_SAMPLES;
    s->avctx = avctx;
    s->bufsize = 0;
    s->last_frame = 0;
    ff_mdct_init(&s->mdct_ctx, 8, 0);
    dsputil_init(&s->dsp, avctx);

    /* Generate overlap window */
    if (!sine_window[0]) {
        ff_sine_window_init(sine_window, 128);
        for (i = 0; i < 128; i++) {
            sine_window[255 - i] = sine_window[i];
        }
    }
    for (i = 0; i < MAX_POW_CACHED; i++)
        pow_table[i] = -pow(2, -i / 2048.0 - 3.0);

    return 0;
}

static av_cold int encode_end(AVCodecContext *avctx)
{
    NellyMoserEncodeContext *s = avctx->priv_data;

    ff_mdct_end(&s->mdct_ctx);
#ifdef LOWPASS
    ff_lowpass_end(&s->lp);
#endif
    return 0;
}

/*
 * Searching index in table with size table_size, where
 * |val-table[best_idx]| is minimal.
 * It assumes that table elements are in increasing order and uses binary search.
 */
#define find_best_value(val, table, table_size, best_idx) \
{ \
    int first=0, last=table_size-1, mid; \
    while(first<=last){ \
        mid=(first+last)/2; \
        if(val > table[mid]){ \
            first = mid + 1; \
        }else{ \
            last = mid - 1; \
        } \
    } \
    if(!first || (first!=table_size && table[first]-val < val-table[last])) \
        best_idx = first; \
    else \
        best_idx = last; \
}

static void encode_block(NellyMoserEncodeContext *s,
                         unsigned char *buf, int buf_size, float *samples)
{
    PutBitContext pb;
    int i, band, block, best_idx, power_idx = 0;
    float power_val, power_candidate, coeff, coeff_sum;
    int band_start, band_end;

    apply_mdct(s, samples, s->mdct_out);
    apply_mdct(s, samples + NELLY_BUF_LEN, s->mdct_out + NELLY_BUF_LEN);

    init_put_bits(&pb, buf, buf_size * 8);

    band_start = 0;
    band_end = ff_nelly_band_sizes_table[0];
    for (band = 0; band < NELLY_BANDS; band++) {
        coeff_sum = 0;
        for (i = band_start; i < band_end; i++) {
            for (block = 0; block < 2; block++) {
                coeff = s->mdct_out[i + block * NELLY_BUF_LEN];
                coeff_sum += coeff * coeff;
            }
        }
        power_candidate =
            (log(FFMAX(64.0, coeff_sum / (ff_nelly_band_sizes_table[band] << 1))) -
             log(64.0)) * 1024.0 / M_LN2;

        if (band) {
            power_candidate -= power_idx;
            find_best_value(power_candidate, ff_nelly_delta_table, 32, best_idx);
            put_bits(&pb, 5, best_idx);
            power_idx += ff_nelly_delta_table[best_idx];
        } else {
            //base exponent
            find_best_value(power_candidate, ff_nelly_init_table, 64, best_idx);
            put_bits(&pb, 6, best_idx);
            power_idx = ff_nelly_init_table[best_idx];
        }

        if (power_idx >= 0) {
            power_val = pow_table[power_idx & 0x7FF] / (1 << (power_idx >> 11));
        } else {
            power_val = -pow(2, -power_idx / 2048.0 - 3.0);
        }
        for (i = band_start; i < band_end; i++) {
            s->mdct_out[i] *= power_val;
            s->mdct_out[i + NELLY_BUF_LEN] *= power_val;
            s->pows[i] = power_idx;
        }
        band_start = band_end;
        if (band != NELLY_BANDS - 1)
            band_end += ff_nelly_band_sizes_table[band + 1];
    }

    ff_nelly_get_sample_bits(s->pows, s->bits);

    for (block = 0; block < 2; block++) {
        for (i = 0; i < NELLY_FILL_LEN; i++) {
            if (s->bits[i] > 0) {
                coeff = s->mdct_out[block * NELLY_BUF_LEN + i];
                find_best_value(coeff,
                                ff_nelly_dequantization_table + (1 << s->bits[i]) - 1,
                                1 << s->bits[i], best_idx);
                put_bits(&pb, s->bits[i], best_idx);
            }
        }
        if (!block)
            put_bits(&pb, NELLY_HEADER_BITS + NELLY_DETAIL_BITS - put_bits_count(&pb), 0);
    }

}

static int encode_tag(AVCodecContext *avctx, uint8_t *frame, int buf_size, void *data)
{
    NellyMoserEncodeContext *s = avctx->priv_data;
    int16_t *samples = data;

    if (s->last_frame)
        return 0;

    if (data) {
#ifdef LOWPASS
        ff_lowpass_filter(&s->lp, samples, s->buf + s->bufsize, avctx->frame_size);
#else
        {
            int i;
            for (i = 0; i < avctx->frame_size; i++) {
                s->buf[i + s->bufsize] = samples[i];
            }
        }
#endif
        s->bufsize += avctx->frame_size;
    } else {
        memset(s->buf + s->bufsize, 0, sizeof(s->buf[0]) * (3 * NELLY_BUF_LEN - s->bufsize));
        s->bufsize = 3 * NELLY_BUF_LEN;
        s->last_frame = 1;
    }

    if (s->bufsize >= 3 * NELLY_BUF_LEN) {
        encode_block(s, frame, buf_size, s->buf);
        memmove(s->buf, s->buf + NELLY_SAMPLES, sizeof(s->buf[0]) * (s->bufsize - NELLY_SAMPLES));
        s->bufsize -= NELLY_SAMPLES;
        return NELLY_BLOCK_LEN;
    }
    return 0;
}

AVCodec nellymoser_encoder = {
    .name = "nellymoser",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_NELLYMOSER,
    .priv_data_size = sizeof(NellyMoserEncodeContext),
    .init = encode_init,
    .encode = encode_tag,
    .close = encode_end,
    .capabilities = CODEC_CAP_SMALL_LAST_FRAME | CODEC_CAP_DELAY,
    .long_name = NULL_IF_CONFIG_SMALL("Nellymoser Asao Codec"),
};
