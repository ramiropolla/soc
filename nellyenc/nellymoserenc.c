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
#include "lowpass.h"

#define MAX_POW_CACHED (1<<11)
#define LOWPASS 1

/*
 * FIXME: Bitstream from vorbis_enc.c (move to seperate file?)
 */
typedef struct {
    int total;
    int total_pos;
    int pos;
    uint8_t * buf_ptr;
} PutBitContext2;

static inline void init_put_bits2(PutBitContext2 * pb, uint8_t * buf, int buffer_len) {
    pb->total = buffer_len * 8;
    pb->total_pos = 0;
    pb->pos = 0;
    pb->buf_ptr = buf;
}

static void put_bits2(PutBitContext2 * pb, int bits, uint64_t val) {
    if ((pb->total_pos += bits) >= pb->total) return;
    if (!bits) return;
    if (pb->pos) {
        if (pb->pos > bits) {
            *pb->buf_ptr |= val << (8 - pb->pos);
            pb->pos -= bits;
            bits = 0;
        } else {
            *pb->buf_ptr++ |= (val << (8 - pb->pos)) & 0xFF;
            val >>= pb->pos;
            bits -= pb->pos;
            pb->pos = 0;
        }
    }
    for (; bits >= 8; bits -= 8) {
        *pb->buf_ptr++ = val & 0xFF;
        val >>= 8;
    }
    if (bits) {
        *pb->buf_ptr = val;
        pb->pos = 8 - bits;
    }
}

static inline int put_bits2_count(PutBitContext2 * pb) {
    return pb->total_pos;
}

typedef struct NellyMoserEncodeContext {
    AVCodecContext* avctx;
    DECLARE_ALIGNED_16(float,float_buf[2*NELLY_SAMPLES]); // NELLY_SAMPLES

    float buf[1024*64]; //FIXME (use any better solution)
    int bufsize;

    DSPContext      dsp;
    MDCTContext     mdct_ctx;
    float pows[NELLY_FILL_LEN];
    DECLARE_ALIGNED_16(float,mdct_tmp[NELLY_BUF_LEN*2]);
    DECLARE_ALIGNED_16(float,mdct_out[NELLY_BUF_LEN*2]);

    LPFilterContext lp;
} NellyMoserEncodeContext;

static DECLARE_ALIGNED_16(float,sine_window[2*NELLY_BUF_LEN]);
static float pow_table[MAX_POW_CACHED];

void apply_mdct(NellyMoserEncodeContext *s, float *in, float *coefs)
{
    DECLARE_ALIGNED_16(float,in_buff[NELLY_SAMPLES]);

    memcpy(&in_buff[0], &in[0], NELLY_SAMPLES*sizeof(float));
    s->dsp.vector_fmul(in_buff,sine_window,NELLY_SAMPLES);
    memset(coefs, 0, NELLY_BUF_LEN*sizeof(float));
    ff_mdct_calc(&s->mdct_ctx, coefs, in_buff, s->mdct_tmp);
}


static av_cold int encode_init(AVCodecContext * avctx) {
    NellyMoserEncodeContext *s = avctx->priv_data;
    int i;

    if(avctx->channels!=1){
        av_log(avctx, AV_LOG_ERROR, "Nellymoser supports only 1 channel\n");
        return -1;
    }

    switch(avctx->sample_rate){
        case 8000:
            ff_lowpass_init(&s->lp, 8000, 2000, 3800, 1, 60);
            break;
        case 11025:
            ff_lowpass_init(&s->lp, 11025, 3000, 5000, 1, 60);
            break;
        case 22050:
            ff_lowpass_init(&s->lp, 22025, 6000, 10000, 1, 60);
            break;
        case 44100:
            ff_lowpass_init(&s->lp, 44100, 12000, 20000, 1, 60);
            break;
        default:
            av_log(avctx, AV_LOG_ERROR,
                    "Nellymoser works only with 8000, 11025, 22050 and 44100 sample rate\n");
        return -1;
    }

    avctx->frame_size = NELLY_SAMPLES;

    s->avctx = avctx;
    ff_mdct_init(&s->mdct_ctx, 8, 0);

    dsputil_init(&s->dsp, avctx);

    /* Generate overlap window */
    if (!sine_window[0]){
        ff_sine_window_init(sine_window, 128);
        for (i=0; i<128; i++) {
            sine_window[255-i] = sine_window[i];
        }
    }
    for(i=0; i<MAX_POW_CACHED; i++)
        pow_table[i] = -pow(2, -i/2048.0 - 3.0);

    s->bufsize = 0;
    return 0;
}

static av_cold int encode_end(AVCodecContext * avctx) {
    NellyMoserEncodeContext *s = avctx->priv_data;

    ff_mdct_end(&s->mdct_ctx);
    ff_lowpass_end(&s->lp);
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
        unsigned char *buf, int buf_size, float *samples){
    PutBitContext2 pb;
    int bits[NELLY_BUF_LEN];
    int i, j, k, l, b;
    int bk;
    int val=0;
    float pval;
    float tmp, stmp;
    int band_start, band_end;

    for(i=0; i<NELLY_BUF_LEN*3; i++){
        s->float_buf[i] = samples[i];
    }
    apply_mdct(s, s->float_buf, s->mdct_out);
    apply_mdct(s, s->float_buf+NELLY_BUF_LEN, s->mdct_out+NELLY_BUF_LEN);

    for(i=0; i<NELLY_SAMPLES; i++){
        av_log(s->avctx, AV_LOG_DEBUG, "%3i: %f (%f)\n", i, s->mdct_out[i],
                s->float_buf[i]);
    }

    init_put_bits2(&pb, buf, buf_size*8);

    band_start = 0;
    band_end = ff_nelly_band_sizes_table[0];
    for(i=0; i<NELLY_BANDS; i++){
        stmp = 0;
        for(l=band_start; l<band_end; l++){
            for(b=0; b<2; b++){
                tmp = s->mdct_out[l+b*NELLY_BUF_LEN];
                stmp += tmp*tmp;
            }
        }
        tmp = ( log(FFMAX(64.0, stmp/(ff_nelly_band_sizes_table[i]<<1))) - log(64.0)) *
            M_LOG2E * 1024.0;

        if(i){
            tmp -= val;
            find_best_value(tmp, ff_nelly_delta_table, 32, bk);
            put_bits2(&pb, 5, bk);
            val += ff_nelly_delta_table[bk];
        }else{
            //base exponent
            find_best_value(tmp, ff_nelly_init_table, 64, bk);
            put_bits2(&pb, 6, bk);
            val = ff_nelly_init_table[bk];
        }

        if(val >= 0){
            pval = pow_table[val&0x7FF] / (1<<(val>>11)) ;
        }else{
            pval = -pow(2, -val/2048.0 - 3.0);
        }
        for (k = band_start; k < band_end; k++) {
            s->mdct_out[k] *= pval;
            s->mdct_out[k+NELLY_BUF_LEN] *= pval;
            s->pows[k] = val;
        }
        band_start = band_end;
        if(i!=NELLY_BANDS-1)
            band_end += ff_nelly_band_sizes_table[i+1];
    }

    ff_nelly_get_sample_bits(s->pows, bits);

    for (i = 0; i < 2; i++) { //2

        for (j = 0; j < NELLY_FILL_LEN; j++) {
            if (bits[j] > 0) {
                tmp = s->mdct_out[i*NELLY_BUF_LEN + j];

                find_best_value(tmp,
                        (ff_nelly_dequantization_table + (1<<bits[j])-1),
                        (1<<bits[j]), bk);
                put_bits2(&pb, bits[j], bk);
            }
        }
        av_log(s->avctx, AV_LOG_DEBUG, "count=%i (%i)\n",
                put_bits2_count(&pb),
                NELLY_HEADER_BITS + NELLY_DETAIL_BITS
                );
        if(!i)
            put_bits2(&pb, NELLY_HEADER_BITS + NELLY_DETAIL_BITS - put_bits2_count(&pb) , 0);

        av_log(s->avctx, AV_LOG_DEBUG, "count=%i (%i)\n",
                put_bits2_count(&pb),
                NELLY_HEADER_BITS + NELLY_DETAIL_BITS
                );
    }

}

static int encode_tag(AVCodecContext *avctx,
        unsigned char *buf, int buf_size, void *data){
    NellyMoserEncodeContext *s = avctx->priv_data;
    int16_t *samples = data;
    int k, i;
    int n = data ? avctx->frame_size : 0;

#if LOWPASS
    ff_lowpass_filter(&s->lp, samples, s->buf+s->bufsize, n);
#else
    for(k=0; k<n; k++){
        s->buf[k+s->bufsize]=samples[k];
    }
#endif
    s->bufsize+=n;

    /* FIXME:
     *  find better method for it...
     */
    for(k=0; (s->bufsize>=3*NELLY_BUF_LEN) && ((k+1)*NELLY_BLOCK_LEN<buf_size); k++){
        encode_block(s, buf+NELLY_BLOCK_LEN*k, buf_size, s->buf);

        //memmove(s->buf, s->buf+NELLY_SAMPLES, s->bufsize-NELLY_SAMPLES);
        for(i=NELLY_SAMPLES; i<s->bufsize; i++){
            s->buf[i-NELLY_SAMPLES] = s->buf[i];
        }
        s->bufsize-=NELLY_SAMPLES;
    }
    return NELLY_BLOCK_LEN*k;
}


AVCodec nellymoser_encoder = {
    .name = "nellymoser",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_NELLYMOSER,
    .priv_data_size = sizeof(NellyMoserEncodeContext),
    .init = encode_init,
    .encode = encode_tag,
    .close = encode_end,
    .long_name = NULL_IF_CONFIG_SMALL("Nellymoser Asao Codec"),
};
