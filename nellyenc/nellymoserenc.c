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
#include "nellymoser.h"
#include "avcodec.h"
#include "libavutil/random.h"
#include "dsputil.h"


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
    float           state[128];

    int16_t buf[1024*64]; //FIXME (use any better solution)
    int bufsize;

    AVRandomState   random_state;
    int             add_bias;
    float           scale_bias;
    DSPContext      dsp;
    MDCTContext     mdct_ctx;
    DECLARE_ALIGNED_16(float,mdct_tmp[NELLY_BUF_LEN*2]);
    DECLARE_ALIGNED_16(float,mdct_out[NELLY_BUF_LEN*2]);
} NellyMoserEncodeContext;

static DECLARE_ALIGNED_16(float,sine_window2[2*NELLY_BUF_LEN]);

void apply_mdct(NellyMoserEncodeContext *s, float *in, float *coefs)
{
    DECLARE_ALIGNED_16(float,in_buff[NELLY_SAMPLES]);

    memcpy(&in_buff[0], &in[0], NELLY_SAMPLES*sizeof(float));
    s->dsp.vector_fmul(in_buff,sine_window2,NELLY_SAMPLES);
    memset(coefs, 0, NELLY_BUF_LEN*sizeof(float));
    ff_mdct_calc(&s->mdct_ctx, coefs, in_buff, s->mdct_tmp);
}


static av_cold int encode_init(AVCodecContext * avctx) {
    NellyMoserEncodeContext *s = avctx->priv_data;
    int i;

    //TODO check bitrate
    if(avctx->channels!=1){
        av_log(avctx, AV_LOG_ERROR, "Nellymoser supports only 1 channel\n");
        return -1;
    }
    av_log(avctx, AV_LOG_DEBUG, "bit_rate=%i\n", avctx->bit_rate);

    avctx->frame_size = NELLY_SAMPLES;

    s->avctx = avctx;
    av_init_random(0, &s->random_state);
    ff_mdct_init(&s->mdct_ctx, 8, 0);

    dsputil_init(&s->dsp, avctx);

    if(s->dsp.float_to_int16 == ff_float_to_int16_c) {
        s->add_bias = 385;
        s->scale_bias = 1.0/(8*32768);
    } else {
        s->add_bias = 0;
        s->scale_bias = 1.0/(1*8);
    }

    /* Generate overlap window */
    if (!sine_window2[0])
        for (i=0 ; i<256; i++) {
            sine_window2[i] = sin((i + 0.5) / 256.0 * M_PI) /8;
        }

    s->bufsize = 0;
    return 0;
}

static av_cold int encode_end(AVCodecContext * avctx) {
    NellyMoserEncodeContext *s = avctx->priv_data;

    ff_mdct_end(&s->mdct_ctx);
    return 0;
}

static void encode_block(NellyMoserEncodeContext *s,
        unsigned char *buf, int buf_size, int16_t *samples){
    PutBitContext2 pb;
    int bits[NELLY_BUF_LEN];
    int i, j, k, l, b;
    int bs, bk;
    float pows[NELLY_FILL_LEN];
    float *bptr, val;
    float pval;
    float best, tmp, stmp;
    int off[NELLY_BANDS], base; // -1
    float ss=1./1.;
    float scale=1./1.;

    for(i=0; i<NELLY_BUF_LEN*3; i++){
        s->float_buf[i] = samples[i]*ss;
    }
    apply_mdct(s, s->float_buf, s->mdct_out);
    apply_mdct(s, s->float_buf+NELLY_BUF_LEN, s->mdct_out+NELLY_BUF_LEN);

    for(i=0; i<NELLY_SAMPLES; i++){
        av_log(s->avctx, AV_LOG_DEBUG, "%3i: %f (%f)\n", i, s->mdct_out[i],
                s->float_buf[i]);
    }

    init_put_bits2(&pb, buf, buf_size*8);

    /* FIXME: use better (fast) algorithm... */
    best = 1e10;
    bk = 0;
    for(k=0; k<64; k++){
        pval = pow(2,(ff_nelly_init_table[k])/2048) * scale;
        tmp = fabs(((s->mdct_out[0]+s->mdct_out[1])/2)/(pval));
        if(tmp<1) tmp = 1/tmp;
        if(best > tmp){
            bk = k;
            best = tmp;
        }
    }

    av_log(s->avctx, AV_LOG_DEBUG, "base_best=%f\n", best  );

    base = bk;
    put_bits2(&pb, 6, base);
    j=0;

    val = ff_nelly_init_table[base];
    for(i=0; i<NELLY_BANDS; i++){
        if (i > 0){
            /* FIXME: use better (fast) algorithm... */
            best = 1e10;
            bk = 0;
            for(k=0; k<32; k++){
                stmp = 0;
                pval = pow(2,(val+ff_nelly_delta_table[k])/2048) * scale;

                for(b=0; b<2; b++){
                    for(l=0; l<ff_nelly_band_sizes_table[i]; l++){
                        tmp = fabs(s->mdct_out[j+l+b*NELLY_BUF_LEN]/(pval));
                        if(tmp<0.90) tmp=0.90;
                        if(tmp<1) tmp = 1/tmp;
                        if(tmp>5) tmp=5;
                    }
                    stmp += tmp*tmp;
                }


                if(best > stmp){
                    bk = k;
                    best = stmp;
                }
            }
            av_log(s->avctx, AV_LOG_DEBUG, "best=%f\n", best  );
            off[i] = bk;
            put_bits2(&pb, 5, off[i]); // [5] -> -1
            j+=ff_nelly_band_sizes_table[i];
            val += ff_nelly_delta_table[off[i]];
        }
    }

    val = ff_nelly_init_table[base];
    bptr = pows;
    for (i=0 ; i<NELLY_BANDS ; i++) {
        if (i > 0)
            val += ff_nelly_delta_table[off[i]];
        for (j = 0; j < ff_nelly_band_sizes_table[i]; j++) {
            *bptr++ = val;
            av_log(s->avctx, AV_LOG_DEBUG, " exp: %f\n",
                -pow(2,val/2048.)*scale);
        }
    }

    ff_nelly_get_sample_bits(pows, bits);

    bs=0;
    for (i = 0; i < 2; i++) { //2

        for (j = 0; j < NELLY_FILL_LEN; j++) {
            bs+=bits[j];
            if (bits[j] > 0) {
                /* FIXME: use better (fast) algorithm... */
                av_log(s->avctx, AV_LOG_DEBUG, "bits=%i\n", bits[j]);
                pval = -pow(2, pows[j]/2048)*scale;
                av_log(s->avctx, AV_LOG_DEBUG, "pval=%f\n", pval);
                best = 1e10;
                bk = 0;
                for(k=0; k<(1<<bits[j]); k++){
                    tmp = s->mdct_out[i*NELLY_BUF_LEN + j]
                        /(ff_nelly_dequantization_table[(1<<bits[j])-1+k]*pval);
                    if(tmp<0) continue;
                    if(tmp<1) tmp = 1/tmp;
                    if(best > tmp){
                        bk = k;
                        best = tmp;
                    }
                }
                av_log(s->avctx, AV_LOG_DEBUG, "bk=%3i best_deq=%f send=%f\n",
                        bk,
                        best,
                        (ff_nelly_dequantization_table[(1<<bits[j])-1+bk]*pval)
                      );
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

    av_log(s->avctx, AV_LOG_DEBUG, "bs=%i\n", bs);
}

static int encode_tag(AVCodecContext *avctx,
        unsigned char *buf, int buf_size, void *data){
    NellyMoserEncodeContext *s = avctx->priv_data;
    int16_t *samples = data;
    int k, i;
    int n = data ? avctx->frame_size : 0;

    for(k=0; k<n; k++){
        s->buf[k+s->bufsize]=samples[k];
    }
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
    .long_name = "Nellymoser Asao Codec"
};
