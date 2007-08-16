/*
 * EAC3 decoder
 * Copyright (c) 2007 Bartlomiej Wolowiec
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

#include "avcodec.h"
#include "eac3.h"
#include "ac3dec.h"
#include "ac3.h"

static void do_imdct_256(EAC3Context *ctx, int ch)
{
    int k;
    float x[128];
    FFTComplex z[2][64];
    float *o_ptr = ctx->tmp_output;
    int i;

    for(i=0; i<2; i++) {
        /* de-interleave coefficients */
        for(k=0; k<128; k++) {
            x[k] = ctx->transform_coeffs[ch][2*k+i];
        }

        /* run standard IMDCT */
        ctx->imdct_256.fft.imdct_calc(&ctx->imdct_256, o_ptr, x, ctx->tmp_imdct);

        /* reverse the post-rotation & reordering from standard IMDCT */
        for(k=0; k<32; k++) {
            z[i][32+k].re = -o_ptr[128+2*k];
            z[i][32+k].im = -o_ptr[2*k];
            z[i][31-k].re =  o_ptr[2*k+1];
            z[i][31-k].im =  o_ptr[128+2*k+1];
        }
    }

    /* apply AC-3 post-rotation & reordering */
    for(k=0; k<64; k++) {
        o_ptr[    2*k  ] = -z[0][   k].im;
        o_ptr[    2*k+1] =  z[0][63-k].re;
        o_ptr[128+2*k  ] = -z[0][   k].re;
        o_ptr[128+2*k+1] =  z[0][63-k].im;
        o_ptr[256+2*k  ] = -z[1][   k].re;
        o_ptr[256+2*k+1] =  z[1][63-k].im;
        o_ptr[384+2*k  ] =  z[1][   k].im;
        o_ptr[384+2*k+1] = -z[1][63-k].re;
    }
}

/**
 * Performs Inverse MDCT transform
 */
void ff_eac3_do_imdct(EAC3Context *ctx)
{
    int ch;

    for(ch=1; ch<=ctx->nfchans+ctx->lfeon; ch++) {
        if(ctx->blksw[ch]) {
            /* 256-point IMDCT */
            do_imdct_256(ctx, ch);
        } else {
            /* 512-point IMDCT */
            ctx->imdct_512.fft.imdct_calc(&ctx->imdct_512, ctx->tmp_output,
                                          ctx->transform_coeffs[ch],
                                          ctx->tmp_imdct);
        }
        /* apply window function, overlap/add output, save delay */
        ctx->dsp.vector_fmul_add_add(ctx->output[ch], ctx->tmp_output,
                                     ctx->window, ctx->delay[ch], 0,
                                     AC3_BLOCK_SIZE, 1);
        ctx->dsp.vector_fmul_reverse(ctx->delay[ch], ctx->tmp_output+256,
                                     ctx->window, AC3_BLOCK_SIZE);
    }
}

static int eac3_decode_frame(AVCodecContext *avctx, void *data, int *data_size,
                                            uint8_t *buf, int buf_size){
    int16_t *out_samples = (int16_t *)data;
    EAC3Context *c = (EAC3Context *)avctx->priv_data;
    int k, i, blk, ch;
    GetBitContext gbc;

    c->gbc = &gbc;
    c->syncword = 0;

    init_get_bits(&gbc, buf, buf_size*8);
    ff_eac3_parse_syncinfo(&gbc, c);

    if(c->syncword != 0x0B77)
        return -1;

    if(ff_eac3_parse_bsi(&gbc, c) || ff_eac3_parse_audfrm(&gbc, c))
        return -1;

    if(c->fscod == 3){
        assert(c->fscod != 3);
        avctx->sample_rate = ff_ac3_freqs[c->fscod2] / 2;
    }else{
        avctx->sample_rate = ff_ac3_freqs[c->fscod];
    }


    avctx->bit_rate = (c->frmsiz * (avctx->sample_rate) * 16 / ( ff_eac3_blocks[c->numblkscod] * 256)) / 1000;
#ifdef DEBUG
    av_log(NULL, AV_LOG_INFO, "bitrate = %i\n", avctx->bit_rate);
#endif
    avctx->channels = c->nfchans + c->lfeon; // TODO lfe

    for(blk = 0; blk < ff_eac3_blocks[c->numblkscod]; blk++)
    {
        for(i=0; i<AC3_MAX_CHANNELS+1; i++){
            c->deltbae[i] = DBA_NONE;
            c->deltnseg[i] = 0;
        }
#ifdef DEBUG
    av_log(NULL, AV_LOG_INFO, "-------START BLK-------\n");
#endif
        if(ff_eac3_parse_audblk(&gbc, c, blk)){
            av_log(c->avctx, AV_LOG_ERROR, "Error in ff_eac3_parse_audblk\n");
            return -1;
        }
#ifdef DEBUG
    av_log(NULL, AV_LOG_INFO, "-------END BLK-------\n");
#endif

    /* recover coefficients if rematrixing is in use */
    if(c->acmod == AC3_ACMOD_STEREO)
        ff_ac3_do_rematrixing(c->transform_coeffs,
                FFMIN(c->endmant[1], c->endmant[2]),
                c->nrematbnds, c->rematflg);
        //TODO downmix_scaling...

        /* apply scaling to coefficients (dialnorm, dynrng) */
        for(ch=1; ch<=c->nfchans + c->lfeon; ch++) {
            float gain=2.0f;
            if(c->acmod == AC3_ACMOD_DUALMONO) {
                gain *= c->dialnorm[ch-1] * ff_ac3_dynrng_tbl[c->dynrng[ch-1]];
            } else {
                gain *= c->dialnorm[0] * ff_ac3_dynrng_tbl[c->dynrng[0]];
            }
            for(i=0; i<c->endmant[ch]; i++) {
                c->transform_coeffs[ch][i] *= gain;
            }
        }

        ff_eac3_do_imdct(c);
        //TODO downmix

#ifdef DEBUG
        av_log(avctx, AV_LOG_INFO, "channels = %i\n", avctx->channels);
#endif

        // set output mode
        c->blkoutput = 0;
        if (avctx->channels == 1) {
            c->blkoutput |= AC3_OUTPUT_MONO;
        } else if (avctx->channels == 2) {
            c->blkoutput |= AC3_OUTPUT_STEREO;
        } else {
            if (avctx->channels && avctx->channels < c->nfchans + c->lfeon )
                av_log(avctx, AV_LOG_INFO, "ac3_decoder: E-AC3 Source Channels Are Less Then Specified %d: Output to %d Channels\n",avctx->channels, c->nfchans + c->lfeon);
            c->blkoutput |= AC3_OUTPUT_UNMODIFIED;
            if (c->lfeon)
                c->blkoutput |= AC3_OUTPUT_LFEON;
            avctx->channels = c->nfchans + c->lfeon;
        }


        // convert float to 16-bit integer
       for(ch = 1; ch<=c->nfchans + c->lfeon; ch++) { // <- out_channels TODO
            for(i=0; i<AC3_BLOCK_SIZE; i++) {
                c->output[ch][i] = c->output[ch][i] * c->mul_bias +
                    c->add_bias;
            }
            c->dsp.float_to_int16(c->int_output[ch], c->output[ch],
                    AC3_BLOCK_SIZE);
        }
        for (k = 0; k < AC3_BLOCK_SIZE; k++) {
            for (i = 1; i <= avctx->channels; i++) {
                *(out_samples++) = c->int_output[i][k];
            }
        }

    }

#ifdef DEBUG
    av_log(NULL, AV_LOG_INFO, "--------------------------------------------------------------------------\n");
#endif

    *data_size = ff_eac3_blocks[c->numblkscod] * 256 * avctx->channels * sizeof (int16_t); // TODO is ok?

    return buf_size;
}

static int eac3_decode_init(AVCodecContext *avctx){
    int ch;
    EAC3Context *ctx = avctx->priv_data;

    ctx->avctx = avctx;
    ac3_common_init();
    ff_ac3_tables_init();
    av_init_random(0, &ctx->dith_state);
    ff_mdct_init(&ctx->imdct_256, 8, 1);
    ff_mdct_init(&ctx->imdct_512, 9, 1);
    dsputil_init(&ctx->dsp, avctx);
    if(ctx->dsp.float_to_int16 == ff_float_to_int16_c) {
        ctx->add_bias = 385.0f;
        ctx->mul_bias = 1.0f;
    } else {
        ctx->add_bias = 0.0f;
        ctx->mul_bias = 32767.0f;
    }
    ff_ac3_window_init(ctx->window);
    for(ch=0; ch<AC3_MAX_CHANNELS; ch++) {
        memset(ctx->delay[ch], 0, sizeof(ctx->delay[ch]));
    }
    memset(ctx->strtmant, 0, sizeof(int)*MAX_CHANNELS);
    return 0;
}

static int eac3_decode_end(AVCodecContext *avctx){
    EAC3Context *ctx = avctx->priv_data;
    ff_mdct_end(&ctx->imdct_512);
    ff_mdct_end(&ctx->imdct_256);

    return 0;
}

AVCodec eac3_decoder = {
    .name = "E-AC3",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_EAC3,
    .priv_data_size = sizeof (EAC3Context),
    .init = eac3_decode_init,
    .close = eac3_decode_end,
    .decode = eac3_decode_frame,

};
