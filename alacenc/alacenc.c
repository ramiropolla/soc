/**
 * ALAC audio encoder
 * Copyright (c) 2008  Jaikrishnan Menon <realityman@gmx.net>
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
#include "bitstream.h"

#define DEFAULT_FRAME_SIZE        4096
#define DEFAULT_SAMPLE_SIZE       16
#define MAX_CHANNELS              8
#define ALAC_EXTRADATA_SIZE       36
#define ALAC_FRAME_HEADER_SIZE    55
#define ALAC_FRAME_FOOTER_SIZE    3

#define ALAC_ESCAPE_CODE          0x1FF

typedef struct RiceContext {
    int history_mult;
    int initial_history;
    int k_modifier;
    int rice_modifier;
} RiceContext;

typedef struct AlacEncodeContext {
    int channels;
    int samplerate;
    int compression_level;
    int max_coded_frame_size;
    int write_sample_size;
    int16_t *sample_buf[MAX_CHANNELS];
    PutBitContext pbctx;
    RiceContext rc;
    AVCodecContext *avctx;
} AlacEncodeContext;

/**
 * put_sbits
 * @param pb PutBitContext pointer
 * @param bits Number of bits to output
 * @param val data Bits
 */
static void put_sbits(PutBitContext *pb, int bits, int32_t val)
{
    put_bits(pb, bits, val & ((1<<bits)-1));
}

static void allocate_sample_buffers(AlacEncodeContext *s)
{
    int i = s->channels;

    while(i) {
        s->sample_buf[i-1] = av_mallocz(s->avctx->frame_size*s->avctx->bits_per_sample>>3);
        i--;
    }
}

static void free_sample_buffers(AlacEncodeContext *s)
{
    int i = s->channels;

    while(i) {
        av_freep(&s->sample_buf[i-1]);
        i--;
    }
}

static void init_sample_buffers(AlacEncodeContext *s, int16_t *input_samples)
{
    int ch, i;

    for(ch=0;ch<s->channels;ch++) {
        int16_t *sptr = input_samples + ch;
        for(i=0;i<s->avctx->frame_size;i++) {
            s->sample_buf[ch][i] = *sptr;
            sptr += s->channels;
        }
    }
}

static void encode_scalar(AlacEncodeContext *s, int x, int k, int write_sample_size)
{
    int divisor, q, r;

    k = FFMIN(k, s->rc.k_modifier);
    divisor = (1<<k) - 1;
    q = x / divisor;
    r = x % divisor;

    if(q > 8) {
        // write escape code and sample value directly
        put_bits(&s->pbctx, 9, ALAC_ESCAPE_CODE);
        put_bits(&s->pbctx, write_sample_size, x);
    } else {
        if(q)
            put_bits(&s->pbctx, q, (1<<q) - 1);
        put_bits(&s->pbctx, 1, 0);

        if(k != 1) {
            if(r > 0)
                put_bits(&s->pbctx, k, r+1);
            else
                put_bits(&s->pbctx, k-1, 0);
        }
    }
}

static void write_frame_header(AlacEncodeContext *s)
{
    put_bits(&s->pbctx, 3,  s->channels-1);                 // No. of channels -1
    put_bits(&s->pbctx, 16, 0);                             // Seems to be zero
    put_bits(&s->pbctx, 1,  1);                             // Sample count is in the header
    put_bits(&s->pbctx, 2,  0);                             // FIXME: Wasted bytes field
    put_bits(&s->pbctx, 1,  !s->compression_level);         // Audio block is verbatim
    put_bits(&s->pbctx, 32, s->avctx->frame_size);          // No. of samples in the frame
}

static void alac_entropy_coder(AlacEncodeContext *s, int16_t *samples)
{
    unsigned int history = s->rc.initial_history;
    int sign_modifier = 0, i = 0, k;

    while(i < s->avctx->frame_size) {
        int x;

        k = av_log2((history >> 9) + 3);

        x = -2*(*samples)-1;
        x ^= (x>>31);

        samples++;
        i++;

        encode_scalar(s, x - sign_modifier, k, s->write_sample_size);

        history += x * s->rc.history_mult
                   - ((history * s->rc.history_mult) >> 9);

        sign_modifier = 0;
        if(x > 0xFFFF)
            history = 0xFFFF;

        if((history < 128) && (i < s->avctx->frame_size)) {
            unsigned int block_size = 0;

            sign_modifier = 1;
            k = 7 - av_log2(history) + ((history + 16) >> 6);

            while((*samples == 0) && (i < s->avctx->frame_size)) {
                samples++;
                i++;
                block_size++;
            }
            encode_scalar(s, block_size, k, 16);

            sign_modifier = (block_size <= 0xFFFF);

            history = 0;
        }

    }
}

static void write_compressed_frame(AlacEncodeContext *s)
{
    int i;

    put_bits(&s->pbctx, 8, 0);      // FIXME: interlacing shift
    put_bits(&s->pbctx, 8, 0);      // FIXME: interlacing leftweight

    for(i=0;i<s->channels;i++) {
        put_bits(&s->pbctx, 4, 0);  // prediction type : currently only type 0 has been RE'd
        put_bits(&s->pbctx, 4, 0);  // FIXME: prediction quantization

        put_bits(&s->pbctx, 3, s->rc.rice_modifier);
        put_bits(&s->pbctx, 5, 0);  // predictor order
        // FIXME: predictor coeff. table goes here
    }

    // apply entropy coding to audio samples

    for(i=0;i<s->channels;i++) {
        alac_entropy_coder(s, s->sample_buf[i]);
    }
}

static av_cold int alac_encode_init(AVCodecContext *avctx)
{
    AlacEncodeContext *s    = avctx->priv_data;
    uint8_t *alac_extradata = av_mallocz(ALAC_EXTRADATA_SIZE+1);

    avctx->frame_size      = DEFAULT_FRAME_SIZE;
    avctx->bits_per_sample = DEFAULT_SAMPLE_SIZE;
    s->channels            = avctx->channels;
    s->samplerate          = avctx->sample_rate;

    if(avctx->sample_fmt != SAMPLE_FMT_S16) {
        av_log(avctx, AV_LOG_ERROR, "only pcm_s16 input samples are supported\n");
        return -1;
    }

    // Set default compression level
    if(avctx->compression_level == FF_COMPRESSION_DEFAULT)
        s->compression_level = 1;
    else
        s->compression_level = av_clip(avctx->compression_level, 0, 1);

    // Initialize default Rice parameters
    s->rc.history_mult    = 40;
    s->rc.initial_history = 10;
    s->rc.k_modifier      = 14;
    s->rc.rice_modifier   = 4;

    s->max_coded_frame_size = (ALAC_FRAME_HEADER_SIZE + ALAC_FRAME_FOOTER_SIZE +
                               avctx->frame_size*s->channels*avctx->bits_per_sample)>>3;

    s->write_sample_size  = avctx->bits_per_sample + s->channels - 1; // FIXME: consider wasted_bytes

    AV_WB32(alac_extradata,    ALAC_EXTRADATA_SIZE);
    AV_WB32(alac_extradata+4,  MKBETAG('a','l','a','c'));
    AV_WB32(alac_extradata+12, avctx->frame_size);
    AV_WB8 (alac_extradata+17, avctx->bits_per_sample);
    AV_WB8 (alac_extradata+21, s->channels);
    AV_WB32(alac_extradata+24, s->max_coded_frame_size);
    AV_WB32(alac_extradata+28, s->samplerate*s->channels*avctx->bits_per_sample); // average bitrate
    AV_WB32(alac_extradata+32, s->samplerate);

    // Set relevant extradata fields
    if(s->compression_level > 0) {
        AV_WB8(alac_extradata+18, s->rc.history_mult);
        AV_WB8(alac_extradata+19, s->rc.initial_history);
        AV_WB8(alac_extradata+20, s->rc.k_modifier);
    }

    avctx->extradata = alac_extradata;
    avctx->extradata_size = ALAC_EXTRADATA_SIZE;

    avctx->coded_frame = avcodec_alloc_frame();
    avctx->coded_frame->key_frame = 1;

    s->avctx = avctx;

    allocate_sample_buffers(s);

    return 0;
}

static int alac_encode_frame(AVCodecContext *avctx, uint8_t *frame,
                             int buf_size, void *data)
{
    AlacEncodeContext *s = avctx->priv_data;
    PutBitContext *pb = &s->pbctx;
    int i, ch;

    if(avctx->frame_size > DEFAULT_FRAME_SIZE) {
        av_log(avctx, AV_LOG_ERROR, "input frame size exceeded\n");
        return -1;
    }

    if(buf_size < s->max_coded_frame_size) {
        av_log(avctx, AV_LOG_ERROR, "buffer size is too small\n");
        return -1;
    }

    init_put_bits(pb, frame, buf_size);
    write_frame_header(s);

    if(s->compression_level == 0) {
        // Verbatim mode
        int16_t *samples = data;
        for(i=0; i<avctx->frame_size*s->channels; i++) {
            put_sbits(pb, 16, *samples++);
        }
    } else {
        init_sample_buffers(s, data);
        write_compressed_frame(s);
    }

    put_bits(pb, 3, 7);
    flush_put_bits(pb);
    return(put_bits_count(pb)>>3);
}

static av_cold int alac_encode_close(AVCodecContext *avctx)
{
    AlacEncodeContext *s = avctx->priv_data;

    av_freep(&avctx->extradata);
    avctx->extradata_size = 0;
    av_freep(&avctx->coded_frame);
    free_sample_buffers(s);
    return 0;
}

AVCodec alac_encoder = {
    "alac",
    CODEC_TYPE_AUDIO,
    CODEC_ID_ALAC,
    sizeof(AlacEncodeContext),
    alac_encode_init,
    alac_encode_frame,
    alac_encode_close,
    .capabilities = CODEC_CAP_SMALL_LAST_FRAME,
    .long_name = "ALAC (Apple Lossless Audio Codec)",
};
