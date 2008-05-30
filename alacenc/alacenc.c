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

#define ALAC_EXTRADATA_SIZE       36
#define ALAC_FRAME_HEADER_SIZE    55
#define ALAC_FRAME_FOOTER_SIZE    3

typedef struct AlacEncodeContext {
    int channels;
    int samplerate;
    int compression_level;
    int max_coded_frame_size;
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

static void write_frame_header(AlacEncodeContext *s, PutBitContext *pbctx)
{
    put_bits(pbctx, 3,  s->channels-1);         // No. of channels -1
    put_bits(pbctx, 16, 0);                     // Seems to be zero
    put_bits(pbctx, 1,  1);                     // Sample count is in the header
    put_bits(pbctx, 2,  0);                     // FIXME: Wasted bytes field
    put_bits(pbctx, 1,  1);                     // Audio block is verbatim
    put_bits(pbctx, 32, s->avctx->frame_size);  // No. of samples in the frame
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

    s->max_coded_frame_size = (ALAC_FRAME_HEADER_SIZE + ALAC_FRAME_FOOTER_SIZE +
                               avctx->frame_size*s->channels*avctx->bits_per_sample)>>3;

    AV_WB32(alac_extradata,    ALAC_EXTRADATA_SIZE);
    AV_WB32(alac_extradata+4,  MKBETAG('a','l','a','c'));
    AV_WB32(alac_extradata+12, avctx->frame_size);
    AV_WB8 (alac_extradata+17, avctx->bits_per_sample);
    AV_WB8 (alac_extradata+21, s->channels);
    AV_WB32(alac_extradata+24, s->max_coded_frame_size);
    AV_WB32(alac_extradata+28, s->samplerate*s->channels*avctx->bits_per_sample); // average bitrate
    AV_WB32(alac_extradata+32, s->samplerate);

    avctx->extradata = alac_extradata;
    avctx->extradata_size = ALAC_EXTRADATA_SIZE;

    avctx->coded_frame = avcodec_alloc_frame();
    avctx->coded_frame->key_frame = 1;

    s->avctx = avctx;
    return 0;
}

static int alac_encode_frame(AVCodecContext *avctx, uint8_t *frame,
                             int buf_size, void *data)
{
    PutBitContext pb;
    AlacEncodeContext *s = avctx->priv_data;
    int16_t *samples;
    int i, ch;

    if(buf_size < s->max_coded_frame_size) {
        av_log(avctx, AV_LOG_ERROR, "buffer size is too small\n");
        return -1;
    }

    init_put_bits(&pb, frame, buf_size);
    write_frame_header(s, &pb);

    for(ch=0; ch<s->channels; ch++) {
        samples = (int16_t *)data + ch;
        for(i=0; i<avctx->frame_size; i++) {
            put_sbits(&pb, 16, *samples);
            samples += s->channels;
        }
    }

    put_bits(&pb, 3, 7);
    flush_put_bits(&pb);
    return(put_bits_count(&pb)>>3);
}

static av_cold int alac_encode_close(AVCodecContext *avctx)
{
    av_freep(&avctx->extradata);
    avctx->extradata_size = 0;
    av_freep(&avctx->coded_frame);
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
