/*
 * DCA encoder
 * Copyright (C) 2008 Alexander E. Patrakov
 *
 * This file is not yet part of FFmpeg.
 *
 * When this file is part of FFmpeg it can be licensed under the LGPL version 2 or later
 */

/* FFmpeg port by Benjamin Larsson */


#include "avcodec.h"
#include "bitstream.h"
#include "dcaenc.h"

#define MAX_CHANNELS (2)
#define DCA_SUBBANDS_32 (32)
#define DCA_MAX_FRAME_SIZE (16383)
#define FIXED_FRAME_SIZE (7167)

typedef struct {
    PutBitContext pb;
    int32_t history[MAX_CHANNELS][512]; /* This is a circular buffer */
    int start[MAX_CHANNELS];

    int32_t pcm[DCA_SUBBANDS_32];
    int32_t subband[64][MAX_CHANNELS][DCA_SUBBANDS_32]; /* [sample][channel][subband] */
    int16_t st_samples[4096];
    int     fill_samples;
} DCAContext;

static int32_t cos_table[128];


static inline int32_t mul32(int32_t a, int32_t b)
{
    /* on >=i686, gcc compiles this into a single "imull" instruction */
    int64_t r = (int64_t)a * b;
    /* round the result before truncating - improves accuracy */
    return (r + 0x80000000) >> 32;
}

/* Integer version of the cosine modulated Pseudo QMF */

void qmf_init(void)
{
    int i;
    int32_t c[17], s[17];
    s[0] = 0;       /* sin(index * PI / 64) * 0x7fffffff */
    c[0] = 0x7fffffff;  /* cos(index * PI / 64) * 0x7fffffff */

    for (i = 1; i <= 16; i++) {
        s[i] = 2 * (mul32(c[i-1], 105372028) + mul32(s[i-1], 2144896908));
        c[i] = 2 * (mul32(c[i-1], 2144896908) - mul32(s[i-1], 105372028));
    }

    for (i = 0; i < 16; i++) {
        cos_table[i] = c[i] >> 3; /* so that the output doesn't overflow */
        cos_table[i+16] = s[16-i] >> 3;
        cos_table[i+32] = -s[i] >> 3;
        cos_table[i+48] = -c[16-i] >> 3;
        cos_table[i+64] = -c[i] >> 3;
        cos_table[i+80] = -s[16-i] >> 3;
        cos_table[i+96] = s[i] >> 3;
        cos_table[i+112] = c[16-i] >> 3;
    }
}

static int32_t band_delta_factor(int band, int sample_num)
{
    int index = band * (2 * sample_num + 1);
    if (band == 0)
        return 0x07ffffff;
    else
        return cos_table[index & 127];
}

static void qmf_decompose(DCAContext *c, int32_t in[32], int32_t out[32], int channel)
{
    int band, i, j, k;
    int32_t resp;
    int32_t accum[DCA_SUBBANDS_32];

    /* Place new samples into the history buffer */
    for (i = 0; i < DCA_SUBBANDS_32; i++)
        c->history[channel][c->start[channel] + i] = in[i];
    c->start[channel] += DCA_SUBBANDS_32;
    if (c->start[channel] == 512)
        c->start[channel] = 0;

    /* Calculate the dot product of the signal with the (possibly inverted)
       reference decoder's response to this vector:
       (0.0, 0.0, ..., 0.0, -1.0, 1.0, 0.0, ..., 0.0)
       so that -1.0 cancels 1.0 from the previous step */

    memset(accum,0,sizeof(int32_t));

    for (k = 48, j = 0, i = c->start[channel]; i < 512; k++, j++, i++)
        accum[(k & 32) ? (31 - (k & 31)) : (k & 31)] += mul32(c->history[channel][i], UnQMF[j]);
    for (i = 0; i < c->start[channel]; k++, j++, i++)
        accum[(k & 32) ? (31 - (k & 31)) : (k & 31)] += mul32(c->history[channel][i], UnQMF[j]);

    resp = 0;
    /* TODO: implement FFT instead of this naive calculation */
    for (band = 0; band < DCA_SUBBANDS_32; band++) {
        for (j = 0; j < 32; j++)
            resp += mul32(accum[j], band_delta_factor(band, j));

        out[band] = (band & 2) ? (-resp) : resp;
    }
}

static void put_frame_header(DCAContext *c)
{
    /* SYNC */
    put_bits(&c->pb, 16, 0x7ffe);
    put_bits(&c->pb, 16, 0x8001);

    /* Frame type: normal */
    put_bits(&c->pb, 1, 1);

    /* Deficit sample count: none */
    put_bits(&c->pb, 5, 31);

    /* CRC is not present */
    put_bits(&c->pb, 1, 0);

    /* Number of PCM sample blocks: 64
       (larger values are unusable with 1:1 compression due to high bitrate
       and frame size limitation) */
    put_bits(&c->pb, 7, 63);

    /* Primary frame byte size: 7168 */
    put_bits(&c->pb, 14, FIXED_FRAME_SIZE);

    /* Audio channel arrangement: L + R (stereo) */
    put_bits(&c->pb, 6, 2);

    /* Core audio sampling frequency: 44100 Hz */
    put_bits(&c->pb, 4, 8);

    /* Transmission bit rate: 1411.2 kbps */
    put_bits(&c->pb, 5, 0x16);

    /* Embedded down mix: disabled */
    put_bits(&c->pb, 1, 0);

    /* Embedded dynamic range flag: not present */
    put_bits(&c->pb, 1, 0);

    /* Embedded time stamp flag: not present */
    put_bits(&c->pb, 1, 0);

    /* Auxiliary data flag: not present */
    put_bits(&c->pb, 1, 0);

    /* HDCD source: no */
    put_bits(&c->pb, 1, 0);

    /* Extension audio ID: N/A */
    put_bits(&c->pb, 3, 0);

    /* Extended audio data: not present */
    put_bits(&c->pb, 1, 0);

    /* Audio sync word insertion flag: after each sub-frame */
    put_bits(&c->pb, 1, 0);

    /* Low frequency effects flag: not present */
    put_bits(&c->pb, 2, 0);

    /* Predictor history switch flag: on */
    put_bits(&c->pb, 1, 1);

    /* No CRC */
    /* Multirate interpolator switch: non-perfect reconstruction */
    put_bits(&c->pb, 1, 0);

    /* Encoder software revision: 7 */
    put_bits(&c->pb, 4, 7);

    /* Copy history: 0 */
    put_bits(&c->pb, 2, 0);

    /* Source PCM resolution: 16 bits, not DTS ES */
    put_bits(&c->pb, 3, 0);

    /* Front sum/difference coding: no */
    put_bits(&c->pb, 1, 0);

    /* Surrounds sum/difference coding: no */
    put_bits(&c->pb, 1, 0);

    /* Dialog normalization: 0 dB */
    put_bits(&c->pb, 4, 0);
}

static void put_primary_audio_header(DCAContext *c)
{
    /* Number of subframes: 2 */
    put_bits(&c->pb, 4, 1);

    /* Number of primary audio channels: 2 */
    put_bits(&c->pb, 3, 1);

    /* Subband activity count: 27 + 27 */
    put_bits(&c->pb, 5, 25);
    put_bits(&c->pb, 5, 25);

    /* High frequency VQ start subband: 27, 27 */
    put_bits(&c->pb, 5, 26);
    put_bits(&c->pb, 5, 26);

    /* Joint intensity coding index: 0, 0 */
    put_bits(&c->pb, 3, 0);
    put_bits(&c->pb, 3, 0);

    /* Transient mode codebook: A4, A4 (arbitrary) */
    put_bits(&c->pb, 2, 0);
    put_bits(&c->pb, 2, 0);

    /* Scale factor code book: 7 bit linear, 7-bit sqrt table (for each channel) */
    put_bits(&c->pb, 3, 6);
    put_bits(&c->pb, 3, 6);

    /* Bit allocation quantizer select: linear 5-bit */
    put_bits(&c->pb, 3, 6);
    put_bits(&c->pb, 3, 6);

    /* Quantization index codebook select: dummy data
       to avoid transmission of scale factor adjustment */
    put_bits(&c->pb, 1, 1); put_bits(&c->pb, 1, 1);
    put_bits(&c->pb, 2, 3); put_bits(&c->pb, 2, 3);
    put_bits(&c->pb, 2, 3); put_bits(&c->pb, 2, 3);
    put_bits(&c->pb, 2, 3); put_bits(&c->pb, 2, 3);
    put_bits(&c->pb, 2, 3); put_bits(&c->pb, 2, 3);
    put_bits(&c->pb, 3, 7); put_bits(&c->pb, 3, 7);
    put_bits(&c->pb, 3, 7); put_bits(&c->pb, 3, 7);
    put_bits(&c->pb, 3, 7); put_bits(&c->pb, 3, 7);
    put_bits(&c->pb, 3, 7); put_bits(&c->pb, 3, 7);
    put_bits(&c->pb, 3, 7); put_bits(&c->pb, 3, 7);

    /* Scale factor adjustment index: not transmitted */
}

/* TODO: don't hardcode 16-bit quantization */
static uint32_t quantize(int32_t d)
{
    d = d >> 16;
    return d & 0xffff;
}


static void put_subframe(DCAContext *c, int32_t subband_data[32][2][32])
{
    int i, sub, ss, ch;

    /* Subsubframes count: 4 */
    put_bits(&c->pb, 2, 3);

    /* Partial subsubframe sample count: dummy */
    put_bits(&c->pb, 3, 0);

    /* Prediction mode: no ADPCM, in each channel and subband */
    for (ch = 0; ch < 2; ch++)
        for (sub = 0; sub < 27; sub++)
            put_bits(&c->pb, 1, 0);

    /* Prediction VQ addres: not transmitted */
    /* Bit allocation index: 19 = "16 bits", for each channel and subband */
    for (ch = 0; ch < 2; ch++)
        for (sub = 0; sub < 27; sub++)
            put_bits(&c->pb, 5, 19);

    /* Transition mode: none for each channel and subband */
    for (ch = 0; ch < 2; ch++)
        for (sub = 0; sub < 27; sub++)
            put_bits(&c->pb, 1, 0); /* according to Huffman codebook A4 */

    /* Scale factors: the same for each channel and subband,
       encoded according to Table D.1.2 */
    for (ch = 0; ch < 2; ch++)
        for (sub = 0; sub < 27; sub++)
            put_bits(&c->pb, 7, 110);

    /* Joint subband scale factor codebook select: not transmitted */
    /* Scale factors for joint subband coding: not transmitted */
    /* Stereo down-mix coefficients: not transmitted */
    /* Dynamic range coefficient: not transmitted */
    /* Stde information CRC check word: not transmitted */
    /* VQ encoded high frequency subbands: not transmitted */
    /* LFE data: none */
    /* Audio data: 4 subsubframes */

    for (ss = 0; ss < 4 ; ss++)
        for (ch = 0; ch < 2; ch++)
            for (sub = 0; sub < 27; sub++)
                for (i = 0; i < 8; i++)
                    put_bits(&c->pb, 16, quantize(subband_data[ss * 8 + i][ch][sub]));
    /* DSYNC */
    align_put_bits(&c->pb);
    put_bits(&c->pb, 0xffff, 16);
}

void put_frame(DCAContext *c, int32_t subband_data[64][2][32], uint8_t *frame)
{
    int channel;
    init_put_bits(&c->pb, frame, DCA_MAX_FRAME_SIZE);

    put_frame_header(c);
    put_primary_audio_header(c);
    for (channel=0 ; channel<2; channel++)
        put_subframe(c, &subband_data[32 * channel]);

    flush_put_bits(&c->pb);
}

static int DCA_encode_frame(AVCodecContext *avctx,
                            uint8_t *frame, int buf_size, void *data)
{
    int i,k,channel;
    DCAContext *c = avctx->priv_data;
//    int16_t *samples = data;

//    if (buf_size < MAX_CHANNELS*2048*sizeof(int16_t))
//        return -1;

    // We always get 2048 16bit samples per call so save then until next call

    if (c->fill_samples) {
        memcpy(c->st_samples, data, 2048*sizeof(int16_t));
        c->fill_samples = 0;
        return 0;
    }
    memcpy(&c->st_samples[2048], data, 2048*sizeof(int16_t));
    c->fill_samples = 0;

    for (i = 0; i < 64; i ++) /* i is the decimated sample number */
        for (channel=0; channel<2 ; channel++) {
            /* Get 32 PCM samples */
            for (k = 0; k < 32; k++) { /* k is the sample number in a 32-sample block */
                c->pcm[k] = c->st_samples[4 * (32*i+k) + 2 * channel] << 16;
            }
        /* Put subband samples into the proper place */
        qmf_decompose(c, c->pcm, &c->subband[i][channel][0], channel);
    }

    put_frame(c, c->subband, frame);

    return 4000;
}

static int DCA_encode_init(AVCodecContext *avctx) {
    DCAContext *c = avctx->priv_data;

    if(avctx->channels != 2 || avctx->sample_rate != 44100) {
        av_log(avctx, AV_LOG_ERROR, "Only 44.1 kHz stereo is supported at the moment!\n");
        return -1;
    }

    // Make sure that the first call to encode is a fill call
    c->fill_samples = 1;

    qmf_init();
    return 0;
}

AVCodec dca_encoder = {
    "dca",
    CODEC_TYPE_AUDIO,
    CODEC_ID_DTS,
    sizeof(DCAContext),
    DCA_encode_init,
    DCA_encode_frame,
    NULL,
    NULL,
};
