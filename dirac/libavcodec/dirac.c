/* -*-  indent-tabs-mode:nil; c-basic-offset:4;  -*- */
/*
 * Copyright (C) 2007 Marco Gerards <marco@gnu.org>
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
 * @file
 * Dirac Decoder
 * @author Marco Gerards <marco@gnu.org>
 */

#include "dirac.h"
#include "avcodec.h"
#include "dsputil.h"
#include "get_bits.h"
#include "bytestream.h"
#include "golomb.h"
#include "dirac_arith.h"
#include "mpeg12data.h"

/* defaults for source parameters */
static const dirac_source_params dirac_source_parameters_defaults[] =
{
    { 640,  480,  2, 0, 0, 1,  1, 640,  480,  0, 0, 1, 0 },
    { 176,  120,  2, 0, 0, 9,  2, 176,  120,  0, 0, 1, 1 },
    { 176,  144,  2, 0, 1, 10, 3, 176,  144,  0, 0, 1, 2 },
    { 352,  240,  2, 0, 0, 9,  2, 352,  240,  0, 0, 1, 1 },
    { 352,  288,  2, 0, 1, 10, 3, 352,  288,  0, 0, 1, 2 },
    { 704,  480,  2, 0, 0, 9,  2, 704,  480,  0, 0, 1, 1 },
    { 704,  576,  2, 0, 1, 10, 3, 704,  576,  0, 0, 1, 2 },
    { 720,  480,  1, 1, 0, 4,  2, 704,  480,  8, 0, 3, 1 },
    { 720,  576,  1, 1, 1, 3,  3, 704,  576,  8, 0, 3, 2 },

    { 1280, 720,  1, 0, 1, 7,  1, 1280, 720,  0, 0, 3, 3 },
    { 1280, 720,  1, 0, 1, 6,  1, 1280, 720,  0, 0, 3, 3 },
    { 1920, 1080, 1, 1, 1, 4,  1, 1920, 1080, 0, 0, 3, 3 },
    { 1920, 1080, 1, 1, 1, 3,  1, 1920, 1080, 0, 0, 3, 3 },
    { 1920, 1080, 1, 0, 1, 7,  1, 1920, 1080, 0, 0, 3, 3 },
    { 1920, 1080, 1, 0, 1, 6,  1, 1920, 1080, 0, 0, 3, 3 },
    { 2048, 1080, 0, 0, 1, 2,  1, 2048, 1080, 0, 0, 4, 4 },
    { 4096, 2160, 0, 0, 1, 2,  1, 4096, 2160, 0, 0, 4, 4 },

    { 3840, 2160, 1, 0, 1, 7,  1, 3840, 2160, 0, 0, 3, 3 },
    { 3840, 2160, 1, 0, 1, 6,  1, 3840, 2160, 0, 0, 3, 3 },
    { 7680, 4320, 1, 0, 1, 7,  1, 3840, 2160, 0, 0, 3, 3 },
    { 7680, 4320, 1, 0, 1, 6,  1, 3840, 2160, 0, 0, 3, 3 },
};

static const AVRational dirac_preset_aspect_ratios[] =
{
    {1, 1},
    {10, 11},
    {12, 11},
    {40, 33},
    {16, 11},
    {4, 3},
};

static const AVRational dirac_frame_rate[] =
{
    {15000, 1001},
    {25, 2},
};

static const dirac_pixel_range dirac_pixel_range_presets[] = {
    { 0,   255,  128,  255  },
    { 16,  219,  128,  224  },
    { 64,  876,  512,  896  },
    { 256, 3504, 2048, 3584 },
};

static const color_specification dirac_color_spec_presets[] = {
    { COLOR_PRIMARY_HDTV,     COLOR_MATRIX_HDTV, TRANSFER_FUNC_TV },
    { COLOR_PRIMARY_SDTV_525, COLOR_MATRIX_SDTV, TRANSFER_FUNC_TV },
    { COLOR_PRIMARY_SDTV_625, COLOR_MATRIX_SDTV, TRANSFER_FUNC_TV },
    { COLOR_PRIMARY_HDTV,     COLOR_MATRIX_HDTV, TRANSFER_FUNC_TV },
    { COLOR_PRIMARY_HDTV,     COLOR_MATRIX_HDTV, TRANSFER_FUNC_DCI_GAMMA },
};

static const enum PixelFormat dirac_pix_fmt[][2] = {
    { PIX_FMT_YUV444P, PIX_FMT_YUVJ444P },
    { PIX_FMT_YUV422P, PIX_FMT_YUVJ422P },
    { PIX_FMT_YUV420P, PIX_FMT_YUVJ420P },
};

/* Quarter pixel interpolation. */
static const uint8_t qpel_weights[][4] = {
    {  4,  0,  0,  0 }, /* rx=0, ry=0 */
    {  2,  0,  2,  0 }, /* rx=0, ry=1 */
    {  2,  2,  0,  0 }, /* rx=1, ry=0 */
    {  1,  1,  1,  1 }, /* rx=1, ry=1 */
};

static const uint8_t eighthpel_weights[][4] = {
    { 16,  0,  0,  0 }, /* rx=0, ry=0 */
    { 12,  0,  4,  0 }, /* rx=0, ry=1 */
    {  8,  0,  8,  0 }, /* rx=0, ry=2 */
    {  4,  0, 12,  0 }, /* rx=0, ry=3 */
    { 12,  4,  0,  0 }, /* rx=1, ry=0 */
    {  9,  3,  3,  1 }, /* rx=1, ry=1 */
    {  6,  2,  6,  2 }, /* rx=1, ry=2 */
    {  3,  1,  9,  3 }, /* rx=1, ry=3 */
    {  8,  8,  0,  0 }, /* rx=2, ry=0 */
    {  6,  6,  2,  2 }, /* rx=2, ry=1 */
    {  4,  4,  4,  4 }, /* rx=2, ry=2 */
    {  2,  2,  6,  6 }, /* rx=2, ry=3 */
    {  4, 12,  0,  0 }, /* rx=3, ry=0 */
    {  3,  9,  1,  3 }, /* rx=3, ry=1 */
    {  2,  6,  2,  6 }, /* rx=3, ry=2 */
    {  1,  3,  3,  9 }, /* rx=3, ry=3 */
};

uint8_t ff_dirac_default_qmat[][4][4] = {
    { { 5,  3,  3,  0}, { 0,  4,  4,  1}, { 0,  5,  5,  2}, { 0,  6,  6,  3} },
    { { 4,  2,  2,  0}, { 0,  4,  4,  2}, { 0,  5,  5,  3}, { 0,  7,  7,  5} },
    { { 5,  3,  3,  0}, { 0,  4,  4,  1}, { 0,  5,  5,  2}, { 0,  6,  6,  3} },
    { { 8,  4,  4,  0}, { 0,  4,  4,  0}, { 0,  4,  4,  0}, { 0,  4,  4,  0} },
    { { 8,  4,  4,  0}, { 0,  4,  4,  0}, { 0,  4,  4,  0}, { 0,  4,  4,  0} },
    { { 0,  4,  4,  8}, { 0,  8,  8, 12}, { 0, 13, 13, 17}, { 0, 17, 17, 21} },
    { { 3,  1,  1,  0}, { 0,  4,  4,  2}, { 0,  6,  6,  5}, { 0,  9,  9,  7} },
};

/**
 * Parse the source parameters in the sequence header.
 */
static int parse_source_parameters(GetBitContext *gb, AVCodecContext *avctx,
                                   dirac_source_params *source)
{
    AVRational frame_rate = (AVRational){0,0};
    unsigned luma_depth, chroma_depth;

    /* Override the luma dimensions. */
    if (get_bits1(gb)) {
        source->width  = svq3_get_ue_golomb(gb);
        source->height = svq3_get_ue_golomb(gb);
    }

    /* Override the chroma format. */
    if (get_bits1(gb))
        source->chroma_format = svq3_get_ue_golomb(gb);
    if (source->chroma_format > 2) {
        av_log(avctx, AV_LOG_ERROR, "Unknown chroma format %d\n",
               source->chroma_format);
        return -1;
    }

    if (get_bits1(gb))
        source->interlaced = svq3_get_ue_golomb(gb);
    if (source->interlaced > 1)
        return -1;

    /* framerate */
    if (get_bits1(gb)) {
        source->frame_rate_index = svq3_get_ue_golomb(gb);

        if (source->frame_rate_index > 10)
            return -1;

        if (!source->frame_rate_index) {
            frame_rate.num = svq3_get_ue_golomb(gb);
            frame_rate.den = svq3_get_ue_golomb(gb);
        }
    }
    if (source->frame_rate_index > 0) {
        if (source->frame_rate_index <= 8)
            frame_rate = ff_frame_rate_tab[source->frame_rate_index];
        else
            frame_rate = dirac_frame_rate[source->frame_rate_index-9];
    }
    av_reduce(&avctx->time_base.num, &avctx->time_base.den,
              frame_rate.den, frame_rate.num, 1<<30);

    /* Override aspect ratio. */
    if (get_bits1(gb)) {
        source->aspect_ratio_index = svq3_get_ue_golomb(gb);

        if (source->aspect_ratio_index > 6)
            return -1;

        if (!source->aspect_ratio_index) {
            avctx->sample_aspect_ratio.num = svq3_get_ue_golomb(gb);
            avctx->sample_aspect_ratio.den = svq3_get_ue_golomb(gb);
        }
    }
    if (source->aspect_ratio_index > 0)
        avctx->sample_aspect_ratio =
                dirac_preset_aspect_ratios[source->aspect_ratio_index-1];

    /* Override clean area. */
    if (get_bits1(gb)) {
        source->clean_width        = svq3_get_ue_golomb(gb);
        source->clean_height       = svq3_get_ue_golomb(gb);
        source->clean_left_offset  = svq3_get_ue_golomb(gb);
        source->clean_right_offset = svq3_get_ue_golomb(gb);
    }

    /* Override signal range. */
    if (get_bits1(gb)) {
        source->pixel_range_index = svq3_get_ue_golomb(gb);

        if (source->pixel_range_index > 4)
            return -1;

        if (!source->pixel_range_index) {
            source->pixel_range.luma_offset      = svq3_get_ue_golomb(gb);
            source->pixel_range.luma_excursion   = svq3_get_ue_golomb(gb);
            source->pixel_range.chroma_offset    = svq3_get_ue_golomb(gb);
            source->pixel_range.chroma_excursion = svq3_get_ue_golomb(gb);
        }
    }
    if (source->pixel_range_index > 0)
        source->pixel_range =
                dirac_pixel_range_presets[source->pixel_range_index-1];

    if (PIXEL_RANGE_EQUAL(source->pixel_range, dirac_pixel_range_presets[0]))
        // full range (JPEG YUV colorspace)
        avctx->pix_fmt = dirac_pix_fmt[source->chroma_format][1];
    else
        // normal YUV colorspace as the default
        avctx->pix_fmt = dirac_pix_fmt[source->chroma_format][0];

    /* color spec */
    source->color_spec = dirac_color_spec_presets[source->color_spec_index];
    if (get_bits1(gb)) {
        source->color_spec_index = svq3_get_ue_golomb(gb);

        if (source->color_spec_index > 4)
            return -1;

        source->color_spec= dirac_color_spec_presets[source->color_spec_index];

        if (!source->color_spec_index) {
            if (get_bits1(gb))
                source->color_spec.primaries = svq3_get_ue_golomb(gb);

            if (get_bits1(gb))
                source->color_spec.matrix = svq3_get_ue_golomb(gb);

            if (get_bits1(gb))
                source->color_spec.transfer_function = svq3_get_ue_golomb(gb);

            if (source->color_spec.primaries > 3 ||
                source->color_spec.matrix > 2 ||
                source->color_spec.transfer_function > 3)
                return -1;
        }
    }

    luma_depth   = av_log2(source->pixel_range.luma_excursion   + 1);
    chroma_depth = av_log2(source->pixel_range.chroma_excursion + 1);
    if (luma_depth > 8 || chroma_depth > 8)
        av_log(avctx, AV_LOG_WARNING, "Bitdepth greater than 8, may not work");

    return 0;
}

/**
 * Parse the sequence header.
 */
int ff_dirac_parse_sequence_header(GetBitContext *gb, AVCodecContext *avctx,
                                   dirac_source_params *source)
{
    unsigned int version_major;
    unsigned int version_minor;
    unsigned int video_format;
    unsigned int picture_coding_mode;

    /* parse parameters */
    version_major  = svq3_get_ue_golomb(gb);
    version_minor  = svq3_get_ue_golomb(gb);
    avctx->profile = svq3_get_ue_golomb(gb);
    avctx->level   = svq3_get_ue_golomb(gb);
    video_format   = svq3_get_ue_golomb(gb);

    if (version_major < 2)
        av_log(avctx, AV_LOG_WARNING, "Stream is old and may not work\n");
    else if (version_major > 2)
        av_log(avctx, AV_LOG_WARNING, "Stream may have unhandled features\n");

    if (video_format > 20)
        return -1;

    /* Fill in defaults for the source parameters. */
    *source = dirac_source_parameters_defaults[video_format];

    /* Override the defaults. */
    if (parse_source_parameters(gb, avctx, source))
        return -1;

    if (avcodec_check_dimensions(avctx, source->width, source->height))
        return -1;

    avcodec_set_dimensions(avctx, source->width, source->height);

    // coded as fields
    picture_coding_mode = svq3_get_ue_golomb(gb);
    if (picture_coding_mode != 0) {
        av_log(avctx, AV_LOG_ERROR, "Unsupported picture coding mode %d",
               picture_coding_mode);
        return -1;
    }

    return 0;
}

const struct dirac_block_params ff_dirac_block_param_defaults[] = {
    {  8,  8,  4,  4 },
    { 12, 12,  8,  8 },
    { 16, 16, 12, 12 },
    { 24, 24, 16, 16 }
};

/**
 * Interpolate a frame
 *
 * @param refframe frame to grab the upconverted pixel from
 * @param width    frame width
 * @param height   frame height
 * @param pixels   buffer to write the interpolated pixels to
 * @param comp     component
 */
static void interpolate_frame_halfpel(AVFrame *refframe, int width, int height,
                                      int8_t *pixels, int comp,
                                      int xpad, int ypad)
{
    int8_t *lineout;
    uint8_t *refdata;
    uint8_t *lineinref;
    int8_t *linein;
    int outwidth = width * 2 + xpad * 4;
    int doutwidth = 2 * outwidth;
    int x, y;
    const int t[4] = { 21, -7, 3, -1 };
    int8_t *pixelsdata = pixels + ypad * doutwidth + 2 * xpad;

    refdata    = refframe->data[comp];

    lineinref  = refdata;
    lineout = pixelsdata;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width; x++) {
            lineout[x * 2] = lineinref[x] - 128;
        }

        lineinref  += refframe->linesize[comp];
        lineout += doutwidth;
    }

    /* Copy top even lines. */
    linein  = pixels + ypad * doutwidth;
    lineout = pixels;
    for (y = 0; y < ypad * 2; y += 2) {
        memcpy(lineout, linein, outwidth);
        lineout += doutwidth;
    }

    /* Copy bottom even lines. */
    linein  = pixels + (ypad + height - 1) * doutwidth;
    lineout = linein + doutwidth;
    for (y = 0; y < ypad * 2; y += 2) {
        memcpy(lineout, linein, outwidth);
        lineout += doutwidth;
    }

    /* Interpolation (vertically). */
    linein  = pixelsdata;
    lineout = pixelsdata + outwidth;
    for (y = 0; y < height; y++) {
        for (x = 0; x < width * 2; x += 2) {
            int val = 16;
            int8_t *li1 = linein;
            int8_t *li2 = linein + doutwidth;

            val += t[0] * (li1[x] + li2[x]);
            li1 -= doutwidth;
            li2 += doutwidth;

            val += t[1] * (li1[x] + li2[x]);
            li1 -= doutwidth;
            li2 += doutwidth;

            val += t[2] * (li1[x] + li2[x]);
            li1 -= doutwidth;
            li2 += doutwidth;

            val += t[3] * (li1[x] + li2[x]);

            val >>= 5;

            lineout[x] = av_clip(val, -128, 127);
        }

        linein += doutwidth;

        /* Skip one line, we are interpolating to odd lines. */
        lineout    += doutwidth;
    }

    /* Add padding on the left and right sides of the frame. */
    lineout = pixels + 2 * ypad * outwidth;
    for (y = 0; y < height * 2; y++) {
        memset(lineout, lineout[2 * xpad], 2 * xpad);
        memset(&lineout[2 * width + xpad * 2],
               lineout[2 * width + xpad * 2 - 2], 2 * xpad);
        lineout += outwidth;
    }

    /* Interpolation (horizontally). */
    lineout = pixelsdata + 1;
    linein  = pixelsdata;
    for (y = 0; y < height * 2; y++) {
        for (x = 0; x < width * 2; x += 2) {
            int8_t *li1 = &linein[x];
            int8_t *li2 = &linein[x + 2];
            int val = 16;

            val += t[0] * (*li1 + *li2);
            li1 -= 2;
            li2 += 2;
            val += t[1] * (*li1 + *li2);
            li1 -= 2;
            li2 += 2;
            val += t[2] * (*li1 + *li2);
            li1 -= 2;
            li2 += 2;
            val += t[3] * (*li1 + *li2);

            val >>= 5;
            lineout[x] = av_clip(val, -128, 127);
        }
        lineout += outwidth;
        linein  += outwidth;
    }

    /* Add padding to right side of the frame. */
    lineout = pixels + 2 * ypad * outwidth;
    for (y = 0; y < height * 2; y++) {
        memset(&lineout[2 * width + xpad * 2],
               lineout[2 * width + xpad * 2 - 1], 2 * xpad);
        lineout += outwidth;
    }

    /* Copy top lines. */
    linein  = pixels + ypad * doutwidth;
    lineout = pixels;
    for (y = 0; y < ypad * 2; y++) {
        memcpy(lineout, linein, outwidth);
        lineout += outwidth;
    }

    /* Copy bottom lines. */
    linein  = pixels + (ypad + height - 1) * doutwidth;
    lineout = linein + outwidth;
    for (y = 0; y < ypad * 2; y++) {
        memcpy(lineout, linein, outwidth);
        lineout += outwidth;
    }
}

/**
 * Calculate WH or WV of the spatial weighting matrix
 *
 * @param i       block position
 * @param x       current pixel
 * @param bsep    block spacing
 * @param blen    block length
 * @param offset  xoffset/yoffset
 * @param blocks  number of blocks
 */
static inline
int spatial_wt(int i, int x, int bsep, int blen, int offset, int blocks)
{
    int pos = x - (i * bsep - offset);
    int max;

    max = 2 * (blen - bsep);
    if (i == 0 && pos < (blen >> 1))
        return max;
    else if (i == blocks - 1 && pos >= (blen >> 1))
        return max;
    else
        return av_clip(blen - FFABS(2*pos - (blen - 1)), 0, max);
}

/**
 * Motion Compensation with two reference frames
 *
 * @param coeffs     coefficients to add the DC to
 * @param i          horizontal position of the MC block
 * @param j          vertical position of the MC block
 * @param xstart     top left coordinate of the MC block
 * @param ystop      top left coordinate of the MC block
 * @param xstop      bottom right coordinate of the MC block
 * @param ystop      bottom right coordinate of the MC block
 * @param ref1       first reference frame
 * @param ref2       second reference frame
 * @param currblock  MC block to use
 * @param comp       component
 * @param border     0 if this is not a border MC block, otherwise 1
 */
static void motion_comp_block2refs(DiracContext *s, int16_t *coeffs,
                                   int i, int j, int xstart, int xstop,
                                   int ystart, int ystop, int8_t *ref1,
                                   int8_t *ref2,
                                   struct dirac_blockmotion *currblock,
                                   int comp, int border)
{
    int x, y;
    int xs, ys;
    int16_t *line;
    int8_t *refline1;
    int8_t *refline2;
    int vect1[2];
    int vect2[2];
    int refxstart1, refystart1;
    int refxstart2, refystart2;
    uint16_t *spatialwt;
    /* Subhalfpixel in qpel/eighthpel interpolated frame. */
    int rx1, ry1, rx2, ry2;
    const uint8_t *w1 = NULL;
    const uint8_t *w2 = NULL;
    int xfix1 = 0, xfix2 = 0;
    Plane *p = &s->plane[comp];

    vect1[0] = currblock->vect[0][0];
    vect1[1] = currblock->vect[0][1];
    vect2[0] = currblock->vect[1][0];
    vect2[1] = currblock->vect[1][1];

    xs = FFMAX(xstart, 0);
    ys = FFMAX(ystart, 0);

    if (comp != 0) {
        vect1[0] >>= s->chroma_hshift;
        vect2[0] >>= s->chroma_hshift;
        vect1[1] >>= s->chroma_vshift;
        vect2[1] >>= s->chroma_vshift;
    }

    switch(s->mv_precision) {
    case 0:
        refxstart1 = (xs + vect1[0]) << 1;
        refystart1 = (ys + vect1[1]) << 1;
        refxstart2 = (xs + vect2[0]) << 1;
        refystart2 = (ys + vect2[1]) << 1;
        break;
    case 1:
        refxstart1   = (xs << 1) + vect1[0];
        refystart1   = (ys << 1) + vect1[1];
        refxstart2   = (xs << 1) + vect2[0];
        refystart2   = (ys << 1) + vect2[1];
        break;
    case 2:
        refxstart1   = ((xs << 2) + vect1[0]) >> 1;
        refystart1   = ((ys << 2) + vect1[1]) >> 1;
        refxstart2   = ((xs << 2) + vect2[0]) >> 1;
        refystart2   = ((ys << 2) + vect2[1]) >> 1;
        rx1 = vect1[0] & 1;
        ry1 = vect1[1] & 1;
        rx2 = vect2[0] & 1;
        ry2 = vect2[1] & 1;
        w1 = qpel_weights[(rx1 << 1) | ry1];
        w2 = qpel_weights[(rx2 << 1) | ry2];
        break;
    case 3:
        refxstart1   = ((xs << 3) + vect1[0]) >> 2;
        refystart1   = ((ys << 3) + vect1[1]) >> 2;
        refxstart2   = ((xs << 3) + vect2[0]) >> 2;
        refystart2   = ((ys << 3) + vect2[1]) >> 2;
        rx1 = vect1[0] & 3;
        ry1 = vect1[1] & 3;
        rx2 = vect2[0] & 3;
        ry2 = vect2[0] & 3;
        w1 = eighthpel_weights[(rx1 << 2) | ry1];
        w2 = eighthpel_weights[(rx2 << 2) | ry2];
        break;
    default:
        /* XXX */
        return;
    }

    spatialwt = &s->spatialwt[p->xblen * (ys - ystart)];

    /* Make sure the vector doesn't point to a block outside the
       padded frame. */
    refystart1 = av_clip(refystart1, -p->yblen, p->height * 2 - 1);
    refystart2 = av_clip(refystart2, -p->yblen, p->height * 2 - 1);
    if (refxstart1 < -p->xblen)
        xfix1 = -p->xblen - refxstart1;
    else if (refxstart1 >= (p->width - 1) * 2)
        xfix1 = (p->width - 1) * 2 - refxstart1;
    if (refxstart2 < -p->xblen * 2)
        xfix2 = -p->xblen * 2 - refxstart2;
    else if (refxstart2 >= (p->width - 1) * 2)
        xfix2 = (p->width - 1) * 2 - refxstart2;

    line = &coeffs[p->width * ys];
    refline1 = &ref1[refystart1 * s->refwidth];
    refline2 = &ref2[refystart2 * s->refwidth];
    for (y = ys; y < ystop; y++) {
        int bx = xs - xstart;
        for (x = xs; x < xstop; x++) {
            int val1;
            int val2;
            int val;

            if (s->mv_precision == 0) {
                /* No interpolation. */
                val1 = refline1[(x + vect1[0]) << 1];
                val2 = refline2[(x + vect2[0]) << 1];
            } else if (s->mv_precision == 1) {
                /* Halfpel interpolation. */
                val1 = refline1[(x << 1) + vect1[0]];
                val2 = refline2[(x << 1) + vect2[0]];
            } else {
                /* Position in halfpel interpolated frame. */
                int hx1, hx2;

                if (s->mv_precision == 2) {
                    /* Do qpel interpolation. */
                    hx1 = ((x << 2) + vect1[0]) >> 1;
                    hx2 = ((x << 2) + vect2[0]) >> 1;
                    val1 = 2;
                    val2 = 2;
                } else {
                    /* Do eighthpel interpolation. */
                    hx1 = ((x << 3) + vect1[0]) >> 2;
                    hx2 = ((x << 3) + vect2[0]) >> 2;
                    val1 = 4;
                    val2 = 4;
                }

                /* Fix the x position on the halfpel interpolated
                   frame so it points to a MC block within the padded
                   region. */
                hx1 += xfix1;
                hx2 += xfix2;

                val1 += w1[0] * refline1[hx1                  ];
                val1 += w1[1] * refline1[hx1               + 1];
                val1 += w1[2] * refline1[hx1 + s->refwidth    ];
                val1 += w1[3] * refline1[hx1 + s->refwidth + 1];

                val2 += w2[0] * refline2[hx2                  ];
                val2 += w2[1] * refline2[hx2               + 1];
                val2 += w2[2] * refline2[hx2 + s->refwidth    ];
                val2 += w2[3] * refline2[hx2 + s->refwidth + 1];

                if (s->mv_precision) {
                    val1 >>= (s->mv_precision-1)<<1;
                    val2 >>= (s->mv_precision-1)<<1;
                }
            }

            val1 *= s->picture_weight_ref1;
            val2 *= s->picture_weight_ref2;
            val = val1 + val2;
            if (border) {
                val *= spatialwt[bx];
            } else {
                val = (val
                       * spatial_wt(i, x, p->xbsep, p->xblen,
                                    p->xoffset, p->current_blwidth)
                       * spatial_wt(j, y, p->ybsep, p->yblen,
                                    p->yoffset, p->current_blheight));
            }

            line[x] += val;
            bx++;
        }
        refline1 += s->refwidth << 1;
        refline2 += s->refwidth << 1;
        line += p->width;
        spatialwt += p->xblen;
    }
}

/**
 * Motion Compensation with one reference frame
 *
 * @param coeffs     coefficients to add the DC to
 * @param i          horizontal position of the MC block
 * @param j          vertical position of the MC block
 * @param xstart     top left coordinate of the MC block
 * @param ystop      top left coordinate of the MC block
 * @param xstop      bottom right coordinate of the MC block
 * @param ystop      bottom right coordinate of the MC block
 * @param refframe   reference frame
 * @param ref        0=first refframe 1=second refframe
 * @param currblock  MC block to use
 * @param comp       component
 * @param border     0 if this is not a border MC block, otherwise 1
 */
static void motion_comp_block1ref(DiracContext *s, int16_t *coeffs,
                                  int i, int j, int xstart, int xstop,
                                  int ystart, int ystop, int8_t *refframe,
                                  int ref,
                                  struct dirac_blockmotion *currblock,
                                  int comp, int border)
{
    int x, y;
    int xs, ys;
    int16_t *line;
    int8_t  *refline;
    int vect[2];
    int refxstart, refystart;
    uint16_t *spatialwt;
    /* Subhalfpixel in qpel/eighthpel interpolated frame. */
    int rx, ry;
    const uint8_t *w = NULL;
    int xfix = 0;
    Plane *p = &s->plane[comp];

    vect[0] = currblock->vect[ref][0];
    vect[1] = currblock->vect[ref][1];

    xs = FFMAX(xstart, 0);
    ys = FFMAX(ystart, 0);

    if (comp != 0) {
        vect[0] >>= s->chroma_hshift;
        vect[1] >>= s->chroma_vshift;
    }

    switch(s->mv_precision) {
    case 0:
        refxstart = (xs + vect[0]) << 1;
        refystart = (ys + vect[1]) << 1;
        break;
    case 1:
        refxstart   = (xs << 1) + vect[0];
        refystart   = (ys << 1) + vect[1];
        break;
    case 2:
        refxstart   = ((xs << 2) + vect[0]) >> 1;
        refystart   = ((ys << 2) + vect[1]) >> 1;
        rx = vect[0] & 1;
        ry = vect[1] & 1;
        w = qpel_weights[(rx << 1) | ry];
        break;
    case 3:
        refxstart   = ((xs << 3) + vect[0]) >> 2;
        refystart   = ((ys << 3) + vect[1]) >> 2;
        rx = vect[0] & 3;
        ry = vect[1] & 3;
        w = eighthpel_weights[(rx << 2) | ry];
        break;
    default:
        /* XXX */
        return;
    }

    /* Make sure the vector doesn't point to a block outside the
       padded frame. */
    refystart = av_clip(refystart, -p->yblen * 2, p->height * 2 - 1);
    if (refxstart < -p->xblen * 2)
        xfix = -p->xblen - refxstart;
    else if (refxstart >= (p->width - 1) * 2)
        xfix = (p->width - 1) * 2 - refxstart;

    spatialwt = &s->spatialwt[p->xblen * (ys - ystart)];

    line = &coeffs[p->width * ys];
    refline = &refframe[refystart * s->refwidth];
    for (y = ys; y < ystop; y++) {
        int bx = xs - xstart;
        for (x = xs; x < xstop; x++) {
            int val;

            if (s->mv_precision == 0) {
                /* No interpolation. */
                val = refline[(x + vect[0]) << 1];
            } else if (s->mv_precision == 1) {
                /* Halfpel interpolation. */
                val = refline[(x << 1) + vect[0]];
            } else {
                /* Position in halfpel interpolated frame. */
                int hx;

                if (s->mv_precision == 2) {
                    /* Do qpel interpolation. */
                    hx = ((x << 2) + vect[0]) >> 1;
                    val = 2;
                } else {
                    /* Do eighthpel interpolation. */
                    hx = ((x << 3) + vect[0]) >> 2;
                    val = 4;
                }

                /* Fix the x position on the halfpel interpolated
                   frame so it points to a MC block within the padded
                   region. */
                hx += xfix;

                val += w[0] * refline[hx                  ];
                val += w[1] * refline[hx               + 1];
                val += w[2] * refline[hx + s->refwidth    ];
                val += w[3] * refline[hx + s->refwidth + 1];
                if (s->mv_precision)
                    val >>= (s->mv_precision-1)<<1;
            }

            val *= s->picture_weight_ref1
                 + s->picture_weight_ref2;

            if (border) {
                val *= spatialwt[bx];
            } else {
                val = (val
                       * spatial_wt(i, x, p->xbsep, p->xblen,
                                    p->xoffset, p->current_blwidth)
                       * spatial_wt(j, y, p->ybsep, p->yblen,
                                    p->yoffset, p->current_blheight));
            }

            line[x] += val;
            bx++;
        }
        line += p->width;
        refline += s->refwidth << 1;
        spatialwt += p->xblen;
    }
}

/**
 * Motion Compensation DC values (no reference frame)
 *
 * @param coeffs coefficients to add the DC to
 * @param i      horizontal position of the MC block
 * @param j      vertical position of the MC block
 * @param xstart top left coordinate of the MC block
 * @param ystop  top left coordinate of the MC block
 * @param xstop  bottom right coordinate of the MC block
 * @param ystop  bottom right coordinate of the MC block
 * @param dcval  DC value to apply to all coefficients in the MC block
 * @param border 0 if this is not a border MC block, otherwise 1
 */
static inline
void motion_comp_dc_block(DiracContext *s, int16_t *coeffs, int i, int j,
                          int xstart, int xstop, int ystart, int ystop,
                          int dcval, int comp, int border)
{
    int x, y;
    int xs, ys;
    int16_t *line;
    int16_t *spatialwt;
    Plane *p = &s->plane[comp];

    ys = FFMAX(ystart, 0);
    xs = FFMAX(xstart, 0);

    dcval <<= s->picture_weight_precision;

    spatialwt = &s->spatialwt[p->xblen * (ys - ystart)];
    line = &coeffs[p->width * ys];
    for (y = ys; y < ystop; y++) {
        int bx = xs - xstart;
        for (x = xs; x < xstop; x++) {
            int val;

            if (border) {
                val = dcval * spatialwt[bx];
            } else {
                val = dcval
                    * spatial_wt(i, x, p->xbsep, p->xblen,
                                 p->xoffset, p->current_blwidth)
                    * spatial_wt(j, y, p->ybsep, p->yblen,
                                 p->yoffset, p->current_blheight);
            }

            line[x] += val;
            bx++;
        }
        line += p->width;
        spatialwt += p->xblen;
    }
}

/**
 * Motion compensation
 *
 * @param coeffs coefficients to which the MC results will be added
 * @param comp component
 * @return returns 0 on succes, otherwise -1
 */
int dirac_motion_compensation(DiracContext *s, int comp)
{
    int i, j;
    int x, y;
    struct dirac_blockmotion *currblock;
    int xstart, ystart;
    int xstop, ystop;
    Plane *p = &s->plane[comp];

    s->refwidth = (p->width + 2 * p->xblen) << 1;
    s->refheight = (p->height + 2 * p->yblen) << 1;

    s->spatialwt = av_malloc(p->xblen * p->yblen * sizeof(int16_t));
    if (!s->spatialwt) {
        av_log(s->avctx, AV_LOG_ERROR, "av_malloc() failed\n");
        return -1;
    }

    /* Set up the spatial weighting matrix. */
    for (x = 0; x < p->xblen; x++) {
        for (y = 0; y < p->yblen; y++) {
            int wh, wv;
            const int xmax = 2 * (p->xblen - p->xbsep);
            const int ymax = 2 * (p->yblen - p->ybsep);

            wh = av_clip(p->xblen - FFABS(2*x - (p->xblen - 1)), 0, xmax);
            wv = av_clip(p->yblen - FFABS(2*y - (p->yblen - 1)), 0, ymax);
            s->spatialwt[x + y * p->xblen] = wh * wv;
        }
    }

    if (avcodec_check_dimensions(s->avctx, s->refwidth, s->refheight))
        return -1;

    for (i = 0; i < s->refs; i++) {
        if (!s->ref_pics[i]) {
            av_log(s->avctx, AV_LOG_ERROR, "Reference frame %d not in buffer\n", i);
            return -1;
        }

            s->refdata[i] = av_malloc(s->refwidth * s->refheight);
            if (!s->refdata[i]) {
                if (i == 1)
                    av_free(s->refdata[0]);
                av_log(s->avctx, AV_LOG_ERROR, "av_malloc() failed\n");
                return -1;
            }
            interpolate_frame_halfpel(s->ref_pics[i], p->width, p->height,
                                      s->refdata[i], comp, p->xblen, p->yblen);
    }

    s->mcpic = av_malloc(p->width * p->height * sizeof(int16_t));
    if (!s->mcpic) {
        for (i = 0; i < s->refs; i++)
            av_free(s->refdata[i]);

        av_log(s->avctx, AV_LOG_ERROR, "av_malloc() failed\n");
        return -1;
    }
    memset(s->mcpic, 0, p->width * p->height * sizeof(int16_t));

    {
        currblock = s->blmotion;
        for (j = 0; j < p->current_blheight; j++) {
            for (i = 0; i < p->current_blwidth; i++) {
                struct dirac_blockmotion *block = &currblock[i];
                int border;
                int padding;

                /* XXX: These calculations do not match those in the
                   Dirac specification, but are correct. */
                xstart  = i * p->xbsep - p->xoffset;
                ystart  = j * p->ybsep - p->yoffset;
                xstop   = FFMIN(xstart + p->xblen, p->width);
                ystop   = FFMIN(ystart + p->yblen, p->height);

                border = (i > 0 && j > 0
                          && i < p->current_blwidth - 1
                          && j < p->current_blheight - 1);

                padding = 2 * ((p->xblen * 2 + p->width) * 2 * p->yblen
                               + p->xblen);

                /* Intra */
                if ((block->use_ref & 3) == 0)
                    motion_comp_dc_block(s, s->mcpic, i, j,
                                         xstart, xstop, ystart, ystop,
                                         block->dc[comp], comp, border);
                /* Reference frame 1 only. */
                else if ((block->use_ref & 3) == DIRAC_REF_MASK_REF1)
                    motion_comp_block1ref(s, s->mcpic, i, j,
                                          xstart, xstop, ystart,
                                          ystop,s->refdata[0] + padding,
                                          0, block, comp, border);
                /* Reference frame 2 only. */
                else if ((block->use_ref & 3) == DIRAC_REF_MASK_REF2)
                    motion_comp_block1ref(s, s->mcpic, i, j,
                                          xstart, xstop, ystart, ystop,
                                          s->refdata[1] + padding,
                                          1, block, comp, border);
                /* Both reference frames. */
                else
                    motion_comp_block2refs(s, s->mcpic, i, j,
                                           xstart, xstop, ystart, ystop,
                                           s->refdata[0] + padding,
                                           s->refdata[1] + padding,
                                           block, comp, border);
            }
            currblock += s->blwidth;
        }
    }

    av_freep(&s->spatialwt);

    for (i = 0; i < s->refs; i++) {
            av_freep(&s->refdata[i]);
    }

    return 0;
}
