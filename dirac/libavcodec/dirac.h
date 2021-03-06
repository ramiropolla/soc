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

#ifndef AVCODEC_DIRAC_H
#define AVCODEC_DIRAC_H

/**
 * @file
 * Interfaces to Dirac Decoder/Encoder
 * @author Marco Gerards <marco@gnu.org>
 */

#include "avcodec.h"
#include "get_bits.h"
#include "dirac_arith.h"
#include "dsputil.h"

typedef enum {
    COLOR_PRIMARY_HDTV,         ///< ITU-R BT. 709, also computer/web/sRGB
    COLOR_PRIMARY_SDTV_525,     ///< SMPTE 170M, 525 primaries
    COLOR_PRIMARY_SDTV_625,     ///< EBU Tech 3213-E, 625 primaries
    COLOR_PRIMARY_DCINEMA,      ///< SMPTE 428.1, CIE XYZ
} dirac_color_primary;

typedef enum {
    COLOR_MATRIX_HDTV,          ///< ITU-R BT.709, also computer/web
    COLOR_MATRIX_SDTV,          ///< ITU-R BT.601
    COLOR_MATRIX_REVERSIBLE,    ///< ITU-T H.264
} dirac_color_matrix;

typedef enum {
    TRANSFER_FUNC_TV,
    TRANSFER_FUNC_EXTENDED_GAMUT,
    TRANSFER_FUNC_LINEAR,
    TRANSFER_FUNC_DCI_GAMMA
} dirac_transfer_func;

typedef struct {
    dirac_color_primary primaries;
    dirac_color_matrix  matrix;
    dirac_transfer_func transfer_function;
} color_specification;

typedef struct {
    uint16_t luma_offset;
    uint16_t luma_excursion;
    uint16_t chroma_offset;
    uint16_t chroma_excursion;
} dirac_pixel_range;

#define PIXEL_RANGE_EQUAL(a,b) \
    (a.luma_offset == b.luma_offset && \
     a.luma_excursion == b.luma_excursion && \
     a.chroma_offset == b.chroma_offset && \
     a.chroma_excursion == b.chroma_excursion)

#define DIRAC_SIGN(x) ((x) > 0 ? ARITH_CONTEXT_SIGN_POS : \
                       (x) < 0 ? ARITH_CONTEXT_SIGN_NEG : \
                                 ARITH_CONTEXT_SIGN_ZERO)
#define DIRAC_PARSE_INFO_PREFIX 0x42424344
#define CALC_PADDING(size, depth) \
         (((size + (1 << depth) - 1) >> depth) << depth)


typedef struct {
    /* information about the frames */
    unsigned int width;                ///< the luma component width
    unsigned int height;               ///< the luma component height
    /** choma format: 0: 4:4:4, 1: 4:2:2, 2: 4:2:0 */
    uint8_t chroma_format;

    /* interlacing */
    uint8_t interlaced;                ///< flag for interlacing
    uint8_t top_field_first;

    uint8_t frame_rate_index;          ///< index into dirac_frame_rate[]
    uint8_t aspect_ratio_index;        ///< index into dirac_aspect_ratio[]

    /* clean area */
    uint16_t clean_width;
    uint16_t clean_height;
    uint16_t clean_left_offset;
    uint16_t clean_right_offset;

    uint8_t pixel_range_index;         ///< index into dirac_pixel_range_presets[]
    uint8_t color_spec_index;          ///< index into dirac_color_spec_presets[]

    dirac_pixel_range pixel_range;
    color_specification color_spec;
} dirac_source_params;

struct dirac_block_params {
    int xblen;
    int yblen;
    int xbsep;
    int ybsep;
};

struct globalmc_parameters {
    unsigned int pan_tilt[2];                   ///< pan/tilt vector
    unsigned int zrs[2][2];                     ///< zoom/rotate/shear matrix
    int perspective[2];                         ///< perspective vector
    unsigned int zrs_exp;
    unsigned int perspective_exp;
};

#define DIRAC_REF_MASK_REF1   1
#define DIRAC_REF_MASK_REF2   2
#define DIRAC_REF_MASK_GLOBAL 4

struct dirac_blockmotion {
    uint8_t use_ref;
    int16_t vect[2][2];
    int16_t dc[3];
};

// Schoedinger limits these to 8
#define MAX_REFERENCE_FRAMES 16
#define MAX_DELAYED_FRAMES 16
#define MAX_FRAMES 32
#define MAX_DECOMPOSITIONS 8

typedef struct SubBand{
    int level;
    int orientation;
    int stride;
    int width;
    int height;
    IDWTELEM *ibuf;
    struct SubBand *parent;
} SubBand;

typedef struct Plane{
    int width;
    int height;
    int padded_width;
    int padded_height;
    SubBand band[MAX_DECOMPOSITIONS][4];

    uint8_t xbsep;
    uint8_t xblen;
    uint8_t ybsep;
    uint8_t yblen;

    uint8_t xoffset;
    uint8_t yoffset;
    uint8_t total_wt_bits;
    uint8_t current_blwidth;
    uint8_t current_blheight;
} Plane;

typedef struct DiracContext {
    AVCodecContext *avctx;
    GetBitContext gb;
    struct dirac_arith_state arith;
    dirac_source_params source;
    Plane plane[3];
    int chroma_hshift;        ///< horizontal bits to shift for choma
    int chroma_vshift;        ///< vertical bits to shift for choma

    int zero_res;             ///< zero residue flag
    int is_arith;             ///< whether coeffs use arith or golomb coding
    int low_delay;            ///< use the low delay syntax
    int globalmc_flag;        ///< use global motion compensation flag
    int refs;                 ///< number of reference pictures

    // wavelet decoding
    uint8_t wavelet_depth;    ///< depth of the IDWT
    unsigned int wavelet_idx;
    unsigned int codeblock_mode;
    unsigned int codeblocksh[MAX_DECOMPOSITIONS+1];
    unsigned int codeblocksv[MAX_DECOMPOSITIONS+1];
    IDWTELEM *spatial_idwt_buffer;

    // low delay
    unsigned int x_slices;
    unsigned int y_slices;
    AVRational slice_bytes;
    uint8_t quant_matrix[MAX_DECOMPOSITIONS][4];

    // motion compensation
    uint8_t mv_precision;
    int16_t picture_weight_ref1;
    int16_t picture_weight_ref2;
    unsigned int picture_weight_precision;

    int blwidth;              ///< number of blocks (horizontally)
    int blheight;             ///< number of blocks (vertically)
    int sbwidth;              ///< number of superblocks (horizontally)
    int sbheight;             ///< number of superblocks (vertically)

    int *sbsplit;     // XXX: int8_t
    struct dirac_blockmotion *blmotion;
    struct globalmc_parameters globalmc;
    int16_t *mcpic;
    int16_t *spatialwt;
    int8_t *refdata[2];
    int refwidth;
    int refheight;


    AVFrame *current_picture;
    AVFrame *ref_pics[2];

    AVFrame *ref_frames[MAX_REFERENCE_FRAMES+1];
    AVFrame *delay_frames[MAX_DELAYED_FRAMES+1];
    AVFrame *all_frames;
} DiracContext;

typedef enum {
    pc_seq_header         = 0x00,
    pc_eos                = 0x10,
    pc_aux_data           = 0x20,
    pc_padding            = 0x60,
    pc_intra_ref          = 0x0c
} dirac_parse_code;

typedef enum {
    subband_ll = 0,
    subband_hl = 1,
    subband_lh = 2,
    subband_hh = 3
} dirac_subband;

// this assumes a max quantizer of 119 (larger would overflow 32 bits),
// which schoedinger and dirac-research also assume
static unsigned int inline coeff_quant_factor(unsigned int quant)
{
    uint64_t base = 1 << (quant >> 2);
    switch(quant & 3) {
    case 0:
        return base << 2;
    case 1:
        return (503829 * base + 52958) / 105917;
    case 2:
        return (665857 * base + 58854) / 117708;
    case 3:
        return (440253 * base + 32722) / 65444;
    }
    assert(0);
    return 0;
}

static unsigned int inline coeff_quant_offset(int is_intra, unsigned int quant)
{
    if (quant == 0)
        return 1;

    if (is_intra) {
        if (quant == 1)
            return 2;
        else
            return (coeff_quant_factor(quant) + 1) >> 1;
    }

    return (coeff_quant_factor(quant) * 3 + 4) >> 3;
}

static inline
int zero_neighbourhood(IDWTELEM *data, int x, int y, int stride)
{
    /* Check if there is a zero to the left and top left of this
       coefficient. */
    if (y > 0 && (data[-stride] || ( x > 0 && data[-stride - 1])))
        return 0;
    else if (x > 0 && data[- 1])
        return 0;

    return 1;
}

/**
 * Determine the most efficient context to use for arithmetic decoding
 * of this coefficient (given by a position in a subband).
 *
 * @param current coefficient
 * @param v vertical position of the coefficient
 * @param h horizontal position of the coefficient
 * @return prediction for the sign: -1 when negative, 1 when positive, 0 when 0
 */
static inline
int sign_predict(IDWTELEM *data, dirac_subband orientation,
                 int x, int y, int stride)
{
    if (orientation == subband_hl && y > 0)
        return DIRAC_SIGN(data[-stride]);
    else if (orientation == subband_lh && x > 0)
        return DIRAC_SIGN(data[-1]);
    else
        return ARITH_CONTEXT_SIGN_ZERO;
}

static inline
int intra_dc_coeff_prediction(IDWTELEM *coeff, int x, int y, int stride)
{
    int pred;
    if (x > 0 && y > 0) {
        pred = coeff[-1] + coeff[-stride] + coeff[-stride - 1];
        if (pred > 0)
            pred = (pred + 1) / 3;
        else
            pred = (pred - 1) / 3;
    } else if (x > 0) {
        /* Just use the coefficient left of this one. */
        pred = coeff[-1];
    } else if (y > 0)
        pred = coeff[-stride];
    else
        pred = 0;

    return pred;
}

static const int avgsplit[7] = { 0, 0, 1, 1, 1, 2, 2 };

static inline int split_prediction(DiracContext *s, int x, int y)
{
    if (x == 0 && y == 0)
        return 0;
    else if (y == 0)
        return s->sbsplit[ y      * s->sbwidth + x - 1];
    else if (x == 0)
        return s->sbsplit[(y - 1) * s->sbwidth + x    ];

    return avgsplit[s->sbsplit[(y - 1) * s->sbwidth + x    ]
                  + s->sbsplit[ y      * s->sbwidth + x - 1]
                  + s->sbsplit[(y - 1) * s->sbwidth + x - 1]];
}

/**
 * Mode prediction
 *
 * @param x    horizontal position of the MC block
 * @param y    vertical position of the MC block
 * @param ref reference frame
 */
static inline
int mode_prediction(DiracContext *s, int x, int y, int refmask, int refshift)
{
    int cnt;

    if (x == 0 && y == 0)
        return 0;
    else if (y == 0)
        return ((s->blmotion[ y      * s->blwidth + x - 1].use_ref & refmask)
                >> refshift);
    else if (x == 0)
        return ((s->blmotion[(y - 1) * s->blwidth + x    ].use_ref & refmask)
                >> refshift);

    /* Return the majority. */
    cnt = (s->blmotion[ y      * s->blwidth + x - 1].use_ref & refmask)
        + (s->blmotion[(y - 1) * s->blwidth + x    ].use_ref & refmask)
        + (s->blmotion[(y - 1) * s->blwidth + x - 1].use_ref & refmask);
    cnt >>= refshift;

    return cnt >> 1;
}

/**
 * Predict the motion vector
 *
 * @param x    horizontal position of the MC block
 * @param y    vertical position of the MC block
 * @param ref reference frame
 * @param dir direction horizontal=0, vertical=1
 */
static inline
int motion_vector_prediction(DiracContext *s, int x, int y, int ref, int dir)
{
    int cnt = 0;
    int left = 0, top = 0, lefttop = 0;
    const int refmask = ref + 1;
    const int mask = refmask | DIRAC_REF_MASK_GLOBAL;
    struct dirac_blockmotion *block = &s->blmotion[y * s->blwidth + x];

    if (x > 0) {
        /* Test if the block to the left has a motion vector for this
           reference frame. */
        if ((block[-1].use_ref & mask) == refmask) {
            left = block[-1].vect[ref][dir];
            cnt++;
        }

        /* This is the only reference, return it. */
        if (y == 0)
            return left;
    }

    if (y > 0) {
        /* Test if the block above the current one has a motion vector
           for this reference frame. */
        if ((block[-s->blwidth].use_ref & mask) == refmask) {
            top = block[-s->blwidth].vect[ref][dir];
            cnt++;
        }

        /* This is the only reference, return it. */
        if (x == 0)
            return top;
        else if (x > 0) {
            /* Test if the block above the current one has a motion vector
               for this reference frame. */
            if ((block[-s->blwidth - 1].use_ref & mask) == refmask) {
                lefttop = block[-s->blwidth - 1].vect[ref][dir];
                cnt++;
            }
        }
    }

    /* No references for the prediction. */
    if (cnt == 0)
        return 0;

    if (cnt == 1)
        return left + top + lefttop;

    /* Return the median of two motion vectors. */
    if (cnt == 2)
        return (left + top + lefttop + 1) >> 1;

    /* Return the median of three motion vectors. */
    return mid_pred(left, top, lefttop);
}

static inline
int block_dc_prediction(DiracContext *s, int x, int y, int comp)
{
    int total = 0;
    int cnt = 0;
    int sign;

    if (x > 0) {
        if (!(s->blmotion[y * s->blwidth + x - 1].use_ref & 3)) {
            total += s->blmotion[y * s->blwidth + x - 1].dc[comp];
            cnt++;
        }
    }

    if (y > 0) {
        if (!(s->blmotion[(y - 1) * s->blwidth + x].use_ref & 3)) {
            total += s->blmotion[(y - 1) * s->blwidth + x].dc[comp];
            cnt++;
        }
    }

    if (x > 0 && y > 0) {
        if (!(s->blmotion[(y - 1) * s->blwidth + x - 1].use_ref & 3)) {
            total += s->blmotion[(y - 1) * s->blwidth + x - 1].dc[comp];
            cnt++;
        }
    }

    if (cnt == 0)
        return 0;

    sign = FFSIGN(total);
    total = FFABS(total);

    /* Return the average of all DC values that were counted. */
    return sign * (total + (cnt >> 1)) / cnt;
}

int dirac_motion_compensation(DiracContext *s, int comp);

int dirac_decode_frame(AVCodecContext *avctx, void *data, int *data_size,
                       const uint8_t *buf, int buf_size);

int ff_dirac_parse_sequence_header(GetBitContext *gb, AVCodecContext *avctx,
                                   dirac_source_params *source);

extern const struct dirac_block_params ff_dirac_block_param_defaults[];
extern uint8_t ff_dirac_default_qmat[][4][4];

#endif /* AVCODEC_DIRAC_H */
