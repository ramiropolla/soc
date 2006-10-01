/*
 * copyright (c) 2006 Oded Shimon <ods15@ods15.dyndns.org>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/**
 * @file vorbis_enc.c
 * Native Vorbis encoder.
 * @author Oded Shimon <ods15@ods15.dyndns.org>
 */

#include "avcodec.h"
#include "dsputil.h"
#include "vorbis.h"

#undef NDEBUG
#include <assert.h>

typedef struct {
    int nentries;
    uint8_t * lens;
    uint32_t * codewords;
    int ndimentions;
    float min;
    float delta;
    int seq_p;
    int lookup;
    int * quantlist;
    float * dimentions;
} codebook_t;

typedef struct {
    int dim;
    int subclass;
    int masterbook;
    int * books;
} floor_class_t;

typedef struct {
    int partitions;
    int * partition_to_class;
    int nclasses;
    floor_class_t * classes;
    int multiplier;
    int rangebits;
    int values;
    floor1_entry_t * list;
} floor_t;

typedef struct {
    int type;
    int begin;
    int end;
    int partition_size;
    int classifications;
    int classbook;
    int8_t (*books)[8];
    float (*maxes)[2];
} residue_t;

typedef struct {
    int submaps;
    int * mux;
    int * floor;
    int * residue;
    int coupling_steps;
    int * magnitude;
    int * angle;
} mapping_t;

typedef struct {
    int blockflag;
    int mapping;
} vorbis_mode_t;

typedef struct {
    int channels;
    int sample_rate;
    int blocksize[2]; // in (1<<n) format
    MDCTContext mdct[2];
    const float * win[2];
    int have_saved;
    float * saved;
    float * samples;
    float * floor; // also used for tmp values for mdct
    float * coeffs; // also used for residue after floor
    float quality;

    int ncodebooks;
    codebook_t * codebooks;

    int nfloors;
    floor_t * floors;

    int nresidues;
    residue_t * residues;

    int nmappings;
    mapping_t * mappings;

    int nmodes;
    vorbis_mode_t * modes;
} venc_context_t;

typedef struct {
    int total;
    int total_pos;
    int pos;
    uint8_t * buf_ptr;
} PutBitContext;

static const uint8_t codebook0[] = {
   2, 10,  8, 14,  7, 12, 11, 14,  1,  5,  3,  7,  4,  9,  7,
  13,
};

static const uint8_t codebook1[] = {
   1,  4,  2,  6,  3,  7,  5,  7,
};

static const uint8_t codebook2[] = {
   1,  5,  7, 21,  5,  8,  9, 21, 10,  9, 12, 20, 20, 16, 20,
  20,  4,  8,  9, 20,  6,  8,  9, 20, 11, 11, 13, 20, 20, 15,
  17, 20,  9, 11, 14, 20,  8, 10, 15, 20, 11, 13, 15, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 13, 20, 20, 20, 18, 18, 20, 20,
  20, 20, 20, 20,  3,  6,  8, 20,  6,  7,  9, 20, 10,  9, 12,
  20, 20, 20, 20, 20,  5,  7,  9, 20,  6,  6,  9, 20, 10,  9,
  12, 20, 20, 20, 20, 20,  8, 10, 13, 20,  8,  9, 12, 20, 11,
  10, 12, 20, 20, 20, 20, 20, 18, 20, 20, 20, 15, 17, 18, 20,
  18, 17, 18, 20, 20, 20, 20, 20,  7, 10, 12, 20,  8,  9, 11,
  20, 14, 13, 14, 20, 20, 20, 20, 20,  6,  9, 12, 20,  7,  8,
  11, 20, 12, 11, 13, 20, 20, 20, 20, 20,  9, 11, 15, 20,  8,
  10, 14, 20, 12, 11, 14, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 11, 16, 18,
  20, 15, 15, 17, 20, 20, 17, 20, 20, 20, 20, 20, 20,  9, 14,
  16, 20, 12, 12, 15, 20, 17, 15, 18, 20, 20, 20, 20, 20, 16,
  19, 18, 20, 15, 16, 20, 20, 17, 17, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20,
};

static const uint8_t codebook3[] = {
   2,  3,  7, 13,  4,  4,  7, 15,  8,  6,  9, 17, 21, 16, 15,
  21,  2,  5,  7, 11,  5,  5,  7, 14,  9,  7, 10, 16, 17, 15,
  16, 21,  4,  7, 10, 17,  7,  7,  9, 15, 11,  9, 11, 16, 21,
  18, 15, 21, 18, 21, 21, 21, 15, 17, 17, 19, 21, 19, 18, 20,
  21, 21, 21, 20,
};

static const uint8_t codebook4[] = {
   5,  5,  5,  5,  6,  5,  6,  5,  6,  5,  6,  5,  6,  5,  6,
   5,  6,  5,  6,  5,  6,  5,  6,  5,  7,  5,  7,  5,  7,  5,
   7,  5,  8,  6,  8,  6,  8,  6,  9,  6,  9,  6, 10,  6, 10,
   6, 11,  6, 11,  7, 11,  7, 12,  7, 12,  7, 12,  7, 12,  7,
  12,  7, 12,  7, 12,  7, 12,  8, 13,  8, 12,  8, 12,  8, 13,
   8, 13,  9, 13,  9, 13,  9, 13,  9, 12, 10, 12, 10, 13, 10,
  14, 11, 14, 12, 14, 13, 14, 13, 14, 14, 15, 16, 15, 15, 15,
  14, 15, 17, 21, 22, 22, 21, 22, 22, 22, 22, 22, 22, 21, 21,
  21, 21, 21, 21, 21, 21, 21, 21,
};

static const uint8_t codebook5[] = {
   2,  5,  5,  4,  5,  4,  5,  4,  5,  4,  6,  5,  6,  5,  6,
   5,  6,  5,  7,  5,  7,  6,  8,  6,  8,  6,  8,  6,  9,  6,
   9,  6,
};

static const uint8_t codebook6[] = {
   8,  5,  8,  4,  9,  4,  9,  4,  9,  4,  9,  4,  9,  4,  9,
   4,  9,  4,  9,  4,  9,  4,  8,  4,  8,  4,  9,  5,  9,  5,
   9,  5,  9,  5,  9,  6, 10,  6, 10,  7, 10,  8, 11,  9, 11,
  11, 12, 13, 12, 14, 13, 15, 13, 15, 14, 16, 14, 17, 15, 17,
  15, 15, 16, 16, 15, 16, 16, 16, 15, 18, 16, 15, 17, 17, 19,
  19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
  19, 19, 19, 19, 19, 19,
};

static const uint8_t codebook7[] = {
   1,  5,  5,  5,  5,  5,  5,  5,  6,  5,  6,  5,  6,  5,  6,
   5,  6,  6,  7,  7,  7,  7,  8,  7,  8,  8,  9,  8, 10,  9,
  10,  9,
};

static const uint8_t codebook8[] = {
   4,  3,  4,  3,  4,  4,  5,  4,  5,  4,  5,  5,  6,  5,  6,
   5,  7,  5,  7,  6,  7,  6,  8,  7,  8,  7,  8,  7,  9,  8,
   9,  9,  9,  9, 10, 10, 10, 11,  9, 12,  9, 12,  9, 15, 10,
  14,  9, 13, 10, 13, 10, 12, 10, 12, 10, 13, 10, 12, 11, 13,
  11, 14, 12, 13, 13, 14, 14, 13, 14, 15, 14, 16, 13, 13, 14,
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16,
  16, 16, 16, 16, 15, 15,
};

static const uint8_t codebook9[] = {
   4,  5,  4,  5,  3,  5,  3,  5,  3,  5,  4,  4,  4,  4,  5,
   5,  5,
};

static const uint8_t codebook10[] = {
   3,  3,  4,  3,  4,  4,  4,  4,  5,  5,  5,  5,  5,  6,  5,
   7,  5,  8,  6,  8,  6,  9,  7, 10,  7, 10,  8, 10,  8, 11,
   9, 11,
};

static const uint8_t codebook11[] = {
   3,  7,  3,  8,  3, 10,  3,  8,  3,  9,  3,  8,  4,  9,  4,
   9,  5,  9,  6, 10,  6,  9,  7, 11,  7, 12,  9, 13, 10, 13,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12,
};

static const uint8_t codebook12[] = {
   4,  5,  4,  5,  4,  5,  4,  5,  3,  5,  3,  5,  3,  5,  4,
   5,  4,
};

static const uint8_t codebook13[] = {
   4,  2,  4,  2,  5,  3,  5,  4,  6,  6,  6,  7,  7,  8,  7,
   8,  7,  8,  7,  9,  8,  9,  8,  9,  8, 10,  8, 11,  9, 12,
   9, 12,
};

static const uint8_t codebook14[] = {
   2,  5,  2,  6,  3,  6,  4,  7,  4,  7,  5,  9,  5, 11,  6,
  11,  6, 11,  7, 11,  6, 11,  6, 11,  9, 11,  8, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 10, 10, 10,
  10, 10, 10,
};

static const uint8_t codebook15[] = {
   5,  6, 11, 11, 11, 11, 10, 10, 12, 11,  5,  2, 11,  5,  6,
   6,  7,  9, 11, 13, 13, 10,  7, 11,  6,  7,  8,  9, 10, 12,
  11,  5, 11,  6,  8,  7,  9, 11, 14, 15, 11,  6,  6,  8,  4,
   5,  7,  8, 10, 13, 10,  5,  7,  7,  5,  5,  6,  8, 10, 11,
  10,  7,  7,  8,  6,  5,  5,  7,  9,  9, 11,  8,  8, 11,  8,
   7,  6,  6,  7,  9, 12, 11, 10, 13,  9,  9,  7,  7,  7,  9,
  11, 13, 12, 15, 12, 11,  9,  8,  8,  8,
};

static const uint8_t codebook16[] = {
   2,  4,  4,  0,  0,  0,  0,  0,  0,  5,  6,  6,  0,  0,  0,
   0,  0,  0,  5,  6,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  5,  7,  7,  0,  0,  0,  0,  0,  0,
   7,  8,  8,  0,  0,  0,  0,  0,  0,  6,  7,  8,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  7,  7,
   0,  0,  0,  0,  0,  0,  6,  8,  7,  0,  0,  0,  0,  0,  0,
   7,  8,  8,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  7,  7,  0,  0,  0,
   0,  0,  0,  7,  8,  8,  0,  0,  0,  0,  0,  0,  7,  8,  8,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   7,  8,  8,  0,  0,  0,  0,  0,  0,  8,  8,  9,  0,  0,  0,
   0,  0,  0,  8,  9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  6,  8,  8,  0,  0,  0,  0,  0,  0,
   7,  9,  8,  0,  0,  0,  0,  0,  0,  8,  9,  9,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  5,  7,  7,  0,  0,  0,  0,  0,  0,  7,  8,  8,
   0,  0,  0,  0,  0,  0,  7,  8,  8,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  8,  8,  0,  0,  0,
   0,  0,  0,  8,  9,  9,  0,  0,  0,  0,  0,  0,  7,  8,  9,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   6,  8,  8,  0,  0,  0,  0,  0,  0,  8,  9,  9,  0,  0,  0,
   0,  0,  0,  8,  9,  8,
};

static const uint8_t codebook17[] = {
   2,  5,  5,  0,  0,  0,  5,  5,  0,  0,  0,  5,  5,  0,  0,
   0,  7,  8,  0,  0,  0,  0,  0,  0,  0,  5,  6,  6,  0,  0,
   0,  7,  7,  0,  0,  0,  7,  7,  0,  0,  0, 10, 10,  0,  0,
   0,  0,  0,  0,  0,  5,  6,  6,  0,  0,  0,  7,  7,  0,  0,
   0,  7,  7,  0,  0,  0, 10, 10,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   5,  7,  7,  0,  0,  0,  7,  7,  0,  0,  0,  7,  7,  0,  0,
   0,  9,  9,  0,  0,  0,  0,  0,  0,  0,  5,  7,  7,  0,  0,
   0,  7,  7,  0,  0,  0,  7,  7,  0,  0,  0,  9,  9,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  5,  7,  7,  0,  0,  0,  7,  7,  0,  0,
   0,  7,  7,  0,  0,  0,  9,  9,  0,  0,  0,  0,  0,  0,  0,
   5,  7,  7,  0,  0,  0,  7,  7,  0,  0,  0,  7,  7,  0,  0,
   0,  9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  8, 10, 10,  0,  0,
   0,  9,  9,  0,  0,  0,  9,  9,  0,  0,  0, 10, 10,  0,  0,
   0,  0,  0,  0,  0,  8, 10, 10,  0,  0,  0,  9,  9,  0,  0,
   0,  9,  9,  0,  0,  0, 10, 10,
};

static const uint8_t codebook18[] = {
   2,  4,  3,  6,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  4,  4,  6,  6,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  4,  4,  4,  6,  6,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   6,  6,  6,  9,  9,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  6,  6,  7,  9,  9,
};

static const uint8_t codebook19[] = {
   2,  3,  3,  6,  6,  0,  0,  0,  0,  0,  4,  4,  6,  6,  0,
   0,  0,  0,  0,  4,  4,  6,  6,  0,  0,  0,  0,  0,  5,  5,
   6,  6,  0,  0,  0,  0,  0,  0,  0,  6,  6,  0,  0,  0,  0,
   0,  0,  0,  7,  8,  0,  0,  0,  0,  0,  0,  0,  7,  7,  0,
   0,  0,  0,  0,  0,  0,  9,  9,
};

static const uint8_t codebook20[] = {
   1,  3,  4,  6,  6,  7,  7,  9,  9,  0,  5,  5,  7,  7,  7,
   8,  9,  9,  0,  5,  5,  7,  7,  8,  8,  9,  9,  0,  7,  7,
   8,  8,  8,  8, 10, 10,  0,  0,  0,  8,  8,  8,  8, 10, 10,
   0,  0,  0,  9,  9,  9,  9, 10, 10,  0,  0,  0,  9,  9,  9,
   9, 10, 10,  0,  0,  0, 10, 10, 10, 10, 11, 11,  0,  0,  0,
   0,  0, 10, 10, 11, 11,
};

static const uint8_t codebook21[] = {
   2,  3,  3,  6,  6,  7,  7,  8,  8,  8,  8,  9,  9, 10, 10,
  11, 10,  0,  5,  5,  7,  7,  8,  8,  9,  9,  9,  9, 10, 10,
  10, 10, 11, 11,  0,  5,  5,  7,  7,  8,  8,  9,  9,  9,  9,
  10, 10, 10, 10, 11, 11,  0,  6,  6,  7,  7,  8,  8,  9,  9,
   9,  9, 10, 10, 11, 11, 11, 11,  0,  0,  0,  7,  7,  8,  8,
   9,  9,  9,  9, 10, 10, 11, 11, 11, 12,  0,  0,  0,  8,  8,
   8,  8,  9,  9,  9,  9, 10, 10, 11, 11, 12, 12,  0,  0,  0,
   8,  8,  8,  8,  9,  9,  9,  9, 10, 10, 11, 11, 12, 12,  0,
   0,  0,  9,  9,  9,  9, 10, 10, 10, 10, 11, 10, 11, 11, 12,
  12,  0,  0,  0,  0,  0,  9,  9, 10, 10, 10, 10, 11, 11, 11,
  11, 12, 12,  0,  0,  0,  0,  0,  9,  8,  9,  9, 10, 10, 11,
  11, 12, 12, 12, 12,  0,  0,  0,  0,  0,  8,  8,  9,  9, 10,
  10, 11, 11, 12, 11, 12, 12,  0,  0,  0,  0,  0,  9, 10, 10,
  10, 11, 11, 11, 11, 12, 12, 13, 13,  0,  0,  0,  0,  0,  0,
   0, 10, 10, 10, 10, 11, 11, 12, 12, 13, 13,  0,  0,  0,  0,
   0,  0,  0, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13,  0,  0,
   0,  0,  0,  0,  0, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13,
   0,  0,  0,  0,  0,  0,  0, 11, 11, 12, 12, 12, 12, 13, 13,
  13, 13,  0,  0,  0,  0,  0,  0,  0,  0,  0, 12, 12, 12, 12,
  13, 13, 13, 13,
};

static const uint8_t codebook22[] = {
   1,  4,  4,  7,  6,  6,  7,  6,  6,  4,  7,  7, 10,  9,  9,
  11,  9,  9,  4,  7,  7, 10,  9,  9, 11,  9,  9,  7, 10, 10,
  11, 11, 10, 12, 11, 11,  6,  9,  9, 11, 10, 10, 11, 10, 10,
   6,  9,  9, 11, 10, 10, 11, 10, 10,  7, 11, 11, 11, 11, 11,
  12, 11, 11,  6,  9,  9, 11, 10, 10, 11, 10, 10,  6,  9,  9,
  11, 10, 10, 11, 10, 10,
};

static const uint8_t codebook23[] = {
   2,  4,  4,  6,  6,  7,  7,  7,  7,  8,  8, 10,  5,  5,  6,
   6,  7,  7,  8,  8,  8,  8, 10,  5,  5,  6,  6,  7,  7,  8,
   8,  8,  8, 10,  6,  6,  7,  7,  8,  8,  8,  8,  8,  8, 10,
  10, 10,  7,  7,  8,  7,  8,  8,  8,  8, 10, 10, 10,  8,  8,
   8,  8,  8,  8,  8,  8, 10, 10, 10,  7,  8,  8,  8,  8,  8,
   8,  8, 10, 10, 10,  8,  8,  8,  8,  8,  8,  8,  8, 10, 10,
  10, 10, 10,  8,  8,  8,  8,  8,  8, 10, 10, 10, 10, 10,  9,
   9,  8,  8,  9,  8, 10, 10, 10, 10, 10,  8,  8,  8,  8,  8,
   8,
};

static const uint8_t codebook24[] = {
   1,  4,  4,  6,  6,  7,  7,  8,  8,  9,  9, 10, 10,  6,  5,
   5,  7,  7,  8,  8,  8,  8,  9,  9, 10, 10,  7,  5,  5,  7,
   7,  8,  8,  8,  8,  9,  9, 11, 10,  0,  8,  8,  8,  8,  9,
   9,  9,  9, 10, 10, 11, 11,  0,  8,  8,  8,  8,  9,  9,  9,
   9, 10, 10, 11, 11,  0, 12, 12,  9,  9, 10, 10, 10, 10, 11,
  11, 11, 12,  0, 13, 13,  9,  9, 10, 10, 10, 10, 11, 11, 12,
  12,  0,  0,  0, 10, 10, 10, 10, 11, 11, 12, 12, 12, 12,  0,
   0,  0, 10, 10, 10, 10, 11, 11, 12, 12, 12, 12,  0,  0,  0,
  14, 14, 11, 11, 11, 11, 12, 12, 13, 13,  0,  0,  0, 14, 14,
  11, 11, 11, 11, 12, 12, 13, 13,  0,  0,  0,  0,  0, 12, 12,
  12, 12, 13, 13, 14, 13,  0,  0,  0,  0,  0, 13, 13, 12, 12,
  13, 12, 14, 13,
};

static const uint8_t codebook25[] = {
   2,  4,  4,  5,  5,  6,  5,  5,  5,  5,  6,  4,  5,  5,  5,
   6,  5,  5,  5,  5,  6,  6,  6,  5,  5,
};

static const uint8_t codebook26[] = {
   1,  4,  4, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,  4,  9,
   8, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,  2,  9,  7, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 11, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12,
  12, 12, 12, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
  11, 11, 11, 11,
};

static const uint8_t codebook27[] = {
   1,  4,  4,  6,  6,  7,  7,  8,  7,  9,  9, 10, 10, 10, 10,
   6,  5,  5,  7,  7,  8,  8, 10,  8, 11, 10, 12, 12, 13, 13,
   6,  5,  5,  7,  7,  8,  8, 10,  9, 11, 11, 12, 12, 13, 12,
  18,  8,  8,  8,  8,  9,  9, 10,  9, 11, 10, 12, 12, 13, 13,
  18,  8,  8,  8,  8,  9,  9, 10, 10, 11, 11, 13, 12, 14, 13,
  18, 11, 11,  9,  9, 10, 10, 11, 11, 11, 12, 13, 12, 13, 14,
  18, 11, 11,  9,  8, 11, 10, 11, 11, 11, 11, 12, 12, 14, 13,
  18, 18, 18, 10, 11, 10, 11, 12, 12, 12, 12, 13, 12, 14, 13,
  18, 18, 18, 10, 11, 11,  9, 12, 11, 12, 12, 12, 13, 13, 13,
  18, 18, 17, 14, 14, 11, 11, 12, 12, 13, 12, 14, 12, 14, 13,
  18, 18, 18, 14, 14, 11, 10, 12,  9, 12, 13, 13, 13, 13, 13,
  18, 18, 17, 16, 18, 13, 13, 12, 12, 13, 11, 14, 12, 14, 14,
  17, 18, 18, 17, 18, 13, 12, 13, 10, 12, 11, 14, 14, 14, 14,
  17, 18, 18, 18, 18, 15, 16, 12, 12, 13, 10, 14, 12, 14, 15,
  18, 18, 18, 16, 17, 16, 14, 12, 11, 13, 10, 13, 13, 14, 15,
};

static const uint8_t codebook28[] = {
   2,  5,  5,  6,  6,  7,  7,  7,  7,  7,  7,  8,  8,  8,  8,
   8,  8, 10,  6,  6,  7,  7,  8,  7,  8,  8,  8,  8,  8,  9,
   9,  9,  9,  9, 10,  6,  6,  7,  7,  7,  7,  8,  8,  8,  8,
   9,  9,  9,  9,  9,  9, 10,  7,  7,  7,  7,  8,  8,  8,  8,
   9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10,  7,  7,  8,  8,
   8,  9,  9,  9,  9,  9,  9,  9,  9,  9, 11, 11, 11,  8,  8,
   8,  8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10, 10, 10,
   8,  8,  8,  8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10,
  10, 10,  8,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9,  9, 10,
   9, 10, 10, 10, 11, 11,  9,  9,  9,  9,  9,  9,  9,  9,  9,
   9,  9,  9, 11, 10, 11, 11, 11,  9,  9,  9,  9,  9,  9, 10,
  10,  9,  9, 10,  9, 11, 10, 11, 11, 11,  9,  9,  9,  9,  9,
   9,  9,  9, 10, 10, 10,  9, 11, 11, 11, 11, 11,  9,  9,  9,
   9, 10, 10,  9,  9,  9,  9, 10,  9, 11, 11, 11, 11, 11, 11,
  11,  9,  9,  9,  9,  9,  9, 10, 10, 10, 10, 11, 11, 11, 11,
  11, 11, 11, 10,  9, 10, 10,  9, 10,  9,  9, 10,  9, 11, 10,
  10, 11, 11, 11, 11,  9, 10,  9,  9,  9,  9, 10, 10, 10, 10,
  11, 11, 11, 11, 11, 11, 10, 10, 10,  9,  9, 10,  9, 10,  9,
  10, 10, 10, 10, 11, 11, 11, 11, 11, 11, 11,  9,  9,  9,  9,
   9, 10, 10, 10,
};

static const struct {
    int dim;
    int len;
    int real_len;
    const uint8_t * clens;
    int lookup;
    float min;
    float delta;
    const uint8_t * quant;
} cvectors[] = {
    { 2,   16,   16, codebook0,  0 },
    { 2,    8,    8, codebook1,  0 },
    { 2,  256,  256, codebook2,  0 },
    { 2,   64,   64, codebook3,  0 },
    { 2,  128,  128, codebook4,  0 },
    { 2,   32,   32, codebook5,  0 },
    { 2,   96,   96, codebook6,  0 },
    { 2,   32,   32, codebook7,  0 },
    { 2,   96,   96, codebook8,  0 },
    { 2,   17,   17, codebook9,  0 },
    { 2,   32,   32, codebook10, 0 },
    { 2,   78,   78, codebook11, 0 },
    { 2,   17,   17, codebook12, 0 },
    { 2,   32,   32, codebook13, 0 },
    { 2,   78,   78, codebook14, 0 },
    { 2,  100,  100, codebook15, 0 },
    { 8, 1641, 6561, codebook16, 1,    -1.0,   1.0, (const uint8_t[]){ 1, 0, 2, } },
    { 4,  443,  625, codebook17, 1,    -2.0,   1.0, (const uint8_t[]){ 2, 1, 3, 0, 4, } },
    { 4,  105,  625, codebook18, 1,    -2.0,   1.0, (const uint8_t[]){ 2, 1, 3, 0, 4, } },
    { 2,   68,   81, codebook19, 1,    -4.0,   1.0, (const uint8_t[]){ 4, 3, 5, 2, 6, 1, 7, 0, 8, } },
    { 2,   81,   81, codebook20, 1,    -4.0,   1.0, (const uint8_t[]){ 4, 3, 5, 2, 6, 1, 7, 0, 8, } },
    { 2,  289,  289, codebook21, 1,    -8.0,   1.0, (const uint8_t[]){ 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15, 0, 16, } },
    { 4,   81,   81, codebook22, 1,   -11.0,  11.0, (const uint8_t[]){ 1, 0, 2, } },
    { 2,  121,  121, codebook23, 1,    -5.0,   1.0, (const uint8_t[]){ 5, 4, 6, 3, 7, 2, 8, 1, 9, 0, 10, } },
    { 2,  169,  169, codebook24, 1,   -30.0,   5.0, (const uint8_t[]){ 6, 5, 7, 4, 8, 3, 9, 2, 10, 1, 11, 0, 12, } },
    { 2,   25,   25, codebook25, 1,    -2.0,   1.0, (const uint8_t[]){ 2, 1, 3, 0, 4, } },
    { 2,  169,  169, codebook26, 1, -1530.0, 255.0, (const uint8_t[]){ 6, 5, 7, 4, 8, 3, 9, 2, 10, 1, 11, 0, 12, } },
    { 2,  225,  225, codebook27, 1,  -119.0,  17.0, (const uint8_t[]){ 7, 6, 8, 5, 9, 4, 10, 3, 11, 2, 12, 1, 13, 0, 14, } },
    { 2,  289,  289, codebook28, 1,    -8.0,   1.0, (const uint8_t[]){ 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15, 0, 16, } },
};

static inline void init_put_bits(PutBitContext * pb, uint8_t * buf, int buffer_len) {
    pb->total = buffer_len * 8;
    pb->total_pos = 0;
    pb->pos = 0;
    pb->buf_ptr = buf;
}

static void put_bits(PutBitContext * pb, int bits, uint64_t val) {
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

static inline void flush_put_bits(PutBitContext * pb) {
}

static inline int put_bits_count(PutBitContext * pb) {
    return pb->total_pos;
}

static inline void put_codeword(PutBitContext * pb, codebook_t * cb, int entry) {
    assert(entry >= 0);
    assert(entry < cb->nentries);
    assert(cb->lens[entry]);
    put_bits(pb, cb->lens[entry], cb->codewords[entry]);
}

static int cb_lookup_vals(int lookup, int dimentions, int entries) {
    if (lookup == 1) return ff_vorbis_nth_root(entries, dimentions);
    else if (lookup == 2) return dimentions * entries;
    return 0;
}

static void ready_codebook(codebook_t * cb) {
    int i;

    ff_vorbis_len2vlc(cb->lens, cb->codewords, cb->nentries);

    if (!cb->lookup) cb->dimentions = NULL;
    else {
        int vals = cb_lookup_vals(cb->lookup, cb->ndimentions, cb->nentries);
        cb->dimentions = av_malloc(sizeof(float) * cb->nentries * cb->ndimentions);
        for (i = 0; i < cb->nentries; i++) {
            float last = 0;
            int j;
            int div = 1;
            for (j = 0; j < cb->ndimentions; j++) {
                int off;
                if (cb->lookup == 1) off = (i / div) % vals; // lookup type 1
                else off = i * cb->ndimentions + j; // lookup type 2

                cb->dimentions[i * cb->ndimentions + j] = last + cb->min + cb->quantlist[off] * cb->delta;
                if (cb->seq_p) last = cb->dimentions[i * cb->ndimentions + j];
                div *= vals;
            }
        }
    }

}

static void ready_residue(residue_t * rc, venc_context_t * venc) {
    int i;
    assert(rc->type == 2);
    rc->maxes = av_mallocz(sizeof(float[2]) * rc->classifications);
    for (i = 0; i < rc->classifications; i++) {
        int j;
        codebook_t * cb;
        for (j = 0; j < 8; j++) if (rc->books[i][j] != -1) break;
        if (j == 8) continue; // zero
        cb = &venc->codebooks[rc->books[i][j]];
        assert(cb->ndimentions >= 2);
        assert(cb->lookup);

        for (j = 0; j < cb->nentries; j++) {
            float a;
            if (!cb->lens[j]) continue;
            a = fabs(cb->dimentions[j * cb->ndimentions]);
            if (a > rc->maxes[i][0]) rc->maxes[i][0] = a;
            a = fabs(cb->dimentions[j * cb->ndimentions + 1]);
            if (a > rc->maxes[i][1]) rc->maxes[i][1] = a;
        }
    }
    // small bias
    for (i = 0; i < rc->classifications; i++) {
        rc->maxes[i][0] += 0.8;
        rc->maxes[i][1] += 0.8;
    }
}

static void create_vorbis_context(venc_context_t * venc, AVCodecContext * avccontext) {
    floor_t * fc;
    residue_t * rc;
    mapping_t * mc;
    int i, book;

    venc->channels = avccontext->channels;
    venc->sample_rate = avccontext->sample_rate;
    venc->blocksize[0] = venc->blocksize[1] = 11;

    venc->ncodebooks = sizeof(cvectors)/sizeof(cvectors[0]);
    venc->codebooks = av_malloc(sizeof(codebook_t) * venc->ncodebooks);

    // codebook 0..14 - floor1 book, values 0..255
    // codebook 15 residue masterbook
    // codebook 16..29 residue
    for (book = 0; book < venc->ncodebooks; book++) {
        codebook_t * cb = &venc->codebooks[book];
        int vals;
        cb->ndimentions = cvectors[book].dim;
        cb->nentries = cvectors[book].real_len;
        cb->min = cvectors[book].min;
        cb->delta = cvectors[book].delta;
        cb->lookup = cvectors[book].lookup;
        cb->seq_p = 0;

        cb->lens = av_malloc(sizeof(uint8_t) * cb->nentries);
        cb->codewords = av_malloc(sizeof(uint32_t) * cb->nentries);
        memcpy(cb->lens, cvectors[book].clens, cvectors[book].len);
        memset(cb->lens + cvectors[book].len, 0, cb->nentries - cvectors[book].len);

        if (cb->lookup) {
            vals = cb_lookup_vals(cb->lookup, cb->ndimentions, cb->nentries);
            cb->quantlist = av_malloc(sizeof(int) * vals);
            for (i = 0; i < vals; i++) cb->quantlist[i] = cvectors[book].quant[i];
        } else {
            cb->quantlist = NULL;
        }
        ready_codebook(cb);
    }

    venc->nfloors = 1;
    venc->floors = av_malloc(sizeof(floor_t) * venc->nfloors);

    // just 1 floor
    fc = &venc->floors[0];
    fc->partitions = 8;
    fc->partition_to_class = av_malloc(sizeof(int) * fc->partitions);
    fc->nclasses = 0;
    for (i = 0; i < fc->partitions; i++) {
        int a[] = {0,1,2,2,3,3,4,4};
        fc->partition_to_class[i] = a[i];
        fc->nclasses = FFMAX(fc->nclasses, fc->partition_to_class[i]);
    }
    fc->nclasses++;
    fc->classes = av_malloc(sizeof(floor_class_t) * fc->nclasses);
    for (i = 0; i < fc->nclasses; i++) {
        floor_class_t * c = &fc->classes[i];
        int j, books;
        int dim[] = {3,4,3,4,3};
        int subclass[] = {0,1,1,2,2};
        int masterbook[] = {0/*none*/,0,1,2,3};
        int * nbooks[] = {
            (int[]){ 4 },
            (int[]){ 5, 6 },
            (int[]){ 7, 8 },
            (int[]){ -1, 9, 10, 11 },
            (int[]){ -1, 12, 13, 14 },
        };
        c->dim = dim[i];
        c->subclass = subclass[i];
        c->masterbook = masterbook[i];
        books = (1 << c->subclass);
        c->books = av_malloc(sizeof(int) * books);
        for (j = 0; j < books; j++) c->books[j] = nbooks[i][j];
    }
    fc->multiplier = 2;
    fc->rangebits = venc->blocksize[0] - 1;

    fc->values = 2;
    for (i = 0; i < fc->partitions; i++)
        fc->values += fc->classes[fc->partition_to_class[i]].dim;

    fc->list = av_malloc(sizeof(floor1_entry_t) * fc->values);
    fc->list[0].x = 0;
    fc->list[1].x = 1 << fc->rangebits;
    for (i = 2; i < fc->values; i++) {
        static const int a[] = {
             93, 23,372,  6, 46,186,750, 14, 33, 65,
            130,260,556,  3, 10, 18, 28, 39, 55, 79,
            111,158,220,312,464,650,850
        };
        fc->list[i].x = a[i - 2];
    }
    ff_vorbis_ready_floor1_list(fc->list, fc->values);

    venc->nresidues = 1;
    venc->residues = av_malloc(sizeof(residue_t) * venc->nresidues);

    // single residue
    rc = &venc->residues[0];
    rc->type = 2;
    rc->begin = 0;
    rc->end = 1600;
    rc->partition_size = 32;
    rc->classifications = 10;
    rc->classbook = 15;
    rc->books = av_malloc(sizeof(*rc->books) * rc->classifications);
    {
        static const int8_t a[10][8] = {
            { -1, -1, -1, -1, -1, -1, -1, -1, },
            { -1, -1, 16, -1, -1, -1, -1, -1, },
            { -1, -1, 17, -1, -1, -1, -1, -1, },
            { -1, -1, 18, -1, -1, -1, -1, -1, },
            { -1, -1, 19, -1, -1, -1, -1, -1, },
            { -1, -1, 20, -1, -1, -1, -1, -1, },
            { -1, -1, 21, -1, -1, -1, -1, -1, },
            { 22, 23, -1, -1, -1, -1, -1, -1, },
            { 24, 25, -1, -1, -1, -1, -1, -1, },
            { 26, 27, 28, -1, -1, -1, -1, -1, },
        };
    	memcpy(rc->books, a, sizeof a);
    }
    ready_residue(rc, venc);

    venc->nmappings = 1;
    venc->mappings = av_malloc(sizeof(mapping_t) * venc->nmappings);

    // single mapping
    mc = &venc->mappings[0];
    mc->submaps = 1;
    mc->mux = av_malloc(sizeof(int) * venc->channels);
    for (i = 0; i < venc->channels; i++) mc->mux[i] = 0;
    mc->floor = av_malloc(sizeof(int) * mc->submaps);
    mc->residue = av_malloc(sizeof(int) * mc->submaps);
    for (i = 0; i < mc->submaps; i++) {
        mc->floor[i] = 0;
        mc->residue[i] = 0;
    }
    mc->coupling_steps = venc->channels == 2 ? 1 : 0;
    mc->magnitude = av_malloc(sizeof(int) * mc->coupling_steps);
    mc->angle = av_malloc(sizeof(int) * mc->coupling_steps);
    if (mc->coupling_steps) {
        mc->magnitude[0] = 0;
        mc->angle[0] = 1;
    }

    venc->nmodes = 1;
    venc->modes = av_malloc(sizeof(vorbis_mode_t) * venc->nmodes);

    // single mode
    venc->modes[0].blockflag = 0;
    venc->modes[0].mapping = 0;

    venc->have_saved = 0;
    venc->saved = av_malloc(sizeof(float) * venc->channels * (1 << venc->blocksize[1]) / 2);
    venc->samples = av_malloc(sizeof(float) * venc->channels * (1 << venc->blocksize[1]));
    venc->floor = av_malloc(sizeof(float) * venc->channels * (1 << venc->blocksize[1]) / 2);
    venc->coeffs = av_malloc(sizeof(float) * venc->channels * (1 << venc->blocksize[1]) / 2);

    venc->win[0] = ff_vorbis_vwin[venc->blocksize[0] - 6];
    venc->win[1] = ff_vorbis_vwin[venc->blocksize[1] - 6];

    ff_mdct_init(&venc->mdct[0], venc->blocksize[0], 0);
    ff_mdct_init(&venc->mdct[1], venc->blocksize[1], 0);
}

static void put_float(PutBitContext * pb, float f) {
    int exp, mant;
    uint32_t res = 0;
    mant = (int)ldexp(frexp(f, &exp), 20);
    exp += 788 - 20;
    if (mant < 0) { res |= (1 << 31); mant = -mant; }
    res |= mant | (exp << 21);
    put_bits(pb, 32, res);
}

static void put_codebook_header(PutBitContext * pb, codebook_t * cb) {
    int i;
    int ordered = 0;

    put_bits(pb, 24, 0x564342); //magic
    put_bits(pb, 16, cb->ndimentions);
    put_bits(pb, 24, cb->nentries);

    for (i = 1; i < cb->nentries; i++) if (cb->lens[i] < cb->lens[i-1]) break;
    if (i == cb->nentries) ordered = 1;

    put_bits(pb, 1, ordered);
    if (ordered) {
        int len = cb->lens[0];
        put_bits(pb, 5, len - 1);
        i = 0;
        while (i < cb->nentries) {
            int j;
            for (j = 0; j+i < cb->nentries; j++) if (cb->lens[j+i] != len) break;
            put_bits(pb, ilog(cb->nentries - i), j);
            i += j;
            len++;
        }
    } else {
        int sparse = 0;
        for (i = 0; i < cb->nentries; i++) if (!cb->lens[i]) break;
        if (i != cb->nentries) sparse = 1;
        put_bits(pb, 1, sparse);

        for (i = 0; i < cb->nentries; i++) {
            if (sparse) put_bits(pb, 1, !!cb->lens[i]);
            if (cb->lens[i]) put_bits(pb, 5, cb->lens[i] - 1);
        }
    }

    put_bits(pb, 4, cb->lookup);
    if (cb->lookup) {
        int tmp = cb_lookup_vals(cb->lookup, cb->ndimentions, cb->nentries);
        int bits = ilog(cb->quantlist[0]);

        for (i = 1; i < tmp; i++) bits = FFMAX(bits, ilog(cb->quantlist[i]));

        put_float(pb, cb->min);
        put_float(pb, cb->delta);

        put_bits(pb, 4, bits - 1);
        put_bits(pb, 1, cb->seq_p);

        for (i = 0; i < tmp; i++) put_bits(pb, bits, cb->quantlist[i]);
    }
}

static void put_floor_header(PutBitContext * pb, floor_t * fc) {
    int i;

    put_bits(pb, 16, 1); // type, only floor1 is supported

    put_bits(pb, 5, fc->partitions);

    for (i = 0; i < fc->partitions; i++) put_bits(pb, 4, fc->partition_to_class[i]);

    for (i = 0; i < fc->nclasses; i++) {
        int j, books;

        put_bits(pb, 3, fc->classes[i].dim - 1);
        put_bits(pb, 2, fc->classes[i].subclass);

        if (fc->classes[i].subclass) put_bits(pb, 8, fc->classes[i].masterbook);

        books = (1 << fc->classes[i].subclass);

        for (j = 0; j < books; j++) put_bits(pb, 8, fc->classes[i].books[j] + 1);
    }

    put_bits(pb, 2, fc->multiplier - 1);
    put_bits(pb, 4, fc->rangebits);

    for (i = 2; i < fc->values; i++) put_bits(pb, fc->rangebits, fc->list[i].x);
}

static void put_residue_header(PutBitContext * pb, residue_t * rc) {
    int i;

    put_bits(pb, 16, rc->type);

    put_bits(pb, 24, rc->begin);
    put_bits(pb, 24, rc->end);
    put_bits(pb, 24, rc->partition_size - 1);
    put_bits(pb, 6, rc->classifications - 1);
    put_bits(pb, 8, rc->classbook);

    for (i = 0; i < rc->classifications; i++) {
        int j, tmp = 0;
        for (j = 0; j < 8; j++) tmp |= (rc->books[i][j] != -1) << j;

        put_bits(pb, 3, tmp & 7);
        put_bits(pb, 1, tmp > 7);

        if (tmp > 7) put_bits(pb, 5, tmp >> 3);
    }

    for (i = 0; i < rc->classifications; i++) {
        int j;
        for (j = 0; j < 8; j++)
            if (rc->books[i][j] != -1)
                put_bits(pb, 8, rc->books[i][j]);
    }
}

static int put_main_header(venc_context_t * venc, uint8_t ** out) {
    int i;
    PutBitContext pb;
    uint8_t buffer[50000] = {0}, * p = buffer;
    int buffer_len = sizeof buffer;
    int len, hlens[3];

    // identification header
    init_put_bits(&pb, p, buffer_len);
    put_bits(&pb, 8, 1); //magic
    for (i = 0; "vorbis"[i]; i++) put_bits(&pb, 8, "vorbis"[i]);
    put_bits(&pb, 32, 0); // version
    put_bits(&pb, 8, venc->channels);
    put_bits(&pb, 32, venc->sample_rate);
    put_bits(&pb, 32, 0); // bitrate
    put_bits(&pb, 32, 0); // bitrate
    put_bits(&pb, 32, 0); // bitrate
    put_bits(&pb, 4, venc->blocksize[0]);
    put_bits(&pb, 4, venc->blocksize[1]);
    put_bits(&pb, 1, 1); // framing

    flush_put_bits(&pb);
    hlens[0] = (put_bits_count(&pb) + 7) / 8;
    buffer_len -= hlens[0];
    p += hlens[0];

    // comment header
    init_put_bits(&pb, p, buffer_len);
    put_bits(&pb, 8, 3); //magic
    for (i = 0; "vorbis"[i]; i++) put_bits(&pb, 8, "vorbis"[i]);
    put_bits(&pb, 32, 0); // vendor length TODO
    put_bits(&pb, 32, 0); // amount of comments
    put_bits(&pb, 1, 1); // framing

    flush_put_bits(&pb);
    hlens[1] = (put_bits_count(&pb) + 7) / 8;
    buffer_len -= hlens[1];
    p += hlens[1];

    // setup header
    init_put_bits(&pb, p, buffer_len);
    put_bits(&pb, 8, 5); //magic
    for (i = 0; "vorbis"[i]; i++) put_bits(&pb, 8, "vorbis"[i]);

    // codebooks
    put_bits(&pb, 8, venc->ncodebooks - 1);
    for (i = 0; i < venc->ncodebooks; i++) put_codebook_header(&pb, &venc->codebooks[i]);

    // time domain, reserved, zero
    put_bits(&pb, 6, 0);
    put_bits(&pb, 16, 0);

    // floors
    put_bits(&pb, 6, venc->nfloors - 1);
    for (i = 0; i < venc->nfloors; i++) put_floor_header(&pb, &venc->floors[i]);

    // residues
    put_bits(&pb, 6, venc->nresidues - 1);
    for (i = 0; i < venc->nresidues; i++) put_residue_header(&pb, &venc->residues[i]);

    // mappings
    put_bits(&pb, 6, venc->nmappings - 1);
    for (i = 0; i < venc->nmappings; i++) {
        mapping_t * mc = &venc->mappings[i];
        int j;
        put_bits(&pb, 16, 0); // mapping type

        put_bits(&pb, 1, mc->submaps > 1);
        if (mc->submaps > 1) put_bits(&pb, 4, mc->submaps - 1);

        put_bits(&pb, 1, !!mc->coupling_steps);
        if (mc->coupling_steps) {
            put_bits(&pb, 8, mc->coupling_steps - 1);
            for (j = 0; j < mc->coupling_steps; j++) {
                put_bits(&pb, ilog(venc->channels - 1), mc->magnitude[j]);
                put_bits(&pb, ilog(venc->channels - 1), mc->angle[j]);
            }
        }

        put_bits(&pb, 2, 0); // reserved

        if (mc->submaps > 1) for (j = 0; j < venc->channels; j++) put_bits(&pb, 4, mc->mux[j]);

        for (j = 0; j < mc->submaps; j++) {
            put_bits(&pb, 8, 0); // reserved time configuration
            put_bits(&pb, 8, mc->floor[j]);
            put_bits(&pb, 8, mc->residue[j]);
        }
    }

    // modes
    put_bits(&pb, 6, venc->nmodes - 1);
    for (i = 0; i < venc->nmodes; i++) {
        put_bits(&pb, 1, venc->modes[i].blockflag);
        put_bits(&pb, 16, 0); // reserved window type
        put_bits(&pb, 16, 0); // reserved transform type
        put_bits(&pb, 8, venc->modes[i].mapping);
    }

    put_bits(&pb, 1, 1); // framing

    flush_put_bits(&pb);
    hlens[2] = (put_bits_count(&pb) + 7) / 8;

    len = hlens[0] + hlens[1] + hlens[2];
    p = *out = av_mallocz(64 + len + len/255);

    *p++ = 2;
    p += av_xiphlacing(p, hlens[0]);
    p += av_xiphlacing(p, hlens[1]);
    buffer_len = 0;
    for (i = 0; i < 3; i++) {
        memcpy(p, buffer + buffer_len, hlens[i]);
        p += hlens[i];
        buffer_len += hlens[i];
    }

    return p - *out;
}

static float get_floor_average(floor_t * fc, float * coeffs, int i) {
    int begin = fc->list[fc->list[FFMAX(i-1, 0)].sort].x;
    int end   = fc->list[fc->list[FFMIN(i+1, fc->values - 1)].sort].x;
    int j;
    float average = 0;

    for (j = begin; j < end; j++) average += fabs(coeffs[j]);
    return average / (end - begin);
}

static void floor_fit(venc_context_t * venc, floor_t * fc, float * coeffs, int * posts, int samples) {
    int range = 255 / fc->multiplier + 1;
    int i;
    float tot_average = 0.;
    float averages[fc->values];
    for (i = 0; i < fc->values; i++) tot_average += averages[i] = get_floor_average(fc, coeffs, i);
    tot_average /= fc->values;
    tot_average /= venc->quality;

    for (i = 0; i < fc->values; i++) {
        int position = fc->list[fc->list[i].sort].x;
        float average = averages[i];
        int j;

        average /= pow(average, 0.5) / tot_average * pow(0.8, position/200.); // MAGIC!
        for (j = 0; j < range - 1; j++) if (ff_vorbis_floor1_inverse_db_table[j * fc->multiplier] > average) break;
        posts[fc->list[i].sort] = j;
    }
}

static int render_point(int x0, int y0, int x1, int y1, int x) {
    return y0 +  (x - x0) * (y1 - y0) / (x1 - x0);
}

static void floor_encode(venc_context_t * venc, floor_t * fc, PutBitContext * pb, int * posts, float * floor, int samples) {
    int range = 255 / fc->multiplier + 1;
    int coded[fc->values]; // first 2 values are unused
    int i, counter;
    int lx, ly;

    put_bits(pb, 1, 1); // non zero
    put_bits(pb, ilog(range - 1), posts[0]);
    put_bits(pb, ilog(range - 1), posts[1]);
    coded[0] = coded[1] = 1;

    for (i = 2; i < fc->values; i++) {
        int predicted = render_point(fc->list[fc->list[i].low].x,
                                     posts[fc->list[i].low],
                                     fc->list[fc->list[i].high].x,
                                     posts[fc->list[i].high],
                                     fc->list[i].x);
        int highroom = range - predicted;
        int lowroom = predicted;
        int room = FFMIN(highroom, lowroom);
        if (predicted == posts[i]) {
            coded[i] = 0; // must be used later as flag!
            continue;
        } else {
            if (!coded[fc->list[i].low]) coded[fc->list[i].low] = -1;
            if (!coded[fc->list[i].high]) coded[fc->list[i].high] = -1;
        }
        if (posts[i] > predicted) {
            if (posts[i] - predicted > room) coded[i] = posts[i] - predicted + lowroom;
            else coded[i] = (posts[i] - predicted) << 1;
        } else {
            if (predicted - posts[i] > room) coded[i] = predicted - posts[i] + highroom - 1;
            else coded[i] = ((predicted - posts[i]) << 1) - 1;
        }
    }

    counter = 2;
    for (i = 0; i < fc->partitions; i++) {
        floor_class_t * c = &fc->classes[fc->partition_to_class[i]];
        int k, cval = 0, csub = 1<<c->subclass;
        if (c->subclass) {
            codebook_t * book = &venc->codebooks[c->masterbook];
            int cshift = 0;
            for (k = 0; k < c->dim; k++) {
                int l;
                for (l = 0; l < csub; l++) {
                    int maxval = 1;
                    if (c->books[l] != -1) maxval = venc->codebooks[c->books[l]].nentries;
                    // coded could be -1, but this still works, cause thats 0
                    if (coded[counter + k] < maxval) break;
                }
                assert(l != csub);
                cval |= l << cshift;
                cshift += c->subclass;
            }
            put_codeword(pb, book, cval);
        }
        for (k = 0; k < c->dim; k++) {
            int book = c->books[cval & (csub-1)];
            int entry = coded[counter++];
            cval >>= c->subclass;
            if (book == -1) continue;
            if (entry == -1) entry = 0;
            put_codeword(pb, &venc->codebooks[book], entry);
        }
    }

    lx = 0;
    ly = posts[0] * fc->multiplier; // sorted 0 is still 0
    for (i = 1; i < fc->values; i++) {
        int pos = fc->list[i].sort;
        if (coded[pos]) {
            render_line(lx, ly, fc->list[pos].x, posts[pos] * fc->multiplier, floor, samples);
            lx = fc->list[pos].x;
            ly = posts[pos] * fc->multiplier;
        }
        if (lx >= samples) break;
    }
    if (lx < samples) render_line(lx, ly, samples, ly, floor, samples);
}

static float * put_vector(codebook_t * book, PutBitContext * pb, float * num) {
    int i;
    int entry = -1;
    float distance = 0;
    assert(book->dimentions);
    for (i = 0; i < book->nentries; i++) {
        float d = 0.;
        int j;
        if (!book->lens[i]) continue;
        for (j = 0; j < book->ndimentions; j++) {
            float a = (book->dimentions[i * book->ndimentions + j] - num[j]);
            d += a*a;
        }
        if (entry == -1 || distance > d) {
            entry = i;
            distance = d;
        }
    }
    put_codeword(pb, book, entry);
    return &book->dimentions[entry * book->ndimentions];
}

static void residue_encode(venc_context_t * venc, residue_t * rc, PutBitContext * pb, float * coeffs, int samples, int real_ch) {
    int pass, i, j, p, k;
    int psize = rc->partition_size;
    int partitions = (rc->end - rc->begin) / psize;
    int channels = (rc->type == 2) ? 1 : real_ch;
    int classes[channels][partitions];
    int classwords = venc->codebooks[rc->classbook].ndimentions;

    assert(rc->type == 2);
    assert(real_ch == 2);
    for (p = 0; p < partitions; p++) {
        float max1 = 0., max2 = 0.;
        int s = rc->begin + p * psize;
        for (k = s; k < s + psize; k += 2) {
            max1 = FFMAX(max1, fabs(coeffs[          k / real_ch]));
            max2 = FFMAX(max2, fabs(coeffs[samples + k / real_ch]));
        }

        for (i = 0; i < rc->classifications - 1; i++) {
            if (max1 < rc->maxes[i][0] && max2 < rc->maxes[i][1]) break;
        }
        classes[0][p] = i;
    }

    for (pass = 0; pass < 8; pass++) {
        p = 0;
        while (p < partitions) {
            if (pass == 0) for (j = 0; j < channels; j++) {
                codebook_t * book = &venc->codebooks[rc->classbook];
                int entry = 0;
                for (i = 0; i < classwords; i++) {
                    entry *= rc->classifications;
                    entry += classes[j][p + i];
                }
                put_codeword(pb, book, entry);
            }
            for (i = 0; i < classwords && p < partitions; i++, p++) {
                for (j = 0; j < channels; j++) {
                    int nbook = rc->books[classes[j][p]][pass];
                    codebook_t * book = &venc->codebooks[nbook];
                    float * buf = coeffs + samples*j + rc->begin + p*psize;
                    if (nbook == -1) continue;

                    assert(rc->type == 0 || rc->type == 2);
                    assert(!(psize % book->ndimentions));

                    if (rc->type == 0) {
                        for (k = 0; k < psize; k += book->ndimentions) {
                            float * a = put_vector(book, pb, &buf[k]);
                            int l;
                            for (l = 0; l < book->ndimentions; l++) buf[k + l] -= a[l];
                        }
                    } else {
                        for (k = 0; k < psize; k += book->ndimentions) {
                            int dim = book->ndimentions, s = rc->begin + p * psize + k, l;
                            float vec[dim], * a = vec;
                            for (l = s; l < s + dim; l++)
                                *a++ = coeffs[(l % real_ch) * samples + l / real_ch];
                            a = put_vector(book, pb, vec);
                            for (l = s; l < s + dim; l++)
                                coeffs[(l % real_ch) * samples + l / real_ch] -= *a++;
                        }
                    }
                }
            }
        }
    }
}

static int apply_window_and_mdct(venc_context_t * venc, signed short * audio, int samples) {
    int i, j, channel;
    const float * win = venc->win[0];
    int window_len = 1 << (venc->blocksize[0] - 1);
    float n = (float)(1 << venc->blocksize[0]) / 4.;
    // FIXME use dsp

    if (!venc->have_saved && !samples) return 0;

    if (venc->have_saved) {
        for (channel = 0; channel < venc->channels; channel++) {
            memcpy(venc->samples + channel*window_len*2, venc->saved + channel*window_len, sizeof(float)*window_len);
        }
    } else {
        for (channel = 0; channel < venc->channels; channel++) {
            memset(venc->samples + channel*window_len*2, 0, sizeof(float)*window_len);
        }
    }

    if (samples) {
        for (channel = 0; channel < venc->channels; channel++) {
            float * offset = venc->samples + channel*window_len*2 + window_len;
            j = channel;
            for (i = 0; i < samples; i++, j += venc->channels)
                offset[i] = audio[j] / 32768. / n * win[window_len - i - 1];
        }
    } else {
        for (channel = 0; channel < venc->channels; channel++) {
            memset(venc->samples + channel*window_len*2 + window_len, 0, sizeof(float)*window_len);
        }
    }

    for (channel = 0; channel < venc->channels; channel++) {
        ff_mdct_calc(&venc->mdct[0], venc->coeffs + channel*window_len, venc->samples + channel*window_len*2, venc->floor/*tmp*/);
    }

    if (samples) {
        for (channel = 0; channel < venc->channels; channel++) {
            float * offset = venc->saved + channel*window_len;
            j = channel;
            for (i = 0; i < samples; i++, j += venc->channels)
                offset[i] = audio[j] / 32768. / n * win[i];
        }
        venc->have_saved = 1;
    } else {
        venc->have_saved = 0;
    }
    return 1;
}

static int vorbis_encode_init(AVCodecContext * avccontext)
{
    venc_context_t * venc = avccontext->priv_data;

    if (avccontext->channels != 2) return -1;

    create_vorbis_context(venc, avccontext);

    if (avccontext->flags & CODEC_FLAG_QSCALE) venc->quality = avccontext->global_quality / (float)FF_QP2LAMBDA / 100.;
    else venc->quality = 0.17;

    avccontext->extradata_size = put_main_header(venc, (uint8_t**)&avccontext->extradata);

    avccontext->frame_size = 1 << (venc->blocksize[0] - 1);

    avccontext->coded_frame = avcodec_alloc_frame();
    avccontext->coded_frame->key_frame = 1;

    return 0;
}

static int vorbis_encode_frame(AVCodecContext * avccontext, unsigned char * packets, int buf_size, void *data)
{
    venc_context_t * venc = avccontext->priv_data;
    signed short * audio = data;
    int samples = data ? avccontext->frame_size : 0;
    vorbis_mode_t * mode;
    mapping_t * mapping;
    PutBitContext pb;
    int i;

    if (!apply_window_and_mdct(venc, audio, samples)) return 0;
    samples = 1 << (venc->blocksize[0] - 1);

    init_put_bits(&pb, packets, buf_size);

    put_bits(&pb, 1, 0); // magic bit

    put_bits(&pb, ilog(venc->nmodes - 1), 0); // 0 bits, the mode

    mode = &venc->modes[0];
    mapping = &venc->mappings[mode->mapping];
    if (mode->blockflag) {
        put_bits(&pb, 1, 0);
        put_bits(&pb, 1, 0);
    }

    for (i = 0; i < venc->channels; i++) {
        floor_t * fc = &venc->floors[mapping->floor[mapping->mux[i]]];
        int posts[fc->values];
        floor_fit(venc, fc, &venc->coeffs[i * samples], posts, samples);
        floor_encode(venc, fc, &pb, posts, &venc->floor[i * samples], samples);
    }

    for (i = 0; i < venc->channels; i++) {
        int j;
        for (j = 0; j < samples; j++) {
            venc->coeffs[i * samples + j] /= venc->floor[i * samples + j];
        }
    }

    for (i = 0; i < mapping->coupling_steps; i++) {
        float * mag = venc->coeffs + mapping->magnitude[i] * samples;
        float * ang = venc->coeffs + mapping->angle[i] * samples;
        int j;
        for (j = 0; j < samples; j++) {
            float a = ang[j];
            ang[j] -= mag[j];
            if (mag[j] > 0) ang[j] = -ang[j];
            if (ang[j] < 0) mag[j] = a;
        }
    }

    residue_encode(venc, &venc->residues[mapping->residue[mapping->mux[0]]], &pb, venc->coeffs, samples, venc->channels);

    flush_put_bits(&pb);
    return (put_bits_count(&pb) + 7) / 8;
}


static int vorbis_encode_close(AVCodecContext * avccontext)
{
    venc_context_t * venc = avccontext->priv_data;
    int i;

    if (venc->codebooks) for (i = 0; i < venc->ncodebooks; i++) {
        av_freep(&venc->codebooks[i].lens);
        av_freep(&venc->codebooks[i].codewords);
        av_freep(&venc->codebooks[i].quantlist);
        av_freep(&venc->codebooks[i].dimentions);
    }
    av_freep(&venc->codebooks);

    if (venc->floors) for (i = 0; i < venc->nfloors; i++) {
        int j;
        av_freep(&venc->floors[i].classes);
        if (venc->floors[i].classes)
            for (j = 0; j < venc->floors[i].nclasses; j++)
                av_freep(&venc->floors[i].classes[j].books);
        av_freep(&venc->floors[i].partition_to_class);
        av_freep(&venc->floors[i].list);
    }
    av_freep(&venc->floors);

    if (venc->residues) for (i = 0; i < venc->nresidues; i++) {
        av_freep(&venc->residues[i].books);
        av_freep(&venc->residues[i].maxes);
    }
    av_freep(&venc->residues);

    if (venc->mappings) for (i = 0; i < venc->nmappings; i++) {
        av_freep(&venc->mappings[i].mux);
        av_freep(&venc->mappings[i].floor);
        av_freep(&venc->mappings[i].residue);
    }
    av_freep(&venc->mappings);

    av_freep(&venc->modes);

    av_freep(&venc->saved);
    av_freep(&venc->samples);
    av_freep(&venc->floor);
    av_freep(&venc->coeffs);

    ff_mdct_end(&venc->mdct[0]);
    ff_mdct_end(&venc->mdct[1]);

    av_freep(&avccontext->coded_frame);
    av_freep(&avccontext->extradata);

    return 0 ;
}

AVCodec vorbis_encoder = {
    "vorbis",
    CODEC_TYPE_AUDIO,
    CODEC_ID_VORBIS,
    sizeof(venc_context_t),
    vorbis_encode_init,
    vorbis_encode_frame,
    vorbis_encode_close,
    .capabilities= CODEC_CAP_DELAY,
};
