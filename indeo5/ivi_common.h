/*
 * common functions for Indeo Video Interactive codecs (indeo4 and indeo5)
 *
 * Copyright (c) 2009 Maxim Poliakovski
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
 * @file libavcodec/ivi_common.h
 * This file contains structures and macros shared by both Indeo4 and
 * Indeo5 decoders.
 */

#ifndef AVCODEC_IVI_COMMON_H
#define AVCODEC_IVI_COMMON_H

#include <stdint.h>

#define IVI_DEBUG

#define IVI_VLC_BITS 13 ///< max number of bits of the ivi's huffman codes

/**
 *  huffman codebook descriptor
 */
typedef struct {
    int32_t     num_rows;
    uint8_t     xbits[16];
} IVIHuffDesc;

/* see ivi_common.c for a definition */
extern const IVIHuffDesc ff_ivi_mb_huff_desc[8];  ///< static macroblock huffman tables
extern const IVIHuffDesc ff_ivi_blk_huff_desc[8]; ///< static block huffman tables


/**
 *  Run-value (RLE) table descriptor
 */
typedef struct {
    uint8_t     eob_sym; ///< huffman symbol for EOB (end of block)
    uint8_t     esc_sym; ///< huffman symbol for ESC (escape)
    uint8_t     runtab[256];
    int8_t      valtab[256];
} RVMapDesc;

extern const RVMapDesc ff_ivi_rvmap_tabs[9]; /* defined in ivi_common.c */


/**
 *  This structure describes an indeo macroblock (16x16, 8x8 or 4x4).
 */
typedef struct {
    int16_t     xpos;
    int16_t     ypos;
    uint32_t    buf_offs; ///< address in the output buffer for this mb
    uint8_t     type;     ///< macroblock type: 0 - INTRA, 1 - INTER
    uint8_t     cbp;      ///< coded block pattern
    uint8_t     q_delta;  ///< quant delta
    int8_t      mv_x;     ///< motion vector (x component)
    int8_t      mv_y;     ///< motion vector (y component)
} IVIMbInfo;


/**
 *  This structure describes an indeo tile.
 */
typedef struct {
    uint32_t    xpos;
    uint32_t    ypos;
    uint32_t    width;
    uint32_t    height;
    int32_t     is_empty;  ///< = 1 if this tile doesn't contain any data
    uint32_t    data_size; ///< size of the data in bytes
    uint32_t    num_MBs;   ///< number of macroblocks in this tile
    IVIMbInfo   *mbs;      ///< array of macroblock descriptors
    IVIMbInfo   *ref_mbs;  ///< ptr to the macroblock descriptors of the reference tile
} IVITile;


/**
 *  This structure describes an indeo wavelet band.
 */
typedef struct {
    uint8_t         plane;          ///< plane number this band belongs to
    uint8_t         band_num;       ///< band number
    uint32_t        width;
    uint32_t        height;
    const uint8_t   *data_ptr;      ///< ptr to the first byte of the band data
    uint32_t        data_size;      ///< size of the band data
    int16_t         *buf;           ///< ptr to the output buffer for this band
    int16_t         *ref_buf;       ///< ptr to the reference frame buffer for motion compensation
    int16_t         *buf1;          ///< primary   band buffer
    int16_t         *buf2;          ///< secondary band buffer
    uint32_t        pitch;          ///< pitch associated with the buffers above
    uint8_t         is_empty;       ///< = 1 if this band doesn't contain any data
    uint8_t         mb_size;        ///< macroblock size
    uint8_t         blk_size;       ///< block size
    uint8_t         mc_resolution;  ///< resolution of the motion compensation: 0 - fullpel, 1 - halfpel
    int8_t          inherit_mv;
    int8_t          inherit_qdelta;
    int8_t          qdelta_present; ///< tells if Qdelta signal is present in the bitstream (indeo5 only)
    uint8_t         quant_mat;      ///< dequant matrix
    uint8_t         glob_quant;     ///< quant base for this band
    const uint8_t   *scan;          ///< ptr to the scan pattern

    uint8_t         huff_sel;       ///< huffman table for this band
    IVIHuffDesc     huff_desc;      ///< table descriptor associated with the selector above
    VLC             *blk_vlc;       ///< ptr to the vlc table for decoding block data
    VLC             blk_vlc_cust;   ///< custom block vlc table

    uint16_t        *dequant_intra; ///< ptr to dequant tables for intra blocks
    uint16_t        *dequant_inter; ///< ptr dequant tables for inter blocks
    uint8_t         num_corr;       ///< number of correction entries
    uint8_t         corr[61*2];     ///< rvmap correction pairs
    uint8_t         rvmap_sel;      ///< rvmap table selector
    RVMapDesc       *rv_map;        ///< ptr to the RLE table for this band
    uint16_t        num_tiles;      ///< number of tiles in this band
    IVITile         *tiles;         ///< array of tile descriptors
    void (*inv_transform)(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags); ///< ptr to inverse transform function
    void (*dc_transform) (int32_t *in, int16_t *out, uint32_t pitch, int blk_size);   ///< ptr to dc transform function or NULL
    uint8_t         is_2d_trans;    ///< 1 indicates that the two-dimensional inverse transform
#ifdef IVI_DEBUG
    uint16_t        checksum; ///< for debug purposes
    int32_t         checksum_present;
    uint32_t        bufsize; ///< band buffer size in bytes
#endif
    const uint8_t   *intra_base;
    const uint8_t   *inter_base;
    const uint8_t   *intra_scale;
    const uint8_t   *inter_scale;
} IVIBandDesc;


/**
 *  This structure describes a color plane (luma or chroma).
 */
typedef struct {
    uint16_t    width;
    uint16_t    height;
   // uint32_t    pitch;
    uint8_t     buf_switch; ///< used to switch between two buffers
    //int16_t     *buf1;      ///< primary buffer to store decoded pixels
    //int16_t     *buf2;      ///< secondary buffer to store decoded pixels
    uint8_t     num_bands;  ///< number of bands this plane subdivided into
    IVIBandDesc *bands;     ///< array of band descriptors
} IVIPlaneDesc;


typedef struct {
    uint16_t    pic_width;
    uint16_t    pic_height;
    uint16_t    chroma_width;
    uint16_t    chroma_height;
    uint16_t    tile_width;
    uint16_t    tile_height;
    uint8_t     luma_bands;
    uint8_t     chroma_bands;
} IVIPicConfig;

static inline int ivi_pic_config_cmp(IVIPicConfig *str1, IVIPicConfig *str2)
{
    return (str1->pic_width    != str2->pic_width    || str1->pic_height    != str2->pic_height    ||
            str1->chroma_width != str2->chroma_width || str1->chroma_height != str2->chroma_height ||
            str1->tile_width   != str2->tile_width   || str1->tile_height   != str2->tile_height   ||
            str1->luma_bands   != str2->luma_bands   || str1->chroma_bands  != str2->chroma_bands);
}

/** calculate number of tiles in a stride */
#define IVI_NUM_TILES(stride, tile_size) (((stride) + (tile_size) - 1) / (tile_size))

/** calculate number of macroblocks in a tile */
#define IVI_MBs_PER_TILE(tile_width, tile_height, mb_size) ((((tile_width) + (mb_size) - 1) / (mb_size)) * (((tile_height) + (mb_size) - 1) / (mb_size)))

/** convert unsigned values into signed ones (the sign is in the LSB) */
/* TODO: find a way to calculate this without the conditional using bit magic */
#define IVI_TOSIGNED(val) (((val) & 1) ? ((val) + 1) >> 1 : -(((val) + 1) >> 1))

/** divide the motion vector mv by 4 */
#define IVI_MV_DIV4(mv) (((mv) + 1 + ((mv) > 0))>>2)

/** divide the motion vector mv by 2 */
#define IVI_MV_DIV2(mv) (((mv) + ((mv) > 0))>>1)

int  ff_ivi_create_huff_from_desc(const IVIHuffDesc *cb, VLC *pOut, int flag);
int  ff_ivi_dec_huff_desc(GetBitContext *gb, IVIHuffDesc *desc);
int  ff_ivi_huff_desc_cmp(const IVIHuffDesc *desc1, const IVIHuffDesc *desc2);
void ff_ivi_huff_desc_copy(IVIHuffDesc *dst, const IVIHuffDesc *src);
int  ff_ivi_init_planes(IVIPlaneDesc *planes, const IVIPicConfig *cfg);
void ff_ivi_free_buffers(IVIPlaneDesc *planes);
int  ff_ivi_init_tiles(IVIPlaneDesc *planes, const int tile_width, const int tile_height);
int  ff_ivi_dec_tile_data_size(GetBitContext *gb);
int  ff_ivi_decode_blocks(GetBitContext *gb, IVIBandDesc *band, IVITile *tile);
void ff_ivi_process_empty_tile(AVCodecContext *avctx, IVIBandDesc *band, IVITile *tile, int32_t mv_scale);
void ff_ivi_output_plane(IVIPlaneDesc *plane, uint8_t *dst, int dst_pitch);
void ff_ivi_put_pixels_8x8(int32_t *in, int16_t *out, uint32_t pitch, uint8_t *flags);
void ff_ivi_put_dc_pixel_8x8(int32_t *in, int16_t *out, uint32_t pitch, int blk_size);

#ifdef IVI_DEBUG
uint16_t ivi_calc_band_checksum (IVIBandDesc *band);
int ivi_check_band (IVIBandDesc *band, uint8_t *ref, int pitch);
#endif

#endif /* AVCODEC_IVI_COMMON_H */
