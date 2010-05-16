/*
 * Indeo Video Interactive v4 compatible decoder
 * Copyright (c) 2009-2010 Maxim Poliakovski
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
 * @file libavcodec/indeo4.c
 * Indeo Video Interactive version 4 decoder
 *
 * Indeo4 data is usually transported within .avi or .mov files.
 * Known FOURCCs: 'IV41'
 */

#define IVI4_STREAM_ANALYSER    1
#define IVI4_DEBUG_COMPARE      1
#define IVI4_DEBUG_CHECKSUM     1

#define ALT_BITSTREAM_READER_LE
#include "avcodec.h"
#include "get_bits.h"
#include "dsputil.h"
#include "ivi_dsp.h"
#include "ivi_common.h"
#include "indeo4data.h"

#if IVI4_DEBUG_COMPARE
#include <dlfcn.h> /* dlsym, dlopen, dlclose */
#endif

#if IVI4_STREAM_ANALYSER
static uint8_t has_b_frames = 0;
static uint8_t has_transp   = 0;
static uint8_t uses_tiling  = 0;
static uint8_t uses_haar    = 0;
static uint8_t uses_fullpel = 0;
#endif

/**
 *  Indeo4 frame types.
 */
enum {
    FRAMETYPE_INTRA       = 0,
    FRAMETYPE_INTER       = 2,  ///< non-droppable P-frame
    FRAMETYPE_BIDIR       = 3,  ///< bidirectional frame
    FRAMETYPE_INTER_NOREF = 4,  ///< droppable P-frame
    FRAMETYPE_NULL        = 5   ///< empty frame with no data
};

#define IVI4_PIC_SIZE_ESC   7


typedef struct {
    GetBitContext   gb;
    AVFrame         frame;
    RVMapDesc       rvmap_tabs[9];   ///< local corrected copy of the static rvmap tables

    uint32_t        frame_num;
    int             frame_type;
    int             prev_frame_type; ///< frame type of the previous frame
    uint32_t        data_size;       ///< size of the frame data in bytes from picture header
    int             is_scalable;
    int             transp_status;   ///< transparency mode status: 1 - enabled

    IVIPicConfig    pic_conf;
    IVIPlaneDesc    planes[3];       ///< color planes

    IVIHuffTab      mb_vlc;          ///< current macroblock table descriptor
    IVIHuffTab      blk_vlc;         ///< current block table descriptor

    uint16_t        checksum;        ///< frame checksum

    uint8_t         rvmap_sel;
    uint8_t         in_imf;
    uint8_t         in_q;            ///< enforce quant
    uint8_t         pic_glob_quant;
    uint8_t         unknown1;
} IVI4DecContext;

/* FIXME: move those declarations into ivi_common.h or ivi_dsp.h */
typedef void (*inv_trans_func)(const int32_t *in, int16_t *out, uint32_t pitch, const uint8_t *flags);
typedef void (*dc_trans_func) (const int32_t *in, int16_t *out, uint32_t pitch, int blk_size);

struct {
    inv_trans_func  inv_trans;
    dc_trans_func   dc_trans;
    int             is_2d_trans;
} transforms[12] = {
    {NULL, NULL, 0}, /* inverse Haar 8x8 */
    {NULL, NULL, 0}, /* inverse Haar 8x1 */
    {NULL, NULL, 0}, /* inverse Haar 1x8 */
    {NULL, NULL, 0}, /* no transform 8x8 */
    {ff_ivi_inverse_slant_8x8, ff_ivi_dc_slant_2d, 1},
    {NULL, NULL, 0}, /* inverse Slant 8x1 */
    {NULL, NULL, 0}, /* inverse Slant 1x8 */
    {NULL, NULL, 0}, /* inverse DCT 8x8 */
    {NULL, NULL, 0}, /* inverse DCT 8x1 */
    {NULL, NULL, 0}, /* inverse DCT 1x8 */
    {NULL, NULL, 0}, /* inverse Haar 4x4 */
    {ff_ivi_inverse_slant_4x4, ff_ivi_dc_slant_2d, 1},
};

#if IVI4_DEBUG_COMPARE

#define INDEO41_DLL "/home/maxim/ffmpegnew/iv41.so"

/**
 *  Define some dummy functions needed to run the XANIM plugin
 */
typedef struct
{
    uint32_t    cmd;            /* decode or query */
    uint32_t    skip_flag;      /* skip_flag */
    uint32_t    imagex;         /* Image Buffer Size */
    uint32_t    imagey;
    uint32_t    imaged;         /* Image depth */
    void        *chdr;          /* Color Map Header = dummy */
    uint32_t    map_flag;       /* remap image? */
    uint32_t    *map;           /* map to use */
    uint32_t    xs,ys;          /* pos of changed area */
    uint32_t    xe,ye;          /* size of change area */
    uint32_t    special;        /* Special Info */
    void        *extra;         /* Decompression specific info */
} XA_DEC_INFO;

typedef struct
{
    uint8_t     *Ybuf;
    uint8_t     *Ubuf;
    uint8_t     *Vbuf;
    uint8_t     *the_buf;
    uint32_t    the_buf_size;
    uint16_t    y_w,y_h;
    uint32_t    uv_w,uv_h;
} YUVBufs;

typedef struct {
    uint32_t    Uskip_mask;
    int32_t     *YUV_Y_tab;
    int32_t     *YUV_UB_tab;
    int32_t     *YUV_VR_tab;
    int32_t     *YUV_UG_tab;
    int32_t     *YUV_VG_tab;
} YUVTabs;

uint8_t *y_buf;
uint8_t *u_buf;
uint8_t *v_buf;

void my_conv_func(uint8_t *image, uint32_t imagex, uint32_t imagey, uint32_t i_x, uint32_t i_y,
                  YUVBufs *yuv_bufs, YUVTabs *yuv_tabs, uint32_t map_flag, uint32_t map, void *chdr);
void my_conv_func(uint8_t *image, uint32_t imagex, uint32_t imagey, uint32_t i_x, uint32_t i_y,
                  YUVBufs *yuv_bufs, YUVTabs *yuv_tabs, uint32_t map_flag, uint32_t map, void *chdr)
{
    y_buf = yuv_bufs->Ybuf;
    u_buf = yuv_bufs->Ubuf;
    v_buf = yuv_bufs->Vbuf;
}

void *XA_YUV1611_Func(int);
void *XA_YUV1611_Func(int a)
{
    return my_conv_func;
}

void XA_Add_Func_To_Free_Chain(void *, void (*function)());
void XA_Add_Func_To_Free_Chain(void *anim_hdr, void (*function)())
{
}

void XA_Gen_YUV_Tabs(void *anim_hdr);
void XA_Gen_YUV_Tabs(void *anim_hdr)
{
}

void JPG_Setup_Samp_Limit_Table(void *anim_hdr);
void JPG_Setup_Samp_Limit_Table(void *anim_hdr)
{
}

void XA_Print(char *fmt, ...);
void XA_Print(char *fmt, ...)
{
    va_list vallist;
    va_start(vallist, fmt);
    //av_log(NULL, AV_LOG_ERROR, "RefLib error:\n");
    vfprintf(stderr,fmt,vallist);
    va_end(vallist);
}

YUVTabs     def_yuv_tabs;
void        *plug_hndl = 0;
void        *what_the  = 0;

static int load_plugin (AVCodecContext *avctx)
{
   /* void    handle, *fptr;*/
    char    *err_msg;
    int     err;
    int     (*start)(uint32_t, uint32_t, uint32_t);

    plug_hndl = dlopen(INDEO41_DLL, RTLD_LAZY);
    if (!plug_hndl) {
        av_log(avctx, AV_LOG_ERROR, "Unable to open reference library!\n");
        err_msg = dlerror();
        if (err_msg != NULL)
            av_log(avctx, AV_LOG_ERROR, "Error: (%s)\n", err_msg);
        return -1;
    }

    what_the = dlsym(plug_hndl, "What_The");
    if (!what_the) {
        av_log(avctx, AV_LOG_ERROR, "Unable to locate What_The!\n");
        dlclose(plug_hndl);
        plug_hndl = 0;
        return -1;
    }

    /* initialize reference decoder */
    start = (int (*)(uint32_t, uint32_t, uint32_t))((int)what_the - 0x520);
    err = (*start)(avctx->width, avctx->height, 1);
    if (err) {
        av_log(avctx, AV_LOG_ERROR, "Unable to initialize ref decoder: %d!\n", err);
        dlclose(plug_hndl);
        plug_hndl = 0;
        return -1;
    }

    return 0;
}

static int plugin_decode (AVCodecContext *avctx, const uint8_t *buf, int data_size)
{
    XA_DEC_INFO info;
    uint32_t    err;
    uint32_t    (*decode)(uint8_t *, uint8_t *, uint32_t, XA_DEC_INFO *);

    info.imagex    = avctx->width;
    info.imagey    = avctx->height;
    info.skip_flag = 0;

    decode = (uint32_t (*)(uint8_t *, uint8_t *, uint32_t, XA_DEC_INFO *))((int)what_the - 0x2A0);
    err = (*decode)(NULL, buf, data_size, &info);
    if (err) {
        av_log(avctx, AV_LOG_ERROR, "Error during decoding: %d!\n", err);
        return -1;
    }

    return 0;
}

static void unload_plugin ()
{
    void    (*end)(void);

    /* close reference decoder */
    end = (void (*)())((int)what_the + 0x140);
    (*end)();

    if (plug_hndl) {
        dlclose(plug_hndl);
        plug_hndl = 0;
    }
}

#endif


/**
 *  Decodes subdivision of a plane.
 *  This is a simplified version that checks for two supported subdivisions:
 *  - 1 wavelet band  per plane, size factor 1:1, code pattern: 3
 *  - 4 wavelet bands per plane, size factor 1:4, code pattern: 2,3,3,3,3
 *  Anything else is either unsupported or corrupt.
 *
 *  @param gb   [in,out] the GetBit context
 *  @return     number of wavelet bands or 0 if error
 */
static int decode_plane_subdivision(GetBitContext *gb)
{
    int i;

    switch (get_bits(gb, 2)) {
    case 3:
        return 1;
    case 2:
        for (i = 4; i && get_bits(gb, 2) == 3; i--);
        return (i ? 0 : 4);
    default:
        return 0;
    }
}

static inline int DEC_TILE_SIZE(int size_factor, int def_size)
{
    return (size_factor == 15 ? def_size : (size_factor + 1) << 5);
}

/**
 *  Decodes Indeo4 picture header.
 *
 *  @param ctx      [in,out] ptr to the decoder context
 *  @param avctx    [in] ptr to the AVCodecContext
 *  @return         result code: 0 = OK, -1 = error
 */
static int decode_pic_hdr(IVI4DecContext *ctx, AVCodecContext *avctx)
{
    int             pic_size_indx, val, i, p;
    IVIPicConfig    pic_conf;

    if (get_bits(&ctx->gb, 18) != 0x3FFF8) {
        av_log(avctx, AV_LOG_ERROR, "Invalid picture start code!\n");
        return -1;
    }

    ctx->prev_frame_type = ctx->frame_type;
    ctx->frame_type      = get_bits(&ctx->gb, 3);
    if (ctx->frame_type == 7) {
        av_log(avctx, AV_LOG_ERROR, "Invalid frame type: %d \n", ctx->frame_type);
        return -1;
    }

#if IVI4_STREAM_ANALYSER
    if (ctx->frame_type == 1 || ctx->frame_type == 3)
        has_b_frames = 1;
#endif

    ctx->transp_status = get_bits1(&ctx->gb);
    if (ctx->transp_status) {
        //av_log(avctx, AV_LOG_ERROR, "Transparency mode is enabled!\n");
#if IVI4_STREAM_ANALYSER
        has_transp = 1;
#endif
    }

    /* Unknown bit: Mac decoder ignores this bit, XANIM one returns error */
    if (get_bits1(&ctx->gb)) {
        av_log(avctx, AV_LOG_ERROR, "Sync bit is set!\n");
        return -1;
    }

    ctx->data_size = (get_bits1(&ctx->gb)) ? get_bits_long(&ctx->gb, 24) : 0;

    /* Null frames don't contain anything else so we just return */
    if (ctx->frame_type >= FRAMETYPE_NULL) {
        //av_log(avctx, AV_LOG_ERROR, "Null frame encountered!\n");
        return 0;
    }

    /* Check key lock status. If enabled - ignore lock word.   */
    /* Usually we have to promt the user for the password, but */
    /* we don't do because indeo4 videos can be decoded anyway */
    if (get_bits1(&ctx->gb)) {
        skip_bits_long(&ctx->gb, 32);
        av_log(avctx, AV_LOG_ERROR, "Password-protected clip!\n");
    }

    pic_size_indx = get_bits(&ctx->gb, 3);
    if (pic_size_indx == IVI4_PIC_SIZE_ESC) {
        pic_conf.pic_height = get_bits(&ctx->gb, 16);
        pic_conf.pic_width  = get_bits(&ctx->gb, 16);
    } else {
        pic_conf.pic_height = ivi4_common_pic_sizes[pic_size_indx * 2 + 1];
        pic_conf.pic_width  = ivi4_common_pic_sizes[pic_size_indx * 2    ];
    }

    /* Decode tile dimensions. */
    if (get_bits1(&ctx->gb)) {
        pic_conf.tile_height = DEC_TILE_SIZE(get_bits(&ctx->gb, 4), pic_conf.pic_height);
        pic_conf.tile_width  = DEC_TILE_SIZE(get_bits(&ctx->gb, 4), pic_conf.pic_width );
#if IVI4_STREAM_ANALYSER
        uses_tiling = 1;
#endif
    } else {
        pic_conf.tile_height = pic_conf.pic_height;
        pic_conf.tile_width  = pic_conf.pic_width;
    }

    /* Decode chroma subsampling. We support only 4:4 aka YVU9. */
    if (get_bits(&ctx->gb, 2)) {
        av_log(avctx, AV_LOG_ERROR, "Only YVU9 picture format is supported!\n");
        return -1;
    }
    pic_conf.chroma_height = (pic_conf.pic_height + 3) >> 2;
    pic_conf.chroma_width  = (pic_conf.pic_width  + 3) >> 2;

    /* decode subdivision of the planes */
    pic_conf.luma_bands = decode_plane_subdivision(&ctx->gb);
    if (pic_conf.luma_bands)
        pic_conf.chroma_bands = decode_plane_subdivision(&ctx->gb);
    ctx->is_scalable = pic_conf.luma_bands != 1 || pic_conf.chroma_bands != 1;
    if (ctx->is_scalable && (pic_conf.luma_bands != 4 || pic_conf.chroma_bands != 1)) {
        av_log(avctx, AV_LOG_ERROR, "Scalability: unsupported subdivision! Luma bands: %d, chroma bands: %d\n",
               pic_conf.luma_bands, pic_conf.chroma_bands);
        return -1;
    }

    /* check if picture layout was changed and reallocate buffers */
    if (ivi_pic_config_cmp(&pic_conf, &ctx->pic_conf)) {
        if (ff_ivi_init_planes(ctx->planes, &pic_conf)) {
            av_log(avctx, AV_LOG_ERROR, "Couldn't reallocate color planes!\n");
            return -1;
        }

        ctx->pic_conf = pic_conf;

        /* set default macroblock/block dimensions */
        for (p = 0; p <= 2; p++) {
            for (i = 0; i < (!p ? pic_conf.luma_bands : pic_conf.chroma_bands); i++) {
                ctx->planes[p].bands[i].mb_size  = !p ? (!ctx->is_scalable ? 16 : 8) : 4;
                ctx->planes[p].bands[i].blk_size = !p ? 8 : 4;
            }
        }

        if (ff_ivi_init_tiles(ctx->planes, ctx->pic_conf.tile_width,
                                   ctx->pic_conf.tile_height)) {
            av_log(avctx, AV_LOG_ERROR,
                   "Couldn't reallocate internal structures!\n");
            return -1;
        }
    }

    ctx->frame_num = (get_bits1(&ctx->gb)) ? get_bits_long(&ctx->gb, 20) : 0;

    /* skip decTimeEst field if present */
    if (get_bits1(&ctx->gb)) skip_bits(&ctx->gb, 8);

    /* decode macroblock and block huffman codebooks */
    if (ff_ivi_dec_huff_desc(&ctx->gb, get_bits1(&ctx->gb), IVI_MB_HUFF,  &ctx->mb_vlc,  avctx) ||
        ff_ivi_dec_huff_desc(&ctx->gb, get_bits1(&ctx->gb), IVI_BLK_HUFF, &ctx->blk_vlc, avctx))
        return -1;

    ctx->rvmap_sel = (get_bits1(&ctx->gb)) ? get_bits(&ctx->gb, 3) : 8;

    ctx->in_imf = get_bits1(&ctx->gb);
    ctx->in_q   = get_bits1(&ctx->gb);

    ctx->pic_glob_quant = get_bits(&ctx->gb, 5);

    /* TODO: ignore this parameter if unused */
    ctx->unknown1 = (get_bits1(&ctx->gb)) ? get_bits(&ctx->gb, 3) : 0;

    ctx->checksum = (get_bits1(&ctx->gb)) ? get_bits(&ctx->gb, 16) : 0;

    /* skip picture header extension if any */
    while (get_bits1(&ctx->gb)) {
        //av_log(avctx, AV_LOG_ERROR, "Pic hdr extension encountered!\n");
        val = get_bits(&ctx->gb, 8);
    }

    if (get_bits1(&ctx->gb)) {
        av_log(avctx, AV_LOG_ERROR, "Bad blocks bits encountered!\n");
    }

    align_get_bits(&ctx->gb);

    return 0;
}


/**
 *  Decodes Indeo4 band header.
 *
 *  @param ctx      [in,out] ptr to the decoder context
 *  @param band     [in,out] ptr to the band descriptor
 *  @param avctx    [in] ptr to the AVCodecContext
 *  @return         result code: 0 = OK, -1 = error
 */
static int decode_band_hdr(IVI4DecContext *ctx, IVIBandDesc *band,
                           AVCodecContext *avctx)
{
    int plane, band_num, hdr_size, indx, transform_id, scan_indx;
    int i;

    plane    = get_bits(&ctx->gb, 2);
    band_num = get_bits(&ctx->gb, 4);
    if (band->plane != plane || band->band_num != band_num) {
        av_log(avctx, AV_LOG_ERROR, "Invalid band header sequence!\n");
        return -1;
    }

    band->is_empty = get_bits1(&ctx->gb);
    if (!band->is_empty) {
        hdr_size = get_bits1(&ctx->gb) ? get_bits(&ctx->gb, 16) : 4;

        band->is_halfpel = get_bits(&ctx->gb, 2);
        if (band->is_halfpel >= 2) {
            av_log(avctx, AV_LOG_ERROR, "Invalid/unsupported mv resolution: %d!\n",
                   band->is_halfpel);
            return -1;
        }
#if IVI4_STREAM_ANALYSER
        if (!band->is_halfpel)
            uses_fullpel = 1;
#endif

        band->checksum_present = get_bits1(&ctx->gb);
        if (band->checksum_present)
            band->checksum = get_bits(&ctx->gb, 16);

        indx = get_bits(&ctx->gb, 2);
        if (indx == 3) {
            av_log(avctx, AV_LOG_ERROR, "Invalid block size!\n");
            return -1;
        }
        band->mb_size  = 16 >> indx;
        band->blk_size = 8 >> (indx >> 1);

        band->inherit_mv     = get_bits1(&ctx->gb);
        band->inherit_qdelta = get_bits1(&ctx->gb);

        band->glob_quant = get_bits(&ctx->gb, 5);

        if (!get_bits1(&ctx->gb) || ctx->frame_type == FRAMETYPE_INTRA) {
            transform_id = get_bits(&ctx->gb, 5);
            if ((transform_id >= 7 && transform_id <= 9) ||
                 transform_id == 17) {
                av_log(avctx, AV_LOG_ERROR, "DCT transform not supported yet!\n");
                return -1;
            }

#if IVI4_STREAM_ANALYSER
            if ((transform_id >= 0 && transform_id <= 2) || transform_id == 10)
                uses_haar = 1;
#endif

            band->inv_transform = transforms[transform_id].inv_trans;
            band->dc_transform  = transforms[transform_id].dc_trans;
            band->is_2d_trans   = transforms[transform_id].is_2d_trans;

            scan_indx = get_bits(&ctx->gb, 4);
            if (scan_indx == 15) {
                av_log(avctx, AV_LOG_ERROR, "Custom scan pattern encountered!\n");
                return -1;
            }
            band->scan = scan_index_to_tab[scan_indx];

            band->quant_mat = get_bits(&ctx->gb, 5);
            if (band->quant_mat == 31) {
                av_log(avctx, AV_LOG_ERROR, "Custom quant matrix encountered!\n");
                return -1;
            }
        }

        /* decode block huffman codebook */
        if (ff_ivi_dec_huff_desc(&ctx->gb, get_bits1(&ctx->gb), IVI_BLK_HUFF, &band->blk_vlc, avctx))
            return -1;

        /* select appropriate rvmap table for this band */
        band->rvmap_sel = (get_bits1(&ctx->gb)) ? get_bits(&ctx->gb, 3) : 8;

        /* decode rvmap probability corrections if any */
        band->num_corr = 0; /* there is no corrections */
        if (get_bits1(&ctx->gb)) {
            band->num_corr = get_bits(&ctx->gb, 8); /* get number of correction pairs */
            if (band->num_corr > 61) {
                av_log(avctx, AV_LOG_ERROR, "Too many corrections: %d\n",
                       band->num_corr);
                return -1;
            }

            /* read correction pairs */
            for (i = 0; i < band->num_corr * 2; i++)
                band->corr[i] = get_bits(&ctx->gb, 8);
        }
    }

    if (band->blk_size == 8) {
        band->intra_base1 = &ivi4_quant_8x8_intra[quant_index_to_tab[band->quant_mat]][0];
        band->inter_base1 = &ivi4_quant_8x8_inter[quant_index_to_tab[band->quant_mat]][0];
    } else {
        band->intra_base1 = &ivi4_quant_4x4_intra[quant_index_to_tab[band->quant_mat]][0];
        band->inter_base1 = &ivi4_quant_4x4_inter[quant_index_to_tab[band->quant_mat]][0];
    }

    /* indeo4 doesn't use scale tables so we set them to the dummy one */
    band->intra_scale1 = dummy_scale_tab;
    band->inter_scale1 = dummy_scale_tab;

    align_get_bits(&ctx->gb);

    return 0;
}


/**
 *  Decodes info (block type, cbp, quant delta, motion vector)
 *  for all macroblocks in the current tile.
 *
 *  @param ctx      [in,out] ptr to the decoder context
 *  @param band     [in,out] ptr to the band descriptor
 *  @param tile     [in,out] ptr to the tile descriptor
 *  @param avctx    [in] ptr to the AVCodecContext
 *  @return         result code: 0 = OK, -1 = error
 */
static int decode_mb_info(IVI4DecContext *ctx, IVIBandDesc *band,
                          IVITile *tile, AVCodecContext *avctx)
{
    int         x, y, offs, mb_offset, blks_per_mb, mb_type_bits;
    IVIMbInfo   *mb, *ref_mb;
    int         row_offset = band->mb_size * band->pitch;

    mb          = tile->mbs;
    ref_mb      = tile->ref_mbs;
    offs        = tile->ypos * band->pitch + tile->xpos;

    blks_per_mb  = band->mb_size != band->blk_size ? 4 : 1;
    mb_type_bits = ctx->frame_type == FRAMETYPE_BIDIR ? 2 : 1;

    for (y = tile->ypos; y < (tile->ypos + tile->height); y += band->mb_size) {
        mb_offset = offs;

        for (x = tile->xpos; x < (tile->xpos + tile->width); x += band->mb_size) {
            mb->xpos     = x;
            mb->ypos     = y;
            mb->buf_offs = mb_offset;

            if (get_bits1(&ctx->gb)) {
                if (ctx->frame_type == FRAMETYPE_INTRA) {
                    av_log(avctx, AV_LOG_ERROR, "Empty macroblock in an INTRA picture!\n");
                    return -1;
                }
                mb->type = 1; /* empty macroblocks are always INTER */
                mb->cbp  = 0; /* all blocks are empty */
            } else {
                if (band->inherit_mv) {
                    mb->type = ref_mb->type; /* copy mb_type from corresponding reference mb */
                } else if (ctx->frame_type == FRAMETYPE_INTRA) {
                    mb->type = 0; /* mb_type is always INTRA for intra-frames */
                } else {
                    mb->type = get_bits(&ctx->gb, mb_type_bits);
                }

                mb->cbp = get_bits(&ctx->gb, blks_per_mb);

                mb->q_delta = 0;
                if (band->inherit_qdelta) {
                    if (ref_mb) mb->q_delta = ref_mb->q_delta;
                } else if (mb->cbp || (!band->plane && !band->band_num &&
                           ctx->in_q)) {
                    mb->q_delta = get_vlc2(&ctx->gb, ctx->mb_vlc.tab->table,
                                           IVI_VLC_BITS, 1);
                    mb->q_delta = IVI_TOSIGNED(mb->q_delta);
                }

                if (!mb->type) {
                    mb->mv_x = mb->mv_y = 0; /* there is no motion vector in intra-macroblocks */
                } else {
                }
            }

            mb++;
            if (ref_mb)
                ref_mb++;
            mb_offset += band->mb_size;
        }

        offs += row_offset;
    }
    align_get_bits(&ctx->gb);

    return 0;
}


/**
 *  Decodes an Indeo4 band.
 *
 *  @param ctx      [in,out] ptr to the decoder context
 *  @param band     [in,out] ptr to the band descriptor
 *  @param avctx    [in] ptr to the AVCodecContext
 *  @return         result code: 0 = OK, -1 = error
 */
static int decode_band(IVI4DecContext *ctx, int plane_num,
                       IVIBandDesc *band, AVCodecContext *avctx)
{
    int         result, i, t, pos, idx1, idx2;
    IVITile     *tile;

    band->buf     = band->bufs[0];//[ctx->dst_buf];
    band->ref_buf = band->bufs[1];//[ctx->ref_buf];

    result = decode_band_hdr(ctx, band, avctx);
    if (result) {
        av_log(avctx, AV_LOG_ERROR, "Error while decoding band header: %d\n",
               result);
        return -1;
    }

    if (band->is_empty) {
        av_log(avctx, AV_LOG_ERROR, "Empty band encountered!\n");
        return -1;
    }

    band->rv_map = &ctx->rvmap_tabs[band->rvmap_sel];

    /* apply corrections to the selected rvmap table if present */
    for (i = 0; i < band->num_corr; i++) {
        idx1 = band->corr[i*2];
        idx2 = band->corr[i*2+1];
        FFSWAP(uint8_t, band->rv_map->runtab[idx1], band->rv_map->runtab[idx2]);
        FFSWAP(int16_t, band->rv_map->valtab[idx1], band->rv_map->valtab[idx2]);
    }

    pos = get_bits_count(&ctx->gb);

    for (t = 0; t < band->num_tiles; t++) {
        tile = &band->tiles[t];

        tile->is_empty = get_bits1(&ctx->gb);
        if (tile->is_empty) {
            //ff_ivi_process_empty_tile(avctx, band, tile,
            //                          (ctx->planes[0].bands[0].mb_size >> 3) - (band->mb_size >> 3));
            av_log(avctx, AV_LOG_ERROR, "Empty tile encountered!\n");
        } else {
            tile->data_size = ff_ivi_dec_tile_data_size(&ctx->gb);
            if (!tile->data_size) {
                av_log(avctx, AV_LOG_ERROR, "Tile data size = NULL!\n");
                return -1;
            }

            result = decode_mb_info(ctx, band, tile, avctx);
            if (result < 0)
                break;

            result = ff_ivi_decode_blocks(&ctx->gb, band, tile);
            if (result < 0 || (get_bits_count(&ctx->gb) - pos) >> 3 != tile->data_size) {
                av_log(avctx, AV_LOG_ERROR, "Corrupted tile data encountered!\n");
                break;
            }

            pos += tile->data_size << 3; // skip to next tile
            //skip_bits_long(&ctx->gb, pos - get_bits_count(&ctx->gb));
        }
    }

    /* restore the selected rvmap table by applying its corrections in reverse order */
    for (i = band->num_corr-1; i >= 0; i--) {
        idx1 = band->corr[i*2];
        idx2 = band->corr[i*2+1];
        FFSWAP(uint8_t, band->rv_map->runtab[idx1], band->rv_map->runtab[idx2]);
        FFSWAP(int16_t, band->rv_map->valtab[idx1], band->rv_map->valtab[idx2]);
    }

#if IVI_DEBUG
    if (band->checksum_present) {
        uint16_t chksum = ivi_calc_band_checksum(band);
        if (chksum != band->checksum) {
            av_log(avctx, AV_LOG_ERROR,
                   "Band checksum mismatch! Plane %d, band %d, received: %x, calculated: %x\n",
                   band->plane, band->band_num, band->checksum, chksum);
        }
    }
#endif

    align_get_bits(&ctx->gb);

    return 0;
}


/**
 *  Initializes Indeo4 decoder.
 */
static av_cold int decode_init(AVCodecContext *avctx)
{
    IVI4DecContext  *ctx = avctx->priv_data;

    ff_ivi_init_static_vlc();

    load_plugin(avctx);

    /* copy rvmap tables in our context so we can apply changes to them */
    memcpy(ctx->rvmap_tabs, ff_ivi_rvmap_tabs, sizeof(ff_ivi_rvmap_tabs));

    /* Force allocation of the internal buffers */
    /* during picture header decoding.          */
    ctx->pic_conf.pic_width     = 0;
    ctx->pic_conf.pic_height    = 0;

    avctx->pix_fmt = PIX_FMT_YUV410P;

    return 0;
}


/**
 *  main decoder function
 */
static int decode_frame(AVCodecContext *avctx, void *data, int *data_size,
                        AVPacket *avpkt)
{
    IVI4DecContext  *ctx = avctx->priv_data;
    const uint8_t   *buf = avpkt->data;
    int             buf_size = avpkt->size;
    int             result, p, b;

    init_get_bits(&ctx->gb, buf, buf_size * 8);

    //plugin_decode (avctx, avpkt->data, avpkt->size);

    result = decode_pic_hdr(ctx, avctx);
    if (result) {
        av_log(avctx, AV_LOG_ERROR,
               "Error while decoding picture header: %d\n", result);
        return -1;
    }

    //START_TIMER;

    if (ctx->frame_type >= FRAMETYPE_NULL) {
        ctx->frame_type = ctx->prev_frame_type;
    } else if (ctx->frame_type == FRAMETYPE_INTRA) {
        for (p = 0; p < 3; p++) {
            for (b = 0; b < ctx->planes[p].num_bands; b++) {
                result = decode_band(ctx, p, &ctx->planes[p].bands[b], avctx);
                if (result) {
                    av_log(avctx, AV_LOG_ERROR,
                           "Error while decoding band: %d, plane: %d\n", b, p);
                    return -1;
                }
            }
        }
    }

    //STOP_TIMER("decode_planes");

    if (ctx->frame.data[0])
        avctx->release_buffer(avctx, &ctx->frame);

    ctx->frame.reference = 0;
    if (avctx->get_buffer(avctx, &ctx->frame) < 0) {
        av_log(avctx, AV_LOG_ERROR, "get_buffer() failed\n");
        return -1;
    }

    ff_ivi_output_plane(&ctx->planes[0], ctx->frame.data[0], ctx->frame.linesize[0]);
    ff_ivi_output_plane(&ctx->planes[2], ctx->frame.data[1], ctx->frame.linesize[1]);
    ff_ivi_output_plane(&ctx->planes[1], ctx->frame.data[2], ctx->frame.linesize[2]);

    *data_size = sizeof(AVFrame);
    *(AVFrame*)data = ctx->frame;

    return buf_size;
}


/**
 *  Closes Indeo4 decoder and cleans up its context.
 */
static av_cold int decode_close(AVCodecContext *avctx)
{
    IVI4DecContext *ctx = avctx->priv_data;

    ff_ivi_free_buffers(&ctx->planes[0]);

    if (ctx->frame.data[0])
        avctx->release_buffer(avctx, &ctx->frame);

#if IVI4_STREAM_ANALYSER
    if (ctx->is_scalable)
        av_log(avctx, AV_LOG_ERROR, "This video uses scalability mode!\n");
    if (uses_tiling)
        av_log(avctx, AV_LOG_ERROR, "This video uses local decoding!\n");
    if (has_b_frames)
        av_log(avctx, AV_LOG_ERROR, "This video contains B-frames!\n");
    if (has_transp)
        av_log(avctx, AV_LOG_ERROR, "Transparency mode is enabled!\n");
    if (uses_haar)
        av_log(avctx, AV_LOG_ERROR, "This video uses Haar transform!\n");
    if (uses_fullpel)
        av_log(avctx, AV_LOG_ERROR, "This video uses fullpel motion vectors!\n");
#endif

    unload_plugin();

    return 0;
}


AVCodec indeo4_decoder = {
    .name           = "indeo4",
    .type           = CODEC_TYPE_VIDEO,
    .id             = CODEC_ID_INDEO4,
    .priv_data_size = sizeof(IVI4DecContext),
    .init           = decode_init,
    .close          = decode_close,
    .decode         = decode_frame,
    .long_name      = NULL_IF_CONFIG_SMALL("Intel Indeo Video Interactive 4"),
};
