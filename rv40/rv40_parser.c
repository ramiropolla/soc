/*
 * RV40 parser
 * Copyright (c) 2007 Konstantin Shishkov
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
 * @file rv40_parser.c
 * RV40 parser that combines slices
 */

#include "parser.h"

/** These bits of slice header should be the same for all slices in one frame */
#define KEY_MASK 0xE0CFFF80

typedef struct RV40ParseContext{
    ParseContext pc;
    int slices;
}RV40ParseContext;

static int rv40_parse_init(AVCodecParserContext *s)
{
    ParseContext *pc  = s->priv_data;

    pc->state = -1;
    return 0;
}

static int rv40_parse(AVCodecParserContext *s,
                      AVCodecContext *avctx,
                      const uint8_t **poutbuf, int *poutbuf_size,
                      const uint8_t *buf, int buf_size)
{
    ParseContext *pc  = s->priv_data;
    RV40ParseContext *rpc  = s->priv_data;
    int hdr, buf_size2;

    hdr = AV_RB32(buf) & KEY_MASK;
    if(pc->state == -1){
        rpc->slices = 0;
        pc->state = hdr;
    }
    if(hdr == pc->state){//add another slice
        rpc->slices++;
        buf_size2 = buf_size;
        /* adjust buffer size for faster slice header searches / bits reading */
        buf_size = (buf_size + 3) & ~3;
        if (ff_combine_frame(pc, END_NOT_FOUND, &buf, &buf_size) < 0) {
            *poutbuf = NULL;
            *poutbuf_size = 0;
            return buf_size2;
        }
        *poutbuf = NULL;
        *poutbuf_size = 0;
        return buf_size;
    }else{ // new frame
        pc->state = -1;
        if (ff_combine_frame(pc, 0, &buf, &buf_size) < 0) {
            *poutbuf = NULL;
            *poutbuf_size = 0;
            return buf_size;
        }
        *poutbuf = buf;
        *poutbuf_size = buf_size;
        return 0;
    }
}

AVCodecParser rv40_parser = {
    { CODEC_ID_RV40 },
    sizeof(RV40ParseContext),
    rv40_parse_init,
    rv40_parse,
    ff_parse_close,
};
