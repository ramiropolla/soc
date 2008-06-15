/*
 * WMA 9/3/PRO compatible decoder
 * Copyright (c) 2007 The FFmpeg Project.
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
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

#include "avcodec.h"
#include "bitstream.h"


/* TODO
   shrink types
*/

typedef struct WMA3DecodeContext {
    AVCodecContext*     avctx;
    GetBitContext       gb;
    int                 buf_bit_size;

    // Packet info
    int                 packet_sequence_number;
    int                 bit5;   // might be that packet contains padding data, last packet in some cbr files have it
    int                 bit6;

    // Stream info
    unsigned int        samples_per_frame;
    unsigned int        log2_block_align;
    unsigned int        log2_block_align_bits;
    unsigned int        log2_frame_size;

    // Extradata
    unsigned int        decode_flags;
    unsigned int        dwChannelMask;
    unsigned int        sample_bit_depth;

    // Packet loss variables
    unsigned int        packet_loss;

    // General frame info
    unsigned int        frame_num;
    int                 len_prefix; //< true if the frame is prefixed with its len
    int                 allow_subframes;
    int                 max_num_subframes;

    // Buffered frame data
    int                 prev_frame_bit_size;
    uint8_t*            prev_frame;

} WMA3DecodeContext;


static void dump_context(WMA3DecodeContext *s)
{
#define PRINT(a,b) av_log(NULL,AV_LOG_ERROR," %s = %d\n", a, b);
#define PRINT_HEX(a,b) av_log(NULL,AV_LOG_ERROR," %s = %x\n", a, b);

    PRINT_HEX("ed channelmask",s->dwChannelMask);
    PRINT("ed sample bit depth",s->sample_bit_depth);
    PRINT_HEX("ed decode flags",s->decode_flags);
    PRINT("samples per frame",s->samples_per_frame);
    PRINT("log2 frame size",s->log2_frame_size);
    PRINT("max num subframes",s->max_num_subframes);
    PRINT("len prefix",s->len_prefix);
}


/**
 *  Get the samples per frame for this stream
 */
static int get_samples_per_frame(int sample_rate, unsigned int decode_flags) {

    int samples_per_frame;
    int tmp;

    if (sample_rate <= 16000)
        samples_per_frame = 512;
    else if (sample_rate <= 22050)
        samples_per_frame = 1024;
    else if (sample_rate <= 48000)
        samples_per_frame = 2048;
    else if (sample_rate <= 96000)
        samples_per_frame = 4096;
    else
        samples_per_frame = 8192;

 /* wma voice code  if (decode_flags & 0x800) {
        tmp = ((decode_flags & 6) >> 1) | ((decode_flags & 0x600) >> 7);
        samples_per_frame = (tmp+1)*160;
    } else { */

    tmp = decode_flags & 0x6;
    if (tmp == 0x2)
        samples_per_frame <<= 1;
    else if (tmp == 0x4)
        samples_per_frame >>= 1;
    else if (tmp == 0x6)
        samples_per_frame >>= 2;


    return samples_per_frame;
}

static av_cold int wma3_decode_init(AVCodecContext *avctx)
{
    WMA3DecodeContext *s = avctx->priv_data;
    uint8_t *edata_ptr = avctx->extradata;
    int i;
    int log2_max_num_subframes;

    s->avctx = avctx;

    if (avctx->extradata_size >= 18) {
        s->decode_flags     = AV_RL16(edata_ptr+14);
        s->dwChannelMask    = AV_RL32(edata_ptr+2);
        s->sample_bit_depth = AV_RL16(edata_ptr);

        /* Dump the extradata */
        for (i=0 ; i<avctx->extradata_size ; i++)
            av_log(avctx, AV_LOG_ERROR, "[%x] ",avctx->extradata[i]);
        av_log(avctx, AV_LOG_ERROR, "\n");

    } else {
        av_log(avctx, AV_LOG_ERROR, "Unknown extradata size %d.\n",avctx->extradata_size);
        return -1;
    }

    s->samples_per_frame = get_samples_per_frame(avctx->sample_rate, s->decode_flags);

    /* Generic init */
    s->packet_loss = 0;
    s->log2_block_align = av_log2(avctx->block_align);
    s->log2_block_align_bits = av_log2(avctx->block_align*8);
    s->log2_frame_size = s->log2_block_align_bits + 1;

    /* frame info */
    s->len_prefix = s->decode_flags & 0x40;

    if(!s->len_prefix){
         av_log(avctx, AV_LOG_ERROR, "file has no len prefix please report\n");
         return -1;
    }

    /* subframe info */
    log2_max_num_subframes = (s->decode_flags & 0x38) >> 3;
    s->max_num_subframes = 1 << log2_max_num_subframes;
    s->allow_subframes = s->max_num_subframes > 1;

    dump_context(s);

    return 0;
}



/* decode one wma frame */
static int wma_decode_frame(WMA3DecodeContext *s,GetBitContext* gb){
    unsigned int gb_start_count = get_bits_count(gb);
    int more_frames = 0;
    /* get frame length */
    int len = 0;

    if(s->len_prefix)
        len = get_bits(gb,s->log2_frame_size);

    av_log(s->avctx,AV_LOG_INFO,"decoding frame with len %x\n",len);

    /* skip the rest of the frame data */
    skip_bits_long(gb,len - (get_bits_count(gb) - gb_start_count) - 1);

    /* decode trailer bit */
    more_frames = get_bits1(gb);

    ++s->frame_num;
    return more_frames;
}


/* calculate remaining bits from the input buffer */
static int remaining_bits(WMA3DecodeContext *s){
    return s->buf_bit_size - get_bits_count(&s->gb);
}

/* store bits from the prev frame into a temporary buffer */
static void save_bits(WMA3DecodeContext *s,int len){
    int buflen = (s->prev_frame_bit_size + len + 8) / 8;
    int bit_offset = s->prev_frame_bit_size % 8;
    int pos = (s->prev_frame_bit_size - bit_offset) / 8;
    s->prev_frame_bit_size += len;

    if(len <= 0)
         return;

    s->prev_frame = av_realloc(s->prev_frame,buflen + FF_INPUT_BUFFER_PADDING_SIZE);

    /* byte align prev_frame buffer */
    if(bit_offset){
        int missing = 8 - bit_offset;
        if(len < missing)
            missing = len;
        s->prev_frame[pos++] |= get_bits(&s->gb, missing) << (8 - bit_offset - missing);
        len -= missing;
    }

    /* copy full bytes */
    while(len > 7){
        s->prev_frame[pos++] = get_bits(&s->gb,8);
        len -= 8;
    }

    /* copy remaining bits */
    if(len > 0)
        s->prev_frame[pos++] = get_bits(&s->gb,len) << (8 - len);
}

static int wma3_decode_packet(AVCodecContext *avctx,
                             void *data, int *data_size,
                             const uint8_t *buf, int buf_size)
{
    WMA3DecodeContext *s = avctx->priv_data;
    int more_frames=1;
    int num_bits_prev_frame;
    s->buf_bit_size = buf_size << 3;

    *data_size = 0;

    //FIXME check minimum buffer size and check for security problems!!!
    if(buf_size < avctx->block_align)
        return 0;

    /* Parse packet header */
    init_get_bits(&s->gb, buf, s->buf_bit_size);
    s->packet_sequence_number = get_bits(&s->gb, 4);
    s->bit5                   = get_bits1(&s->gb);
    s->bit6                   = get_bits1(&s->gb);

    /* get number of bits that need to be added to the previous frame */
    num_bits_prev_frame = get_bits(&s->gb, s->log2_frame_size);
    av_log(avctx, AV_LOG_INFO, "[%d]: nbpf %x\n", avctx->frame_number, num_bits_prev_frame);

    /* check for packet loss */
    if (s->packet_sequence_number != (avctx->frame_number&0xF)) {
        s->packet_loss = 1;
        av_log(avctx, AV_LOG_ERROR, "!!Packet loss detected! seq %x vs %x\n",s->packet_sequence_number,avctx->frame_number&0xF);
    }

    if (num_bits_prev_frame > 0) {
        /* append the prev frame data to the remaining data from the previous packet to create a full frame */
        save_bits(s,num_bits_prev_frame);
        av_log(avctx, AV_LOG_INFO, "accumulated %x bits of frame data\n",s->prev_frame_bit_size);

        /* decode the cross packet frame if it is valid */
        if(!s->packet_loss){
            GetBitContext gb_prev;
            init_get_bits(&gb_prev, s->prev_frame, s->prev_frame_bit_size);
            wma_decode_frame(s,&gb_prev);
        }

        /* reset prev frame buffer */
        s->prev_frame_bit_size = 0;
    }

    /* decode the rest of the packet */
    while(more_frames && remaining_bits(s) > s->log2_frame_size){
        int frame_size = show_bits(&s->gb, s->log2_frame_size);

        /* there is enough data for a full frame */
        if(remaining_bits(s) >= frame_size){
            /* decode the frame */
            more_frames = wma_decode_frame(s,&s->gb);

            if(!more_frames){
                av_log(avctx, AV_LOG_ERROR, "no more frames\n");
            }
        }else
            more_frames = 0;
    }

    /* save the rest of the data so that it can be decoded with the next packet */
    save_bits(s,remaining_bits(s));

    return avctx->block_align;
}

static av_cold int wma3_decode_end(AVCodecContext *avctx)
{
    WMA3DecodeContext *s = avctx->priv_data;
    if(s->prev_frame)
        av_free(s->prev_frame);
    return 0;
}

AVCodec wmav3pro_decoder =
{
    "wmav3Pro",
    CODEC_TYPE_AUDIO,
    CODEC_ID_WMAPRO,
    sizeof(WMA3DecodeContext),
    wma3_decode_init,
    NULL,
    wma3_decode_end,
    wma3_decode_packet,
    .long_name = "Windows Media Audio 9 Professional",
};
