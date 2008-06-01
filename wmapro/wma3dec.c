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

typedef struct WMA3DecodeContext {
    AVCodecContext*     avctx;
    GetBitContext       gb;

    // Packet info
    int                 packet_sequence_number;
    int                 bit5;   // might be that packet contains padding data, last packet in some cbr files have it
    int                 bit6;

    // Stream info
    unsigned int        samples_per_frame;
    unsigned int        log2_block_align;
    unsigned int        log2_block_align_bits;
    unsigned int        log2_frame_size;
    unsigned int        num_bits_prev_frame;
    unsigned int        num_bits_curr_frame;
    unsigned int        cross_frame_needed_bits;
    unsigned int        cross_length_prev_bits;
    unsigned int        cross_length_prev;

    // Extradata
    unsigned int        decode_flags;
    unsigned int        dwChannelMask;
    unsigned int        sample_bit_depth;

    // Packet loss variables
    unsigned int        packet_loss;
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

static int wma3_decode_init(AVCodecContext *avctx)
{
    WMA3DecodeContext *s = avctx->priv_data;
    uint8_t *edata_ptr = avctx->extradata;
    int i;

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

    dump_context(s);

    return 0;
}

static int wma3_decode_packet(AVCodecContext *avctx,
                             void *data, int *data_size,
                             uint8_t *buf, int buf_size)
{
    WMA3DecodeContext *s = avctx->priv_data;
    int more_frames=1;
    int num_bits_curr_frame, num_bits_prev_frame;
    int i=0, sum=0, cross_frame = 0, buf_bit_size = buf_size << 3;

    /* Parse packet header */
    init_get_bits(&s->gb, buf, buf_bit_size);
    s->packet_sequence_number = get_bits(&s->gb, 4);
    s->bit5                   = get_bits1(&s->gb);
    s->bit6                   = get_bits1(&s->gb);

    /* Parse frames */
    num_bits_prev_frame = get_bits(&s->gb, s->log2_frame_size);
    av_log(avctx, AV_LOG_ERROR, "[%d]: nbpf %x\n", avctx->frame_number, num_bits_prev_frame);
    if (num_bits_prev_frame) {
        if (s->cross_length_prev_bits) {
            int num_bits_cross_frame = s->cross_length_prev | get_bits(&s->gb, s->log2_frame_size - s->cross_length_prev_bits);
            av_log(avctx, AV_LOG_ERROR, "[%d]: cross length field nbpf %x, %x\n", avctx->frame_number, num_bits_prev_frame + s->cross_length_prev_bits, num_bits_cross_frame);
            skip_bits_long(&s->gb, num_bits_cross_frame - s->log2_frame_size);
        }
        else
            skip_bits_long(&s->gb, num_bits_prev_frame);
        more_frames = 1;
    }
    s->cross_length_prev_bits = 0;
    s->cross_length_prev = 0;
    i = 0;
    do {
        num_bits_curr_frame = get_bits(&s->gb, s->log2_frame_size);
        av_log(avctx, AV_LOG_ERROR, "[%d]:%d nbcf %x ", avctx->frame_number,i, num_bits_curr_frame);
        if (num_bits_curr_frame >= s->log2_frame_size) {
            //Check if the next frame will need data from the next packet
            if (num_bits_curr_frame - s->log2_frame_size <= buf_bit_size - get_bits_count(&s->gb)) {
                int left_bits = buf_bit_size - get_bits_count(&s->gb) - (num_bits_curr_frame - s->log2_frame_size);
                if (left_bits == 0) {
                    skip_bits_long(&s->gb, num_bits_curr_frame - s->log2_frame_size);
                    more_frames = 0;
                    av_log(avctx, AV_LOG_ERROR, "\tfull end\n");
                }
                else {
                    skip_bits_long(&s->gb, num_bits_curr_frame - s->log2_frame_size - 1);
                    more_frames = get_bits1(&s->gb);
                    av_log(avctx, AV_LOG_ERROR, "\tlast bit %d \t lastbit align pos %d\n", more_frames,get_bits_count(&s->gb)%8);
                    // Need to save the bits to be handled in next packet?
                    if (more_frames && left_bits < s->log2_frame_size) {
                        s->cross_length_prev_bits = left_bits;
                        s->cross_length_prev = get_bits(&s->gb, left_bits);
                        s->cross_length_prev <<= s->log2_frame_size - left_bits;
                        more_frames = 0;
                        av_log(avctx, AV_LOG_ERROR, "cross length field found: %x/%x, %x\n", left_bits, s->log2_frame_size, s->cross_length_prev);
                    }
                }
            } else {
                more_frames = 0;
                cross_frame = 1;
                av_log(avctx, AV_LOG_ERROR, "\tlast bit: cross_frame %x \tleft_bit_size %x\n", num_bits_curr_frame, buf_bit_size - get_bits_count(&s->gb) + s->log2_frame_size);
                s->cross_frame_needed_bits = num_bits_curr_frame - s->log2_frame_size - (buf_bit_size - get_bits_count(&s->gb));
                av_log(avctx, AV_LOG_ERROR, "    Not available needed bits %x\n",s->cross_frame_needed_bits);
            }
        } else
            av_log(avctx, AV_LOG_ERROR, "\n");
        i++;
    } while(more_frames);


    if (!cross_frame) {
        av_log(avctx, AV_LOG_ERROR, "    Available bits %x - Consumed bits %x \t diff %x\n",buf_bit_size, get_bits_count(&s->gb),buf_bit_size-get_bits_count(&s->gb));
        s->cross_frame_needed_bits = buf_bit_size-get_bits_count(&s->gb);
    }

    //Check amount of non zero bits in non crossing frames
    if (!cross_frame) {
        for (i=0 ; i<buf_bit_size-get_bits_count(&s->gb) ; i++) {
            sum+=get_bits1(&s->gb);
        }
        if (sum)
            av_log(avctx, AV_LOG_ERROR, "!!Non crossing frame contains %x non zero bits!\n",sum);
    }

    if (s->packet_sequence_number != (avctx->frame_number&0xF)) {
        s->packet_loss = 1;
        av_log(avctx, AV_LOG_ERROR, "!!Packet loss detected! seq %x vs %x\n",s->packet_sequence_number,avctx->frame_number&0xF);
    }

    return avctx->block_align;
}

static int wma3_decode_end(AVCodecContext *avctx)
{
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
