/*
 * WMA 9/3/PRO compatible decoder
 * Copyright (c) 2007 Baptiste Coudurier, Benjamin Larsson, Ulion
 * Copyright (c) 2008 Sascha Sommer
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

#define MAX_CHANNELS 6               //< max number of handled channels
#define MAX_SUBFRAMES 32             //< max number of subframes per channel

/**
 *@brief decoder context for a single channel
 */
typedef struct {
    uint8_t num_subframes;                //< number of subframes for the current channel
    uint16_t subframe_len[MAX_SUBFRAMES]; //< subframe len in samples
    uint16_t channel_len;                 //< channel len in samples
} wma_channel_t;


/**
 *@brief main decoder context
 */
typedef struct WMA3DecodeContext {
    AVCodecContext*  avctx;                    //< codec context for av_log
    GetBitContext    gb;                       //< getbitcontext for the packet
    int              buf_bit_size;             //< buffer size in bits

    /** Packet info */
    uint8_t          packet_sequence_number;   //< current packet nr
    uint8_t          bit5;                     //< padding bit? (cbr files)
    uint8_t          bit6;
    uint8_t          packet_loss;              //< set in case of bitstream error

    /** Stream info */
    uint16_t         samples_per_frame;        //< nr of outputed samples
    uint16_t         log2_frame_size;          //< frame size
    uint8_t          lossless;                 //< lossless mode
    uint8_t          no_tiling;                //< frames are split in subframes
    int8_t           nb_channels;              //< nr of channels
    wma_channel_t    channel[MAX_CHANNELS];    //< per channel data

    /** Extradata */
    unsigned int     decode_flags;             //< used compression features
    unsigned int     dwChannelMask;
    uint8_t          sample_bit_depth;         //< bits per sample

    /** General frame info */
    unsigned int     frame_num;                //< current frame number
    uint8_t          len_prefix;               //< frame is prefixed with its len
    uint8_t          allow_subframes;          //< frames may contain subframes
    uint8_t          max_num_subframes;        //< maximum number of subframes
    uint16_t         min_samples_per_subframe; //< minimum samples per subframe
    uint8_t          dynamic_range_compression;//< frame contains drc data
    uint8_t          drc_gain;                 //< gain for the drc tool
    uint8_t          update_samples_per_frame; //< recalculate output size


    /** Buffered frame data */
    int              prev_frame_bit_size;      //< saved number of bits
    uint8_t*         prev_frame;               //< prev frame data

} WMA3DecodeContext;

/**
 *@brief helper function to print the most important members of the context
 *@param s context
 */
static void dump_context(WMA3DecodeContext *s)
{
#define PRINT(a,b) av_log(s->avctx,AV_LOG_ERROR," %s = %d\n", a, b);
#define PRINT_HEX(a,b) av_log(s->avctx,AV_LOG_ERROR," %s = %x\n", a, b);

    PRINT_HEX("ed channelmask",s->dwChannelMask);
    PRINT("ed sample bit depth",s->sample_bit_depth);
    PRINT_HEX("ed decode flags",s->decode_flags);
    PRINT("samples per frame",s->samples_per_frame);
    PRINT("log2 frame size",s->log2_frame_size);
    PRINT("max num subframes",s->max_num_subframes);
    PRINT("len prefix",s->len_prefix);
    PRINT("nb channels",s->nb_channels);
    PRINT("lossless",s->lossless);
}


/**
 *@brief Get the samples per frame for this stream
 *@param sample_rate output sample_rate
 *@param decode_flags codec compression features
 *@return number of output samples per frame
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

/**
 *@brief initialize the decoder
 *@param avctx codec context
 *@return 0 on success, -1 otherwise
 */
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

        /** Dump the extradata */
        for (i=0 ; i<avctx->extradata_size ; i++)
            av_log(avctx, AV_LOG_ERROR, "[%x] ",avctx->extradata[i]);
        av_log(avctx, AV_LOG_ERROR, "\n");

    } else {
        av_log(avctx, AV_LOG_ERROR, "Unknown extradata size %d.\n",
                      avctx->extradata_size);
        return -1;
    }

    s->samples_per_frame = get_samples_per_frame(avctx->sample_rate,
                                                 s->decode_flags);

    /** Generic init */
    s->packet_loss = 0;
    s->log2_frame_size = av_log2(avctx->block_align*8)+1;

    /** frame info */
    s->len_prefix = s->decode_flags & 0x40;

    if(!s->len_prefix){
         av_log(avctx, AV_LOG_ERROR, "file has no len prefix please report\n");
         return -1;
    }

    /** subframe info */
    log2_max_num_subframes = (s->decode_flags & 0x38) >> 3;
    s->max_num_subframes = 1 << log2_max_num_subframes;
    s->allow_subframes = s->max_num_subframes > 1;
    s->min_samples_per_subframe = s->samples_per_frame / s->max_num_subframes;
    s->dynamic_range_compression = s->decode_flags & 0x80;

    if(s->max_num_subframes > MAX_SUBFRAMES){
        av_log(avctx, AV_LOG_ERROR, "invalid number of subframes %i\n",
                      s->max_num_subframes);
        return -1;
    }

    s->nb_channels = avctx->channels;

    if(s->nb_channels < 0 || s->nb_channels > MAX_CHANNELS){
        av_log(avctx, AV_LOG_ERROR, "invalid number of channels %i\n",
                      s->nb_channels);
        return -1;
    }

    dump_context(s);

    return 0;
}



/**
 *@brief decode how the data in the frame is split into subframes
 *       every wma frame contains the encoded data for a fixed number of
 *       samples per channel the data for every channel might be split
 *       into several subframes this function will reconstruct the list of
 *       subframes for every channel
 *
 *       If the subframes are not evenly split the algorithm estimates the
 *       channels with the lowest number of total samples.
 *       Afterwards for every of these channels a bit is read from the
 *       bitstream that indicates if the channel contains a frame with the
 *       next subframesize that is going to be read from the bitstream or not.
 *       If a channel contains such a subframe the subframesize gets added to
 *       the channel's subframelist.
 *       The algorithm repeats these steps until the frame is properly divided
 *       between the individual channels.
 *
 *
 *@param s context
 *@param gb current get bit context
 *@return 0 on success < 0 in case of an error
 */
static int wma_decode_tilehdr(WMA3DecodeContext *s, GetBitContext* gb){
    int c;
    int missing_samples = s->nb_channels * s->samples_per_frame;

    /** reset tiling information */
    for(c=0;c<s->nb_channels;c++){
        s->channel[c].num_subframes = 0;
        s->channel[c].channel_len = 0;
    }

    /** handle the easy case whith one constant sized subframe per channel */
    if(s->max_num_subframes == 1){
        for(c=0;c<s->nb_channels;c++){
            s->channel[c].num_subframes = 1;
            s->channel[c].subframe_len[0] = s->samples_per_frame;
            s->channel[c].channel_len = 0;
        }
    }else{ /** subframe len and number of subframes is not constant */
        int subframe_len_bits = 0;     /** bits needed for the subframe len */
        int subframe_len_zero_bit = 0; /** how many of the len bits indicate
                                           if the subframe is zero */

        /** calculate subframe len bits */
        if(s->lossless)
            subframe_len_bits = av_log2(s->max_num_subframes - 1) + 1;
        else if(s->max_num_subframes == 16){
            subframe_len_zero_bit = 1;
            subframe_len_bits = 3;
        }else
            subframe_len_bits = av_log2(s->max_num_subframes) + 1;

        /** loop until the frame data is split between the subframes */
        while(missing_samples > 0){
             int64_t tileinfo = -1;
             int min_channel_len = s->samples_per_frame;
             int num_subframes_per_channel = 0;
             int num_channels = 0;
             int subframe_len = s->samples_per_frame / s->max_num_subframes;

             /** find channel with the smallest overall len */
             for(c=0;c<s->nb_channels;c++){
                 if(min_channel_len > s->channel[c].channel_len)
                     min_channel_len = s->channel[c].channel_len;
             }

             /** check if this is the start of a new frame */
             if(missing_samples == s->nb_channels * s->samples_per_frame){
                 s->no_tiling = get_bits1(gb);
             }

             if(s->no_tiling){
                 num_subframes_per_channel = 1;
                 num_channels = s->nb_channels;
             }else{
                  /** count how many channels have the minimum len */
                  for(c=0;c<s->nb_channels;c++){
                      if(min_channel_len == s->channel[c].channel_len){
                          ++num_channels;
                      }
                  }
                  if(num_channels <= 1)
                      num_subframes_per_channel = 1;
             }

             /** maximum number of subframes, evenly split */
             if(subframe_len == missing_samples / num_channels){
                 num_subframes_per_channel = 1;
             }

             /** if there might be multiple subframes per channel */
             if(num_subframes_per_channel != 1){
                 int total_num_bits = num_channels;
                 tileinfo = 0;
                 /** for every channel with the minimum len 1 bit is
                     transmitted that informs us if the channel
                     contains a subframe with the next subframe_len */
                 while(total_num_bits){
                    int num_bits = total_num_bits;
                    if(num_bits > 32)
                        num_bits = 32;
                    tileinfo |= get_bits_long(gb,num_bits);
                    total_num_bits -= num_bits;
                    num_bits = total_num_bits;
                    tileinfo <<= (num_bits > 32)? 32 : num_bits;
                 }
             }


             /** if the frames are not evenly split get the next subframe len
                 from the bitstream */
             if(subframe_len != missing_samples / num_channels){
                  int log2_subframe_len = 0;
                  /* 1 bit indicates if the subframe len is zero */
                  if(subframe_len_zero_bit){
                      if(get_bits1(gb)){
                          log2_subframe_len = get_bits(gb,subframe_len_bits-1);
                          ++log2_subframe_len;
                      }
                  }else
                      log2_subframe_len = get_bits(gb,subframe_len_bits);

                  if(s->lossless){
                      subframe_len =
                          s->samples_per_frame / s->max_num_subframes;
                      subframe_len *= log2_subframe_len + 1;
                  }else
                      subframe_len =
                          s->samples_per_frame / (1 << log2_subframe_len);

             }

             /** sanity check the len */
             if(subframe_len < s->min_samples_per_subframe
                  || subframe_len > s->samples_per_frame){
                  av_log(s->avctx, AV_LOG_ERROR,
                         "broken frame: subframe_len %i\n", subframe_len);
                  return -1;
             }
             for(c=0; c<s->nb_channels;c++){
                  wma_channel_t* chan = &s->channel[c];
                  if(chan->num_subframes > 32){
                      av_log(s->avctx, AV_LOG_ERROR,
                             "broken frame: num subframes %i\n",
                             chan->num_subframes);
                      return -1;
                  }

                  /** add subframes to the individual channels */
                  if(min_channel_len == chan->channel_len){
                       --num_channels;
                       if(tileinfo & (1<<num_channels)){
                            if(chan->num_subframes > 31){
                               av_log(s->avctx, AV_LOG_ERROR,
                                      "broken frame: num subframes > 31\n");
                               return -1;
                            }
                            chan->subframe_len[chan->num_subframes] = subframe_len;
                            chan->channel_len += subframe_len;
                            missing_samples -= subframe_len;
                            ++chan->num_subframes;
                            if(missing_samples < 0
                               || chan->channel_len > s->samples_per_frame){
                                av_log(s->avctx, AV_LOG_ERROR,"broken frame: "
                                       "channel len > samples_per_frame\n");
                                return -1;
                            }
                        }
                  }
             }
        }

    }

    for(c=0;c<s->nb_channels;c++){
        int i;
        for(i=0;i<s->channel[c].num_subframes;i++){
            av_log(s->avctx, AV_LOG_INFO,"frame[%i] channel[%i] subframe[%i]"
                   " len %i\n",s->frame_num,c,i,s->channel[c].subframe_len[i]);
        }
    }

    return 0;
}


/**
 *@brief decode one wma frame
 *@param s context
 *@param gb current get bit context
 *@return 0 if the trailer bit indicates that this is the last frame
 *        1 if there are more frames
 */
static int wma_decode_frame(WMA3DecodeContext *s,GetBitContext* gb){
    unsigned int gb_start_count = get_bits_count(gb);
    int more_frames = 0;
    /** get frame length */
    int len = 0;

    if(s->len_prefix)
        len = get_bits(gb,s->log2_frame_size);

    av_log(s->avctx,AV_LOG_INFO,"decoding frame with len %x\n",len);

    /** decode tile information */
    if(wma_decode_tilehdr(s,gb)){
        s->packet_loss = 1;
        return 0;
    }


    /** read postproc transform */
    if(get_bits1(gb)){
        av_log(s->avctx,AV_LOG_ERROR,"Unsupported postproc transform found\n");
        s->packet_loss = 1;
        return 0;
    }

    /** read drc info */
    if(s->dynamic_range_compression){
        s->drc_gain = get_bits(gb,8);
        av_log(s->avctx,AV_LOG_INFO,"drc_gain %i\n",s->drc_gain);
    }

    s->update_samples_per_frame = 0;

    /** check if num output samples might change */
    if(get_bits(gb,1)){
        s->update_samples_per_frame = get_bits1(gb);

        if(s->update_samples_per_frame){
            get_bits(gb,av_log2(s->samples_per_frame * 2));
        }
     }
    /** skip the rest of the frame data */
    skip_bits_long(gb,len - (get_bits_count(gb) - gb_start_count) - 1);

    /** decode trailer bit */
    more_frames = get_bits1(gb);

    ++s->frame_num;
    return more_frames;
}


/**
 *@brief calculate remaining input buffer len
 *@param s codec context
 *@return remaining size in bits
 */
static int remaining_bits(WMA3DecodeContext *s){
    return s->buf_bit_size - get_bits_count(&s->gb);
}

/**
 *@brief fill the bit reservoir with a partial frame
 *@param s codec context
 *@param len length of the partial frame
 */
static void save_bits(WMA3DecodeContext *s,int len){
    int buflen = (s->prev_frame_bit_size + len + 8) / 8;
    int bit_offset = s->prev_frame_bit_size % 8;
    int pos = (s->prev_frame_bit_size - bit_offset) / 8;
    s->prev_frame_bit_size += len;

    if(len <= 0)
         return;

    /** increase length if needed */
    s->prev_frame = av_realloc(s->prev_frame,buflen +
                               FF_INPUT_BUFFER_PADDING_SIZE);

    /** byte align prev_frame buffer */
    if(bit_offset){
        int missing = 8 - bit_offset;
        if(len < missing)
            missing = len;
        s->prev_frame[pos++] |=
            get_bits(&s->gb, missing) << (8 - bit_offset - missing);
        len -= missing;
    }

    /** copy full bytes */
    while(len > 7){
        s->prev_frame[pos++] = get_bits(&s->gb,8);
        len -= 8;
    }

    /** copy remaining bits */
    if(len > 0)
        s->prev_frame[pos++] = get_bits(&s->gb,len) << (8 - len);
}

/**
 *@brief decode a single wma packet
 *@param avctx codec context
 *@param data the output buffer
 *@param data_size number of bytes that were written to the output buffer
 *@param buf input buffer
 *@param buf_size input buffer length
 *@return number of bytes that were read from the input buffer
 */
static int wma3_decode_packet(AVCodecContext *avctx,
                             void *data, int *data_size,
                             const uint8_t *buf, int buf_size)
{
    WMA3DecodeContext *s = avctx->priv_data;
    int more_frames=1;
    int num_bits_prev_frame;
    s->buf_bit_size = buf_size << 3;

    *data_size = 0;

    /** sanity check for the buffer length */
    if(buf_size < avctx->block_align)
        return 0;

    /** Parse packet header */
    init_get_bits(&s->gb, buf, s->buf_bit_size);
    s->packet_sequence_number = get_bits(&s->gb, 4);
    s->bit5                   = get_bits1(&s->gb);
    s->bit6                   = get_bits1(&s->gb);

    /** get number of bits that need to be added to the previous frame */
    num_bits_prev_frame = get_bits(&s->gb, s->log2_frame_size);
    av_log(avctx, AV_LOG_INFO, "[%d]: nbpf %x\n", avctx->frame_number,
                  num_bits_prev_frame);

    /** check for packet loss */
    if (s->packet_sequence_number != (avctx->frame_number&0xF)) {
        s->packet_loss = 1;
        av_log(avctx, AV_LOG_ERROR, "!!Packet loss detected! seq %x vs %x\n",
                      s->packet_sequence_number,avctx->frame_number&0xF);
    }

    if (num_bits_prev_frame > 0) {
        /** append the prev frame data to the remaining data from the
            previous packet to create a full frame */
        save_bits(s,num_bits_prev_frame);
        av_log(avctx, AV_LOG_INFO, "accumulated %x bits of frame data\n",
                      s->prev_frame_bit_size);

        /** decode the cross packet frame if it is valid */
        if(!s->packet_loss){
            GetBitContext gb_prev;
            init_get_bits(&gb_prev, s->prev_frame, s->prev_frame_bit_size);
            wma_decode_frame(s,&gb_prev);
        }
    }else if(s->prev_frame_bit_size){
        av_log(avctx, AV_LOG_ERROR, "ignoring %x previously saved bits\n",
                      s->prev_frame_bit_size);
    }

    /** reset prev frame buffer */
    s->prev_frame_bit_size = 0;
    s->packet_loss = 0;
    /** decode the rest of the packet */
    while(more_frames && remaining_bits(s) > s->log2_frame_size){
        int frame_size = show_bits(&s->gb, s->log2_frame_size);

        /** there is enough data for a full frame */
        if(remaining_bits(s) >= frame_size){
            /** decode the frame */
            more_frames = wma_decode_frame(s,&s->gb);

            if(!more_frames){
                av_log(avctx, AV_LOG_ERROR, "no more frames\n");
            }
        }else
            more_frames = 0;
    }

    /** save the rest of the data so that it can be decoded
       with the next packet */
    save_bits(s,remaining_bits(s));

    return avctx->block_align;
}

/**
 *@brief uninitialize the decoder and free all ressources
 *@param avctx codec context
 *@return 0 on success, < 0 otherwise
 */
static av_cold int wma3_decode_end(AVCodecContext *avctx)
{
    WMA3DecodeContext *s = avctx->priv_data;
    if(s->prev_frame)
        av_free(s->prev_frame);
    return 0;
}

/**
 *@brief WMA9 decoder
 */
AVCodec wmav3pro_decoder =
{
    "wmav3pro",
    CODEC_TYPE_AUDIO,
    CODEC_ID_WMAPRO,
    sizeof(WMA3DecodeContext),
    wma3_decode_init,
    NULL,
    wma3_decode_end,
    wma3_decode_packet,
    .long_name = "Windows Media Audio 9 Professional",
};
