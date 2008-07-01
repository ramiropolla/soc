/*
 * MLP decoder
 * Copyright (c) 2007-2008 Ian Caulfield
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
 * @file mlpdec.c
 * MLP decoder
 */

#include "avcodec.h"
#include "libavutil/intreadwrite.h"
#include "bitstream.h"
#include "libavutil/crc.h"
#include "parser.h"
#include "mlp_parser.h"

/** Maximum number of channels that can be decoded. */
#define MAX_CHANNELS        16

/** Maximum number of matrices used in decoding. Most streams have one matrix
 *  per output channel, but some rematrix a channel (usually 0) more than once.
 */

#define MAX_MATRICES        15

/** Maximum number of substreams that can be decoded. This could also be set
 *  higher, but again I haven't seen any examples with more than two. */
#define MAX_SUBSTREAMS      2

/** Maximum sample frequency supported. */
#define MAX_SAMPLERATE      192000

/** The maximum number of audio samples within one access unit. */
#define MAX_BLOCKSIZE       (40 * (MAX_SAMPLERATE / 48000))
/** The next power of two greater than MAX_BLOCKSIZE. */
#define MAX_BLOCKSIZE_POW2  (64 * (MAX_SAMPLERATE / 48000))

/** Number of allowed filters. */
#define NUM_FILTERS         2

/** The maximum number of taps in either the IIR or FIR filter.
 *  I believe MLP actually specifies the maximum order for IIR filters is four,
 *  and that the sum of the orders of both filters must be <= 8. */
#define MAX_FILTER_ORDER    8

/** Number of bits used for VLC lookup - longest huffman code is 9. */
#define VLC_BITS            9


static const char* sample_message =
    "Please file a bug report following the instructions at "
    "http://ffmpeg.mplayerhq.hu/bugreports.html and include "
    "a sample of this file.";

typedef struct SubStream {
    //! For each substream, whether a restart header has been read
    uint8_t     restart_seen;

    //@{
    /** Restart header data */
    //! The sync word used at the start of the last restart header
    uint16_t    restart_sync_word;

    //! The index of the first channel coded in this substream
    uint8_t     min_channel;
    //! The index of the last channel coded in this substream
    uint8_t     max_channel;
    //! The number of channels input into the rematrix stage
    uint8_t     max_matrix_channel;

    //! The left shift applied to random noise in 0x31ea substreams
    uint8_t     noise_shift;
    //! The current seed value for the pseudorandom noise generator(s)
    uint32_t    noisegen_seed;

    //! Does this substream contain extra info to check the size of VLC blocks?
    uint8_t     data_check_present;

    //! Bitmask of which parameter sets are conveyed in a decoding parameter block
    uint8_t     param_presence_flags;
#define PARAM_BLOCKSIZE     (1 << 7)
#define PARAM_MATRIX        (1 << 6)
#define PARAM_OUTSHIFT      (1 << 5)
#define PARAM_QUANTSTEP     (1 << 4)
#define PARAM_FIR           (1 << 3)
#define PARAM_IIR           (1 << 2)
#define PARAM_HUFFOFFSET    (1 << 1)
    //@}

    //@{
    /** Matrix data */

    //! Number of matrices to be applied
    uint8_t     num_primitive_matrices;

    //! Output channel of matrix
    uint8_t     matrix_ch[MAX_MATRICES];

    //! Whether the LSBs of the matrix output are encoded in the bitstream
    uint8_t     lsb_bypass[MAX_MATRICES];
    //! Matrix coefficients, stored as 2.14 fixed point
    int32_t     matrix_coeff[MAX_MATRICES][MAX_CHANNELS+2];
    //! Left shift to apply to noise values in 0x31eb substreams
    uint8_t     matrix_noise_shift[MAX_MATRICES];
    //@}

    //! Left shift to apply to huffman-decoded residuals
    uint8_t     quant_step_size[MAX_CHANNELS];

    //! Number of PCM samples in current audio block
    uint16_t    blocksize;
    //! Number of PCM samples decoded so far in this frame
    uint16_t    blockpos;

    //! Left shift to apply to decoded PCM values to get final 24-bit output
    int8_t      output_shift[MAX_CHANNELS];

    //! Running XOR of all output samples
    int32_t     lossless_check_data;

} SubStream;

typedef struct MLPDecodeContext {
    AVCodecContext *avctx;

    //! Do we have valid stream data read from a major sync block?
    uint8_t     params_valid;

    //! Number of substreams contained within this stream
    uint8_t     num_substreams;

    //! Index of the last substream to decode - further substreams are skipped
    uint8_t     max_decoded_substream;

    //! Number of PCM samples contained in each frame
    int         access_unit_size;
    //! Next power of two above the number of samples in each frame
    int         access_unit_size_pow2;

    SubStream   substream[MAX_SUBSTREAMS];

    //@{
    /** Filter data. */
#define FIR 0
#define IIR 1
    //! Number of taps in filter
    uint8_t     filter_order[MAX_CHANNELS][NUM_FILTERS];
    //! Right shift to apply to output of filter
    uint8_t     filter_coeff_q[MAX_CHANNELS][NUM_FILTERS];

    int32_t     filter_coeff[MAX_CHANNELS][NUM_FILTERS][MAX_FILTER_ORDER];
    int32_t     filter_state[MAX_CHANNELS][NUM_FILTERS][MAX_FILTER_ORDER];
    int32_t     filter_state_buffer       [NUM_FILTERS][MAX_BLOCKSIZE + MAX_FILTER_ORDER];
    //@}

    //@{
    /** Sample data coding infomation */
    //! Offset to apply to residual values
    int16_t     huff_offset[MAX_CHANNELS];
    //! Sign/rounding corrected version of huff_offset
    int32_t     sign_huff_offset[MAX_CHANNELS];
    //! Which VLC codebook to use to read residuals
    uint8_t     codebook[MAX_CHANNELS];
    //! Size of residual suffix not encoded using VLC
    uint8_t     huff_lsbs[MAX_CHANNELS];
    //@}

    int8_t      noise_buffer[MAX_BLOCKSIZE_POW2];
    int8_t      bypassed_lsbs[MAX_BLOCKSIZE][MAX_CHANNELS];
    int32_t     sample_buffer[MAX_BLOCKSIZE][MAX_CHANNELS+2];
} MLPDecodeContext;

/** Tables defining the huffman codes.
 *  There are three entropy coding methods used in MLP (four if you count "none"
 *  as a method). These use the same sequences for codes starting 00... or 01...
 *  but have different codes starting 1....
 */

static const uint8_t huffman_tables[3][18][2] = {
    { /* Huffman table 0, -7 - +10 */
        {0x01, 9}, {0x01, 8}, {0x01, 7}, {0x01, 6}, {0x01, 5}, {0x01, 4}, {0x01, 3},
        {0x04, 3}, {0x05, 3}, {0x06, 3}, {0x07, 3},
        {0x03, 3}, {0x05, 4}, {0x09, 5}, {0x11, 6}, {0x21, 7}, {0x41, 8}, {0x81, 9},
    }, { /* Huffman table 1, -7 - +8 */
        {0x01, 9}, {0x01, 8}, {0x01, 7}, {0x01, 6}, {0x01, 5}, {0x01, 4}, {0x01, 3},
        {0x02, 2}, {0x03, 2},
        {0x03, 3}, {0x05, 4}, {0x09, 5}, {0x11, 6}, {0x21, 7}, {0x41, 8}, {0x81, 9},
    }, { /* Huffman table 2, -7 - +7 */
        {0x01, 9}, {0x01, 8}, {0x01, 7}, {0x01, 6}, {0x01, 5}, {0x01, 4}, {0x01, 3},
        {0x01, 1},
        {0x03, 3}, {0x05, 4}, {0x09, 5}, {0x11, 6}, {0x21, 7}, {0x41, 8}, {0x81, 9},
    }
};

static VLC huff_vlc[3];

static AVCRC crc_63[1024];
static AVCRC crc_1D[1024];


/** Initialize static data, constant between all invocations of the codec. */

static void init_static()
{
    if (!huff_vlc[0].bits) {
        INIT_VLC_STATIC(&huff_vlc[0], VLC_BITS, 18,
                 &huffman_tables[0][0][1], 2, 1,
                 &huffman_tables[0][0][0], 2, 1, 512);
        INIT_VLC_STATIC(&huff_vlc[1], VLC_BITS, 16,
                 &huffman_tables[1][0][1], 2, 1,
                 &huffman_tables[1][0][0], 2, 1, 512);
        INIT_VLC_STATIC(&huff_vlc[2], VLC_BITS, 15,
                 &huffman_tables[2][0][1], 2, 1,
                 &huffman_tables[2][0][0], 2, 1, 512);

        av_crc_init(crc_63, 0,  8,   0x63, sizeof(crc_63));
        av_crc_init(crc_1D, 0,  8,   0x1D, sizeof(crc_1D));
    }
}


/** MLP uses checksums that seem to be based on the standard CRC algorithm,
 *  but not (in implementation terms, the table lookup and XOR are reversed).
 *  We can implement this behaviour using a standard av_crc on all but the
 *  last element, then XOR that with the last element.
 */

static uint8_t mlp_checksum8(const uint8_t *buf, unsigned int buf_size)
{
    uint8_t checksum = av_crc(crc_63, 0x3c, buf, buf_size - 1); // crc_63[0xa2] == 0x3c
    checksum ^= buf[buf_size-1];
    return checksum;
}

/** Calculate an 8-bit checksum over a restart header -- a non-multiple-of-8
 *  number of bits, starting two bits into the first byte of buf.
 */

static uint8_t mlp_restart_checksum(const uint8_t *buf, unsigned int bit_size)
{
    int i;
    int num_bytes = (bit_size + 2) / 8;

    int crc = crc_1D[buf[0] & 0x3f];
    crc = av_crc(crc_1D, crc, buf + 1, num_bytes - 2);
    crc ^= buf[num_bytes - 1];

    for (i = 0; i < ((bit_size + 2) & 7); i++) {
        crc <<= 1;
        if (crc & 0x100)
            crc ^= 0x11D;
        crc ^= (buf[num_bytes] >> (7 - i)) & 1;
    }

    return crc;
}

static inline void calculate_sign_huff(MLPDecodeContext *m, unsigned int substr,
                                       unsigned int ch)
{
    SubStream *s = &m->substream[substr];
    int lsb_bits = m->huff_lsbs[ch] - s->quant_step_size[ch];
    int sign_shift = lsb_bits + (m->codebook[ch] ? 2 - m->codebook[ch] : -1);

    m->sign_huff_offset[ch] = m->huff_offset[ch];

    if (m->codebook[ch] > 0)
        m->sign_huff_offset[ch] -= 7 << lsb_bits;

    if (sign_shift >= 0)
        m->sign_huff_offset[ch] -= 1 << sign_shift;
}

/** Read a sample, consisting of either, both or neither of entropy-coded MSBs
 *  and plain LSBs.
 */

static inline int read_huff_channels(MLPDecodeContext *m, GetBitContext *gbp,
                                     unsigned int substr, unsigned int pos)
{
    SubStream *s = &m->substream[substr];
    unsigned int mat, channel;

    for (mat = 0; mat < s->num_primitive_matrices; mat++)
        if (s->lsb_bypass[mat])
            m->bypassed_lsbs[pos + s->blockpos][mat] = get_bits1(gbp);

    for (channel = s->min_channel; channel <= s->max_channel; channel++) {
        int codebook = m->codebook[channel];
        int quant_step_size = s->quant_step_size[channel];
        int lsb_bits = m->huff_lsbs[channel] - quant_step_size;
        int result = 0;

        if (codebook > 0)
            result = get_vlc2(gbp, huff_vlc[codebook-1].table,
                            VLC_BITS, (9 + VLC_BITS - 1) / VLC_BITS);

        if (result < 0)
            return -1;

        if (lsb_bits > 0)
            result = (result << lsb_bits) + get_bits(gbp, lsb_bits);

        result  += m->sign_huff_offset[channel];
        result <<= quant_step_size;

        m->sample_buffer[pos + s->blockpos][channel] = result;
    }

    return 0;
}

static int mlp_decode_init(AVCodecContext *avctx)
{
    MLPDecodeContext *m = avctx->priv_data;
    int substr;

    init_static();
    m->avctx = avctx;
    for (substr = 0; substr < MAX_SUBSTREAMS; substr++)
        m->substream[substr].lossless_check_data = 0xffffffff;
    return 0;
}

/** Read a major sync info header - contains high level information about
 *  the stream - sample rate, channel arrangement etc. Most of this
 *  information is not actually necessary for decoding, only for playback.
 */

static int read_major_sync(MLPDecodeContext *m, GetBitContext *gb)
{
    MLPHeaderInfo mh;
    int substr;

    if (ff_mlp_read_major_sync(m->avctx, &mh, gb) != 0)
        return -1;

    if (mh.group1_bits == 0) {
        av_log(m->avctx, AV_LOG_ERROR, "Invalid/unknown bits per sample\n");
        return -1;
    }
    if (mh.group2_bits > mh.group1_bits) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Channel group 2 cannot have more bits per sample than group 1\n");
        return -1;
    }

    if (mh.group2_samplerate && mh.group2_samplerate != mh.group1_samplerate) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Channel groups with differing sample rates not currently supported\n");
        return -1;
    }

    if (mh.group1_samplerate == 0) {
        av_log(m->avctx, AV_LOG_ERROR, "Invalid/unknown sampling rate\n");
        return -1;
    }
    if (mh.group1_samplerate > MAX_SAMPLERATE) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Sampling rate %d is greater than maximum supported (%d)\n",
               mh.group1_samplerate, MAX_SAMPLERATE);
        return -1;
    }
    if (mh.access_unit_size > MAX_BLOCKSIZE) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Block size %d is greater than maximum supported (%d)\n",
               mh.access_unit_size, MAX_BLOCKSIZE);
        return -1;
    }
    if (mh.access_unit_size_pow2 > MAX_BLOCKSIZE_POW2) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Block size pow2 %d is greater than maximum supported (%d)\n",
               mh.access_unit_size_pow2, MAX_BLOCKSIZE_POW2);
        return -1;
    }

    if (mh.num_substreams == 0)
        return -1;
    if (mh.num_substreams > MAX_SUBSTREAMS) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Number of substreams %d is more than maximum supported by "
               "decoder. %s\n", mh.num_substreams, sample_message);
        return -1;
    }

    m->access_unit_size      = mh.access_unit_size;
    m->access_unit_size_pow2 = mh.access_unit_size_pow2;

    m->num_substreams        = mh.num_substreams;
    m->max_decoded_substream = m->num_substreams - 1;

    m->avctx->sample_rate    = mh.group1_samplerate;
    m->avctx->frame_size     = mh.access_unit_size;

#ifdef CONFIG_AUDIO_NONSHORT
    m->avctx->bits_per_sample = mh.group1_bits;
    if (mh.group1_bits > 16) {
        m->avctx->sample_fmt = SAMPLE_FMT_S32;
    }
#endif

    m->params_valid = 1;
    for (substr = 0; substr < MAX_SUBSTREAMS; substr++)
        m->substream[substr].restart_seen = 0;

    return 0;
}

/** Read a restart header from a block in a substream. This contains parameters
 *  required to decode the audio that do not change very often. Generally
 *  (always) present only in blocks following a major sync.
 */

static int read_restart_header(MLPDecodeContext *m, GetBitContext *gbp,
                               const uint8_t *buf, unsigned int substr)
{
    SubStream *s = &m->substream[substr];
    unsigned int ch;
    int sync_word, tmp;
    uint8_t checksum;
    uint8_t lossless_check;
    int start_count = get_bits_count(gbp);

    sync_word = get_bits(gbp, 14);

    if ((sync_word & 0x3ffe) != 0x31ea) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Restart header sync incorrect (got 0x%04x)\n", sync_word);
        return -1;
    }
    s->restart_sync_word = sync_word;

    skip_bits(gbp, 16); /* Output timestamp */

    s->min_channel        = get_bits(gbp, 4);
    s->max_channel        = get_bits(gbp, 4);
    s->max_matrix_channel = get_bits(gbp, 4);

    if (s->min_channel > s->max_channel) {
        av_log(m->avctx, AV_LOG_ERROR,
               "Substream min channel cannot be greater than max channel.\n");
        return -1;
    }

    if (m->avctx->request_channels > 0
        && s->max_channel + 1 >= m->avctx->request_channels
        && substr < m->max_decoded_substream) {
        av_log(m->avctx, AV_LOG_INFO,
               "Extracting %d channel downmix from substream %d. "
               "Further substreams will be skipped.\n",
               s->max_channel + 1, substr);
        m->max_decoded_substream = substr;
    }

    s->noise_shift   = get_bits(gbp,  4);
    s->noisegen_seed = get_bits(gbp, 23);

    skip_bits(gbp, 19);

    s->data_check_present = get_bits1(gbp);
    lossless_check = get_bits(gbp, 8);
    if (substr == m->max_decoded_substream
        && s->lossless_check_data != 0xffffffff) {
        tmp = s->lossless_check_data;
        tmp ^= tmp >> 16;
        tmp ^= tmp >> 8;
        tmp &= 0xff;
        if (tmp != lossless_check)
            av_log(m->avctx, AV_LOG_WARNING,
                   "Lossless check failed - expected %02x, calculated %02x\n",
                   lossless_check, tmp);
        else
            dprintf(m->avctx, "Lossless check passed for substream %d (%x)\n",
                    substr, tmp);
    }

    skip_bits(gbp, 16);

    for (ch = 0; ch <= s->max_matrix_channel; ch++) {
        int ch_assign = get_bits(gbp, 6);
        dprintf(m->avctx, "ch_assign[%d][%d] = %d\n", substr, ch,
                ch_assign);
        if (ch_assign != ch) {
            av_log(m->avctx, AV_LOG_ERROR,
                   "Non 1:1 channel assignments are used in this stream. %s\n",
                   sample_message);
            return -1;
        }
    }

    checksum = mlp_restart_checksum(buf, get_bits_count(gbp) - start_count);

    if (checksum != get_bits(gbp, 8))
        av_log(m->avctx, AV_LOG_ERROR, "Restart header checksum error\n");

    /* Set default decoding parameters */
    s->param_presence_flags   = 0xff;
    s->num_primitive_matrices = 0;
    s->blocksize              = 8;
    s->lossless_check_data    = 0;

    memset(s->output_shift   , 0, sizeof(s->output_shift   ));
    memset(s->quant_step_size, 0, sizeof(s->quant_step_size));

    for (ch = s->min_channel; ch <= s->max_channel; ch++) {
        m->filter_order  [ch][FIR] = 0;
        m->filter_order  [ch][IIR] = 0;
        m->filter_coeff_q[ch][FIR] = 0;
        m->filter_coeff_q[ch][IIR] = 0;

        memset(m->filter_coeff[ch], 0, sizeof(m->filter_coeff[ch]));
        memset(m->filter_state[ch], 0, sizeof(m->filter_state[ch]));

        /* Default audio coding is 24-bit raw PCM */
        m->huff_offset     [ch] = 0;
        m->sign_huff_offset[ch] = (-1) << 23;
        m->codebook        [ch] = 0;
        m->huff_lsbs       [ch] = 24;
    }

    if (substr == m->max_decoded_substream) {
        m->avctx->channels = s->max_channel + 1;
    }

    return 0;
}

/** Read parameters for one of the prediction filters.
 */

static int read_filter_params(MLPDecodeContext *m, GetBitContext *gbp,
                              unsigned int channel, unsigned int filter)
{
    const char fchar = filter ? 'I' : 'F';
    int i, order;

    // filter is 0 for FIR, 1 for IIR
    assert(filter < 2);

    order = get_bits(gbp, 4);
    if (order > MAX_FILTER_ORDER) {
        av_log(m->avctx, AV_LOG_ERROR,
               "%cIR filter order %d is greater than maximum %d\n",
               fchar, order, MAX_FILTER_ORDER);
        return -1;
    }
    m->filter_order[channel][filter] = order;

    if (order > 0) {
        int coeff_bits, coeff_shift;

        m->filter_coeff_q[channel][filter] = get_bits(gbp, 4);

        coeff_bits  = get_bits(gbp, 5);
        coeff_shift = get_bits(gbp, 3);
        if (coeff_bits < 1 || coeff_bits > 16) {
            av_log(m->avctx, AV_LOG_ERROR,
                   "%cIR filter coeff_bits must be between 1 and 16\n",
                   fchar);
            return -1;
        }
        if (coeff_bits + coeff_shift > 16) {
            av_log(m->avctx, AV_LOG_ERROR,
                   "Sum of coeff_bits and coeff_shift for %cIR filter must be 16 or less\n",
                   fchar);
            return -1;
        }

        for (i = 0; i < order; i++)
            m->filter_coeff[channel][filter][i] =
                    get_sbits(gbp, coeff_bits) << coeff_shift;

        if (get_bits1(gbp)) {
            int state_bits, state_shift;

            if (filter == FIR) {
                av_log(m->avctx, AV_LOG_ERROR,
                       "FIR filter has state data specified\n");
                return -1;
            }

            state_bits  = get_bits(gbp, 4);
            state_shift = get_bits(gbp, 4);

            /* TODO: check validity of state data */

            for (i = 0; i < order; i++)
                m->filter_state[channel][filter][i] =
                    get_sbits(gbp, state_bits) << state_shift;
        }
    }

    return 0;
}

/** Read decoding parameters that change more often than those in the restart
 *  header.
 */

static int read_decoding_params(MLPDecodeContext *m, GetBitContext *gbp,
                                unsigned int substr)
{
    SubStream *s = &m->substream[substr];
    unsigned int mat, ch;

    if (get_bits1(gbp))
        s->param_presence_flags = get_bits(gbp, 8);

    if (s->param_presence_flags & PARAM_BLOCKSIZE)
        if (get_bits1(gbp)) {
            s->blocksize = get_bits(gbp, 9);
            if (s->blocksize > MAX_BLOCKSIZE) {
                av_log(m->avctx, AV_LOG_ERROR, "Block size too large\n");
                s->blocksize = 0;
                return -1;
            }
        }

    if (s->param_presence_flags & PARAM_MATRIX)
        if (get_bits1(gbp)) {
            s->num_primitive_matrices = get_bits(gbp, 4);

            for (mat = 0; mat < s->num_primitive_matrices; mat++) {
                int frac_bits, max_chan;
                s->matrix_ch[mat] = get_bits(gbp, 4);
                frac_bits = get_bits(gbp, 4);
                s->lsb_bypass[mat] = get_bits1(gbp);

                if (s->matrix_ch[mat] > s->max_channel) {
                    av_log(m->avctx, AV_LOG_ERROR,
                           "Invalid channel %d specified as output from matrix\n",
                           s->matrix_ch[mat]);
                    return -1;
                }
                if (frac_bits > 14) {
                    av_log(m->avctx, AV_LOG_ERROR,
                           "Too many fractional bits specified\n");
                    return -1;
                }

                max_chan = s->max_matrix_channel;
                if (s->restart_sync_word == 0x31ea)
                    max_chan+=2;

                for (ch = 0; ch <= max_chan; ch++) {
                    int coeff_val = 0;
                    if (get_bits1(gbp))
                        coeff_val = get_sbits(gbp, frac_bits + 2);

                    s->matrix_coeff[mat][ch] = coeff_val << (14 - frac_bits);
                }

                if (s->restart_sync_word == 0x31eb)
                    s->matrix_noise_shift[mat] = get_bits(gbp, 4);
                else
                    s->matrix_noise_shift[mat] = 0;
            }
        }

    if (s->param_presence_flags & PARAM_OUTSHIFT)
        if (get_bits1(gbp))
            for (ch = 0; ch <= s->max_matrix_channel; ch++) {
                s->output_shift[ch] = get_bits(gbp, 4);
                dprintf(m->avctx, "output shift[%d] = %d\n",
                        ch, s->output_shift[ch]);
                /* TODO: validate */
            }

    if (s->param_presence_flags & PARAM_QUANTSTEP)
        if (get_bits1(gbp))
            for (ch = 0; ch <= s->max_channel; ch++) {
                s->quant_step_size[ch] = get_bits(gbp, 4);
                /* TODO: validate */

                calculate_sign_huff(m, substr, ch);
            }

    for (ch = s->min_channel; ch <= s->max_channel; ch++)
        if (get_bits1(gbp)) {
            if (s->param_presence_flags & PARAM_FIR)
                if (get_bits1(gbp))
                    if (read_filter_params(m, gbp, ch, FIR) < 0)
                        return -1;

            if (s->param_presence_flags & PARAM_IIR)
                if (get_bits1(gbp))
                    if (read_filter_params(m, gbp, ch, IIR) < 0)
                        return -1;

            if (m->filter_order  [ch][FIR] && m->filter_order  [ch][IIR] &&
                m->filter_coeff_q[ch][FIR] != m->filter_coeff_q[ch][IIR]) {
                av_log(m->avctx, AV_LOG_ERROR,
                       "FIR and IIR filters must use same precision\n");
                return -1;
            }
            /* Both filters must have the same precision, so the filtering
             * code always use the FIR precision. If only IIR is used, we copy
             * its precision to FIR. */
            if (!m->filter_order[ch][FIR] && m->filter_order[ch][IIR])
                m->filter_coeff_q[ch][FIR] = m->filter_coeff_q[ch][IIR];

            if (s->param_presence_flags & PARAM_HUFFOFFSET)
                if (get_bits1(gbp))
                    m->huff_offset[ch] = get_sbits(gbp, 15);

            m->codebook [ch] = get_bits(gbp, 2);
            m->huff_lsbs[ch] = get_bits(gbp, 5);

            calculate_sign_huff(m, substr, ch);

            /* TODO: validate */
        }

    return 0;
}

/** Generate PCM samples using the prediction filters and residual values
 *  read from the data stream, and update the filter state.
 */

static void filter_channel(MLPDecodeContext *m, unsigned int substr,
                           unsigned int channel)
{
    SubStream *s = &m->substream[substr];
    unsigned int quant_step_size = s->quant_step_size[channel];
    unsigned int filter_coeff_q = m->filter_coeff_q[channel][FIR];
    int index = MAX_BLOCKSIZE;
    int j, i;

    for (j = 0; j < NUM_FILTERS; j++) {
        memcpy(&m->filter_state_buffer  [j][MAX_BLOCKSIZE],
               &m->filter_state[channel][j][0],
               MAX_FILTER_ORDER * sizeof(int32_t));
    }

    for (i = 0; i < s->blocksize; i++) {
        int32_t residual = m->sample_buffer[i + s->blockpos][channel];
        unsigned int order;
        int64_t accum = 0;
        int32_t result;

        /* TODO: Move this code to DSPContext? */

        for (j = 0; j < NUM_FILTERS; j++)
            for (order = 0; order < m->filter_order[channel][j]; order++)
                accum += (int64_t)m->filter_state_buffer[j][index + order] *
                        m->filter_coeff[channel][j][order];

        accum  = accum >> filter_coeff_q;
        result = (accum + residual) & ~((1 << quant_step_size) - 1);

        --index;

        m->filter_state_buffer[FIR][index] = result;
        m->filter_state_buffer[IIR][index] = result - accum;

        m->sample_buffer[i + s->blockpos][channel] = result;
    }

    for (j = 0; j < NUM_FILTERS; j++) {
        memcpy(&m->filter_state[channel][j][0],
               &m->filter_state_buffer  [j][index],
               MAX_FILTER_ORDER * sizeof(int32_t));
    }
}

/** Read a block of PCM residual (or actual if no filtering active) data.
 */

static int read_block_data(MLPDecodeContext *m, GetBitContext *gbp,
                           unsigned int substr)
{
    SubStream *s = &m->substream[substr];
    unsigned int i, ch, expected_stream_pos = 0;

    if (s->data_check_present) {
        expected_stream_pos = get_bits_count(gbp) + get_bits(gbp, 16);
        av_log(m->avctx, AV_LOG_WARNING, "This file contains some features "
               "we have not tested yet. %s\n", sample_message);
    }

    if (s->blockpos + s->blocksize > m->access_unit_size) {
        av_log(m->avctx, AV_LOG_ERROR, "Too many audio samples in frame\n");
        return -1;
    }

    memset(&m->bypassed_lsbs[s->blockpos][0], 0,
           s->blocksize * sizeof(m->bypassed_lsbs[0]));

    for (i = 0; i < s->blocksize; i++) {
        if (read_huff_channels(m, gbp, substr, i) < 0)
            return -1;
    }

    for (ch = s->min_channel; ch <= s->max_channel; ch++) {
        filter_channel(m, substr, ch);
    }

    s->blockpos += s->blocksize;

    if (s->data_check_present) {
        if (get_bits_count(gbp) != expected_stream_pos)
            av_log(m->avctx, AV_LOG_ERROR, "Block data length mismatch\n");
        skip_bits(gbp, 8);
    }

    return 0;
}

/** Data table used for TrueHD noise generation function */

static const int8_t noise_table[256] = {
     30,  51,  22,  54,   3,   7,  -4,  38,  14,  55,  46,  81,  22,  58,  -3,   2,
     52,  31,  -7,  51,  15,  44,  74,  30,  85, -17,  10,  33,  18,  80,  28,  62,
     10,  32,  23,  69,  72,  26,  35,  17,  73,  60,   8,  56,   2,   6,  -2,  -5,
     51,   4,  11,  50,  66,  76,  21,  44,  33,  47,   1,  26,  64,  48,  57,  40,
     38,  16, -10, -28,  92,  22, -18,  29, -10,   5, -13,  49,  19,  24,  70,  34,
     61,  48,  30,  14,  -6,  25,  58,  33,  42,  60,  67,  17,  54,  17,  22,  30,
     67,  44,  -9,  50, -11,  43,  40,  32,  59,  82,  13,  49, -14,  55,  60,  36,
     48,  49,  31,  47,  15,  12,   4,  65,   1,  23,  29,  39,  45,  -2,  84,  69,
      0,  72,  37,  57,  27,  41, -15, -16,  35,  31,  14,  61,  24,   0,  27,  24,
     16,  41,  55,  34,  53,   9,  56,  12,  25,  29,  53,   5,  20, -20,  -8,  20,
     13,  28,  -3,  78,  38,  16,  11,  62,  46,  29,  21,  24,  46,  65,  43, -23,
     89,  18,  74,  21,  38, -12,  19,  12, -19,   8,  15,  33,   4,  57,   9,  -8,
     36,  35,  26,  28,   7,  83,  63,  79,  75,  11,   3,  87,  37,  47,  34,  40,
     39,  19,  20,  42,  27,  34,  39,  77,  13,  42,  59,  64,  45,  -1,  32,  37,
     45,  -5,  53,  -6,   7,  36,  50,  23,   6,  32,   9, -21,  18,  71,  27,  52,
    -25,  31,  35,  42,  -1,  68,  63,  52,  26,  43,  66,  37,  41,  25,  40,  70,
};

/** Noise generation functions.
 *  I'm not sure what these are for - they seem to be some kind of pseudorandom
 *  sequence generators, used to generate noise data which is used when the
 *  channels are rematrixed. I'm not sure if they provide a practical benefit
 *  to compression, or just obfuscate the decoder. Are they for some kind of
 *  dithering?
 */

/** Generate two channels of noise, used in the matrix when restart_sync_word == 0x31ea. */

static void generate_noise_1(MLPDecodeContext *m, unsigned int substr)
{
    SubStream *s = &m->substream[substr];
    unsigned int i;
    uint32_t seed = s->noisegen_seed;
    unsigned int maxchan = s->max_matrix_channel;

    for (i = 0; i < s->blockpos; i++) {
        uint16_t seed_shr7 = seed >> 7;
        m->sample_buffer[i][maxchan+1] = ((int8_t)(seed >> 15)) << s->noise_shift;
        m->sample_buffer[i][maxchan+2] = ((int8_t) seed_shr7)   << s->noise_shift;

        seed = (seed << 16) ^ seed_shr7 ^ (seed_shr7 << 5);
    }

    s->noisegen_seed = seed;
}

/** Generate a block of noise, used when restart_sync_word == 0x31eb. */

static void generate_noise_2(MLPDecodeContext *m, unsigned int substr)
{
    SubStream *s = &m->substream[substr];
    unsigned int i;
    uint32_t seed = s->noisegen_seed;

    for (i = 0; i < m->access_unit_size_pow2; i++) {
        uint8_t seed_shr15 = seed >> 15;
        m->noise_buffer[i] = noise_table[seed_shr15];
        seed = (seed << 8) ^ seed_shr15 ^ (seed_shr15 << 5);
    }

    s->noisegen_seed = seed;
}


/** Apply the channel matrices in turn to reconstruct the original audio samples.
 */

static void rematrix_channels(MLPDecodeContext *m, unsigned int substr)
{
    SubStream *s = &m->substream[substr];
    unsigned int mat, src_ch, i;
    unsigned int maxchan;

    maxchan = s->max_matrix_channel;
    if (s->restart_sync_word == 0x31ea) {
        generate_noise_1(m, substr);
        maxchan += 2;
    } else {
        generate_noise_2(m, substr);
    }

    for (mat = 0; mat < s->num_primitive_matrices; mat++) {
        unsigned int dest_ch = s->matrix_ch[mat];
        int quant_step_size = s->quant_step_size[dest_ch];

        /* TODO: DSPContext? */

        for (i = 0; i < s->blockpos; i++) {
            int64_t accum = 0;
            for (src_ch = 0; src_ch <= maxchan; src_ch++) {
                accum += (int64_t)m->sample_buffer[i][src_ch]
                                  * s->matrix_coeff[mat][src_ch];
            }
            if (s->matrix_noise_shift[mat]) {
                uint32_t index = s->num_primitive_matrices - mat;
                index = (i * (index * 2 + 1) + index) & (m->access_unit_size_pow2 - 1);
                accum += m->noise_buffer[index] << (s->matrix_noise_shift[mat] + 7);
            }
            m->sample_buffer[i][dest_ch] = ((accum >> 14) & ~((1 << quant_step_size) - 1))
                                             + m->bypassed_lsbs[i][mat];
        }
    }
}

/** Write the audio data into the output buffer.
 */

static int output_data_internal(MLPDecodeContext *m, unsigned int substr,
                                uint8_t *data, unsigned int *data_size, int is32)
{
    SubStream *s = &m->substream[substr];
    unsigned int i, ch = 0;
    int32_t *data_32 = (int32_t*) data;
    int16_t *data_16 = (int16_t*) data;

    if (*data_size < (s->max_channel + 1) * s->blockpos * (is32 ? 4 : 2))
        return -1;

    for (i = 0; i < s->blockpos; i++) {
        for (ch = 0; ch <= s->max_channel; ch++) {
            int32_t sample = m->sample_buffer[i][ch] << s->output_shift[ch];
            s->lossless_check_data ^= (sample & 0xffffff) << ch;
            if (is32) *data_32++ = sample << 8;
            else      *data_16++ = sample >> 8;
        }
    }

    *data_size = i * ch * (is32 ? 4 : 2);

    return 0;
}

static int output_data(MLPDecodeContext *m, unsigned int substr,
                       uint8_t *data, unsigned int *data_size)
{
    if (m->avctx->sample_fmt == SAMPLE_FMT_S32)
        return output_data_internal(m, substr, data, data_size, 1);
    else
        return output_data_internal(m, substr, data, data_size, 0);
}


/** XOR together all the bytes of a buffer.
 *  Does this belong in dspcontext?
 */

static uint8_t calculate_parity(const uint8_t *buf, unsigned int buf_size)
{
    uint32_t scratch = 0;
    const uint8_t *buf_end = buf + buf_size;

    for (; buf < buf_end - 3; buf += 4)
        scratch ^= *((const uint32_t*)buf);

    scratch ^= scratch >> 16;
    scratch ^= scratch >> 8;

    for (; buf < buf_end; buf++)
        scratch ^= *buf;

    return scratch;
}

/**
 * Read an access unit from the stream.
 * Returns -1 on error, 0 if not enough data is present in the input stream
 * otherwise returns the number of bytes consumed.
 */

static int read_access_unit(AVCodecContext *avctx, void* data, int *data_size,
                            const uint8_t *buf, int buf_size)
{
    MLPDecodeContext *m = avctx->priv_data;
    GetBitContext gb;
    unsigned int length, substr;
    unsigned int substream_start;
    unsigned int header_size = 4;
    unsigned int substr_header_size = 0;
    uint8_t substream_parity_present[MAX_SUBSTREAMS];
    uint16_t substream_data_len[MAX_SUBSTREAMS];
    uint8_t parity_bits;

    if (buf_size < 4)
        return 0;

    length = (AV_RB16(buf) & 0xfff) * 2;

    if (length > buf_size)
        return -1;

    init_get_bits(&gb, (buf + 4), (length - 4) * 8);

    if (show_bits_long(&gb, 31) == (0xf8726fba >> 1)) {
        dprintf(m->avctx, "Found major sync\n");
        if (read_major_sync(m, &gb) < 0)
            goto error;
        header_size += 28;
    }

    if (!m->params_valid) {
        av_log(m->avctx, AV_LOG_WARNING,
               "Stream parameters not seen; skipping frame\n");
        *data_size = 0;
        return length;
    }

    substream_start = 0;

    for (substr = 0; substr < m->num_substreams; substr++) {
        int extraword_present, checkdata_present, end;

        extraword_present = get_bits1(&gb);
        skip_bits1(&gb);
        checkdata_present = get_bits1(&gb);
        skip_bits1(&gb);

        end = get_bits(&gb, 12) * 2;

        substr_header_size += 2;

        if (extraword_present) {
            skip_bits(&gb, 16);
            substr_header_size += 2;
        }

        if (end + header_size + substr_header_size > length) {
            av_log(m->avctx, AV_LOG_ERROR,
                   "Substream %d data indicated length goes off end of packet.\n",
                   substr);
            end = length - header_size - substr_header_size;
        }

        if (end < substream_start) {
            av_log(avctx, AV_LOG_ERROR,
                   "Substream %d data indicated end offset "
                   "is before calculated start offset.\n",
                   substr);
            goto error;
        }

        if (substr > m->max_decoded_substream)
            continue;

        substream_parity_present[substr] = checkdata_present;
        substream_data_len[substr] = end - substream_start;
        substream_start = end;
    }

    parity_bits  = calculate_parity(buf, 4);
    parity_bits ^= calculate_parity(buf + header_size, substr_header_size);

    if ((((parity_bits >> 4) ^ parity_bits) & 0xF) != 0xF) {
        av_log(avctx, AV_LOG_ERROR, "Parity check failed.\n");
        goto error;
    }

    buf += header_size + substr_header_size;

    for (substr = 0; substr <= m->max_decoded_substream; substr++) {
        SubStream *s = &m->substream[substr];
        init_get_bits(&gb, buf, substream_data_len[substr] * 8);

        s->blockpos = 0;
        do {
            if (get_bits1(&gb)) {
                if (get_bits1(&gb)) {
                    /* A restart header should be present */
                    if (read_restart_header(m, &gb, buf, substr) < 0)
                        goto next_substr;
                    s->restart_seen = 1;
                }

                if (!s->restart_seen) {
                    av_log(m->avctx, AV_LOG_ERROR,
                           "No restart header present in substream %d.\n",
                           substr);
                    goto next_substr;
                }

                if (read_decoding_params(m, &gb, substr) < 0)
                    goto next_substr;
            }

            if (!s->restart_seen) {
                av_log(m->avctx, AV_LOG_ERROR,
                       "No restart header present in substream %d.\n",
                       substr);
                goto next_substr;
            }

            if (read_block_data(m, &gb, substr) < 0)
                return -1;

        } while ((get_bits_count(&gb) < substream_data_len[substr] * 8)
                 && get_bits1(&gb) == 0);

        skip_bits(&gb, (-get_bits_count(&gb)) & 15);
        if (substream_data_len[substr] * 8 - get_bits_count(&gb) >= 48 &&
            (show_bits_long(&gb, 32) == 0xd234d234 ||
             show_bits_long(&gb, 20) == 0xd234e)) {
            skip_bits(&gb, 18);
            if (substr == m->max_decoded_substream)
                av_log(m->avctx, AV_LOG_INFO, "End of stream indicated\n");

            if (get_bits1(&gb)) {
                int shorten_by = get_bits(&gb, 13);
                shorten_by = FFMIN(shorten_by, s->blockpos);
                s->blockpos -= shorten_by;
            } else
                skip_bits(&gb, 13);
        }
        if (substream_parity_present[substr]) {
            uint8_t parity, checksum;

            parity = calculate_parity(buf, substream_data_len[substr] - 2);
            if ((parity ^ get_bits(&gb, 8)) != 0xa9)
                av_log(m->avctx, AV_LOG_ERROR,
                       "Substream %d parity check failed\n", substr);

            checksum = mlp_checksum8(buf, substream_data_len[substr] - 2);
            if (checksum != get_bits(&gb, 8))
                av_log(m->avctx, AV_LOG_ERROR, "Substream %d checksum failed\n",
                       substr);
        }
        if (substream_data_len[substr] * 8 != get_bits_count(&gb)) {
            av_log(m->avctx, AV_LOG_ERROR, "Substream %d length mismatch.\n",
                   substr);
            return -1;
        }

next_substr:
        buf += substream_data_len[substr];
    }

    rematrix_channels(m, substr - 1);

    if (output_data(m, substr - 1, data, data_size) < 0)
        return -1;

    return length;

error:
    m->params_valid = 0;
    return -1;
}

AVCodec mlp_decoder = {
    "mlp",
    CODEC_TYPE_AUDIO,
    CODEC_ID_MLP,
    sizeof(MLPDecodeContext),
    mlp_decode_init,
    NULL,
    NULL,
    read_access_unit,
    .long_name = NULL_IF_CONFIG_SMALL("Meridian Lossless Packing"),
};

