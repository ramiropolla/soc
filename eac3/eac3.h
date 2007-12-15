/*
 * E-AC3 parser
 * Copyright (c) 2007 Bartlomiej Wolowiec <bartek.wolowiec@gmail.com>
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
#ifndef FFMPEG_EAC3_H
#define FFMPEG_EAC3_H

#include "dsputil.h"
#include "avcodec.h"
#include "ac3.h"
#include "random.h"
#include "bitstream.h"

/* override ac3.h to include coupling channel */
#undef AC3_MAX_CHANNELS
#define AC3_MAX_CHANNELS 7

#define AC3_MAX_COEFS   256

#define AC3_BLOCK_SIZE  256
#define MAX_BLOCKS 6
#define MAX_SPX_CODES 18

typedef struct EAC3Context{
    AVCodecContext *avctx;           ///< Parent context
    int syncword;
///@name Bit stream information
///@{
    int stream_type;                 ///< Stream type (strmtyp)
    int substreamid;                 ///< Substream identification
    int frame_size;                  ///< Frame size, in bytes
    int sr_code;                     ///< Sample rate code (fscod)
    int sr_code2;                    ///< Sample rate code 2 (fscod2)
    int num_blocks;                  ///< Number of audio blocks
    int channel_mode;                ///< Channel mode (acmod)
    int lfe_on;                      ///< Low frequency effect channel on (lfeon)
    int bitstream_id;                ///< Bit stream identification (bsid)
    float dialog_norm[2];            ///< Dialogue normalization (dialnorm)
///@}
///@name Audio Frame
///@{
    int snr_offset_strategy;         ///< SNR offset strategy (snroffststr)
    int snr_offset[AC3_MAX_CHANNELS];///< SNR offset (snroffst)
    int block_switch_syntax;         ///< Block switch syntax enabled (blkswe)
    int dither_flag_syntax;          ///< Dither flag syntax enabled (dithflage)
    int bit_allocation_syntax;       ///< Bit allocation model syntax enabled (bamode)
    int fast_gain_syntax;            ///< Fast gain codes enabled (frmfgaincode)
    int dba_syntax;                  ///< Delta bit allocation syntax enabled (dbaflde)
    int skip_syntax;                 ///< Skip Filed syntax enabled (skipflde)
    int cpl_in_use[MAX_BLOCKS];      ///< Coupling in use (cplinu)
    int cpl_stratety_exists[MAX_BLOCKS];            ///< Coupling strategy exists (cplstre)
    int exp_strategy[MAX_BLOCKS][AC3_MAX_CHANNELS]; ///< Channel exponent strategy (chexpstr)
    int channel_uses_aht[AC3_MAX_CHANNELS];         ///< Channel AHT in use (chahtinu)
    int channel_uses_spx[AC3_MAX_CHANNELS];         ///< Channel in spectral extension attenuation process (chinspxatten)
    int spx_atten_code[AC3_MAX_CHANNELS];           ///< spectral extension attenuation code (spxattencod)
///@}
///@name Audio block
///@{
    int blksw[AC3_MAX_CHANNELS];     ///< Block switch flag
    int dithflag[AC3_MAX_CHANNELS];  ///< Dither flag
    float dynamic_range[2];          ///< Dynamic range gain (dynrng)
    float downmix_coeffs[AC3_MAX_CHANNELS][2];  ///< stereo downmix coefficients
    int spxinu;                      ///< spectral extension in use
    int chinspx[AC3_MAX_CHANNELS];   ///< Channel in spectral extension
    int spxstrtf;                    ///< Spectral extension start copy frequency code
    int spxbegf;                     ///< Spectral extension begin frequency code
    int spxendf;                     ///< Spectral extension end frequency code
    int spxbndstrc[MAX_SPX_CODES];   ///< Spectral extension band structure
    int spxcoe[AC3_MAX_CHANNELS];    ///< Spectral extension coordinates exists
    int spxblnd[AC3_MAX_CHANNELS];   ///< Spectral extension blend
    int ecpl_in_use;                 ///< Enhanced coupling in use
    int channel_in_cpl[AC3_MAX_CHANNELS];   ///< Channel in coupling (chincpl)
    int phase_flags_in_use;                 ///< Phase flag in use (phsflginu)
    int cpl_band_struct[18];         ///< Coupling band structure (cplbndstrc)
    int ecplbegf;                    ///< Enhanced coupling begin frequency code
    int ecplendf;                    ///< Enhanced coupling end frequency code
    int ecplbndstrc[23];             ///< Enhanced coupling band structure
    int phase_flags[18];             ///< Phase flag
    int ecplangleintrp;              ///< Enhanced coupling angle interpolation flag
    int ecplparam1e[AC3_MAX_CHANNELS];   ///< Enhanced coupling parameters 1 exists
    int ecplparam2e[AC3_MAX_CHANNELS];   ///< Enhanced coupling parameters 2 exists
    int ecplamp[AC3_MAX_CHANNELS][23];   ///< Enhanced coupling amplitude scaling
    int ecplangle[AC3_MAX_CHANNELS][23]; ///< Enhanced coupling angle
    int ecplchaos[AC3_MAX_CHANNELS][23]; ///< Enhanced coupling chaos
    int ecpltrans[AC3_MAX_CHANNELS];     ///< Enhanced coupling transient present
    int rematflg[4];                 ///< Rematrixing flag
    int cplabsexp;                   ///< Coupling absolute exponent

    int fgain[AC3_MAX_CHANNELS];     ///< Channel fast gain
    int cplleake;                    ///< Coupling leak initialization exists
    uint8_t deltbae[AC3_MAX_CHANNELS];   ///< Delta bit allocation exists
    uint8_t deltnseg[AC3_MAX_CHANNELS];  ///< Channel delta bit allocation number of segments
    uint8_t deltoffst[AC3_MAX_CHANNELS][9]; ///< Channel delta bit allocation offset
    uint8_t deltlen[AC3_MAX_CHANNELS][9];   ///< Channel delta bit allocation length
    uint8_t deltba[AC3_MAX_CHANNELS][9];    ///< Channel delta bit allocation

    int got_cplchan;
    int chgaqgain[256];                             ///< Channel gain adaptive quantization gain
    float pre_chmant[6][AC3_MAX_CHANNELS][256];     ///< Pre channel mantissas

    int firstspxcos[AC3_MAX_CHANNELS];              ///< First spectral extension coordinates states
    int first_cpl_coords[AC3_MAX_CHANNELS];         ///< First coupling coordinates states (firstcplcos)
    int first_cpl_leak;                             ///< First coupling leak state (firstcplleak)
///@}

    int nrematbnds;                    ///< Number of rematrixing bands
    int num_cpl_subbands;              ///< Number of coupling sub bands (ncplsubnd)
    int num_cpl_bands;                 ///< Number of coupling bands (ncplbnd)

    int nchgrps[AC3_MAX_CHANNELS];                  ///< Number of fbw channel exponent groups
    uint8_t dexps[AC3_MAX_CHANNELS][AC3_MAX_COEFS]; ///< Differential exponents

    int start_freq[AC3_MAX_CHANNELS];   ///< Start frequency bin (strtmant)
    int end_freq[AC3_MAX_CHANNELS];     ///< End frequency bin (endmant)
    int firstchincpl;
    int ecpl_start_subbnd;             ///< Enhanced coupling begin frequency
    int ecpl_end_subbnd;               ///< Enhanced coupling end frequency

    int necplbnd;                      ///< Number of structured enhanced coupling bands
    int nspxbnds;                      ///< Number of structured spectral extension bands
    int spxbndsztab[MAX_SPX_CODES];    ///< Sizes of spectral extension bands
    int fbw_channels;                  ///< Number of fbw channels

    uint8_t bap[AC3_MAX_CHANNELS][AC3_MAX_COEFS];   ///< bit allocation pointers
    uint8_t hebap[AC3_MAX_CHANNELS][AC3_MAX_COEFS]; ///< high-efficiency bit allocation pointers for AHT
    int16_t psd[AC3_MAX_CHANNELS][AC3_MAX_COEFS];   ///< scaled exponents
    int16_t bndpsd[AC3_MAX_CHANNELS][50];           ///< interpolated exponents
    int16_t mask[AC3_MAX_CHANNELS][50];             ///< masking values

    float   cpl_coords[AC3_MAX_CHANNELS][18];       ///< coupling coordinates (cplco)
    float   spxco[AC3_MAX_CHANNELS][18];            ///< Spectral extension coordinates

    AC3BitAllocParameters bit_alloc_params;         ///< Bit allocation parameters

    AVRandomState dith_state;        ///< for dither generation

    int num_channels;                ///< Total of all channels
    int lfe_channel;                 ///< Index of LFE channel

    GetBitContext *gbc;              ///< Bitstream reader

    MDCTContext imdct_512;           ///< for 512 sample imdct transform
    MDCTContext imdct_256;           ///< for 256 sample imdct transform
    DSPContext  dsp;                 ///< for optimization

    DECLARE_ALIGNED_16(float, transform_coeffs[AC3_MAX_CHANNELS][AC3_MAX_COEFS]);
    DECLARE_ALIGNED_16(float, delay[AC3_MAX_CHANNELS][AC3_BLOCK_SIZE]);  ///< delay - added to the next block
    DECLARE_ALIGNED_16(float, window[AC3_BLOCK_SIZE]);               ///< window coefficients
    DECLARE_ALIGNED_16(float, tmp_output[AC3_BLOCK_SIZE * 24]);      ///< temp storage for output before windowing
    DECLARE_ALIGNED_16(float, tmp_imdct[AC3_BLOCK_SIZE * 24]);       ///< temp storage for imdct transform
    DECLARE_ALIGNED_16(float, output[AC3_MAX_CHANNELS][AC3_BLOCK_SIZE]); ///< output after imdct transform and windowing
    DECLARE_ALIGNED_16(int16_t, int_output[AC3_MAX_CHANNELS][AC3_BLOCK_SIZE]);///< final 16-bit integer output


    float add_bias;                  ///< offset for float_to_int16 conversion
    float mul_bias;                  ///< scaling for float_to_int16 conversion
}EAC3Context;

/** Channel gain adaptive quantization mode */
typedef enum {
    EAC3_GAQ_NO =0,
    EAC3_GAQ_12,
    EAC3_GAQ_14,
    EAC3_GAQ_124
} EAC3GaqMode;

/** Stream Type */
typedef enum {
    EAC3_STREAM_TYPE_INDEPENDENT = 0,
    EAC3_STREAM_TYPE_DEPENDENT,
    EAC3_STREAM_TYPE_AC3_CONVERT,
    EAC3_STREAM_TYPE_RESERVED
} EAC3StreamType;

#define EAC3_SR_CODE_REDUCED  3

#endif
