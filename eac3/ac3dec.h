/*
 * Common code between AC3 and E-AC3 decoder
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

/**
 * @file ac3.h
 * Common code between AC3 encoder and decoder.
 */

#ifndef AC3DEC_H
#define AC3DEC_H

#include "libavutil/random.h"
#include "ac3tab.h"
#include "ac3.h"
#include "ac3dec_data.h"
#include "bitstream.h"
#include "dsputil.h"

#define TEST_SPX 0
#define TEST_ECPL 0

/* override ac3.h to include coupling channel */
#undef AC3_MAX_CHANNELS
#define AC3_MAX_CHANNELS 7

#define CPL_CH 0

#define AC3_OUTPUT_LFEON  8

#define AC3_MAX_COEFS   256
#define AC3_BLOCK_SIZE  256
#define MAX_BLOCKS 6
#define MAX_SPX_CODES 18

typedef struct AC3DecodeContext {
    AVCodecContext *avctx;  ///< Parent context
    GetBitContext gbc;      ///< Bitstream reader
    uint8_t *input_buffer;  ///< temp buffer to prevent overread

///@defgroup bsi Bit Stream Information
///@{
    int frame_type;         ///< Frame type (strmtyp)
    int substreamid;        ///< Substream identification
    int frame_size;         ///< Frame size, in bytes
    int bit_rate;           ///< Bitrate, in bits-per-second
    int sample_rate;        ///< Sample rate, in Hz
    int num_blocks;         ///< Number of audio blocks
    int channel_mode;       ///< Channel mode (acmod)
    int lfe_on;             ///< Low frequency effect channel on (lfeon)
    int channel_map;        ///< Custom channel map
    int bitstream_id;       ///< Bit stream identification (bsid)
    int center_mix_level;   ///< Center mix level index
    int surround_mix_level; ///< Surround mix level index
    int eac3;               ///< indicates if current frame is E-AC3
///@}

///@defgroup audfrm Frame Syntax Parameters
    int snr_offset_strategy;    ///< SNR offset strategy (snroffststr)
    int block_switch_syntax;    ///< Block switch syntax enabled (blkswe)
    int dither_flag_syntax;     ///< Dither flag syntax enabled (dithflage)
    int bit_allocation_syntax;  ///< Bit allocation model syntax enabled (bamode)
    int fast_gain_syntax;       ///< Fast gain codes enabled (frmfgaincode)
    int dba_syntax;             ///< Delta bit allocation syntax enabled (dbaflde)
    int skip_syntax;            ///< Skip Filed syntax enabled (skipflde)
///@}

///@defgroup cpl Standard Coupling
    int cpl_in_use[MAX_BLOCKS];                 ///< Coupling in use (cplinu)
    int cpl_strategy_exists[MAX_BLOCKS];        ///< Coupling strategy exists (cplstre)
    int channel_in_cpl[AC3_MAX_CHANNELS];       ///< Channel in coupling (chincpl)
    int phase_flags_in_use;                     ///< Phase flag in use (phsflginu)
    int phase_flags[18];                        ///< Phase flag
    int num_cpl_subbands;                       ///< Number of coupling sub bands (ncplsubnd)
    int num_cpl_bands;                          ///< Number of coupling bands (ncplbnd)
    int cpl_band_struct[18];                    ///< Coupling band structure (cplbndstrc)
    int firstchincpl;                           ///< First channel in coupling
    int first_cpl_coords[AC3_MAX_CHANNELS];     ///< First coupling coordinates states (firstcplcos)
    int cpl_coords[AC3_MAX_CHANNELS][18];       ///< coupling coordinates (cplco)
///@}

///@defgroup aht Adaptive Hybrid Transform
    int channel_uses_aht[AC3_MAX_CHANNELS];     ///< Channel AHT in use (chahtinu)
    int pre_mantissa[6][AC3_MAX_CHANNELS][256]; ///< Pre-IDCT mantissas
///@}

#if TEST_SPX
///@defgroup spx Spectral Extension
    int channel_uses_spx[AC3_MAX_CHANNELS]; ///< Channel in spectral extension attenuation process (chinspxatten)
    int spx_atten_code[AC3_MAX_CHANNELS];   ///< spectral extension attenuation code (spxattencod)
    int spx_in_use;                         ///< spectral extension in use
    int chinspx[AC3_MAX_CHANNELS];          ///< Channel in spectral extension
    int spxstrtf;                           ///< Spectral extension start copy frequency code
    int spxbegf;                            ///< Spectral extension begin frequency code
    int spxendf;                            ///< Spectral extension end frequency code
    int nspxbnds;                           ///< Number of structured spectral extension bands
    int spxbndsztab[MAX_SPX_CODES];         ///< Sizes of spectral extension bands
    int spxbndstrc[MAX_SPX_CODES];          ///< Spectral extension band structure
    int spxcoe[AC3_MAX_CHANNELS];           ///< Spectral extension coordinates exists
    int spxblnd[AC3_MAX_CHANNELS];          ///< Spectral extension blend
    int firstspxcos[AC3_MAX_CHANNELS];      ///< First spectral extension coordinates states
    float spxco[AC3_MAX_CHANNELS][18];      ///< Spectral extension coordinates
///@}
#endif

#if TEST_ECPL
///@defgroup ecpl Enhanced Coupling
    int ecpl_in_use;                        ///< Enhanced coupling in use
    int ecplbegf;                           ///< Enhanced coupling begin frequency code
    int ecplendf;                           ///< Enhanced coupling end frequency code
    int ecpl_start_subbnd;                  ///< Enhanced coupling begin frequency
    int ecpl_end_subbnd;                    ///< Enhanced coupling end frequency
    int necplbnd;                           ///< Number of structured enhanced coupling bands
    int ecplbndstrc[23];                    ///< Enhanced coupling band structure
    int ecplangleintrp;                     ///< Enhanced coupling angle interpolation flag
    int ecplparam1e[AC3_MAX_CHANNELS];      ///< Enhanced coupling parameters 1 exists
    int ecplparam2e[AC3_MAX_CHANNELS];      ///< Enhanced coupling parameters 2 exists
    int ecplamp[AC3_MAX_CHANNELS][23];      ///< Enhanced coupling amplitude scaling
    int ecplangle[AC3_MAX_CHANNELS][23];    ///< Enhanced coupling angle
    int ecplchaos[AC3_MAX_CHANNELS][23];    ///< Enhanced coupling chaos
    int ecpltrans[AC3_MAX_CHANNELS];        ///< Enhanced coupling transient present
///@}
#endif

///@defgroup channel Channel
    int fbw_channels;                           ///< Number of fbw channels
    int channels;                               ///< Total of all channels
    int lfe_ch;                                 ///< Index of LFE channel
    float downmix_coeffs[AC3_MAX_CHANNELS][2];  ///< stereo downmix coefficients
    float downmix_coeff_adjust[2];              ///< adjustment needed for each output channel when downmixing
    int downmixed;                              ///< indicates if transform coeffs are currently downmixed
    int output_mode;                            ///< output channel configuration
    int out_channels;                           ///< number of output channels
///@}

///@defgroup dynrng Dynamic Range
    float dynamic_range[2]; ///< Dynamic range gain (dynrng)
///@}

///@defgroup bandwidth Bandwidth
    int start_freq[AC3_MAX_CHANNELS];   ///< Start frequency bin (strtmant)
    int end_freq[AC3_MAX_CHANNELS];     ///< End frequency bin (endmant)
///@}

///@defgroup rematrixing Rematrixing
    int num_rematrixing_bands;  ///< Number of rematrixing bands (nrematbnds)
    int rematrixing_flags[4];   ///< Rematrixing flags (rematflg)
///@}

///@defgroup exponents Exponents
    int num_exp_groups[AC3_MAX_CHANNELS];           ///< Number of exponent groups
    uint8_t dexps[AC3_MAX_CHANNELS][AC3_MAX_COEFS]; ///< Differential exponents
    int exp_strategy[MAX_BLOCKS][AC3_MAX_CHANNELS]; ///< Channel exponent strategy (chexpstr)
///@}

///@defgroup bitalloc Bit Allocation
    AC3BitAllocParameters bit_alloc_params;         ///< Bit allocation parameters
    int first_cpl_leak;                             ///< First coupling leak state (firstcplleak)
    int snr_offset[AC3_MAX_CHANNELS];               ///< SNR offset (snroffst)
    int fast_gain[AC3_MAX_CHANNELS];                ///< Channel fast gain (fgain)
    uint8_t bap[AC3_MAX_CHANNELS][AC3_MAX_COEFS];   ///< bit allocation pointers
    int16_t psd[AC3_MAX_CHANNELS][AC3_MAX_COEFS];   ///< scaled exponents
    int16_t band_psd[AC3_MAX_CHANNELS][50];         ///< interpolated exponents (bndpsd)
    int16_t mask[AC3_MAX_CHANNELS][50];             ///< masking values
    uint8_t dba_mode[AC3_MAX_CHANNELS];             ///< Delta bit allocation mode (deltbae)
    uint8_t dba_nsegs[AC3_MAX_CHANNELS];            ///< Number of delta segments (deltnseg)
    uint8_t dba_offsets[AC3_MAX_CHANNELS][9];       ///< Delta segment offsets (deltoffst)
    uint8_t dba_lengths[AC3_MAX_CHANNELS][9];       ///< Delta segment lengths (deltlen)
    uint8_t dba_values[AC3_MAX_CHANNELS][9];        ///< Delta values for each segment (deltba)
///@}

///@defgroup dithering Zero-Mantissa Dithering
    int dither_all;                     ///< Indicates whether all channels use dithering
    int dither_flag[AC3_MAX_CHANNELS];  ///< Dither flag (dithflag)
    AVRandomState dith_state;           ///< for dither generation
///@}

///@defgroup imdct IMDCT
    int block_switch[AC3_MAX_CHANNELS]; ///< Block switch flag (blksw)
    MDCTContext imdct_512;              ///< for 512 sample imdct transform
    MDCTContext imdct_256;              ///< for 256 sample imdct transform
///@}

///@defgroup opt Optimization
    DSPContext  dsp;    ///< for optimization
    float add_bias;     ///< offset for float_to_int16 conversion
    float mul_bias;     ///< scaling for float_to_int16 conversion
///@}

    int fixed_coeffs[AC3_MAX_CHANNELS][256];    ///< fixed-point transform coefficients

///@defgroup arrays Aligned Arrays
    DECLARE_ALIGNED_16(float, transform_coeffs[AC3_MAX_CHANNELS][AC3_MAX_COEFS]);   ///< Frequency Coefficients
    DECLARE_ALIGNED_16(float, delay[AC3_MAX_CHANNELS][AC3_BLOCK_SIZE]);             ///< delay - added to the next block
    DECLARE_ALIGNED_16(float, window[AC3_BLOCK_SIZE]);                              ///< window coefficients
    DECLARE_ALIGNED_16(float, tmp_output[AC3_BLOCK_SIZE * 24]);                     ///< temp storage for output before windowing
    DECLARE_ALIGNED_16(float, tmp_imdct[AC3_BLOCK_SIZE * 24]);                      ///< temp storage for imdct transform
    DECLARE_ALIGNED_16(float, output[AC3_MAX_CHANNELS][AC3_BLOCK_SIZE]);            ///< output after imdct transform and windowing
    DECLARE_ALIGNED_16(int16_t, int_output[AC3_MAX_CHANNELS][AC3_BLOCK_SIZE]);      ///< final 16-bit integer output
///@}
} AC3DecodeContext;

/* TEMPORARY SOLUTION */
int ff_eac3_parse_header(AC3DecodeContext *s);
void ff_eac3_get_transform_coeffs_aht_ch(AC3DecodeContext *s, int ch);
void ff_eac3_idct_transform_coeffs_ch(AC3DecodeContext *s, int ch, int blk);
void ff_eac3_tables_init(void);
void ff_eac3_log_missing_feature(AVCodecContext *avctx, const char *log);

#endif /* AC3DEC_H */
