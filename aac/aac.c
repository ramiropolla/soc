/*
 * AAC decoder
 * Copyright (c) 2005-2006 Oded Shimon ( ods15 ods15 dyndns org )
 * Copyright (c) 2006-2007 Maxim Gavrilov ( maxim.gavrilov gmail com )
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
 * @file aac.c
 * AAC decoder
 * @author Oded Shimon  ( ods15 ods15 dyndns org )
 * @author Maxim Gavrilov ( maxim.gavrilov gmail com )
 */

/**
 * AAC SSR (Scalable Sample Rate) is currently not working, and therefore
 * not compiled in. SSR files play without crashing but produce audible
 * artifacts that seem to be related to EIGHT_SHORT_SEQUENCE windows.
 */
//#define AAC_SSR

/**
 * AAC LTP (Long Term Prediction) is currently not working, and therefore
 * not compiled in. Playing LTP files with LTP support compiled in results
 * in crashes due to SSE alignment issues. Also, there are major audible
 * artifacts.
 */
//#define AAC_LTP

/*
 * supported tools
 *
 * Support?             Name
 * N (code in SoC repo) gain control
 * Y                    block switching
 * Y                    window shapes - standard
 * N                    window shapes - Low Delay
 * Y                    filterbank - standard
 * N (code in SoC repo) filterbank - Scalable Sample Rate
 * Y                    Temporal Noise Shaping
 * N (code in SoC repo) Long Term Prediction
 * Y                    intensity stereo
 * Y                    channel coupling
 * N                    frequency domain prediction
 * Y                    Perceptual Noise Substitution
 * Y                    Mid/Side stereo
 * N                    Scalable Inverse AAC Quantization
 * N                    Frequency Selective Switch
 * N                    upsampling filter
 * Y                    quantization & coding - AAC
 * N                    quantization & coding - TwinVQ
 * N                    quantization & coding - BSAC
 * N                    AAC Error Resilience tools
 * N                    Error Resilience payload syntax
 * N                    Error Protection tool
 * N                    CELP
 * N                    Silence Compression
 * N                    HVXC
 * N                    HVXC 4kbits/s VR
 * N                    Structured Audio tools
 * N                    Structured Audio Sample Bank Format
 * N                    MIDI
 * N                    Harmonic and Individual Lines plus Noise
 * N                    Text-To-Speech Interface
 * N (in progress)      Spectral Band Replication
 * Y (not in this code) Layer-1
 * Y (not in this code) Layer-2
 * Y (not in this code) Layer-3
 * N                    SinuSoidal Coding (Transient, Sinusoid, Noise)
 * N (planned)          Parametric Stereo
 * N                    Direct Stream Transfer
 *
 * Note: - HE AAC v1 comprises LC AAC with Spectral Band Replication.
 *       - HE AAC v2 comprises LC AAC with Spectral Band Replication and
           Parametric Stereo.
 */


#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"
#include "libavutil/random.h"

#include "aactab.h"
#include "mpeg4audio.h"

#include <assert.h>
#include <math.h>
#include <string.h>

#define MAX_CHANNELS 64
#define MAX_TAGID 16

#ifndef CONFIG_HARDCODED_TABLES
    static float ivquant_tab[IVQUANT_SIZE];
    static float pow2sf_tab[316];
#endif /* CONFIG_HARDCODED_TABLES */
DECLARE_ALIGNED_16(static float,  kbd_long_1024[1024]);
DECLARE_ALIGNED_16(static float,  kbd_short_128[128]);
DECLARE_ALIGNED_16(static float, sine_long_1024[1024]);
DECLARE_ALIGNED_16(static float, sine_short_128[128]);

static VLC mainvlc;
static VLC books[11];

enum AudioObjectType {
    AOT_NULL,
                               // Support?                Name
    AOT_AAC_MAIN,              ///< Y                       Main
    AOT_AAC_LC,                ///< Y                       Low Complexity
    AOT_AAC_SSR,               ///< N (code in SoC repo)    Scalable Sample Rate
    AOT_AAC_LTP,               ///< N (code in SoC repo)    Long Term Prediction
    AOT_SBR,                   ///< N (in progress)         Spectral Band Replication
    AOT_AAC_SCALABLE,          ///< N                       Scalable
    AOT_TWINVQ,                ///< N                       Twin Vector Quantizer
    AOT_CELP,                  ///< N                       Code Excited Linear Prediction
    AOT_HVXC,                  ///< N                       Harmonic Vector eXcitation Coding
    AOT_TTSI             = 12, ///< N                       Text-To-Speech Interface
    AOT_MAINSYNTH,             ///< N                       Main Synthesis
    AOT_WAVESYNTH,             ///< N                       Wavetable Synthesis
    AOT_MIDI,                  ///< N                       General MIDI
    AOT_SAFX,                  ///< N                       Algorithmic Synthesis and Audio Effects
    AOT_ER_AAC_LC,             ///< N                       Error Resilient Low Complexity
    AOT_ER_AAC_LTP       = 19, ///< N                       Error Resilient Long Term Prediction
    AOT_ER_AAC_SCALABLE,       ///< N                       Error Resilient Scalable
    AOT_ER_TWINVQ,             ///< N                       Error Resilient Twin Vector Quantizer
    AOT_ER_BSAC,               ///< N                       Error Resilient Bit-Sliced Arithmetic Coding
    AOT_ER_AAC_LD,             ///< N                       Error Resilient Low Delay
    AOT_ER_CELP,               ///< N                       Error Resilient Code Excited Linear Prediction
    AOT_ER_HVXC,               ///< N                       Error Resilient Harmonic Vector eXcitation Coding
    AOT_ER_HILN,               ///< N                       Error Resilient Harmonic and Individual Lines plus Noise
    AOT_ER_PARAM,              ///< N                       Error Resilient Parametric
    AOT_SSC,                   ///< N                       SinuSoidal Coding
};

enum RawDataBlockID {
    ID_SCE,
    ID_CPE,
    ID_CCE,
    ID_LFE,
    ID_DSE,
    ID_PCE,
    ID_FIL,
    ID_END,
};

enum ExtensionPayloadID {
    EXT_FILL,
    EXT_FILL_DATA,
    EXT_DATA_ELEMENT,
    EXT_DYNAMIC_RANGE = 0xb,
    EXT_SBR_DATA      = 0xd,
    EXT_SBR_DATA_CRC  = 0xe,
};

enum WindowSequence {
    ONLY_LONG_SEQUENCE,
    LONG_START_SEQUENCE,
    EIGHT_SHORT_SEQUENCE,
    LONG_STOP_SEQUENCE,
};

enum BandType {
    ZERO_BT        = 0,     ///< Scalefactors and spectral data are all zero.
    FIRST_PAIR_BT  = 5,     ///< This and later band types encode two values (rather than four) with one code word.
    ESC_BT         = 11,    ///< Spectral data are coded with an escape sequence.
    NOISE_BT       = 13,    ///< Spectral data are scaled white noise not coded in the bitstream.
    INTENSITY_BT2  = 14,    ///< Scalefactor data are intensity stereo positions.
    INTENSITY_BT   = 15,    ///< Scalefactor data are intensity stereo positions.
};

#define IS_CODEBOOK_UNSIGNED(x) ((x - 1) & 10)

enum ChannelType {
    AAC_CHANNEL_FRONT = 1,
    AAC_CHANNEL_SIDE  = 2,
    AAC_CHANNEL_BACK  = 3,
    AAC_CHANNEL_LFE   = 4,
    AAC_CHANNEL_CC    = 5,
};

/**
 * mix-down channel types
 * MIXDOWN_CENTER is the index into the mix-down arrays for a Single Channel Element with AAC_CHANNEL_FRONT.
 * MIXDOWN_(BACK|FRONT) are the indices for Channel Pair Elements with AAC_CHANNEL_(BACK|FRONT).
 */
enum {
    MIXDOWN_CENTER,
    MIXDOWN_FRONT,
    MIXDOWN_BACK,
};

/**
 * Program configuration - describes how channels are arranged. Either read from
 * stream (ID_PCE) or created based on a default fixed channel arrangement.
 */
typedef struct {
    enum ChannelType che_type[4][MAX_TAGID]; ///< channel element type with the first index as the first 4 raw_data_block IDs
    int mono_mixdown_tag;                    ///< The SCE tag to use if user requests mono   output, -1 if not available.
    int stereo_mixdown_tag;                  ///< The CPE tag to use if user requests stereo output, -1 if not available.
    int mixdown_coeff_index;                 ///< 0-3
    int pseudo_surround;                     ///< Mix surround channels out of phase.
} ProgramConfig;

#ifdef AAC_LTP
/**
 * Long Term Prediction
 */
#define MAX_LTP_LONG_SFB 40
typedef struct {
    int present;
    int lag;
    float coef;
    int used[MAX_LTP_LONG_SFB];
} LongTermPrediction;
#endif /* AAC_LTP */

/**
 * Individual Channel Stream
 */
typedef struct {
    int intensity_present;
    uint8_t max_sfb;            ///< number of scalefactor bands per group
    enum WindowSequence window_sequence[2];
    uint8_t use_kb_window[2];   ///< If set, use Kaiser-Bessel window, otherwise use a sinus window.
    int num_window_groups;
    uint8_t group_len[8];
#ifdef AAC_LTP
    LongTermPrediction ltp;
    LongTermPrediction ltp2;
#endif /* AAC_LTP */
    const uint16_t *swb_offset; ///< table of offsets to the lowest spectral coefficient of a scalefactor band, sfb, for a particular window
    int num_swb;                ///< number of scalefactor window bands
    int num_windows;
    int tns_max_bands;
} IndividualChannelStream;

/**
 * Temporal Noise Shaping
 */
typedef struct {
    int present;
    int n_filt[8];
    int length[8][4];
    int direction[8][4];
    int order[8][4];
    const float *tmp2_map[8][4];
    int coef[8][4][TNS_MAX_ORDER];
} TemporalNoiseShaping;

/**
 * M/S joint channel coding
 */
typedef struct {
    int present;
    uint8_t mask[8][64];
} MidSideStereo;

/**
 * Dynamic Range Control - decoded from the bitstream but not processed further.
 */
typedef struct {
    int pce_instance_tag;                           ///< Indicates with which program the DRC info is associated.
    int dyn_rng_sgn[17];                            ///< DRC sign information; 0 - positive, 1 - negative
    int dyn_rng_ctl[17];                            ///< DRC magnitude information
    int exclude_mask[MAX_CHANNELS];                 ///< Channels to be excluded from DRC processing.
    int additional_excluded_chns[MAX_CHANNELS];     /**< The exclude_mask bits are
                                                        coded in groups of 7 with 1 bit preceeding each group (except the first)
                                                        indicating that 7 more mask bits are coded. */
    int band_incr;                                  ///< Number of DRC bands greater than 1 having DRC info.
    int interpolation_scheme;                       ///< Indicates the interpolation scheme used in the SBR QMF domain.
    int band_top[17];                               ///< Indicates the top of the i-th DRC band in units of 4 spectral lines.
    int prog_ref_level;                             /**< A reference level for the long-term program audio level for all
                                                        channels combined. */
} DynamicRangeControl;

/**
 * pulse tool
 */
typedef struct {
    int num_pulse;
    int start;
    int offset[4];
    int amp[4];
} Pulse;

#ifdef AAC_SSR
/**
 * parameters for the SSR Inverse Polyphase Quadrature Filter
 */
typedef struct {
    float q[4][4];
    float t0[4][12];
    float t1[4][12];
} ssr_context;

/**
 * per-element gain control for SSR
 */
typedef struct {
    int max_band;
    int adjust_num[4][8];
    int alev[4][8][8];
    int aloc[4][8][8];
    float buf[4][24];
} ScalableSamplingRate;
#endif /* AAC_SSR */

/**
 * coupling parameters
 */
typedef struct {
    int is_indep_coup;     ///< Set if independent coupling (i.e. after IMDCT).
    int domain;            ///< Controls if coupling is performed before (0) or after (1) the TNS decoding of the target channels.
    int num_coupled;       ///< number of target elements
    int is_cpe[9];         ///< Set if target is an CPE (otherwise it's an SCE).
    int tag_select[9];     ///< element tag index
    int l[9];              ///< Apply gain to left channel of a CPE.
    int r[9];              ///< Apply gain to right channel of a CPE.
    float gain[18][8][64];
} ChannelCoupling;


/**
 * Single Channel Element - used for both SCE and LFE elements.
 */
typedef struct {
    float mixing_gain;                        /**< Channel gain (not used by AAC bitstream).
                                               *   Note that this is applied before joint stereo decoding.
                                               *   Thus, when used inside CPE elements, both channels must have equal gain.
                                               */
    IndividualChannelStream ics;
    TemporalNoiseShaping tns;
    enum BandType band_type[8][64];           ///< band types
    int band_type_run_end[8][64];             ///< band type run end points
    float sf[8][64];                          ///< scalefactors
    DECLARE_ALIGNED_16(float, coeffs[1024]);  ///< coefficients for IMDCT
    DECLARE_ALIGNED_16(float, saved[1024]);   ///< overlap
    DECLARE_ALIGNED_16(float, ret[1024]);     ///< PCM output
#ifdef AAC_LTP
    int16_t *ltp_state;
#endif /* AAC_LTP */
#ifdef AAC_SSR
    ScalableSamplingRate *ssr;
#endif /* AAC_SSR */
} SingleChannelElement;

/**
 * channel element - generic struct for SCE/CPE/CCE/LFE
 */
typedef struct {
    // CPE specific
    MidSideStereo ms;
    // shared
    SingleChannelElement ch[2];
    // CCE specific
    ChannelCoupling coup;
} ChannelElement;

/**
 * main AAC context
 */
typedef struct {
    AVCodecContext * avccontext;

    MPEG4AudioConfig m4ac;

    int is_saved;                 ///< Set if elements have stored overlap from previous frame.
    DynamicRangeControl che_drc;

    /**
     * @defgroup elements
     * @{
     */
    ProgramConfig pcs;
    ChannelElement * che[4][MAX_TAGID];
    /** @} */

    /**
     * @defgroup temporary aligned temporary buffers (We do not want to have these on the stack.)
     * @{
     */
    DECLARE_ALIGNED_16(float, buf_mdct[2048]);
    DECLARE_ALIGNED_16(float, revers[1024]);
    /** @} */

    /**
     * @defgroup tables   Computed / set up during initialization.
     * @{
     */
    MDCTContext mdct;
    MDCTContext mdct_small;
#ifdef AAC_LTP
    MDCTContext mdct_ltp;
#endif /* AAC_LTP */
    DSPContext dsp;
#ifdef AAC_SSR
    ssr_context ssrctx;
#endif /* AAC_SSR */
    AVRandomState random_state;
    /** @} */

    /**
     * @defgroup output   Members used for output interleaving and down-mixing.
     * @{
     */
    float *interleaved_output;                        ///< Interim buffer for interleaving PCM samples.
    float *output_data[MAX_CHANNELS];                 ///< Points to each element's 'ret' buffer (PCM output).
    ChannelElement *mm[3];                            ///< Center/Front/Back channel elements to use for matrix mix-down.
    float add_bias;                                   ///< offset for dsp.float_to_int16
    float sf_scale;                                   ///< Pre-scale for correct IMDCT and dsp.float_to_int16.
    int sf_offset;                                    ///< offset into pow2sf_tab as appropriate for dsp.float_to_int16
    /** @} */

} AACContext;


// TODO: Maybe add to dsputil?!
#if defined(AAC_LTP) || defined(AAC_SSR)
static void vector_fmul_dst(AACContext * ac, float * dst, const float * src0, const float * src1, int len) {
    memcpy(dst, src0, len * sizeof(float));
    ac->dsp.vector_fmul(dst, src1, len);
}
#endif

#if 0
static void vector_fmul_add_add_add(AACContext * ac, float * dst, const float * src0, const float * src1,
        const float * src2, const float * src3, float src4, int len) {
    int i;
    ac->dsp.vector_fmul_add_add(dst, src0, src1, src2, src4, len, 1);
    for (i = 0; i < len; i++)
        dst[i] += src3[i];
}
#endif

#ifdef AAC_SSR
static void ssr_context_init(ssr_context * ctx) {
    int b, i;
    for (b = 0; b < 2; b++) {
        for (i = 0; i < 4; i++) {
            // 2
            ctx->q[b][i] = cos((2*i+1)*(2*b+1-4)*M_PI/16);
            ctx->q[b+2][i] = cos((2*i+1)*(2*(b+4)+1-4)*M_PI/16);
        }
    }
    for (b = 0; b < 4; b++) {
        for (i = 0; i < 12; i++) {
            float sgn = 1 - 2 * (i&1); //4
            if (i < 6) {
                ctx->t0[b][i] = sgn * ssr_q_table[8*i + b];
                ctx->t1[b][i] = sgn * ssr_q_table[8*i + b + 4];
            } else {
                ctx->t0[b][i] = sgn * ssr_q_table[95 - (8*i + b)];
                ctx->t1[b][i] = sgn * ssr_q_table[95 - (8*i + b + 4)];
            }
        }
    }
}
#endif /* AAC_SSR */

/**
 * Free a channel element.
 */
static void che_freep(ChannelElement **s) {
#ifdef AAC_SSR
    av_freep(&(*s)->ch[0].ssr);
    av_freep(&(*s)->ch[1].ssr);
#endif /* AAC_SSR */
#ifdef AAC_LTP
    av_freep(&(*s)->ch[0].ltp_state);
    av_freep(&(*s)->ch[1].ltp_state);
#endif /* AAC_LTP */
    av_freep(s);
}

/**
 * Configure output channel order and optional mixing based on the current
 * program configuration element and user requested channels.
 *
 * \param newpcs New program configuration struct - we only do something if it differs from the current one.
 */
static int output_configure(AACContext *ac, ProgramConfig *newpcs) {
    AVCodecContext *avctx = ac->avccontext;
    ProgramConfig * pcs = &ac->pcs;
    int i, j, channels = 0;
    float a, b;
    ChannelElement *mixdown[3] = { NULL, NULL, NULL };

    static const float mixdowncoeff[4] = {
        /* matrix mix-down coefficient, table 4.70 */
        1. / M_SQRT2,
        1. / 2.,
        1. / (2 * M_SQRT2),
        0
    };

    if(!memcmp(&ac->pcs, newpcs, sizeof(ProgramConfig)))
        return 0; /* no change */

    *pcs = *newpcs;

    /* Allocate or free elements depending on if they are in the
     * current program configuration struct.
     *
     * Set up default 1:1 output mapping.
     *
     * For a 5.1 stream the output order will be:
     *    [ Front Left ] [ Front Right ] [ Center ] [ LFE ] [ Surround Left ] [ Surround Right ]
     *
     * Locate front, center and back channels for later matrix mix-down.
     */

    for(i = 0; i < MAX_TAGID; i++) {
        for(j = 0; j < 4; j++) {
            if(pcs->che_type[j][i]) {
                if(!ac->che[j][i] && !(ac->che[j][i] = av_mallocz(sizeof(ChannelElement))))
                    return AVERROR(ENOMEM);
                if(j != ID_CCE) {
                    ac->output_data[channels++] = ac->che[j][i]->ch[0].ret;
                    ac->che[j][i]->ch[0].mixing_gain = 1.0f;
                    if(j == ID_CPE) {
                        ac->output_data[channels++] = ac->che[j][i]->ch[1].ret;
                        ac->che[j][i]->ch[1].mixing_gain = 1.0f;
                        if(!mixdown[MIXDOWN_FRONT] && pcs->che_type[j][i] == AAC_CHANNEL_FRONT)
                            mixdown[MIXDOWN_FRONT] = ac->che[j][i];
                        if(!mixdown[MIXDOWN_BACK ] && pcs->che_type[j][i] == AAC_CHANNEL_BACK)
                            mixdown[MIXDOWN_BACK ] = ac->che[j][i];
                    }
                    if(j == ID_SCE && !mixdown[MIXDOWN_CENTER] && pcs->che_type[j][i] == AAC_CHANNEL_FRONT)
                        mixdown[MIXDOWN_CENTER] = ac->che[j][i];
                }
            } else
                che_freep(&ac->che[j][i]);
        }
    }

    // allocate appropriately aligned buffer for interleaved output
    if(channels > avctx->channels)
        av_freep(&ac->interleaved_output);
    if(!ac->interleaved_output && !(ac->interleaved_output = av_malloc(channels * 1024 * sizeof(float))))
        return AVERROR(ENOMEM);

    ac->mm[MIXDOWN_FRONT] = ac->mm[MIXDOWN_BACK] = ac->mm[MIXDOWN_CENTER] = NULL;

    /* Check for matrix mix-down to mono or stereo. */

    if(avctx->request_channels && avctx->request_channels <= 2 &&
       avctx->request_channels != channels) {

        if((avctx->request_channels == 1 && pcs->mono_mixdown_tag   != -1) ||
           (avctx->request_channels == 2 && pcs->stereo_mixdown_tag != -1)) {
            /* Add support for this as soon as we get a sample so we can figure out
               exactly how this is supposed to work. */
            av_log(avctx, AV_LOG_ERROR,
                   "Mix-down using pre-mixed elements is not supported, please file a bug. "
                   "Reverting to matrix mix-down.\n");
        }

        /* We need 'center + L + R + sL + sR' for matrix mix-down. */
        if(mixdown[MIXDOWN_CENTER] && mixdown[MIXDOWN_FRONT] && mixdown[MIXDOWN_BACK]) {
            a = mixdowncoeff[pcs->mixdown_coeff_index];

            if(avctx->request_channels == 2) {
                b = 1. / (1. + (1. / M_SQRT2) + a * (pcs->pseudo_surround ? 2. : 1.));
                mixdown[MIXDOWN_CENTER]->ch[0].mixing_gain      = b / M_SQRT2;
            } else {
                b = 1. / (3. + 2. * a);
                mixdown[MIXDOWN_CENTER]->ch[0].mixing_gain      = b;
            }
            mixdown[MIXDOWN_FRONT]->ch[0].mixing_gain = b;
            mixdown[MIXDOWN_FRONT]->ch[1].mixing_gain = b;
            mixdown[MIXDOWN_BACK ]->ch[0].mixing_gain = b * a;
            mixdown[MIXDOWN_BACK ]->ch[1].mixing_gain = b * a;
            for(i = 0; i < 3; i++) ac->mm[i] = mixdown[i];

            channels = avctx->request_channels;
        } else {
            av_log(avctx, AV_LOG_WARNING, "Matrix mixing from %d to %d channels is not supported.\n",
                   channels, avctx->request_channels);
        }
    }

    avctx->channels = channels;
    return 0;
}


/**
 * Decode an array of 4 bit tag IDs, optionally interleaved with a stereo/mono switching bit.
 *
 * @param cpe_map Stereo (Channel Pair Element) map, NULL if stereo bit is not present.
 * @param sce_map mono (Single Channel Element) map
 * @param type speaker type/position for these channels
 */
static void program_config_element_parse_tags(GetBitContext * gb, enum ChannelType *cpe_map,
                                              enum ChannelType *sce_map, int n, enum ChannelType type) {
    while(n--) {
        enum ChannelType *map = cpe_map && get_bits1(gb) ? cpe_map : sce_map; // stereo or mono map
        map[get_bits(gb, 4)] = type;
    }
}


/**
 * Parse program configuration element; reference: table 4.2.
 */
static int program_config_element(AACContext * ac, GetBitContext * gb) {
    ProgramConfig pcs;
    int num_front, num_side, num_back, num_lfe, num_assoc_data, num_cc;

    memset(&pcs, 0, sizeof(pcs));

    skip_bits(gb, 2);  // object_type

    ac->m4ac.sampling_index = get_bits(gb, 4);
    if(ac->m4ac.sampling_index > 12) {
        av_log(ac->avccontext, AV_LOG_ERROR, "invalid sampling rate index %d\n", ac->m4ac.sampling_index);
        return -1;
    }
    ac->m4ac.sample_rate = ff_mpeg4audio_sample_rates[ac->m4ac.sampling_index];
    num_front       = get_bits(gb, 4);
    num_side        = get_bits(gb, 4);
    num_back        = get_bits(gb, 4);
    num_lfe         = get_bits(gb, 2);
    num_assoc_data  = get_bits(gb, 3);
    num_cc          = get_bits(gb, 4);

    pcs.mono_mixdown_tag   = get_bits1(gb) ? get_bits(gb, 4) : -1;
    pcs.stereo_mixdown_tag = get_bits1(gb) ? get_bits(gb, 4) : -1;

    if (get_bits1(gb)) {
        pcs.mixdown_coeff_index = get_bits(gb, 2);
        pcs.pseudo_surround     = get_bits1(gb);
    }

    program_config_element_parse_tags(gb, pcs.che_type[ID_CPE], pcs.che_type[ID_SCE], num_front, AAC_CHANNEL_FRONT);
    program_config_element_parse_tags(gb, pcs.che_type[ID_CPE], pcs.che_type[ID_SCE], num_side,  AAC_CHANNEL_SIDE );
    program_config_element_parse_tags(gb, pcs.che_type[ID_CPE], pcs.che_type[ID_SCE], num_back,  AAC_CHANNEL_BACK );
    program_config_element_parse_tags(gb, NULL,                 pcs.che_type[ID_LFE], num_lfe,   AAC_CHANNEL_LFE  );

    skip_bits_long(gb, 4 * num_assoc_data);

    program_config_element_parse_tags(gb, pcs.che_type[ID_CCE], pcs.che_type[ID_CCE], num_cc,    AAC_CHANNEL_CC   );

    align_get_bits(gb);

    /* comment field, first byte is length */
    skip_bits_long(gb, 8 * get_bits(gb, 8));
    return output_configure(ac, &pcs);
}

/**
 * Set up ProgramConfig, but based on a default channel configuration
 * as specified in table 1.17.
 */
static int program_config_element_default(AACContext *ac, int channels)
{
    ProgramConfig pcs;

    memset(&pcs, 0, sizeof(ProgramConfig));

    /* Pre-mixed down-mix outputs are not available. */
    pcs.mono_mixdown_tag   = -1;
    pcs.stereo_mixdown_tag = -1;

    if(channels < 1 || channels > 7) {
        av_log(ac->avccontext, AV_LOG_ERROR, "invalid default channel configuration (%d channels)\n",
               channels);
        return -1;
    }

    /* default channel configurations:
     *
     * 1ch : front center (mono)
     * 2ch : L + R (stereo)
     * 3ch : front center + L + R
     * 4ch : front center + L + R + back center
     * 5ch : front center + L + R + back stereo
     * 6ch : front center + L + R + back stereo + LFE
     * 7ch : front center + L + R + outer front left + outer front right + back stereo + LFE
     */

    if(channels != 2)
        pcs.che_type[ID_SCE][0] = AAC_CHANNEL_FRONT; // front center (or mono)
    if(channels > 1)
        pcs.che_type[ID_CPE][0] = AAC_CHANNEL_FRONT; // L + R (or stereo)
    if(channels == 4)
        pcs.che_type[ID_SCE][1] = AAC_CHANNEL_BACK;  // back center
    if(channels > 4)
        pcs.che_type[ID_CPE][(channels == 7) + 1]
                                = AAC_CHANNEL_BACK;  // back stereo
    if(channels > 5)
        pcs.che_type[ID_LFE][0] = AAC_CHANNEL_LFE;   // LFE
    if(channels == 7)
        pcs.che_type[ID_CPE][1] = AAC_CHANNEL_FRONT; // outer front left + outer front right

    return output_configure(ac, &pcs);
}


/**
 * Parse GA "General Audio" specific configuration; reference: table 4.1.
 */
static int ga_specific_config(AACContext * ac, GetBitContext * gb, int channels) {
    int extension_flag, ret;

    if(get_bits1(gb)) {  // frameLengthFlag
        av_log(ac->avccontext, AV_LOG_ERROR, "960/120 MDCT window is not supported.\n");
        return -1;
    }

    if (get_bits1(gb))       // dependsOnCoreCoder
        skip_bits(gb, 14);   // coreCoderDelay
    extension_flag = get_bits1(gb);

    if(ac->m4ac.object_type == AOT_AAC_SCALABLE ||
       ac->m4ac.object_type == AOT_ER_AAC_SCALABLE)
        skip_bits(gb, 3);     // layerNr

    if (channels == 0) {
        skip_bits(gb, 4);  // element_instance_tag
        if((ret = program_config_element(ac, gb)))
            return ret;
    } else {
        if((ret = program_config_element_default(ac, channels)))
            return ret;
    }

    if (extension_flag) {
        switch (ac->m4ac.object_type) {
            case AOT_ER_BSAC:
                skip_bits(gb, 5);    // numOfSubFrame
                skip_bits(gb, 11);   // layer_length
                break;
            case AOT_ER_AAC_LC:
            case AOT_ER_AAC_LTP:
            case AOT_ER_AAC_SCALABLE:
            case AOT_ER_AAC_LD:
                skip_bits(gb, 3);  /* aacSectionDataResilienceFlag
                                    * aacScalefactorDataResilienceFlag
                                    * aacSpectralDataResilienceFlag
                                    */
                break;
        }
        skip_bits1(gb);    // extensionFlag3 (TBD in version 3)
    }
    return 0;
}


/**
 * Parse audio specific configuration; reference: table 1.13.
 *
 * @param   data        pointer to AVCodecContext extradata
 * @param   data_size   size of AVCCodecContext extradata
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int audio_specific_config(AACContext * ac, void *data, int data_size) {
    GetBitContext gb;
    int i;

    init_get_bits(&gb, data, data_size * 8);

    if((i = ff_mpeg4audio_get_config(&ac->m4ac, data, data_size)) < 0)
        return -1;

    skip_bits_long(&gb, i);

    switch (ac->m4ac.object_type) {
    case AOT_AAC_LC:
#ifdef AAC_SSR
    case AOT_AAC_SSR:
#endif /* AAC_SSR */
#ifdef AAC_LTP
    case AOT_AAC_LTP:
#endif /* AAC_LTP */
        if (ga_specific_config(ac, &gb, ac->m4ac.chan_config))
            return -1;
        break;
    default:
        av_log(ac->avccontext, AV_LOG_ERROR, "Audio object type %s%d is not supported.\n",
               ac->m4ac.sbr == 1? "SBR+" : "", ac->m4ac.object_type);
        return -1;
    }
    return 0;
}

static av_cold int aac_decode_init(AVCodecContext * avccontext) {
    static const struct {
        const uint16_t     *a_code;
        const unsigned int s;
        const uint8_t      *a_bits;
    } tmp[] = {
        { code1 , sizeof code1 , bits1  },
        { code2 , sizeof code2 , bits2  },
        { code3 , sizeof code3 , bits3  },
        { code4 , sizeof code4 , bits4  },
        { code5 , sizeof code5 , bits5  },
        { code6 , sizeof code6 , bits6  },
        { code7 , sizeof code7 , bits7  },
        { code8 , sizeof code8 , bits8  },
        { code9 , sizeof code9 , bits9  },
        { code10, sizeof code10, bits10 },
        { code11, sizeof code11, bits11 },
    };
    AACContext * ac = avccontext->priv_data;
    int i;

    ac->avccontext = avccontext;

    if (avccontext->extradata_size <= 0 ||
        audio_specific_config(ac, avccontext->extradata, avccontext->extradata_size))
        return -1;

    avccontext->sample_rate = ac->m4ac.sample_rate;
    avccontext->frame_size  = 1024;

    INIT_VLC_STATIC(&books[0], 6, tmp[0].s/sizeof(tmp[0].a_code[0]),
        tmp[0].a_bits, sizeof(tmp[0].a_bits[0]), sizeof(tmp[0].a_bits[0]),
        tmp[0].a_code, sizeof(tmp[0].a_code[0]), sizeof(tmp[0].a_code[0]),
        144);
    INIT_VLC_STATIC(&books[1], 6, tmp[1].s/sizeof(tmp[1].a_code[0]),
        tmp[1].a_bits, sizeof(tmp[1].a_bits[0]), sizeof(tmp[1].a_bits[0]),
        tmp[1].a_code, sizeof(tmp[1].a_code[0]), sizeof(tmp[1].a_code[0]),
        114);
    INIT_VLC_STATIC(&books[2], 6, tmp[2].s/sizeof(tmp[2].a_code[0]),
        tmp[2].a_bits, sizeof(tmp[2].a_bits[0]), sizeof(tmp[2].a_bits[0]),
        tmp[2].a_code, sizeof(tmp[2].a_code[0]), sizeof(tmp[2].a_code[0]),
        188);
    INIT_VLC_STATIC(&books[3], 6, tmp[3].s/sizeof(tmp[3].a_code[0]),
        tmp[3].a_bits, sizeof(tmp[3].a_bits[0]), sizeof(tmp[3].a_bits[0]),
        tmp[3].a_code, sizeof(tmp[3].a_code[0]), sizeof(tmp[3].a_code[0]),
        180);
    INIT_VLC_STATIC(&books[4], 6, tmp[4].s/sizeof(tmp[4].a_code[0]),
        tmp[4].a_bits, sizeof(tmp[4].a_bits[0]), sizeof(tmp[4].a_bits[0]),
        tmp[4].a_code, sizeof(tmp[4].a_code[0]), sizeof(tmp[4].a_code[0]),
        172);
    INIT_VLC_STATIC(&books[5], 6, tmp[5].s/sizeof(tmp[5].a_code[0]),
        tmp[5].a_bits, sizeof(tmp[5].a_bits[0]), sizeof(tmp[5].a_bits[0]),
        tmp[5].a_code, sizeof(tmp[5].a_code[0]), sizeof(tmp[5].a_code[0]),
        140);
    INIT_VLC_STATIC(&books[6], 6, tmp[6].s/sizeof(tmp[6].a_code[0]),
        tmp[6].a_bits, sizeof(tmp[6].a_bits[0]), sizeof(tmp[6].a_bits[0]),
        tmp[6].a_code, sizeof(tmp[6].a_code[0]), sizeof(tmp[6].a_code[0]),
        168);
    INIT_VLC_STATIC(&books[7], 6, tmp[7].s/sizeof(tmp[7].a_code[0]),
        tmp[7].a_bits, sizeof(tmp[7].a_bits[0]), sizeof(tmp[7].a_bits[0]),
        tmp[7].a_code, sizeof(tmp[7].a_code[0]), sizeof(tmp[7].a_code[0]),
        114);
    INIT_VLC_STATIC(&books[8], 6, tmp[8].s/sizeof(tmp[8].a_code[0]),
        tmp[8].a_bits, sizeof(tmp[8].a_bits[0]), sizeof(tmp[8].a_bits[0]),
        tmp[8].a_code, sizeof(tmp[8].a_code[0]), sizeof(tmp[8].a_code[0]),
        262);
    INIT_VLC_STATIC(&books[9], 6, tmp[9].s/sizeof(tmp[9].a_code[0]),
        tmp[9].a_bits, sizeof(tmp[9].a_bits[0]), sizeof(tmp[9].a_bits[0]),
        tmp[9].a_code, sizeof(tmp[9].a_code[0]), sizeof(tmp[9].a_code[0]),
        248);
    INIT_VLC_STATIC(&books[10], 6, tmp[10].s/sizeof(tmp[10].a_code[0]),
        tmp[10].a_bits, sizeof(tmp[10].a_bits[0]), sizeof(tmp[10].a_bits[0]),
        tmp[10].a_code, sizeof(tmp[10].a_code[0]), sizeof(tmp[10].a_code[0]),
        384);

    dsputil_init(&ac->dsp, avccontext);

    av_init_random(0x1f2e3d4c, &ac->random_state);

    // -1024 - Compensate wrong IMDCT method.
    // 32768 - Required to scale values to the correct range for the bias method
    //         for float to int16 conversion.

    if(ac->dsp.float_to_int16 == ff_float_to_int16_c) {
        ac->add_bias = 385.0f;
        ac->sf_scale = 1. / (-1024. * 32768.);
        ac->sf_offset = 0;
    } else {
        ac->add_bias = 0.0f;
        ac->sf_scale = 1. / -1024.;
        ac->sf_offset = 60;
    }

#ifndef CONFIG_HARDCODED_TABLES
    for (i = 1 - IVQUANT_SIZE/2; i < IVQUANT_SIZE/2; i++)
        ivquant_tab[i + IVQUANT_SIZE/2 - 1] =  cbrt(fabs(i)) * i;
    for (i = 0; i < 316; i++)
        pow2sf_tab[i] = pow(2, (i - 200)/4.);
#endif /* CONFIG_HARDCODED_TABLES */

    INIT_VLC_STATIC(&mainvlc, 7, sizeof(code)/sizeof(code[0]),
        bits, sizeof(bits[0]), sizeof(bits[0]),
        code, sizeof(code[0]), sizeof(code[0]),
        352);

#ifdef AAC_LTP
    ff_mdct_init(&ac->mdct_ltp, 11, 0);
#endif /* AAC_LTP */
#ifdef AAC_SSR
    if (ac->m4ac.object_type == AOT_AAC_SSR) {
        ff_mdct_init(&ac->mdct, 9, 1);
        ff_mdct_init(&ac->mdct_small, 6, 1);
        // window initialization
        ff_kbd_window_init(kbd_long_1024, 4.0, 256);
        ff_kbd_window_init(kbd_short_128, 6.0, 32);
        ff_sine_window_init(sine_long_1024, 256);
        ff_sine_window_init(sine_short_128, 32);
        ssr_context_init(&ac->ssrctx);
    } else {
#endif /* AAC_SSR */
        ff_mdct_init(&ac->mdct, 11, 1);
        ff_mdct_init(&ac->mdct_small, 8, 1);
        // window initialization
        ff_kbd_window_init(kbd_long_1024, 4.0, 1024);
        ff_kbd_window_init(kbd_short_128, 6.0, 128);
        ff_sine_window_init(sine_long_1024, 1024);
        ff_sine_window_init(sine_short_128, 128);
#ifdef AAC_SSR
    }
#endif /* AAC_SSR */
    return 0;
}

/**
 * Decode a data_stream_element; reference: table 4.10.
 */
static void data_stream_element(AACContext * ac, GetBitContext * gb) {
    int byte_align = get_bits1(gb);
    int count = get_bits(gb, 8);
    if (count == 255)
        count += get_bits(gb, 8);
    if (byte_align)
        align_get_bits(gb);
    skip_bits_long(gb, 8 * count);
}

#ifdef AAC_LTP
static void decode_ltp(AACContext * ac, GetBitContext * gb, uint8_t max_sfb, LongTermPrediction * ltp) {
    int sfb;
    if (ac->m4ac.object_type == AOT_ER_AAC_LD) {
        assert(0);
    } else {
        ltp->lag = get_bits(gb, 11);
        ltp->coef = ltp_coef[get_bits(gb, 3)] * ac->sf_scale;
        for (sfb = 0; sfb < FFMIN(max_sfb, MAX_LTP_LONG_SFB); sfb++)
            ltp->used[sfb] = get_bits1(gb);
    }
}
#endif /* AAC_LTP */

/**
 * Decode Individual Channel Stream info; reference: table 4.6.
 *
 * @param   common_window   Channels have independent [0], or shared [1], Individual Channel Stream information.
 */
static int decode_ics_info(AACContext * ac, GetBitContext * gb, int common_window, IndividualChannelStream * ics) {
    uint8_t grouping;
    if (get_bits1(gb)) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Reserved bit set.\n");
        return -1;
    }
    ics->window_sequence[1] = ics->window_sequence[0];
    ics->window_sequence[0] = get_bits(gb, 2);
    ics->use_kb_window[1] = ics->use_kb_window[0];
    ics->use_kb_window[0] = get_bits1(gb);
    ics->num_window_groups = 1;
    ics->group_len[0] = 1;
    if (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE) {
        int i;
        ics->max_sfb = get_bits(gb, 4);
        grouping = get_bits(gb, 7);
        for (i = 0; i < 7; i++) {
            if (grouping & (1<<(6-i))) {
                ics->group_len[ics->num_window_groups-1]++;
            } else {
                ics->num_window_groups++;
                ics->group_len[ics->num_window_groups-1] = 1;
            }
        }
        ics->swb_offset    =    swb_offset_128[ac->m4ac.sampling_index];
        ics->num_swb       =       num_swb_128[ac->m4ac.sampling_index];
    } else {
        ics->max_sfb = get_bits(gb, 6);
        ics->swb_offset    =    swb_offset_1024[ac->m4ac.sampling_index];
        ics->num_swb       =       num_swb_1024[ac->m4ac.sampling_index];
    }

    if(ics->max_sfb > ics->num_swb) {
        av_log(ac->avccontext, AV_LOG_ERROR,
            "Number of scalefactor bands in group (%d) exceeds limit (%d).\n",
            ics->max_sfb, ics->num_swb);
        return -1;
    }

    if (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE) {
        ics->num_windows   = 8;
        ics->tns_max_bands = tns_max_bands_128[ac->m4ac.sampling_index];
    } else {
        ics->num_windows   = 1;
        ics->tns_max_bands = tns_max_bands_1024[ac->m4ac.sampling_index];
        if (get_bits1(gb)) {
#ifdef AAC_LTP
            if (ac->m4ac.object_type == AOT_AAC_MAIN) {
                assert(0);
            } else {
                if ((ics->ltp.present = get_bits(gb, 1))) {
                    decode_ltp(ac, gb, ics->max_sfb, &ics->ltp);
                }
                if (common_window) {
                    if ((ics->ltp2.present = get_bits(gb, 1))) {
                        decode_ltp(ac, gb, ics->max_sfb, &ics->ltp2);
                    }
                }
            }
#else /* AAC_LTP */
            av_log(ac->avccontext, AV_LOG_ERROR,
                   "Predictor bit set but LTP is not supported.\n");
            return -1;
#endif /* AAC_LTP */
#ifdef AAC_LTP
        } else {
            ics->ltp.present = 0;
            ics->ltp2.present = 0;
#endif /* AAC_LTP */
        }
    }
    return 0;
}

/**
 * inverse quantization
 *
 * @param   a   quantized value to be dequantized
 * @return  Returns dequantized value.
 */
static inline float ivquant(AACContext * ac, int a) {
    if (a + (unsigned int)IVQUANT_SIZE/2 - 1 < (unsigned int)IVQUANT_SIZE - 1)
        return ivquant_tab[a + IVQUANT_SIZE/2 - 1];
    else
        return cbrtf(fabsf(a)) * a;
}

/**
 * Decode section_data payload; reference: table 4.46.
 *
 * @param   band_type           array of the band type used for a window group's scalefactor band
 * @param   band_type_run_end   array of the last scalefactor band of a band type run for a window group's scalefactor band
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_band_types(AACContext * ac, GetBitContext * gb, IndividualChannelStream * ics, enum BandType band_type[][64], int band_type_run_end[][64]) {
    int g;
    const int bits = (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE) ? 3 : 5;
    for (g = 0; g < ics->num_window_groups; g++) {
        int k = 0;
        while (k < ics->max_sfb) {
            uint8_t sect_len = k;
            int sect_len_incr;
            int sect_band_type = get_bits(gb, 4);
            if (sect_band_type == 12) {
                av_log(ac->avccontext, AV_LOG_ERROR, "invalid band type\n");
                return -1;
            }
            while ((sect_len_incr = get_bits(gb, bits)) == (1 << bits)-1)
                sect_len += sect_len_incr;
            sect_len += sect_len_incr;
            if (sect_len > ics->max_sfb) {
                av_log(ac->avccontext, AV_LOG_ERROR,
                    "Number of bands (%d) exceeds limit (%d).\n",
                    sect_len, ics->max_sfb);
                return -1;
            }
            for (; k < sect_len; k++) {
                band_type        [g][k] = sect_band_type;
                band_type_run_end[g][k] = sect_len;
            }
        }
    }
    return 0;
}

/**
 * Decode scale_factor_data; reference: table 4.47.
 *
 * @param   mix_gain            channel gain (Not used by AAC bitstream.)
 * @param   global_gain         first scalefactor value as scalefactors are differentially coded
 * @param   band_type           array of the band type used for a window group's scalefactor band
 * @param   band_type_run_end   array of the last scalefactor band of a band type run for a window group's scalefactor band
 * @param   sf                  array of scalefactors or intensity stereo positions used for a window group's scalefactor band
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_scalefactors(AACContext * ac, GetBitContext * gb, float mix_gain, unsigned int global_gain,
        IndividualChannelStream * ics, const enum BandType band_type[][64], const int band_type_run_end[][64], float sf[][64]) {
    const int sf_offset = ac->sf_offset + (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE ? 12 : 0);
    int g, i;
    int offset[3] = { global_gain, global_gain - 90, 100 };
    int noise_flag = 1;
    static const char *sf_str[3] = { "Global gain", "Noise gain", "Intensity stereo position" };
    ics->intensity_present = 0;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb;) {
            int run_end = band_type_run_end[g][i];
            if (band_type[g][i] == ZERO_BT) {
                for(; i < run_end; i++)
                    sf[g][i] = 0.;
            }else if((band_type[g][i] == INTENSITY_BT) || (band_type[g][i] == INTENSITY_BT2)) {
                ics->intensity_present = 1;
                for(; i < run_end; i++) {
                    offset[2] += get_vlc2(gb, mainvlc.table, 7, 3) - 60;
                    if(offset[2] > 255) {
                        av_log(ac->avccontext, AV_LOG_ERROR,
                            "%s (%d) out of range.\n", sf_str[2], offset[2]);
                        return -1;
                    }
                    sf[g][i] =  pow2sf_tab[-offset[2] + 300];
                    sf[g][i] *= mix_gain;
                }
            }else if(band_type[g][i] == NOISE_BT) {
                for(; i < run_end; i++) {
                    if(noise_flag-- > 0)
                        offset[1] += get_bits(gb, 9) - 256;
                    else
                        offset[1] += get_vlc2(gb, mainvlc.table, 7, 3) - 60;
                    if(offset[1] > 255) {
                        av_log(ac->avccontext, AV_LOG_ERROR,
                            "%s (%d) out of range.\n", sf_str[1], offset[1]);
                        return -1;
                    }
                    sf[g][i] = -pow2sf_tab[ offset[1] + sf_offset];
                    sf[g][i] *= mix_gain;
                }
            }else {
                for(; i < run_end; i++) {
                    offset[0] += get_vlc2(gb, mainvlc.table, 7, 3) - 60;
                    if(offset[0] > 255) {
                        av_log(ac->avccontext, AV_LOG_ERROR,
                            "%s (%d) out of range.\n", sf_str[0], offset[0]);
                        return -1;
                    }
                    sf[g][i] = -pow2sf_tab[ offset[0] + sf_offset];
                    sf[g][i] *= mix_gain;
                }
            }
        }
    }
    return 0;
}

/**
 * Decode pulse data; reference: table 4.7.
 */
static void decode_pulses(AACContext * ac, GetBitContext * gb, Pulse * pulse) {
    int i;
    pulse->num_pulse = get_bits(gb, 2) + 1;
    pulse->start = get_bits(gb, 6);
    for (i = 0; i < pulse->num_pulse; i++) {
        pulse->offset[i] = get_bits(gb, 5);
        pulse->amp   [i] = get_bits(gb, 4);
    }
}

/**
 * Decode Temporal Noise Shaping data; reference: table 4.48.
 */
static void decode_tns(AACContext * ac, GetBitContext * gb, const IndividualChannelStream * ics, TemporalNoiseShaping * tns) {
    int w, filt, i, coef_len, coef_res = 0, coef_compress;
    const int is8 = ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE;
    for (w = 0; w < ics->num_windows; w++) {
        tns->n_filt[w] = get_bits(gb, 2 - is8);

        if (tns->n_filt[w])
            coef_res = get_bits1(gb) + 3;

        for (filt = 0; filt < tns->n_filt[w]; filt++) {
            tns->length[w][filt] = get_bits(gb, 6 - 2*is8);

            if ((tns->order[w][filt] = get_bits(gb, 5 - 2*is8))) {
                tns->direction[w][filt] = get_bits1(gb);
                coef_compress = get_bits1(gb);
                coef_len = coef_res - coef_compress;
                tns->tmp2_map[w][filt] = tns_tmp2_map[2*coef_compress + coef_res - 3];

                for (i = 0; i < tns->order[w][filt]; i++)
                    tns->coef[w][filt][i] = get_bits(gb, coef_len);
            }
        }
    }
}

#ifdef AAC_SSR
static int decode_gain_control(AACContext * ac, GetBitContext * gb, SingleChannelElement * sce) {
    // wd_num wd_test aloc_size
    static const int gain_mode[4][3] = {
        {1, 0, 5}, //ONLY_LONG_SEQUENCE = 0,
        {2, 1, 2}, //LONG_START_SEQUENCE,
        {8, 0, 2}, //EIGHT_SHORT_SEQUENCE,
        {2, 1, 5}, //LONG_STOP_SEQUENCE
    };
    const int mode = sce->ics.window_sequence[0];
    int bd, wd, ad;
    ScalableSamplingRate * ssr = sce->ssr;
    if (!ssr && !(ssr = sce->ssr = av_mallocz(sizeof(ScalableSamplingRate))))
        return AVERROR(ENOMEM);
    ssr->max_band = get_bits(gb, 2);
    for (bd = 0; bd < ssr->max_band; bd++) {
        for (wd = 0; wd < gain_mode[mode][0]; wd++) {
            ssr->adjust_num[bd][wd] = get_bits(gb, 3);
            for (ad = 0; ad < ssr->adjust_num[bd][wd]; ad++) {
                ssr->alev[bd][wd][ad] = get_bits(gb, 4);
                if (gain_mode[mode][1] && (wd == 0))
                    ssr->aloc[bd][wd][ad] = get_bits(gb, 4);
                else
                    ssr->aloc[bd][wd][ad] = get_bits(gb, gain_mode[mode][2]);
            }
        }
    }
    return 0;
}
#endif /* AAC_SSR */

/**
 * Decode Mid/Side data; reference: table 4.54.
 */
static void decode_mid_side_stereo(AACContext * ac, GetBitContext * gb, ChannelElement * cpe) {
    MidSideStereo * ms = &cpe->ms;
    int g, i;
    ms->present = get_bits(gb, 2);
    if (ms->present == 1) {
        for (g = 0; g < cpe->ch[0].ics.num_window_groups; g++)
            for (i = 0; i < cpe->ch[0].ics.max_sfb; i++)
                ms->mask[g][i] = get_bits1(gb);// << i;
    } else if (ms->present == 2) {
        for (g = 0; g < cpe->ch[0].ics.num_window_groups; g++)
            memset(ms->mask[g], 1, cpe->ch[0].ics.max_sfb * sizeof(ms->mask[g][0]));
    }
}

/**
 * Decode spectral data; reference: table 4.50.
 *
 * @param   band_type   array of the band type used for a window group's scalefactor band
 * @param   icoef       array of quantized spectral data
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_spectrum(AACContext * ac, GetBitContext * gb, const IndividualChannelStream * ics, const enum BandType band_type[][64], int icoef[1024]) {
    int i, k, g;
    const uint16_t * offsets = ics->swb_offset;

    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            const int cur_band_type = band_type[g][i];
            const int dim = cur_band_type >= FIRST_PAIR_BT ? 2 : 4;
            const int is_cb_unsigned = IS_CODEBOOK_UNSIGNED(cur_band_type);
            int group;
            if (cur_band_type == ZERO_BT) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    memset(icoef + group * 128 + offsets[i], 0, (offsets[i+1] - offsets[i])*sizeof(int));
                }
            }else if (cur_band_type != NOISE_BT && cur_band_type != INTENSITY_BT2 && cur_band_type != INTENSITY_BT) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k += dim) {
                        const int index = get_vlc2(gb, books[cur_band_type - 1].table, 6, 3);
                        const int coef_idx = (group << 7) + k;
                        const int8_t *vq_ptr = &codebook_vectors[cur_band_type - 1][index * dim];
                        int j;
                        if (is_cb_unsigned) {
                            for (j = 0; j < dim; j++)
                                if (vq_ptr[j])
                                    icoef[coef_idx + j] = 1 - 2*get_bits1(gb);
                        }else {
                            for (j = 0; j < dim; j++)
                                icoef[coef_idx + j] = 1;
                        }
                        if (cur_band_type == ESC_BT) {
                            for (j = 0; j < 2; j++) {
                                if (vq_ptr[j] == 16) {
                                    int n = 4;
                                    /* The total length of escape_sequence must be < 22 bits according
                                       to the specification (i.e. max is 11111111110xxxxxxxxxx). */
                                    while (get_bits1(gb) && n < 15) n++;
                                    if(n == 15) {
                                        av_log(ac->avccontext, AV_LOG_ERROR, "error in spectral data, ESC overflow\n");
                                        return -1;
                                    }
                                    icoef[coef_idx + j] *= (1<<n) + get_bits(gb, n);
                                }else
                                    icoef[coef_idx + j] *= vq_ptr[j];
                            }
                        }else
                            for (j = 0; j < dim; j++)
                                icoef[coef_idx + j] *= vq_ptr[j];
                    }
                }
            }
        }
        icoef += ics->group_len[g]<<7;
    }
    return 0;
}

/**
 * Add pulses with particular amplitudes to the quantized spectral data; reference: 4.6.3.3.
 *
 * @param   pulse   pointer to pulse data struct
 * @param   icoef   array of quantized spectral data
 */
static void add_pulses(AACContext * ac, const IndividualChannelStream * ics, const Pulse * pulse, int icoef[1024]) {
    int i, off = ics->swb_offset[pulse->start];
    for (i = 0; i < pulse->num_pulse; i++) {
        int ic;
        off += pulse->offset[i];
        ic = (icoef[off] - 1)>>31;
        icoef[off] += (pulse->amp[i]^ic) - ic;
    }
}

/**
 * Dequantize and scale spectral data; reference: 4.6.3.3.
 *
 * @param   icoef       array of quantized spectral data
 * @param   band_type   array of the band type used for a window group's scalefactor band
 * @param   sf          array of scalefactors or intensity stereo positions used for a window group's scalefactor band
 * @param   coef        array of dequantized, scaled spectral data
 */
static void dequant(AACContext * ac, const IndividualChannelStream * ics, const int icoef[1024],
        const enum BandType band_type[][64], const float sf[][64], float coef[1024]) {
    const uint16_t * offsets = ics->swb_offset;
    const int c = 1024/ics->num_window_groups;
    int g, i, group, k;

    for (g = 0; g < ics->num_window_groups; g++) {
        memset(coef + g * 128 + offsets[ics->max_sfb], 0, sizeof(float)*(c - offsets[ics->max_sfb]));
        for (i = 0; i < ics->max_sfb; i++) {
            if (band_type[g][i] == NOISE_BT) {
                const float scale = sf[g][i] / ((offsets[i+1] - offsets[i]) * PNS_MEAN_ENERGY);
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++)
                        coef[group*128+k] = (int32_t)av_random(&ac->random_state) * scale;
                }
            } else if (band_type[g][i] != INTENSITY_BT && band_type[g][i] != INTENSITY_BT2) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++) {
                        coef[group*128+k] = ivquant(ac, icoef[group*128+k]) * sf[g][i];
                    }
                }
            }
        }
        coef  += ics->group_len[g]*128;
        icoef += ics->group_len[g]*128;
    }
}

/**
 * Decode an individual_channel_stream payload; reference: table 4.44.
 *
 * @param   common_window   Channels have independent [0], or shared [1], Individual Channel Stream information.
 * @param   scale_flag      scalable [1] or non-scalable [0] AAC (Unused until scalable AAC is implemented.)
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_ics(AACContext * ac, GetBitContext * gb, int common_window, int scale_flag, SingleChannelElement * sce) {
    int icoeffs[1024];
    Pulse pulse;
    TemporalNoiseShaping * tns = &sce->tns;
    IndividualChannelStream * ics = &sce->ics;
    float * out = sce->coeffs;
    int global_gain, pulse_present = 0;

    global_gain = get_bits(gb, 8);

    if (!common_window && !scale_flag) {
        if (decode_ics_info(ac, gb, 0, ics) < 0)
            return -1;
    }

    if (decode_band_types(ac, gb, ics, sce->band_type, sce->band_type_run_end) < 0)
        return -1;
    if (decode_scalefactors(ac, gb, sce->mixing_gain, global_gain, ics, sce->band_type, sce->band_type_run_end, sce->sf) < 0)
        return -1;

    pulse_present = 0;
    if (!scale_flag) {
        if ((pulse_present = get_bits1(gb))) {
            if (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE) {
                av_log(ac->avccontext, AV_LOG_ERROR, "Pulse tool not allowed in eight short sequence.\n");
                return -1;
            }
            decode_pulses(ac, gb, &pulse);
        }
        if ((tns->present = get_bits1(gb)))
            decode_tns(ac, gb, ics, tns);
        if (get_bits1(gb)) {
#ifdef AAC_SSR
            int ret;
            if ((ret = decode_gain_control(ac, gb, sce))) return ret;
#else
            av_log(ac->avccontext, AV_LOG_ERROR, "SSR not supported.\n");
            return -1;
#endif
        }
    }

    if (decode_spectrum(ac, gb, ics, sce->band_type, icoeffs) < 0)
        return -1;
    if (pulse_present)
        add_pulses(ac, ics, &pulse, icoeffs);
    dequant(ac, ics, icoeffs, sce->band_type, sce->sf, out);
    return 0;
}

/**
 * Mid/Side stereo decoding; reference: 4.6.8.1.3.
 */
static void apply_mid_side_stereo(AACContext * ac, ChannelElement * cpe) {
    const MidSideStereo * ms = &cpe->ms;
    const IndividualChannelStream * ics = &cpe->ch[0].ics;
    float *ch0 = cpe->ch[0].coeffs;
    float *ch1 = cpe->ch[1].coeffs;
    if (ms->present) {
        int g, i, k, gp;
        const uint16_t * offsets = ics->swb_offset;
        for (g = 0; g < ics->num_window_groups; g++) {
            for (gp = 0; gp < ics->group_len[g]; gp++) {
                for (i = 0; i < ics->max_sfb; i++) {
                    if (ms->mask[g][i] &&
                        cpe->ch[0].band_type[g][i] < NOISE_BT && cpe->ch[1].band_type[g][i] < NOISE_BT) {
                        for (k = offsets[i]; k < offsets[i+1]; k++) {
                            float tmp = ch0[k] - ch1[k];
                            ch0[k] += ch1[k];
                            ch1[k] = tmp;
                        }
                    }
                }
                ch0 += 128;
                ch1 += 128;
            }
        }
    }
}

/**
 * intensity stereo decoding; reference: 4.6.8.2.3
 */
static void apply_intensity_stereo(AACContext * ac, ChannelElement * cpe) {
    const IndividualChannelStream * ics = &cpe->ch[1].ics;
    SingleChannelElement * sce1 = &cpe->ch[1];
    float *coef0 = cpe->ch[0].coeffs, *coef1 = cpe->ch[1].coeffs;
    const uint16_t * offsets = ics->swb_offset;
    int g, gp, i, k;
    int c;
    float scale;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (gp = 0; gp < ics->group_len[g]; gp++) {
            for (i = 0; i < ics->max_sfb;) {
                if (sce1->band_type[g][i] == INTENSITY_BT || sce1->band_type[g][i] == INTENSITY_BT2) {
                    const int bt_run_end = sce1->band_type_run_end[g][i];
                    while (i < bt_run_end) {
                        c = -1 + 2 * (sce1->band_type[g][i] - 14);
                        if (cpe->ms.present)
                            c *= 1 - 2 * cpe->ms.mask[g][i];
                        scale = c * sce1->sf[g][i];
                        for (k = offsets[i]; k < offsets[i+1]; k++)
                            coef1[k] = scale * coef0[k];
                        i++;
                    }
                } else
                    i = sce1->band_type_run_end[g][i];
            }
            coef0 += 128;
            coef1 += 128;
        }
    }
}

/**
 * Decode a channel_pair_element; reference: table 4.4.
 *
 * @param   tag Identifies the instance of a syntax element.
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_cpe(AACContext * ac, GetBitContext * gb, int tag) {
    int i, ret, common_window;
    ChannelElement * cpe;

    cpe = ac->che[ID_CPE][tag];
    common_window = get_bits1(gb);
    if (common_window) {
        if (decode_ics_info(ac, gb, 1, &cpe->ch[0].ics))
            return -1;
        i = cpe->ch[1].ics.use_kb_window[0];
        cpe->ch[1].ics = cpe->ch[0].ics;
        cpe->ch[1].ics.use_kb_window[1] = i;
#ifdef AAC_LTP
        cpe->ch[1].ics.ltp = cpe->ch[0].ics.ltp2;
#endif /* AAC_LTP */
        decode_mid_side_stereo(ac, gb, cpe);
    } else {
        cpe->ms.present = 0;
    }
    if ((ret = decode_ics(ac, gb, common_window, 0, &cpe->ch[0])))
        return ret;
    if ((ret = decode_ics(ac, gb, common_window, 0, &cpe->ch[1])))
        return ret;

    if (common_window)
        apply_mid_side_stereo(ac, cpe);

    if (cpe->ch[1].ics.intensity_present)
        apply_intensity_stereo(ac, cpe);
    return 0;
}

/**
 * Decode coupling_channel_element; reference: table 4.8.
 *
 * @param   tag Identifies the instance of a syntax element.
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_cce(AACContext * ac, GetBitContext * gb, int tag) {
    int num_gain = 0;
    int c, g, sfb, ret;
    int sign;
    float scale;
    SingleChannelElement * sce;
    ChannelCoupling * coup;

    sce = &ac->che[ID_CCE][tag]->ch[0];
    sce->mixing_gain = 1.0;

    coup = &ac->che[ID_CCE][tag]->coup;

    coup->is_indep_coup = get_bits1(gb);
    coup->num_coupled = get_bits(gb, 3);
    for (c = 0; c <= coup->num_coupled; c++) {
        num_gain++;
        coup->is_cpe[c] = get_bits1(gb);
        coup->tag_select[c] = get_bits(gb, 4);
        if (coup->is_cpe[c]) {
            coup->l[c] = get_bits1(gb);
            coup->r[c] = get_bits1(gb);
            if (coup->l[c] && coup->r[c])
                num_gain++;
        }
    }
    coup->domain = get_bits1(gb);
    sign = get_bits(gb, 1);
    scale = pow(2., pow(2., get_bits(gb, 2) - 3));

    if ((ret = decode_ics(ac, gb, 0, 0, sce)))
        return ret;

    for (c = 0; c < num_gain; c++) {
        int cge = 1;
        int gain = 0;
        float gain_cache = 1.;
        if (c) {
            cge = coup->is_indep_coup ? 1 : get_bits1(gb);
            gain = cge ? get_vlc2(gb, mainvlc.table, 7, 3) - 60: 0;
            gain_cache = pow(scale, gain);
        }
        for (g = 0; g < sce->ics.num_window_groups; g++)
            for (sfb = 0; sfb < sce->ics.max_sfb; sfb++)
                if (sce->band_type[g][sfb] == ZERO_BT) {
                    coup->gain[c][g][sfb] = 0;
                } else {
                    if (cge) {
                        coup->gain[c][g][sfb] = gain_cache;
                    } else {
                        int s, t = get_vlc2(gb, mainvlc.table, 7, 3) - 60;
                        if (sign) {
                            s = 1 - 2 * (t & 0x1);
                            t >>= 1;
                        } else
                            s = 1;
                        gain += t;
                        coup->gain[c][g][sfb] = pow(scale, gain) * s;
                    }
                }
    }
    return 0;
}

/**
 * Parse Spectral Band Replication extension data; reference: table 4.55.
 *
 * @param   crc flag indicating the presence of CRC checksum
 * @param   cnt length of ID_FIL syntactic element in bytes
 * @return  Returns number of bytes consumed from the ID_FIL element.
 */
static int sbr_extension_data(AACContext * ac, GetBitContext * gb, int crc, int cnt) {
    // TODO : sbr_extension implementation
    av_log(ac->avccontext, AV_LOG_DEBUG, "aac: SBR not yet supported.\n");
    skip_bits_long(gb, 8*cnt - 4); // -4 due to reading extension type
    return cnt;
}

/**
 * Parse whether channels are to be excluded from Dynamic Range Compression; reference: table 4.53.
 *
 * @return  Returns number of bytes consumed.
 */
static int excluded_channels(AACContext * ac, GetBitContext * gb) {
    int i;
    int n = 1;
    int num_excl_chan = 7;

    for (i = 0; i < 7; i++)
         ac->che_drc.exclude_mask[i] = get_bits1(gb);

    while (n <= MAX_CHANNELS && num_excl_chan < MAX_CHANNELS - 7 && get_bits1(gb)) {
        ac->che_drc.additional_excluded_chns[n-1]=1;
        for (i = num_excl_chan; i < num_excl_chan+7; i++)
            ac->che_drc.exclude_mask[i] = get_bits1(gb);
        n++;
        num_excl_chan += 7;
    }
    return n;
}

/**
 * Decode dynamic range information; reference: table 4.52.
 *
 * @param   cnt length of ID_FIL syntactic element in bytes
 * @return  Returns number of bytes consumed.
 */
static int dynamic_range_info(AACContext * ac, GetBitContext * gb, int cnt) {
    int n = 1;
    int drc_num_bands = 1;
    int i;

    /* pce_tag_present? */
    if(get_bits1(gb)) {
        ac->che_drc.pce_instance_tag  = get_bits(gb, 4);
        skip_bits(gb, 4); // tag_reserved_bits
        n++;
    }

    /* excluded_chns_present? */
    if(get_bits1(gb)) {
        n += excluded_channels(ac, gb);
    }

    /* drc_bands_present? */
    if (get_bits1(gb)) {
        ac->che_drc.band_incr            = get_bits(gb, 4);
        ac->che_drc.interpolation_scheme = get_bits(gb, 4);
        n++;
        drc_num_bands += ac->che_drc.band_incr;
        for (i = 0; i < drc_num_bands; i++) {
            ac->che_drc.band_top[i] = get_bits(gb, 8);
            n++;
        }
    }

    /* prog_ref_level_present? */
    if (get_bits1(gb)) {
        ac->che_drc.prog_ref_level = get_bits(gb, 7);
        skip_bits1(gb); // prog_ref_level_reserved_bits
        n++;
    }

    for (i = 0; i < drc_num_bands; i++) {
        ac->che_drc.dyn_rng_sgn[i] = get_bits1(gb);
        ac->che_drc.dyn_rng_ctl[i] = get_bits(gb, 7);
        n++;
    }

    return n;
}

/**
 * Parse extension data (incomplete); reference: table 4.51.
 *
 * @param   cnt length of ID_FIL syntactic element in bytes
 */
static int extension_payload(AACContext * ac, GetBitContext * gb, int cnt) {
    int crc_flag = 0;
    int res = cnt;
    switch (get_bits(gb, 4)) { // extension type
        case EXT_SBR_DATA_CRC:
            crc_flag++;
        case EXT_SBR_DATA:
            res = sbr_extension_data(ac, gb, crc_flag, cnt);
            break;
        case EXT_DYNAMIC_RANGE:
            res = dynamic_range_info(ac, gb, cnt);
            break;
        case EXT_FILL:
        case EXT_FILL_DATA:
        case EXT_DATA_ELEMENT:
        default:
            skip_bits_long(gb, 8*cnt - 4);
            break;
    };
    return res;
}

/**
 * Decode Temporal Noise Shaping filter coefficients and apply all-pole filters; reference: 4.6.9.3.
 *
 * @param   decode  1 if tool is used normally, 0 if tool is used in LTP.
 * @param   coef    spectral coefficients
 */
static void apply_tns(AACContext * ac, int decode, SingleChannelElement * sce, float coef[1024]) {
    const IndividualChannelStream * ics = &sce->ics;
    const TemporalNoiseShaping * tns = &sce->tns;
    const int mmm = FFMIN(ics->tns_max_bands,  ics->max_sfb);
    int w, filt, m, i, ib;
    int bottom, top, order, start, end, size, inc;
    float tmp;
    float lpc[TNS_MAX_ORDER + 1], b[2 * TNS_MAX_ORDER];

    for (w = 0; w < ics->num_windows; w++) {
        bottom = ics->num_swb;
        for (filt = 0; filt < tns->n_filt[w]; filt++) {
            top = bottom;
            bottom = FFMAX(                  0, top - tns->length[w][filt]);
            order  = FFMIN(tns->order[w][filt], TNS_MAX_ORDER);
            if (order == 0)
                continue;

            // tns_decode_coef
            lpc[0] = 1;
            for (m = 1; m <= order; m++) {
                lpc[m] = tns->tmp2_map[w][filt][tns->coef[w][filt][m - 1]];
                for (i = 1; i < m; i++)
                    b[i] = lpc[i] + lpc[m] * lpc[m-i];
                for (i = 1; i < m; i++)
                    lpc[i] = b[i];
            }

            start = ics->swb_offset[FFMIN(bottom, mmm)];
            end   = ics->swb_offset[FFMIN(   top, mmm)];
            if ((size = end - start) <= 0)
                continue;
            if (tns->direction[w][filt]) {
                inc = -1; start = end - 1;
            } else {
                inc = 1;
            }
            start += w * 128;

            // ar filter
            memset(b, 0, sizeof(b));
            ib = 0;
            for (m = 0; m < size; m++) {
                tmp = coef[start];
                if (decode) {
                    for (i = 0; i < order; i++)
                        tmp -= b[ib + i] * lpc[i + 1];
                } else { // encode
                    for (i = 0; i < order; i++)
                        tmp += b[i]      * lpc[i + 1];
                }
                if (--ib < 0)
                    ib = order - 1;
                b[ib] = b[ib + order] = tmp;
                coef[start] = tmp;
                start += inc;
            }
        }
    }
}

#ifdef AAC_LTP
static void windowing_and_mdct_ltp(AACContext * ac, SingleChannelElement * sce, float * in, float * out) {
    IndividualChannelStream * ics = &sce->ics;
    const float * lwindow      = ics->use_kb_window[0] ? kbd_long_1024 : sine_long_1024;
    const float * swindow      = ics->use_kb_window[0] ? kbd_short_128 : sine_short_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? kbd_long_1024 : sine_long_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? kbd_short_128 : sine_short_128;
    float * buf = ac->buf_mdct;
    assert(ics->window_sequence[0] != EIGHT_SHORT_SEQUENCE);
    if (ics->window_sequence[0] != LONG_STOP_SEQUENCE) {
        vector_fmul_dst(ac, buf, in, lwindow_prev, 1024);
    } else {
        memset(buf, 0, 448 * sizeof(float));
        vector_fmul_dst(ac, buf + 448, in + 448, swindow_prev, 128);
        memcpy(buf + 576, in + 576, 448 * sizeof(float));
    }
    if (ics->window_sequence[0] != LONG_START_SEQUENCE) {
        ac->dsp.vector_fmul_reverse(buf + 1024, in + 1024, lwindow, 1024);
    } else {
        memcpy(buf + 1024, in + 1024, 448 * sizeof(float));
        ac->dsp.vector_fmul_reverse(buf + 1024 + 448, in + 1024 + 448, swindow, 128);
        memset(buf + 1024 + 576, 0, 448 * sizeof(float));
    }
    ff_mdct_calc(&ac->mdct_ltp, out, buf, in); // Using in as buffer for MDCT.
}

static int apply_ltp(AACContext * ac, SingleChannelElement * sce) {
    const LongTermPrediction * ltp = &sce->ics.ltp;
    const uint16_t * offsets = sce->ics.swb_offset;
    int i, sfb;
    if (!ltp->present)
        return 0;
    if (!sce->ltp_state && !(sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t))))
        return AVERROR(ENOMEM);
    if (sce->ics.window_sequence[0] != EIGHT_SHORT_SEQUENCE && ac->is_saved) {
        float x_est[2 * 1024], X_est[2 * 1024];
        for (i = 0; i < 2 * 1024; i++)
            x_est[i] = sce->ltp_state[i + 2 * 1024 - ltp->lag] * ltp->coef;

        windowing_and_mdct_ltp(ac, sce, x_est, X_est);
        if(sce->tns.present) apply_tns(ac, 0, sce, X_est);

        for (sfb = 0; sfb < FFMIN(sce->ics.max_sfb, MAX_LTP_LONG_SFB); sfb++)
            if (ltp->used[sfb])
                for (i = offsets[sfb]; i < offsets[sfb + 1]; i++)
                    sce->coeffs[i] += X_est[i];
    }
    return 0;
}


/**
 * @todo: Replace this with float_to_int16().
 */
static inline int16_t ltp_round(float x) {
    if (x >= 0)
    {
        if (x >= 1.0f)
            return 32767;
    } else {
        if (x <= -1.0f)
            return -32768;
    }

    return lrintf(32768 * x);
}


static int update_ltp(AACContext * ac, SingleChannelElement * sce) {
    int i;
    if (!sce->ltp_state && !(sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t))))
        return AVERROR(ENOMEM);
    if (ac->is_saved) {
        for (i = 0; i < 1024; i++) {
            sce->ltp_state[i] = sce->ltp_state[i + 1024];
            sce->ltp_state[i + 1024] = ltp_round(sce->ret[i] - ac->add_bias);
            sce->ltp_state[i + 2 * 1024] = ltp_round(sce->saved[i]);
            //sce->ltp_state[i + 3 * 1024] = 0;
        }
    }
    return 0;
}
#endif /* AAC_LTP */

/**
 * Conduct IMDCT and windowing.
 */
static void imdct_and_windowing(AACContext * ac, SingleChannelElement * sce) {
    IndividualChannelStream * ics = &sce->ics;
    float * in = sce->coeffs;
    float * out = sce->ret;
    float * saved = sce->saved;
    const float * lwindow      = ics->use_kb_window[0] ? kbd_long_1024 : sine_long_1024;
    const float * swindow      = ics->use_kb_window[0] ? kbd_short_128 : sine_short_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? kbd_long_1024 : sine_long_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? kbd_short_128 : sine_short_128;
    float * buf = ac->buf_mdct;
    int i;

    if (ics->window_sequence[0] != EIGHT_SHORT_SEQUENCE) {
        ff_imdct_calc(&ac->mdct, buf, in, out); // Out can be abused, for now, as a temp buffer.
        if (ics->window_sequence[0] != LONG_STOP_SEQUENCE) {
            ac->dsp.vector_fmul_add_add(out, buf, lwindow_prev, saved, ac->add_bias, 1024, 1);
        } else {
            for (i = 0;   i < 448;  i++)   out[i] =          saved[i] + ac->add_bias;
            ac->dsp.vector_fmul_add_add(out + 448, buf + 448, swindow_prev, saved + 448, ac->add_bias, 128, 1);
            for (i = 576; i < 1024; i++)   out[i] = buf[i] + saved[i] + ac->add_bias;
        }
        if (ics->window_sequence[0] != LONG_START_SEQUENCE) {
            ac->dsp.vector_fmul_reverse(saved, buf + 1024, lwindow, 1024);
        } else {
            memcpy(saved, buf + 1024, 448 * sizeof(float));
            ac->dsp.vector_fmul_reverse(saved + 448, buf + 1024 + 448, swindow, 128);
            memset(saved + 576, 0, 448 * sizeof(float));
        }
    } else {
        if (ics->window_sequence[1] == ONLY_LONG_SEQUENCE || ics->window_sequence[1] == LONG_STOP_SEQUENCE)
            av_log(ac->avccontext, AV_LOG_WARNING,
                   "Transition from an ONLY_LONG or LONG_STOP to an EIGHT_SHORT sequence detected. "
                   "If you heard an audible artifact, please submit the sample to the FFmpeg developers.\n");
        for (i = 0; i < 2048; i += 256) {
            ff_imdct_calc(&ac->mdct_small, buf + i, in + i/2, out);
            ac->dsp.vector_fmul_reverse(ac->revers + i/2, buf + i + 128, swindow, 128);
        }
        for (i = 0; i < 448; i++)   out[i] = saved[i] + ac->add_bias;

        ac->dsp.vector_fmul_add_add(out + 448 + 0*128, buf + 0*128, swindow_prev, saved + 448 ,       ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 1*128, buf + 2*128, swindow,      ac->revers + 0*128, ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 2*128, buf + 4*128, swindow,      ac->revers + 1*128, ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 3*128, buf + 6*128, swindow,      ac->revers + 2*128, ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 4*128, buf + 8*128, swindow,      ac->revers + 3*128, ac->add_bias,  64, 1);

#if 0
        vector_fmul_add_add_add(ac, out + 448 + 1*128, buf + 2*128, swindow,      saved + 448 + 1*128, ac->revers + 0*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 448 + 2*128, buf + 4*128, swindow,      saved + 448 + 2*128, ac->revers + 1*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 448 + 3*128, buf + 6*128, swindow,      saved + 448 + 3*128, ac->revers + 2*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 448 + 4*128, buf + 8*128, swindow,      saved + 448 + 4*128, ac->revers + 3*128, ac->add_bias, 64);
#endif

        ac->dsp.vector_fmul_add_add(saved,       buf + 1024 + 64,    swindow + 64, ac->revers + 3*128+64,  0, 64, 1);
        ac->dsp.vector_fmul_add_add(saved + 64,  buf + 1024 + 2*128, swindow,      ac->revers + 4*128,     0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 192, buf + 1024 + 4*128, swindow,      ac->revers + 5*128,     0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 320, buf + 1024 + 6*128, swindow,      ac->revers + 6*128,     0, 128, 1);
        memcpy(                     saved + 448, ac->revers + 7*128, 128 * sizeof(float));
        memset(                     saved + 576, 0,                  448 * sizeof(float));
    }
}

#ifdef AAC_SSR
static void windowing_and_imdct_ssr(AACContext * ac, SingleChannelElement * sce, float * in, float * out) {
    IndividualChannelStream * ics = &sce->ics;
    const float * lwindow      = ics->use_kb_window[0] ? kbd_long_1024 : sine_long_1024;
    const float * swindow      = ics->use_kb_window[0] ? kbd_short_128 : sine_short_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? kbd_long_1024 : sine_long_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? kbd_short_128 : sine_short_128;
    float * buf = ac->buf_mdct;
    if (ics->window_sequence[0] != EIGHT_SHORT_SEQUENCE) {
        ff_imdct_calc(&ac->mdct, buf, in, out);
        if (ics->window_sequence[0] != LONG_STOP_SEQUENCE) {
            vector_fmul_dst(ac, out, buf, lwindow_prev, 256);
        } else {
            memset(out, 0, 112 * sizeof(float));
            vector_fmul_dst(ac, out + 112, buf + 112, swindow_prev, 32);
            memcpy(out + 144, buf + 144, 112 * sizeof(float));
        }
        if (ics->window_sequence[0] != LONG_START_SEQUENCE) {
            ac->dsp.vector_fmul_reverse(out + 256, buf + 256, lwindow, 256);
        } else {
            memcpy(out + 256, buf + 256, 112 * sizeof(float));
            ac->dsp.vector_fmul_reverse(out + 256 + 112, buf + 256 + 112, swindow, 32);
            memset(out + 144, 0, 112 * sizeof(float));
        }
    } else {
        int i;
        for (i = 0; i < 8; i++) {
            ff_imdct_calc(&ac->mdct_small, buf, in + i * 32, out);
            vector_fmul_dst(ac, out + 64 * i, buf, (i == 0) ? swindow_prev : swindow, 32);
            ac->dsp.vector_fmul_reverse(out + 64 * i + 32, buf + 32, swindow, 32);
        }
    }
}

static void vector_add_dst(AACContext * ac, float * dst, const float * src0, const float * src1, int len) {
    int i;
    for (i = 0; i < len; i++)
        dst[i] = src0[i] + src1[i];
}

static void apply_ssr_gains(AACContext * ac, SingleChannelElement * sce, int band, float * in, float * preret, float * saved) {
    // TODO: 'in' buffer gain normalization
    if (sce->ics.window_sequence[0] != EIGHT_SHORT_SEQUENCE) {
        vector_add_dst(ac, preret, in, saved, 256);
        memcpy(saved, in + 256, 256 * sizeof(float));
    } else {
        memcpy(preret, saved, 112 * sizeof(float));
        preret += 112; saved += 112;
        vector_add_dst(ac, preret       , in            , saved    , 32);
        vector_add_dst(ac, preret + 1*32, in + 0*64 + 32, in + 1*64, 32);
        vector_add_dst(ac, preret + 2*32, in + 1*64 + 32, in + 2*64, 32);
        vector_add_dst(ac, preret + 3*32, in + 2*64 + 32, in + 3*64, 32);
        vector_add_dst(ac, preret + 4*32, in + 3*64 + 32, in + 4*64, 16);

        vector_add_dst(ac, saved            , in + 3*64 + 32 + 16, in + 4*64 + 16, 16);
        vector_add_dst(ac, saved        + 16, in + 4*64 + 32     , in + 5*64     , 32);
        vector_add_dst(ac, saved + 1*32 + 16, in + 5*64 + 32     , in + 6*64     , 32);
        vector_add_dst(ac, saved + 2*32 + 16, in + 6*64 + 32     , in + 7*64     , 32);
        memcpy(saved + 3*32 + 16, in + 7*64 + 32, 32 * sizeof(float));
        memset(saved + 144, 0, 112 * sizeof(float));
    }
}

static void inverse_polyphase_quadrature_filter(AACContext * ac, SingleChannelElement * sce, float * preret) {
    ssr_context * ctx = &ac->ssrctx;
    ScalableSamplingRate * ssr = sce->ssr;
    int i, b, j;
    float x;
    for (i = 0; i < 256; i++) {
        memcpy(&ssr->buf[0][0], &ssr->buf[0][1], 23 * sizeof(float));
        memcpy(&ssr->buf[1][0], &ssr->buf[1][1], 23 * sizeof(float));
        memcpy(&ssr->buf[2][0], &ssr->buf[2][1], 23 * sizeof(float));
        memcpy(&ssr->buf[3][0], &ssr->buf[3][1], 23 * sizeof(float));

        ssr->buf[0][23] = ctx->q[0][0] * preret[0*256+i] + ctx->q[0][1] * preret[1*256+i] +
            ctx->q[0][2] * preret[2*256+i] + ctx->q[0][3] * preret[3*256+i];
        ssr->buf[1][23] = ctx->q[1][0] * preret[0*256+i] + ctx->q[1][1] * preret[1*256+i] +
            ctx->q[1][2] * preret[2*256+i] + ctx->q[1][3] * preret[3*256+i];
        ssr->buf[2][23] = ctx->q[2][0] * preret[0*256+i] + ctx->q[2][1] * preret[1*256+i] +
            ctx->q[2][2] * preret[2*256+i] + ctx->q[2][3] * preret[3*256+i];
        ssr->buf[3][23] = ctx->q[3][0] * preret[0*256+i] + ctx->q[3][1] * preret[1*256+i] +
            ctx->q[3][2] * preret[2*256+i] + ctx->q[3][3] * preret[3*256+i];

        for (b = 0; b < 2; b++) {
            x = 0.0;
            for (j = 0; j < 12; j++)
                x += ctx->t0[b][j] * ssr->buf[b][23-2*j] + ctx->t1[b][j] * ssr->buf[b+2][22-2*j];
            sce->ret[4*i + b] = x + ac->add_bias;
            x = 0.0;
            for (j = 0; j < 12; j++)
                x += ctx->t0[3-b][j] * ssr->buf[b][23-2*j] - ctx->t1[3-b][j] * ssr->buf[b+2][22-2*j];
            sce->ret[4*i + 3-b] = x + ac->add_bias;
        }
    }
}

static void vector_reverse(AACContext * ac, float * dst, const float * src, int len) {
    int i;
    for (i = 0; i < len; i++)
        dst[i] = src[len - i];
}

static void apply_ssr(AACContext * ac, SingleChannelElement * sce) {
    float * in = sce->coeffs;
    DECLARE_ALIGNED_16(float, tmp_buf[512]);
    DECLARE_ALIGNED_16(float, tmp_ret[1024]);
    int b;
    for (b = 0; b < 4; b++) {
        if (b & 1) { // spectral reverse
            vector_reverse(ac, tmp_buf, in + 256 * b, 256);
            memcpy(in + 256 * b, tmp_buf, 256 * sizeof(float));
        }
        windowing_and_imdct_ssr(ac, sce, in + 256 * b, tmp_buf);
        apply_ssr_gains(ac, sce, b, tmp_buf, tmp_ret + 256 * b, sce->saved + 256 * b);
    }
    inverse_polyphase_quadrature_filter(ac, sce, tmp_ret);
}
#endif /* AAC_SSR */

/**
 * Apply dependent channel coupling (applied before IMDCT).
 *
 * @param   index   which gain to use for coupling
 */
static void dependent_coupling(AACContext * ac, ChannelElement * cc, SingleChannelElement * sce, int index) {
    IndividualChannelStream * ics = &cc->ch[0].ics;
    const uint16_t * offsets = ics->swb_offset;
    float * dest = sce->coeffs;
    const float * src = cc->ch[0].coeffs;
    int g, i, group, k;
    if(ac->m4ac.object_type == AOT_AAC_LTP) {
        av_log(ac->avccontext, AV_LOG_ERROR,
               "Dependent coupling is not supported together with LTP\n");
        return;
    }
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            if (cc->ch[0].band_type[g][i] != ZERO_BT) {
                float gain = cc->coup.gain[index][g][i] * sce->mixing_gain;
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++) {
                        // XXX dsputil-ize
                        dest[group*128+k] += gain * src[group*128+k];
                    }
                }
            }
        }
        dest += ics->group_len[g]*128;
        src  += ics->group_len[g]*128;
    }
}

/**
 * Apply independent channel coupling (applied after IMDCT).
 *
 * @param   index   which gain to use for coupling
 */
static void independent_coupling(AACContext * ac, ChannelElement * cc, SingleChannelElement * sce, int index) {
    int i;
    float gain = cc->coup.gain[index][0][0] * sce->mixing_gain;
    for (i = 0; i < 1024; i++)
        sce->ret[i] += gain * (cc->ch[0].ret[i] - ac->add_bias);
}

/**
 * channel coupling transformation interface
 *
 * @param   index   which gain to use for coupling
 */
static void apply_channel_coupling(AACContext * ac, ChannelElement * cc,
        void (*apply_coupling_method)(AACContext * ac, ChannelElement * cc, SingleChannelElement * sce, int index))
{
    int c;
    int index = 0;
    ChannelCoupling * coup = &cc->coup;
    for (c = 0; c <= coup->num_coupled; c++) {
        if (     !coup->is_cpe[c] && ac->che[ID_SCE][coup->tag_select[c]]) {
            apply_coupling_method(ac, cc, &ac->che[ID_SCE][coup->tag_select[c]]->ch[0], index++);
        } else if(coup->is_cpe[c] && ac->che[ID_CPE][coup->tag_select[c]]) {
            if (!coup->l[c] && !coup->r[c]) {
                apply_coupling_method(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[0], index);
                apply_coupling_method(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[1], index++);
            }
            if (coup->l[c])
                apply_coupling_method(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[0], index++);
            if (coup->r[c])
                apply_coupling_method(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[1], index++);
        } else {
            av_log(ac->avccontext, AV_LOG_ERROR,
                   "coupling target %sE[%d] not available\n",
                   coup->is_cpe[c] ? "CP" : "SC", coup->tag_select[c]);
            break;
        }
    }
}

/**
 * Convert spectral data to float samples, applying all supported tools as appropriate.
 */
static int spectral_to_sample(AACContext * ac) {
    int i, j;
    for (i = 0; i < MAX_TAGID; i++) {
        for(j = 0; j < 4; j++) {
            ChannelElement *che = ac->che[j][i];
            if(che) {
                if(j == ID_CCE && !che->coup.is_indep_coup && (che->coup.domain == 0))
                    apply_channel_coupling(ac, che, dependent_coupling);
#ifdef AAC_LTP
                if (ac->m4ac.object_type == AOT_AAC_LTP) {
                    int ret;
                    if((ret = apply_ltp(ac, &che->ch[0])))
                        return ret;
                    if(j == ID_CPE && (ret = apply_ltp(ac, &che->ch[1])))
                        return ret;
                }
#endif /* AAC_LTP */
                if(che->ch[0].tns.present)
                    apply_tns(ac, 1, &che->ch[0], che->ch[0].coeffs);
                if(che->ch[1].tns.present)
                    apply_tns(ac, 1, &che->ch[1], che->ch[1].coeffs);
                if(j == ID_CCE && !che->coup.is_indep_coup && (che->coup.domain == 1))
                    apply_channel_coupling(ac, che, dependent_coupling);
#ifdef AAC_SSR
                if (ac->m4ac.object_type == AOT_AAC_SSR) {
                    apply_ssr(ac, &che->ch[0]);
                    if(j == ID_CPE)
                        apply_ssr(ac, &che->ch[1]);
                } else {
#endif /* AAC_SSR */
                    imdct_and_windowing(ac, &che->ch[0]);
                    if(j == ID_CPE)
                        imdct_and_windowing(ac, &che->ch[1]);
#ifdef AAC_SSR
                }
#endif /* AAC_SSR */
                if(j == ID_CCE && che->coup.is_indep_coup && (che->coup.domain == 1))
                    apply_channel_coupling(ac, che, independent_coupling);
#ifdef AAC_LTP
                if (ac->m4ac.object_type == AOT_AAC_LTP) {
                    int ret;
                    if((ret = update_ltp(ac, &che->ch[0])))
                        return ret;
                    if(j == ID_CPE && (ret = update_ltp(ac, &che->ch[1])))
                        return ret;
                }
#endif /* AAC_LTP */
            }
        }
    }
    return 0;
}

/**
 * Conduct matrix mix-down and float to int16 conversion.
 *
 * @param   data        pointer to output data
 * @param   data_size   output data size in bytes
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int output_samples(AVCodecContext * avccontext, uint16_t * data, int * data_size) {
    AACContext * ac = avccontext->priv_data;
    int i;
    float *c, *l, *r, *sl, *sr, *out;

    if (!ac->is_saved) {
        ac->is_saved = 1;
        *data_size = 0;
        return 0;
    }

    i = 1024 * avccontext->channels * sizeof(int16_t);
    if(*data_size < i) {
        av_log(avccontext, AV_LOG_ERROR,
               "Output buffer too small (%d) or trying to output too many samples (%d) for this frame.\n",
               *data_size, i);
        return -1;
    }
    *data_size = i;

    if(ac->mm[MIXDOWN_CENTER]) {
        /* matrix mix-down */
        l   = ac->mm[MIXDOWN_FRONT ]->ch[0].ret;
        r   = ac->mm[MIXDOWN_FRONT ]->ch[1].ret;
        c   = ac->mm[MIXDOWN_CENTER]->ch[0].ret;
        sl  = ac->mm[MIXDOWN_BACK  ]->ch[0].ret;
        sr  = ac->mm[MIXDOWN_BACK  ]->ch[1].ret;
        out = ac->interleaved_output;

        // XXX dsputil-ize
        if(avccontext->channels == 2) {
            if(ac->pcs.pseudo_surround) {
                for(i = 0; i < 1024; i++) {
                    *out++ = *l++ + *c   - *sl   - *sr   + ac->add_bias;
                    *out++ = *r++ + *c++ + *sl++ + *sr++ - ac->add_bias * 3;
                }
            } else {
                for(i = 0; i < 1024; i++) {
                    *out++ = *l++ + *c   + *sl++ - ac->add_bias * 2;
                    *out++ = *r++ + *c++ + *sr++ - ac->add_bias * 2;
                }
            }

        } else {
            assert(avccontext->channels == 1);
            for(i = 0; i < 1024; i++) {
                *out++ = *l++ + *r++ + *c++ + *sl++ + *sr++ - ac->add_bias * 4;
            }
        }

        ac->dsp.float_to_int16(data, ac->interleaved_output, 1024 * avccontext->channels);
    } else {
        ac->dsp.float_to_int16_interleave(data, ac->output_data, 1024, avccontext->channels);
    }

    return 0;
}


static int aac_decode_frame(AVCodecContext * avccontext, void * data, int * data_size, const uint8_t * buf, int buf_size) {
    AACContext * ac = avccontext->priv_data;
    GetBitContext gb;
    enum RawDataBlockID id;
    int err, tag;

    init_get_bits(&gb, buf, buf_size*8);

    // parse
    while ((id = get_bits(&gb, 3)) != ID_END) {
        tag = get_bits(&gb, 4);
        err = -1;
        switch (id) {

        case ID_SCE:
            if(!ac->che[ID_SCE][tag]) {
                if(tag == 1 && ac->che[ID_LFE][0]) {
                    /* Some streams incorrectly code 5.1 audio as SCE[0] CPE[0] CPE[1] SCE[1]
                       instead of SCE[0] CPE[0] CPE[0] LFE[0].
                       If we seem to have encountered such a stream,
                       transfer the LFE[0] element to SCE[1] */
                    ac->che[ID_SCE][tag] = ac->che[ID_LFE][0];
                    ac->che[ID_LFE][0] = NULL;
                } else
                    break;
            }
            err = decode_ics(ac, &gb, 0, 0, &ac->che[ID_SCE][tag]->ch[0]);
            break;

        case ID_CPE:
            if (ac->che[ID_CPE][tag])
                err = decode_cpe(ac, &gb, tag);
            break;

        case ID_FIL:
            if (tag == 15)
                tag += get_bits(&gb, 8) - 1;
            while (tag > 0)
                tag -= extension_payload(ac, &gb, tag);
            err = 0; /* FIXME */
            break;

        case ID_PCE:
            err = program_config_element(ac, &gb);
            break;

        case ID_DSE:
            data_stream_element(ac, &gb);
            err = 0;
            break;

        case ID_CCE:
            if (ac->che[ID_CCE][tag])
                err = decode_cce(ac, &gb, tag);
            break;

        case ID_LFE:
            if (ac->che[ID_LFE][tag])
                err = decode_ics(ac, &gb, 0, 0, &ac->che[ID_LFE][tag]->ch[0]);
            break;

        default:
            err = -1; /* should not happen, but keeps compiler happy */
            break;
        }

        if(err)
            return err;
    }

    if((err = spectral_to_sample(ac)))
        return err;
    output_samples(avccontext, data, data_size);

    return buf_size;
}

static av_cold int aac_decode_close(AVCodecContext * avccontext) {
    AACContext * ac = avccontext->priv_data;
    int i, j;

    for (i = 0; i < MAX_TAGID; i++) {
        for(j = 0; j < 4; j++)
            che_freep(&ac->che[j][i]);
    }

    ff_mdct_end(&ac->mdct);
    ff_mdct_end(&ac->mdct_small);
#ifdef AAC_LTP
    ff_mdct_end(&ac->mdct_ltp);
#endif /* AAC_LTP */
    av_freep(&ac->interleaved_output);
    return 0 ;
}

AVCodec aac_decoder = {
    "aac",
    CODEC_TYPE_AUDIO,
    CODEC_ID_AAC,
    sizeof(AACContext),
    aac_decode_init,
    NULL,
    aac_decode_close,
    aac_decode_frame,
    .long_name = NULL_IF_CONFIG_SMALL("Advanced Audio Coding"),
};

