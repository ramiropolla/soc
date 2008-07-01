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


#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"
#include "libavutil/random.h"

#include "aactab.h"
#include "mpeg4audio.h"

#include <assert.h>

#define MAX_CHANNELS 64
#define MAX_TAGID 16

#ifndef CONFIG_HARDCODED_TABLES
    static float ivquant_tab[IVQUANT_SIZE];
    static float pow2sf_tab[316];
#endif /* CONFIG_HARDCODED_TABLES */
DECLARE_ALIGNED_16(static float, kbd_long_1024[1024]);
DECLARE_ALIGNED_16(static float, kbd_short_128[128]);
DECLARE_ALIGNED_16(static float, sine_long_1024[1024]);
DECLARE_ALIGNED_16(static float, sine_short_128[128]);

static VLC mainvlc;
static VLC books[11];

/**
 * Audio Object Types
 */
enum {
    AOT_NULL,
    AOT_AAC_MAIN,
    AOT_AAC_LC,
    AOT_AAC_SSR,
    AOT_AAC_LTP,
    AOT_SBR,
    AOT_AAC_SCALABLE,
    AOT_TWINVQ,
    AOT_CELP,
    AOT_HVXC,
    AOT_TTSI             = 12,
    AOT_MAINSYNTH,
    AOT_WAVESYNTH,
    AOT_MIDI,
    AOT_SAFX,
    AOT_ER_AAC_LC,
    AOT_ER_AAC_LTP       = 19,
    AOT_ER_AAC_SCALABLE,
    AOT_ER_TWINVQ,
    AOT_ER_BSAC,
    AOT_ER_AAC_LD,
    AOT_ER_CELP,
    AOT_ER_HVXC,
    AOT_ER_HILN,
    AOT_ER_PARAM,
    AOT_SSC,
};

/**
 * IDs for raw_data_block
 */
enum {
    ID_SCE,
    ID_CPE,
    ID_CCE,
    ID_LFE,
    ID_DSE,
    ID_PCE,
    ID_FIL,
    ID_END,
};

/**
 * IDs for extension_payload
 */
enum {
    EXT_FILL,
    EXT_FILL_DATA,
    EXT_DATA_ELEMENT,
    EXT_DYNAMIC_RANGE = 0xb,
    EXT_SBR_DATA      = 0xd,
    EXT_SBR_DATA_CRC  = 0xe,
};

/**
 * window sequences
 */
enum WindowSequence {
    ONLY_LONG_SEQUENCE,
    LONG_START_SEQUENCE,
    EIGHT_SHORT_SEQUENCE,
    LONG_STOP_SEQUENCE,
};

/**
 * special codebooks
 */
enum {
    ZERO_HCB       = 0,
    FIRST_PAIR_HCB = 5,
    ESC_HCB        = 11,
    NOISE_HCB      = 13,
    INTENSITY_HCB2 = 14,
    INTENSITY_HCB  = 15,
    ESC_FLAG       = 16,
};

#define IS_CODEBOOK_UNSIGNED(x) ((x - 1) & 10)

/**
 * channel types
 */
enum {
    AAC_CHANNEL_FRONT = 1,
    AAC_CHANNEL_SIDE  = 2,
    AAC_CHANNEL_BACK  = 3,
    AAC_CHANNEL_LFE   = 4,
    AAC_CHANNEL_CC    = 5,
};

/**
 * mix-down channel types
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
    int che_type[4][MAX_TAGID]; ///< channel element type with the first index as the first 4 raw_data_block IDs
    int mono_mixdown;           ///< The SCE tag to use if user requests mono   output, -1 if not available.
    int stereo_mixdown;         ///< The CPE tag to use if user requests stereo output, -1 if not available.
    int mixdown_coeff_index;    ///< 0-3
    int pseudo_surround;        ///< Mix surround channels out of phase.
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
    int max_sfb;                ///< number of scalefactor bands per group
    enum WindowSequence window_sequence;
    enum WindowSequence window_sequence_prev;
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
    int pce_instance_tag;
    int tag_reserved_bits;
    int dyn_rng_sgn[17];
    int dyn_rng_ctl[17];
    int exclude_mask[MAX_CHANNELS];
    int additional_excluded_chns[MAX_CHANNELS];
    int band_incr;
    int interpolation_scheme;
    int band_top[17];
    int prog_ref_level;
    int prog_ref_level_reserved_bits;
} DynamicRangeControl;

/**
 * pulse tool
 */
typedef struct {
    int present;
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
    int ind_sw;            ///< Set if independent coupling (i.e. after IMDCT).
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
    int cb[8][64];                            ///< codebooks
    int cb_run_end[8][64];                    ///< codebook run end points
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
    int common_window;     ///< Set if channels share a common 'IndividualChannelStream' in bitstream.
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
    MDCTContext *mdct_ltp;
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
    float* interleaved_output;                        ///< Interim buffer for interleaving PCM samples.
    float *output_data[MAX_CHANNELS];                 ///< Points to each element's 'ret' buffer (PCM output).
    ChannelElement *mm[3];                            ///< Center/Front/Back channel elements to use for matrix mix-down.
    float add_bias;                                   ///< offset for dsp.float_to_int16
    float sf_scale;                                   ///< Pre-scale for correct IMDCT and dsp.float_to_int16.
    int sf_offset;                                    ///< offset into pow2sf_tab as appropriate for dsp.float_to_int16
    /** @} */

} AACContext;


//aux
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
    if(!*s)
        return;
#ifdef AAC_SSR
    av_free((*s)->ch[0].ssr);
    av_free((*s)->ch[1].ssr);
#endif /* AAC_SSR */
#ifdef AAC_LTP
    av_free((*s)->ch[0].ltp_state);
    av_free((*s)->ch[1].ltp_state);
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
            if(pcs->che_type[j][i] && !ac->che[j][i]) {
                ac->che[j][i] = av_mallocz(sizeof(ChannelElement));
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

    ac->mm[MIXDOWN_FRONT] = ac->mm[MIXDOWN_BACK] = ac->mm[MIXDOWN_CENTER] = NULL;

    /* Check for matrix mix-down to mono or stereo. */

    if(avctx->request_channels && avctx->request_channels <= 2 &&
       avctx->request_channels != channels) {

        if((avctx->request_channels == 1 && pcs->mono_mixdown   != -1) ||
           (avctx->request_channels == 2 && pcs->stereo_mixdown != -1)) {
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
    ac->interleaved_output = av_realloc(ac->interleaved_output, channels * 1024 * sizeof(float));
    return ac->interleaved_output ? 0 : -1;
}


/**
 * Decode an array of 4 bit tag IDs, optionally interleaved with a stereo/mono switching bit.
 *
 * @param cpe_map Stereo (Channel Pair Element) map, NULL if stereo bit is not present.
 * @param sce_map mono (Single Channel Element) map
 * @param type speaker type/position for these channels
 */
static void program_config_element_parse_tags(GetBitContext * gb, int *cpe_map,
                                              int *sce_map, int n, int type) {
    int *map;
    while(n--) {
        map = cpe_map && get_bits1(gb) ? cpe_map : sce_map; // stereo or mono map
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

    pcs.mono_mixdown   = get_bits1(gb) ? get_bits(gb, 4) : -1;
    pcs.stereo_mixdown = get_bits1(gb) ? get_bits(gb, 4) : -1;

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
    pcs.mono_mixdown   = -1;
    pcs.stereo_mixdown = -1;

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
static int GASpecificConfig(AACContext * ac, GetBitContext * gb, int channels) {
    int ext;

    if(get_bits1(gb)) {  // frameLengthFlag
        av_log(ac->avccontext, AV_LOG_ERROR, "960/120 MDCT window is not supported.\n");
        return -1;
    }

    if (get_bits1(gb))       // dependsOnCoreCoder
        skip_bits(gb, 14);   // coreCoderDelay
    ext = get_bits1(gb);

    if(ac->m4ac.object_type == AOT_AAC_SCALABLE ||
       ac->m4ac.object_type == AOT_ER_AAC_SCALABLE)
        skip_bits(gb, 3);     // layerNr

    if (channels == 0) {
        skip_bits(gb, 4);  // element_instance_tag
        if(program_config_element(ac, gb) < 0)
            return -1;
    } else {
        if(program_config_element_default(ac, channels) < 0)
            return -1;
    }

    if (ext) {
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
 */
static int AudioSpecificConfig(AACContext * ac, void *data, int data_size) {
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
        if (GASpecificConfig(ac, &gb, ac->m4ac.chan_config))
            return -1;
        break;
    default:
        av_log(ac->avccontext, AV_LOG_ERROR, "Audio object type %s%d is not supported.\n",
               ac->m4ac.sbr == 1? "SBR+" : "", ac->m4ac.object_type);
        return -1;
    }
    return 0;
}

static int aac_decode_init(AVCodecContext * avccontext) {
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

    if (AudioSpecificConfig(ac, avccontext->extradata, avccontext->extradata_size))
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

    /* Initialize RNG dither. */
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

#ifdef AAC_SSR
    if (ac->audioObjectType == AOT_AAC_SSR) {
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

// parser implementation

/**
 * Decode a data_stream_element; reference: table 4.10.
 */
static void data_stream_element(AACContext * ac, GetBitContext * gb, int id) {
    int byte_align = get_bits1(gb);
    int count = get_bits(gb, 8);
    if (count == 255)
        count += get_bits(gb, 8);
    if (byte_align)
        align_get_bits(gb);
    skip_bits_long(gb, 8 * count);
}

#ifdef AAC_LTP
static void decode_ltp_data(AACContext * ac, GetBitContext * gb, int max_sfb, LongTermPrediction * ltp) {
    int sfb;
    if (ac->audioObjectType == AOT_ER_AAC_LD) {
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
 */
static int decode_ics_info(AACContext * ac, GetBitContext * gb, int common_window, IndividualChannelStream * ics) {
    uint8_t grouping;
    if (get_bits1(gb)) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Reserved bit set.\n");
        return -1;
    }
    ics->window_sequence_prev = ics->window_sequence;
    ics->window_sequence = get_bits(gb, 2);
    ics->use_kb_window[1] = ics->use_kb_window[0];
    ics->use_kb_window[0] = get_bits1(gb);
    ics->num_window_groups = 1;
    ics->group_len[0] = 1;
    if (ics->window_sequence == EIGHT_SHORT_SEQUENCE) {
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
        if(ics->max_sfb > ics->num_swb) {
            av_log(ac->avccontext, AV_LOG_ERROR,
                "Number of scalefactor bands in group (%d) exceeds limit (%d).\n",
                ics->max_sfb, ics->num_swb);
            return -1;
        }

        ics->num_windows   = 8;
        ics->tns_max_bands = tns_max_bands_128[ac->m4ac.sampling_index];
    } else {
        ics->max_sfb = get_bits(gb, 6);
        ics->swb_offset    =    swb_offset_1024[ac->m4ac.sampling_index];
        ics->num_swb       =       num_swb_1024[ac->m4ac.sampling_index];
        if(ics->max_sfb > ics->num_swb) {
            av_log(ac->avccontext, AV_LOG_ERROR,
                "Number of scalefactor bands in group (%d) exceeds limit (%d).\n",
                ics->max_sfb, ics->num_swb);
            return -1;
        }

        ics->num_windows   = 1;
        ics->tns_max_bands = tns_max_bands_1024[ac->m4ac.sampling_index];
        if (get_bits1(gb)) {
#ifdef AAC_LTP
            if (ac->audioObjectType == AOT_AAC_MAIN) {
                assert(0);
            } else {
                if ((ics->ltp.present = get_bits(gb, 1))) {
                    decode_ltp_data(ac, gb, ics->max_sfb, &ics->ltp);
                }
                if (common_window) {
                    if ((ics->ltp2.present = get_bits(gb, 1))) {
                        decode_ltp_data(ac, gb, ics->max_sfb, &ics->ltp2);
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

static inline float ivquant(AACContext * ac, int a) {
    if (a + (unsigned int)IVQUANT_SIZE/2 - 1 < (unsigned int)IVQUANT_SIZE - 1)
        return ivquant_tab[a + IVQUANT_SIZE/2 - 1];
    else
        return cbrtf(fabsf(a)) * a;
}

/**
 * Decode section_data payload; reference: table 4.46.
 */
static int decode_section_data(AACContext * ac, GetBitContext * gb, IndividualChannelStream * ics, int cb[][64], int cb_run_end[][64]) {
    int g;
    for (g = 0; g < ics->num_window_groups; g++) {
        int bits = (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 3 : 5;
        int k = 0;
        while (k < ics->max_sfb) {
            int sect_len = k;
            int sect_len_incr;
            int sect_cb = get_bits(gb, 4);
            if (sect_cb == 12) {
                av_log(ac->avccontext, AV_LOG_ERROR, "invalid codebook\n");
                return -1;
            }
            while ((sect_len_incr = get_bits(gb, bits)) == (1 << bits)-1)
                sect_len += sect_len_incr;
            sect_len += sect_len_incr;
            if (sect_len > ics->max_sfb) {
                av_log(ac->avccontext, AV_LOG_ERROR,
                    "Number of codebooks (%d) exceeds limit (%d).\n",
                    sect_len, ics->max_sfb);
                return -1;
            }
            for (; k < sect_len; k++) {
                cb[g][k] = sect_cb;
                cb_run_end[g][k] = sect_len;
            }
        }
    }
    return 0;
}

/**
 * Decode scale_factor_data; reference: table 4.47.
 */
static int decode_scale_factor_data(AACContext * ac, GetBitContext * gb, float mix_gain, unsigned int global_gain,
        IndividualChannelStream * ics, const int cb[][64], const int cb_run_end[][64], float sf[][64]) {
    const int sf_offset = ac->sf_offset + (ics->window_sequence == EIGHT_SHORT_SEQUENCE ? 12 : 0);
    int g, i;
    int offset[3] = { global_gain, global_gain - 90, 100 };
    int noise_flag = 1;
    static const char *sf_str[3] = { "Global gain", "Noise gain", "Intensity stereo position" };
    ics->intensity_present = 0;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb;) {
            if (cb[g][i] == ZERO_HCB) {
                int run_end = cb_run_end[g][i];
                for(; i < run_end; i++)
                    sf[g][i] = 0.;
                continue;
            }else if((cb[g][i] == INTENSITY_HCB) || (cb[g][i] == INTENSITY_HCB2)) {
                int run_end = cb_run_end[g][i];
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
            }else if(cb[g][i] == NOISE_HCB) {
                int run_end = cb_run_end[g][i];
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
                int run_end = cb_run_end[g][i];
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

static void decode_pulse_data(AACContext * ac, GetBitContext * gb, Pulse * pulse) {
    int i;
    pulse->num_pulse = get_bits(gb, 2) + 1;
    pulse->start = get_bits(gb, 6);
    for (i = 0; i < pulse->num_pulse; i++) {
        pulse->offset[i] = get_bits(gb, 5);
        pulse->amp[i] = get_bits(gb, 4);
    }
}

static void decode_tns_data(AACContext * ac, GetBitContext * gb, const IndividualChannelStream * ics, TemporalNoiseShaping * tns) {
    int w, filt, i, coef_len, coef_res = 0, coef_compress;
    for (w = 0; w < ics->num_windows; w++) {
        tns->n_filt[w] = get_bits(gb, ics->window_sequence == EIGHT_SHORT_SEQUENCE ? 1 : 2);

        if (tns->n_filt[w])
            coef_res = get_bits1(gb) + 3;

        for (filt = 0; filt < tns->n_filt[w]; filt++) {
            tns->length[w][filt] = get_bits(gb, ics->window_sequence == EIGHT_SHORT_SEQUENCE ? 4 : 6);

            if ((tns->order[w][filt] = get_bits(gb, ics->window_sequence == EIGHT_SHORT_SEQUENCE ? 3 : 5))) {
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
static int decode_gain_control_data(AACContext * ac, GetBitContext * gb, SingleChannelElement * sce) {
    // wd_num wd_test aloc_size
    static const int gain_mode[4][3] = {
        {1, 0, 5}, //ONLY_LONG_SEQUENCE = 0,
        {2, 1, 2}, //LONG_START_SEQUENCE,
        {8, 0, 2}, //EIGHT_SHORT_SEQUENCE,
        {2, 1, 5}, //LONG_STOP_SEQUENCE
    };
    const int mode = sce->ics.window_sequence;
    int bd, wd, ad;
    ScalableSamplingRate * ssr = sce->ssr;
    if (!ssr)
        ssr = sce->ssr = av_mallocz(sizeof(ScalableSamplingRate));
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

static void decode_ms_data(AACContext * ac, GetBitContext * gb, ChannelElement * cpe) {
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
 */
static int decode_spectral_data(AACContext * ac, GetBitContext * gb, const IndividualChannelStream * ics, const int cb[][64], int icoef[1024]) {
    int i, k, g;
    const uint16_t * offsets = ics->swb_offset;

    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            const int cur_cb = cb[g][i];
            const int dim = cur_cb >= FIRST_PAIR_HCB ? 2 : 4;
            int group;
            if (cur_cb == ZERO_HCB) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    memset(icoef + group * 128 + offsets[i], 0, (offsets[i+1] - offsets[i])*sizeof(int));
                }
            }else if (cur_cb != INTENSITY_HCB2 && cur_cb != INTENSITY_HCB) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k += dim) {
                        const int index = get_vlc2(gb, books[cur_cb - 1].table, 6, 3);
                        const int coef_idx = (group << 7) + k;
                        const int8_t *vq_ptr = &codebook_vectors[cur_cb - 1][index * dim];
                        int j;
                        for (j = 0; j < dim; j++)
                            icoef[coef_idx + j] = 1;
                        if (IS_CODEBOOK_UNSIGNED(cur_cb)) {
                            for (j = 0; j < dim; j++)
                                if (vq_ptr[j] && get_bits1(gb))
                                    icoef[coef_idx + j] = -1;
                        }
                        if (cur_cb == ESC_HCB) {
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

static void pulse_tool(AACContext * ac, const IndividualChannelStream * ics, const Pulse * pulse, int * icoef) {
    int i, off = ics->swb_offset[pulse->start];
    for (i = 0; i < pulse->num_pulse; i++) {
        off += pulse->offset[i];
        if (icoef[off] > 0)
            icoef[off] += pulse->amp[i];
        else
            icoef[off] -= pulse->amp[i];
    }
}

static void quant_to_spec_tool(AACContext * ac, const IndividualChannelStream * ics, const int * icoef,
        const int cb[][64], const float sf[][64], float * coef) {
    const uint16_t * offsets = ics->swb_offset;
    const int c = 1024/ics->num_window_groups;
    int g, i, group, k;

    for(g = 0; g < ics->num_window_groups; g++)
        memset(coef + g * 128 + offsets[ics->max_sfb], 0, sizeof(float)*(c - offsets[ics->max_sfb]));

    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            if (cb[g][i] == NOISE_HCB) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    float scale = sf[g][i] / ((offsets[i+1] - offsets[i]) * PNS_MEAN_ENERGY);
                    for (k = offsets[i]; k < offsets[i+1]; k++)
                        coef[group*128+k] = (int32_t)av_random(&ac->random_state) * scale;
                }
            } else if (cb[g][i] != INTENSITY_HCB && cb[g][i] != INTENSITY_HCB2) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++) {
                        coef[group*128+k] = ivquant(ac, icoef[group*128+k]) * sf[g][i];
                    }
                }
            }
        }
        coef += ics->group_len[g]*128;
        icoef += ics->group_len[g]*128;
    }
}

/**
 * Decode an individual_channel_stream payload; reference: table 4.44.
 */
static int decode_ics(AACContext * ac, GetBitContext * gb, int common_window, int scale_flag, SingleChannelElement * sce) {
    int icoeffs[1024];
    Pulse pulse;
    TemporalNoiseShaping * tns = &sce->tns;
    IndividualChannelStream * ics = &sce->ics;
    float * out = sce->coeffs;
    int global_gain;

    pulse.present = 0;
    global_gain = get_bits(gb, 8);

    if (!common_window && !scale_flag) {
        if (decode_ics_info(ac, gb, 0, ics) < 0)
            return -1;
    }

    if (decode_section_data(ac, gb, ics, sce->cb, sce->cb_run_end) < 0)
        return -1;
    if (decode_scale_factor_data(ac, gb, sce->mixing_gain, global_gain, ics, sce->cb, sce->cb_run_end, sce->sf) < 0)
        return -1;

    if (!scale_flag) {
        if ((pulse.present = get_bits1(gb))) {
            if (ics->window_sequence == EIGHT_SHORT_SEQUENCE) {
                av_log(ac->avccontext, AV_LOG_ERROR, "Pulse tool not allowed in eight short sequence.\n");
                return -1;
            }
            decode_pulse_data(ac, gb, &pulse);
        }
        if ((tns->present = get_bits1(gb)))
            decode_tns_data(ac, gb, ics, tns);
        if (get_bits1(gb)) {
#ifdef AAC_SSR
            if (decode_gain_control_data(ac, gb, sce)) return -1;
#else
            av_log(ac->avccontext, AV_LOG_ERROR, "SSR not supported.\n");
            return -1;
#endif
        }
    }

    if (decode_spectral_data(ac, gb, ics, sce->cb, icoeffs) < 0)
        return -1;
    if (pulse.present)
        pulse_tool(ac, ics, &pulse, icoeffs);
    quant_to_spec_tool(ac, ics, icoeffs, sce->cb, sce->sf, out);
    return 0;
}

static void ms_tool(AACContext * ac, ChannelElement * cpe) {
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
                        cpe->ch[0].cb[g][i] < NOISE_HCB && cpe->ch[1].cb[g][i] < NOISE_HCB) {
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


static void intensity_tool(AACContext * ac, ChannelElement * cpe) {
    const IndividualChannelStream * ics = &cpe->ch[1].ics;
    SingleChannelElement * sce1 = &cpe->ch[1];
    float *coef0 = cpe->ch[0].coeffs, *coef1 = cpe->ch[1].coeffs;
    const uint16_t * offsets = ics->swb_offset;
    int g, gp, i, k;
    int c;
    float scale;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (gp = 0; gp < ics->group_len[g]; gp++) {
            for (i = 0; i < ics->max_sfb; i++) {
                if (sce1->cb[g][i] == INTENSITY_HCB || sce1->cb[g][i] == INTENSITY_HCB2) {
                    c = -1 + 2 * (sce1->cb[g][i] - 14);
                    if (cpe->ms.present)
                        c *= 1 - 2 * cpe->ms.mask[g][i];
                    scale = c * sce1->sf[g][i];
                    for (k = offsets[i]; k < offsets[i+1]; k++)
                        coef1[k] = scale * coef0[k];
                }
            }
            coef0 += 128;
            coef1 += 128;
        }
    }
}

/**
 * Decode a channel_pair_element; reference: table 4.4.
 */
static int decode_cpe(AACContext * ac, GetBitContext * gb, int id) {
    int i;
    ChannelElement * cpe;

    cpe = ac->che[ID_CPE][id];
    cpe->common_window = get_bits1(gb);
    if (cpe->common_window) {
        if (decode_ics_info(ac, gb, 1, &cpe->ch[0].ics))
            return -1;
        i = cpe->ch[1].ics.use_kb_window[0];
        cpe->ch[1].ics = cpe->ch[0].ics;
        cpe->ch[1].ics.use_kb_window[1] = i;
#ifdef AAC_LTP
        cpe->ch[1].ics.ltp = cpe->ch[0].ics.ltp2;
#endif /* AAC_LTP */
        decode_ms_data(ac, gb, cpe);
    } else {
        cpe->ms.present = 0;
    }
    if (decode_ics(ac, gb, cpe->common_window, 0, &cpe->ch[0]))
        return -1;
    if (decode_ics(ac, gb, cpe->common_window, 0, &cpe->ch[1]))
        return -1;

    if (cpe->common_window)
        ms_tool(ac, cpe);

    if (cpe->ch[1].ics.intensity_present)
        intensity_tool(ac, cpe);
    return 0;
}

static int decode_cce(AACContext * ac, GetBitContext * gb, int id) {
    int num_gain = 0;
    int c, g, sfb;
    int sign;
    float scale;
    SingleChannelElement * sce;
    ChannelCoupling * coup;

    sce = &ac->che[ID_CCE][id]->ch[0];
    sce->mixing_gain = 1.0;

    coup = &ac->che[ID_CCE][id]->coup;

    coup->ind_sw = get_bits1(gb);
    //if (coup->ind_sw)
    //    av_log(ac->avccontext, AV_LOG_ERROR, "aac: independently switched coupling\n");
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

    if (decode_ics(ac, gb, 0, 0, sce))
        return -1;

    for (c = 0; c < num_gain; c++) {
        int cge = 1;
        int gain = 0;
        float gain_cache = 1.;
        if (c) {
            cge = coup->ind_sw ? 1 : get_bits1(gb);
            gain = cge ? get_vlc2(gb, mainvlc.table, 7, 3) - 60: 0;
            gain_cache = pow(scale, gain);
        }
        for (g = 0; g < sce->ics.num_window_groups; g++)
            for (sfb = 0; sfb < sce->ics.max_sfb; sfb++)
                if (sce->cb[g][sfb] == ZERO_HCB) {
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

static int sbr_extension_data(AACContext * ac, GetBitContext * gb, int crc, int cnt) {
    // TODO : sbr_extension implementation
    av_log(ac->avccontext, AV_LOG_DEBUG, "aac: SBR not yet supported.\n");
    skip_bits_long(gb, 8*cnt - 4);
    return cnt;
}


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



static int dynamic_range_info(AACContext * ac, GetBitContext * gb, int cnt) {
    int n = 1;
    int drc_num_bands = 1;
    int i;

    /* pce_tag_present? */
    if(get_bits1(gb)) {
        ac->che_drc.pce_instance_tag = get_bits(gb, 4);
        ac->che_drc.tag_reserved_bits = get_bits(gb, 4);
        n++;
    }

    /* excluded_chns_present? */
    if(get_bits1(gb)) {
        n += excluded_channels(ac, gb);
    }

    /* drc_bands_present? */
    if (get_bits1(gb)) {
        ac->che_drc.band_incr = get_bits(gb, 4);
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
        ac->che_drc.prog_ref_level_reserved_bits = get_bits1(gb);
        n++;
    }

    for (i = 0; i < drc_num_bands; i++) {
        ac->che_drc.dyn_rng_sgn[i] = get_bits1(gb);
        ac->che_drc.dyn_rng_ctl[i] = get_bits(gb, 7);
        n++;
    }

    return n;
}

/** Parse extension data (incomplete).
 */
static int extension_payload(AACContext * ac, GetBitContext * gb, int cnt) {
    int i = 0;
    int res = cnt;
    switch (get_bits(gb, 4)) { // extension type
        case EXT_SBR_DATA_CRC:
            i++;
        case EXT_SBR_DATA:
            res = sbr_extension_data(ac, gb, i, cnt);
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

static void tns_filter_tool(AACContext * ac, int decode, SingleChannelElement * sce, float * coef) {
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

static void tns_trans(AACContext * ac, SingleChannelElement * sce) {
    if(sce->tns.present) tns_filter_tool(ac, 1, sce, sce->coeffs);
}

#ifdef AAC_LTP
static void window_ltp_tool(AACContext * ac, SingleChannelElement * sce, float * in, float * out) {
    IndividualChannelStream * ics = &sce->ics;
    const float * lwindow      = ics->use_kb_window[0] ? kbd_long_1024 : sine_long_1024;
    const float * swindow      = ics->use_kb_window[0] ? kbd_short_128 : sine_short_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? kbd_long_1024 : sine_long_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? kbd_short_128 : sine_short_128;
    float * buf = ac->buf_mdct;
    int i;
    assert(ics->window_sequence != EIGHT_SHORT_SEQUENCE);
    if (!ac->mdct_ltp) {
        ac->mdct_ltp = av_malloc(sizeof(MDCTContext));
        ff_mdct_init(ac->mdct_ltp, 11, 0);
    }
    if (ics->window_sequence != LONG_STOP_SEQUENCE) {
        vector_fmul_dst(ac, buf, in, lwindow_prev, 1024);
    } else {
        memset(buf, 0, 448 * sizeof(float));
        vector_fmul_dst(ac, buf + 448, in + 448, swindow_prev, 128);
        memcpy(buf + 576, in + 576, 448 * sizeof(float));
    }
    if (ics->window_sequence != LONG_START_SEQUENCE) {
        ac->dsp.vector_fmul_reverse(buf + 1024, in + 1024, lwindow, 1024);
    } else {
        memcpy(buf + 1024, in + 1024, 448 * sizeof(float));
        ac->dsp.vector_fmul_reverse(buf + 1024 + 448, in + 1024 + 448, swindow, 128);
        memset(buf + 1024 + 576, 0, 448 * sizeof(float));
    }
    ff_mdct_calc(ac->mdct_ltp, out, buf, in); // Using in as buffer for MDCT.
}

static void ltp_trans(AACContext * ac, SingleChannelElement * sce) {
    const LongTermPrediction * ltp = &sce->ics.ltp;
    const uint16_t * offsets = sce->ics.swb_offset;
    int i, sfb;
    if (!ltp->present)
        return;
    if (!sce->ltp_state)
        sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t));
    if (sce->ics.window_sequence != EIGHT_SHORT_SEQUENCE && ac->is_saved) {
        float x_est[2 * 1024], X_est[2 * 1024];
        for (i = 0; i < 2 * 1024; i++)
            x_est[i] = sce->ltp_state[i + 2 * 1024 - ltp->lag] * ltp->coef;

        window_ltp_tool(ac, sce, x_est, X_est);
        if(sce->tns.present) tns_filter_tool(ac, 0, sce, X_est);

        for (sfb = 0; sfb < FFMIN(sce->ics.max_sfb, MAX_LTP_LONG_SFB); sfb++)
            if (ltp->used[sfb])
                for (i = offsets[sfb]; i < offsets[sfb + 1]; i++)
                    sce->coeffs[i] += X_est[i];
    }
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


static void ltp_update_trans(AACContext * ac, SingleChannelElement * sce) {
    int i;
    if (!sce->ltp_state)
        sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t));
    if (ac->is_saved) {
        for (i = 0; i < 1024; i++) {
            sce->ltp_state[i] = sce->ltp_state[i + 1024];
            sce->ltp_state[i + 1024] = ltp_round(sce->ret[i] - ac->add_bias);
            sce->ltp_state[i + 2 * 1024] = ltp_round(sce->saved[i]);
            //sce->ltp_state[i + 3 * 1024] = 0;
        }
    }
}
#endif /* AAC_LTP */

static void window_trans(AACContext * ac, SingleChannelElement * sce) {
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

    if (ics->window_sequence != EIGHT_SHORT_SEQUENCE) {
        ff_imdct_calc(&ac->mdct, buf, in, out); // Out can be abused, for now, as a temp buffer.
        if (ics->window_sequence != LONG_STOP_SEQUENCE) {
            ac->dsp.vector_fmul_add_add(out, buf, lwindow_prev, saved, ac->add_bias, 1024, 1);
        } else {
            for (i = 0;   i < 448;  i++)   out[i] =          saved[i] + ac->add_bias;
            ac->dsp.vector_fmul_add_add(out + 448, buf + 448, swindow_prev, saved + 448, ac->add_bias, 128, 1);
            for (i = 576; i < 1024; i++)   out[i] = buf[i] + saved[i] + ac->add_bias;
        }
        if (ics->window_sequence != LONG_START_SEQUENCE) {
            ac->dsp.vector_fmul_reverse(saved, buf + 1024, lwindow, 1024);
        } else {
            memcpy(saved, buf + 1024, 448 * sizeof(float));
            ac->dsp.vector_fmul_reverse(saved + 448, buf + 1024 + 448, swindow, 128);
            memset(saved + 576, 0, 448 * sizeof(float));
        }
    } else {
        if (ics->window_sequence_prev == ONLY_LONG_SEQUENCE || ics->window_sequence_prev == LONG_STOP_SEQUENCE)
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
static void window_ssr_tool(AACContext * ac, SingleChannelElement * sce, float * in, float * out) {
    IndividualChannelStream * ics = &sce->ics;
    const float * lwindow      = ics->use_kb_window[0] ? kbd_long_1024 : sine_long_1024;
    const float * swindow      = ics->use_kb_window[0] ? kbd_short_128 : sine_short_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? kbd_long_1024 : sine_long_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? kbd_short_128 : sine_short_128;
    float * buf = ac->buf_mdct;
    if (ics->window_sequence != EIGHT_SHORT_SEQUENCE) {
        int i;
        ff_imdct_calc(&ac->mdct, buf, in, out);
        if (ics->window_sequence != LONG_STOP_SEQUENCE) {
            vector_fmul_dst(ac, out, buf, lwindow_prev, 256);
        } else {
            memset(out, 0, 112 * sizeof(float));
            vector_fmul_dst(ac, out + 112, buf + 112, swindow_prev, 32);
            memcpy(out + 144, buf + 144, 112 * sizeof(float));
        }
        if (ics->window_sequence != LONG_START_SEQUENCE) {
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

static void ssr_gain_tool(AACContext * ac, SingleChannelElement * sce, int band, float * in, float * preret, float * saved) {
    // TODO: 'in' buffer gain normalization
    if (sce->ics.window_sequence != EIGHT_SHORT_SEQUENCE) {
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

static void ssr_ipqf_tool(AACContext * ac, SingleChannelElement * sce, float * preret) {
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

static void ssr_trans(AACContext * ac, SingleChannelElement * sce) {
    float * in = sce->coeffs;
    DECLARE_ALIGNED_16(float, tmp_buf[512]);
    DECLARE_ALIGNED_16(float, tmp_ret[1024]);
    int b;
    for (b = 0; b < 4; b++) {
        if (b & 1) { // spectral reverse
            vector_reverse(ac, tmp_buf, in + 256 * b, 256);
            memcpy(in + 256 * b, tmp_buf, 256 * sizeof(float));
        }
        window_ssr_tool(ac, sce, in + 256 * b, tmp_buf);
        ssr_gain_tool(ac, sce, b, tmp_buf, tmp_ret + 256 * b, sce->saved + 256 * b);
    }
    ssr_ipqf_tool(ac, sce, tmp_ret);
}
#endif /* AAC_SSR */

static void coupling_dependent_trans(AACContext * ac, ChannelElement * cc, SingleChannelElement * sce, int index) {
    IndividualChannelStream * ics = &cc->ch[0].ics;
    const uint16_t * offsets = ics->swb_offset;
    float * dest = sce->coeffs;
    float * src = cc->ch[0].coeffs;
    int g, i, group, k;
    if(ac->m4ac.object_type == AOT_AAC_LTP) {
        av_log(ac->avccontext, AV_LOG_ERROR,
               "Dependent coupling is not supported together with LTP\n");
        return;
    }
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            if (cc->ch[0].cb[g][i] != ZERO_HCB) {
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

static void coupling_independent_trans(AACContext * ac, ChannelElement * cc, SingleChannelElement * sce, int index) {
    int i;
    float gain = cc->coup.gain[index][0][0] * sce->mixing_gain;
    for (i = 0; i < 1024; i++)
        sce->ret[i] += gain * (cc->ch[0].ret[i] - ac->add_bias);
}

static void transform_coupling_tool(AACContext * ac, ChannelElement * cc,
        void (*cc_trans)(AACContext * ac, ChannelElement * cc, SingleChannelElement * sce, int index))
{
    int c;
    int index = 0;
    ChannelCoupling * coup = &cc->coup;
    for (c = 0; c <= coup->num_coupled; c++) {
        if (     !coup->is_cpe[c] && ac->che[ID_SCE][coup->tag_select[c]]) {
            cc_trans(ac, cc, &ac->che[ID_SCE][coup->tag_select[c]]->ch[0], index++);
        } else if(coup->is_cpe[c] && ac->che[ID_CPE][coup->tag_select[c]]) {
            if (!coup->l[c] && !coup->r[c]) {
                cc_trans(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[0], index);
                cc_trans(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[1], index++);
            }
            if (coup->l[c])
                cc_trans(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[0], index++);
            if (coup->r[c])
                cc_trans(ac, cc, &ac->che[ID_CPE][coup->tag_select[c]]->ch[1], index++);
        } else {
            av_log(ac->avccontext, AV_LOG_ERROR,
                   "coupling target %sE[%d] not available\n",
                   coup->is_cpe[c] ? "CP" : "SC", coup->tag_select[c]);
            break;
        }
    }
}

static void coupling_tool(AACContext * ac, int independent, int domain) {
    int i;
    for (i = 0; i < MAX_TAGID; i++) {
        ChannelElement * cc = ac->che[ID_CCE][i];
        if (cc) {
            if (cc->coup.ind_sw && independent) {
                transform_coupling_tool(ac, cc, coupling_independent_trans);
            } else if (!cc->coup.ind_sw && !independent && (cc->coup.domain == domain)) {
                transform_coupling_tool(ac, cc, coupling_dependent_trans);
            }
        }
    }
}

static void transform_sce_tool(AACContext * ac, void (*sce_trans)(AACContext * ac, SingleChannelElement * sce)) {
    int i, j;
    for (i = 0; i < MAX_TAGID; i++) {
        for(j = 0; j < 4; j++) {
            if(ac->che[j][i]) {
                sce_trans(ac, &ac->che[j][i]->ch[0]);
                if(j == ID_CPE)
                    sce_trans(ac, &ac->che[j][i]->ch[1]);
            }
        }
    }
}

static void spec_to_sample(AACContext * ac) {
    coupling_tool(ac, 0, 0);
#ifdef AAC_LTP
    if (ac->audioObjectType == AOT_AAC_LTP)
        transform_sce_tool(ac, ltp_trans);
#endif /* AAC_LTP */
    transform_sce_tool(ac, tns_trans);
    coupling_tool(ac, 0, 1);
#ifdef AAC_SSR
    if (ac->audioObjectType == AOT_AAC_SSR)
        transform_sce_tool(ac, ssr_trans);
    else
#endif /* AAC_SSR */
        transform_sce_tool(ac, window_trans);
    coupling_tool(ac, 1, 1);
#ifdef AAC_LTP
    if (ac->audioObjectType == AOT_AAC_LTP)
        transform_sce_tool(ac, ltp_update_trans);
#endif /* AAC_LTP */
}


static int output_samples(AVCodecContext * avccontext, uint16_t * data, int * data_size) {
    AACContext * ac = avccontext->priv_data;
    int i, ch;
    float *c, *l, *r, *sl, *sr, *out;

    if (!ac->is_saved) {
        ac->is_saved = 1;
        *data_size = 0;
        return 0;
    }

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
                    *out++ = *l++ + *c   - *sl   - *sr   - ac->add_bias * 3;
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

    } else {
        for(i = 0; i < 1024; i++)
            for(ch = 0; ch < avccontext->channels; ch++) {
                ac->interleaved_output[i * avccontext->channels + ch] = ac->output_data[ch][i];
            }
    }

    i = 1024 * avccontext->channels * sizeof(uint16_t);
    if(*data_size < i)
        return -1;
    *data_size = i;

    ac->dsp.float_to_int16(data, ac->interleaved_output, 1024 * avccontext->channels);
    return 0;
}


static int aac_decode_frame(AVCodecContext * avccontext, void * data, int * data_size, const uint8_t * buf, int buf_size) {
    AACContext * ac = avccontext->priv_data;
    GetBitContext gb;
    int id, err, tag;

    init_get_bits(&gb, buf, buf_size*8);

    // parse
    while ((id = get_bits(&gb, 3)) != ID_END) {
        tag = get_bits(&gb, 4);
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
                } else {
                    err = 1;
                    break;
                }
            }
            err = decode_ics(ac, &gb, 0, 0, &ac->che[ID_SCE][tag]->ch[0]);
            break;

        case ID_CPE:
            if (ac->che[ID_CPE][tag]) {
                err = decode_cpe(ac, &gb, tag);
            } else
                err = -1;
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
            data_stream_element(ac, &gb, tag);
            err = 0;
            break;

        case ID_CCE:
            if (ac->che[ID_CCE][tag]) {
                err = decode_cce(ac, &gb, tag);
            } else
                err = -1;
            break;

        case ID_LFE:
            if (ac->che[ID_LFE][tag]) {
                err = decode_ics(ac, &gb, 0, 0, &ac->che[ID_LFE][tag]->ch[0]);
            } else
                err = -1;
            break;

        default:
            err = -1; /* should not happen, but keeps compiler happy */
            break;
        }

        if(err)
            return -1;
    }

    spec_to_sample(ac);
    output_samples(avccontext, data, data_size);

    return buf_size;
}

static int aac_decode_close(AVCodecContext * avccontext) {
    AACContext * ac = avccontext->priv_data;
    int i, j;

    for (i = 0; i < MAX_TAGID; i++) {
        for(j = 0; j < 4; j++)
            che_freep(&ac->che[j][i]);
    }

    ff_mdct_end(&ac->mdct);
    ff_mdct_end(&ac->mdct_small);
#ifdef AAC_LTP
    if (ac->mdct_ltp) {
        ff_mdct_end(ac->mdct_ltp);
        av_free(ac->mdct_ltp);
    }
#endif /* AAC_LTP */
    av_free(ac->interleaved_output);
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

