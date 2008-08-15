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

#include "aac.h"
#include "aactab.h"
#include "aacdectab.h"
#include "mpeg4audio.h"

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <string.h>

static VLC vlc_scalefactors;
static VLC vlc_spectral[11];


// TODO: Maybe add to dsputil?!
#if defined(AAC_LTP) || defined(AAC_SSR)
static void vector_fmul_dst(DSPContext * dsp, float * dst, const float * src0, const float * src1, int len) {
    memcpy(dst, src0, len * sizeof(float));
    dsp->vector_fmul(dst, src1, len);
}
#endif

#if 0
static void vector_fmul_add_add_add(DSPContext * dsp, float * dst, const float * src0, const float * src1,
        const float * src2, const float * src3, float src4, int len) {
    int i;
    dsp->vector_fmul_add_add(dst, src0, src1, src2, src4, len, 1);
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
 * Configure output channel order based on the current program configuration element.
 *
 * @param   che_pos current channel position configuration
 * @param   new_che_pos New channel position configuration - we only do something if it differs from the current one.
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int output_configure(AACContext *ac, enum ChannelPosition che_pos[4][MAX_ELEM_ID],
        enum ChannelPosition new_che_pos[4][MAX_ELEM_ID]) {
    AVCodecContext *avctx = ac->avccontext;
    int i, type, channels = 0;

    if(!memcmp(che_pos, new_che_pos, 4 * MAX_ELEM_ID * sizeof(new_che_pos[0][0])))
        return 0; /* no change */

    memcpy(che_pos, new_che_pos, 4 * MAX_ELEM_ID * sizeof(new_che_pos[0][0]));

    /* Allocate or free elements depending on if they are in the
     * current program configuration.
     *
     * Set up default 1:1 output mapping.
     *
     * For a 5.1 stream the output order will be:
     *    [ Front Left ] [ Front Right ] [ Center ] [ LFE ] [ Surround Left ] [ Surround Right ]
     */

    for(i = 0; i < MAX_ELEM_ID; i++) {
        for(type = 0; type < 4; type++) {
            if(che_pos[type][i]) {
                if(!ac->che[type][i] && !(ac->che[type][i] = av_mallocz(sizeof(ChannelElement))))
                    return AVERROR(ENOMEM);
                if(type != TYPE_CCE) {
                    ac->output_data[channels++] = ac->che[type][i]->ch[0].ret;
                    if(type == TYPE_CPE) {
                        ac->output_data[channels++] = ac->che[type][i]->ch[1].ret;
                    }
                }
            } else
                che_freep(&ac->che[type][i]);
        }
    }

    avctx->channels = channels;
    return 0;
}

/**
 * Decode an array of 4 bit element IDs, optionally interleaved with a stereo/mono switching bit.
 *
 * @param cpe_map Stereo (Channel Pair Element) map, NULL if stereo bit is not present.
 * @param sce_map mono (Single Channel Element) map
 * @param type speaker type/position for these channels
 */
static void decode_channel_map(enum ChannelPosition *cpe_map,
        enum ChannelPosition *sce_map, enum ChannelPosition type, GetBitContext * gb, int n) {
    while(n--) {
        enum ChannelPosition *map = cpe_map && get_bits1(gb) ? cpe_map : sce_map; // stereo or mono map
        map[get_bits(gb, 4)] = type;
    }
}

/**
 * Decode program configuration element; reference: table 4.2.
 *
 * @param   new_che_pos New channel position configuration - we only do something if it differs from the current one.
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_pce(AACContext * ac, enum ChannelPosition new_che_pos[4][MAX_ELEM_ID],
        GetBitContext * gb) {
    int num_front, num_side, num_back, num_lfe, num_assoc_data, num_cc;

    skip_bits(gb, 2);  // object_type

    ac->m4ac.sampling_index = get_bits(gb, 4);
    if(ac->m4ac.sampling_index > 11) {
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

    if (get_bits1(gb))
        skip_bits(gb, 4); // mono_mixdown_tag
    if (get_bits1(gb))
        skip_bits(gb, 4); // stereo_mixdown_tag

    if (get_bits1(gb))
        skip_bits(gb, 3); // mixdown_coeff_index and pseudo_surround

    decode_channel_map(new_che_pos[TYPE_CPE], new_che_pos[TYPE_SCE], AAC_CHANNEL_FRONT, gb, num_front);
    decode_channel_map(new_che_pos[TYPE_CPE], new_che_pos[TYPE_SCE], AAC_CHANNEL_SIDE,  gb, num_side );
    decode_channel_map(new_che_pos[TYPE_CPE], new_che_pos[TYPE_SCE], AAC_CHANNEL_BACK,  gb, num_back );
    decode_channel_map(NULL,                  new_che_pos[TYPE_LFE], AAC_CHANNEL_LFE,   gb, num_lfe  );

    skip_bits_long(gb, 4 * num_assoc_data);

    decode_channel_map(new_che_pos[TYPE_CCE], new_che_pos[TYPE_CCE], AAC_CHANNEL_CC,    gb, num_cc   );

    align_get_bits(gb);

    /* comment field, first byte is length */
    skip_bits_long(gb, 8 * get_bits(gb, 8));
    return 0;
}

/**
 * Set up channel positions based on a default channel configuration
 * as specified in table 1.17.
 *
 * @param   new_che_pos New channel position configuration - we only do something if it differs from the current one.
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int set_default_channel_config(AACContext *ac, enum ChannelPosition new_che_pos[4][MAX_ELEM_ID],
        int channel_config)
{
    if(channel_config < 1 || channel_config > 7) {
        av_log(ac->avccontext, AV_LOG_ERROR, "invalid default channel configuration (%d)\n",
               channel_config);
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

    if(channel_config != 2)
        new_che_pos[TYPE_SCE][0] = AAC_CHANNEL_FRONT; // front center (or mono)
    if(channel_config > 1)
        new_che_pos[TYPE_CPE][0] = AAC_CHANNEL_FRONT; // L + R (or stereo)
    if(channel_config == 4)
        new_che_pos[TYPE_SCE][1] = AAC_CHANNEL_BACK;  // back center
    if(channel_config > 4)
        new_che_pos[TYPE_CPE][(channel_config == 7) + 1]
                                 = AAC_CHANNEL_BACK;  // back stereo
    if(channel_config > 5)
        new_che_pos[TYPE_LFE][0] = AAC_CHANNEL_LFE;   // LFE
    if(channel_config == 7)
        new_che_pos[TYPE_CPE][1] = AAC_CHANNEL_FRONT; // outer front left + outer front right

    return 0;
}

/**
 * Decode GA "General Audio" specific configuration; reference: table 4.1.
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_ga_specific_config(AACContext * ac, GetBitContext * gb, int channel_config) {
    enum ChannelPosition new_che_pos[4][MAX_ELEM_ID];
    int extension_flag, ret;

    if(get_bits1(gb)) {  // frameLengthFlag
        av_log_missing_feature(ac->avccontext, "960/120 MDCT window is", 1);
        return -1;
    }

    if (get_bits1(gb))       // dependsOnCoreCoder
        skip_bits(gb, 14);   // coreCoderDelay
    extension_flag = get_bits1(gb);

    if(ac->m4ac.object_type == AOT_AAC_SCALABLE ||
       ac->m4ac.object_type == AOT_ER_AAC_SCALABLE)
        skip_bits(gb, 3);     // layerNr

    memset(new_che_pos, 0, 4 * MAX_ELEM_ID * sizeof(new_che_pos[0][0]));
    if (channel_config == 0) {
        skip_bits(gb, 4);  // element_instance_tag
        if((ret = decode_pce(ac, new_che_pos, gb)))
            return ret;
    } else {
        if((ret = set_default_channel_config(ac, new_che_pos, channel_config)))
            return ret;
    }
    if((ret = output_configure(ac, ac->che_pos, new_che_pos)))
        return ret;

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
 * Decode audio specific configuration; reference: table 1.13.
 *
 * @param   data        pointer to AVCodecContext extradata
 * @param   data_size   size of AVCCodecContext extradata
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_audio_specific_config(AACContext * ac, void *data, int data_size) {
    GetBitContext gb;
    int i;

    init_get_bits(&gb, data, data_size * 8);

    if((i = ff_mpeg4audio_get_config(&ac->m4ac, data, data_size)) < 0)
        return -1;
    if(ac->m4ac.sampling_index > 11) {
        av_log(ac->avccontext, AV_LOG_ERROR, "invalid sampling rate index %d\n", ac->m4ac.sampling_index);
        return -1;
    }

    skip_bits_long(&gb, i);

    switch (ac->m4ac.object_type) {
    case AOT_AAC_LC:
#ifdef AAC_SSR
    case AOT_AAC_SSR:
#endif /* AAC_SSR */
#ifdef AAC_LTP
    case AOT_AAC_LTP:
#endif /* AAC_LTP */
        if (decode_ga_specific_config(ac, &gb, ac->m4ac.chan_config))
            return -1;
        break;
    default:
        av_log(ac->avccontext, AV_LOG_ERROR, "Audio object type %s%d is not supported.\n",
               ac->m4ac.sbr == 1? "SBR+" : "", ac->m4ac.object_type);
        return -1;
    }
    return 0;
}

/**
 * linear congruential pseudorandom number generator
 *
 * @param   previous_val    pointer to the current state of the generator
 *
 * @return  Returns a 32-bit pseudorandom integer
 */
static av_always_inline int lcg_random(int previous_val) {
    return previous_val * 1664525 + 1013904223;
}

static av_cold int aac_decode_init(AVCodecContext * avccontext) {
    AACContext * ac = avccontext->priv_data;
    int i;

    ac->avccontext = avccontext;

    if (avccontext->extradata_size <= 0 ||
        decode_audio_specific_config(ac, avccontext->extradata, avccontext->extradata_size))
        return -1;

    avccontext->sample_fmt  = SAMPLE_FMT_S16;
    avccontext->sample_rate = ac->m4ac.sample_rate;
    avccontext->frame_size  = 1024;

    AAC_INIT_VLC_STATIC( 0, 144);
    AAC_INIT_VLC_STATIC( 1, 114);
    AAC_INIT_VLC_STATIC( 2, 188);
    AAC_INIT_VLC_STATIC( 3, 180);
    AAC_INIT_VLC_STATIC( 4, 172);
    AAC_INIT_VLC_STATIC( 5, 140);
    AAC_INIT_VLC_STATIC( 6, 168);
    AAC_INIT_VLC_STATIC( 7, 114);
    AAC_INIT_VLC_STATIC( 8, 262);
    AAC_INIT_VLC_STATIC( 9, 248);
    AAC_INIT_VLC_STATIC(10, 384);

    dsputil_init(&ac->dsp, avccontext);

    ac->random_state = 0x1f2e3d4c;

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
    for (i = 0; i < 316; i++)
        ff_aac_pow2sf_tab[i] = pow(2, (i - 200)/4.);
#endif /* CONFIG_HARDCODED_TABLES */

    INIT_VLC_STATIC(&vlc_scalefactors, 7, sizeof(ff_aac_scalefactor_code)/sizeof(ff_aac_scalefactor_code[0]),
        ff_aac_scalefactor_bits, sizeof(ff_aac_scalefactor_bits[0]), sizeof(ff_aac_scalefactor_bits[0]),
        ff_aac_scalefactor_code, sizeof(ff_aac_scalefactor_code[0]), sizeof(ff_aac_scalefactor_code[0]),
        352);

#ifdef AAC_LTP
    ff_mdct_init(&ac->mdct_ltp, 11, 0);
#endif /* AAC_LTP */
#ifdef AAC_SSR
    if (ac->m4ac.object_type == AOT_AAC_SSR) {
        ff_mdct_init(&ac->mdct, 9, 1);
        ff_mdct_init(&ac->mdct_small, 6, 1);
        // window initialization
        ff_kbd_window_init(ff_aac_kbd_long_1024, 4.0, 256);
        ff_kbd_window_init(ff_aac_kbd_short_128, 6.0, 32);
        ff_sine_window_init(ff_sine_1024, 256);
        ff_sine_window_init(ff_sine_128, 32);
        ssr_context_init(&ac->ssrctx);
    } else {
#endif /* AAC_SSR */
        ff_mdct_init(&ac->mdct, 11, 1);
        ff_mdct_init(&ac->mdct_small, 8, 1);
        // window initialization
        ff_kbd_window_init(ff_aac_kbd_long_1024, 4.0, 1024);
        ff_kbd_window_init(ff_aac_kbd_short_128, 6.0, 128);
        ff_sine_window_init(ff_sine_1024, 1024);
        ff_sine_window_init(ff_sine_128, 128);
#ifdef AAC_SSR
    }
#endif /* AAC_SSR */
    return 0;
}

/**
 * Skip data_stream_element; reference: table 4.10.
 */
static void skip_data_stream_element(GetBitContext * gb) {
    int byte_align = get_bits1(gb);
    int count = get_bits(gb, 8);
    if (count == 255)
        count += get_bits(gb, 8);
    if (byte_align)
        align_get_bits(gb);
    skip_bits_long(gb, 8 * count);
}

#ifdef AAC_LTP
static void decode_ltp(AACContext * ac, LongTermPrediction * ltp, GetBitContext * gb, uint8_t max_sfb) {
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
static int decode_ics_info(AACContext * ac, IndividualChannelStream * ics, GetBitContext * gb, int common_window) {
    if (get_bits1(gb)) {
        av_log(ac->avccontext, AV_LOG_ERROR, "Reserved bit set.\n");
        memset(ics, 0, sizeof(IndividualChannelStream));
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
        for (i = 0; i < 7; i++) {
            if (get_bits1(gb)) {
                ics->group_len[ics->num_window_groups-1]++;
            } else {
                ics->num_window_groups++;
                ics->group_len[ics->num_window_groups-1] = 1;
            }
        }
        ics->num_windows   = 8;
        ics->swb_offset    =      swb_offset_128[ac->m4ac.sampling_index];
        ics->num_swb       =  ff_aac_num_swb_128[ac->m4ac.sampling_index];
        ics->tns_max_bands =   tns_max_bands_128[ac->m4ac.sampling_index];
    } else {
        ics->max_sfb       = get_bits(gb, 6);
        ics->num_windows   = 1;
        ics->swb_offset    =     swb_offset_1024[ac->m4ac.sampling_index];
        ics->num_swb       = ff_aac_num_swb_1024[ac->m4ac.sampling_index];
        ics->tns_max_bands =  tns_max_bands_1024[ac->m4ac.sampling_index];
        if (get_bits1(gb)) {
#ifdef AAC_LTP
            if (ac->m4ac.object_type == AOT_AAC_MAIN) {
                assert(0);
            } else {
                if ((ics->ltp.present = get_bits(gb, 1))) {
                    decode_ltp(ac, &ics->ltp, gb, ics->max_sfb);
                }
                if (common_window) {
                    if ((ics->ltp2.present = get_bits(gb, 1))) {
                        decode_ltp(ac, &ics->ltp2, gb, ics->max_sfb);
                    }
                }
            }
#else /* AAC_LTP */
            av_log_missing_feature(ac->avccontext, "Predictor bit set but LTP is", 1);
            memset(ics, 0, sizeof(IndividualChannelStream));
            return -1;
#endif /* AAC_LTP */
#ifdef AAC_LTP
        } else {
            ics->ltp.present = 0;
            ics->ltp2.present = 0;
#endif /* AAC_LTP */
        }
    }

    if(ics->max_sfb > ics->num_swb) {
        av_log(ac->avccontext, AV_LOG_ERROR,
            "Number of scalefactor bands in group (%d) exceeds limit (%d).\n",
            ics->max_sfb, ics->num_swb);
        memset(ics, 0, sizeof(IndividualChannelStream));
        return -1;
    }

    return 0;
}

/**
 * Decode band types (section_data payload); reference: table 4.46.
 *
 * @param   band_type           array of the used band type
 * @param   band_type_run_end   array of the last scalefactor band of a band type run
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_band_types(AACContext * ac, enum BandType band_type[120],
        int band_type_run_end[120], GetBitContext * gb, IndividualChannelStream * ics) {
    int g, idx = 0;
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
                band_type        [idx]   = sect_band_type;
                band_type_run_end[idx++] = sect_len;
            }
        }
    }
    return 0;
}

/**
 * Decode scalefactors; reference: table 4.47.
 *
 * @param   global_gain         first scalefactor value as scalefactors are differentially coded
 * @param   band_type           array of the used band type
 * @param   band_type_run_end   array of the last scalefactor band of a band type run
 * @param   sf                  array of scalefactors or intensity stereo positions
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_scalefactors(AACContext * ac, float sf[120], GetBitContext * gb,
        unsigned int global_gain, IndividualChannelStream * ics,
        enum BandType band_type[120], int band_type_run_end[120]) {
    const int sf_offset = ac->sf_offset + (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE ? 12 : 0);
    int g, i, idx = 0;
    int offset[3] = { global_gain, global_gain - 90, 100 };
    int noise_flag = 1;
    static const char *sf_str[3] = { "Global gain", "Noise gain", "Intensity stereo position" };
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb;) {
            int run_end = band_type_run_end[idx];
            if (band_type[idx] == ZERO_BT) {
                for(; i < run_end; i++, idx++)
                    sf[idx] = 0.;
            }else if((band_type[idx] == INTENSITY_BT) || (band_type[idx] == INTENSITY_BT2)) {
                for(; i < run_end; i++, idx++) {
                    offset[2] += get_vlc2(gb, vlc_scalefactors.table, 7, 3) - 60;
                    if(offset[2] > 255U) {
                        av_log(ac->avccontext, AV_LOG_ERROR,
                            "%s (%d) out of range.\n", sf_str[2], offset[2]);
                        return -1;
                    }
                    sf[idx]  = ff_aac_pow2sf_tab[-offset[2] + 300];
                }
            }else if(band_type[idx] == NOISE_BT) {
                for(; i < run_end; i++, idx++) {
                    if(noise_flag-- > 0)
                        offset[1] += get_bits(gb, 9) - 256;
                    else
                        offset[1] += get_vlc2(gb, vlc_scalefactors.table, 7, 3) - 60;
                    if(offset[1] > 255U) {
                        av_log(ac->avccontext, AV_LOG_ERROR,
                            "%s (%d) out of range.\n", sf_str[1], offset[1]);
                        return -1;
                    }
                    sf[idx]  = -ff_aac_pow2sf_tab[ offset[1] + sf_offset];
                }
            }else {
                for(; i < run_end; i++, idx++) {
                    offset[0] += get_vlc2(gb, vlc_scalefactors.table, 7, 3) - 60;
                    if(offset[0] > 255U) {
                        av_log(ac->avccontext, AV_LOG_ERROR,
                            "%s (%d) out of range.\n", sf_str[0], offset[0]);
                        return -1;
                    }
                    sf[idx] = -ff_aac_pow2sf_tab[ offset[0] + sf_offset];
                }
            }
        }
    }
    return 0;
}

/**
 * Decode pulse data; reference: table 4.7.
 */
static void decode_pulses(Pulse * pulse, GetBitContext * gb, const uint16_t * swb_offset) {
    int i;
    pulse->num_pulse = get_bits(gb, 2) + 1;
    pulse->pos[0]    = get_bits(gb, 5) + swb_offset[get_bits(gb, 6)];
    pulse->amp[0]    = get_bits(gb, 4);
    for (i = 1; i < pulse->num_pulse; i++) {
        pulse->pos[i] = get_bits(gb, 5) + pulse->pos[i-1];
        pulse->amp[i] = get_bits(gb, 4);
    }
}

/**
 * Decode Temporal Noise Shaping data; reference: table 4.48.
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_tns(AACContext * ac, TemporalNoiseShaping * tns,
        GetBitContext * gb, const IndividualChannelStream * ics) {
    int w, filt, i, coef_len, coef_res, coef_compress;
    const int is8 = ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE;
    const int tns_max_order = is8 ? 7 : ac->m4ac.object_type == AOT_AAC_MAIN ? 20 : 12;
    for (w = 0; w < ics->num_windows; w++) {
        tns->n_filt[w] = get_bits(gb, 2 - is8);

        if (tns->n_filt[w])
            coef_res = get_bits1(gb) + 3;

        for (filt = 0; filt < tns->n_filt[w]; filt++) {
            tns->length[w][filt] = get_bits(gb, 6 - 2*is8);

            if ((tns->order[w][filt] = get_bits(gb, 5 - 2*is8)) > tns_max_order) {
                av_log(ac->avccontext, AV_LOG_ERROR, "TNS filter order %d is greater than maximum %d.",
                       tns->order[w][filt], tns_max_order);
                tns->order[w][filt] = 0;
                return -1;
            } else {
                tns->direction[w][filt] = get_bits1(gb);
                coef_compress = get_bits1(gb);
                coef_len = coef_res - coef_compress;
                tns->tmp2_map[w][filt] = tns_tmp2_map[2*coef_compress + coef_res - 3];

                for (i = 0; i < tns->order[w][filt]; i++)
                    tns->coef[w][filt][i] = get_bits(gb, coef_len);
            }
        }
    }
    return 0;
}

#ifdef AAC_SSR
static int decode_gain_control(SingleChannelElement * sce, GetBitContext * gb) {
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
 *
 * @param   ms_present  Indicates mid/side stereo presence. [0] mask is all 0s;
 *                      [1] mask is decoded from bitstream; [2] mask is all 1s;
 *                      [3] reserved for scalable AAC
 */
static void decode_mid_side_stereo(ChannelElement * cpe, GetBitContext * gb,
        int ms_present) {
    int idx;
    if (ms_present == 1) {
        for (idx = 0; idx < cpe->ch[0].ics.num_window_groups * cpe->ch[0].ics.max_sfb; idx++)
            cpe->ms_mask[idx] = get_bits1(gb);
    } else if (ms_present == 2) {
        memset(cpe->ms_mask, 1, cpe->ch[0].ics.num_window_groups * cpe->ch[0].ics.max_sfb * sizeof(cpe->ms_mask[0]));
    }
}

/**
 * Decode spectral data; reference: table 4.50.
 * Dequantize and scale spectral data; reference: 4.6.3.3.
 *
 * @param   coef            array of dequantized, scaled spectral data
 * @param   sf              array of scalefactors or intensity stereo positions
 * @param   pulse_present   set if pulses are present
 * @param   pulse           pointer to pulse data struct
 * @param   band_type       array of the used band type
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_spectrum_and_dequant(AACContext * ac, float coef[1024], GetBitContext * gb, float sf[120],
        int pulse_present, const Pulse * pulse, const IndividualChannelStream * ics, enum BandType band_type[120]) {
    int i, k, g, idx = 0;
    const int c = 1024/ics->num_windows;
    const uint16_t * offsets = ics->swb_offset;
    float *coef_base = coef;

    for (g = 0; g < ics->num_windows; g++)
        memset(coef + g * 128 + offsets[ics->max_sfb], 0, sizeof(float)*(c - offsets[ics->max_sfb]));

    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++, idx++) {
            const int cur_band_type = band_type[idx];
            const int dim = cur_band_type >= FIRST_PAIR_BT ? 2 : 4;
            const int is_cb_unsigned = IS_CODEBOOK_UNSIGNED(cur_band_type);
            int group;
            if (cur_band_type == ZERO_BT) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    memset(coef + group * 128 + offsets[i], 0, (offsets[i+1] - offsets[i])*sizeof(float));
                }
            }else if (cur_band_type == NOISE_BT) {
                const float scale = sf[idx] / ((offsets[i+1] - offsets[i]) * PNS_MEAN_ENERGY);
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++) {
                        ac->random_state  = lcg_random(ac->random_state);
                        coef[group*128+k] = ac->random_state * scale;
                    }
                }
            }else if (cur_band_type != INTENSITY_BT2 && cur_band_type != INTENSITY_BT) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k += dim) {
                        const int index = get_vlc2(gb, vlc_spectral[cur_band_type - 1].table, 6, 3);
                        const int coef_tmp_idx = (group << 7) + k;
                        const float *vq_ptr;
                        int j;
                        if(index >= ff_aac_spectral_sizes[cur_band_type - 1]) {
                            av_log(ac->avccontext, AV_LOG_ERROR,
                                "Read beyond end of ff_aac_codebook_vectors[%d][]. index %d >= %d\n",
                                cur_band_type - 1, index, ff_aac_spectral_sizes[cur_band_type - 1]);
                            return -1;
                        }
                        vq_ptr = &ff_aac_codebook_vectors[cur_band_type - 1][index * dim];
                        if (is_cb_unsigned) {
                            for (j = 0; j < dim; j++)
                                if (vq_ptr[j])
                                    coef[coef_tmp_idx + j] = 1 - 2*(int)get_bits1(gb);
                        }else {
                            for (j = 0; j < dim; j++)
                                coef[coef_tmp_idx + j] = 1.0f;
                        }
                        if (cur_band_type == ESC_BT) {
                            for (j = 0; j < 2; j++) {
                                if (vq_ptr[j] == 64.0f) {
                                    int n = 4;
                                    /* The total length of escape_sequence must be < 22 bits according
                                       to the specification (i.e. max is 11111111110xxxxxxxxxx). */
                                    while (get_bits1(gb) && n < 15) n++;
                                    if(n == 15) {
                                        av_log(ac->avccontext, AV_LOG_ERROR, "error in spectral data, ESC overflow\n");
                                        return -1;
                                    }
                                    n = (1<<n) + get_bits(gb, n);
                                    coef[coef_tmp_idx + j] *= cbrtf(fabsf(n)) * n;
                                }else
                                    coef[coef_tmp_idx + j] *= vq_ptr[j];
                            }
                        }else
                            for (j = 0; j < dim; j++)
                                coef[coef_tmp_idx + j] *= vq_ptr[j];
                        for (j = 0; j < dim; j++)
                            coef[coef_tmp_idx + j] *= sf[idx];
                    }
                }
            }
        }
        coef += ics->group_len[g]<<7;
    }

    if (pulse_present) {
        for(i = 0; i < pulse->num_pulse; i++){
            float co  = coef_base[ pulse->pos[i] ];
            float ico = co / sqrtf(sqrtf(fabsf(co))) + pulse->amp[i];
            coef_base[ pulse->pos[i] ] = cbrtf(fabsf(ico)) * ico;
        }
    }
    return 0;
}

/**
 * Decode an individual_channel_stream payload; reference: table 4.44.
 *
 * @param   common_window   Channels have independent [0], or shared [1], Individual Channel Stream information.
 * @param   scale_flag      scalable [1] or non-scalable [0] AAC (Unused until scalable AAC is implemented.)
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_ics(AACContext * ac, SingleChannelElement * sce, GetBitContext * gb, int common_window, int scale_flag) {
    Pulse pulse;
    TemporalNoiseShaping * tns = &sce->tns;
    IndividualChannelStream * ics = &sce->ics;
    float * out = sce->coeffs;
    int global_gain, pulse_present = 0;

    /* This assignment is to silence a GCC warning about the variable being used
     * uninitialized when in fact it always is.
     */
    pulse.num_pulse = 0;

    global_gain = get_bits(gb, 8);

    if (!common_window && !scale_flag) {
        if (decode_ics_info(ac, ics, gb, 0) < 0)
            return -1;
    }

    if (decode_band_types(ac, sce->band_type, sce->band_type_run_end, gb, ics) < 0)
        return -1;
    if (decode_scalefactors(ac, sce->sf, gb, global_gain, ics, sce->band_type, sce->band_type_run_end) < 0)
        return -1;

    pulse_present = 0;
    if (!scale_flag) {
        if ((pulse_present = get_bits1(gb))) {
            if (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE) {
                av_log(ac->avccontext, AV_LOG_ERROR, "Pulse tool not allowed in eight short sequence.\n");
                return -1;
            }
            decode_pulses(&pulse, gb, ics->swb_offset);
        }
        if ((tns->present = get_bits1(gb)) && decode_tns(ac, tns, gb, ics))
            return -1;
        if (get_bits1(gb)) {
#ifdef AAC_SSR
            int ret;
            if ((ret = decode_gain_control(sce, gb))) return ret;
#else
            av_log_missing_feature(ac->avccontext, "SSR", 1);
            return -1;
#endif
        }
    }

    if (decode_spectrum_and_dequant(ac, out, gb, sce->sf, pulse_present, &pulse, ics, sce->band_type) < 0)
        return -1;
    return 0;
}

/**
 * Mid/Side stereo decoding; reference: 4.6.8.1.3.
 */
static void apply_mid_side_stereo(ChannelElement * cpe) {
    const IndividualChannelStream * ics = &cpe->ch[0].ics;
    float *ch0 = cpe->ch[0].coeffs;
    float *ch1 = cpe->ch[1].coeffs;
    int g, i, k, group, idx = 0;
    const uint16_t * offsets = ics->swb_offset;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++, idx++) {
            if (cpe->ms_mask[idx] &&
                cpe->ch[0].band_type[idx] < NOISE_BT && cpe->ch[1].band_type[idx] < NOISE_BT) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++) {
                        float tmp = ch0[group*128 + k] - ch1[group*128 + k];
                        ch0[group*128 + k] += ch1[group*128 + k];
                        ch1[group*128 + k] = tmp;
                    }
                }
            }
        }
        ch0 += ics->group_len[g]*128;
        ch1 += ics->group_len[g]*128;
    }
}

/**
 * intensity stereo decoding; reference: 4.6.8.2.3
 *
 * @param   ms_present  Indicates mid/side stereo presence. [0] mask is all 0s;
 *                      [1] mask is decoded from bitstream; [2] mask is all 1s;
 *                      [3] reserved for scalable AAC
 */
static void apply_intensity_stereo(ChannelElement * cpe, int ms_present) {
    const IndividualChannelStream * ics = &cpe->ch[1].ics;
    SingleChannelElement * sce1 = &cpe->ch[1];
    float *coef0 = cpe->ch[0].coeffs, *coef1 = cpe->ch[1].coeffs;
    const uint16_t * offsets = ics->swb_offset;
    int g, group, i, k, idx = 0;
    int c;
    float scale;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb;) {
            if (sce1->band_type[idx] == INTENSITY_BT || sce1->band_type[idx] == INTENSITY_BT2) {
                const int bt_run_end = sce1->band_type_run_end[idx];
                for (; i < bt_run_end; i++, idx++) {
                    c = -1 + 2 * (sce1->band_type[idx] - 14);
                    if (ms_present)
                        c *= 1 - 2 * cpe->ms_mask[idx];
                    scale = c * sce1->sf[idx];
                    for (group = 0; group < ics->group_len[g]; group++)
                        for (k = offsets[i]; k < offsets[i+1]; k++)
                            coef1[group*128 + k] = scale * coef0[group*128 + k];
                }
            } else {
                int bt_run_end = sce1->band_type_run_end[idx];
                idx += bt_run_end - i;
                i    = bt_run_end;
            }
        }
        coef0 += ics->group_len[g]*128;
        coef1 += ics->group_len[g]*128;
    }
}

/**
 * Decode a channel_pair_element; reference: table 4.4.
 *
 * @param   elem_id Identifies the instance of a syntax element.
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_cpe(AACContext * ac, GetBitContext * gb, int elem_id) {
    int i, ret, common_window, ms_present = 0;
    ChannelElement * cpe;

    cpe = ac->che[TYPE_CPE][elem_id];
    common_window = get_bits1(gb);
    if (common_window) {
        if (decode_ics_info(ac, &cpe->ch[0].ics, gb, 1))
            return -1;
        i = cpe->ch[1].ics.use_kb_window[0];
        cpe->ch[1].ics = cpe->ch[0].ics;
        cpe->ch[1].ics.use_kb_window[1] = i;
#ifdef AAC_LTP
        cpe->ch[1].ics.ltp = cpe->ch[0].ics.ltp2;
#endif /* AAC_LTP */
        ms_present = get_bits(gb, 2);
        if(ms_present == 3) {
            av_log(ac->avccontext, AV_LOG_ERROR, "ms_present = 3 is reserved.\n");
            return -1;
        } else if(ms_present)
            decode_mid_side_stereo(cpe, gb, ms_present);
    }
    if ((ret = decode_ics(ac, &cpe->ch[0], gb, common_window, 0)))
        return ret;
    if ((ret = decode_ics(ac, &cpe->ch[1], gb, common_window, 0)))
        return ret;

    if (common_window && ms_present)
        apply_mid_side_stereo(cpe);

    apply_intensity_stereo(cpe, ms_present);
    return 0;
}

/**
 * Decode coupling_channel_element; reference: table 4.8.
 *
 * @param   elem_id Identifies the instance of a syntax element.
 *
 * @return  Returns error status. 0 - OK, !0 - error
 */
static int decode_cce(AACContext * ac, GetBitContext * gb, ChannelElement * che) {
    int num_gain = 0;
    int c, g, sfb, ret, idx = 0;
    int sign;
    float scale;
    SingleChannelElement * sce = &che->ch[0];
    ChannelCoupling * coup     = &che->coup;

    coup->coupling_point = 2*get_bits1(gb);
    coup->num_coupled = get_bits(gb, 3);
    for (c = 0; c <= coup->num_coupled; c++) {
        num_gain++;
        coup->type[c] = get_bits1(gb) ? TYPE_CPE : TYPE_SCE;
        coup->id_select[c] = get_bits(gb, 4);
        if (coup->type[c] == TYPE_CPE) {
            coup->ch_select[c] = get_bits(gb, 2);
            if (coup->ch_select[c] == 3)
                num_gain++;
        } else
            coup->ch_select[c] = 1;
    }
    coup->coupling_point += get_bits1(gb);

    if (coup->coupling_point == 2) {
        av_log(ac->avccontext, AV_LOG_ERROR,
            "Independently switched CCE with 'invalid' domain signalled.\n");
        memset(coup, 0, sizeof(ChannelCoupling));
        return -1;
    }

    sign = get_bits(gb, 1);
    scale = pow(2., pow(2., get_bits(gb, 2) - 3));

    if ((ret = decode_ics(ac, sce, gb, 0, 0)))
        return ret;

    for (c = 0; c < num_gain; c++) {
        int cge = 1;
        int gain = 0;
        float gain_cache = 1.;
        if (c) {
            cge = coup->coupling_point == AFTER_IMDCT ? 1 : get_bits1(gb);
            gain = cge ? get_vlc2(gb, vlc_scalefactors.table, 7, 3) - 60: 0;
            gain_cache = pow(scale, gain);
        }
        for (g = 0; g < sce->ics.num_window_groups; g++)
            for (sfb = 0; sfb < sce->ics.max_sfb; sfb++, idx++)
                if (sce->band_type[idx] != ZERO_BT) {
                    if (!cge) {
                        int t = get_vlc2(gb, vlc_scalefactors.table, 7, 3) - 60;
                        if (t) {
                            int s = 1;
                            if (sign) {
                                s  -= 2 * (t & 0x1);
                                t >>= 1;
                            }
                            gain += t;
                            gain_cache = pow(scale, gain) * s;
                        }
                    }
                    coup->gain[c][idx] = gain_cache;
                }
    }
    return 0;
}

/**
 * Decode Spectral Band Replication extension data; reference: table 4.55.
 *
 * @param   crc flag indicating the presence of CRC checksum
 * @param   cnt length of TYPE_FIL syntactic element in bytes
 *
 * @return  Returns number of bytes consumed from the TYPE_FIL element.
 */
static int decode_sbr_extension(AACContext * ac, GetBitContext * gb, int crc, int cnt) {
    // TODO : sbr_extension implementation
    av_log_missing_feature(ac->avccontext, "SBR", 0);
    skip_bits_long(gb, 8*cnt - 4); // -4 due to reading extension type
    return cnt;
}

/**
 * Parse whether channels are to be excluded from Dynamic Range Compression; reference: table 4.53.
 *
 * @return  Returns number of bytes consumed.
 */
static int decode_drc_channel_exclusions(DynamicRangeControl *che_drc, GetBitContext * gb) {
    int i;
    int num_excl_chan = 0;

    do {
        for (i = 0; i < 7; i++)
            che_drc->exclude_mask[num_excl_chan++] = get_bits1(gb);
    } while (num_excl_chan < MAX_CHANNELS - 7 && get_bits1(gb));

    return num_excl_chan / 7;
}

/**
 * Decode dynamic range information; reference: table 4.52.
 *
 * @param   cnt length of TYPE_FIL syntactic element in bytes
 *
 * @return  Returns number of bytes consumed.
 */
static int decode_dynamic_range(DynamicRangeControl *che_drc, GetBitContext * gb, int cnt) {
    int n = 1;
    int drc_num_bands = 1;
    int i;

    /* pce_tag_present? */
    if(get_bits1(gb)) {
        che_drc->pce_instance_tag  = get_bits(gb, 4);
        skip_bits(gb, 4); // tag_reserved_bits
        n++;
    }

    /* excluded_chns_present? */
    if(get_bits1(gb)) {
        n += decode_drc_channel_exclusions(che_drc, gb);
    }

    /* drc_bands_present? */
    if (get_bits1(gb)) {
        che_drc->band_incr            = get_bits(gb, 4);
        che_drc->interpolation_scheme = get_bits(gb, 4);
        n++;
        drc_num_bands += che_drc->band_incr;
        for (i = 0; i < drc_num_bands; i++) {
            che_drc->band_top[i] = get_bits(gb, 8);
            n++;
        }
    }

    /* prog_ref_level_present? */
    if (get_bits1(gb)) {
        che_drc->prog_ref_level = get_bits(gb, 7);
        skip_bits1(gb); // prog_ref_level_reserved_bits
        n++;
    }

    for (i = 0; i < drc_num_bands; i++) {
        che_drc->dyn_rng_sgn[i] = get_bits1(gb);
        che_drc->dyn_rng_ctl[i] = get_bits(gb, 7);
        n++;
    }

    return n;
}

/**
 * Decode extension data (incomplete); reference: table 4.51.
 *
 * @param   cnt length of TYPE_FIL syntactic element in bytes
 *
 * @return Returns number of bytes consumed
 */
static int decode_extension_payload(AACContext * ac, GetBitContext * gb, int cnt) {
    int crc_flag = 0;
    int res = cnt;
    switch (get_bits(gb, 4)) { // extension type
        case EXT_SBR_DATA_CRC:
            crc_flag++;
        case EXT_SBR_DATA:
            res = decode_sbr_extension(ac, gb, crc_flag, cnt);
            break;
        case EXT_DYNAMIC_RANGE:
            res = decode_dynamic_range(&ac->che_drc, gb, cnt);
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
static void apply_tns(float coef[1024], TemporalNoiseShaping * tns, IndividualChannelStream * ics, int decode) {
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
            order  = tns->order[w][filt];
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
static void windowing_and_mdct_ltp(AACContext * ac, float * out, float * in, IndividualChannelStream * ics) {
    const float * lwindow      = ics->use_kb_window[0] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow      = ics->use_kb_window[0] ? ff_aac_kbd_short_128 : ff_sine_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? ff_aac_kbd_short_128 : ff_sine_128;
    float * buf = ac->buf_mdct;
    assert(ics->window_sequence[0] != EIGHT_SHORT_SEQUENCE);
    if (ics->window_sequence[0] != LONG_STOP_SEQUENCE) {
        vector_fmul_dst(&ac->dsp, buf, in, lwindow_prev, 1024);
    } else {
        memset(buf, 0, 448 * sizeof(float));
        vector_fmul_dst(&ac->dsp, buf + 448, in + 448, swindow_prev, 128);
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
    if (!sce->ltp_state && !(sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t))))
        return AVERROR(ENOMEM);
    if (sce->ics.window_sequence[0] != EIGHT_SHORT_SEQUENCE && ac->is_saved) {
        float x_est[2 * 1024], X_est[2 * 1024];
        for (i = 0; i < 2 * 1024; i++)
            x_est[i] = sce->ltp_state[i + 2 * 1024 - ltp->lag] * ltp->coef;

        windowing_and_mdct_ltp(ac, X_est, x_est, &sce->ics);
        if(sce->tns.present) apply_tns(X_est, &sce->tns, &sce->ics, 0);

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

static int update_ltp(SingleChannelElement * sce, int is_saved) {
    int i;
    if (!sce->ltp_state && !(sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t))))
        return AVERROR(ENOMEM);
    if (is_saved) {
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
    const float * lwindow      = ics->use_kb_window[0] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow      = ics->use_kb_window[0] ? ff_aac_kbd_short_128 : ff_sine_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? ff_aac_kbd_short_128 : ff_sine_128;
    float * buf = ac->buf_mdct;
    int i;

    if (ics->window_sequence[0] == EIGHT_SHORT_SEQUENCE) {
        if (ics->window_sequence[1] == ONLY_LONG_SEQUENCE || ics->window_sequence[1] == LONG_STOP_SEQUENCE)
            av_log(ac->avccontext, AV_LOG_WARNING,
                   "Transition from an ONLY_LONG or LONG_STOP to an EIGHT_SHORT sequence detected. "
                   "If you heard an audible artifact, please submit the sample to the FFmpeg developers.\n");
        for (i = 0; i < 2048; i += 256) {
            ff_imdct_calc(&ac->mdct_small, buf + i, in + i/2);
            ac->dsp.vector_fmul_reverse(ac->revers + i/2, buf + i + 128, swindow, 128);
        }
        for (i = 0; i < 448; i++)   out[i] = saved[i] + ac->add_bias;

        ac->dsp.vector_fmul_add_add(out + 448 + 0*128, buf + 0*128, swindow_prev, saved + 448 ,       ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 1*128, buf + 2*128, swindow,      ac->revers + 0*128, ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 2*128, buf + 4*128, swindow,      ac->revers + 1*128, ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 3*128, buf + 6*128, swindow,      ac->revers + 2*128, ac->add_bias, 128, 1);
        ac->dsp.vector_fmul_add_add(out + 448 + 4*128, buf + 8*128, swindow,      ac->revers + 3*128, ac->add_bias,  64, 1);

#if 0
        vector_fmul_add_add_add(&ac->dsp, out + 448 + 1*128, buf + 2*128, swindow,      saved + 448 + 1*128, ac->revers + 0*128, ac->add_bias, 128);
        vector_fmul_add_add_add(&ac->dsp, out + 448 + 2*128, buf + 4*128, swindow,      saved + 448 + 2*128, ac->revers + 1*128, ac->add_bias, 128);
        vector_fmul_add_add_add(&ac->dsp, out + 448 + 3*128, buf + 6*128, swindow,      saved + 448 + 3*128, ac->revers + 2*128, ac->add_bias, 128);
        vector_fmul_add_add_add(&ac->dsp, out + 448 + 4*128, buf + 8*128, swindow,      saved + 448 + 4*128, ac->revers + 3*128, ac->add_bias, 64);
#endif

        ac->dsp.vector_fmul_add_add(saved,       buf + 1024 + 64,    swindow + 64, ac->revers + 3*128+64,  0, 64, 1);
        ac->dsp.vector_fmul_add_add(saved + 64,  buf + 1024 + 2*128, swindow,      ac->revers + 4*128,     0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 192, buf + 1024 + 4*128, swindow,      ac->revers + 5*128,     0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 320, buf + 1024 + 6*128, swindow,      ac->revers + 6*128,     0, 128, 1);
        memcpy(                     saved + 448, ac->revers + 7*128, 128 * sizeof(float));
        memset(                     saved + 576, 0,                  448 * sizeof(float));
    } else {
        ff_imdct_calc(&ac->mdct, buf, in);
        if (ics->window_sequence[0] == LONG_STOP_SEQUENCE) {
            for (i = 0;   i < 448;  i++)   out[i] =          saved[i] + ac->add_bias;
            ac->dsp.vector_fmul_add_add(out + 448, buf + 448, swindow_prev, saved + 448, ac->add_bias, 128, 1);
            for (i = 576; i < 1024; i++)   out[i] = buf[i] + saved[i] + ac->add_bias;
        } else {
            ac->dsp.vector_fmul_add_add(out, buf, lwindow_prev, saved, ac->add_bias, 1024, 1);
        }
        if (ics->window_sequence[0] == LONG_START_SEQUENCE) {
            memcpy(saved, buf + 1024, 448 * sizeof(float));
            ac->dsp.vector_fmul_reverse(saved + 448, buf + 1024 + 448, swindow, 128);
            memset(saved + 576, 0, 448 * sizeof(float));
        } else {
            ac->dsp.vector_fmul_reverse(saved, buf + 1024, lwindow, 1024);
        }
    }
}

#ifdef AAC_SSR
static void windowing_and_imdct_ssr(AACContext * ac, float * out, float * in, IndividualChannelStream * ics) {
    const float * lwindow      = ics->use_kb_window[0] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow      = ics->use_kb_window[0] ? ff_aac_kbd_short_128 : ff_sine_128;
    const float * lwindow_prev = ics->use_kb_window[1] ? ff_aac_kbd_long_1024 : ff_sine_1024;
    const float * swindow_prev = ics->use_kb_window[1] ? ff_aac_kbd_short_128 : ff_sine_128;
    float * buf = ac->buf_mdct;
    if (ics->window_sequence[0] != EIGHT_SHORT_SEQUENCE) {
        ff_imdct_calc(&ac->mdct, buf, in);
        if (ics->window_sequence[0] != LONG_STOP_SEQUENCE) {
            vector_fmul_dst(&ac->dsp, out, buf, lwindow_prev, 256);
        } else {
            memset(out, 0, 112 * sizeof(float));
            vector_fmul_dst(&ac->dsp, out + 112, buf + 112, swindow_prev, 32);
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
            ff_imdct_calc(&ac->mdct_small, buf, in + i * 32);
            vector_fmul_dst(&ac->dsp, out + 64 * i, buf, (i == 0) ? swindow_prev : swindow, 32);
            ac->dsp.vector_fmul_reverse(out + 64 * i + 32, buf + 32, swindow, 32);
        }
    }
}

static void vector_add_dst(float * dst, const float * src0, const float * src1, int len) {
    int i;
    for (i = 0; i < len; i++)
        dst[i] = src0[i] + src1[i];
}

static void apply_ssr_gains(SingleChannelElement * sce, float * preret, float * saved, float * in, int band) {
    // TODO: 'in' buffer gain normalization
    if (sce->ics.window_sequence[0] != EIGHT_SHORT_SEQUENCE) {
        vector_add_dst(preret, in, saved, 256);
        memcpy(saved, in + 256, 256 * sizeof(float));
    } else {
        memcpy(preret, saved, 112 * sizeof(float));
        preret += 112; saved += 112;
        vector_add_dst(preret       , in            , saved    , 32);
        vector_add_dst(preret + 1*32, in + 0*64 + 32, in + 1*64, 32);
        vector_add_dst(preret + 2*32, in + 1*64 + 32, in + 2*64, 32);
        vector_add_dst(preret + 3*32, in + 2*64 + 32, in + 3*64, 32);
        vector_add_dst(preret + 4*32, in + 3*64 + 32, in + 4*64, 16);

        vector_add_dst(saved            , in + 3*64 + 32 + 16, in + 4*64 + 16, 16);
        vector_add_dst(saved        + 16, in + 4*64 + 32     , in + 5*64     , 32);
        vector_add_dst(saved + 1*32 + 16, in + 5*64 + 32     , in + 6*64     , 32);
        vector_add_dst(saved + 2*32 + 16, in + 6*64 + 32     , in + 7*64     , 32);
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

static void vector_reverse(float * dst, const float * src, int len) {
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
            vector_reverse(tmp_buf, in + 256 * b, 256);
            memcpy(in + 256 * b, tmp_buf, 256 * sizeof(float));
        }
        windowing_and_imdct_ssr(ac, tmp_buf, in + 256 * b, &sce->ics);
        apply_ssr_gains(sce, tmp_ret + 256 * b, sce->saved + 256 * b, tmp_buf, b);
    }
    inverse_polyphase_quadrature_filter(ac, sce, tmp_ret);
}
#endif /* AAC_SSR */

/**
 * Apply dependent channel coupling (applied before IMDCT).
 *
 * @param   index   index into coupling gain array
 */
static void apply_dependent_coupling(AACContext * ac, SingleChannelElement * sce, ChannelElement * cc, int index) {
    IndividualChannelStream * ics = &cc->ch[0].ics;
    const uint16_t * offsets = ics->swb_offset;
    float * dest = sce->coeffs;
    const float * src = cc->ch[0].coeffs;
    int g, i, group, k, idx = 0;
    if(ac->m4ac.object_type == AOT_AAC_LTP) {
        av_log(ac->avccontext, AV_LOG_ERROR,
               "Dependent coupling is not supported together with LTP\n");
        return;
    }
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++, idx++) {
            if (cc->ch[0].band_type[idx] != ZERO_BT) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++) {
                        // XXX dsputil-ize
                        dest[group*128+k] += cc->coup.gain[index][idx] * src[group*128+k];
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
 * @param   index   index into coupling gain array
 */
static void apply_independent_coupling(AACContext * ac, SingleChannelElement * sce, ChannelElement * cc, int index) {
    int i;
    for (i = 0; i < 1024; i++)
        sce->ret[i] += cc->coup.gain[index][0] * (cc->ch[0].ret[i] - ac->add_bias);
}

/**
 * channel coupling transformation interface
 *
 * @param   index   index into coupling gain array
 * @param   apply_coupling_method   pointer to (in)dependent coupling function
 */
static void apply_channel_coupling(AACContext * ac, ChannelElement * cc,
        void (*apply_coupling_method)(AACContext * ac, SingleChannelElement * sce, ChannelElement * cc, int index))
{
    int c;
    int index = 0;
    ChannelCoupling * coup = &cc->coup;
    for (c = 0; c <= coup->num_coupled; c++) {
        if (ac->che[coup->type[c]][coup->id_select[c]]) {
            if (coup->ch_select[c] != 2) {
                apply_coupling_method(ac, &ac->che[coup->type[c]][coup->id_select[c]]->ch[0], cc, index);
                if (coup->ch_select[c] != 0)
                    index++;
            }
            if (coup->ch_select[c] != 1)
                apply_coupling_method(ac, &ac->che[coup->type[c]][coup->id_select[c]]->ch[1], cc, index++);
        } else {
            av_log(ac->avccontext, AV_LOG_ERROR,
                   "coupling target %sE[%d] not available\n",
                   coup->type[c] == TYPE_CPE ? "CP" : "SC", coup->id_select[c]);
            break;
        }
    }
}

/**
 * Convert spectral data to float samples, applying all supported tools as appropriate.
 */
static int spectral_to_sample(AACContext * ac) {
    int i, type;
    for (i = 0; i < MAX_ELEM_ID; i++) {
        for(type = 0; type < 4; type++) {
            ChannelElement *che = ac->che[type][i];
            if(che) {
                if(che->coup.coupling_point == BEFORE_TNS)
                    apply_channel_coupling(ac, che, apply_dependent_coupling);
#ifdef AAC_LTP
                if (ac->m4ac.object_type == AOT_AAC_LTP) {
                    int ret;
                    if(che->ch[0].ics.ltp.present && (ret = apply_ltp(ac, &che->ch[0])))
                        return ret;
                    if(type == TYPE_CPE &&
                       che->ch[1].ics.ltp.present && (ret = apply_ltp(ac, &che->ch[1])))
                        return ret;
                }
#endif /* AAC_LTP */
                if(che->ch[0].tns.present)
                    apply_tns(che->ch[0].coeffs, &che->ch[0].tns, &che->ch[0].ics, 1);
                if(che->ch[1].tns.present)
                    apply_tns(che->ch[1].coeffs, &che->ch[1].tns, &che->ch[1].ics, 1);
                if(che->coup.coupling_point == BETWEEN_TNS_AND_IMDCT)
                    apply_channel_coupling(ac, che, apply_dependent_coupling);
#ifdef AAC_SSR
                if (ac->m4ac.object_type == AOT_AAC_SSR) {
                    apply_ssr(ac, &che->ch[0]);
                    if(type == TYPE_CPE)
                        apply_ssr(ac, &che->ch[1]);
                } else {
#endif /* AAC_SSR */
                    imdct_and_windowing(ac, &che->ch[0]);
                    if(type == TYPE_CPE)
                        imdct_and_windowing(ac, &che->ch[1]);
#ifdef AAC_SSR
                }
#endif /* AAC_SSR */
                if(che->coup.coupling_point == AFTER_IMDCT)
                    apply_channel_coupling(ac, che, apply_independent_coupling);
#ifdef AAC_LTP
                if (ac->m4ac.object_type == AOT_AAC_LTP) {
                    int ret;
                    if((ret = update_ltp(&che->ch[0], ac->is_saved)))
                        return ret;
                    if(type == TYPE_CPE && (ret = update_ltp(&che->ch[1], ac->is_saved)))
                        return ret;
                }
#endif /* AAC_LTP */
            }
        }
    }
    return 0;
}

static int aac_decode_frame(AVCodecContext * avccontext, void * data, int * data_size, const uint8_t * buf, int buf_size) {
    AACContext * ac = avccontext->priv_data;
    GetBitContext gb;
    enum RawDataBlockType elem_type;
    int err, elem_id, data_size_tmp;

    init_get_bits(&gb, buf, buf_size*8);

    // parse
    while ((elem_type = get_bits(&gb, 3)) != TYPE_END) {
        elem_id = get_bits(&gb, 4);
        err = -1;

        if(elem_type == TYPE_SCE && elem_id == 1 &&
                !ac->che[TYPE_SCE][elem_id] && ac->che[TYPE_LFE][0]) {
            /* Some streams incorrectly code 5.1 audio as SCE[0] CPE[0] CPE[1] SCE[1]
               instead of SCE[0] CPE[0] CPE[0] LFE[0]. If we seem to have
               encountered such a stream, transfer the LFE[0] element to SCE[1] */
            ac->che[TYPE_SCE][elem_id] = ac->che[TYPE_LFE][0];
            ac->che[TYPE_LFE][0] = NULL;
        }
        if(elem_type && elem_type < TYPE_DSE) {
            if(!ac->che[elem_type][elem_id])
                return -1;
            if(elem_type != TYPE_CCE)
                ac->che[elem_type][elem_id]->coup.coupling_point = 4;
        }

        switch (elem_type) {

        case TYPE_SCE:
            err = decode_ics(ac, &ac->che[TYPE_SCE][elem_id]->ch[0], &gb, 0, 0);
            break;

        case TYPE_CPE:
            err = decode_cpe(ac, &gb, elem_id);
            break;

        case TYPE_CCE:
            err = decode_cce(ac, &gb, ac->che[TYPE_SCE][elem_id]);
            break;

        case TYPE_LFE:
            err = decode_ics(ac, &ac->che[TYPE_LFE][elem_id]->ch[0], &gb, 0, 0);
            break;

        case TYPE_DSE:
            skip_data_stream_element(&gb);
            err = 0;
            break;

        case TYPE_PCE:
        {
            enum ChannelPosition new_che_pos[4][MAX_ELEM_ID];
            memset(new_che_pos, 0, 4 * MAX_ELEM_ID * sizeof(new_che_pos[0][0]));
            if((err = decode_pce(ac, new_che_pos, &gb)))
                break;
            err = output_configure(ac, ac->che_pos, new_che_pos);
            break;
        }

        case TYPE_FIL:
            if (elem_id == 15)
                elem_id += get_bits(&gb, 8) - 1;
            while (elem_id > 0)
                elem_id -= decode_extension_payload(ac, &gb, elem_id);
            err = 0; /* FIXME */
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

    if (!ac->is_saved) {
        ac->is_saved = 1;
        *data_size = 0;
        return buf_size;
    }

    data_size_tmp = 1024 * avccontext->channels * sizeof(int16_t);
    if(*data_size < data_size_tmp) {
        av_log(avccontext, AV_LOG_ERROR,
               "Output buffer too small (%d) or trying to output too many samples (%d) for this frame.\n",
               *data_size, data_size_tmp);
        return -1;
    }
    *data_size = data_size_tmp;

    ac->dsp.float_to_int16_interleave(data, (const float **)ac->output_data, 1024, avccontext->channels);

    return buf_size;
}

static av_cold int aac_decode_close(AVCodecContext * avccontext) {
    AACContext * ac = avccontext->priv_data;
    int i, type;

    for (i = 0; i < MAX_ELEM_ID; i++) {
        for(type = 0; type < 4; type++)
            che_freep(&ac->che[type][i]);
    }

    ff_mdct_end(&ac->mdct);
    ff_mdct_end(&ac->mdct_small);
#ifdef AAC_LTP
    ff_mdct_end(&ac->mdct_ltp);
#endif /* AAC_LTP */
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
    .sample_fmts = (enum SampleFormat[]){SAMPLE_FMT_S16,SAMPLE_FMT_NONE},
};
