/*
 * AAC decoder
 * Copyright (c) 2005-2006 Oded Shimon ( ods15 ods15 dyndns org )
 * Copyright (c) 2006-2007 Maxim Gavrilov ( maxim.gavrilov gmail com )
 *
 * Kaiser-Bessel Derived Window by Justin Ruggles
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

/** TODO
 *  To reduce the memory required, reorder the ret buffer so all buffers are
 *  after each other. And then try to do the interleaving inplace.
 */

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"
#include "random.h"

#include "aactab.h"

#include <assert.h>

#define MAX_CHANNELS 64
#define MAX_TAGID 16

// Audio Object Types
enum {
    AOT_NULL = 0x0,
    AOT_AAC_MAIN,
    AOT_AAC_LC,
    AOT_AAC_SSR,
    AOT_AAC_LTP,
    AOT_SBR,
    AOT_AAC_SCALABLE,
    AOT_TWINVQ,
    AOT_CELP,
    AOT_HVXC,
    AOT_TTSI = 12,
    AOT_MAINSYNTH,
    AOT_WAVESYNTH,
    AOT_MIDI,
    AOT_SAFX,
    AOT_ER_AAC_LC,
    AOT_ER_AAC_LTP = 19,
    AOT_ER_AAC_SCALABLE,
    AOT_ER_TWINVQ,
    AOT_ER_BSAC,
    AOT_ER_AAC_LD,
    AOT_ER_CELP,
    AOT_ER_HVXC,
    AOT_ER_HILN,
    AOT_ER_PARAM,
    AOT_SSC
};

// IDs for raw_data_block
enum {
    ID_SCE = 0x0,
    ID_CPE,
    ID_CCE,
    ID_LFE,
    ID_DSE,
    ID_PCE,
    ID_FIL,
    ID_END
};

// IDs for extension_payload
enum {
    EXT_FILL = 0x0,
    EXT_FILL_DATA,
    EXT_DATA_ELEMENT,
    EXT_DYNAMIC_RANGE = 0xb,
    EXT_SBR_DATA = 0xd,
    EXT_SBR_DATA_CRC = 0xe
};

// window sequences
enum {
    ONLY_LONG_SEQUENCE = 0,
    LONG_START_SEQUENCE,
    EIGHT_SHORT_SEQUENCE,
    LONG_STOP_SEQUENCE
};

// special codebooks
#define ZERO_HCB 0
#define FIRST_PAIR_HCB 5
#define ESC_HCB 11
#define NOISE_HCB 13
#define INTENSITY_HCB2 14
#define INTENSITY_HCB 15
#define ESC_FLAG 16

//tns
#define TNS_MAX_ORDER 20

//ltp
#define MAX_LTP_LONG_SFB 40

typedef struct {
    int present;
    int generated;

    int num_channels;

    int num_front;
    int front_cpe;
    int front_tag[MAX_TAGID];

    int num_side;
    int side_cpe;
    int side_tag[MAX_TAGID];

    int num_back;
    int back_cpe;
    int back_tag[MAX_TAGID];

    int num_lfe;
    int lfe_tag[MAX_TAGID];

    int num_assoc_data;
    int assoc_data_tag[MAX_TAGID];

    int num_cc;
    int cc_ind_sw;
    int cc_tag[MAX_TAGID];

    int mono_mixdown;
    int stereo_mixdown;
    int matrix_mixdown;
    int pseudo_surround;
} program_config_struct;

enum {
    MIXMODE_DEFAULT = 0,
    MIXMODE_1TO1,
    MIXMODE_2TO1,
    MIXMODE_1TO2,
    MIXMODE_2TO2,
    MIXMODE_MATRIX1,
    MIXMODE_MATRIX2,
    MIXMODE_UNKNOWN
};

typedef struct {
    int mode;
    int c_tag;
    int lr_tag;
    int sur_tag;
    float sce_gain[MAX_TAGID];
    float cpe_gain[MAX_TAGID][2];
    float lfe_gain[MAX_TAGID];
} mix_config_struct;

typedef struct {
    int present;
    int lag;
    float coef;
    int used[MAX_LTP_LONG_SFB];
} ltp_struct;

typedef struct {
    int intensity_present;
    int noise_present;

    int max_sfb;
    int window_sequence;
    int window_shape;
    int window_shape_prev;
    int predictor;
    int num_window_groups;
    uint8_t grouping;
    uint8_t group_len[8];
    // ltp
    ltp_struct ltp;
    ltp_struct ltp2;
    // calculated
    const uint16_t *swb_offset;
    int num_swb;
    int num_windows;
    int tns_max_bands;
} ics_struct;

typedef struct {
    int present;
    int n_filt[8];
    int length[8][4];
    int direction[8][4];
    int order[8][4];
    float *tmp2_map[8];
    int coef[8][4][TNS_MAX_ORDER];
} tns_struct;

typedef struct {
    int present;
    int mask[8][64];
} ms_struct;

// dynamic range compression
typedef struct {
    int pce_instance_tag;
    int drc_tag_reserved_bits;
    int dyn_rng_sgn[17];
    int dyn_rng_ctl[17];
    int exclude_mask[MAX_CHANNELS];
    int additional_excluded_chns[MAX_CHANNELS];
    int drc_band_incr;
    int drc_interpolation_scheme;
    int drc_band_top[17];
    int prog_ref_level;
    int prog_ref_level_reserved_bits;
} drc_struct;

typedef struct {
    int present;
    int num_pulse;
    int start;
    int offset[4];
    int amp[4];
} pulse_struct;

typedef struct {
    int max_band;
    int adjust_num[4][8];
    int alev[4][8][8];
    int aloc[4][8][8];
    float buf[4][24];
} ssr_struct;

typedef struct {
    int ind_sw;
    int domain;

    int num_coupled;
    int is_cpe[9];
    int tag_select[9];
    int l[9];
    int r[9];

    float gain[18][8][64];
} coupling_struct;

// individual channel element
typedef struct {
    int global_gain;
    float mixing_gain;
    ics_struct ics;
    tns_struct tns;
    int cb[8][64];   // codebooks
    float sf[8][64];
    DECLARE_ALIGNED_16(float, coeffs[1024]);
    DECLARE_ALIGNED_16(float, saved[1024]);
    DECLARE_ALIGNED_16(float, ret[1024]);
    int16_t *ltp_state;
    ssr_struct *ssr;
} sce_struct;

// channel element
typedef struct {
    int common_window;
    ms_struct ms;
    sce_struct ch[2];
} cpe_struct;

typedef struct {
    coupling_struct coup;
    sce_struct ch;
} cc_struct;

typedef struct {
    float q[4][4];
    float t0[4][12];
    float t1[4][12];
} ssr_context;

typedef struct {
    // objects
    AVCodecContext * avccontext;
    GetBitContext gb;
    VLC mainvlc;
    VLC books[11];

    // main config
    int audioObjectType;
    int ext_audioObjectType;
    int sbr_present;
    int sampling_index;
    int ext_sampling_index;
    int sample_rate;
    int ext_sample_rate;
    int channels;
    int frame_length;

    // decoder param
    program_config_struct pcs;
    mix_config_struct mix;
    sce_struct * che_sce[MAX_TAGID];
    cpe_struct * che_cpe[MAX_TAGID];
    sce_struct * che_lfe[MAX_TAGID];
    cc_struct * che_cc[MAX_TAGID];
    drc_struct * che_drc;

    DECLARE_ALIGNED_16(float, buf_mdct[2048]);
    int is_saved;

    //cashes
    const uint16_t *swb_offset_1024;
    const uint16_t *swb_offset_128;
    int num_swb_1024;
    int num_swb_128;
    int tns_max_bands_1024;
    int tns_max_bands_128;

    // tables
    DECLARE_ALIGNED_16(float, kbd_long_1024[1024]);
    DECLARE_ALIGNED_16(float, kbd_short_128[128]);
    DECLARE_ALIGNED_16(float, sine_long_1024[1024]);
    DECLARE_ALIGNED_16(float, sine_short_128[128]);
    DECLARE_ALIGNED_16(float, pow2sf_tab[256]);
    DECLARE_ALIGNED_16(float, intensity_tab[256]);
    DECLARE_ALIGNED_16(float, ivquant_tab[256]);
    DECLARE_ALIGNED_16(float, revers[1024]);
    float* interleaved_output;
    float* iop;

    MDCTContext mdct;
    MDCTContext mdct_small;
    MDCTContext *mdct_ltp;
    DSPContext dsp;
    int * vq[11];
    ssr_context * ssrctx;
    AVRandomState random_state;

    //bias values
    float add_bias;
    float scale_bias;

    // statistics
    int num_frame;
} AACContext;


//aux
// TODO: Maybe add to dsputil?!
static void vector_add_dst(AACContext * ac, float * dst, const float * src0, const float * src1, int len) {
    int i;
    for (i = 0; i < len; i++)
        dst[i] = src0[i] + src1[i];
}

static void vector_fmul_dst(AACContext * ac, float * dst, const float * src0, const float * src1, int len) {
    memcpy(dst, src0, len * sizeof(float));
    ac->dsp.vector_fmul(dst, src1, len);
}

static void vector_fmul_add_add_add(AACContext * ac, float * dst, const float * src0, const float * src1, const float * src2, const float * src3, float src4, int len) {
    int i;
    ac->dsp.vector_fmul_add_add(dst, src0, src1, src2, src4, len, 1);
    for (i = 0; i < len; i++)
        dst[i] += src3[i];
}

static void vector_reverse(AACContext * ac, float * dst, const float * src, int len) {
    int i;
    for (i = 0; i < len; i++)
        dst[i] = src[len - i];
}

static inline int16_t LTP_ROUND(float x) {
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


// aux
/**
 * Generate a Kaiser Window.
 */
static void k_window_init(int alpha, float *window, int n, int iter) {
    int j, k;
    float a, x;
    a = alpha * M_PI / n;
    a = a*a;
    for(k=0; k<n; k++) {
        x = k * (n - k) * a;
        window[k] = 1.0;
        for(j=iter; j>0; j--) {
            window[k] = (window[k] * x / (j*j)) + 1.0;
        }
    }
}

/**
 * Generate a Kaiser-Bessel Derived Window.
 * @param alpha  determines window shape
 * @param window array to fill with window values
 * @param n      length of the window
 * @param iter   number of iterations to use in BesselI0
 */
static void kbd_window_init(int alpha, float *window, int n, int iter) {
    int k, n2;
    float *kwindow;

    n2 = n >> 1;
    kwindow = &window[n2];
    k_window_init(alpha, kwindow, n2, iter);
    window[0] = kwindow[0];
    for(k=1; k<n2; k++) {
        window[k] = window[k-1] + kwindow[k];
    }
    for(k=0; k<n2; k++) {
        window[k] = sqrt(window[k] / (window[n2-1]+1));
        //window[n-1-k] = window[k];
    }
}

/**
 * Generate a sine Window.
 */
static void sine_window_init(float *window, int n) {
    const float alpha = M_PI / n;
    int i;
    for(i = 0; i < n/2; i++)
        window[i] = sin((i + 0.5) * alpha);
}

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

// General functions
#define TAG_MASK 0x00f
#define FLAG_SCE 0x100
#define FLAG_CPE 0x200
#define FLAG_LFE 0x400
#define FLAG_CCE 0x800

static int program_config_element_add_channel(AACContext * ac, int flag_tag) {
    program_config_struct * pcs = &ac->pcs;
    if (pcs->present)
        return 0;
    pcs->generated = 1;
    switch (ac->channels) {
        case 8:
        case 7:
            if ((pcs->num_channels == 3) && (FLAG_CPE & flag_tag)) {
                pcs->num_side = 1;
                pcs->side_cpe = 1;
                pcs->side_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 5;
                return 0;
            }
            if ((pcs->num_channels == 5) && (FLAG_CPE & flag_tag)) {
                pcs->num_back = 1;
                pcs->back_cpe = 1;
                pcs->back_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 7;
                return 0;
            }
            if ((pcs->num_channels == 7) && (FLAG_LFE & flag_tag)) {
                pcs->num_lfe = 1;
                pcs->lfe_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 8;
                return 0;
            }
            goto lab3;
        case 6:
            if ((pcs->num_channels == 5) && (FLAG_LFE & flag_tag)) {
                pcs->num_lfe = 1;
                pcs->lfe_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 6;
                return 0;
            }
        case 5:
            if ((pcs->num_channels == 3) && (FLAG_CPE & flag_tag)) {
                pcs->num_back = 1;
                pcs->back_cpe = 1;
                pcs->back_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 5;
                return 0;
            }
            goto lab3;
        case 4:
            if ((pcs->num_channels == 3) && (FLAG_SCE & flag_tag)) {
                pcs->num_back = 1;
                pcs->back_cpe = 0;
                pcs->back_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 4;
                return 0;
            }
lab3:
        case 3:
            if ((pcs->num_channels == 1) && (FLAG_CPE & flag_tag)) {
                pcs->num_front = 2;
                pcs->front_cpe = 2;
                pcs->front_tag[1] = flag_tag & TAG_MASK;
                pcs->num_channels = 3;
                return 0;
            }
        case 1:
            if ((pcs->num_channels == 0) && (FLAG_SCE & flag_tag)) {
                pcs->num_front = 1;
                pcs->front_cpe = 0;
                pcs->front_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 1;
                return 0;
            }
            break;
        case 2:
            if ((pcs->num_channels == 0) && (FLAG_CPE & flag_tag)) {
                pcs->num_front = 1;
                pcs->front_cpe = 1;
                pcs->front_tag[0] = flag_tag & TAG_MASK;
                pcs->num_channels = 2;
                return 0;
            }
            break;
    }
    pcs->generated = 0;
    return 1;
}

static int program_config_element(AACContext * ac, GetBitContext * gb) {
    program_config_struct * pcs = &ac->pcs;
    int id, object_type, i;
    assert(ac->channels == 0);
    pcs->present = 1;
    id = get_bits(gb, 4);
    object_type = get_bits(gb, 2);

    ac->sampling_index = get_bits(gb, 4);
    assert(ac->sampling_index <= 12);
    ac->sample_rate = sampling_table[ac->sampling_index];
    pcs->num_front = get_bits(gb, 4);
    pcs->num_side = get_bits(gb, 4);
    pcs->num_back = get_bits(gb, 4);
    pcs->num_lfe = get_bits(gb, 2);
    pcs->num_assoc_data = get_bits(gb, 3);
    pcs->num_cc = get_bits(gb, 4);

    pcs->mono_mixdown = get_bits1(gb) ? get_bits(gb, 4) + 1: 0;
    pcs->stereo_mixdown = get_bits1(gb) ? get_bits(gb, 4) + 1: 0;
    assert(pcs->mono_mixdown == 0 && pcs->stereo_mixdown == 0);

    if (get_bits1(gb)) {
        pcs->matrix_mixdown = get_bits(gb, 2) + 1;
        pcs->pseudo_surround = get_bits1(gb);
    } else {
        pcs->matrix_mixdown = 0;
        pcs->pseudo_surround = 0;
    }

    pcs->front_cpe = 0;
    ac->channels += pcs->num_front;
    for (i = 0; i < pcs->num_front; i++) {
        if (get_bits1(gb)) {
            pcs->front_cpe |= (1 << i);
            ac->channels++;
        }
        pcs->front_tag[i] = get_bits(gb, 4);
    }
    pcs->side_cpe = 0;
    ac->channels += pcs->num_side;
    for (i = 0; i < pcs->num_side; i++) {
        if (get_bits1(gb)) {
            pcs->side_cpe |= (1 << i);
            ac->channels++;
        }
        pcs->side_tag[i] = get_bits(gb, 4);
    }
    pcs->back_cpe = 0;
    ac->channels += pcs->num_back;
    for (i = 0; i < pcs->num_back; i++) {
        if (get_bits1(gb)) {
            pcs->back_cpe |= (1 << i);
            ac->channels++;
        }
        pcs->back_tag[i] = get_bits(gb, 4);
    }
    ac->channels += pcs->num_lfe;
    for (i = 0; i < pcs->num_lfe; i++)
        pcs->lfe_tag[i] = get_bits(gb, 4);

    pcs->num_channels = ac->channels;
    // not a real audio channel
    for (i = 0; i < pcs->num_assoc_data; i++)
        pcs->assoc_data_tag[i] = get_bits(gb, 4);
    pcs->cc_ind_sw = 0;
    for (i = 0; i < pcs->num_cc; i++) {
        pcs->cc_ind_sw |= (get_bits1(gb) << i);
        pcs->cc_tag[i] = get_bits(gb, 4);
    }
    align_get_bits(gb);
    skip_bits(gb, 8 * get_bits(gb, 8));
    return 0;
}

static int GASpecificConfig(AACContext * ac, GetBitContext * gb) {
    int ext = 0;
    ac->frame_length = get_bits1(gb);
    if (get_bits1(gb))
        get_bits(gb, 14);
    ext = get_bits1(gb);
    assert(ac->frame_length == 0);
    assert(ext == 0);
    if (ac->channels == 0)
        program_config_element(ac, gb);
    if (ext) {
        switch (ac->audioObjectType) {
            case 22:
                get_bits(gb, 5);
                get_bits(gb, 11);
                break;
            case 17:
            case 19:
            case 20:
            case 23:
                get_bits(gb, 3);
                break;
        }
        if (get_bits1(gb)) ;
    }
    return 0;
}

static inline int GetAudioObjectType(GetBitContext * gb) {
    int result = get_bits(gb, 5);
    if (result == 31)
        result = 32 + get_bits(gb, 6);
    return result;
}

static inline int GetSampleRate(GetBitContext * gb, int *index, int *rate) {
    *index = get_bits(gb, 4);
    if (*index == 0xf) {
        *index = -1;
        *rate = get_bits(gb, 24);
    } else {
        assert(*index <= 12);
        *rate = sampling_table[*index];
    }
    return 0;
}

static int AudioSpecificConfig(AACContext * ac, void *data, int data_size) {
    GetBitContext * gb = &ac->gb;

    init_get_bits(gb, data, data_size * 8);

    memset(&ac->pcs, 0, sizeof(ac->pcs));

    ac->audioObjectType = GetAudioObjectType(gb);
    assert(ac->audioObjectType == AOT_AAC_LC || //ac->audioObjectType == AOT_AAC_MAIN ||
            ac->audioObjectType == AOT_AAC_LTP || ac->audioObjectType == AOT_AAC_SSR);
    if (GetSampleRate(gb, &ac->sampling_index, &ac->sample_rate)) return 1;
    ac->channels = get_bits(gb, 4);
    //assert(ac->channels == 2);

    ac->sbr_present = 0;
    if (ac->audioObjectType == AOT_SBR) {
        ac->ext_audioObjectType = ac->audioObjectType;
        ac->sbr_present = 1;
        if (GetSampleRate(gb, &ac->ext_sampling_index, &ac->ext_sample_rate)) return 1;
        ac->audioObjectType = GetAudioObjectType(gb);
    } else {
        ac->ext_audioObjectType = 0;
    }

    switch (ac->audioObjectType) {
        case AOT_AAC_MAIN:
        case AOT_AAC_LC:
        case AOT_AAC_SSR:
        case AOT_AAC_LTP:
        case AOT_AAC_SCALABLE:
        case AOT_TWINVQ:
            if (GASpecificConfig(ac, gb))
                return 1;
            break;
        case AOT_SBR:
            return 1;
        case AOT_CELP:
        case AOT_HVXC:
            assert(0);
            break;
    };
    if ((ac->ext_audioObjectType != 5) && (8 * data_size - get_bits_count(gb) >= 16)) {
        if (get_bits(gb, 11) == 0x2b7) { // syncExtensionType
            ac->ext_audioObjectType = GetAudioObjectType(gb);
            if (ac->ext_audioObjectType == AOT_SBR) {
                ac->sbr_present = get_bits1(gb);
                if (ac->sbr_present) {
                    if (GetSampleRate(gb, &ac->ext_sampling_index, &ac->ext_sample_rate)) return 1;
                }
            }
        }
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
    ac->is_saved = 0;

    if (AudioSpecificConfig(ac, avccontext->extradata, avccontext->extradata_size))
        return 1;

    if (avccontext->channels == 0) {
        avccontext->channels = ac->channels;
    } else if ((avccontext->channels == 1) && ((ac->channels == 7) || (ac->channels == 6) || (ac->channels == 5))) {
        av_log(avccontext, AV_LOG_INFO, "aac: Downmix from %d to mono\n", ac->channels);
        if (!ac->pcs.matrix_mixdown) {
            ac->pcs.matrix_mixdown = 1;
            ac->pcs.pseudo_surround = 0;
        }
    } else if ((avccontext->channels == 2) && ((ac->channels == 7) || (ac->channels == 6) || (ac->channels == 5))){
        av_log(avccontext, AV_LOG_INFO, "aac: Downmix from %d to stereo\n", ac->channels);
        if (!ac->pcs.matrix_mixdown) {
            ac->pcs.matrix_mixdown = 1;
            ac->pcs.pseudo_surround = 0;
        }
    } else if (((avccontext->channels == 2) || (avccontext->channels == 1)) && ((ac->channels == 2) || (ac->channels == 1))){
        // it's ok stereo <-> mono
    } else {
        if (avccontext->channels < ac->channels)
            av_log(avccontext, AV_LOG_INFO, "aac: AAC stream has %d channels but output to %d channels\n",
                    ac->channels, avccontext->channels);
        avccontext->channels = ac->channels;
    }
    avccontext->sample_rate = ac->sample_rate;

    /* Allocate aligned reorder buffer */
    ac->interleaved_output = av_malloc(ac->channels * 1024 * sizeof(float));

    for (i = 0; i < 11; i++) {
        static const int mod_cb[11] = { 3, 3, 3, 3, 9, 9, 8, 8, 13, 13, 17 };
        static const int off_cb[11] = { 1, 1, 0, 0, 4, 4, 0, 0,  0,  0,  0 };

        int a_bits_size = sizeof(tmp[i].a_bits[0]);
        int a_code_size = sizeof(tmp[i].a_code[0]);
        int j, values = tmp[i].s/a_code_size;
        int dim = (i >= 4 ? 2 : 4);
        int mod = mod_cb[i], off = off_cb[i], index = 0;
        int ret;
        ret = init_vlc(&ac->books[i], 6, values,
                tmp[i].a_bits, a_bits_size, a_bits_size,
                tmp[i].a_code, a_code_size, a_code_size,
                0);
        assert(!ret);
        ac->vq[i] = av_malloc(dim * values * sizeof(int));
        if (dim == 2) {
            for (j = 0; j < values * dim; j += dim) {
                index = j/dim;
                ac->vq[i][j  ] = (index / (mod            ) - off); index %= mod;
                ac->vq[i][j+1] = (index                     - off);
            }
        } else {
            for (j = 0; j < values * dim; j += dim) {
                index = j/dim;
                ac->vq[i][j  ] = (index / (mod * mod * mod) - off); index %= mod*mod*mod;
                ac->vq[i][j+1] = (index / (mod * mod      ) - off); index %= mod*mod;
                ac->vq[i][j+2] = (index / (mod            ) - off); index %= mod;
                ac->vq[i][j+3] = (index                     - off);
            }
        }
    }

    /* Initialize RNG dither */
    av_init_random(0x1f2e3d4c, &ac->random_state);

    // 1024  - compensate wrong IMDCT method
    // 32768 - values in AAC build for ready float->int 16 bit audio, using
    // BIAS method instead needs values -1<x<1
    for (i = 0; i < 256; i++)
        ac->intensity_tab[i] = pow(0.5, (i - 100) / 4.);
    for (i = 0; i < sizeof(ac->ivquant_tab)/sizeof(ac->ivquant_tab[0]); i++)
        ac->ivquant_tab[i] = pow(i, 4./3);

    if(ac->dsp.float_to_int16 == ff_float_to_int16_c) {
        ac->add_bias = 385.0f;
        ac->scale_bias = 32768.0f;
    } else {
        ac->add_bias = 0.0f;
        ac->scale_bias = 1.0f;
    }
    for (i = 0; i < 256; i++)
        ac->pow2sf_tab[i] = pow(2, (i - 100)/4.) /1024./ac->scale_bias;


    // general init
    memset(ac->che_sce, 0, sizeof(ac->che_sce));
    memset(ac->che_cpe, 0, sizeof(ac->che_cpe));
    memset(ac->che_lfe, 0, sizeof(ac->che_lfe));
    memset(ac->che_cc, 0, sizeof(ac->che_cc));
    ac->num_frame = -1;

    ac->swb_offset_1024 = swb_offset_1024[ac->sampling_index];
    ac->num_swb_1024 = num_swb_1024[ac->sampling_index];
    ac->tns_max_bands_1024 = tns_max_bands_1024[ac->sampling_index];
    ac->swb_offset_128 = swb_offset_128[ac->sampling_index];
    ac->num_swb_128 = num_swb_128[ac->sampling_index];
    ac->tns_max_bands_128 = tns_max_bands_128[ac->sampling_index];

    init_vlc(&ac->mainvlc, 7, sizeof(code)/sizeof(code[0]),
            bits, sizeof(bits[0]), sizeof(bits[0]),
            code, sizeof(code[0]), sizeof(code[0]),
            0);

    if (ac->audioObjectType == AOT_AAC_SSR) {
        ff_mdct_init(&ac->mdct, 9, 1);
        ff_mdct_init(&ac->mdct_small, 6, 1);
        // windows init
        kbd_window_init(4, ac->kbd_long_1024, 512, 50);
        kbd_window_init(6, ac->kbd_short_128, 64, 50);
        sine_window_init(ac->sine_long_1024, 512);
        sine_window_init(ac->sine_short_128, 64);
        ac->ssrctx = av_malloc(sizeof(ssr_context));
        ssr_context_init(ac->ssrctx);
    } else {
        ff_mdct_init(&ac->mdct, 11, 1);
        ff_mdct_init(&ac->mdct_small, 8, 1);
        // windows init
        kbd_window_init(4, ac->kbd_long_1024, 2048, 50);
        kbd_window_init(6, ac->kbd_short_128, 256, 50);
        sine_window_init(ac->sine_long_1024, 2048);
        sine_window_init(ac->sine_short_128, 256);
        ac->ssrctx = NULL;
    }
    for (i = 0; i < 128; i++) {
        ac->sine_short_128[i] *= 8.;
        ac->kbd_short_128[i] *= 8.;
    }
    ac->mdct_ltp = NULL;
    dsputil_init(&ac->dsp, avccontext);

    return 0;
}

// Parsers implementation
static int data_stream_element(AACContext * ac, GetBitContext * gb) {
    int id, byte_align;
    int count;
    id = get_bits(gb, 4);
    byte_align = get_bits1(gb);
    count = get_bits(gb, 8);
    if (count == 255)
        count += get_bits(gb, 8);
    if (byte_align)
        align_get_bits(gb);
    skip_bits(gb, 8 * count);
    return 0;
}

static void ltp_data(AACContext * ac, GetBitContext * gb, int max_sfb, ltp_struct * ltp) {
    int sfb;
    if (ac->audioObjectType == AOT_ER_AAC_LD) {
        assert(0);
    } else {
        ltp->lag = get_bits(gb, 11);
        ltp->coef = ltp_coef[get_bits(gb, 3)] * (-2./ac->scale_bias/1024.); // wrong mdct method
        for (sfb = 0; sfb < FFMIN(max_sfb, MAX_LTP_LONG_SFB); sfb++)
            ltp->used[sfb] = get_bits1(gb);
    }
}

static void ics_info(AACContext * ac, GetBitContext * gb, int common_window, ics_struct * ics) {
    int reserved;
    reserved = get_bits1(gb);
    assert(reserved == 0);
    ics->window_sequence = get_bits(gb, 2);
    ics->window_shape_prev = ics->window_shape;
    ics->window_shape = get_bits1(gb);
    if (ics->window_shape_prev == -1)
        ics->window_shape_prev = ics->window_shape;
    ics->num_window_groups = 1;
    ics->group_len[0] = 1;
    if (ics->window_sequence == EIGHT_SHORT_SEQUENCE) {
        int i;
        ics->max_sfb = get_bits(gb, 4);
        ics->grouping = get_bits(gb, 7);
        for (i = 0; i < 7; i++) {
            if (ics->grouping & (1<<(6-i))) {
                ics->group_len[ics->num_window_groups-1]++;
            } else {
                ics->num_window_groups++;
                ics->group_len[ics->num_window_groups-1] = 1;
            }
        }
        ics->swb_offset = ac->swb_offset_128;
        ics->num_swb = ac->num_swb_128;
        ics->num_windows = 8;
        ics->tns_max_bands = ac->tns_max_bands_128;
        //av_log(ac->avccontext, AV_LOG_INFO, " %d groups for %d\n", ics->num_window_groups, ics->grouping);
    } else {
        ics->max_sfb = get_bits(gb, 6);
        ics->swb_offset = ac->swb_offset_1024;
        ics->num_swb = ac->num_swb_1024;
        ics->num_windows = 1;
        ics->tns_max_bands = ac->tns_max_bands_1024;
        ics->predictor = get_bits1(gb);
        if (ics->predictor) {
            if (ac->audioObjectType == AOT_AAC_MAIN) {
                assert(0);
            } else {
                if ((ics->ltp.present = get_bits(gb, 1))) {
                    ltp_data(ac, gb, ics->max_sfb, &ics->ltp);
                }
                if (common_window) {
                    if ((ics->ltp2.present = get_bits(gb, 1))) {
                        ltp_data(ac, gb, ics->max_sfb, &ics->ltp2);
                    }
                }
            }
        } else {
            ics->ltp.present = 0;
            ics->ltp2.present = 0;
        }
    }
}

static inline float ivquant(AACContext * ac, int a) {
    static const float sign[2] = { -1., 1. };
    int tmp = (a>>31);
    int abs_a = (a^tmp)-tmp;
    if (abs_a < sizeof(ac->ivquant_tab)/sizeof(ac->ivquant_tab[0]))
        return sign[tmp+1] * ac->ivquant_tab[abs_a];
    else
        return sign[tmp+1] * pow(abs_a, 4./3);
}

static void section_data(AACContext * ac, GetBitContext * gb, ics_struct * ics, int cb[][64]) {
    int g;
    for (g = 0; g < ics->num_window_groups; g++) {
        int bits = (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 3 : 5;
        int k = 0;
        while (k < ics->max_sfb) {
            int sect_len = 0;
            int sect_len_incr = 1;
            int sect_cb = get_bits(gb, 4);
            assert(sect_cb < 12 || sect_cb == INTENSITY_HCB || sect_cb == INTENSITY_HCB2 || sect_cb == NOISE_HCB);
            while ((sect_len_incr = get_bits(gb, bits)) == (1 << bits)-1)
                sect_len += sect_len_incr;
            sect_len += sect_len_incr;
            sect_len += k;
            for (; k < sect_len && k < ics->max_sfb; k++)
                cb[g][k] = sect_cb;
            assert(k == sect_len);
        }
    }
}

static void scale_factor_data(AACContext * ac, GetBitContext * gb, float mix_gain, int global_gain, ics_struct * ics, const int cb[][64], float sf[][64]) {
    int g, i;
    int intensity = 100; // normalization for intensity_tab lookup table
    int noise = global_gain - 90;
    int noise_flag = 1;
    ics->intensity_present = 0;
    ics->noise_present = 0;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            if (cb[g][i] == ZERO_HCB) {
                sf[g][i] = 0.;
            } else if ((cb[g][i] == INTENSITY_HCB) || (cb[g][i] == INTENSITY_HCB2)) {
                ics->intensity_present = 1;
                intensity += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                assert(!(intensity & (~255)));
                sf[g][i] = ac->intensity_tab[intensity];
            } else if (cb[g][i] == NOISE_HCB) {
                ics->noise_present = 1;
                if (noise_flag) {
                    noise_flag = 0;
                    noise += get_bits(gb, 9) - 256;
                } else {
                    noise += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                }
                sf[g][i] = pow(2.0, 0.25 * noise)/1024./ac->scale_bias;
            } else {
                global_gain += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                assert(!(global_gain & (~255)));
                sf[g][i] = ac->pow2sf_tab[global_gain];
            }
            sf[g][i] *= mix_gain;
        }
    }
}

static void pulse_data(AACContext * ac, GetBitContext * gb, pulse_struct * pulse) {
    int i;
    pulse->num_pulse = get_bits(gb, 2);
    pulse->start = get_bits(gb, 6);
    for (i = 0; i <= pulse->num_pulse; i++) {
        pulse->offset[i] = get_bits(gb, 5);
        pulse->amp[i] = get_bits(gb, 4);
    }
}

static void tns_data(AACContext * ac, GetBitContext * gb, const ics_struct * ics, tns_struct * tns) {
    int w, filt, i, coef_len, coef_res = 0, coef_compress;
    for (w = 0; w < ics->num_windows; w++) {
        tns->n_filt[w] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 1 : 2);
        if (tns->n_filt[w])
            coef_res = get_bits1(gb) + 3;
        for (filt = 0; filt < tns->n_filt[w]; filt++) {
            tns->length[w][filt] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 4 : 6);
            if ((tns->order[w][filt] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 3 : 5))) {
                tns->direction[w][filt] = get_bits1(gb);
                assert(coef_res == 3 || coef_res == 4);
                coef_compress = get_bits1(gb);
                coef_len = coef_res - coef_compress;
                tns->tmp2_map[w] = tns_tmp2_map[(coef_compress << 1) + (coef_res - 3)];
                for (i = 0; i < tns->order[w][filt]; i++)
                    tns->coef[w][filt][i] = get_bits(gb, coef_len);
            }
        }
    }
}

static int gain_control_data(AACContext * ac, GetBitContext * gb, sce_struct * sce) {
    // wd_num wd_test aloc_size
    static const int gain_mode[4][3] = {
        {1, 0, 5}, //ONLY_LONG_SEQUENCE = 0,
        {2, 1, 2}, //LONG_START_SEQUENCE,
        {8, 0, 2}, //EIGHT_SHORT_SEQUENCE,
        {2, 1, 5} //LONG_STOP_SEQUENCE
    };
    const int mode = sce->ics.window_sequence;
    int bd, wd, ad;
    ssr_struct * ssr = sce->ssr;
    if (ssr == NULL)
        ssr = sce->ssr = av_mallocz(sizeof(ssr_struct));
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

static int ms_data(AACContext * ac, GetBitContext * gb, cpe_struct * cpe) {
    ms_struct * ms = &cpe->ms;
    ms->present = get_bits(gb, 2);
    if (ms->present == 1) {
        int g, i;
        for (g = 0; g < cpe->ch[0].ics.num_window_groups; g++)
            for (i = 0; i < cpe->ch[0].ics.max_sfb; i++)
                ms->mask[g][i] = get_bits1(gb);// << i;
    } else if (ms->present == 2) {
        int g, i;
        for (g = 0; g < cpe->ch[0].ics.num_window_groups; g++)
            for (i = 0; i < cpe->ch[0].ics.max_sfb; i++)
                ms->mask[g][i] = 1;// = 0xFFFFFFFFFFFFFFFFULL;
    }
    return 0;
}

static void spectral_data(AACContext * ac, GetBitContext * gb, const ics_struct * ics, const int cb[][64], const float sf[][64], int * icoef) {
    static const int unsigned_cb[] = { 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1 };
    int i, k, g;
    const uint16_t * offsets = ics->swb_offset;

    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            const int cur_cb = cb[g][i];
            const int dim = (cur_cb >= FIRST_PAIR_HCB) ? 2 : 4;
            int group;
            if ((cur_cb == INTENSITY_HCB2) || (cur_cb == INTENSITY_HCB)) {
                continue;
            }
            if (cur_cb == NOISE_HCB) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++)
                        icoef[group*128+k] = av_random(&ac->random_state) & 0x0000FFFF;
                }
                continue;
            }
            if (cur_cb == ZERO_HCB) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    memset(icoef + group * 128 + offsets[i], 0, (offsets[i+1] - offsets[i])*sizeof(int));
                }
                continue;
            }
            for (group = 0; group < ics->group_len[g]; group++) {
                for (k = offsets[i]; k < offsets[i+1]; k += dim) {
                    int index = get_vlc2(gb, ac->books[cur_cb - 1].table, 6, 3);
                    int j;
                    int sign[4] = {1,1,1,1};
                    int ptr[4];
                    memcpy(ptr, &ac->vq[cur_cb - 1][index * dim], dim*sizeof(int));
                    //if (index == -1) av_log(ac->avccontext, AV_LOG_INFO, " tried book %d, at pos %d\n", cb[g][i]-1, before);
                    assert(index != -1);
                    if (unsigned_cb[cur_cb - 1]) {
                        for (j = 0; j < dim; j++)
                            if (ptr[j] && get_bits1(gb))
                                sign[j] = -1;
                    }
                    if (cur_cb == 11) {
                        for (j = 0; j < 2; j++) {
                            if (ptr[j] == 16) {
                                int n = 4;
                                while (get_bits1(gb)) n++;
                                ptr[j] = (1<<n) + get_bits(gb, n);
                            }
                        }
                    }
                    for (j = 0; j < dim; j++)
                        icoef[group*128+k+j] = sign[j] * ptr[j];
                    //out[group*128+j+k] = ivquant(ac, icoef[group*128+j+k]) * sf[g][i];
                    //for (j = 0; j < dim; j++) av_log(ac->avccontext, AV_LOG_INFO, " %4d: %5d %10.3lf => %10.3lf\n", j+k, ptr[j]*sign[j], sf[g][i]*1024*32768, out[group*128+j+k]*1024*32768);
                }
                //av_log(ac->avccontext, AV_LOG_INFO, " checking escape %d[%d] %d\n", ptr[j], j, index);
                assert(k == offsets[i+1]);
            }
        }
        icoef += ics->group_len[g]*128;
    }
}

static void pulse_tool(AACContext * ac, const ics_struct * ics, const pulse_struct * pulse, int * icoef) {
    int i, off;
    if (pulse->present) {
        assert(ics->window_sequence != EIGHT_SHORT_SEQUENCE);
        off = ics->swb_offset[pulse->start];
        for (i = 0; i <= pulse->num_pulse; i++) {
            off += pulse->offset[i];
            if (icoef[off] > 0)
                icoef[off] += pulse->amp[i];
            else
                icoef[off] -= pulse->amp[i];
        }
    }
}

// Tools implementation
static void quant_to_spec_tool(AACContext * ac, const ics_struct * ics, const int * icoef, const int cb[][64], const float sf[][64], float * coef) {
    const uint16_t * offsets = ics->swb_offset;
    int g, i, group, k;

    if(ics->window_sequence == EIGHT_SHORT_SEQUENCE)
      for(g = 0; g < 8; g++)
        memset(coef + g * 128 + offsets[ics->max_sfb], 0, sizeof(float)*(128 - offsets[ics->max_sfb]));
    else
      memset(coef + offsets[ics->max_sfb], 0, sizeof(float)*(1024 - offsets[ics->max_sfb]));

    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            if (cb[g][i] == NOISE_HCB) {
                for (group = 0; group < ics->group_len[g]; group++) {
                    float energy = 0;
                    float scale = 1.;// / (float)(offsets[i+1] - offsets[i]);
                    for (k = offsets[i]; k < offsets[i+1]; k++)
                        energy += (float)icoef[group*128+k] * icoef[group*128+k];
                    scale *= sf[g][i] / sqrt(energy);
                    for (k = offsets[i]; k < offsets[i+1]; k++)
                        coef[group*128+k] = icoef[group*128+k] * scale;
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

static int individual_channel_stream(AACContext * ac, GetBitContext * gb, int common_window, int scale_flag, sce_struct * sce) {
    int icoeffs[1024];
    pulse_struct pulse;
    tns_struct * tns = &sce->tns;
    ics_struct * ics = &sce->ics;
    float * out = sce->coeffs;

    //memset(sf, 0, sizeof(sf));
    memset(&pulse, 0, sizeof(pulse));
    sce->global_gain = get_bits(gb, 8);

    if (!common_window && !scale_flag) {
        ics_info(ac, gb, 0, ics);
    }

    //av_log(ac->avccontext, AV_LOG_INFO, " global_gain: %d, groups: %d\n", global_gain, ics->window_sequence);
    section_data(ac, gb, ics, sce->cb);
    scale_factor_data(ac, gb, sce->mixing_gain, sce->global_gain, ics, sce->cb, sce->sf);

    if (!scale_flag) {
        if ((pulse.present = get_bits1(gb)))
            pulse_data(ac, gb, &pulse);
        if ((tns->present = get_bits1(gb)))
            tns_data(ac, gb, ics, tns);
        if (get_bits1(gb))
            if (gain_control_data(ac, gb, sce)) return 1;
    }

    spectral_data(ac, gb, ics, sce->cb, sce->sf, icoeffs);
    pulse_tool(ac, ics, &pulse, icoeffs);
    quant_to_spec_tool(ac, ics, icoeffs, sce->cb, sce->sf, out);
    return 0;
}

static void ms_tool(AACContext * ac, cpe_struct * cpe) {
    const ms_struct * ms = &cpe->ms;
    const ics_struct * ics = &cpe->ch[0].ics;
    float *ch0 = cpe->ch[0].coeffs;
    float *ch1 = cpe->ch[1].coeffs;
    if (ms->present) {
        int g, i, k, gp;
        const uint16_t * offsets = ics->swb_offset;
        for (g = 0; g < ics->num_window_groups; g++) {
            //av_log(ac->avccontext, AV_LOG_INFO, " masking[%d]: ", g);
            for (gp = 0; gp < ics->group_len[g]; gp++) {
                for (i = 0; i < ics->max_sfb; i++) {
                    if (ms->mask[g][i]) {
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
            //av_log(ac->avccontext, AV_LOG_INFO, "\n");
        }
    }
}

static int single_channel_struct(AACContext * ac, GetBitContext * gb) {
    sce_struct * sce;
    int id = get_bits(gb, 4);
    if (ac->che_sce[id] == NULL) {
        ac->che_sce[id] = av_mallocz(sizeof(sce_struct));
        program_config_element_add_channel(ac, FLAG_SCE | id);
    }
    sce = ac->che_sce[id];
    sce->mixing_gain = ac->mix.sce_gain[id];
    if (individual_channel_stream(ac, gb, 0, 0, sce))
        return 1;
    return 0;
}

static void intensity_tool(AACContext * ac, cpe_struct * cpe) {
    const ics_struct * ics = &cpe->ch[1].ics;
    sce_struct * sce0 = &cpe->ch[0];
    sce_struct * sce1 = &cpe->ch[1];
    if (ics->intensity_present) {
        const uint16_t * offsets = ics->swb_offset;
        int g, gp, i, k;
        int c;
        float scale;
        for (g = 0; g < ics->num_window_groups; g++) {
            for (gp = 0; gp < ics->group_len[g]; gp++) {
                for (i = 0; i < ics->max_sfb; i++) {
                    if ((sce1->cb[g][i] == INTENSITY_HCB) || (sce1->cb[g][i] == INTENSITY_HCB2)) {
                        c = (-1 + 2 * (sce1->cb[g][i] - 14));
                        if (cpe->ms.present)
                            c *= (1 - 2 * cpe->ms.mask[g][i]);
                        scale = c * sce1->sf[g][i];
                        for (k = offsets[i] + gp*128; k < offsets[i+1] + gp*128; k++) {
                            sce1->coeffs[k] = scale * sce0->coeffs[k];
                        }
                    }
                }
            }
        }
    }
}

static int channel_pair_element(AACContext * ac, GetBitContext * gb) {
    int i;
    cpe_struct * cpe;
    int id = get_bits(gb, 4);
    if (ac->che_cpe[id] == NULL) {
        ac->che_cpe[id] = av_mallocz(sizeof(cpe_struct));
        program_config_element_add_channel(ac, FLAG_CPE | id);
    }
    cpe = ac->che_cpe[id];
    cpe->ch[0].mixing_gain = ac->mix.cpe_gain[id][0];
    cpe->ch[1].mixing_gain = ac->mix.cpe_gain[id][1];
    cpe->common_window = get_bits1(gb);
    if (cpe->common_window) {
        ics_info(ac, gb, 1, &cpe->ch[0].ics);
        i = cpe->ch[1].ics.window_shape_prev;
        cpe->ch[1].ics = cpe->ch[0].ics;
        cpe->ch[1].ics.window_shape_prev = i;
        cpe->ch[1].ics.ltp = cpe->ch[0].ics.ltp2;
        ms_data(ac, gb, cpe);
    } else {
        cpe->ms.present = 0;
    }
    if (individual_channel_stream(ac, gb, cpe->common_window, 0, &cpe->ch[0]))
        return 1;
    if (individual_channel_stream(ac, gb, cpe->common_window, 0, &cpe->ch[1]))
        return 1;

    // M/S tool
    if (cpe->common_window) {
        ms_tool(ac, cpe);
    }
    intensity_tool(ac, cpe);
    return 0;
}

static int coupling_channel_element(AACContext * ac, GetBitContext * gb) {
    float cc_scale[] = {
        pow(2, 1/8.), pow(2, 1/4.), pow(2, 0.5), 2.
    };
    int num_gain = 0;
    int c, g, sfb;
    int sign;
    float scale;
    sce_struct * sce;
    coupling_struct * coup;
    int id = get_bits(gb, 4);
    if (ac->che_cc[id] == NULL) {
        ac->che_cc[id] = av_mallocz(sizeof(cc_struct));
        program_config_element_add_channel(ac, FLAG_CCE | id);
    }
    sce = &ac->che_cc[id]->ch;
    sce->mixing_gain = 1.0;

    coup = &ac->che_cc[id]->coup;

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
    scale = cc_scale[get_bits(gb, 2)];

    if (individual_channel_stream(ac, gb, 0, 0, sce))
        return 1;

    for (c = 0; c < num_gain; c++) {
        int cge = 1;
        int gain = 0;
        float gain_cash = 1.;
        if (c != 0) {
            cge = coup->ind_sw ? 1 : get_bits1(gb);
            gain = cge ? get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60: 0;
            gain_cash = pow(scale, (float)gain);
        }
        for (g = 0; g < sce->ics.num_window_groups; g++)
            for (sfb = 0; sfb < sce->ics.max_sfb; sfb++)
                if (sce->cb[g][sfb] == ZERO_HCB) {
                    coup->gain[c][g][sfb] = 0;
                } else {
                    if (cge) {
                        coup->gain[c][g][sfb] = gain_cash;
                    } else {
                        if (sign) {
                            int t = get_vlc2(gb, ac->mainvlc.table, 7, 3);
                            int s = 1 - 2 * (t & 0x1);
                            gain += (t >> 1) - 30;
                            coup->gain[c][g][sfb] = s * pow(scale, (float)gain);
                        } else {
                            gain += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                            coup->gain[c][g][sfb] = pow(scale, (float)gain);
                        }
                    }
                }
    }
    return 0;
}

static int lfe_channel_struct(AACContext * ac, GetBitContext * gb) {
    sce_struct * sce;
    int id = get_bits(gb, 4);
    if (ac->che_lfe[id] == NULL) {
        ac->che_lfe[id] = av_mallocz(sizeof(sce_struct));
        program_config_element_add_channel(ac, FLAG_LFE | id);
    }
    sce = ac->che_lfe[id];
    sce->mixing_gain = ac->mix.lfe_gain[id];
    if (individual_channel_stream(ac, gb, 0, 0, sce))
        return 1;
    return 0;
}

static int sbr_extension_data(AACContext * ac, GetBitContext * gb, int crc, int cnt) {
    // TODO : sbr_extension implementation
    av_log(ac->avccontext, AV_LOG_DEBUG, "aac: SBR not yet supported.\n");
    skip_bits(gb, 8*cnt - 4);
    return cnt;
}


static int excluded_channels(AACContext * ac, GetBitContext * gb) {
    int i;
    int n = 0;
    int num_excl_chan = 7;

    for (i = 0; i < 7; i++)
         ac->che_drc->exclude_mask[i] = get_bits1(gb);
    n++;

    while (get_bits1(gb)) {
        ac->che_drc->additional_excluded_chns[n-1]=1;
        for (i = num_excl_chan; i < num_excl_chan+7; i++)
            ac->che_drc->exclude_mask[i] = get_bits1(gb);
        n++;
        num_excl_chan += 7;
    }
    return n;
}



static int dynamic_range_info(AACContext * ac, GetBitContext * gb, int cnt) {
    int n = 1;
    int drc_num_bands = 1;
    int i;

    if (ac->che_drc == NULL)
        ac->che_drc = av_mallocz(sizeof(drc_struct));

    /* pce_tag_present? */
    if(get_bits1(gb)) {
        ac->che_drc->pce_instance_tag = get_bits(gb, 4);
        ac->che_drc->drc_tag_reserved_bits = get_bits(gb, 4);
        n++;
    }

    /* excluded_chns_present? */
    if(get_bits1(gb)) {
        n += excluded_channels(ac, gb);
    }

    /* drc_bands_present? */
    if (get_bits1(gb)) {
        ac->che_drc->drc_band_incr = get_bits(gb, 4);
        ac->che_drc->drc_interpolation_scheme = get_bits(gb, 4);
        n++;
        drc_num_bands += ac->che_drc->drc_band_incr;
        for (i = 0; i < drc_num_bands; i++) {
            ac->che_drc->drc_band_top[i] = get_bits(gb, 8);
            n++;
        }
    }

    /* prog_ref_level_present? */
    if (get_bits1(gb)) {
        ac->che_drc->prog_ref_level = get_bits(gb, 7);
        ac->che_drc->prog_ref_level_reserved_bits = get_bits1(gb);
        n++;
    }

    for (i = 0; i < drc_num_bands; i++) {
        ac->che_drc->dyn_rng_sgn[i] = get_bits1(gb);
        ac->che_drc->dyn_rng_ctl[i] = get_bits(gb, 7);
        n++;
    }

    return n;
}

/** Parse extension data (incomplete)
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
            skip_bits(gb, 8*cnt - 4);
            break;
    };
    return res;
}

static void tns_filter_tool(AACContext * ac, int decode, sce_struct * sce, float * coef) {
    const ics_struct * ics = &sce->ics;
    const tns_struct * tns = &sce->tns;
    const int mmm = FFMIN(ics->tns_max_bands,  ics->max_sfb);
    int w, filt, m, i, ib;
    int bottom, top, order, start, end, size, inc;
    float tmp;
    float lpc[TNS_MAX_ORDER + 1], b[2 * TNS_MAX_ORDER];
    if (!tns->present) return;
    //av_log(ac->avccontext, AV_LOG_INFO, "%d ", ac->num_frame);
    for (w = 0; w < ics->num_windows; w++) {
        bottom = ics->num_swb;
        for (filt = 0; filt < tns->n_filt[w]; filt++) {
            top = bottom;
            bottom = top - tns->length[w][filt];
            if (bottom < 0)
                bottom = 0;
            order = FFMIN(tns->order[w][filt], TNS_MAX_ORDER);
            if (order == 0)
                continue;

            // tns_decode_coef
            lpc[0] = 1;
            for (m = 1; m <= order; m++) {
                lpc[m] = tns->tmp2_map[w][tns->coef[w][filt][m - 1]];
                for (i = 1; i < m; i++)
                    b[i] = lpc[i] + lpc[m] * lpc[m-i];
                for (i = 1; i < m; i++)
                    lpc[i] = b[i];
            }

            start = ics->swb_offset[FFMIN(bottom, mmm)];
            end = ics->swb_offset[FFMIN(top, mmm)];
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
            if (decode) {
                for (m = 0; m < size; m++) {
                    tmp = coef[start];
                    for (i = 0; i < order; i++)
                        tmp -= b[ib + i] * lpc[i + 1];
                    if (--ib < 0)
                        ib = order - 1;
                    b[ib] = b[ib + order] = tmp;
                    coef[start] = tmp;
                    start += inc;
                }
            } else { // encode
                for (m = 0; m < size; m++) {
                    tmp = coef[start];
                    for (i = 0; i < order; i++)
                        tmp += b[i] * lpc[i + 1];
                    if (--ib < 0)
                        ib = order - 1;
                    b[ib] = b[ib + order] = tmp;
                    coef[start] = tmp;
                    start += inc;
                }
            }
        }
    }
}

static void tns_trans(AACContext * ac, sce_struct * sce) {
    tns_filter_tool(ac, 1, sce, sce->coeffs);
}

static void window_ltp_tool(AACContext * ac, sce_struct * sce, float * in, float * out) {
    ics_struct * ics = &sce->ics;
    const float * lwindow = (ics->window_shape) ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float * swindow = (ics->window_shape) ? ac->kbd_short_128 : ac->sine_short_128;
    const float * lwindow_prev = (ics->window_shape_prev) ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float * swindow_prev = (ics->window_shape_prev) ? ac->kbd_short_128 : ac->sine_short_128;
    float * buf = ac->buf_mdct;
    int i;
    assert(ics->window_sequence != EIGHT_SHORT_SEQUENCE);
    if (ac->mdct_ltp == NULL) {
        ac->mdct_ltp = av_malloc(sizeof(MDCTContext));
        ff_mdct_init(ac->mdct_ltp, 11, 0);
    }
    if (ics->window_sequence != LONG_STOP_SEQUENCE) {
        vector_fmul_dst(ac, buf, in, lwindow_prev, 1024);
    } else {
        memset(buf, 0, 448 * sizeof(float));
        for (i = 448; i < 576; i++) in[i] *= 0.125; // normalize
        vector_fmul_dst(ac, buf + 448, in + 448, swindow_prev, 128);
        memcpy(buf + 576, in + 576, 448 * sizeof(float));
    }
    if (ics->window_sequence != LONG_START_SEQUENCE) {
        ac->dsp.vector_fmul_reverse(buf + 1024, in + 1024, lwindow, 1024);
    } else {
        memcpy(buf + 1024, in + 1024, 448 * sizeof(float));
        for (i = 448; i < 576; i++) in[i + 1024] *= 0.125; // normalize
        ac->dsp.vector_fmul_reverse(buf + 1024 + 448, in + 1024 + 448, swindow, 128);
        memset(buf + 1024 + 576, 0, 448 * sizeof(float));
    }
     ff_mdct_calc(ac->mdct_ltp, out, buf, in); // using in as buffer for mdct
}

static void ltp_trans(AACContext * ac, sce_struct * sce) {
    const ltp_struct * ltp = &sce->ics.ltp;
    const uint16_t * offsets = sce->ics.swb_offset;
    int i, sfb;
    if (!ltp->present)
        return;
    if (sce->ltp_state == NULL)
        sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t));
    if ((sce->ics.window_sequence != EIGHT_SHORT_SEQUENCE) && (ac->is_saved)) {
        float x_est[2 * 1024], X_est[2 * 1024];
        for (i = 0; i < 2 * 1024; i++)
            x_est[i] = (float)sce->ltp_state[i + 2 * 1024 - ltp->lag] * ltp->coef;

        window_ltp_tool(ac, sce, x_est, X_est);
        tns_filter_tool(ac, 0, sce, X_est);

        for (sfb = 0; sfb < FFMIN(sce->ics.max_sfb, MAX_LTP_LONG_SFB); sfb++)
            if (ltp->used[sfb])
                for (i = offsets[sfb]; i < offsets[sfb + 1]; i++)
                    sce->coeffs[i] += X_est[i];
    }
}

static void ltp_update_trans(AACContext * ac, sce_struct * sce) {
    int i;
    if (sce->ltp_state == NULL)
        sce->ltp_state = av_mallocz(4 * 1024 * sizeof(int16_t));
    if (ac->is_saved) {
        for (i = 0; i < 1024; i++) {
            sce->ltp_state[i] = sce->ltp_state[i + 1024];
            sce->ltp_state[i + 1024] = LTP_ROUND(sce->ret[i] - ac->add_bias);
            sce->ltp_state[i + 2 * 1024] = LTP_ROUND(sce->saved[i]);
            //sce->ltp_state[i + 3 * 1024] = 0;
        }
    }
}

static void window_trans(AACContext * ac, sce_struct * sce) {
    ics_struct * ics = &sce->ics;
    float * in = sce->coeffs;
    float * out = sce->ret;
    float * saved = sce->saved;
    const float * lwindow = (ics->window_shape) ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float * swindow = (ics->window_shape) ? ac->kbd_short_128 : ac->sine_short_128;
    const float * lwindow_prev = (ics->window_shape_prev) ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float * swindow_prev = (ics->window_shape_prev) ? ac->kbd_short_128 : ac->sine_short_128;
    float * buf = ac->buf_mdct;
    if (ics->window_sequence != EIGHT_SHORT_SEQUENCE) {
        int i;
        ff_imdct_calc(&ac->mdct, buf, in, out); // out can be abused for now as a temp buffer
        if (ac->is_saved) {
            if (ics->window_sequence != LONG_STOP_SEQUENCE) {
                ac->dsp.vector_fmul_add_add(out, buf, lwindow_prev, saved, ac->add_bias, 1024, 1);
            } else {
                for (i = 0; i < 448; i++) out[i] = saved[i] + ac->add_bias;
                for (i = 448; i < 576; i++) buf[i] *= 0.125; // normalize
                ac->dsp.vector_fmul_add_add(out + 448, buf + 448, swindow_prev, saved + 448, ac->add_bias, 128, 1);
                for (i = 576; i < 1024; i++)   out[i] = buf[i] + saved[i] + ac->add_bias;
            }
        }
        if (ics->window_sequence != LONG_START_SEQUENCE) {
            ac->dsp.vector_fmul_reverse(saved, buf + 1024, lwindow, 1024);
        } else {
            memcpy(saved, buf + 1024, 448 * sizeof(float));
            for (i = 448; i < 576; i++) buf[i + 1024] *= 0.125; // normalize
            ac->dsp.vector_fmul_reverse(saved + 448, buf + 1024 + 448, swindow, 128);
            for (i = 576; i < 1024; i++)   saved[i] = 0.0;
        }
    } else {
        int i;

        for (i = 0; i < 2048; i += 256) {
            ff_imdct_calc(&ac->mdct_small, buf + i, in + i/2, out);
            ac->dsp.vector_fmul_reverse(ac->revers + i/2, buf + i + 128, swindow, 128);
        }
        for (i = 0; i < 448; i++)   out[i] = saved[i] + ac->add_bias;

        ac->dsp.vector_fmul_add_add(out + 448 + 0*128, buf + 0*128, swindow_prev, saved + 448 ,                            ac->add_bias, 128, 1);
        vector_fmul_add_add_add(ac, out + 448 + 1*128, buf + 2*128, swindow,      saved + 448 + 1*128, ac->revers + 0*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 448 + 2*128, buf + 4*128, swindow,      saved + 448 + 2*128, ac->revers + 1*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 448 + 3*128, buf + 6*128, swindow,      saved + 448 + 3*128, ac->revers + 2*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 448 + 4*128, buf + 8*128, swindow,      saved + 448 + 4*128, ac->revers + 3*128, ac->add_bias, 64);

        ac->dsp.vector_fmul_add_add(saved,       buf + 1024 + 64,    swindow + 64, ac->revers + 3*128+64,  0, 64, 1);
        ac->dsp.vector_fmul_add_add(saved + 64,  buf + 1024 + 2*128, swindow,      ac->revers + 4*128,     0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 192, buf + 1024 + 4*128, swindow,      ac->revers + 5*128,     0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 320, buf + 1024 + 6*128, swindow,      ac->revers + 6*128,     0, 128, 1);
        memcpy(                     saved + 448, ac->revers + 7*128, 128 * sizeof(float));
        for (i = 576; i < 1024; i++) saved[i] = 0.0;
    }
}

static void window_ssr_tool(AACContext * ac, sce_struct * sce, float * in, float * out) {
    ics_struct * ics = &sce->ics;
    const float * lwindow = (ics->window_shape) ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float * swindow = (ics->window_shape) ? ac->kbd_short_128 : ac->sine_short_128;
    const float * lwindow_prev = (ics->window_shape_prev) ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float * swindow_prev = (ics->window_shape_prev) ? ac->kbd_short_128 : ac->sine_short_128;
    float * buf = ac->buf_mdct;
    if (ics->window_sequence != EIGHT_SHORT_SEQUENCE) {
        int i;
        ff_imdct_calc(&ac->mdct, buf, in, out);
        if (ics->window_sequence != LONG_STOP_SEQUENCE) {
            vector_fmul_dst(ac, out, buf, lwindow_prev, 256);
        } else {
            memset(out, 0, 112 * sizeof(float));
            for (i = 112; i < 144; i++) buf[i] *= 0.125; // normalize
            vector_fmul_dst(ac, out + 112, buf + 112, swindow_prev, 32);
            memcpy(out + 144, buf + 144, 112 * sizeof(float));
        }
        if (ics->window_sequence != LONG_START_SEQUENCE) {
            ac->dsp.vector_fmul_reverse(out + 256, buf + 256, lwindow, 256);
        } else {
            memcpy(out + 256, buf + 256, 112 * sizeof(float));
            for (i = 112; i < 144; i++) buf[i + 256] *= 0.125; // normalize
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

static void ssr_gain_tool(AACContext * ac, sce_struct * sce, int band, float * in, float * preret, float * saved) {
    // TODO: 'in' buffer gain normalization
    if (sce->ics.window_sequence != EIGHT_SHORT_SEQUENCE) {
        vector_add_dst(ac, preret, in, saved, 256);
        memcpy(saved, in + 256, 256 * sizeof(float));
    } else {
        memcpy(preret, saved, 112 * sizeof(float));
        preret += 112; saved += 112;
        vector_add_dst(ac, preret, in, saved, 32);
        vector_add_dst(ac, preret + 1*32, in + 0*64 + 32, in + 1*64, 32);
        vector_add_dst(ac, preret + 2*32, in + 1*64 + 32, in + 2*64, 32);
        vector_add_dst(ac, preret + 3*32, in + 2*64 + 32, in + 3*64, 32);
        vector_add_dst(ac, preret + 4*32, in + 3*64 + 32, in + 4*64, 16);

        vector_add_dst(ac, saved, in + 3*64 + 32 + 16, in + 4*64 + 16, 16);
        vector_add_dst(ac, saved + 16, in + 4*64 + 32, in + 5*64, 32);
        vector_add_dst(ac, saved + 1*32 + 16, in + 5*64 + 32, in + 6*64, 32);
        vector_add_dst(ac, saved + 2*32 + 16, in + 6*64 + 32, in + 7*64, 32);
        memcpy(saved + 3*32 + 16, in + 7*64 + 32, 32 * sizeof(float));
        memset(saved + 144, 0, 112 * sizeof(float));
    }
}

static void ssr_ipqf_tool(AACContext * ac, sce_struct * sce, float * preret) {
    ssr_context * ctx = ac->ssrctx;
    ssr_struct * ssr = sce->ssr;
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

static void ssr_trans(AACContext * ac, sce_struct * sce) {
    ics_struct * ics = &sce->ics;
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

static void coupling_dependent_trans(AACContext * ac, cc_struct * cc, sce_struct * sce, int index) {
    ics_struct * ics = &cc->ch.ics;
    const uint16_t * offsets = ics->swb_offset;
    float * dest = sce->coeffs;
    float * src = cc->ch.coeffs;
    int g, i, group, k;
    assert(ac->audioObjectType != AOT_AAC_LTP);
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            if (cc->ch.cb[g][i] != ZERO_HCB) {
                float gain = cc->coup.gain[index][g][i] * sce->mixing_gain;
                for (group = 0; group < ics->group_len[g]; group++) {
                    for (k = offsets[i]; k < offsets[i+1]; k++) {
                        dest[group*128+k] += gain * src[group*128+k];
                    }
                }
            }
        }
        dest += ics->group_len[g]*128;
        src += ics->group_len[g]*128;
    }
}

static void coupling_independent_trans(AACContext * ac, cc_struct * cc, sce_struct * sce, int index) {
    int i;
    float gain = cc->coup.gain[index][0][0] * sce->mixing_gain;
    for (i = 0; i < 1024; i++)
        sce->ret[i] += gain * (cc->ch.ret[i] - ac->add_bias);
}

static void transform_coupling_tool(AACContext * ac, cc_struct * cc,
        void (*cc_trans)(AACContext * ac, cc_struct * cc, sce_struct * sce, int index))
{
    int c;
    int index = 0;
    coupling_struct * coup = &cc->coup;
    for (c = 0; c <= coup->num_coupled; c++) {
        if (!coup->is_cpe[c]) {
            assert(ac->che_sce[coup->tag_select[c]] != NULL);
            cc_trans(ac, cc, ac->che_sce[coup->tag_select[c]], index++);
        } else {
            assert(ac->che_cpe[coup->tag_select[c]] != NULL);
            if (!coup->l[c] && !coup->r[c]) {
                cc_trans(ac, cc, &ac->che_cpe[coup->tag_select[c]]->ch[0], index);
                cc_trans(ac, cc, &ac->che_cpe[coup->tag_select[c]]->ch[1], index++);
            }
            if (coup->l[c])
                cc_trans(ac, cc, &ac->che_cpe[coup->tag_select[c]]->ch[0], index++);
            if (coup->r[c])
                cc_trans(ac, cc, &ac->che_cpe[coup->tag_select[c]]->ch[1], index++);
        }
    }
}

static void coupling_tool(AACContext * ac, int independent, int domain) {
    int i;
    for (i = 0; i < MAX_TAGID; i++) {
        cc_struct * cc = ac->che_cc[i];
        if (cc != NULL) {
            if (cc->coup.ind_sw && independent) {
                transform_coupling_tool(ac, cc, coupling_independent_trans);
            } else if (!cc->coup.ind_sw && !independent && (cc->coup.domain == domain)) {
                transform_coupling_tool(ac, cc, coupling_dependent_trans);
            }
        }
    }
}

static void transform_sce_tool(AACContext * ac, void (*sce_trans)(AACContext * ac, sce_struct * sce)) {
    int i;
    for (i = 0; i < MAX_TAGID; i++) {
        if (ac->che_sce[i] != NULL)
            sce_trans(ac, ac->che_sce[i]);
        if (ac->che_cpe[i] != NULL) {
            sce_trans(ac, &ac->che_cpe[i]->ch[0]);
            sce_trans(ac, &ac->che_cpe[i]->ch[1]);
        }
        if (ac->che_lfe[i] != NULL)
            sce_trans(ac, ac->che_lfe[i]);
        if (ac->che_cc[i] != NULL)
            sce_trans(ac, &ac->che_cc[i]->ch);
    }
}

static void spec_to_sample(AACContext * ac) {
    coupling_tool(ac, 0, 0);
    if (ac->audioObjectType == AOT_AAC_LTP)
        transform_sce_tool(ac, ltp_trans);
    transform_sce_tool(ac, tns_trans);
    coupling_tool(ac, 0, 1);
    if (ac->audioObjectType != AOT_AAC_SSR)
        transform_sce_tool(ac, window_trans);
    else
        transform_sce_tool(ac, ssr_trans);
    coupling_tool(ac, 1, 1);
    if (ac->audioObjectType == AOT_AAC_LTP)
        transform_sce_tool(ac, ltp_update_trans);
}

static int output_coefs(AVCodecContext * avccontext) {
    AACContext * ac = avccontext->priv_data;
    program_config_struct * pcs = &ac->pcs;
    mix_config_struct * mix = &ac->mix;
    int ichannels = ac->channels;
    int ochannels = avccontext->channels;
    int i;
    for (i = 0; i < MAX_TAGID; i++) {
        mix->sce_gain[i] = 1.;
        mix->cpe_gain[i][0] = 1.;
        mix->cpe_gain[i][1] = 1.;
        mix->lfe_gain[i] = 1.;
    }
    mix->c_tag = 0;
    mix->lr_tag = 0;
    mix->sur_tag = 0;
    mix->mode = MIXMODE_UNKNOWN;
    if ((ochannels == 1) && ((ichannels == 2) || (ichannels == 1) || (pcs->mono_mixdown))) {
        int tag = pcs->mono_mixdown ? pcs->mono_mixdown - 1 : ((pcs->num_front == 1) ? pcs->front_tag[0] : 0);
        if (ichannels == 2) {
            mix->mode = MIXMODE_2TO1;
            mix->lr_tag = tag;
            mix->cpe_gain[tag][0] = mix->cpe_gain[tag][1] = 0.5;
        } else {
            mix->mode = MIXMODE_1TO1;
            mix->c_tag = tag;
        }
    } else if ((ochannels == 2) && ((ichannels == 1) || (ichannels == 2) || (pcs->stereo_mixdown))) {
        int tag = pcs->stereo_mixdown ? pcs->stereo_mixdown - 1 :  ((pcs->num_front == 1) ? pcs->front_tag[0] : 0);
        if (ichannels == 1) {
            mix->mode = MIXMODE_1TO2;
            mix->c_tag = tag;
        } else {
            mix->mode = MIXMODE_2TO2;
            mix->lr_tag = tag;
        }
    } else if (((ochannels == 1) || (ochannels == 2)) && (pcs->matrix_mixdown)) {
        float alpha_tab[] = {sqrt(2)/2, 1./2, sqrt(2)/4., 0};
        float alpha = alpha_tab[pcs->matrix_mixdown - 1];
        float ialpha = 0;
        if (ochannels == 1) {
            mix->mode = MIXMODE_MATRIX1;
            ialpha = 1. / (3 + 2 * alpha);
        } else {
            mix->mode = MIXMODE_MATRIX2;
            ialpha = pcs->pseudo_surround ? 1. / (1. + sqrt(2) / 2 + 2 * alpha) : 1. / (1. + sqrt(2) / 2 + alpha);
        }
        for (i = 0; i < pcs->num_front; i++) {
            if (pcs->front_cpe & (1 << i)) {
                mix->lr_tag = pcs->front_tag[i];
                mix->cpe_gain[mix->lr_tag][0] = mix->cpe_gain[mix->lr_tag][1] = ialpha;
                break;
            }
        }
        for (i = 0; i < pcs->num_front; i++) {
            if (!(pcs->front_cpe & (1 << i))) {
                mix->c_tag = pcs->front_tag[i];
                mix->sce_gain[mix->c_tag] = ialpha * sqrt(2) / 2.;
                break;
            }
        }
        mix->sur_tag = -1;
        for (i = 0; i < pcs->num_back; i++) {
            if (pcs->back_cpe & (1 << i)) {
                mix->sur_tag = pcs->back_tag[i];
                break;
            }
        }
        if (mix->sur_tag == -1) {
            for (i = 0; i < pcs->num_side; i++) {
                if (pcs->side_cpe & (1 << i)) {
                    mix->sur_tag = pcs->back_tag[i];
                    break;
                }
            }
        }
        if (mix->sur_tag != -1) {
            mix->cpe_gain[mix->sur_tag][0] = ialpha * alpha;
            mix->cpe_gain[mix->sur_tag][1] = ialpha * alpha;
        }
    } else if (ochannels >= ichannels) {
        mix->mode = MIXMODE_DEFAULT;
    }
    return 0;
}

static int output_samples(AVCodecContext * avccontext, uint16_t * data, int * data_size) {
    AACContext * ac = avccontext->priv_data;
    program_config_struct * pcs = &ac->pcs;
    mix_config_struct * mix = &ac->mix;
    int ichannels = ac->channels;
    int ochannels = avccontext->channels;
    int size = ochannels * 1024 * sizeof(uint16_t);
    int i;

    /* set a default float2int16 buffer */
    ac->iop = ac->interleaved_output;

    if (!ac->is_saved) {
        ac->is_saved = 1;
        *data_size = 0;
        return 0;
    }
    *data_size = size;

    /* the matrixmix modes are probably broken and they shouldn't be here anyway */
    switch (mix->mode) {
        case MIXMODE_DEFAULT:
            break;
        case MIXMODE_1TO1:
            ac->iop = ac->che_sce[mix->c_tag]->ret;
            break;
        case MIXMODE_2TO1:
            for (i = 0; i < 1024; i++)
                ac->interleaved_output[i] = ac->che_cpe[0]->ch[mix->lr_tag].ret[i] + ac->che_cpe[mix->lr_tag]->ch[1].ret[i];
            break;
        case MIXMODE_1TO2:
            for (i = 0; i < 1024; i++)
                ac->interleaved_output[i*2] = ac->interleaved_output[i*2+1] = ac->che_sce[mix->c_tag]->ret[i];
            break;
        case MIXMODE_2TO2:
            for (i = 0; i < 1024; i++) {
                ac->interleaved_output[i*2]    = ac->che_cpe[mix->lr_tag]->ch[0].ret[i];
                ac->interleaved_output[i*2+1]  = ac->che_cpe[mix->lr_tag]->ch[1].ret[i];
            }
            break;
        case MIXMODE_MATRIX1:
            {
                cpe_struct *ch_lr = ac->che_cpe[mix->lr_tag];
                sce_struct *ch_c = ac->che_sce[mix->c_tag];
                cpe_struct *ch_sur = ac->che_cpe[mix->sur_tag];
                float cBIAS = - ac->add_bias;
                float out[1024];
                if (ch_c) {
                    cBIAS += ac->add_bias;
                    for (i = 0; i < 1024; i++)
                        out[i] = ch_c->ret[i];
                } else {
                    memset(out, 0, sizeof(out));
                }
                if (ch_lr) {
                    cBIAS += 2 * ac->add_bias;
                    for (i = 0; i < 1024; i++)
                        out[i] += ch_lr->ch[0].ret[i] + ch_lr->ch[1].ret[i];
                }
                if (ch_sur) {
                    cBIAS += 2 * ac->add_bias;
                    for (i = 0; i < 1024; i++)
                        out[i] += ch_sur->ch[0].ret[i] + ch_sur->ch[1].ret[i];
                }
                for (i = 0; i < 1024; i++)
                    ac->interleaved_output[i] = out[i] - cBIAS;
            }
            break;
        case MIXMODE_MATRIX2:
            {
                cpe_struct *ch_lr = ac->che_cpe[mix->lr_tag];
                sce_struct *ch_c = ac->che_sce[mix->c_tag];
                cpe_struct *ch_sur = ac->che_cpe[mix->sur_tag];
                float lBIAS = -ac->add_bias, rBIAS = -ac->add_bias;
                float out[1024][2];
                if (ch_c) {
                    lBIAS += ac->add_bias; rBIAS += ac->add_bias;
                    for (i = 0; i < 1024; i++) {
                        out[i][0] = out[i][1] = ch_c->ret[i];
                    }
                } else {
                    memset(out, 0, sizeof(out));
                }
                if (ch_lr) {
                    lBIAS += ac->add_bias; rBIAS += ac->add_bias;
                    for (i = 0; i < 1024; i++) {
                        out[i][0] += ch_lr->ch[0].ret[i];
                        out[i][1] += ch_lr->ch[1].ret[i];
                    }
                }
                if (ch_sur) {
                    if (pcs->pseudo_surround) {
                        lBIAS -= 2 * ac->add_bias; rBIAS += 2 * ac->add_bias;
                        for (i = 0; i < 1024; i++) {
                            out[i][0] -= (ch_sur->ch[0].ret[i] + ch_sur->ch[1].ret[i]);
                            out[i][1] += (ch_sur->ch[1].ret[i] + ch_sur->ch[0].ret[i]);
                        }
                    } else {
                        lBIAS += ac->add_bias; rBIAS += ac->add_bias;
                        for (i = 0; i < 1024; i++) {
                            out[i][0] += ch_sur->ch[0].ret[i];
                            out[i][1] += ch_sur->ch[1].ret[i];
                        }
                    }
                }
                for (i = 0; i < 1024; i++) {
                    ac->interleaved_output[i*2]   = out[i][0] - lBIAS;
                    ac->interleaved_output[i*2+1] = out[i][1] - rBIAS;
                }
            }
            break;
        case MIXMODE_UNKNOWN:
        default:
            *data_size = 0;
            return 1;
    }
    if (mix->mode == MIXMODE_DEFAULT) {
        float *order[MAX_CHANNELS];
        int i, j = 0;
        if (pcs->present || pcs->generated) {
            for (i = 0; i < pcs->num_front; i++)
                if (!(pcs->front_cpe & (1 << i)))
                    order[j++] = ac->che_sce[pcs->front_tag[i]]->ret;
            for (i = 0; i < pcs->num_front; i++)
                if (pcs->front_cpe & (1 << i)) {
                    order[j++] = ac->che_cpe[pcs->front_tag[i]]->ch[0].ret;
                    order[j++] = ac->che_cpe[pcs->front_tag[i]]->ch[1].ret;
                }
            for (i = 0; i < pcs->num_side; i++)
                if (!(pcs->side_cpe & (1 << i))) {
                    order[j++] = ac->che_sce[pcs->side_tag[i]]->ret;
                } else {
                    order[j++] = ac->che_cpe[pcs->side_tag[i]]->ch[0].ret;
                    order[j++] = ac->che_cpe[pcs->side_tag[i]]->ch[1].ret;
                }
            for (i = 0; i < pcs->num_back; i++)
                if (pcs->back_cpe & (1 << i)) {
                    order[j++] = ac->che_cpe[pcs->back_tag[i]]->ch[0].ret;
                    order[j++] = ac->che_cpe[pcs->back_tag[i]]->ch[1].ret;
                }
            for (i = 0; i < pcs->num_back; i++)
                if (!(pcs->back_cpe & (1 << i)))
                    order[j++] = ac->che_sce[pcs->back_tag[i]]->ret;
            for (i = 0; i < pcs->num_lfe; i++)
                order[j++] = ac->che_lfe[pcs->lfe_tag[i]]->ret;
        } else {
            for (i = 0; i < MAX_TAGID; i++)
                if (ac->che_sce[i] != NULL)
                    order[j++] = ac->che_sce[i]->ret;
            for (i = 0; i < MAX_TAGID; i++)
                if (ac->che_cpe[i] != NULL) {
                    order[j++] = ac->che_cpe[i]->ch[0].ret;
                    order[j++] = ac->che_cpe[i]->ch[1].ret;
                }
            for (i = 0; i < MAX_TAGID; i++)
                if (ac->che_lfe[i] != NULL)
                    order[j++] = ac->che_lfe[i]->ret;
        }
        assert(j == ichannels);
        for (i = 0; i < ochannels; i++) {
            if (i < ichannels) {
                for (j = 0; j < 1024; j++)
                    ac->interleaved_output[j * ochannels + i] = order[i][j];
            } else {
                for (j = 0; j < 1024; j++)
                    ac->interleaved_output[j * ochannels + i] = 0.0;
            }
        }
    }

    /* Convert from float to int16 */
    ac->dsp.float_to_int16(data, ac->iop, 1024*ochannels);

    return 0;
}

static int aac_decode_frame(AVCodecContext * avccontext, void * data, int * data_size, uint8_t * buf, int buf_size) {
    AACContext * ac = avccontext->priv_data;
    GetBitContext * gb = &ac->gb;
    int id;
    int num_decoded = 0;

    ac->num_frame++;
    //if (ac->num_frame == 40)
    //    __asm int 3;

    init_get_bits(gb, buf, buf_size*8);
    //av_log(avccontext, AV_LOG_INFO, "%d ", buf_size);

    if (!ac->is_saved) {
        output_coefs(avccontext);
    }
    // parse
    while ((id = get_bits(gb, 3)) != ID_END) {
        switch (id) {
            case ID_SCE: {
                         if (!single_channel_struct(ac, gb))
                             num_decoded += 1;
                         break;
                     }
            case ID_CPE: {
                         if (!channel_pair_element(ac, gb))
                             num_decoded += 2;
                         break;
                     }
            case ID_FIL: {
                         int cnt = get_bits(gb, 4);
                         if (cnt == 15) cnt += get_bits(gb, 8) - 1;
                         while (cnt > 0)
                            cnt -= extension_payload(ac, gb, cnt);
                         break;
                     }
            case ID_PCE:
                     program_config_element(ac, gb);
                     break;
            case ID_DSE:
                     data_stream_element(ac, gb);
                     break;
            case ID_CCE:
                     coupling_channel_element(ac, gb);
                     break;
            case ID_LFE:
                     if (!lfe_channel_struct(ac, gb))
                         num_decoded += 1;
                     break;
            default:
                     assert(0 && 0);
                     break;
        }
    }

    spec_to_sample(ac);
    output_samples(avccontext, data, data_size);

    return buf_size;
}

static int aac_decode_close(AVCodecContext * avccontext) {
    AACContext * ac = avccontext->priv_data;
    int i;

    for (i = 0; i < MAX_TAGID; i++) {
        if (ac->che_sce[i]) {
            av_free(ac->che_sce[i]->ssr);
            av_free(ac->che_sce[i]->ltp_state);
            av_free(ac->che_sce[i]);
        }
        if (ac->che_cpe[i]) {
            av_free(ac->che_cpe[i]->ch[0].ssr);
            av_free(ac->che_cpe[i]->ch[1].ssr);
            av_free(ac->che_cpe[i]->ch[0].ltp_state);
            av_free(ac->che_cpe[i]->ch[1].ltp_state);
            av_free(ac->che_cpe[i]);
        }
        if (ac->che_lfe[i]) {
            av_free(ac->che_lfe[i]->ssr);
            av_free(ac->che_lfe[i]->ltp_state);
            av_free(ac->che_lfe[i]);
        }
        if (ac->che_cc[i]) {
            av_free(ac->che_cc[i]->ch.ssr);
            // ltp never used in cc
            av_free(ac->che_cc[i]);
        }
    }

    for (i = 0; i < 11; i++) {
        free_vlc(&ac->books[i]);
        av_free(ac->vq[i]);
    }
    av_free(ac->ssrctx);
    free_vlc(&ac->mainvlc);
    ff_mdct_end(&ac->mdct);
    ff_mdct_end(&ac->mdct_small);
    if (ac->mdct_ltp) {
        ff_mdct_end(ac->mdct_ltp);
        av_free(ac->mdct_ltp);
    }
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
};

