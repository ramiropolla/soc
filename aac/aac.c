/**
 * @file aac.c
 * AAC decoder
 * @author Oded Shimon  ( ods15 ods15 dyndns org )
 * @author Maxim Gavrilov ( maxim.gavrilov gmail com )

 * Kaiser-Bessel Derived Window by Justin Ruggles
 * Mersenne Twister by Kartikey Mahendra

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
 *
 */

#define V_DEBUG

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"

#include "aac.h"

#ifndef V_DEBUG
#define AV_DEBUG(...)
#endif

#undef NDEBUG
#include <assert.h>

// aux memory
static sce_struct * create_sce_struct() {
    return av_mallocz(sizeof(sce_struct));
}

static void free_sce_struct(sce_struct * sce) {
    if (sce != NULL)
        av_free(sce);
}

static cpe_struct * create_cpe_struct() {
    cpe_struct * cpe = av_mallocz(sizeof(cpe_struct));
    cpe->ch[0] = create_sce_struct();
    cpe->ch[1] = create_sce_struct();
    return cpe;
}

static void free_cpe_struct(cpe_struct * cpe) {
    if (cpe != NULL) {
        if (cpe->ch[0] != NULL)
            free_sce_struct(cpe->ch[0]);
        if (cpe->ch[1] != NULL)
            free_sce_struct(cpe->ch[1]);
        av_free(cpe);
    }
}

static cc_struct * create_cc_struct() {
    cc_struct * cc = av_malloc(sizeof(cc_struct));
    cc->ch = create_sce_struct();
    return cc;
}

static void free_cc_struct(cc_struct * cc) {
    if (cc != NULL) {
        if (cc->ch != NULL)
            free_sce_struct(cc->ch);
        av_free(cc);
    }
}


// Parsing predefines
static int program_config_element_add_channel(aac_context_t * ac, int flag_tag);
static int program_config_element(aac_context_t * ac, GetBitContext * gb);
static void ics_info(aac_context_t * ac, GetBitContext * gb, ics_struct * out);
static void section_data(aac_context_t * ac, GetBitContext * gb, ics_struct * ics, int cb[][64]);
static void scale_factor_data(aac_context_t * ac, GetBitContext * gb, int global_gain, ics_struct * ics, const int cb[][64], float sf[][64]);
static void pulse_data(aac_context_t * ac, GetBitContext * gb, pulse_struct * pulse);
static void tns_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, tns_struct * tns);
static int gain_control_data(aac_context_t * ac, GetBitContext * gb);
static int ms_data(aac_context_t * ac, GetBitContext * gb, cpe_struct * cpe);
static void spectral_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, const int cb[][64], const float sf[][64], int * icoef);
static int individual_channel_stream(aac_context_t * ac, GetBitContext * gb, int common_window, int scale_flag, sce_struct * sce);
static int single_channel_struct(aac_context_t * ac, GetBitContext * gb);
static int channel_pair_element(aac_context_t * ac, GetBitContext * gb);
static int coupling_channel_element(aac_context_t * ac, GetBitContext * gb);
static int lfe_channel_struct(aac_context_t * ac, GetBitContext * gb);

// Tools predefines
static void quant_to_spec_tool(aac_context_t * ac, const ics_struct * ics, const int * icoef, const int cb[][64], const float sf[][64], float * coef);
static void pulse_tool(aac_context_t * ac, const ics_struct * ics, const pulse_struct * pulse, int * icoef);
static void ms_tool(aac_context_t * ac, cpe_struct * cpe);
static void intensity_tool(aac_context_t * ac, cpe_struct * cpe);
static void coupling_tool(aac_context_t * ac, int independent, int domain);
static void spec_to_sample(aac_context_t * ac);

// Transformations
static void coupling_dependent_trans(aac_context_t * ac, cc_struct * cc, sce_struct * sce, int index);
static void coupling_independent_trans(aac_context_t * ac, cc_struct * cc, sce_struct * sce, int index);
static void tns_trans(aac_context_t * ac, sce_struct * sce);
static void window_trans(aac_context_t * ac, sce_struct * sce);

// Generic transformations
static void transform_sce_tool(aac_context_t * ac, void (*f)(aac_context_t * ac, sce_struct * sce));
static void transform_coupling_tool(aac_context_t * ac, cc_struct * cc, void (*f)(aac_context_t * ac, cc_struct * cc, sce_struct * sce, int index));

// Output
static int output_samples(AVCodecContext * avccontext, void * data, int * data_size);


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

/* BEGIN Mersenne Twister Code. */
static void dither_seed(dither_state *state, uint32_t seed) {
    if (seed == 0)
        seed = 0x1f2e3d4c;

    state->mt[0] = seed;
    for (state->mti = 1; state->mti < N; state->mti++)
        state->mt[state->mti] = ((69069 * state->mt[state->mti - 1]) + 1);
}

static uint32_t dither_uint32(dither_state *state) {
    uint32_t y;
    static const uint32_t mag01[2] = { 0x00, MATRIX_A };
    int kk;

    if (state->mti >= N) {
        for (kk = 0; kk < N - M; kk++) {
            y = (state->mt[kk] & UPPER_MASK) | (state->mt[kk + 1] & LOWER_MASK);
            state->mt[kk] = state->mt[kk + M] ^ (y >> 1) ^ mag01[y & 0x01];
        }
        for (;kk < N - 1; kk++) {
            y = (state->mt[kk] & UPPER_MASK) | (state->mt[kk + 1] & LOWER_MASK);
            state->mt[kk] = state->mt[kk + (M - N)] ^ (y >> 1) ^ mag01[y & 0x01];
        }
        y = (state->mt[N - 1] & UPPER_MASK) | (state->mt[0] & LOWER_MASK);
        state->mt[N - 1] = state->mt[M - 1] ^ (y >> 1) ^ mag01[y & 0x01];

        state->mti = 0;
    }

    y = state->mt[state->mti++];
    y ^= (y >> 11);
    y ^= ((y << 7) & 0x9d2c5680);
    y ^= ((y << 15) & 0xefc60000);
    y ^= (y >> 18);

    return y;
}

static inline int16_t dither_int16(dither_state *state) {
    return ((dither_uint32(state) << 16) >> 16);
}

//static inline int32_t dither_int32(dither_state *state) {
//    return (int32_t)dither_uint32(state);
//}

/* END Mersenne Twister */

// General functions
#define TAG_MASK 0x00f
#define FLAG_SCE 0x100
#define FLAG_CPE 0x200
#define FLAG_LFE 0x400
#define FLAG_CCE 0x800

static int program_config_element_add_channel(aac_context_t * ac, int flag_tag) {
    program_config_struct * pcs = &ac->pcs;
    if (pcs->present)
        return 0;
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
    return 1;
}

static int program_config_element(aac_context_t * ac, GetBitContext * gb) {
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

static int GASpecificConfig(aac_context_t * ac, GetBitContext * gb) {
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
    *index= get_bits(gb, 4);
    if (*index == 0xf) {
        *index = -1;
        *rate = get_bits(gb, 24);
    } else {
        assert(*index <= 12);
        *rate = sampling_table[*index];
    }
    return 0;
}

static int AudioSpecificConfig(aac_context_t * ac, void *data, int data_size) {
    GetBitContext * gb = &ac->gb;

    init_get_bits(gb, data, data_size * 8);

    memset(&ac->pcs, 0, sizeof(ac->pcs));

    ac->audioObjectType = GetAudioObjectType(gb);
    assert(ac->audioObjectType == AOT_AAC_LC);
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
        const uint16_t (*a)[2];
    const unsigned int s;
    } tmp[] = {
        { codebook1 , sizeof codebook1  },
        { codebook2 , sizeof codebook2  },
        { codebook3 , sizeof codebook3  },
        { codebook4 , sizeof codebook4  },
        { codebook5 , sizeof codebook5  },
        { codebook6 , sizeof codebook6  },
        { codebook7 , sizeof codebook7  },
        { codebook8 , sizeof codebook8  },
        { codebook9 , sizeof codebook9  },
        { codebook10, sizeof codebook10 },
        { codebook11, sizeof codebook11 },
    };
    aac_context_t * ac = avccontext->priv_data;
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

    for (i = 0; i < 11; i++) {
        static const int mod_cb[11] = { 3, 3, 3, 3, 9, 9, 8, 8, 13, 13, 17 };
        static const int off_cb[11] = { 1, 1, 0, 0, 4, 4, 0, 0,  0,  0,  0 };

        int j, values = tmp[i].s/sizeof(tmp[i].a[0]);
        int dim = (i >= 4 ? 2 : 4);
        int mod = mod_cb[i], off = off_cb[i], index = 0;
        int ret;
        ret = init_vlc(&ac->books[i], 6, values,
                 &tmp[i].a[0][1], sizeof(tmp[i].a[0]), sizeof(tmp[i].a[0][1]),
                 &tmp[i].a[0][0], sizeof(tmp[i].a[0]), sizeof(tmp[i].a[0][0]),
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
    // dither init
    dither_seed(&ac->dither, 0);
    // windows init
    kbd_window_init(4, ac->kbd_long_1024, 2048, 50);
    kbd_window_init(6, ac->kbd_short_128, 256, 50);
    sine_window_init(ac->sine_long_1024, 2048);
    sine_window_init(ac->sine_short_128, 256);
    // 1024  - compensate wrong imdct method
    // 32768 - values in AAC build for ready float->int 16 bit audio, using
    // BIAS method instead needs values -1<x<1
    for (i = 0; i < 256; i++)
        ac->pow2sf_tab[i] = pow(2, (i - 100)/4.) /1024./32768.;
    for (i = 0; i < 256; i++)
        ac->intensity_tab[i] = pow(0.5, (i - 100) / 4.);
    for (i = 0; i < sizeof(ac->ivquant_tab)/sizeof(ac->ivquant_tab[0]); i++)
        ac->ivquant_tab[i] = pow(i, 4./3);
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

    init_vlc(&ac->mainvlc, 7, sizeof(aac_scalefactor_huffman_table)/sizeof(aac_scalefactor_huffman_table[0]),
             &aac_scalefactor_huffman_table[0][1], sizeof(aac_scalefactor_huffman_table[0]), sizeof(aac_scalefactor_huffman_table[0][1]),
             &aac_scalefactor_huffman_table[0][0], sizeof(aac_scalefactor_huffman_table[0]), sizeof(aac_scalefactor_huffman_table[0][0]),
             0);

    ff_mdct_init(&ac->mdct, 11, 1);
    ff_mdct_init(&ac->mdct_small, 8, 1);

    return 0;
}

// Parsers implementation
static int data_stream_element(aac_context_t * ac, GetBitContext * gb) {
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

static void ics_info(aac_context_t * ac, GetBitContext * gb, ics_struct * ics) {
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
        ics->predictor = get_bits1(gb);
        assert(ics->predictor == 0);
        ics->swb_offset = ac->swb_offset_1024;
        ics->num_swb = ac->num_swb_1024;
        ics->num_windows = 1;
        ics->tns_max_bands = ac->tns_max_bands_1024;
    }
}

static inline float ivquant(aac_context_t * ac, int a) {
    static const float sign[2] = { -1., 1. };
    int tmp = (a>>31);
    int abs_a = (a^tmp)-tmp;
    if (abs_a < sizeof(ac->ivquant_tab)/sizeof(ac->ivquant_tab[0]))
        return sign[tmp+1] * ac->ivquant_tab[abs_a];
    else
        return sign[tmp+1] * pow(abs_a, 4./3);
}

static void section_data(aac_context_t * ac, GetBitContext * gb, ics_struct * ics, int cb[][64]) {
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

static void scale_factor_data(aac_context_t * ac, GetBitContext * gb, int global_gain, ics_struct * ics, const int cb[][64], float sf[][64]) {
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
                sf[g][i] = pow(2.0, 0.25 * noise)/1024./32768.;
            } else {
                global_gain += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                assert(!(global_gain & (~255)));
                sf[g][i] = ac->pow2sf_tab[global_gain];
            }
        }
    }
}

static void pulse_data(aac_context_t * ac, GetBitContext * gb, pulse_struct * pulse) {
    int i;
    pulse->num_pulse = get_bits(gb, 2);
    pulse->start = get_bits(gb, 6);
    for (i = 0; i <= pulse->num_pulse; i++) {
        pulse->offset[i] = get_bits(gb, 5);
        pulse->amp[i] = get_bits(gb, 4);
    }
}

static void tns_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, tns_struct * tns) {
    int w, filt, i, coef_len, coef_res = 0, coef_compress;
    for (w = 0; w < ics->num_windows; w++) {
        tns->n_filt[w] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 1 : 2);
        if (tns->n_filt[w])
            coef_res = get_bits1(gb) + 3;
        for (filt = 0; filt < tns->n_filt[w]; filt++) {
            tns->length[w][filt] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 4 : 6);
            if (tns->order[w][filt] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 3 : 5)) {
                tns->direction[w][filt] = get_bits1(gb);
                assert(coef_res == 3 || coef_res == 4);
                coef_compress = get_bits1(gb);
                coef_len = coef_res - coef_compress;
                tns->tmp2_map = tns_tmp2_map[(coef_compress << 1) + (coef_res - 3)];
                for (i = 0; i < tns->order[w][filt]; i++)
                    tns->coef[w][filt][i] = get_bits(gb, coef_len);
            }
        }
    }
}

static int gain_control_data(aac_context_t * ac, GetBitContext * gb) {
    // FIXME ignored
    av_log(ac->avccontext, AV_LOG_INFO, " gain data ignored\n");
    return 0;
}

static int ms_data(aac_context_t * ac, GetBitContext * gb, cpe_struct * cpe) {
    ms_struct * ms = &cpe->ms;
    ms->present = get_bits(gb, 2);
    if (ms->present == 1) {
        int g, i;
        for (g = 0; g < cpe->ch[0]->ics.num_window_groups; g++)
            for (i = 0; i < cpe->ch[0]->ics.max_sfb; i++)
                ms->mask[g][i] = get_bits1(gb);// << i;
    } else if (ms->present == 2) {
        int g, i;
        for (g = 0; g < cpe->ch[0]->ics.num_window_groups; g++)
            for (i = 0; i < cpe->ch[0]->ics.max_sfb; i++)
                ms->mask[g][i] = 1;// = 0xFFFFFFFFFFFFFFFFULL;
    }
    return 0;
}

static void spectral_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, const int cb[][64], const float sf[][64], int * icoef) {
    static const int unsigned_cb[] = { 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1 };
    int i, k, g;
    const uint16_t * offsets = ics->swb_offset;
    for (g = 0; g < ics->num_window_groups; g++) {
        int total = (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 128 : 1024;
        memset(icoef + g*total + offsets[ics->max_sfb], 0, sizeof(int)*(total - offsets[ics->max_sfb]));
    }
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
                        icoef[group*128+k] = dither_int16(&ac->dither);
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

static int individual_channel_stream(aac_context_t * ac, GetBitContext * gb, int common_window, int scale_flag, sce_struct * sce) {
    int icoeffs[1024];
    pulse_struct pulse;
    tns_struct * tns = &sce->tns;
    ics_struct * ics = &sce->ics;
    float * out = sce->coeffs;

    //memset(sf, 0, sizeof(sf));

    sce->global_gain = get_bits(gb, 8);

    if (!common_window && !scale_flag) {
        ics_info(ac, gb, ics);
    }

    //av_log(ac->avccontext, AV_LOG_INFO, " global_gain: %d, groups: %d\n", global_gain, ics->window_sequence);
    section_data(ac, gb, ics, sce->cb);
    scale_factor_data(ac, gb, sce->global_gain, ics, sce->cb, sce->sf);

    if (!scale_flag) {
        if (pulse.present = get_bits1(gb))
            pulse_data(ac, gb, &pulse);
        if (tns->present = get_bits1(gb))
            tns_data(ac, gb, ics, tns);
        if (get_bits1(gb))
            if (gain_control_data(ac, gb)) return 1;
    }

    spectral_data(ac, gb, ics, sce->cb, sce->sf, icoeffs);
    pulse_tool(ac, ics, &pulse, icoeffs);
    quant_to_spec_tool(ac, ics, icoeffs, sce->cb, sce->sf, out);
    return 0;
}

static int single_channel_struct(aac_context_t * ac, GetBitContext * gb) {
    sce_struct * sce;
    int id = get_bits(gb, 4);
    program_config_element_add_channel(ac, FLAG_SCE | id);
    if (ac->che_sce[id] == NULL)
        ac->che_sce[id] = create_sce_struct();
    sce = ac->che_sce[id];
    if (individual_channel_stream(ac, gb, 0, 0, sce))
        return 1;
    return 0;
}

static int channel_pair_element(aac_context_t * ac, GetBitContext * gb) {
    int i, common_window;
    cpe_struct * cpe;
    int id = get_bits(gb, 4);
    program_config_element_add_channel(ac, FLAG_CPE | id);
    if (ac->che_cpe[id] == NULL)
        ac->che_cpe[id] = create_cpe_struct();
    cpe = ac->che_cpe[id];
    //assert(id == 0);
    common_window = get_bits1(gb);
    if (common_window) {
        ics_info(ac, gb, &cpe->ch[0]->ics);
        i = cpe->ch[1]->ics.window_shape_prev;
        cpe->ch[1]->ics = cpe->ch[0]->ics;
        cpe->ch[1]->ics.window_shape_prev = i;
        ms_data(ac, gb, cpe);
    } else {
        cpe->ms.present = 0;
    }
    if (individual_channel_stream(ac, gb, common_window, 0, cpe->ch[0]))
        return 1;
    if (individual_channel_stream(ac, gb, common_window, 0, cpe->ch[1]))
        return 1;

    // M/S tool
    if (common_window) {
        ms_tool(ac, cpe);
    }
    intensity_tool(ac, cpe);
    return 0;
}

static int coupling_channel_element(aac_context_t * ac, GetBitContext * gb) {
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
    program_config_element_add_channel(ac, FLAG_CCE | id);
    if (ac->che_cc[id] == NULL)
        ac->che_cc[id] = create_cc_struct();
    sce = ac->che_cc[id]->ch;
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

static int lfe_channel_struct(aac_context_t * ac, GetBitContext * gb) {
    sce_struct * sce;
    int id = get_bits(gb, 4);
    program_config_element_add_channel(ac, FLAG_LFE | id);
    if (ac->che_lfe[id] == NULL)
        ac->che_lfe[id] = create_sce_struct();
    sce = ac->che_lfe[id];
    if (individual_channel_stream(ac, gb, 0, 0, sce))
        return 1;
    return 0;
}

// Tools implementation
static void quant_to_spec_tool(aac_context_t * ac, const ics_struct * ics, const int * icoef, const int cb[][64], const float sf[][64], float * coef) {
    const uint16_t * offsets = ics->swb_offset;
    int g, i, group, k;
    for (g = 0; g < ics->num_window_groups; g++) {
        int total = (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 128 : 1024;
        memset(coef + g*total + offsets[ics->max_sfb], 0, sizeof(float)*(total - offsets[ics->max_sfb]));
    }
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

static void pulse_tool(aac_context_t * ac, const ics_struct * ics, const pulse_struct * pulse, int * icoef) {
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

static void ms_tool(aac_context_t * ac, cpe_struct * cpe) {
    const ms_struct * ms = &cpe->ms;
    const ics_struct * ics = &cpe->ch[0]->ics;
    float *ch0 = cpe->ch[0]->coeffs;
    float *ch1 = cpe->ch[1]->coeffs;
    if (ms->present) {
        int g, i, k, start = 0, gp;
        const uint16_t * offsets = ics->swb_offset;
        for (g = 0; g < ics->num_window_groups; g++) {
            //av_log(ac->avccontext, AV_LOG_INFO, " masking[%d]: ", g);
            for (gp = start; gp < start + ics->group_len[g]; gp++) {
                for (i = 0; i < ics->max_sfb; i++) {
                    if (ms->mask[g][i]) {
                        for (k = offsets[i] + gp*128; k < offsets[i+1] + gp*128; k++) {
                            float tmp = ch0[k] - ch1[k];
                            ch0[k] += ch1[k];
                            ch1[k] = tmp;
                        }
                    }
                }
            }
            start += ics->group_len[g];
            //av_log(ac->avccontext, AV_LOG_INFO, "\n");
        }
    }
}

static void intensity_tool(aac_context_t * ac, cpe_struct * cpe) {
    const ics_struct * ics = &cpe->ch[1]->ics;
    sce_struct * sce0 = cpe->ch[0];
    sce_struct * sce1 = cpe->ch[1];
    if (ics->intensity_present) {
        const uint16_t * offsets = ics->swb_offset;
        int g, gp, i, k, start = 0;
        int c;
        float scale;
        for (g = 0; g < ics->num_window_groups; g++) {
            for (gp = start; gp < start + ics->group_len[g]; gp++) {
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
            start += ics->group_len[g];
        }
    }
}

static void tns_trans(aac_context_t * ac, sce_struct * sce) {
    const ics_struct * ics = &sce->ics;
    const tns_struct * tns = &sce->tns;
    float * coef = sce->coeffs;
    const int mmm = (ics->tns_max_bands > ics->max_sfb) ? ics->max_sfb : ics->tns_max_bands;
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
            order = tns->order[w][filt];
            if (order > TNS_MAX_ORDER)
                order = TNS_MAX_ORDER;
            if (order == 0)
                continue;

            // tns_decode_coef
            lpc[0] = 1;
            for (m = 1; m <= order; m++) {
                lpc[m] = tns->tmp2_map[tns->coef[w][filt][m - 1]];
                for (i = 1; i < m; i++)
                    b[i] = lpc[i] + lpc[m] * lpc[m-i];
                for (i = 1; i < m; i++)
                    lpc[i] = b[i];
            }

            start = ics->swb_offset[(bottom > mmm) ? mmm : bottom];
            end = ics->swb_offset[(top > mmm) ? mmm : top];
            if ((size = end - start) <= 0)
                continue;
            if (tns->direction[w][filt]) {
                inc = -1; start = end - 1;
            } else {
                inc = 1;
            }

            // ar filter
            memset(b, 0, sizeof(b));
            ib = 0;
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
        }
    }
}

#define BIAS 385

static void window_trans(aac_context_t * ac, sce_struct * sce) {
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
                for (i = 0; i < 1024; i++)     out[i] = saved[i] + buf[i] * lwindow_prev[i           ] + 0;//BIAS;
            } else {
                for (i = 0; i < 448; i++)      out[i] = saved[i] +                                  0;//BIAS;
                for (i = 448; i < 576; i++)    out[i] = saved[i] + buf[i] * swindow_prev[i - 448] + 0;//BIAS;
                for (i = 576; i < 1024; i++)   out[i] =            buf[i]                         + 0;//BIAS;
            }
        }
        if (ics->window_sequence != LONG_START_SEQUENCE) {
            for (i = 0; i < 1024; i++)     saved[i] = buf[i+1024] * lwindow[1023-i];
        } else {
            for (i = 0; i < 448; i++)      saved[i] = buf[i+1024];
            for (i = 448; i < 576; i++)    saved[i] = buf[i+1024] * swindow[127 - (i - 448)];
            //for (i = 576; i < 1024; i++)   saved[i] = 0.0;
        }
    } else {
        float * ptr[8];
        int i;
        for (i = 0; i < 8; i++) {
            int k;
            ptr[i] = buf + i*256;
            ff_imdct_calc(&ac->mdct_small, ptr[i], in + i*128, out);
            if (i == 0) {
                for (k = 0; k < 128; k++)   ptr[i][k] *= swindow_prev[k] * 8.;
                for (k = 128; k < 256; k++) ptr[i][k] *= swindow[255 - k] * 8.;
            } else {
                for (k = 0; k < 128; k++)   ptr[i][k] *= swindow[k] * 8.;
                for (k = 128; k < 256; k++)   ptr[i][k] *= swindow[255 - k] * 8.;
            }
        }
        for (i = 0; i < 448; i++)    out[i] = saved[i]                           + 0;//BIAS;
        for (; i < 576; i++)         out[i] = saved[i]         + ptr[0][i - 448] + 0;//BIAS;
        for (; i < 704; i++)         out[i] = ptr[0][i - 448]  + ptr[1][i - 576] + 0;//BIAS;
        for (; i < 832; i++)         out[i] = ptr[1][i - 576]  + ptr[2][i - 704] + 0;//BIAS;
        for (; i < 960; i++)         out[i] = ptr[2][i - 704]  + ptr[3][i - 832] + 0;//BIAS;
        for (; i < 1024; i++)        out[i] = ptr[3][i - 832]  + ptr[4][i - 960] + 0;//BIAS;

        for (; i < 1088; i++) saved[i-1024] = ptr[3][i - 832]  + ptr[4][i - 960];
        for (; i < 1216; i++) saved[i-1024] = ptr[4][i - 960]  + ptr[5][i - 1088];
        for (; i < 1344; i++) saved[i-1024] = ptr[5][i - 1088] + ptr[6][i - 1216];
        for (; i < 1472; i++) saved[i-1024] = ptr[6][i - 1216] + ptr[7][i - 1344];
        for (; i < 1600; i++) saved[i-1024] = ptr[7][i - 1344]; // memcpy?
        for (; i < 2048; i++) saved[i-1024] = 0.0;
    }
}


static void coupling_dependent_trans(aac_context_t * ac, cc_struct * cc, sce_struct * sce, int index) {
    ics_struct * ics = &cc->ch->ics;
    const uint16_t * offsets = ics->swb_offset;
    float * dest = sce->coeffs;
    float * src = cc->ch->coeffs;
    int g, i, group, k;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            if (cc->ch->cb[g][i] != ZERO_HCB) {
                float gain = cc->coup.gain[index][g][i];
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

static void coupling_independent_trans(aac_context_t * ac, cc_struct * cc, sce_struct * sce, int index) {
    int i;
    float gain = cc->coup.gain[index][0][0];
    for (i = 0; i < 1024; i++)
        sce->ret[i] += gain * (cc->ch->ret[i] - BIAS);
}

static void coupling_tool(aac_context_t * ac, int independent, int domain) {
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

static void spec_to_sample(aac_context_t * ac) {
    coupling_tool(ac, 0, 0);
    transform_sce_tool(ac, tns_trans);
    coupling_tool(ac, 0, 1);
    transform_sce_tool(ac, window_trans);
    coupling_tool(ac, 1, 1);
}

static void transform_sce_tool(aac_context_t * ac, void (*f)(aac_context_t * ac, sce_struct * sce)) {
    int i;
    for (i = 0; i < MAX_TAGID; i++) {
        if (ac->che_sce[i] != NULL)
            f(ac, ac->che_sce[i]);
        if (ac->che_cpe[i] != NULL) {
            f(ac, ac->che_cpe[i]->ch[0]);
            f(ac, ac->che_cpe[i]->ch[1]);
        }
        if (ac->che_lfe[i] != NULL)
            f(ac, ac->che_lfe[i]);
        if (ac->che_cc[i] != NULL)
            f(ac, ac->che_cc[i]->ch);
    }
}

static void transform_coupling_tool(aac_context_t * ac, cc_struct * cc,
                                    void (*f)(aac_context_t * ac, cc_struct * cc, sce_struct * sce, int index))
{
    int c;
    int index = 0;
    coupling_struct * coup = &cc->coup;
    for (c = 0; c <= coup->num_coupled; c++) {
        if (!coup->is_cpe[c]) {
            assert(ac->che_sce[coup->tag_select[c]] != NULL);
            f(ac, cc, ac->che_sce[coup->tag_select[c]], index++);
        } else {
            assert(ac->che_cpe[coup->tag_select[c]] != NULL);
            if (!coup->l[c] && !coup->r[c]) {
                f(ac, cc, ac->che_cpe[coup->tag_select[c]]->ch[0], index);
                f(ac, cc, ac->che_cpe[coup->tag_select[c]]->ch[1], index++);
            }
            if (coup->l[c])
                f(ac, cc, ac->che_cpe[coup->tag_select[c]]->ch[0], index++);
            if (coup->r[c])
                f(ac, cc, ac->che_cpe[coup->tag_select[c]]->ch[1], index++);
        }
    }
}

static inline uint16_t F2U16(float x) {
    int32_t tmp = 0;
    x += BIAS;
    tmp = ((int32_t*)&x)[0];
    if (tmp & 0xf0000) {
        if (tmp > 0x43c0ffff) tmp = 0xFFFF;
        else                  tmp = 0;
    }
    return (uint16_t)(tmp - 0x8000);
}

static int output_samples(AVCodecContext * avccontext, void * data, int * data_size) {
    aac_context_t * ac = avccontext->priv_data;
    program_config_struct * pcs = &ac->pcs;
    int ichannels = ac->channels;
    int ochannels = avccontext->channels;
    int size = ochannels * 1024 * sizeof(uint16_t);
    int i;

    if (!ac->is_saved) {
        ac->is_saved = 1;
        *data_size = 0;
        return 0;
    }
    *data_size = size;
    if ((ochannels == 1) && ((ichannels == 2) || (ichannels == 1) || (pcs->mono_mixdown))) {
        int tag = pcs->mono_mixdown ? pcs->mono_mixdown - 1 : 0;
        if (ichannels == 2) {
            for (i = 0; i < 1024; i++)
                ((uint16_t *)data)[i] = F2U16(0.5 * (ac->che_cpe[0]->ch[0]->ret[i] + ac->che_cpe[0]->ch[1]->ret[i]));
        } else {
            for (i = 0; i < 1024; i++)
                ((uint16_t *)data)[i] = F2U16(ac->che_sce[tag]->ret[i]);
        }
    } else if ((ochannels == 2) && ((ichannels == 1) || (ichannels == 2) || (pcs->stereo_mixdown))) {
        if (ichannels == 1) {
            for (i = 0; i < 1024; i++) {
                ((uint16_t(*)[2])data)[i][0] = ((uint16_t(*)[2])data)[i][1] = F2U16(ac->che_sce[0]->ret[i]);
            }
        } else {
            int tag = pcs->stereo_mixdown ? pcs->stereo_mixdown - 1 : 0;
            for (i = 0; i < 1024; i++) {
                ((uint16_t(*)[2])data)[i][0] = F2U16(ac->che_cpe[tag]->ch[0]->ret[i]);
                ((uint16_t(*)[2])data)[i][1] = F2U16(ac->che_cpe[tag]->ch[1]->ret[i]);
            }
        }
    } else if (((ochannels == 1) || (ochannels == 2)) && (pcs->matrix_mixdown)) {
        float alpha_tab[] = {sqrt(2)/2, 1./2, sqrt(2)/4., 0};
        float alpha = alpha_tab[pcs->matrix_mixdown - 1];
        sce_struct *ch_l = NULL, *ch_r = NULL, *ch_c = NULL, *ch_ls = NULL, *ch_rs = NULL;
        for (i = 0; i < pcs->num_front; i++) {
            if (pcs->front_cpe & (1 << i)) {
                ch_l = ac->che_cpe[pcs->front_tag[i]]->ch[0];
                ch_r = ac->che_cpe[pcs->front_tag[i]]->ch[1];
                break;
            }
        }
        for (i = 0; i < pcs->num_front; i++) {
            if (!(pcs->front_cpe & (1 << i)) && (ch_c == NULL)) {
                ch_c = ac->che_sce[pcs->front_tag[0]];
                break;
            }
        }
        for (i = 0; i < pcs->num_back; i++) {
            if (pcs->back_cpe & (1 << i)) {
                ch_ls = ac->che_cpe[pcs->back_tag[i]]->ch[0];
                ch_rs = ac->che_cpe[pcs->back_tag[i]]->ch[1];
                break;
            }
        }
        if (ch_ls == NULL) {
            for (i = 0; i < pcs->num_side; i++) {
                if (pcs->side_cpe & (1 << i)) {
                    ch_ls = ac->che_cpe[pcs->side_tag[i]]->ch[0];
                    ch_rs = ac->che_cpe[pcs->side_tag[i]]->ch[1];
                    break;
                }
            }
        }
        if (ochannels == 1) {
            float ialpha = 1. / (3 + 2 * alpha);
            float out[1024];
            memset(out, 0, sizeof(out));
            if (ch_c) for (i = 0; i < 1024; i++)
                out[i] = ch_c->ret[i];
            if (ch_l) for (i = 0; i < 1024; i++)
                out[i] += ch_l->ret[i] + ch_r->ret[i];
            if (ch_ls) for (i = 0; i < 1024; i++)
                out[i] += alpha * (ch_ls->ret[i] + ch_rs->ret[i]);
            for (i = 0; i < 1024; i++)
                ((uint16_t *)data)[i] = F2U16(ialpha * out[i]);
        }
        if (ochannels == 2) {
            float ialpha = pcs->pseudo_surround ? 1. / (1. + sqrt(2) / 2 + 2 * alpha) : 1. / (1. + sqrt(2) / 2 + alpha);
            float c = sqrt(2) / 2.;
            float out[1024][2];
            memset(out, 0, sizeof(out));
            if (ch_c) for (i = 0; i < 1024; i++) {
                out[i][0] = out[i][1] = c * ch_c->ret[i];
            }
            if (ch_l) for (i = 0; i < 1024; i++) {
                out[i][0] += ch_l->ret[i];
                out[i][1] += ch_r->ret[i];
            }
            if (ch_ls) {
                if (pcs->pseudo_surround)
                    for (i = 0; i < 1024; i++) {
                        out[i][0] -= alpha * (ch_ls->ret[i] + ch_rs->ret[i]);
                        out[i][1] += alpha * (ch_ls->ret[i] + ch_rs->ret[i]);
                    }
                else
                    for (i = 0; i < 1024; i++) {
                        out[i][0] += alpha * ch_ls->ret[i];
                        out[i][1] += alpha * ch_rs->ret[i];
                    }
            }
            for (i = 0; i < 1024; i++) {
                ((uint16_t(*)[2])data)[i][0] = F2U16(ialpha * out[i][0]);
                ((uint16_t(*)[2])data)[i][1] = F2U16(ialpha * out[i][1]);
            }
        }
    } else if (ichannels <= ochannels) {
        float *order[MAX_CHANNELS];
        int i, j = 0;
        if (pcs->present) {
            for (i = 0; i < pcs->num_front; i++)
                if (!(pcs->front_cpe & (1 << i)))
                    order[j++] = ac->che_sce[pcs->front_tag[i]]->ret;
            for (i = 0; i < pcs->num_front; i++)
                if (pcs->front_cpe & (1 << i)) {
                    order[j++] = ac->che_cpe[pcs->front_tag[i]]->ch[0]->ret;
                    order[j++] = ac->che_cpe[pcs->front_tag[i]]->ch[1]->ret;
                }
            for (i = 0; i < pcs->num_side; i++)
                if (!(pcs->side_cpe & (1 << i))) {
                    order[j++] = ac->che_sce[pcs->side_tag[i]]->ret;
                } else {
                    order[j++] = ac->che_cpe[pcs->side_tag[i]]->ch[0]->ret;
                    order[j++] = ac->che_cpe[pcs->side_tag[i]]->ch[1]->ret;
                }
            for (i = 0; i < pcs->num_back; i++)
                if (pcs->back_cpe & (1 << i)) {
                    order[j++] = ac->che_cpe[pcs->back_tag[i]]->ch[0]->ret;
                    order[j++] = ac->che_cpe[pcs->back_tag[i]]->ch[1]->ret;
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
                    order[j++] = ac->che_cpe[i]->ch[0]->ret;
                    order[j++] = ac->che_cpe[i]->ch[1]->ret;
                }
            for (i = 0; i < MAX_TAGID; i++)
                if (ac->che_lfe[i] != NULL)
                    order[j++] = ac->che_lfe[i]->ret;
        }
        assert(j == ichannels);
        for (i = 0; i < ochannels; i++) {
            if (i < ichannels) {
                for (j = 0; j < 1024; j++)
                    ((uint16_t *)data)[j * ochannels + i] = F2U16(order[i][j]);
            } else {
                for (j = 0; j < 1024; j++)
                    ((uint16_t *)data)[j * ochannels + i] = 0;
            }
        }
    } else {
        *data_size = 0;
        return 1;
    }
    return 0;
}

static int aac_decode_frame(AVCodecContext * avccontext, void * data, int * data_size, uint8_t * buf, int buf_size) {
    aac_context_t * ac = avccontext->priv_data;
    GetBitContext * gb = &ac->gb;
    int id;
    int num_decoded = 0;

    ac->num_frame++;
    //if (ac->num_frame == 40)
    //    __asm int 3;

    init_get_bits(gb, buf, buf_size*8);
    //av_log(avccontext, AV_LOG_INFO, "%d ", buf_size);

    // parse
    while ((id = get_bits(gb, 3)) != 7) {
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
                skip_bits(gb, 8*cnt);
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
    aac_context_t * ac = avccontext->priv_data;
    int i;

    for (i = 0; i < MAX_TAGID; i++) {
        if (ac->che_sce[i])
            free_sce_struct(ac->che_sce[i]);
        if (ac->che_cpe[i])
            free_cpe_struct(ac->che_cpe[i]);
        if (ac->che_lfe[i])
            free_sce_struct(ac->che_lfe[i]);
        if (ac->che_cc[i])
            free_cc_struct(ac->che_cc[i]);
    }

    for (i = 0; i < 11; i++) {
        free_vlc(&ac->books[i]);
        av_free(ac->vq[i]);
    }
    free_vlc(&ac->mainvlc);
    ff_mdct_end(&ac->mdct);
    ff_mdct_end(&ac->mdct_small);

    return 0 ;
}

AVCodec native_aac_decoder = {
    "native_aac",
    CODEC_TYPE_AUDIO,
    CODEC_ID_NATIVE_AAC,
    sizeof(aac_context_t),
    aac_decode_init,
    NULL,
    aac_decode_close,
    aac_decode_frame,
};

