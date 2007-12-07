/*
 * AAC (LC) decoder
 * This code is developed as part of Google Summer of Code 2006 Program.
 *
 * Copyright (c) 2005 Oded Shimon( ods15 ods15 dyndns org )
 * Copyright (c) 2005-2006 Maxim Gavrilov ( maxim.gavrilov gmail com )
 * Copyright (c) 2007 Andreas Öman ( andreas lonelycoder com)
 *
 * This file is part of FFmpeg.
 *
 * FFmpeg is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * FFmpeg is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with FFmpeg; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 */

/*
 * The 'reference:'-statements refers to the ISO/IEC 14496-3 specification
 */

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"
#include "random.h"
#include "aac.h"
#include "aactab.h"

#define MAX_TAGID 16

/**
 * Program config
 */
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
} aac_program_config;

/**
 * Individual Channel Stream
 */
typedef struct {
    int intensity_present;
    int noise_present;

    int max_sfb;
    int window_sequence;
    int window_sequence_prev;
    int window_shape;
    int window_shape_prev;
    int predictor;
    int num_window_groups;
    uint8_t grouping;
    uint8_t group_len[8];
    // calculated
    const uint16_t *swb_offset;
    int num_swb;
    int num_windows;
    int tns_max_bands;
} aac_ics;

#define ICS_ESS(ics) ((ics)->window_sequence == EIGHT_SHORT_SEQUENCE)


/**
 * Temporal Noise Shaping
 */
typedef struct {
    int present;
    int n_filt[8];
    int length[8][4];
    int direction[8][4];
    int order[8][4];
    const float *tmp2_map[8];
    int coef[8][4][TNS_MAX_ORDER];
} aac_tns;


/**
 * M/S tool
 */
typedef struct {
    int present;
    int mask[8][64];
} aac_ms;


/**
 * Pulse tool
 */
typedef struct {
    int present;
    int num_pulse;
    int start;
    int offset[4];
    int amp[4];
} aac_pulse;


/**
 * Coupling struct
 */
typedef struct {
    int ind_sw;
    int domain;

    int num_coupled;
    int is_cpe[9];
    int tag_select[9];
    int l[9];
    int r[9];

    float gain[18][8][64];
} aac_coupling;


/**
 * Single Channel Element
 */
typedef struct {
    int global_gain;
    aac_ics ics;
    aac_tns tns;
    int cb[8][64];
    float sf[8][64];
    DECLARE_ALIGNED_16(float, coeffs[1024]);
    DECLARE_ALIGNED_16(float, saved[1024]);
    DECLARE_ALIGNED_16(float, ret[1024]);
} aac_sce;


/**
 * Channel Pair Element
 */
typedef struct {
    int common_window;
    aac_ms ms;
    aac_sce ch[2];
} aac_cpe;


/**
 * Channel Coupling Element
 */
typedef struct {
    aac_coupling coup;
    aac_sce ch;
} aac_cce;


/**
 * AAC context
 */
typedef struct {
    AVCodecContext *avctx;
    GetBitContext gb;
    VLC mainvlc;
    VLC books[11];

    // main config
    int audioObjectType;
    int ext_audioObjectType;
    int sbr_present;
    int sample_rate;
    int sample_rate_index;
    int channels;
    int ochannels;
    int frame_length;

    // decoder param
    aac_program_config pcs;
    aac_sce *che_sce[MAX_TAGID];
    aac_cpe *che_cpe[MAX_TAGID];
    aac_sce *che_lfe[MAX_TAGID];
    aac_cce *che_cce[MAX_TAGID];

    DECLARE_ALIGNED_16(float, buf_mdct[2048]);
    int is_saved;

    const uint16_t *swb_offset_1024;
    const uint16_t *swb_offset_128;
    int num_swb_1024;
    int num_swb_128;
    int tns_max_bands_1024;
    int tns_max_bands_128;

    DECLARE_ALIGNED_16(float, kbd_long_1024[1024]);
    DECLARE_ALIGNED_16(float, kbd_short_128[128]);
    DECLARE_ALIGNED_16(float, sine_long_1024[1024]);
    DECLARE_ALIGNED_16(float, sine_short_128[128]);
    DECLARE_ALIGNED_16(float, pow2sf_tab[256]);
    DECLARE_ALIGNED_16(float, intensity_tab[256]);
    DECLARE_ALIGNED_16(float, ivquant_tab[256]);
    DECLARE_ALIGNED_16(float, revers[1024]);
    float* iop;

    MDCTContext mdct;
    MDCTContext mdct_small;
    DSPContext dsp;
    int * vq[11];
    AVRandomState random_state;

    int add_bias;
    int scale_bias;

    int num_frame;
} AACContext;


static int is_intensity(int cb)
{
    if(cb == INTENSITY_HCB)
        return 1;
    if(cb  == INTENSITY_HCB2)
        return -1;
    return 0;
}

#define is_noise(cb) (cb == NOISE_HCB)

static void vector_fmul_add_add_add(AACContext *ac, float *dst,
                                    const float *src0, const float *src1,
                                    const float *src2, const float *src3,
                                    float src4, int len) {
    int i;
    ac->dsp.vector_fmul_add_add(dst, src0, src1, src2, src4, len, 1);
    for(i = 0; i < len; i++)
        dst[i] += src3[i];
}


#define TAG_MASK 0x00f
#define FLAG_SCE 0x100
#define FLAG_CPE 0x200
#define FLAG_LFE 0x400
#define FLAG_CCE 0x800

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


/**
 * Parse audio object type
 * reference: Table 1.14
 */
static int GetAudioObjectType(GetBitContext *gb)
{
    int result = get_bits(gb, 5);
    if(result == 31)
        result = 32 + get_bits(gb, 6);
    return result;
}

/**
 * Parse program config element
 * reference: Table 4.2
 *
 * XXX: Needs fixup
 */
static int program_config_element(AACContext *ac, GetBitContext *gb)
{
    aac_program_config *pcs = &ac->pcs;
    int id, object_type, i;

    pcs->present = 1;
    id = get_bits(gb, 4);
    object_type = get_bits(gb, 2);

    ac->sample_rate_index = get_bits(gb, 4);

    pcs->num_front      = get_bits(gb, 4);
    pcs->num_side       = get_bits(gb, 4);
    pcs->num_back       = get_bits(gb, 4);
    pcs->num_lfe        = get_bits(gb, 2);
    pcs->num_assoc_data = get_bits(gb, 3);
    pcs->num_cc         = get_bits(gb, 4);

    pcs->mono_mixdown   = get_bits1(gb) ? get_bits(gb, 4) + 1 : 0;
    pcs->stereo_mixdown = get_bits1(gb) ? get_bits(gb, 4) + 1 : 0;

    if(get_bits1(gb)) {
        pcs->matrix_mixdown = get_bits(gb, 2) + 1;
        pcs->pseudo_surround = get_bits1(gb);
    } else {
        pcs->matrix_mixdown = 0;
        pcs->pseudo_surround = 0;
    }

    pcs->front_cpe = 0;
    ac->channels += pcs->num_front;
    for(i = 0; i < pcs->num_front; i++) {
        if(get_bits1(gb)) {
            pcs->front_cpe |= (1 << i);
            ac->channels++;
        }
        pcs->front_tag[i] = get_bits(gb, 4);
    }
    pcs->side_cpe = 0;
    ac->channels += pcs->num_side;
    for(i = 0; i < pcs->num_side; i++) {
        if(get_bits1(gb)) {
            pcs->side_cpe |= (1 << i);
            ac->channels++;
        }
        pcs->side_tag[i] = get_bits(gb, 4);
    }
    pcs->back_cpe = 0;
    ac->channels += pcs->num_back;
    for(i = 0; i < pcs->num_back; i++) {
        if(get_bits1(gb)) {
            pcs->back_cpe |= (1 << i);
            ac->channels++;
        }
        pcs->back_tag[i] = get_bits(gb, 4);
    }
    ac->channels += pcs->num_lfe;
    for(i = 0; i < pcs->num_lfe; i++)
        pcs->lfe_tag[i] = get_bits(gb, 4);

    pcs->num_channels = ac->channels;
    // not a real audio channel
    for(i = 0; i < pcs->num_assoc_data; i++)
        pcs->assoc_data_tag[i] = get_bits(gb, 4);
    pcs->cc_ind_sw = 0;
    for(i = 0; i < pcs->num_cc; i++) {
        pcs->cc_ind_sw |= (get_bits1(gb) << i);
        pcs->cc_tag[i] = get_bits(gb, 4);
    }
    align_get_bits(gb);
    skip_bits(gb, 8 * get_bits(gb, 8));
    return 0;
}

/**
 * Implicit channel configuration
 * reference: Table 1.17
 *
 * XXX: Needs fixup
 */
static int implicit_channel_config(AACContext *ac)
{
    switch(ac->channels) {
    case 1:        /* C */
        ac->che_sce[0] = av_mallocz(sizeof(aac_sce));
        ac->ochannels = 1;
        break;

    case 2:        /* L + R */
        ac->che_cpe[0] = av_mallocz(sizeof(aac_cpe));
        ac->ochannels = 2;
        break;


    case 6:
        ac->che_sce[0] = av_mallocz(sizeof(aac_sce));
        ac->che_cpe[0] = av_mallocz(sizeof(aac_cpe));
        ac->che_cpe[1] = av_mallocz(sizeof(aac_cpe));
        ac->che_lfe[0] = av_mallocz(sizeof(aac_sce));
        ac->ochannels = 6;
        break;

    default:
        return -1;
    }
    return 0;
}

/**
 * Parse GA specific configuration
 * reference: Table 4.1
 */
static int GASpecificConfig(AACContext *ac, GetBitContext *gb)
{
    if(get_bits1(gb))     // frameLengthFlag
        return -1;        // We only support 1024 / 128*8 MDCT windows

    if(get_bits1(gb))
        get_bits(gb, 14); // coreCoderDelay

    if(get_bits1(gb))     // extensionFlag, should be 0 for AAC_LC
        return -1;

    if(ac->channels == 0)
        return program_config_element(ac, gb);
    else
        return implicit_channel_config(ac);
}


/**
 * Parse audio specific configuration
 * reference: Table 1.13
 */
static int AudioSpecificConfig(AACContext *ac, void *data, int data_size)
{
    GetBitContext * gb = &ac->gb;

    init_get_bits(gb, data, data_size * 8);

    ac->audioObjectType = GetAudioObjectType(gb);
    if(ac->audioObjectType != AOT_AAC_LC)
        return -1;

    ac->sample_rate_index = get_bits(gb, 4);
    if(ac->sample_rate_index == 0xf)
        ac->sample_rate = get_bits(gb, 24);

    ac->channels = get_bits(gb, 4);

    if(GASpecificConfig(ac, gb) < 0)
        return -1;

    if(ac->sample_rate_index == 15) {
        /* Explicit rate configured, XXX: Use table 4.68 */
    } else if(ac->sample_rate_index > 12) {
        return -1; // Reserved
    } else {
        ac->sample_rate = aac_sample_rates[ac->sample_rate_index];
    }
    return 0;
}

/**
 * Top level init
 */
static int aac_decode_init(AVCodecContext *avctx)
{
    AACContext *ac = avctx->priv_data;
    int i;

    ac->avctx = avctx;

    if(AudioSpecificConfig(ac, avctx->extradata, avctx->extradata_size))
        return -1;

    avctx->sample_rate = ac->sample_rate;
    avctx->channels = ac->ochannels;

    for(i = 0; i < 11; i++) {
        static const int mod_cb[11] = { 3, 3, 3, 3, 9, 9, 8, 8, 13, 13, 17 };
        static const int off_cb[11] = { 1, 1, 0, 0, 4, 4, 0, 0,  0,  0,  0 };
        const aac_codebook *aco = &aac_codebooks[i];

        int j, values = aco->s/sizeof(aco->a[0]);
        int dim = (i >= 4 ? 2 : 4);
        int mod = mod_cb[i], off = off_cb[i], index = 0;
        int ret;
        ret = init_vlc(&ac->books[i], 6, values,
                       &aco->a[0][1], sizeof(aco->a[0]),
                       sizeof(aco->a[0][1]),
                       &aco->a[0][0], sizeof(aco->a[0]),
                       sizeof(aco->a[0][0]),
                       0);
        assert(!ret);
        ac->vq[i] = av_malloc(dim * values * sizeof(int));
        if(dim == 2) {
            for(j = 0; j < values * dim; j += dim) {
                index = j/dim;
                ac->vq[i][j  ] = (index / (mod            ) - off);
                index %= mod;
                ac->vq[i][j+1] = (index                     - off);
            }
        } else {
            for(j = 0; j < values * dim; j += dim) {
                index = j/dim;
                ac->vq[i][j  ] = (index / (mod * mod * mod) - off);
                index %= mod*mod*mod;
                ac->vq[i][j+1] = (index / (mod * mod      ) - off);
                index %= mod*mod;
                ac->vq[i][j+2] = (index / (mod            ) - off);
                index %= mod;
                ac->vq[i][j+3] = (index                     - off);
            }
        }
    }

    /* Initialize RNG dither */
    av_init_random(0x1f2e3d4c, &ac->random_state);

    /* Speedup tables */
    for(i = 0; i < 256; i++)
        ac->intensity_tab[i] = pow(0.5, (i - 100) / 4.);
    for(i = 0; i < sizeof(ac->ivquant_tab)/sizeof(ac->ivquant_tab[0]); i++)
        ac->ivquant_tab[i] = pow(i, 4./3);

    if(ac->dsp.float_to_int16 == ff_float_to_int16_c) {
        ac->add_bias = 385;
        ac->scale_bias = 32768;
    } else {
        ac->add_bias = 0;
        ac->scale_bias = 1;
    }
    for(i = 0; i < 256; i++)
        ac->pow2sf_tab[i] = pow(2, (i - 100)/4.) /1024./ac->scale_bias;

    ac->num_frame = -1;

    ac->swb_offset_1024    = aac_swb_offset_1024[ac->sample_rate_index];
    ac->num_swb_1024       = aac_num_swb_1024[ac->sample_rate_index];
    ac->tns_max_bands_1024 = aac_tns_max_bands_1024[ac->sample_rate_index];
    ac->swb_offset_128     = aac_swb_offset_128[ac->sample_rate_index];
    ac->num_swb_128        = aac_num_swb_128[ac->sample_rate_index];
    ac->tns_max_bands_128  = aac_tns_max_bands_128[ac->sample_rate_index];

    init_vlc(&ac->mainvlc, 7,
             sizeof(aac_scalefactor_huffman_table)/
             sizeof(aac_scalefactor_huffman_table[0]),
             &aac_scalefactor_huffman_table[0][1],
             sizeof(aac_scalefactor_huffman_table[0]),
             sizeof(aac_scalefactor_huffman_table[0][1]),
             &aac_scalefactor_huffman_table[0][0],
             sizeof(aac_scalefactor_huffman_table[0]),
             sizeof(aac_scalefactor_huffman_table[0][0]),
             0);

    ff_mdct_init(&ac->mdct, 11, 1);
    ff_mdct_init(&ac->mdct_small, 8, 1);
    kbd_window_init(4, ac->kbd_long_1024, 2048, 50);
    kbd_window_init(6, ac->kbd_short_128, 256,  50);
    sine_window_init(ac->sine_long_1024, 2048);
    sine_window_init(ac->sine_short_128, 256);

    for(i = 0; i < 128; i++) {
        ac->sine_short_128[i] *= 8.;
        ac->kbd_short_128[i] *= 8.;
    }
    dsputil_init(&ac->dsp, avctx);
    return 0;
}

/**
 * Joint coding, M/S Stereo parser
 */
static int ms_parser(AACContext *ac, GetBitContext *gb, aac_cpe *cpe)
{
    aac_ms *ms = &cpe->ms;
    int g, i;

    ms->present = get_bits(gb, 2);
    switch(ms->present) {
    case 1:
    case 2:
        for(g = 0; g < cpe->ch[0].ics.num_window_groups; g++) {
            for(i = 0; i < cpe->ch[0].ics.max_sfb; i++) {
                ms->mask[g][i] = ms->present == 1 ? get_bits1(gb) : 1;
            }
        }
        break;

    case 3:
        return -1; /* Reserved */
    }
    return 0;
}

/**
 * Joint coding, M/S Stereo
 */
static void ms_tool(AACContext *ac, aac_cpe *cpe) {
    const aac_ms *ms = &cpe->ms;
    const aac_ics *ics = &cpe->ch[0].ics;
    aac_sce *sce0 = &cpe->ch[0];
    aac_sce *sce1 = &cpe->ch[1];
    float *ch0 = sce0->coeffs;
    float *ch1 = sce1->coeffs;
    int g, i, k, gp;
    const uint16_t * offsets = ics->swb_offset;

    if(!ms->present)
        return;

    for(g = 0; g < ics->num_window_groups; g++) {
        for(gp = 0; gp < ics->group_len[g]; gp++) {
            for(i = 0; i < ics->max_sfb; i++) {
                if(ms->mask[g][i] && !is_intensity(sce1->cb[g][i]) &&
                    !is_noise(sce1->cb[g][i])) {
                    for(k = offsets[i]; k < offsets[i+1]; k++) {
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


/**
 * Pulse tool parser
 */
static void pulse_parser(AACContext *ac, GetBitContext *gb, aac_pulse *pulse)
{
    int i;

    if(!(pulse->present = get_bits1(gb)))
       return;
    pulse->num_pulse = get_bits(gb, 2);
    pulse->start = get_bits(gb, 6);
    for(i = 0; i <= pulse->num_pulse; i++) {
        pulse->offset[i] = get_bits(gb, 5);
        pulse->amp[i] = get_bits(gb, 4);
    }
}


/**
 * Pulse tool
 */
static void pulse_tool(AACContext *ac, aac_ics * const ics,
                       aac_pulse * const pulse, int *icoef)
{
    int i, off;
    if(!pulse->present)
        return;

    off = ics->swb_offset[pulse->start];
    for(i = 0; i <= pulse->num_pulse; i++) {
        off += pulse->offset[i];
        if(icoef[off] > 0)
            icoef[off] += pulse->amp[i];
        else
            icoef[off] -= pulse->amp[i];
    }
}


/**
 * Intensity stereo tool
 */
static void intensity_tool(AACContext *ac, aac_cpe *cpe)
{
    const aac_ics *ics = &cpe->ch[1].ics;
    aac_sce *sce0 = &cpe->ch[0];
    aac_sce *sce1 = &cpe->ch[1];
    const uint16_t *offsets = ics->swb_offset;
    int g, gp, i, k, c;
    float scale;

    if(!ics->intensity_present)
        return;

    for(g = 0; g < ics->num_window_groups; g++) {
        for(gp = 0; gp < ics->group_len[g]; gp++) {
            for(i = 0; i < ics->max_sfb; i++) {
                if((c = is_intensity(sce1->cb[g][i]))) {
                    if(cpe->ms.present == 1)
                        c *= (1 - 2 * cpe->ms.mask[g][i]);

                    scale = c * sce1->sf[g][i];

                    for(k = offsets[i] + gp * 128;
                         k < offsets[i+1] + gp * 128; k++) {
                        sce1->coeffs[k] = scale * sce0->coeffs[k];
                    }
                }
            }
        }
    }
}


/**
 * Temporal Noise Shaping tool parser
 */
static void tns_parser(AACContext *ac, GetBitContext *gb, const aac_ics *ics,
                       aac_tns *tns)
{
    int w, filt, i, coef_len, coef_res = 0, coef_compress;

    if(!(tns->present = get_bits1(gb)))
        return;

    for(w = 0; w < ics->num_windows; w++) {
        tns->n_filt[w] = get_bits(gb, ICS_ESS(ics) ? 1 : 2);
        if(tns->n_filt[w])
            coef_res = get_bits1(gb) + 3;
        for(filt = 0; filt < tns->n_filt[w]; filt++) {
            tns->length[w][filt] = get_bits(gb, ICS_ESS(ics) ? 4 : 6);
            if((tns->order[w][filt] = get_bits(gb, ICS_ESS(ics) ? 3 : 5))) {
                tns->direction[w][filt] = get_bits1(gb);
                coef_compress = get_bits1(gb);
                coef_len = coef_res - coef_compress;
                tns->tmp2_map[w] = aac_tns_coeffs_table[(coef_compress << 1) +
                                                        (coef_res - 3)];
                for(i = 0; i < tns->order[w][filt]; i++)
                    tns->coef[w][filt][i] = get_bits(gb, coef_len);
            }
        }
    }
}


/**
 * Temporal Noise Shaping tool
 */
static void tns_filter_tool(AACContext * ac, aac_sce *sce, float *coef)
{
    const aac_ics *ics = &sce->ics;
    const aac_tns *tns = &sce->tns;
    const int mmm = FFMIN(ics->tns_max_bands,  ics->max_sfb);
    int w, filt, m, i, ib, bottom, top, order, start, end, size, inc;
    float tmp;
    float lpc[TNS_MAX_ORDER + 1], b[2 * TNS_MAX_ORDER];
    if(!tns->present)
        return;

    for(w = 0; w < ics->num_windows; w++) {
        bottom = ics->num_swb;
        for(filt = 0; filt < tns->n_filt[w]; filt++) {
            top = bottom;
            bottom = FFMAX(top - tns->length[w][filt], 0);
            order  = FFMIN(tns->order[w][filt], TNS_MAX_ORDER);
            if(order == 0)
                continue;

            // tns_decode_coef
            lpc[0] = 1;
            for(m = 1; m <= order; m++) {
                lpc[m] = tns->tmp2_map[w][tns->coef[w][filt][m - 1]];
                for(i = 1; i < m; i++)
                    b[i] = lpc[i] + lpc[m] * lpc[m-i];
                for(i = 1; i < m; i++)
                    lpc[i] = b[i];
            }

            start = ics->swb_offset[FFMIN(bottom, mmm)];
            end = ics->swb_offset[FFMIN(top, mmm)];
            if((size = end - start) <= 0)
                continue;
            if(tns->direction[w][filt]) {
                inc = -1; start = end - 1;
            } else {
                inc = 1;
            }
            start += w * 128;

            // ar filter
            memset(b, 0, sizeof(b));
            ib = 0;

            for(m = 0; m < size; m++) {
                tmp = coef[start];
                for(i = 0; i < order; i++)
                    tmp -= b[ib + i] * lpc[i + 1];
                if(--ib < 0)
                    ib = order - 1;
                b[ib] = b[ib + order] = tmp;
                coef[start] = tmp;
                start += inc;
            }
        }
    }
}


/**
 * Wrapper for use together with transform_sce_tool()
 */
static void tns_trans(AACContext *ac, aac_sce *sce)
{
    tns_filter_tool(ac, sce, sce->coeffs);
}


/**
 *
 */
static float ivquant(AACContext *ac, int a) {
    const float sign = FFSIGN(a);
    a = FFABS(a);
    if(a < sizeof(ac->ivquant_tab) / sizeof(ac->ivquant_tab[0]))
        return sign * ac->ivquant_tab[a];
    else
        return sign * pow(a, 4./3);
}


/**
 * Quantization tool
 */
static void quant_to_spec_tool(AACContext *ac, aac_ics *ics, const int *icoef,
                               int cb[][64], float sf[][64], float *coef)
{
    const uint16_t * offsets = ics->swb_offset;
    int g, i, group, k;
    int total = ICS_ESS(ics) ? 128 : 1024;
    float energy, scale;
    for(g = 0; g < ics->num_window_groups; g++) {
        memset(coef + g*total + offsets[ics->max_sfb], 0,
               sizeof(float) * (total - offsets[ics->max_sfb]));
    }

    for(g = 0; g < ics->num_window_groups; g++) {
        for(i = 0; i < ics->max_sfb; i++) {
            if(cb[g][i] == NOISE_HCB) {
                for(group = 0; group < ics->group_len[g]; group++) {
                    energy = 0;
                    scale = 1;
                    for(k = offsets[i]; k < offsets[i+1]; k++)
                        energy += (float)icoef[group*128+k] *
                            icoef[group*128+k];

                    scale *= sf[g][i] / sqrt(energy);
                    for(k = offsets[i]; k < offsets[i+1]; k++)
                        coef[group*128+k] = icoef[group*128+k] * scale;
                }
            } else if(cb[g][i] != INTENSITY_HCB &&
                      cb[g][i] != INTENSITY_HCB2) {

                for(group = 0; group < ics->group_len[g]; group++) {
                    for(k = offsets[i]; k < offsets[i+1]; k++) {
                        coef[group*128+k] = ivquant(ac, icoef[group*128+k]) *
                            sf[g][i];
                    }
                }
            } else {
                assert(0);
            }
        }
        coef  += ics->group_len[g]*128;
        icoef += ics->group_len[g]*128;
    }
}


/**
 * Inverse MDCT Filterbank tool
 *
 * The imdct transform uses 'out' as a temporary buffer
 */
static void window_trans(AACContext *ac, aac_sce *sce) {
    aac_ics *ics = &sce->ics;
    float *in = sce->coeffs;
    float *out = sce->ret;
    float *saved = sce->saved;

    const float *lwindow = ics->window_shape ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float *swindow = ics->window_shape ? ac->kbd_short_128 : ac->sine_short_128;
    const float *lwindow_prev = ics->window_shape_prev ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float *swindow_prev = ics->window_shape_prev ? ac->kbd_short_128 : ac->sine_short_128;
    float * buf = ac->buf_mdct;
    int i;

    if(ics->window_sequence != EIGHT_SHORT_SEQUENCE) {
        ff_imdct_calc(&ac->mdct, buf, in, out);

        if(ac->is_saved) {
            if(ics->window_sequence != LONG_STOP_SEQUENCE) {
                ac->dsp.vector_fmul_add_add(out, buf, lwindow_prev, saved,
                                            ac->add_bias, 1024, 1);
            } else {
                for(i = 0; i < 448; i++) out[i] = saved[i] + ac->add_bias;
                for(i = 448; i < 576; i++) buf[i] *= 0.125; // normalize
                ac->dsp.vector_fmul_add_add(out + 448, buf + 448, swindow_prev,
                                            saved + 448, ac->add_bias, 128, 1);
                for(i = 576; i < 1024; i++)   out[i] = buf[i] + ac->add_bias;
            }
        }
        if(ics->window_sequence != LONG_START_SEQUENCE) {
            ac->dsp.vector_fmul_reverse(saved, buf + 1024, lwindow, 1024);
        } else {
            memcpy(saved, buf + 1024, 448 * sizeof(float));
            for(i = 448; i < 576; i++) buf[i + 1024] *= 0.125; // normalize
            ac->dsp.vector_fmul_reverse(saved + 448, buf + 1024 + 448, swindow,
                                        128);
            memset(saved + 576, 0, 448 * sizeof(float));
        }
    } else {
        int i;

        for(i = 0; i < 2048; i += 256) {
            ff_imdct_calc(&ac->mdct_small, buf + i, in + i/2, out);
            ac->dsp.vector_fmul_reverse(ac->revers + i/2, buf + i + 128,
                                        swindow, 128);
        }

        for(i = 0; i < 448; i++)   out[i] = saved[i] + ac->add_bias;
        out += 448; saved += 448;
        ac->dsp.vector_fmul_add_add(out + 0*128, buf + 0*128, swindow_prev, saved, ac->add_bias, 128, 1);
        vector_fmul_add_add_add(ac, out + 1*128, buf + 2*128, swindow, saved + 1*128, ac->revers + 0*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 2*128, buf + 4*128, swindow, saved + 2*128, ac->revers + 1*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 3*128, buf + 6*128, swindow, saved + 3*128, ac->revers + 2*128, ac->add_bias, 128);
        vector_fmul_add_add_add(ac, out + 4*128, buf + 8*128, swindow, saved + 4*128, ac->revers + 3*128, ac->add_bias, 64);

        saved -= 448;
        buf += 1024;
        ac->dsp.vector_fmul_add_add(saved,       buf + 64, swindow, ac->revers + 3*128+64,  0, 64, 1);
        ac->dsp.vector_fmul_add_add(saved + 64,  buf + 2*128, swindow, ac->revers + 4*128, 0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 192, buf + 4*128, swindow, ac->revers + 5*128, 0, 128, 1);
        ac->dsp.vector_fmul_add_add(saved + 320, buf + 6*128, swindow, ac->revers + 6*128, 0, 128, 1);
        memcpy(                     saved + 448, ac->revers + 7*128, 128 * sizeof(float));

        memset(saved + 576, 0, 448 * sizeof(float));
    }
}


/**
 * Decode Individual Channel Stream info
 * reference: table 4.6
 */
static int ics_info(AACContext *ac, GetBitContext *gb, int common_window,
                     aac_ics *ics)
{
    int i;

    if(get_bits1(gb))   // ics_reserved_bit
        return -1;
    ics->window_sequence = get_bits(gb, 2);
    ics->window_shape_prev = ics->window_shape;
    ics->window_shape = get_bits1(gb);
    if(ics->window_shape_prev == -1)
        ics->window_shape_prev = ics->window_shape;
    ics->num_window_groups = 1;
    ics->group_len[0] = 1;
    if(ICS_ESS(ics)) {
        ics->max_sfb = get_bits(gb, 4);
        ics->grouping = get_bits(gb, 7);
        for(i = 0; i < 7; i++) {
            if(ics->grouping & (1<<(6-i))) {
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
    } else {
        ics->max_sfb = get_bits(gb, 6);
        ics->swb_offset = ac->swb_offset_1024;
        ics->num_swb = ac->num_swb_1024;
        ics->num_windows = 1;
        ics->tns_max_bands = ac->tns_max_bands_1024;
        if(get_bits1(gb)) { // predictor_data_present
            av_log(ac->avctx, AV_LOG_ERROR, "Prediction not supported\n");
            return -1;
        }
    }
    return 0;
}


/**
 * Decode scale_factor_data
 * reference: Table 4.47
 */
static int scale_factor_data(AACContext *ac, GetBitContext *gb,
                             aac_sce *sce, aac_ics *ics)
{
    int g, i;
    unsigned int intensity = 100;
    unsigned int global_gain = sce->global_gain;
    int noise = sce->global_gain - 90;
    int noise_flag = 1;
    float sf;

    ics->intensity_present = 0;
    ics->noise_present = 0;
    for(g = 0; g < ics->num_window_groups; g++) {
        for(i = 0; i < ics->max_sfb; i++) {
            switch(sce->cb[g][i]) {
            case ZERO_HCB:
                sce->sf[g][i] = 0;
                continue;

            case INTENSITY_HCB:
            case INTENSITY_HCB2:
                ics->intensity_present = 1;
                intensity += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                if(intensity > 255)
                    return -1;
                sf = ac->intensity_tab[intensity];
                break;

            case NOISE_HCB:
                ics->noise_present = 1;
                if(noise_flag) {
                    noise_flag = 0;
                    noise += get_bits(gb, 9) - 256;
                } else {
                    noise += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                }
                sf = pow(2.0, 0.25 * noise)/1024./ac->scale_bias;
                break;

            default:
                global_gain += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                if(global_gain > 255 || global_gain < 0)
                    return -1;
                sf = ac->pow2sf_tab[global_gain];
            }
            sce->sf[g][i] = sf;
        }
    }
    return 0;
}


/**
 * Decode section_data payload
 * reference: Table 4.46
 */
static void section_data(AACContext *ac, GetBitContext *gb, aac_ics *ics,
                         int cb[][64]) {
    int g, k, sect_len, sect_len_incr, sect_cb;
    const int bits = ICS_ESS(ics) ? 3 : 5;
    const int sect_esc_val = (1 << bits) - 1;

    for(g = 0; g < ics->num_window_groups; g++) {
        for(k = 0; k < ics->max_sfb;) {
            sect_len = 0;
            sect_cb = get_bits(gb, 4);

            while((sect_len_incr = get_bits(gb, bits)) == sect_esc_val)
                sect_len += sect_esc_val;
            sect_len += sect_len_incr;
            sect_len += k;
            for(; k < sect_len && k < ics->max_sfb; k++)
                cb[g][k] = sect_cb;
        }
    }
}


/**
 * Decode spectral data for VLC codebooks
 */
static int spectral_data_vlc(AACContext *ac, GetBitContext *gb,
                             const int cur_cb, const int group,
                             const int k, int *icoef, int d)
{
    static const int unsigned_cb[] = { 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1 };
    int sign[4] = {1,1,1,1};
    int j, ptr[4], index;

    if((index = get_vlc2(gb, ac->books[cur_cb - 1].table, 6, 3)) < 0)
        return -1;

    memcpy(ptr, &ac->vq[cur_cb - 1][index * d], d * sizeof(int));

    if(unsigned_cb[cur_cb - 1]) {
        for(j = 0; j < d; j++)
            if(ptr[j] && get_bits1(gb))
                sign[j] = -1;
    }
    if(cur_cb == ESC_HCB) {
        for(j = 0; j < 2; j++) {
            if(ptr[j] == 16) {
                int n = 4;
                while(get_bits1(gb)) n++;
                ptr[j] = (1<<n) + get_bits(gb, n);
            }
        }
    }
    for(j = 0; j < d; j++)
        icoef[group*128+k+j] = sign[j] * ptr[j];
    return 0;
}


/**
 * Decode spectral data
 * reference: Table 4.50
 */
static int spectral_data(AACContext *ac, GetBitContext *gb,
                         const aac_ics *ics, int cb[8][64],
                         int *icoef)
{
    int i, k, g, cur_cb, d, group;
    const uint16_t * offsets = ics->swb_offset;
    const int total = ICS_ESS(ics) ? 128 : 1024;

    for(g = 0; g < ics->num_window_groups; g++)
        memset(icoef + g * total + offsets[ics->max_sfb], 0,
               sizeof(int) * (total - offsets[ics->max_sfb]));

    for(g = 0; g < ics->num_window_groups; g++) {
        for(i = 0; i < ics->max_sfb; i++) {
            cur_cb = cb[g][i];
            d = cur_cb >= FIRST_PAIR_HCB ? 2 : 4;


            switch(cur_cb) {
            case INTENSITY_HCB:
            case INTENSITY_HCB2:
                continue;

            case NOISE_HCB:
                for(group = 0; group < ics->group_len[g]; group++)
                    for(k = offsets[i]; k < offsets[i+1]; k++)
                        icoef[group*128+k] =
                            av_random(&ac->random_state) & 0x0000FFFF;
                continue;

            case ZERO_HCB:
                for(group = 0; group < ics->group_len[g]; group++)
                    memset(icoef + group * 128 + offsets[i], 0,
                           (offsets[i+1] - offsets[i])*sizeof(int));
                continue;

            default:
                for(group = 0; group < ics->group_len[g]; group++)
                    for(k = offsets[i]; k < offsets[i+1]; k += d)
                        if(spectral_data_vlc(ac, gb, cur_cb, group, k, icoef, d))
                            return -1;
            }
        }
        icoef += ics->group_len[g]*128;
    }
    return 0;
}


/**
 * Decode an individual_channel_stream payload
 * reference: Table 4.44
 */
static int individual_channel_stream(AACContext *ac, GetBitContext *gb,
                                     int common_window, int scale_flag,
                                     aac_sce *sce)
{
    int icoeffs[1024];
    aac_pulse pulse;
    aac_tns *tns = &sce->tns;
    aac_ics *ics = &sce->ics;
    float *out = sce->coeffs;

    memset(&pulse, 0, sizeof(pulse));
    sce->global_gain = get_bits(gb, 8);

    if(!common_window && !scale_flag)
        ics_info(ac, gb, 0, ics);

    section_data(ac, gb, ics, sce->cb);
    if(scale_factor_data(ac, gb, sce, ics))
        return -1;

    if(!scale_flag) {
        pulse_parser(ac, gb, &pulse);
        tns_parser(ac, gb, ics, tns);

        if(get_bits1(gb)) {
            av_log(ac->avctx, AV_LOG_INFO, "gain control not supported\n");
            return -1;
        }
    }
    if(spectral_data(ac, gb, ics, sce->cb, icoeffs))
        return -1;
    pulse_tool(ac, ics, &pulse, icoeffs);
    quant_to_spec_tool(ac, ics, icoeffs, sce->cb, sce->sf, out);
    return 0;
}


/**
 * Decode a single_channel_element
 * reference: Table 4.4
 */
static int single_channel_element(AACContext *ac, GetBitContext *gb)
{
    aac_sce *sce;
    int id = get_bits(gb, 4);
    if((sce = ac->che_sce[id]) == NULL) {
        av_log(ac->avctx, AV_LOG_ERROR,
               "Single channel element %d not configured\n", id);
        return -1;
    }
    return individual_channel_stream(ac, gb, 0, 0, sce);
}


/**
 * Decode a channel_pair_element
 * reference: Table 4.4
 */
static int channel_pair_element(AACContext *ac, GetBitContext *gb) {
    int i;
    aac_cpe *cpe;
    int id = get_bits(gb, 4);
    if((cpe = ac->che_cpe[id]) == NULL) {
        av_log(ac->avctx, AV_LOG_ERROR,
               "Channel pair element %d not configured\n", id);
        return -1;
    }

    cpe->common_window = get_bits1(gb);
    if(cpe->common_window) {
        ics_info(ac, gb, 1, &cpe->ch[0].ics);
        i = cpe->ch[1].ics.window_shape_prev;
        cpe->ch[1].ics = cpe->ch[0].ics;
        cpe->ch[1].ics.window_shape_prev = i;
        if(ms_parser(ac, gb, cpe))
            return -1;
    } else {
        cpe->ms.present = 0;
    }

    if(individual_channel_stream(ac, gb, cpe->common_window, 0, &cpe->ch[0]))
        return -1;
    if(individual_channel_stream(ac, gb, cpe->common_window, 0, &cpe->ch[1]))
        return -1;

    if(cpe->common_window)
        ms_tool(ac, cpe);

    intensity_tool(ac, cpe);
    return 0;
}


/**
 * Decode a lfe_channel_element
 * reference: Table 4.4
 */
static int lfe_channel_element(AACContext *ac, GetBitContext *gb)
{
    aac_sce *sce;
    int id = get_bits(gb, 4);
    if((sce = ac->che_lfe[id]) == NULL) {
        av_log(ac->avctx, AV_LOG_ERROR,
               "Low frequency element %d not configured\n", id);
        return -1;
    }

    return individual_channel_stream(ac, gb, 0, 0, sce);
}


/**
 * Decode a data_stream_element
 * reference: Table 4.10
 */
static int data_stream_element(AACContext *ac, GetBitContext *gb)
{
    int id, byte_align;
    int count;
    id = get_bits(gb, 4);
    byte_align = get_bits1(gb);
    count = get_bits(gb, 8);
    if(count == 255)
        count += get_bits(gb, 8);
    if(byte_align)
        align_get_bits(gb);
    skip_bits(gb, 8 * count);
    return 0;
}


/**
 * Perform function \p f on all active channel elements
 */
static void transform_sce_tool(AACContext * ac,
                               void (*f)(AACContext *ac, aac_sce *sce))
{
    int i;

    for(i = 0; i < MAX_TAGID; i++) {
        if(ac->che_sce[i] != NULL)
            f(ac, ac->che_sce[i]);
        if(ac->che_cpe[i] != NULL) {
            f(ac, &ac->che_cpe[i]->ch[0]);
            f(ac, &ac->che_cpe[i]->ch[1]);
        }
        if(ac->che_lfe[i] != NULL)
            f(ac, ac->che_lfe[i]);
        if(ac->che_cce[i] != NULL)
            f(ac, &ac->che_cce[i]->ch);
    }
}


/**
 * Write out samples. XXX: This needs to be rewritten
 */
static void output_samples(AACContext *ac, AVCodecContext *avctx,
                           void *data, int *data_size)
{
    int size = ac->ochannels * 1024 * sizeof(uint16_t);
    float *mixbuf = NULL, *src = NULL;
    int i;

    *data_size = 0;

    if(!ac->is_saved) {
        ac->is_saved = 1;
        return;
    }

    switch(ac->ochannels) {
    case 1:
        src = ac->che_sce[0]->ret;
        break;

    case 2:
        src = mixbuf = av_malloc(size * 2);
        for(i = 0; i < 1024; i++) {
            mixbuf[i * 2 + 0] = ac->che_cpe[0]->ch[0].ret[i];
            mixbuf[i * 2 + 1] = ac->che_cpe[0]->ch[1].ret[i];
        }
        break;

    case 6:
        src = mixbuf = av_malloc(size * 6);
        for(i = 0; i < 1024; i++) {
            mixbuf[i * 6 + 0] = ac->che_cpe[0]->ch[0].ret[i];
            mixbuf[i * 6 + 1] = ac->che_cpe[0]->ch[1].ret[i];
            mixbuf[i * 6 + 2] = ac->che_sce[0]->ret[i];
            mixbuf[i * 6 + 3] = ac->che_lfe[0]->ret[i];
            mixbuf[i * 6 + 4] = ac->che_cpe[1]->ch[0].ret[i];
            mixbuf[i * 6 + 5] = ac->che_cpe[1]->ch[1].ret[i];
        }
        break;
    }


    ac->dsp.float_to_int16(data, src, 1024 * ac->ochannels);
    av_free(mixbuf);
    *data_size = size;
}


/**
 * Decode an AAC frame
 */
static int aac_decode_frame(AVCodecContext *avctx, void *data, int *data_size,
                            uint8_t * buf, int buf_size)
{
    AACContext * ac = avctx->priv_data;
    GetBitContext * gb = &ac->gb;
    int id, cnt, err;

    init_get_bits(gb, buf, buf_size * 8);

    while((id = get_bits(gb, 3)) != ID_END) {
        err = 0;
        switch (id) {
        case ID_SCE:
            err = single_channel_element(ac, gb);
            break;

        case ID_CPE:
            err = channel_pair_element(ac, gb);
            break;

        case ID_FIL:
            cnt = get_bits(gb, 4);
            if(cnt == 15)
                cnt += get_bits(gb, 8) - 1;
            skip_bits(gb, 8 * cnt);
            break;

        case ID_DSE:
            err = data_stream_element(ac, gb);
            break;

        case ID_LFE:
            err = lfe_channel_element(ac, gb);
            break;

        default:
            av_log(avctx, AV_LOG_ERROR, "Unhandled id %d\n", id);
            err = 1;
            break;
        }
        if(err)
            return buf_size;
    }

    transform_sce_tool(ac, tns_trans);
    transform_sce_tool(ac, window_trans);

    ac->num_frame++;
    output_samples(ac, avctx, data, data_size);
    return buf_size;
}


/**
 * Decoder close
 */
static int aac_decode_close(AVCodecContext *avctx)
{
    AACContext *ac = avctx->priv_data;
    int i;

    for(i = 0; i < MAX_TAGID; i++) {
        av_free(ac->che_sce[i]);
        av_free(ac->che_cpe[i]);
        av_free(ac->che_lfe[i]);
        av_free(ac->che_cce[i]);
    }

    for(i = 0; i < 11; i++) {
        free_vlc(&ac->books[i]);
        av_free(ac->vq[i]);
    }
    free_vlc(&ac->mainvlc);
    ff_mdct_end(&ac->mdct);
    ff_mdct_end(&ac->mdct_small);
    return 0;
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
