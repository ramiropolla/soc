/**
 * @file aac.c
 * AAC decoder
 * @author Oded Shimon  ( ods15 ods15 dyndns org )
 * @author Maxim Gavrilov ( maxim.gavrilov gmail com )

 * Kaiser-Bessel Derived Window by Justin Ruggles

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

// Parsing predefines
static void ics_info(aac_context_t * ac, GetBitContext * gb, ics_struct * out);
static void section_data(aac_context_t * ac, GetBitContext * gb, ics_struct * ics, int cb[][64]);
static void scale_factor_data(aac_context_t * ac, GetBitContext * gb, int global_gain, const ics_struct * ics, const int cb[][64], float sf[8][64]);
static void pulse_data(aac_context_t * ac, GetBitContext * gb);
static void tns_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, tns_struct * tns);
static int gain_control_data(aac_context_t * ac, GetBitContext * gb);
static int ms_data(aac_context_t * ac, GetBitContext * gb, ms_struct * ms);
static void spectral_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, const int cb[][64], const float sf[8][64], float * out);
static int individual_channel_stream(aac_context_t * ac, GetBitContext * gb, int common_window, int scale_flag, ics_struct * ics, tns_struct * tns, float * out);
static int single_channel_element(aac_context_t * ac, GetBitContext * gb);
static int channel_pair_element(aac_context_t * ac, GetBitContext * gb);

// Tools predefines
static void ms_tool(aac_context_t * ac, ms_struct * ms);
static void tns_tool(aac_context_t * ac, const ics_struct * ics, const tns_struct * tns, float * coef);
static void window(aac_context_t * ac, ics_struct * ics, float * in, float * out, float * saved);


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
        window[n-1-k] = window[k];
    }
}

/**
* Generate a sine Window.
*/
static void sine_window_init(float *window, int n) {
    const float alpha = M_PI / (2.0 * n);
    int i;
    for(i = 0; i < n; i++)
        window[i] = sin((i + 0.5) * alpha);
}

// General functions
static int GASpecificConfig(aac_context_t * ac, GetBitContext * gb) {
    int ext = 0;
    ac->frame_length = get_bits(gb, 1);
    if (get_bits(gb, 1))
        get_bits(gb, 14);
    ext = get_bits(gb, 1);
    assert(ac->frame_length == 0);
    assert(ext == 0);
    assert(ac->channels != 0); //FIX: program_config_element
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
        if (get_bits(gb, 1)) ;
    }
    return 0;
}

static inline int GetAudioObjectType(GetBitContext * gb) {
    int result = get_bits(gb, 5);
    if (result == 31)
        result = 32 + get_bits(gb, 6);
}

static inline int GetSampleRate(GetBitContext * gb, int *index, int *rate) {
    static const int sampling_table[] = { 96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050, 16000, 12000, 11025, 8000, 7350 };
    *index= get_bits(gb, 4);
    if (*index == 0xf) {
        *index = -1;
        *rate = get_bits(gb, 24);
    } else if (*index <= 12) {
        *rate = sampling_table[*index];
    } else {
        return 1;
    }
    return 0;
}

static int AudioSpecificConfig(aac_context_t * ac, void *data, int data_size) {
    GetBitContext * gb = &ac->gb;

    init_get_bits(gb, data, data_size * 8);

    ac->audioObjectType = GetAudioObjectType(gb);
    assert(ac->audioObjectType == AOT_AAC_LC);
    if (GetSampleRate(gb, &ac->sampling_index, &ac->sample_rate)) return 1;
    ac->channels = get_bits(gb, 4);
    assert(ac->channels == 2);

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
                ac->sbr_present = get_bits(gb, 1);
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
    avccontext->channels = ac->channels;
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

    // windows init
    kbd_window_init(4, ac->kbd_long_1024, 2048, 50);
    kbd_window_init(6, ac->kbd_short_128, 256, 50);
    sine_window_init(ac->sine_long_1024, 1024);
    sine_window_init(ac->sine_short_128, 128);

    // general init
    memset(ac->saved, 0, sizeof(ac->saved));
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
static void ics_info(aac_context_t * ac, GetBitContext * gb, ics_struct * ics) {
    int reserved;
    reserved = get_bits(gb, 1);
    assert(reserved == 0);
    ics->window_sequence = get_bits(gb, 2);
    ics->window_shape = get_bits(gb, 1);
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
        ics->predictor = get_bits(gb, 1);
        assert(ics->predictor == 0);
        ics->swb_offset = ac->swb_offset_1024;
        ics->num_swb = ac->num_swb_1024;
        ics->num_windows = 1;
        ics->tns_max_bands = ac->tns_max_bands_1024;
    }
}

static inline float ivquant(int abs_a) {
    static const float tab[21] = {
        0,
        1,
        2.5198420997897464,
        4.3267487109222245,
        6.3496042078727974,
        8.5498797333834844,
        10.902723556992836,
        13.390518279406722,
        16.000000000000000, // =)
        18.720754407467133,
        21.544346900318832,
        24.463780996262464,
        27.47314182127996,
        30.567350940369842,
        33.741991698453212,
        36.993181114957046,
        40.317473596635935,
        43.711787041189993,
        47.173345095760126,
        50.699631325716943,
        54.288352331898118,
    };
    if (abs_a <= 20) return tab[abs_a];
    else             return pow(abs_a, 4./3);
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
            assert(sect_cb < 12);
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

static void scale_factor_data(aac_context_t * ac, GetBitContext * gb, int global_gain, const ics_struct * ics, const int cb[][64], float sf[8][64]) {
    int g, i;
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            static const float pow2sf_tab[] = {
                2.9802322387695313E-008, 5.9604644775390625E-008, 1.1920928955078125E-007,
                2.384185791015625E-007, 4.76837158203125E-007, 9.5367431640625E-007,
                1.9073486328125E-006, 3.814697265625E-006, 7.62939453125E-006,
                1.52587890625E-005, 3.0517578125E-005, 6.103515625E-005,
                0.0001220703125, 0.000244140625, 0.00048828125,
                0.0009765625, 0.001953125, 0.00390625,
                0.0078125, 0.015625, 0.03125,
                0.0625, 0.125, 0.25,
                0.5, 1.0, 2.0,
                4.0, 8.0, 16.0, 32.0,
                64.0, 128.0, 256.0,
                512.0, 1024.0, 2048.0,
                4096.0, 8192.0, 16384.0,
                32768.0, 65536.0, 131072.0,
                262144.0, 524288.0, 1048576.0,
                2097152.0, 4194304.0, 8388608.0,
                16777216.0, 33554432.0, 67108864.0,
                134217728.0, 268435456.0, 536870912.0,
                1073741824.0, 2147483648.0, 4294967296.0,
                8589934592.0, 17179869184.0, 34359738368.0,
                68719476736.0, 137438953472.0, 274877906944.0
            };
            static const float pow2_table[] = {
                1.0,
                1.1892071150027210667174999705605, /* 2^0.25 */
                1.4142135623730950488016887242097, /* 2^0.5 */
                1.6817928305074290860622509524664 /* 2^0.75 */
            };
            if ((cb[g][i] == ZERO_HCB) || (cb[g][i] == INTENSITY_HCB) || (cb[g][i] == INTENSITY_HCB2)) {
                sf[g][i] = 0.;
            } else {
                global_gain += get_vlc2(gb, ac->mainvlc.table, 7, 3) - 60;
                assert(!(global_gain & (~255)));
                sf[g][i] = pow2sf_tab[global_gain >> 2] * pow2_table[global_gain & 3]/1024./32768.;
            }
        }
    }
}

static void pulse_data(aac_context_t * ac, GetBitContext * gb) {
    int i, number_pulse, pulse_start_sfb;
    // FIXME ignored
    if (get_bits(gb, 1)) { // pulse data
        av_log(ac->avccontext, AV_LOG_INFO, " pulse skipped\n");
        number_pulse = get_bits(gb, 2);
        pulse_start_sfb = get_bits(gb, 6);
        for (i = 0; i < number_pulse + 1; i++) {
            int pulse_offset = get_bits(gb, 5);
            int pulse_amp = get_bits(gb, 4);
        }
    }
}

static void tns_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, tns_struct * tns) {
    int w, filt, i, coef_len, coef_res = 0, coef_compress;
    if (tns->present = get_bits(gb, 1)) {
        for (w = 0; w < ics->num_windows; w++) {
            tns->n_filt[w] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 1 : 2);
            if (tns->n_filt[w])
                coef_res = get_bits(gb, 1) + 3;
            for (filt = 0; filt < tns->n_filt[w]; filt++) {
                tns->length[w][filt] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 4 : 6);
                if (tns->order[w][filt] = get_bits(gb, (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 3 : 5)) {
                    tns->direction[w][filt] = get_bits(gb, 1);
                    assert(coef_res == 3 || coef_res == 4);
                    coef_compress = get_bits(gb, 1);
                    coef_len = coef_res - coef_compress;
                    tns->tmp2_map = tns_tmp2_map[(coef_compress << 1) + (coef_res - 3)];
                    for (i = 0; i < tns->order[w][filt]; i++)
                        tns->coef[w][filt][i] = get_bits(gb, coef_len);
                }
            }
        }
    }
}

static int gain_control_data(aac_context_t * ac, GetBitContext * gb) {
    // FIXME ignored
    if (get_bits(gb, 1)) { // gain control data
    av_log(ac->avccontext, AV_LOG_INFO, " gain data ignored\n");
        return 1;
    }
    return 0;
}

static int ms_data(aac_context_t * ac, GetBitContext * gb, ms_struct * ms) {
    ms->present = get_bits(gb, 2);
    if (ms->present == 1) {
        int g, i;
        for (g = 0; g < ac->ics[0].num_window_groups; g++)
            for (i = 0; i < ac->ics[0].max_sfb; i++)
                ms->mask[g][i] = get_bits(gb, 1);// << i;
    } else if (ms->present == 2) {
        int g, i;
        for (g = 0; g < ac->ics[0].num_window_groups; g++)
            for (i = 0; i < ac->ics[0].max_sfb; i++)
                ms->mask[g][i] = 1;// = 0xFFFFFFFFFFFFFFFFULL;
    }
    return 0;
}

static void spectral_data(aac_context_t * ac, GetBitContext * gb, const ics_struct * ics, const int cb[][64], const float sf[8][64], float * out) {
    static const int unsigned_cb[] = { 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1 };
    int i, k, g;
    const uint16_t * offsets = ics->swb_offset;

    for (g = 0; g < ics->num_window_groups; g++) {
        int total = (ics->window_sequence == EIGHT_SHORT_SEQUENCE) ? 128 : 1024;
        memset(out + g*total + offsets[ics->max_sfb], 0, sizeof(float)*(total - offsets[ics->max_sfb]));
    }
    for (g = 0; g < ics->num_window_groups; g++) {
        for (i = 0; i < ics->max_sfb; i++) {
            const int cur_cb = cb[g][i];
            const int dim = (cur_cb >= FIRST_PAIR_HCB) ? 2 : 4;
            int group;
            if ((cur_cb == ZERO_HCB) || (cur_cb == NOISE_HCB) ||
                (cur_cb == INTENSITY_HCB2) || (cur_cb == INTENSITY_HCB)) {
                int j;
                for (j = 0; j < ics->group_len[g]; j++) {
                    memset(out + j * 128 + offsets[i], 0, (offsets[i+1] - offsets[i])*sizeof(float));
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
                            if (ptr[j] && get_bits(gb, 1))
                                sign[j] = -1;
                    }
                    if (cur_cb == 11) {
                        for (j = 0; j < 2; j++) {
                            if (ptr[j] == 16) {
                                int N = 4;
                                while (get_bits(gb, 1)) N++;
                                ptr[j] = (1<<N) + get_bits(gb, N);
                            }
                        }
                    }
                    for (j = 0; j < dim; j++)
                        out[group*128+j+k] = sign[j] * ivquant(ptr[j]) * sf[g][i];
                    //for (j = 0; j < dim; j++) av_log(ac->avccontext, AV_LOG_INFO, " %4d: %5d %10.3lf => %10.3lf\n", j+k, ptr[j]*sign[j], sf[g][i]*1024*32768, out[group*128+j+k]*1024*32768);
                }
                //av_log(ac->avccontext, AV_LOG_INFO, " checking escape %d[%d] %d\n", ptr[j], j, index);
                assert(k == offsets[i+1]);
            }
        }
        out += ics->group_len[g]*128;
    }
}

static int individual_channel_stream(aac_context_t * ac, GetBitContext * gb, int common_window, int scale_flag, ics_struct * ics, tns_struct * tns, float * out) {
    int global_gain;
    int cb[8][64];   // codebooks
    float sf[8][64]; // scale factors

    memset(sf, 0, sizeof(sf));

    global_gain = get_bits(gb, 8);

    if (!common_window && !scale_flag) {
        ics_info(ac, gb, ics);
    }

    //av_log(ac->avccontext, AV_LOG_INFO, " global_gain: %d, groups: %d\n", global_gain, ics->window_sequence);
    section_data(ac, gb, ics, cb);
    scale_factor_data(ac, gb, global_gain, ics, cb, sf);

    if (!scale_flag) {
        pulse_data(ac, gb);
        tns_data(ac, gb, ics, tns);
        if (gain_control_data(ac, gb))
            return 1;
    }

    //av_log(ac->avccontext, AV_LOG_INFO, " read: %d ->", get_bits_count(gb));
    spectral_data(ac, gb, ics, cb, sf, out);
    //av_log(ac->avccontext, AV_LOG_INFO, " done: %d\n", get_bits_count(gb));
    return 0;
}

static int single_channel_element(aac_context_t * ac, GetBitContext * gb) {
    int id;
    id = get_bits(gb, 4);
    return individual_channel_stream(ac, gb, 0, 0, &ac->ics[0], &ac->tns[0], ac->coeffs[0]);
}

static int channel_pair_element(aac_context_t * ac, GetBitContext * gb) {
    int common_window;
    ms_struct ms = {0};
    int id;

    id = get_bits(gb, 4); // element instance tag
    assert(id == 0);
    common_window = get_bits(gb, 1);
    if (common_window) {
        ics_info(ac, gb, &ac->ics[0]);
        ac->ics[1] = ac->ics[0];
        ms_data(ac, gb, &ms);
    } else {
        ms.present = 0;
    }
    if (individual_channel_stream(ac, gb, common_window, 0, &ac->ics[0], &ac->tns[0], ac->coeffs[0]))
        return 1;
    if (individual_channel_stream(ac, gb, common_window, 0, &ac->ics[1], &ac->tns[1], ac->coeffs[1]))
        return 1;

    // M/S tool
    if (common_window) {
        ms_tool(ac, &ms);
    }
    return 0;
}

// Tools implementation
static void ms_tool(aac_context_t * ac, ms_struct * ms) {
    if (ms->present) {
        int g, i, k, start = 0, gp;
        const uint16_t * offsets = ac->ics[0].swb_offset;
        for (g = 0; g < ac->ics[0].num_window_groups; g++) {
            //av_log(ac->avccontext, AV_LOG_INFO, " masking[%d]: ", g);
            for (gp = start; gp < start + ac->ics[0].group_len[g]; gp++) {
                for (i = 0; i < ac->ics[0].max_sfb; i++)
                    if (ms->mask[g][i]) {
                        for (k = offsets[i] + gp*128; k < offsets[i+1] + gp*128; k++) {
                            float tmp = ac->coeffs[0][k] - ac->coeffs[1][k];
                            ac->coeffs[0][k] += ac->coeffs[1][k];
                            ac->coeffs[1][k] = tmp;
                        }
                    }
            }
            start += ac->ics[0].group_len[g];
            //av_log(ac->avccontext, AV_LOG_INFO, "\n");
        }
    }
}

static void tns_tool(aac_context_t * ac, const ics_struct * ics, const tns_struct * tns, float * coef) {
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
                b[ib] = b[ib + order] = coef[start];
                coef[start] = tmp;
                start += inc;
            }
        }
    }
}

#define BIAS 385

static void window(aac_context_t * ac, ics_struct * ics, float * in, float * out, float * saved) {
    const float * lwindow = (ics->window_shape) ? ac->kbd_long_1024 : ac->sine_long_1024;
    const float * swindow = (ics->window_shape) ? ac->kbd_short_128 : ac->sine_long_1024;
    float * buf = ac->buf_mdct;
    if (ics->window_sequence != EIGHT_SHORT_SEQUENCE) {
        int i;
        ff_imdct_calc(&ac->mdct, buf, in, out); // out can be abused for now as a temp buffer
        if (ac->is_saved) {
            if (ics->window_sequence != LONG_STOP_SEQUENCE) {
                for (i = 0; i < 1024; i++)        out[i] = saved[i] + buf[i] * lwindow[i           ] + BIAS;
            } else {
                for (i = 0; i < 512-64; i++)      out[i] = saved[i] +                                  BIAS;
                for (i = 512-64; i < 512+64; i++) out[i] = saved[i] + buf[i] * swindow[i - (512-64)] + BIAS;
                for (i = 512+64; i < 1024; i++)   out[i] =            buf[i]                         + BIAS;
            }
        }
        if (ics->window_sequence != LONG_START_SEQUENCE) {
            for (i = 0; i < 1024; i++)        saved[i] = buf[i+1024] * lwindow[1023-i];
        } else {
            for (i = 0; i < 512-64; i++)      saved[i] = buf[i+1024];
            for (i = 512-64; i < 512+64; i++) saved[i] = buf[i+1024] * swindow[127-(i - (512-64))];
        }
    } else {
        float * ptr[8];
        int i;
        for (i = 0; i < 8; i++) {
            int k;
            ptr[i] = buf + i*256;
            ff_imdct_calc(&ac->mdct_small, ptr[i], in + i*128, out);
            for (k = 0; k < 128; k++)   ptr[i][k] *= swindow[k] * 8.;
            for (k = 128; k < 256; k++) ptr[i][k] *= swindow[255-k] * 8.;
        }
        for (i = 0; i < 512-64; i++) out[i] = saved[i]                           + BIAS;
        for (; i < 576; i++)         out[i] = saved[i]         + ptr[0][i - 448] + BIAS;
        for (; i < 704; i++)         out[i] = ptr[0][i - 448]  + ptr[1][i - 576] + BIAS;
        for (; i < 832; i++)         out[i] = ptr[1][i - 576]  + ptr[2][i - 704] + BIAS;
        for (; i < 960; i++)         out[i] = ptr[2][i - 704]  + ptr[3][i - 832] + BIAS;
        for (; i < 1024; i++)        out[i] = ptr[3][i - 832]  + ptr[4][i - 960] + BIAS;

        for (; i < 1088; i++) saved[i-1024] = ptr[3][i - 832]  + ptr[4][i - 960];
        for (; i < 1216; i++) saved[i-1024] = ptr[4][i - 960]  + ptr[5][i - 1088];
        for (; i < 1344; i++) saved[i-1024] = ptr[5][i - 1088] + ptr[6][i - 1216];
        for (; i < 1472; i++) saved[i-1024] = ptr[6][i - 1216] + ptr[7][i - 1344];
        for (; i < 1600; i++) saved[i-1024] = ptr[7][i - 1344]; // memcpy?
    }
}

static int aac_decode_frame(AVCodecContext * avccontext, void * data, int * data_size, uint8_t * buf, int buf_size) {
    aac_context_t * ac = avccontext->priv_data;
    GetBitContext * gb = &ac->gb;
    int id;

    //ac->num_frame++;
    //if (ac->num_frame == 7)
    //    __asm int 3;

    init_get_bits(gb, buf, buf_size*8);

    // parse
    while ((id = get_bits(gb, 3)) != 7) {
        switch (id) {
            case ID_SCE: {
                single_channel_element(ac, gb);
                break;
            }
            case ID_CPE: {
                channel_pair_element(ac, gb);
                break;
            }
            case ID_FIL: {
                int cnt = get_bits(gb, 4);
                if (cnt == 15) cnt += get_bits(gb, 8) - 1;
                skip_bits(gb, 8*cnt);
                break;
            }
            case ID_CCE:
            case ID_LFE:
            case ID_DSE:
            case ID_PCE:
                assert(0);
                break;
            default:
                assert(0 && 0);
                break;
        }
        break;
    }
    //av_log(avccontext, AV_LOG_INFO, " %d %d %d\n", get_bits_count(gb), buf_size*8 - get_bits_count(gb), id);
    // Processing
    tns_tool(ac, &ac->ics[0], &ac->tns[0], ac->coeffs[0]);
    tns_tool(ac, &ac->ics[1], &ac->tns[1], ac->coeffs[1]);

    for (id = 0; id < 2; id++) {
        int i;
        window(ac, &ac->ics[id], ac->coeffs[id], ac->ret, ac->saved[id]);
        if (ac->is_saved)
            for (i = 0; i < 1024; i++) {
                int32_t tmp = ((int32_t*)ac->ret)[i];
                if (tmp & 0xf0000) {
                    if (tmp > 0x43c0ffff) tmp = 0xFFFF;
                    else                  tmp = 0;
                }
                ((uint16_t(*)[2])data)[i][id] = tmp - 0x8000;
        }
    }
    *data_size = ac->is_saved ? 4096 : 0;
    ac->is_saved = 1;

    return buf_size;
}

static int aac_decode_close(AVCodecContext * avccontext) {
    aac_context_t * ac = avccontext->priv_data;
    int i;

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

