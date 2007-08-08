/*
 * Common code between AC3 encoder and decoder
 * Copyright (c) 2000 Fabrice Bellard.
 * Copyright (c) 2007 Bartlomiej Wolowiec
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
 * @file ac3.c
 * Common code between AC3 encoder and decoder.
 */
#define DEBUG

#include "avcodec.h"
#include "ac3.h"
#include "bitstream.h"
#include "random.h"

uint8_t bndtab[51];
uint8_t masktab[253];

/** tables for ungrouping mantissas */
float ff_ac3_b1_mantissas[ 32][3];
float ff_ac3_b2_mantissas[128][3];
float ff_ac3_b3_mantissas[8];
float ff_ac3_b4_mantissas[128][2];
float ff_ac3_b5_mantissas[16];

/** dynamic range table. converts codes to scale factors. */
float ff_ac3_dynrng_tbl[256];

/** dialogue normalization table */
float ff_ac3_dialnorm_tbl[32];

/** table for exponent to scale_factor mapping */
float ff_ac3_scale_factors[25];

/** table for grouping exponents */
uint8_t ff_ac3_exp_ungroup_tbl[128][3];

/**
 * Table for number of exponent groups
 * format: nexpgrp_tbl[cplinu][expstr-1][nmant]
 */
uint8_t ff_ac3_nexpgrp_tbl[2][3][256];

static inline int calc_lowcomp1(int a, int b0, int b1, int c)
{
    if ((b0 + 256) == b1) {
        a = c;
    } else if (b0 > b1) {
        a = FFMAX(a - 64, 0);
    }
    return a;
}

static inline int calc_lowcomp(int a, int b0, int b1, int bin)
{
    if (bin < 7) {
        return calc_lowcomp1(a, b0, b1, 384);
    } else if (bin < 20) {
        return calc_lowcomp1(a, b0, b1, 320);
    } else {
        return FFMAX(a - 128, 0);
    }
}

void ff_ac3_bit_alloc_calc_psd(int8_t *exp, int start, int end, int16_t *psd,
                               int16_t *bndpsd)
{
    int bin, i, j, k, end1, v;

    /* exponent mapping to PSD */
    for(bin=start;bin<end;bin++) {
        psd[bin]=(3072 - (exp[bin] << 7));
    }

    /* PSD integration */
    j=start;
    k=masktab[start];
    do {
        v=psd[j];
        j++;
        end1 = FFMIN(bndtab[k+1], end);
        for(i=j;i<end1;i++) {
            /* logadd */
            int adr = FFMIN(FFABS(v - psd[j]) >> 1, 255);
            v = FFMAX(v, psd[j]) + ff_ac3_latab[adr];
            j++;
        }
        bndpsd[k]=v;
        k++;
    } while (end > bndtab[k]);
}

void ff_ac3_bit_alloc_calc_mask(AC3BitAllocParameters *s, int16_t *bndpsd,
                                int start, int end, int fgain, int is_lfe,
                                int deltbae, int deltnseg, uint8_t *deltoffst,
                                uint8_t *deltlen, uint8_t *deltba,
                                int16_t *mask)
{
    int16_t excite[50]; /* excitation */
    int bin, k;
    int bndstrt, bndend, begin, end1, tmp;
    int lowcomp, fastleak, slowleak;

    /* excitation function */
    bndstrt = masktab[start];
    bndend = masktab[end-1] + 1;

    if (bndstrt == 0) {
        lowcomp = 0;
        lowcomp = calc_lowcomp1(lowcomp, bndpsd[0], bndpsd[1], 384);
        excite[0] = bndpsd[0] - fgain - lowcomp;
        lowcomp = calc_lowcomp1(lowcomp, bndpsd[1], bndpsd[2], 384);
        excite[1] = bndpsd[1] - fgain - lowcomp;
        begin = 7;
        for (bin = 2; bin < 7; bin++) {
            if (!(is_lfe && bin == 6))
                lowcomp = calc_lowcomp1(lowcomp, bndpsd[bin], bndpsd[bin+1], 384);
            fastleak = bndpsd[bin] - fgain;
            slowleak = bndpsd[bin] - s->sgain;
            excite[bin] = fastleak - lowcomp;
            if (!(is_lfe && bin == 6)) {
                if (bndpsd[bin] <= bndpsd[bin+1]) {
                    begin = bin + 1;
                    break;
                }
            }
        }

        end1=bndend;
        if (end1 > 22) end1=22;

        for (bin = begin; bin < end1; bin++) {
            if (!(is_lfe && bin == 6))
                lowcomp = calc_lowcomp(lowcomp, bndpsd[bin], bndpsd[bin+1], bin);

            fastleak = FFMAX(fastleak - s->fdecay, bndpsd[bin] - fgain);
            slowleak = FFMAX(slowleak - s->sdecay, bndpsd[bin] - s->sgain);
            excite[bin] = FFMAX(fastleak - lowcomp, slowleak);
        }
        begin = 22;
    } else {
        /* coupling channel */
        begin = bndstrt;

        fastleak = (s->cplfleak << 8) + 768;
        slowleak = (s->cplsleak << 8) + 768;
    }

    for (bin = begin; bin < bndend; bin++) {
        fastleak = FFMAX(fastleak - s->fdecay, bndpsd[bin] - fgain);
        slowleak = FFMAX(slowleak - s->sdecay, bndpsd[bin] - s->sgain);
        excite[bin] = FFMAX(fastleak, slowleak);
    }

    /* compute masking curve */

    for (bin = bndstrt; bin < bndend; bin++) {
        tmp = s->dbknee - bndpsd[bin];
        if (tmp > 0) {
            excite[bin] += tmp >> 2;
        }
        mask[bin] = FFMAX(ff_ac3_hth[bin >> s->halfratecod][s->fscod], excite[bin]);
    }

    /* delta bit allocation */

    if (deltbae == 0 || deltbae == 1) {
        int band, seg, delta;
        band = 0;
        for (seg = 0; seg < deltnseg; seg++) {
            band += deltoffst[seg];
            if (deltba[seg] >= 4) {
                delta = (deltba[seg] - 3) << 7;
            } else {
                delta = (deltba[seg] - 4) << 7;
            }
            for (k = 0; k < deltlen[seg]; k++) {
                mask[band] += delta;
                band++;
            }
        }
    }
}

void ff_ac3_bit_alloc_calc_bap(int16_t *mask, int16_t *psd, int start, int end,
                               int snroffset, int floor, uint8_t *bap)
{
    int i, j, k, end1, v, address;

    /* special case, if snroffset is -960, set all bap's to zero */
    if(snroffset == -960) {
        memset(bap, 0, 256);
        return;
    }
        memset(bap, 0, 256);

    i = start;
    j = masktab[start];
    do {
        v = (FFMAX(mask[j] - snroffset - floor, 0) & 0x1FE0) + floor;
        end1 = FFMIN(bndtab[j] + ff_ac3_bndsz[j], end);
        for (k = i; k < end1; k++) {
            address = av_clip((psd[i] - v) >> 5, 0, 63);
            bap[i] = ff_ac3_baptab[address];
            assert(bap[i]<16);
            i++;
        }
    } while (end > bndtab[j++]);
}

/* AC3 bit allocation. The algorithm is the one described in the AC3
   spec. */
void ac3_parametric_bit_allocation(AC3BitAllocParameters *s, uint8_t *bap,
                                   int8_t *exp, int start, int end,
                                   int snroffset, int fgain, int is_lfe,
                                   int deltbae,int deltnseg,
                                   uint8_t *deltoffst, uint8_t *deltlen,
                                   uint8_t *deltba)
{
    int16_t psd[256];   /* scaled exponents */
    int16_t bndpsd[50]; /* interpolated exponents */
    int16_t mask[50];   /* masking value */

    ff_ac3_bit_alloc_calc_psd(exp, start, end, psd, bndpsd);

    ff_ac3_bit_alloc_calc_mask(s, bndpsd, start, end, fgain, is_lfe,
                               deltbae, deltnseg, deltoffst, deltlen, deltba,
                               mask);

    ff_ac3_bit_alloc_calc_bap(mask, psd, start, end, snroffset, s->floor, bap);
}

/**
 * Initializes some tables.
 * note: This function must remain thread safe because it is called by the
 *       AVParser init code.
 */
void ac3_common_init(void)
{
    int i, j, k, l, v;
    /* compute bndtab and masktab from bandsz */
    k = 0;
    l = 0;
    for(i=0;i<50;i++) {
        bndtab[i] = l;
        v = ff_ac3_bndsz[i];
        for(j=0;j<v;j++) masktab[k++]=i;
        l += v;
    }
    bndtab[50] = l;
}

/**
 * Generates a Kaiser-Bessel Derived Window.
 */
void ff_ac3_window_init(float *window)
{
   int i, j;
   double sum = 0.0, bessel, tmp;
   double local_window[256];
   double alpha2 = (5.0 * M_PI / 256.0) * (5.0 * M_PI / 256.0);

   for(i=0; i<256; i++) {
       tmp = i * (256 - i) * alpha2;
       bessel = 1.0;
       for(j=100; j>0; j--) /* default to 100 iterations */
           bessel = bessel * tmp / (j * j) + 1;
       sum += bessel;
       local_window[i] = sum;
   }

   sum++;
   for(i=0; i<256; i++)
       window[i] = sqrt(local_window[i] / sum);
}

static inline float
symmetric_dequant(int code, int levels)
{
    return (code - (levels >> 1)) * (2.0f / levels);
}

void ff_ac3_decoder_tables_init(void)
{
    int i, expstr, cplinu;

    /* generate grouped mantissa tables
       reference: Section 7.3.5 Ungrouping of Mantissas */
    for(i=0; i<32; i++) {
        /* bap=1 mantissas */
        ff_ac3_b1_mantissas[i][0] = symmetric_dequant( i / 9     , 3);
        ff_ac3_b1_mantissas[i][1] = symmetric_dequant((i % 9) / 3, 3);
        ff_ac3_b1_mantissas[i][2] = symmetric_dequant((i % 9) % 3, 3);
    }
    for(i=0; i<128; i++) {
        /* bap=2 mantissas */
        ff_ac3_b2_mantissas[i][0] = symmetric_dequant( i / 25     , 5);
        ff_ac3_b2_mantissas[i][1] = symmetric_dequant((i % 25) / 5, 5);
        ff_ac3_b2_mantissas[i][2] = symmetric_dequant((i % 25) % 5, 5);

        /* bap=4 mantissas */
        ff_ac3_b4_mantissas[i][0] = symmetric_dequant(i / 11, 11);
        ff_ac3_b4_mantissas[i][1] = symmetric_dequant(i % 11, 11);
    }
    /* generate ungrouped mantissa tables
       reference: Tables 7.21 and 7.23 */
    for(i=0; i<7; i++) {
        /* bap=3 mantissas */
        ff_ac3_b3_mantissas[i] = symmetric_dequant(i, 7);
    }
    for(i=0; i<15; i++) {
        /* bap=5 mantissas */
        ff_ac3_b5_mantissas[i] = symmetric_dequant(i, 15);
    }

    /* generate dynamic range table
       reference: Section 7.7.1 Dynamic Range Control */
    for(i=0; i<256; i++) {
        int v = (i >> 5) - ((i >> 7) << 3) - 5;
        ff_ac3_dynrng_tbl[i] = powf(2.0f, v) * ((i & 0x1F) | 0x20);
    }

    /* generate dialogue normalization table
       references: Section 5.4.2.8 dialnorm
                   Section 7.6 Dialogue Normalization */
    for(i=1; i<32; i++) {
        ff_ac3_dialnorm_tbl[i] = expf((i-31) * M_LN10 / 20.0f);
    }
    ff_ac3_dialnorm_tbl[0] = ff_ac3_dialnorm_tbl[31];

    /* generate scale factors */
    for(i=0; i<25; i++)
        ff_ac3_scale_factors[i] = pow(2.0, -i);

    /* generate exponent tables
       reference: Section 7.1.3 Exponent Decoding */
    for(i=0; i<128; i++) {
        ff_ac3_exp_ungroup_tbl[i][0] =  i / 25;
        ff_ac3_exp_ungroup_tbl[i][1] = (i % 25) / 5;
        ff_ac3_exp_ungroup_tbl[i][2] = (i % 25) % 5;
    }

    /** generate table for number of exponent groups
        reference: Section 7.1.3 Exponent Decoding */
    for(cplinu=0; cplinu<=1; cplinu++) {
        for(expstr=EXP_D15; expstr<=EXP_D45; expstr++) {
            for(i=0; i<256; i++) {
                int grpsize = expstr + (expstr == EXP_D45);
                int ngrps = 0;
                if(cplinu) {
                    ngrps = i / (3 * grpsize);
                } else {
                    if(i == 7)
                        ngrps = 2;
                    else
                        ngrps = (i + (grpsize * 3) - 4) / (3 * grpsize);
                }
                ff_ac3_nexpgrp_tbl[cplinu][expstr-1][i] = ngrps;
            }
        }
    }
}


/** ungrouped mantissas */
typedef struct {
    float b1_mant[3];
    float b2_mant[3];
    float b4_mant[2];
    int b1ptr;
    int b2ptr;
    int b4ptr;
} mant_groups;

/** Gets the transform coefficients for particular channel */
//TODO add @param
static void get_transform_coeffs_ch(
        GetBitContext *gb, uint8_t *bap, uint8_t *dexps, mant_groups *m,
        int dithflag, float *transform_coeffs, int strmant, int endmant
        )
{
//    GetBitContext *gb = &ctx->gb;
    int i, gcode ;

    for(i=strmant; i<endmant; i++) {
       // int tbap = ctx->bap[ch][i];
        int tbap = bap[i];
        assert(tbap>=0 && tbap<=15);
        switch(tbap) {
            case 0:
                if(!dithflag) {
                    transform_coeffs[i] = 0.0f;
                } else {
                    transform_coeffs[i] = LEVEL_MINUS_3DB;
                }
                break;
            case 1:
                if(m->b1ptr > 2) {
                    gcode = get_bits(gb, 5);
                    m->b1_mant[0] = ff_ac3_b1_mantissas[gcode][0];
                    m->b1_mant[1] = ff_ac3_b1_mantissas[gcode][1];
                    m->b1_mant[2] = ff_ac3_b1_mantissas[gcode][2];
                    m->b1ptr = 0;
                }
                transform_coeffs[i] = m->b1_mant[m->b1ptr++];
                break;
            case 2:
                if(m->b2ptr > 2) {
                    gcode = get_bits(gb, 7);
                    m->b2_mant[0] = ff_ac3_b2_mantissas[gcode][0];
                    m->b2_mant[1] = ff_ac3_b2_mantissas[gcode][1];
                    m->b2_mant[2] = ff_ac3_b2_mantissas[gcode][2];
                    m->b2ptr = 0;
                }
                transform_coeffs[i] = m->b2_mant[m->b2ptr++];
                break;
            case 3:
                transform_coeffs[i] = ff_ac3_b3_mantissas[get_bits(gb, 3)];
                break;
            case 4:
                if(m->b4ptr > 1) {
                    gcode = get_bits(gb, 7);
                    m->b4_mant[0] = ff_ac3_b4_mantissas[gcode][0];
                    m->b4_mant[1] = ff_ac3_b4_mantissas[gcode][1];
                    m->b4ptr = 0;
                }
                transform_coeffs[i] = m->b4_mant[m->b4ptr++];
                break;
            case 5:
                transform_coeffs[i] = ff_ac3_b5_mantissas[get_bits(gb, 4)];
                break;
            default:
                /* asymmetric dequantization */
                transform_coeffs[i] = get_sbits(gb, ff_qntztab[tbap]) * ff_ac3_scale_factors[ff_qntztab[tbap]-1];
                break;
        }
        transform_coeffs[i] *= ff_ac3_scale_factors[dexps[i]];
        //av_log(NULL, AV_LOG_INFO, "dexps[%i] = %i transform_coeffs[%i] = %f\n", i, dexps[i], i, transform_coeffs[i]);
    }
}

/**
 * Applies random dithering to coefficients
 * reference: Section 7.3.4 Dither for Zero Bit Mantissas (bap=0)
 */
static void apply_dithering(int nchans, int *endmant, int *dithflag, float (*transform_coeffs)[AC3_MAX_COEFS],int *chincpl, AVRandomState *dith_state, uint8_t (*bap)[AC3_MAX_COEFS]) {
    int ch, i;
    int end=0;
    float *coeffs;

    for(ch=1; ch<=nchans; ch++) {
        coeffs = transform_coeffs[ch];
        if(chincpl[ch])
            end = endmant[CPL_CH];
        else
            end = endmant[ch];
        if(dithflag[ch]) {
            for(i=0; i<end; i++) {
                if(bap[ch][i] == 0) {
                    coeffs[i] *= (av_random(dith_state) & 0xFFFF) / 32768.0f;
                }
            }
        }
    }
}

/**
 * Gets the transform coefficients.
 * This function extracts the tranform coefficients from the AC-3 bitstream
 * using the previously decoded bit allocation pointers.  If coupling is in
 * use, coupled coefficients are also reconstructed here.
 */

void ff_ac3_get_transform_coeffs(GetBitContext *gb, uint8_t (*bap)[AC3_MAX_COEFS], uint8_t (*dexps)[AC3_MAX_COEFS], int nchans, int *chincpl, int *dithflag, float (*transform_coeffs)[AC3_MAX_COEFS], int *strtmant, int *endmant, AVRandomState *dith_state, int ncplbnd, int *cplbndstrc, float (*cplco)[18])
{
    int ch, end;
    int got_cplch = 0;
    mant_groups m;

    m.b1ptr = m.b2ptr = m.b4ptr = 3;

    for(ch=1; ch<=nchans; ch++) {
        /* transform coefficients for individual channel */
        get_transform_coeffs_ch(gb, bap[ch], dexps[ch], &m, dithflag[ch], transform_coeffs[ch], strtmant[ch],endmant[ch]);

        /* tranform coefficients for coupling channel */
        if(chincpl[ch]) {
            if(!got_cplch) {
                int i;
                // TODO
                for(i=0; i<AC3_MAX_COEFS; i++) transform_coeffs[0][i]=0;

                get_transform_coeffs_ch(gb, bap[CPL_CH], dexps[CPL_CH], &m, dithflag[CPL_CH], transform_coeffs[CPL_CH], strtmant[CPL_CH],endmant[CPL_CH]);
               // uncouple_channels(ctx);
                {
                    //TODO (form ac3)
                    int i, j, ch, bnd, subbnd;

                    subbnd = 0;
                    i = strtmant[CPL_CH];
                    /*av_log(NULL, AV_LOG_INFO, "strtmant=%i endmant=%i\n", strtmant[CPL_CH], endmant[CPL_CH]);
                    for(bnd=0; bnd<256; bnd++){
                        av_log(NULL, AV_LOG_INFO, "%i: %f\n", bnd, transform_coeffs[CPL_CH][bnd]);
                    }*/
                    for(bnd=0; bnd<ncplbnd; bnd++) {
                        do {
                            for(j=0; j<12; j++) {
                                for(ch=1; ch<=nchans; ch++) {// TODO lfe?
                                    if(chincpl[ch]) {
                                        transform_coeffs[ch][i] =
                                            transform_coeffs[CPL_CH][i] *
                                            cplco[ch][bnd] * 8.0f;
                                    }
                                }
                                //av_log(NULL, AV_LOG_INFO, "%i ", i);
                                i++;
                            }
                        } while(cplbndstrc[subbnd++]);
                    }
                    //av_log(NULL, AV_LOG_INFO, "\n");
                }
               got_cplch = 1;
            }
            end = endmant[CPL_CH];
        } else {
            end = endmant[ch];
        }
        memset(&transform_coeffs[ch][end], 0, (256-end)*sizeof(float));
    }
    // TODO
//static void apply_dithering(int nfchans, int *endmant, int *dithflag, float (*transform_coeffs)[AC3_MAX_COEFS],int *chincpl, AVRandomState *dith_state, uint8_t *bap) {
    apply_dithering(nchans, endmant, dithflag, transform_coeffs,chincpl, dith_state, bap);
}


/**
 * Decodes the grouped exponents.
 * This function decodes the coded exponents according to exponent strategy
 * and stores them in the decoded exponents buffer.
 */
void ff_ac3_decode_exponents(GetBitContext *gb, int expstr, int ngrps,
                             uint8_t absexp, uint8_t *dexps)
{
    int i, j, grp, grpsize;
    int dexp[256];
    int expacc, prevexp;

    /* unpack groups */
    grpsize = expstr + (expstr == EXP_D45);
    for(grp=0,i=0; grp<ngrps; grp++) {
        expacc = get_bits(gb, 7);
        dexp[i++] = ff_ac3_exp_ungroup_tbl[expacc][0];
        dexp[i++] = ff_ac3_exp_ungroup_tbl[expacc][1];
        dexp[i++] = ff_ac3_exp_ungroup_tbl[expacc][2];
    }

    /* convert to absolute exps and expand groups */
    prevexp = absexp;
    for(i=0; i<ngrps*3; i++) {
        prevexp = av_clip(prevexp + dexp[i]-2, 0, 24);
        for(j=0; j<grpsize; j++) {
            dexps[(i*grpsize)+j] = prevexp;
        }
    }
}


