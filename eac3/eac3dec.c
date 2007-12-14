/*
 * EAC3 decoder
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

#include "avcodec.h"
#include "eac3.h"
#include "ac3dec.h"
#include "ac3.h"

/**
 * Table for default stereo downmixing coefficients
 * reference: Section 7.8.2 Downmixing Into Two Channels
 */
static const uint8_t eac3_default_coeffs[8][5][2] = {
    { { 2, 7 }, { 7, 2 },                               },
    { { 4, 4 },                                         },
    { { 2, 7 }, { 7, 2 },                               },
    { { 2, 7 }, { 5, 5 }, { 7, 2 },                     },
    { { 2, 7 }, { 7, 2 }, { 6, 6 },                     },
    { { 2, 7 }, { 5, 5 }, { 7, 2 }, { 8, 8 },           },
    { { 2, 7 }, { 7, 2 }, { 6, 7 }, { 7, 6 },           },
    { { 2, 7 }, { 5, 5 }, { 7, 2 }, { 6, 7 }, { 7, 6 }, },
};

static const float mixlevels[9] = {
    LEVEL_PLUS_3DB,
    LEVEL_PLUS_1POINT5DB,
    LEVEL_ONE,
    LEVEL_MINUS_1POINT5DB,
    LEVEL_MINUS_3DB,
    LEVEL_MINUS_4POINT5DB,
    LEVEL_MINUS_6DB,
    LEVEL_ZERO,
    LEVEL_MINUS_9DB
};

static float idct_cos_tab[6][5];

static void log_missing_feature(AVCodecContext *avctx, const char *log){
    av_log(avctx, AV_LOG_ERROR, "%s is not implemented. If you want to help, "
            "update your FFmpeg version to the newest one from SVN. If the "
            "problem still occurs, it means that your file has extension "
            "which has not been tested due to a lack of samples exhibiting "
            "this feature. Upload a sample of the audio from this file to "
            "ftp://upload.mplayerhq.hu/incoming and contact the ffmpeg-devel "
            "mailing list.\n", log);
}

static void uncouple_channels(EAC3Context *s){
    int i, j, ch, bnd, subbnd;

    subbnd = s->cplbegf+1;
    i = s->strtmant[CPL_CH];
    for (bnd = 0; bnd < s->ncplbnd; bnd++) {
        do {
            for (j = 0; j < 12; j++) {
                for (ch = 1; ch <= s->nfchans; ch++) {
                    if (s->chincpl[ch]) {
                        s->transform_coeffs[ch][i] =
                            s->transform_coeffs[CPL_CH][i] *
                            s->cplco[ch][bnd] * 8.0f;
                    }
                }
                i++;
            }
        } while(s->cplbndstrc[subbnd++] && subbnd<=s->cplendf);
    }
}

static void spectral_extension(EAC3Context *s){
    //Now turned off, because there are no samples for testing it.
#if 0
    int copystartmant, copyendmant, copyindex, insertindex;
    int wrapflag[18];
    int bandsize, bnd, bin, spxmant, filtbin, ch;
    float nratio, accum, nscale, sscale, spxcotemp;
    float noffset[AC3_MAX_CHANNELS], nblendfact[AC3_MAX_CHANNELS][18], sblendfact[AC3_MAX_CHANNELS][18];
    float rmsenergy[AC3_MAX_CHANNELS][18];

    //XXX spxbandtable[bnd] = 25 + 12 * bnd ?

    copystartmant = spxbandtable[s->spxstrtf];
    copyendmant = spxbandtable[s->spxbegf];

    for (ch = 1; ch <= s->nfchans; ch++) {
        if (!s->chinspx[ch])
            continue;

        copyindex = copystartmant;
        insertindex = copyendmant;

        for (bnd = 0; bnd < s->nspxbnds; bnd++){
            bandsize = s->spxbndsztab[bnd];
            if ((copyindex + bandsize) > copyendmant) {
                copyindex = copystartmant;
                wrapflag[bnd] = 1;
            } else
                wrapflag[bnd] = 0;
            for (bin = 0; bin < bandsize; bin++){
                if (copyindex == copyendmant)
                    copyindex = copystartmant;
                s->transform_coeffs[ch][insertindex++] = s->transform_coeffs[ch][copyindex++];
            }
        }

        noffset[ch] = s->spxblnd[ch] / 32.0;
        spxmant = spxbandtable[s->spxbegf];
        if (s->spxcoe[ch]) {
            for (bnd = 0; bnd < s->nspxbnds; bnd++){
                bandsize = s->spxbndsztab[bnd];
                nratio = ((spxmant + 0.5*bandsize) / spxbandtable[s->spxendf]) - noffset[ch];
                nratio = FFMAX(FFMIN(nratio, 1.0), 0.0);
                nblendfact[ch][bnd] = sqrt(nratio);
                sblendfact[ch][bnd] = sqrt(1 - nratio);
                spxmant += bandsize;
            }
        }

        spxmant = spxbandtable[s->spxbegf];
        for (bnd = 0; bnd < s->nspxbnds; bnd++){
            bandsize = s->spxbndsztab[bnd];
            accum = 0;
            for (bin = 0; bin < bandsize; bin++){
                accum += (s->transform_coeffs[ch][spxmant] * s->transform_coeffs[ch][spxmant]);
                spxmant++;
            }
            rmsenergy[ch][bnd] = sqrt(accum / bandsize);
        }

        if (s->chinspxatten[ch]) {
            /* apply notch filter at baseband / extension region border */
            filtbin = spxbandtable[s->spxbegf] - 2;
            for (bin = 0; bin < 3; bin++){
                s->transform_coeffs[ch][filtbin] *= ff_eac3_spxattentab[s->spxattencod[ch]][bin];
                filtbin++;
            }
            for (bin = 1; bin >= 0; bin--){
                s->transform_coeffs[ch][filtbin] *= ff_eac3_spxattentab[s->spxattencod[ch]][bin];
                filtbin++;
            }
            filtbin += s->spxbndsztab[0];
            /* apply notch at all other wrap points */
            for (bnd = 1; bnd < s->nspxbnds; bnd++){
                if (wrapflag[bnd]) {
                    filtbin -= 5;
                    for (bin = 0; bin < 3; bin++){
                        s->transform_coeffs[ch][filtbin] *= ff_eac3_spxattentab[s->spxattencod[ch]][bin];
                        filtbin++;
                    }
                    for (bin = 1; bin >= 0; bin--){
                        s->transform_coeffs[ch][filtbin] *= ff_eac3_spxattentab[s->spxattencod[ch]][bin];
                        filtbin++;
                    }
                }
                filtbin += s->spxbndsztab[bnd];
            }
        }

        spxmant = spxbandtable[s->spxbegf];
        for (bnd = 0; bnd < s->nspxbnds; bnd++){
            nscale = rmsenergy[ch][bnd] * nblendfact[ch][bnd];
            sscale = sblendfact[ch][bnd];
            for (bin = 0; bin < s->spxbndsztab[bnd]; bin++){
                //TODO generate noise()
                s->transform_coeffs[ch][spxmant] =
                    s->transform_coeffs[ch][spxmant] * sscale + noise() * nscale;
                spxmant++;
            }
        }

        spxmant = spxbandtable[s->spxbegf];
        for (bnd = 0; bnd < s->nspxbnds; bnd++){
            spxcotemp = s->spxco[ch][bnd];
            for (bin = 0; bin < s->spxbndsztab[bnd]; bin++){
                s->transform_coeffs[ch][spxmant] *= spxcotemp * 32;
                spxmant++;
            }
        }
    }
#endif
}

static void get_transform_coeffs_aht_ch(GetBitContext *gbc, EAC3Context *s, int ch){
    int endbap, bin, n, m;
    int bg, g, bits, pre_chmant, remap, chgaqsections, chgaqmod;
    float mant;

    chgaqmod = get_bits(gbc, 2);

    endbap = chgaqmod<2?12:17;

    chgaqsections = 0;
    for (bin = 0; bin < s->endmant[ch]; bin++) {
        if (s->hebap[ch][bin] > 7 && s->hebap[ch][bin] < endbap)
            chgaqsections++;
    }

    if (chgaqmod == EAC3_GAQ_12 || chgaqmod == EAC3_GAQ_14) {
        for (n = 0; n < chgaqsections; n++) {
            s->chgaqgain[n] = get_bits1(gbc);
        }
    } else {
        if (chgaqmod == EAC3_GAQ_124) {
            int grpgain;
            chgaqsections = (chgaqsections+2)/3;
            for (n = 0; n < chgaqsections; n++) {
                grpgain = get_bits(gbc, 5);
                s->chgaqgain[3*n]   = grpgain/9;
                s->chgaqgain[3*n+1] = (grpgain%9)/3;
                s->chgaqgain[3*n+2] = grpgain%3;
            }
        }
    }

    m=0;
    for (bin = s->strtmant[ch]; bin < s->endmant[ch]; bin++) {
        if (s->hebap[ch][bin] > 7) {
            // GAQ (E3.3.4.2)
            // XXX what about gaqmod = 0 ?
            // difference between Gk=1 and gaqmod=0 ?
            if (s->hebap[ch][bin] < endbap) {
                // hebap in active range
                // Gk = 1<<bg
                bg = ff_gaq_gk[chgaqmod][s->chgaqgain[m++]];
            } else {
                bg = 0;
            }
            bits = ff_bits_vs_hebap[s->hebap[ch][bin]];

            for (n = 0; n < 6; n++) {
                // pre_chmant[n][ch][bin]
                pre_chmant = get_sbits(gbc, bits-bg);
                if (bg && pre_chmant == -(1 << (bits - bg - 1))) {
                    // large mantissa
                    pre_chmant = get_sbits(gbc, bits - ((bg==1)?1:0));
                    if (bg == 1)
                        //Gk = 2
                        mant = (float)pre_chmant/((1<<(bits-1))-1);
                    else
                        //Gk = 4
                        mant = (float)pre_chmant*3.0f/((1<<(bits+1))-2);

                    g = 0;
                    remap = 1;
                } else {
                    // small mantissa
                    if (bg)
                        //Gk = 2 or 4
                        mant = (float)pre_chmant/((1<<(bits-1))-1);
                    else
                        //Gk = 1
                        mant = (float)pre_chmant*2.0f/((1<<bits)-1); ///XXX

                    g = bg;
                    remap = (!bg) && (s->hebap[ch][bin] < endbap);
                }

                //TODO when remap needed ?
                if (remap) {
                    mant = (float)
                        (ff_eac3_gaq_remap[s->hebap[ch][bin]-8][0][g][0]/32768.0f + 1.0f)
                        * mant / (1<<g) +
                        (ff_eac3_gaq_remap[s->hebap[ch][bin]-8][mant<0][g][1]) / 32768.0f;
                }
                s->pre_chmant[n][ch][bin] = mant;
            }
        } else {
            // hebap = 0 or VQ
            if (s->hebap[ch][bin]) {
                pre_chmant = get_bits(gbc, ff_bits_vs_hebap[s->hebap[ch][bin]]);
                for (n = 0; n < 6; n++) {
                    s->pre_chmant[n][ch][bin] =
                        ff_vq_hebap[s->hebap[ch][bin]][pre_chmant][n] / 32768.0f;
                }
            } else {
                for (n = 0; n < 6; n++) {
                    s->pre_chmant[n][ch][bin] = 0;
                }
            }
        }
    }
}

static void idct_transform_coeffs_ch(EAC3Context *s, int ch, int blk){
    // TODO fast IDCT
    int bin, i;
    float tmp;
    for (bin = s->strtmant[ch]; bin < s->endmant[ch]; bin++) {
        tmp = s->pre_chmant[0][ch][bin];
        for (i = 1; i < 6; i++) {
            tmp += idct_cos_tab[blk][i-1] * s->pre_chmant[i][ch][bin];
        }
        s->transform_coeffs[ch][bin] = tmp * ff_ac3_scale_factors[s->dexps[ch][bin]];
    }
}

static void get_eac3_transform_coeffs_ch(GetBitContext *gbc, EAC3Context *s, int blk,
        int ch, mant_groups *m){
    if (s->chahtinu[ch] == 0) {
        ff_ac3_get_transform_coeffs_ch(m, gbc, s->dexps[ch], s->bap[ch],
                s->transform_coeffs[ch], s->strtmant[ch], s->endmant[ch],
                &s->dith_state);
    } else {
        if (s->chahtinu[ch] == 1) {
            get_transform_coeffs_aht_ch(gbc, s, ch);
            s->chahtinu[ch] = -1; /* AHT info for this frame has been read - do not read again */
        }
    }
    if (s->chahtinu[ch] != 0) {
        idct_transform_coeffs_ch(s, ch, blk);
    }

    memset(s->transform_coeffs[ch]+s->endmant[ch], 0,
           sizeof(s->transform_coeffs[ch]) -
           s->endmant[ch] * sizeof(*s->transform_coeffs[ch]));
}

static int parse_bsi(GetBitContext *gbc, EAC3Context *s){
    int i, blk;

    s->stream_type = get_bits(gbc, 2);
    if (s->stream_type == EAC3_STREAM_TYPE_DEPENDENT) {
        log_missing_feature(s->avctx, "Dependent substream");
        return -1;
    } else if (s->stream_type == EAC3_STREAM_TYPE_RESERVED) {
        av_log(s->avctx, AV_LOG_ERROR, "Reserved stream type\n");
        return -1;
    }
    s->substreamid = get_bits(gbc, 3);
    s->frame_size = (get_bits(gbc, 11) + 1) * 2;
    s->sr_code = get_bits(gbc, 2);
    if (s->sr_code == EAC3_SR_CODE_REDUCED) {
        /* The E-AC3 specification does not tell how to handle reduced sample
           rates in bit allocation.  The best assumption would be that it is
           handled like AC3 DolbyNet, but we cannot be sure until we have a
           sample which utilizes this feature. */
        log_missing_feature(s->avctx, "Reduced Sampling Rates");
        return -1;
#if 0
        s->sr_code2 = get_bits(gbc, 2);
        s->num_blocks = 6;
#endif
    } else {
        s->num_blocks = ff_eac3_blocks[get_bits(gbc, 2)];
    }
    s->channel_mode = get_bits(gbc, 3);
    s->lfe_on = get_bits1(gbc);

    // calculate number of channels
    s->nfchans = ff_ac3_channels_tab[s->channel_mode];
    s->num_channels = s->nfchans;
    s->lfe_channel = s->num_channels+1;
    if (s->lfe_on) {
        s->strtmant[s->lfe_channel] = 0;
        s->endmant [s->lfe_channel] = 7;
        s->nchgrps [s->lfe_channel] = 2;
        s->chincpl [s->lfe_channel] = 0;
        s->num_channels++;
    }

    s->bitstream_id = get_bits(gbc, 5);
    if (s->bitstream_id < 11 || s->bitstream_id > 16) {
        av_log(s->avctx, AV_LOG_ERROR, "bitstream id is not within 11 and 16\n");
        return -1;
    }

    for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
        s->dialog_norm[i] = ff_ac3_dialog_norm_tab[get_bits(gbc, 5)];
        if (get_bits1(gbc)) {
            skip_bits(gbc, 8); //skip Compression gain word
        }
    }
#if 0
    /* TODO: Add support for dependent streams */
    if (s->stream_type == EAC3_STREAM_TYPE_DEPENDENT) {
        if (get_bits1(gbc)) {
            skip_bits(gbc, 16); // skip custom channel map
        } else {
            //TODO default channel map based on acmod and lfeon
        }
    }
#endif

    /* set stereo downmixing coefficients
       reference: Section 7.8.2 Downmixing Into Two Channels */
    for (i = 0; i < s->nfchans; i++) {
        s->downmix_coeffs[i][0] = mixlevels[eac3_default_coeffs[s->channel_mode][i][0]];
        s->downmix_coeffs[i][1] = mixlevels[eac3_default_coeffs[s->channel_mode][i][1]];
    }

    if (get_bits1(gbc)) {
        /* Mixing metadata */
        if (s->channel_mode > 2) {
            /* if more than 2 channels */
            skip_bits(gbc, 2);  // skip preferred stereo downmix mode

            if (s->channel_mode & 1) {
                /* if three front channels exist */
                skip_bits(gbc, 3); //skip Lt/Rt center mix level
                s->downmix_coeffs[1][0] = s->downmix_coeffs[1][1] = mixlevels[get_bits(gbc, 3)];
            }
            if (s->channel_mode & 4) {
                /* if a surround channel exists */
                float surmixlev;
                skip_bits(gbc, 3); //skip Lt/Rt surround mix level
                surmixlev = mixlevels[get_bits(gbc, 3)];
                if (s->channel_mode & 2) {
                    //two surround channels
                    s->downmix_coeffs[s->channel_mode-4][0] = s->downmix_coeffs[s->channel_mode-3][1] =
                        surmixlev;
                } else {
                    s->downmix_coeffs[s->channel_mode-2][0] = s->downmix_coeffs[s->channel_mode-2][1] =
                        surmixlev * LEVEL_MINUS_3DB;
                }
            }
        }
        if (s->lfe_on) {
            /* if the LFE channel exists */
            if (get_bits1(gbc)) {
                // TODO: use LFE mix level
                skip_bits(gbc, 5); // skip LFE mix level code
            }
        }
        if (s->stream_type == EAC3_STREAM_TYPE_INDEPENDENT) {
            for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
                // TODO: apply program scale factor
                if (get_bits1(gbc)) {
                    skip_bits(gbc, 6);  // skip program scale factor
                }
            }
            if (get_bits1(gbc)) {
                skip_bits(gbc, 6);  // skip external program scale factor
            }
            /* skip mixing parameter data */
            switch(get_bits(gbc, 2)) {
                case 1: skip_bits(gbc, 5);                     break;
                case 2: skip_bits(gbc, 12);                    break;
                case 3: skip_bits(gbc, 8*get_bits(gbc, 5)+16); break;
            }
            /* skip pan information for mono or dual mono source */
            if (s->channel_mode < 2) {
                for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
                    if (get_bits1(gbc)) {
                        /* note: this is not in the ATSC A/52B specification
                           reference: ETSI TS 102 366 V1.1.1
                                      section: E.1.3.1.25 */
                        skip_bits(gbc, 8);  // skip Pan mean direction index
                        skip_bits(gbc, 6);  // skip reserved paninfo bits
                    }
                }
            }
            if (get_bits1(gbc)) {
                /* mixing configuration information */
                if (s->num_blocks == 1) {
                    s->blkmixcfginfo[0] = get_bits(gbc, 5);
                } else {
                    for (blk = 0; blk < s->num_blocks; blk++) {
                        if (get_bits1(gbc)) {
                            s->blkmixcfginfo[blk] = get_bits(gbc, 5);
                        }
                    }
                }
            }
        }
    }
    if (get_bits1(gbc)) {
        /* Informational metadata */
        skip_bits(gbc, 3); //skip Bit stream mode
        skip_bits(gbc, 2); //skip copyright bit and original bitstream bit
        if (s->channel_mode == AC3_CHMODE_STEREO) { /* if in 2/0 mode */
            skip_bits(gbc, 4); //skip Dolby surround and headphone mode
        }
        if (s->channel_mode >= 6) {
            /* if both surround channels exist */
            skip_bits(gbc, 2); //skip Dolby surround EX mode
        }
        for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
            if (get_bits1(gbc)) {
                skip_bits(gbc, 8); //skip Mix level, Room type and A/D converter type
            }
        }
        if (s->sr_code != EAC3_SR_CODE_REDUCED) {
            skip_bits1(gbc); //skip Source sample rate code
        }
    }
    if (s->stream_type == EAC3_STREAM_TYPE_INDEPENDENT && s->num_blocks != 6) {
        skip_bits1(gbc); //converter synchronization flag
    }
    if (s->stream_type == EAC3_STREAM_TYPE_AC3_CONVERT) {
        if (s->num_blocks == 6 || get_bits1(gbc)) {
            skip_bits(gbc, 6); // skip Frame size code
        }
    }
    if (get_bits1(gbc)) {
        int addbsil = get_bits(gbc, 6);
        for (i = 0; i < addbsil + 1; i++) {
            skip_bits(gbc, 8); // Additional bit stream information
        }
    }

    return 0;
} /* end of bsi */

static int parse_audfrm(GetBitContext *gbc, EAC3Context *s){
    int blk, ch;

    /* Audio frame exist flags and strategy data */
    if (s->num_blocks == 6) {
        /* LUT-based exponent strategy syntax */
        s->expstre = get_bits1(gbc);
        s->ahte = get_bits1(gbc);
    } else {
        /* AC-3 style exponent strategy syntax */
        s->expstre = 1;
        s->ahte = 0;
    }
    s->snroffststr = get_bits(gbc, 2);
    s->transproce = get_bits1(gbc);
    s->blkswe = get_bits1(gbc);
    if (!s->blkswe) {
        for (ch = 1; ch <= s->nfchans; ch++)
            s->blksw[ch] = 0;
    }
    s->dithflage = get_bits1(gbc);
    if (!s->dithflage) {
        for (ch = 1; ch <= s->nfchans; ch++)
            s->dithflag[ch] = 1; /* dither on */
    }
    s->dithflag[CPL_CH] = s->dithflag[s->lfe_channel] = 0;

    /* frame-based syntax flags */
    s->bamode = get_bits1(gbc);
    if (!s->bamode) {
        /* set default bit allocation parameters */
        s->bit_alloc_params.slow_decay = ff_ac3_slow_decay_tab[2];  /* Table 7.6 */
        s->bit_alloc_params.fast_decay = ff_ac3_fast_decay_tab[1];  /* Table 7.7 */
        s->bit_alloc_params.slow_gain  = ff_ac3_slow_gain_tab [1];  /* Table 7.8 */
        s->bit_alloc_params.db_per_bit = ff_ac3_db_per_bit_tab[2];  /* Table 7.9 */
        s->bit_alloc_params.floor      = ff_ac3_floor_tab     [7];  /* Table 7.10 */
    }
    s->frmfgaincode = get_bits1(gbc);
    s->dbaflde = get_bits1(gbc);
    s->skipflde = get_bits1(gbc);
    s->spxattene = get_bits1(gbc);
    /* Coupling data */
    if (s->channel_mode > 1) {
        s->cplstre[0] = 1;
        s->cplinu[0] = get_bits1(gbc);
        s->ncplblks = s->cplinu[0];
        for (blk = 1; blk < s->num_blocks; blk++) {
            s->cplstre[blk] = get_bits1(gbc);

            if (s->cplstre[blk]) {
                s->cplinu[blk] = get_bits1(gbc);
            } else {
                s->cplinu[blk] = s->cplinu[blk-1];
            }
            s->ncplblks += s->cplinu[blk];
        }
    } else {
        memset(s->cplinu, 0, sizeof(*s->cplinu) * s->num_blocks);
        s->ncplblks = 0;
    }

    /* Exponent strategy data */
    if (s->expstre) {
        /* AC-3 style exponent strategy syntax */
        for (blk = 0; blk < s->num_blocks; blk++) {
            for (ch = !s->cplinu[blk]; ch <= s->nfchans; ch++) {
                s->chexpstr[blk][ch] = get_bits(gbc, 2);
            }
        }
    } else {
        /* LUT-based exponent strategy syntax */
        int frmchexpstr;
        /* cplexpstr[blk] and chexpstr[blk][ch] derived from table lookups. see Table E2.14 */
        for (ch = !((s->channel_mode > 1) && (s->ncplblks)); ch <= s->nfchans; ch++) {
            frmchexpstr = get_bits(gbc, 5);
            for (blk = 0; blk < 6; blk++) {
                s->chexpstr[blk][ch] = ff_eac3_frm_expstr[frmchexpstr][blk];
            }
        }
    }
    /* LFE exponent strategy */
    if (s->lfe_on) {
        for (blk = 0; blk < s->num_blocks; blk++) {
            s->chexpstr[blk][s->lfe_channel] = get_bits1(gbc);
        }
    }
    /* Converter exponent strategy data */
    if (s->stream_type == EAC3_STREAM_TYPE_INDEPENDENT) {
        if (s->num_blocks == 6 || get_bits1(gbc)) {
            for (ch = 1; ch <= s->nfchans; ch++) {
                skip_bits(gbc, 5); //skip Converter channel exponent strategy
            }
        }
    }
    /* AHT data */
    if (s->ahte) {
        /* AHT is only available in 6 block mode (numblkscod ==3) */
        /* coupling can use AHT only when coupling in use for all blocks */
        /* ncplregs derived from cplstre and cplexpstr - see Section E3.3.2 */
        int nchregs;
        s->chahtinu[CPL_CH]=0;
        for (ch = (s->ncplblks != 6); ch <= s->num_channels; ch++) {
            nchregs = 0;
            for (blk = 0; blk < 6; blk++)
                nchregs += (s->chexpstr[blk][ch] != EXP_REUSE);
            s->chahtinu[ch] = (nchregs == 1) && get_bits1(gbc);
        }
    } else {
        memset(s->chahtinu, 0, sizeof(s->chahtinu));
    }
    /* Audio frame SNR offset data */
    if (!s->snroffststr) {
        int csnroffst = (get_bits(gbc, 6) - 15) << 4;
        int snroffst = (csnroffst + get_bits(gbc, 4)) << 2;
        for (ch = 0; ch <= s->num_channels; ch++)
            s->snroffst[ch] = snroffst;
    }
    /* Audio frame transient pre-noise processing data */
    if (s->transproce) {
        for (ch = 1; ch <= s->nfchans; ch++) {
            if (get_bits1(gbc)) { // channel in transient processing
                skip_bits(gbc, 10); // skip transient processing location
                skip_bits(gbc, 8);  // skip transient processing length
            }
        }
    }
    /* Spectral extension attenuation data */
    if (s->spxattene) {
        for (ch = 1; ch <= s->nfchans; ch++) {
            s->chinspxatten[ch] = get_bits1(gbc);
            if (s->chinspxatten[ch]) {
                s->spxattencod[ch] = get_bits(gbc, 5);
            }
        }
    } else {
        for (ch = 1; ch <= s->nfchans; ch++)
            s->chinspxatten[ch]=0;
    }
    /* Block start information */
    if (s->num_blocks > 1 && get_bits1(gbc)) {
        /* reference: Section E2.3.2.27
           nblkstrtbits = (numblks - 1) * (4 + ceiling(log2(words_per_frame)))
           Great that the spec tells how to parse. Unfortunately it doesn't say
           what this data is or what it's used for. */
        int block_start_bits = (s->num_blocks-1) * (4 + av_log2(s->frame_size-2));
        skip_bits(gbc, block_start_bits);
    }
    /* Syntax state initialization */
    for (ch = 1; ch <= s->nfchans; ch++) {
        s->firstspxcos[ch] = 1;
        s->firstcplcos[ch] = 1;
    }
    s->first_cpl_leak = 1;

    return 0;
} /* end of audfrm */

static int parse_audblk(GetBitContext *gbc, EAC3Context *s, const int blk){
    //int grp, sbnd, n, bin;
    int seg, bnd, ch, i, chbwcod, grpsize;
    int got_cplchan;
    mant_groups m;

    m.b1ptr = m.b2ptr = m.b4ptr = 3;

    /* Block switch and dither flags */
    if (s->blkswe) {
        for (ch = 1; ch <= s->nfchans; ch++) {
            s->blksw[ch] = get_bits1(gbc);
        }
    }
    if (s->dithflage) {
        for (ch = 1; ch <= s->nfchans; ch++) {
            s->dithflag[ch] = get_bits1(gbc);
        }
    }

    /* Dynamic range control */
    for (i = 0; i < (s->channel_mode ? 1 : 2); i++) {
        if (get_bits1(gbc)) {
            s->dynamic_range[i] = ff_ac3_dynamic_range_tab[get_bits(gbc, 8)];
        } else {
            if (!blk) {
                s->dynamic_range[i] = 1.0f;
            }
        }
    }
    /* Spectral extension strategy information */
    if ((!blk) || get_bits1(gbc)) {
        s->spxinu = get_bits1(gbc);
        if (s->spxinu) {
            log_missing_feature(s->avctx, "Spectral extension");
            return -1;
#if 0
            if (s->channel_mode == AC3_CHMODE_MONO) {
                s->chinspx[1] = 1;
            } else {
                for (ch = 1; ch <= s->nfchans; ch++) {
                    s->chinspx[ch] = get_bits1(gbc);
                }
            }
#if 0
            {
                int nspx=0;
                for (ch = 1; ch <= s->nfchans; ch++) {
                    nspx+=s->chinspx[ch];
                }
                if (!nspx)
                    av_log(s->avctx, AV_LOG_INFO, "No channels in spectral extension\n");
            }
#endif
            s->spxstrtf = get_bits(gbc, 2);
            s->spxbegf = get_bits(gbc, 3);
            s->spxendf = get_bits(gbc, 3);
            if (s->spxbegf < 6) {
                s->spxbegf += 2;
            } else {
                s->spxbegf = s->spxbegf * 2 - 3;
            }
            if (s->spxendf < 3) {
                s->spxendf += 5;
            } else {
                s->spxendf = s->spxendf * 2 + 3;
            }
            for (ch = 1; ch <= s->nfchans; ch++) {
                if (s->chinspx[ch])
                    s->endmant[ch] = 25 + 12 * s->spxbegf;
            }
            if (get_bits1(gbc)) {
                for (bnd = s->spxbegf + 1; bnd < s->spxendf; bnd++) {
                    s->spxbndstrc[bnd] = get_bits1(gbc);
                }
            } else {
                if (!blk) {
                    for (bnd = 0; bnd < 17; bnd++)
                        s->spxbndstrc[bnd] = ff_eac3_defspxbndstrc[bnd];
                }
            }
            // calculate number of spectral extension bands
            s->nspxbnds = 1;
            s->spxbndsztab[0] = 12;
            for (bnd = s->spxbegf+1; bnd < s->spxendf; bnd ++){
                if (!s->spxbndstrc[bnd]) {
                    s->spxbndsztab[s->nspxbnds] = 12;
                    s->nspxbnds++;
                } else {
                    s->spxbndsztab[s->nspxbnds - 1] += 12;
                }
            }
#endif
        } else {
            /* !spxinu */
            for (ch = 1; ch <= s->nfchans; ch++) {
                s->chinspx[ch] = 0;
                s->firstspxcos[ch] = 1;
            }
        }
    }

#if 0
    /* Spectral extension coordinates */
    if (s->spxinu) {
        for (ch = 1; ch <= s->nfchans; ch++) {
            if (s->chinspx[ch]) {
                if (s->firstspxcos[ch]) {
                    s->spxcoe[ch] = 1;
                    s->firstspxcos[ch] = 0;
                } else {
                    /* !firstspxcos[ch] */
                    s->spxcoe[ch] = get_bits1(gbc);
                }
                if (!blk && !s->spxcoe[ch]) {
                    av_log(s->avctx, AV_LOG_ERROR, "no spectral extension coordinates in first block\n");
                    return -1;
                }

                if (s->spxcoe[ch]) {
                    int spxcoexp, spxcomant, mstrspxco;
                    s->spxblnd[ch] = get_bits(gbc, 5);
                    mstrspxco = get_bits(gbc, 2);
                    mstrspxco*=3;
                    /* nspxbnds determined from spxbegf, spxendf, and spxbndstrc[ ] */
                    for (bnd = 0; bnd < s->nspxbnds; bnd++) {
                        spxcoexp = get_bits(gbc, 4);
                        spxcomant = get_bits(gbc, 2);
                        if (spxcoexp == 15)
                            s->spxco[ch][bnd] = spxcomant / 4.0f;
                        else
                            s->spxco[ch][bnd] = (spxcomant+4) / 8.0f;
                        s->spxco[ch][bnd] *= ff_ac3_scale_factors[spxcoexp + mstrspxco];
                    }
                }
            } else {
                /* !chinspx[ch] */
                s->firstspxcos[ch] = 1;
            }
        }
    }
#endif

    /* Coupling strategy and enhanced coupling strategy information */
    if (s->cplstre[blk]) {
        if (s->cplinu[blk]) {
            s->ecplinu = get_bits1(gbc);
            if (s->channel_mode == AC3_CHMODE_STEREO) {
                s->chincpl[1] = 1;
                s->chincpl[2] = 1;
            } else {
                for (ch = 1; ch <= s->nfchans; ch++) {
                    s->chincpl[ch] = get_bits1(gbc);
                }
            }
            if (!s->ecplinu) {
                /* standard coupling in use */
                if (s->channel_mode == AC3_CHMODE_STEREO) { /* if in 2/0 mode */
                    s->phsflginu = get_bits1(gbc);
                }
                s->cplbegf = get_bits(gbc, 4);
                if (!s->spxinu) {
                    /* if SPX not in use */
                    s->cplendf = get_bits(gbc, 4);
                    s->cplendf += 3;
                } else {
                    /* SPX in use */
                    s->cplendf = s->spxbegf - 1;
                }

                s->strtmant[CPL_CH] = 37 + (12 * s->cplbegf);
                s->endmant[CPL_CH] = 37 + (12 * s->cplendf);
                if (s->strtmant[CPL_CH] > s->endmant[CPL_CH]) {
                    av_log(s->avctx, AV_LOG_ERROR, "cplstrtmant > cplendmant [blk=%i]\n", blk);
                    return -1;
                }
                for (ch = 1; ch <= s->nfchans; ch++) {
                    if (s->chincpl[ch])
                        s->endmant[ch] = s->strtmant[CPL_CH];
                }
                if (get_bits1(gbc)) {
                    for (bnd = s->cplbegf + 1; bnd < s->cplendf; bnd++) {
                        s->cplbndstrc[bnd] = get_bits1(gbc);
                    }
                } else {
                    if (!blk) {
                        for (bnd = 0; bnd < 18; bnd++)
                            s->cplbndstrc[bnd] = ff_eac3_defcplbndstrc[bnd];
                    }
                }
                s->ncplsubnd =  s->cplendf - s->cplbegf;
                s->ncplbnd = s->ncplsubnd;
                for (bnd = s->cplbegf + 1; bnd < s->cplendf; bnd++) {
                    s->ncplbnd -= s->cplbndstrc[bnd];
                }
            } else {
                /* enhanced coupling in use */
                log_missing_feature(s->avctx, "Enhanced coupling");
                return -1;
#if 0
                s->ecplbegf = get_bits(gbc, 4);
                if (s->ecplbegf < 3) {
                    s->ecpl_start_subbnd = s->ecplbegf * 2;
                } else {
                    if (s->ecplbegf < 13) {
                        s->ecpl_start_subbnd = s->ecplbegf + 2;
                    } else {
                        s->ecpl_start_subbnd = s->ecplbegf * 2 - 10;
                    }
                }
                if (!s->spxinu) {
                    /* if SPX not in use */
                    s->ecplendf = get_bits(gbc, 4);
                    s->ecpl_end_subbnd = s->ecplendf + 7;
                } else {
                    /* SPX in use */
                    if (s->spxbegf < 6) {
                        s->ecpl_end_subbnd = s->spxbegf + 5;
                    } else {
                        s->ecpl_end_subbnd = s->spxbegf * 2;
                    }
                }
                if (get_bits1(gbc)) {
                    for (sbnd = FFMAX(9, s->ecpl_start_subbnd + 1);
                            sbnd < s->ecpl_end_subbnd; sbnd++){
                        s->ecplbndstrc[sbnd] = get_bits1(gbc);
                    }
                } else {
                    if (!blk) {
                        for (sbnd = 0; sbnd < 22; sbnd++)
                            s->ecplbndstrc[sbnd] = ff_eac3_defecplbndstrc[sbnd];
                    }
                }
                //necplbnd = ecpl_end_subbnd - ecpl_start_subbnd;
                //necplbnd -= ecplbndstrc[ecpl_start_subbnd] + ... + ecplbndstrc[ecpl_end_subbnd -1]
                s->necplbnd = s->ecpl_end_subbnd - s->ecpl_start_subbnd;
                for (bnd = s->ecpl_start_subbnd; bnd < s->ecpl_end_subbnd; bnd++) {
                    s->necplbnd -= s->ecplbndstrc[bnd];
                }
#endif
            } /* ecplinu[blk] */
        } else {
            /* !cplinu[blk] */
            for (ch = 1; ch <= s->nfchans; ch++) {
                s->chincpl[ch] = 0;
                s->firstcplcos[ch] = 1;
            }
            s->first_cpl_leak = 1;
            s->phsflginu = 0;
            s->ecplinu = 0;
        }
    } /* cplstre[blk] */
    /* Coupling coordinates */
    if (s->cplinu[blk]) {
        if (!s->ecplinu) {
            /* standard coupling in use */
            for (ch = 1; ch <= s->nfchans; ch++) {
                if (s->chincpl[ch]) {
                    if (s->firstcplcos[ch]) {
                        s->cplcoe[ch] = 1;
                        s->firstcplcos[ch] = 0;
                    } else {
                        /* !firstcplcos[ch] */
                        s->cplcoe[ch] = get_bits1(gbc);
                    }
                    if (s->cplcoe[ch]) {
                        int cplcoexp, cplcomant, mstrcplco;
                        mstrcplco = get_bits(gbc, 2);
                        mstrcplco = 3 * mstrcplco;
                        /* ncplbnd derived from cplbegf, cplendf, and cplbndstrc */
                        for (bnd = 0; bnd < s->ncplbnd; bnd++) {
                            cplcoexp = get_bits(gbc, 4);
                            cplcomant = get_bits(gbc, 4);
                            if (cplcoexp == 15)
                                s->cplco[ch][bnd] = cplcomant / 16.0f;
                            else
                                s->cplco[ch][bnd] = (cplcomant + 16.0f) / 32.0f;
                            s->cplco[ch][bnd] *=  ff_ac3_scale_factors[cplcoexp + mstrcplco];
                        }
                    } /* cplcoe[ch] */
                    else {
                        if (!blk) {
                            av_log(s->avctx, AV_LOG_ERROR,  "no coupling coordinates in first block\n");
                            return -1;
                        }
                    }
                } else {
                    /* ! chincpl[ch] */
                    s->firstcplcos[ch] = 1;
                }
            } /* ch */
            if ((s->channel_mode == AC3_CHMODE_STEREO) && s->phsflginu
                    && (s->cplcoe[1] || s->cplcoe[2])) {
                for (bnd = 0; bnd < s->ncplbnd; bnd++) {
                    s->phsflg[bnd] = get_bits1(gbc);
                }
            }
            s->nchgrps[CPL_CH] = (s->endmant[CPL_CH] - s->strtmant[CPL_CH]) /
                (3 << (s->chexpstr[blk][CPL_CH] - 1));
        } else {
            /* enhanced coupling in use */
            //TODO calc nchgrps[CPL_CH]
#if 0
            s->firstchincpl = -1;
            s->ecplangleintrp = get_bits1(gbc);
            for (ch = 1; ch <= s->nfchans; ch++) {
                if (s->chincpl[ch]) {
                    if (s->firstchincpl == -1) {
                        s->firstchincpl = ch;
                    }
                    if (s->firstcplcos[ch]) {
                        s->ecplparam1e[ch] = 1;
                        if (ch > s->firstchincpl) {
                            s->ecplparam2e[ch] = 1;
                        } else {
                            s->ecplparam2e[ch] = 0;
                        }
                        s->firstcplcos[ch] = 0;
                    } else {
                        /* !firstcplcos[ch] */
                        s->ecplparam1e[ch] = get_bits1(gbc);
                        if (ch > s->firstchincpl) {
                            s->ecplparam2e[ch] = get_bits1(gbc);
                        } else {
                            s->ecplparam2e[ch] = 0;
                        }
                    }
                    if (s->ecplparam1e[ch]) {
                        /* necplbnd derived from ecpl_start_subbnd, ecpl_end_subbnd, and ecplbndstrc */
                        for (bnd = 0; bnd < s->necplbnd; bnd++) {
                            s->ecplamp[ch][bnd] = get_bits(gbc, 5);
                        }
                    }
                    if (s->ecplparam2e[ch]) {
                        /* necplbnd derived from ecpl_start_subbnd, ecpl_end_subbnd, and ecplbndstrc */
                        for (bnd = 0; bnd < s->necplbnd; bnd++) {
                            s->ecplangle[ch][bnd] = get_bits(gbc, 6);
                            s->ecplchaos[ch][bnd] = get_bits(gbc, 3);
                        }
                    }
                    if (ch > s->firstchincpl) {
                        s->ecpltrans[ch] = get_bits1(gbc);
                    }
                } else {
                    /* !chincpl[ch] */
                    s->firstcplcos[ch] = 1;
                }
            } /* ch */
#endif
        } /* ecplinu[blk] */
    } /* cplinu[blk] */
    /* Rematrixing operation in the 2/0 mode */
    if (s->channel_mode == AC3_CHMODE_STEREO) { /* if in 2/0 mode */
        if (!blk || get_bits1(gbc)) {
            /* nrematbnds determined from cplinu, ecplinu, spxinu, cplbegf, ecplbegf and spxbegf */
            // TODO spx in one channel
            int end = (s->cplinu[blk] || s->spxinu) ?
                FFMIN(s->endmant[1], s->endmant[2]) : (ff_ac3_rematrix_band_tab[4]-1);
            for (bnd = 0; ff_ac3_rematrix_band_tab[bnd] <= end; bnd++) {
                s->rematflg[bnd] = get_bits1(gbc);
            }
            s->nrematbnds = bnd;
        }
    }
    /* Channel bandwidth code */
    for (ch = 1; ch <= s->nfchans; ch++) {
        if (!blk && s->chexpstr[blk][ch] == EXP_REUSE) {
            av_log(s->avctx, AV_LOG_ERROR,  "no channel exponent strategy in first block\n");
            return -1;
        }
        if (s->chexpstr[blk][ch] != EXP_REUSE) {
            grpsize = 3 << (s->chexpstr[blk][ch] - 1);
            s->strtmant[ch] = 0;
            if ((!s->chincpl[ch]) && (!s->chinspx[ch])) {
                chbwcod = get_bits(gbc, 6);
                if (chbwcod > 60) {
                    av_log(s->avctx, AV_LOG_ERROR, "chbwcod > 60\n");
                    return -1;
                }
                s->endmant[ch] = ((chbwcod + 12) * 3) + 37; /* (ch is not coupled) */
            }
            grpsize = 3 << (s->chexpstr[blk][ch] - 1);
            s->nchgrps[ch] = (s->endmant[ch] + grpsize - 4) / grpsize;
        }
    }
    /* Exponents */
    for (ch = !s->cplinu[blk]; ch <= s->num_channels; ch++) {
        if (s->chexpstr[blk][ch] != EXP_REUSE) {
            s->dexps[ch][0] = get_bits(gbc, 4) << !ch;
            ff_ac3_decode_exponents(gbc, s->chexpstr[blk][ch], s->nchgrps[ch],
                    s->dexps[ch][0], s->dexps[ch]+s->strtmant[ch]+!!ch);
            if (ch != CPL_CH && ch != s->lfe_channel)
                skip_bits(gbc, 2); /* skip gainrng */
        }
    }

    /* Bit-allocation parametric information */
    if (s->bamode) {
        if (get_bits1(gbc)) {
            s->bit_alloc_params.slow_decay = ff_ac3_slow_decay_tab[get_bits(gbc, 2)];   /* Table 7.6 */
            s->bit_alloc_params.fast_decay = ff_ac3_fast_decay_tab[get_bits(gbc, 2)];   /* Table 7.7 */
            s->bit_alloc_params.slow_gain  = ff_ac3_slow_gain_tab [get_bits(gbc, 2)];   /* Table 7.8 */
            s->bit_alloc_params.db_per_bit = ff_ac3_db_per_bit_tab[get_bits(gbc, 2)];   /* Table 7.9 */
            s->bit_alloc_params.floor      = ff_ac3_floor_tab     [get_bits(gbc, 3)];   /* Table 7.10 */
        } else {
            if (!blk) {
                av_log(s->avctx, AV_LOG_ERROR, "no bit allocation information in first block\n");
                return -1;
            }
        }
    }

    if (s->snroffststr) {
        if (!blk || get_bits1(gbc)) {
            int csnroffst = (get_bits(gbc, 6) - 15) << 4;
            int snroffst = 0;
            for (i = !s->cplinu[blk]; ch <= s->num_channels; ch++){
                if (ch == !s->cplinu[blk] || s->snroffststr == 2)
                    snroffst = (csnroffst + get_bits(gbc, 4)) << 2;
                s->snroffst[ch] = snroffst;
            }
        }
    }

    if (s->frmfgaincode && get_bits1(gbc)) {
        for (ch = !s->cplinu[blk]; ch <= s->num_channels; ch++)
            s->fgain[ch] = ff_ac3_fast_gain_tab[get_bits(gbc, 3)];
    } else {
        if (!blk) {
            for (ch = !s->cplinu[blk]; ch <= s->num_channels; ch++)
                s->fgain[ch] = ff_ac3_fast_gain_tab[4];
        }
    }
    if (s->stream_type == EAC3_STREAM_TYPE_INDEPENDENT) {
        if (get_bits1(gbc)) {
            skip_bits(gbc, 10); //Converter SNR offset
        }
    }
    if (s->cplinu[blk]) {
        if (s->first_cpl_leak || get_bits1(gbc)) {
            s->bit_alloc_params.cpl_fast_leak = get_bits(gbc, 3);
            s->bit_alloc_params.cpl_slow_leak = get_bits(gbc, 3);
        }
        if(s->first_cpl_leak)
            s->first_cpl_leak = 0;
    }
    /* Delta bit allocation information */
    if (s->dbaflde && get_bits1(gbc)) {
        for (ch = !s->cplinu[blk]; ch <= s->nfchans; ch++) {
            s->deltbae[ch] = get_bits(gbc, 2);
        }
        for (ch = !s->cplinu[blk]; ch <= s->nfchans; ch++) {
            if (s->deltbae[ch] == DBA_NEW) {
                s->deltnseg[ch] = get_bits(gbc, 3);
                for (seg = 0; seg <= s->deltnseg[ch]; seg++) {
                    s->deltoffst[ch][seg] = get_bits(gbc, 5);
                    s->deltlen[ch][seg] = get_bits(gbc, 4);
                    s->deltba[ch][seg] = get_bits(gbc, 3);
                }
            }
        }
    } else {
        if (!blk) {
            for (ch = 0; ch <= s->num_channels; ch++) {
                s->deltbae[ch] = DBA_NONE;
            }
        }
    }

    /* Inclusion of unused dummy data */
    if (s->skipflde) {
        if (get_bits1(gbc)) {
            int skipl = get_bits(gbc, 9);
            while(skipl--) skip_bits(gbc, 8);
        }
    }

    /* run bit allocation */
    for (ch = !s->cplinu[blk]; ch <= s->num_channels; ch++) {
        ff_ac3_bit_alloc_calc_psd((int8_t *)s->dexps[ch], s->strtmant[ch],
                s->endmant[ch], s->psd[ch], s->bndpsd[ch]);

        s->bit_alloc_params.sr_code = s->sr_code;
        s->bit_alloc_params.sr_shift = 0;

        ff_ac3_bit_alloc_calc_mask(&s->bit_alloc_params,
                s->bndpsd[ch], s->strtmant[ch], s->endmant[ch], s->fgain[ch],
                (ch == s->lfe_channel),
                s->deltbae[ch], s->deltnseg[ch],
                s->deltoffst[ch], s->deltlen[ch],
                s->deltba[ch], s->mask[ch]);

        if (s->chahtinu[ch] == 0)
            ff_ac3_bit_alloc_calc_bap(s->mask[ch], s->psd[ch], s->strtmant[ch],
                    s->endmant[ch], s->snroffst[ch], s->bit_alloc_params.floor, ff_ac3_bap_tab,
                    s->bap[ch]);
        else
            if (s->chahtinu[ch] == 1)
                ff_ac3_bit_alloc_calc_bap(s->mask[ch], s->psd[ch], s->strtmant[ch], s->endmant[ch],
                        s->snroffst[ch], s->bit_alloc_params.floor, ff_ac3_hebaptab,
                        s->hebap[ch]);
    }

    got_cplchan = 0;

    /* Quantized mantissa values */
    for (ch = 1; ch <= s->num_channels; ch++) {
        get_eac3_transform_coeffs_ch(gbc, s, blk, ch, &m);
        if (s->cplinu[blk] && s->chincpl[ch] && !got_cplchan) {
            get_eac3_transform_coeffs_ch(gbc, s, blk, CPL_CH, &m);
            got_cplchan = 1;
        }
    }

    if (s->cplinu[blk])
        uncouple_channels(s);

    //apply spectral extension
    if (s->spxinu)
        spectral_extension(s);

    return 0;
}

/**
 * Performs Inverse MDCT transform
 */
static void do_imdct(EAC3Context *ctx){
    int ch;

    for (ch = 1; ch <= ctx->nfchans + ctx->lfe_on; ch++) {
        if (ctx->blksw[ch]) {
            /* 256-point IMDCT */
            ff_ac3_do_imdct_256(ctx->tmp_output, ctx->transform_coeffs[ch],
                    &ctx->imdct_256, ctx->tmp_imdct);
        } else {
            /* 512-point IMDCT */
            ctx->imdct_512.fft.imdct_calc(&ctx->imdct_512, ctx->tmp_output,
                    ctx->transform_coeffs[ch],
                    ctx->tmp_imdct);
        }
        /* apply window function, overlap/add output, save delay */
        ctx->dsp.vector_fmul_add_add(ctx->output[ch-1], ctx->tmp_output,
                ctx->window, ctx->delay[ch-1], 0,
                AC3_BLOCK_SIZE, 1);
        ctx->dsp.vector_fmul_reverse(ctx->delay[ch-1], ctx->tmp_output+256,
                ctx->window, AC3_BLOCK_SIZE);
    }
}

static int eac3_decode_frame(AVCodecContext *avctx, void *data, int *data_size,
        uint8_t *buf, int buf_size){
    int16_t *out_samples = (int16_t *)data;
    EAC3Context *c = (EAC3Context *)avctx->priv_data;
    int k, i, blk, ch;
    GetBitContext gbc;

    *data_size = 0;
    c->gbc = &gbc;
    c->syncword = 0;
    init_get_bits(&gbc, buf, buf_size*8);
    c->syncword = get_bits(&gbc, 16);

    if (c->syncword != 0x0B77)
        return -1;

    if (parse_bsi(&gbc, c) || parse_audfrm(&gbc, c))
        return -1;

    if (c->substreamid) {
        // TODO: allow user to select which substream to decode
        av_log(avctx, AV_LOG_INFO, "Skipping additional substream #%d\n",
               c->substreamid);
        return -1;
    }

    if (c->sr_code == EAC3_SR_CODE_REDUCED) {
        avctx->sample_rate = ff_ac3_sample_rate_tab[c->sr_code2] / 2;
    } else {
        avctx->sample_rate = ff_ac3_sample_rate_tab[c->sr_code];
    }

    avctx->bit_rate = (c->frame_size * avctx->sample_rate * 8 / (c->num_blocks * 256)) / 1000;

    /* channel config */
    if (!avctx->request_channels) {
        if (!avctx->channels)
            avctx->channels = c->num_channels;
    } else {
        if (c->num_channels < avctx->request_channels) {
            av_log(avctx, AV_LOG_ERROR, "Cannot upmix EAC3 from %d to %d channels.\n",
                    c->num_channels, avctx->request_channels);
            return -1;
        } else {
            if (avctx->request_channels > 2
                    && avctx->request_channels != c->num_channels) {
                av_log(avctx, AV_LOG_ERROR, "Cannot downmix EAC3 from %d to %d channels.\n",
                        c->num_channels, avctx->request_channels);
                return -1;
            }
            avctx->channels = avctx->request_channels;
        }
    }

    for (blk = 0; blk < c->num_blocks; blk++) {
        if (parse_audblk(&gbc, c, blk)) {
            av_log(c->avctx, AV_LOG_ERROR, "Error in parse_audblk\n");
            return -1;
        }

        /* recover coefficients if rematrixing is in use */
        if (c->channel_mode == AC3_CHMODE_STEREO)
            ff_ac3_do_rematrixing(c->transform_coeffs,
                    FFMIN(c->endmant[1], c->endmant[2]),
                    c->nrematbnds, c->rematflg);

        /* apply scaling to coefficients (dialnorm, dynrng) */
        for (ch = 1; ch <= c->nfchans + c->lfe_on; ch++) {
            float gain=2.0f;
            if (c->channel_mode == AC3_CHMODE_DUALMONO) {
                gain *= c->dialog_norm[ch-1] * c->dynamic_range[ch-1];
            } else {
                gain *= c->dialog_norm[0] * c->dynamic_range[0];
            }
            for (i = 0; i < c->endmant[ch]; i++) {
                c->transform_coeffs[ch][i] *= gain;
            }
        }

        do_imdct(c);

        // TODO: Transient Pre-Noise Cross-Fading

        if (avctx->channels != c->num_channels) {
            ff_ac3_downmix(c->output, c->nfchans, avctx->channels, c->downmix_coeffs);
        }

        // convert float to 16-bit integer
        for (ch = 0; ch < avctx->channels; ch++) {
            for (i = 0; i < AC3_BLOCK_SIZE; i++) {
                c->output[ch][i] = c->output[ch][i] * c->mul_bias +
                    c->add_bias;
            }
            c->dsp.float_to_int16(c->int_output[ch], c->output[ch],
                    AC3_BLOCK_SIZE);
        }
        for (k = 0; k < AC3_BLOCK_SIZE; k++) {
            for (i = 0; i < avctx->channels; i++) {
                *(out_samples++) = c->int_output[i][k];
            }
        }
    }

    *data_size = c->num_blocks * 256 * avctx->channels * sizeof (int16_t); // TODO is ok?

    return c->frame_size;
}

static void eac3_tables_init(void) {
    int blk, i;

    // initialize IDCT cosine table for use with AHT
    for(blk=0; blk<6; blk++) {
        for(i=1; i<6; i++) {
            idct_cos_tab[blk][i-1] = M_SQRT2 * cos(M_PI*i*(2*blk + 1)/12);
        }
    }
}

static int eac3_decode_init(AVCodecContext *avctx){
    EAC3Context *ctx = avctx->priv_data;

    ctx->avctx = avctx;
    ac3_common_init();
    ff_ac3_tables_init();
    eac3_tables_init();
    av_init_random(0, &ctx->dith_state);
    ff_mdct_init(&ctx->imdct_256, 8, 1);
    ff_mdct_init(&ctx->imdct_512, 9, 1);
    dsputil_init(&ctx->dsp, avctx);
    if (ctx->dsp.float_to_int16 == ff_float_to_int16_c) {
        ctx->add_bias = 385.0f;
        ctx->mul_bias = 1.0f;
    } else {
        ctx->add_bias = 0.0f;
        ctx->mul_bias = 32767.0f;
    }
    ff_ac3_window_init(ctx->window);
    return 0;
}

static int eac3_decode_end(AVCodecContext *avctx){
    EAC3Context *ctx = avctx->priv_data;
    ff_mdct_end(&ctx->imdct_512);
    ff_mdct_end(&ctx->imdct_256);

    return 0;
}

AVCodec eac3_decoder = {
    .name = "E-AC3",
    .type = CODEC_TYPE_AUDIO,
    .id = CODEC_ID_EAC3,
    .priv_data_size = sizeof (EAC3Context),
    .init = eac3_decode_init,
    .close = eac3_decode_end,
    .decode = eac3_decode_frame,

};
