/*
 * E-AC3 parser
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

//#define DEBUG

#include "avcodec.h"
#include "bitstream.h"
#include "ac3.h"
#include "random.h"

#ifdef DEBUG
#define GET_BITS(a, gbc, n) a = get_bits(gbc, n); av_log(NULL, AV_LOG_INFO, "%s: %i\n", __STRING(a), a)
#else
#define GET_BITS(a, gbc, n) a = get_bits(gbc, n)
#endif //DEBUG

#include "eac3.h"

int eac3_parse_syncinfo(GetBitContext *gbc, EAC3Context *s){
    GET_BITS(s->syncword, gbc, 16);
    return 0;
}

int eac3_parse_bsi(GetBitContext *gbc, EAC3Context *s){
    int i, blk;

    GET_BITS(s->strmtyp, gbc, 2);
    GET_BITS(s->substreamid, gbc, 3);
    GET_BITS(s->frmsiz, gbc, 11);
    GET_BITS(s->fscod, gbc, 2);
    if(s->fscod == 0x3)
    {
        av_log(s, AV_LOG_ERROR, "NOT IMPLEMENTED");
        return -1;
#if 0
        GET_BITS(s->fscod2, gbc, 2);
        s->numblkscod = 0x3; /* six blocks per frame */
#endif
    }
    else
    {
        GET_BITS(s->numblkscod, gbc, 2);
    }
    GET_BITS(s->acmod, gbc, 3);
    GET_BITS(s->lfeon, gbc, 1);

    // calculate number of channels
    s->nfchans = ff_ac3_channels[s->acmod];
    s->ntchans = s->nfchans;
    if(s->lfeon){
        s->lfe_channel = s->ntchans+1;

        s->nlfemant = 7;
        s->endmant[s->lfe_channel] = 7;
        s->nchgrps[s->lfe_channel] = 2;

        s->ntchans ++ ;
    }

    GET_BITS(s->bsid, gbc, 5);
    if(s->bsid < 11 || s->bsid > 16){
        av_log(s, AV_LOG_ERROR, "bsid should be between 11 and 16");
        return -1;
    }

    GET_BITS(s->dialnorm, gbc, 5);
    GET_BITS(s->compre, gbc, 1);
    if(s->compre) {
        GET_BITS(s->compr, gbc, 8);
    }
    if(s->acmod == 0x0) /* if 1+1 mode (dual mono, so some items need a second value) */
    {
        GET_BITS(s->dialnorm2, gbc, 5);
        GET_BITS(s->compr2e, gbc, 1);
        if(s->compr2e) {
            GET_BITS(s->compr2, gbc, 8);
        }

    }
    if(s->strmtyp == 0x1) /* if dependent stream */
    {
        GET_BITS(s->chanmape, gbc, 1);
        if(s->chanmape) {
            GET_BITS(s->chanmap, gbc, 16);
        }
    }
    GET_BITS(s->mixmdate, gbc, 1);
    if(s->mixmdate) /* Mixing metadata */
    {
        if(s->acmod > 0x2) /* if more than 2 channels */ {
            GET_BITS(s->dmixmod, gbc, 2);
        }
        if((s->acmod & 0x1) && (s->acmod > 0x2)) /* if three front channels exist */
        {
            GET_BITS(s->ltrtcmixlev, gbc, 3);
            GET_BITS(s->lorocmixlev, gbc, 3);
        }
        if(s->acmod & 0x4) /* if a surround channel exists */
        {

            GET_BITS(s->ltrtsurmixlev, gbc, 3);
            GET_BITS(s->lorosurmixlev, gbc, 3);
        }
        if(s->lfeon) /* if the LFE channel exists */
        {
            GET_BITS(s->lfemixlevcode, gbc, 1);
            if(s->lfemixlevcode) {
                GET_BITS(s->lfemixlevcod, gbc, 5);
            }
        }
        if(s->strmtyp == 0x0) /* if independent stream */
        {
            GET_BITS(s->pgmscle, gbc, 1);
            if(s->pgmscle) {
                GET_BITS(s->pgmscl, gbc, 6);
            }
            if(s->acmod == 0x0) /* if 1+1 mode (dual mono, so some items need a second value) */
            {
                GET_BITS(s->pgmscl2e, gbc, 1);
                if(s->pgmscl2e) {
                    GET_BITS(s->pgmscl2, gbc, 6);
                }
            }
            GET_BITS(s->extpgmscle, gbc, 1);
            if(s->extpgmscle) {
                GET_BITS(s->extpgmscl, gbc, 6);
            }
            GET_BITS(s->mixdef, gbc, 2);
            if(s->mixdef == 0x1) /* mixing option 2 */ {
                GET_BITS(s->mixdata, gbc, 5);
            }
            else if(s->mixdef == 0x2) /* mixing option 3 */ {
                GET_BITS(s->mixdata, gbc, 12);
            }
            else if(s->mixdef == 0x3) /* mixing option 4 */
            {
                GET_BITS(s->mixdeflen, gbc, 5);
                av_log(s, AV_LOG_ERROR, "NOT IMPLEMENTED");
                return -1;
//                GET_BITS(s->mixdata, gbc, 8*(mixdeflen+2));
            }
            if(s->acmod < 0x2) /* if mono or dual mono source */
            {
                GET_BITS(s->paninfoe, gbc, 1);
                if(s->paninfoe) {
                    GET_BITS(s->paninfo, gbc, 14);
                }
                if(s->acmod == 0x0) /* if 1+1 mode (dual mono, so some items need a second value) */
                {
                    GET_BITS(s->paninfo2e, gbc, 1);
                    if(s->paninfo2e) {
                        GET_BITS(s->paninfo2, gbc, 14);
                    }
                }
            }
            GET_BITS(s->frmmixcfginfoe, gbc, 1);
            if(s->frmmixcfginfoe) /* mixing configuration information */
            {
                if(s->numblkscod == 0x0) {
                    GET_BITS(s->blkmixcfginfo0, gbc, 5);
                }
                else
                {
                    for(blk = 0; blk < ff_eac3_blocks[s->numblkscod]; blk++)
                    {
                        GET_BITS(s->blkmixcfginfoe, gbc, 1);
                        if(s->blkmixcfginfoe){
                            GET_BITS(s->blkmixcfginfoblk, gbc, 5);
                        }
                    }
                }
            }
        }
    }
    GET_BITS(s->infomdate, gbc, 1);
    if(s->infomdate) /* Informational metadata */
    {
        GET_BITS(s->bsmod, gbc, 3);

        GET_BITS(s->copyrightb, gbc, 1);
        GET_BITS(s->origbs, gbc, 1);
        if(s->acmod == 0x2) /* if in 2/0 mode */
        {
            GET_BITS(s->dsurmod, gbc, 2);
            GET_BITS(s->dheadphonmod, gbc, 2);
        }
        if(s->acmod >= 0x6) /* if both surround channels exist */ {
            GET_BITS(s->dsurexmod, gbc, 2);
        }
        GET_BITS(s->audprodie, gbc, 1);
        if(s->audprodie)
        {
            GET_BITS(s->mixlevel, gbc, 5);
            GET_BITS(s->roomtyp, gbc, 2);
            GET_BITS(s->adconvtyp, gbc, 1);
        }
        if(s->acmod == 0x0) /* if 1+1 mode (dual mono, so some items need a second value) */
        {
            GET_BITS(s->audprodi2e, gbc, 1);
            if(s->audprodi2e)
            {
                GET_BITS(s->mixlevel2, gbc, 5);
                GET_BITS(s->roomtyp2, gbc, 2);
                GET_BITS(s->adconvtyp2, gbc, 1);
            }
        }
        if(s->fscod < 0x3) /* if not half sample rate */ {
            GET_BITS(s->sourcefscod, gbc, 1); // TODO
        }
    }
    if((s->strmtyp == 0x0) && (s->numblkscod != 0x3) ) {
        GET_BITS(s->convsync, gbc, 1);
    }
    if(s->strmtyp == 0x2) /* if bit stream converted from AC-3 */
    {
        if(s->numblkscod == 0x3) /* 6 blocks per frame */ {
            s->blkid = 1;
        }
        else {
            GET_BITS(s->blkid, gbc, 1);
        }
        if(s->blkid) {
            GET_BITS(s->frmsizecod, gbc, 6);
        }
    }
    GET_BITS(s->addbsie, gbc, 1);
    if(s->addbsie)
    {
        GET_BITS(s->addbsil, gbc, 6);
        for(i=0; i<s->addbsil+1; i++){
            GET_BITS(s->addbsi[i], gbc, 8);
        }
    }

    return 0;
} /* end of bsi */


int eac3_parse_audfrm(GetBitContext *gbc, EAC3Context *s){
    int blk, ch;

    /* These fields for audio frame exist flags and strategy data */
    if(s->numblkscod == 0x3) /* six blocks per frame */
    {
        GET_BITS(s->expstre, gbc, 1);
        GET_BITS(s->ahte, gbc, 1);
    }
    else
    {
        s->expstre = 1;
        s->ahte = 0;
    }
    GET_BITS(s->snroffststr, gbc, 2);
    GET_BITS(s->transproce, gbc, 1);
    GET_BITS(s->blkswe, gbc, 1);
    GET_BITS(s->dithflage, gbc, 1);
    GET_BITS(s->bamode, gbc, 1);
    GET_BITS(s->frmfgaincode, gbc, 1);
    GET_BITS(s->dbaflde, gbc, 1);
    GET_BITS(s->skipflde, gbc, 1);
    GET_BITS(s->spxattene, gbc, 1);
    /* These fields for coupling data */
    if(s->acmod > 0x1)
    {
        s->cplstre[0] = 1;
        GET_BITS(s->cplinu[0], gbc, 1);
        for(blk = 1; blk < ff_eac3_blocks[s->numblkscod]; blk++)
        {
            GET_BITS(s->cplstre[blk], gbc, 1);

            if(s->cplstre[blk] == 1) {
                GET_BITS(s->cplinu[blk], gbc, 1);
            }
            else {
                s->cplinu[blk] = s->cplinu[blk-1];
            }
        }
    }
    else
    {
        for(blk = 0; blk < ff_eac3_blocks[s->numblkscod]; blk++) {
            s->cplinu[blk] = 0;
        }
    }

    // calculate number of coupling blocks
    s->ncplblks = 0;
    for(blk = 0; blk < ff_eac3_blocks[s->numblkscod]; blk++) {
        s->ncplblks += s->cplinu[blk];
    }

    /* These fields for exponent strategy data */
    if(s->expstre)
    {
        for(blk = 0; blk < ff_eac3_blocks[s->numblkscod]; blk++)
        {
            if(s->cplinu[blk] == 1) {
                GET_BITS(s->cplexpstr[blk], gbc, 2);
            }

            for(ch = 1; ch <= s->nfchans; ch++) {
                GET_BITS(s->chexpstr[blk][ch], gbc, 2); // TODO
            }
        }
    }
    else
    {
        if( (s->acmod > 0x1) && (s->ncplblks > 0) ) {
            GET_BITS(s->frmcplexpstr, gbc, 5);
        }
        for(ch = 1; ch <= s->nfchans; ch++) {
            GET_BITS(s->frmchexpstr[ch], gbc, 5);
        }
        /* cplexpstr[blk] and chexpstr[blk][ch] derived from table lookups ? see Table E2.14 */
        av_log(s, AV_LOG_ERROR, "NOT IMPLEMENTED");
        return -1;
    }
    if(s->lfeon)
    {
        for(blk = 0; blk < ff_eac3_blocks[s->numblkscod]; blk++) {
            GET_BITS(s->lfeexpstr[blk], gbc, 1);
            s->chexpstr[blk][s->lfe_channel] = s->lfeexpstr[blk];
        }
    }
    /* These fields for converter exponent strategy data */
    if(s->strmtyp == 0x0)
    {
        if(s->numblkscod != 0x3) {
            GET_BITS(s->convexpstre, gbc, 1);
        }
        else {
            s->convexpstre = 1;
        }
        if(s->convexpstre == 1)
        {
            for(ch = 1; ch <= s->nfchans; ch++) {
                GET_BITS(s->convexpstr[ch], gbc, 5);
            }
        }
    }
    /* These fields for AHT data */
    if(s->ahte)
    {
        /* coupling can use AHT only when coupling in use for all blocks */
        /* ncplregs derived from cplstre and cplexpstr ? see Section E3.3.2 */
        av_log(s, AV_LOG_ERROR, "AHT NOT IMPLEMENTED");
        return -1;
#if 0
        /*
        if( (s->ncplblks == 6) && (s->ncplregs ==1) ) {
            GET_BITS(s->cplahtinu, gbc, 1);
        }
        else {
            cplahtinu = 0
        }

        for(ch = 1; ch <= s->nfchans; ch++)
        {
            // nchregs derived from chexpstr ? see Section E3.3.2
            if(s->nchregs[ch] == 1) {
                GET_BITS(s->chahtinu[ch], gbc, 1);
            }
            else {
                chahtinu[ch] = 0
            }
        }
        if(lfeon)
        {
            // nlferegs derived from lfeexpstr ? see Section E3.3.2
            if(nlferegs == 1) {
                GET_BITS(s->lfeahtinu, gbc, 1);
            }
            else {
                lfeahtinu = 0
            }
        }
        */
#endif
    }
    /* These fields for audio frame SNR offset data */
    if(s->snroffststr == 0x0)
    {
        GET_BITS(s->frmcsnroffst, gbc, 6);
        GET_BITS(s->frmfsnroffst, gbc, 4);
    }
    /* These fields for audio frame transient pre-noise processing data */
    if(s->transproce)
    {
        for(ch = 1; ch <= s->nfchans; ch++)
        {
            GET_BITS(s->chintransproc[ch], gbc, 1);
            if(s->chintransproc[ch])
            {
                GET_BITS(s->transprocloc[ch], gbc, 10);
                GET_BITS(s->transproclen[ch], gbc, 8);
            }
        }
    }
    /* These fields for spectral extension attenuation data */
    if(s->spxattene)
    {
        for(ch = 1; ch <= s->nfchans; ch++)
        {
            GET_BITS(s->chinspxatten[ch], gbc, 1);
            if(s->chinspxatten[ch])
            {
                GET_BITS(s->spxattencod[ch], gbc, 5);
            }
        }
    }
    /* These fields for block start information */
    if (s->numblkscod != 0x0) {
        GET_BITS(s->blkstrtinfoe, gbc, 1);
    }
    else {
        s->blkstrtinfoe = 0;
    }
    if(s->blkstrtinfoe)
    {
        /* nblkstrtbits determined from frmsiz (see Section E2.3.2.27) */
        // nblkstrtbits = (numblks - 1) * (4 + ceiling (log2 (words_per_frame)))
        // where numblks is derived from the numblkscod in Table E2.15 <- TODO ???
        // words_per_frame = frmsiz + 1
        int nblkstrtbits = (ff_eac3_blocks[s->numblkscod]-1) * (4 + (av_log2(s->frmsiz-1)+1) );
        av_log(NULL, AV_LOG_INFO, "nblkstrtbits = %i\n", nblkstrtbits);
        GET_BITS(s->blkstrtinfo, gbc, nblkstrtbits);
    }
    /* These fields for syntax state initialization */
    for(ch = 1; ch <= s->nfchans; ch++)
    {
        s->firstspxcos[ch] = 1;
        s->firstcplcos[ch] = 1;
    }
    s->firstcplleak = 1;

    return 0;
} /* end of audfrm */

int eac3_parse_audblk(GetBitContext *gbc, EAC3Context *s, const int blk){
    //int grp, sbnd, n, bin;
    int seg, bnd, ch, i;
    int got_cplchan;
    /* These fields for block switch and dither flags */
    if(s->blkswe)
    {
        for(ch = 1; ch <= s->nfchans; ch++) {
            GET_BITS(s->blksw[ch], gbc, 1);
        }
    }
    else
    {
        for(ch = 1; ch <= s->nfchans; ch++) {
            s->blksw[ch] = 0;
        }
    }
    if(s->dithflage)
    {
        for(ch = 1; ch <= s->nfchans; ch++) {
            GET_BITS(s->dithflag[ch], gbc, 1);
        }
    }
    else
    {
        for(ch = 1; ch <= s->nfchans; ch++) {
            s->dithflag[ch] = 1; /* dither on */
        }
    }
    /* These fields for dynamic range control */
    GET_BITS(s->dynrnge, gbc, 1);
    if(s->dynrnge) {
        GET_BITS(s->dynrng, gbc, 8);
    }
    if(s->acmod == 0x0) /* if 1+1 mode */
    {
        GET_BITS(s->dynrng2e, gbc, 1);
        if(s->dynrng2e) {
            GET_BITS(s->dynrng2, gbc, 8);
        }
    }
    /* These fields for spectral extension strategy information */
    if(blk == 0) {
        s->spxstre = 1;
    }
    else {
        GET_BITS(s->spxstre, gbc, 1);
    }
    if(!blk && !s->spxstre){
        av_log(s, AV_LOG_ERROR, "no spectral extension strategy in first block");
        return -1;
    }
    if(s->spxstre)
    {
        GET_BITS(s->spxinu, gbc, 1);
        if(s->spxinu)
        {
            if(s->acmod == 0x1)
            {
                s->chinspx[0] = 1;
            }
            else
            {
                for(ch = 1; ch <= s->nfchans; ch++) {
                    GET_BITS(s->chinspx[ch], gbc, 1);
                }
            }
            GET_BITS(s->spxstrtf, gbc, 2);
            GET_BITS(s->spxbegf, gbc, 3);
            GET_BITS(s->spxendf, gbc, 3);
            if(s->spxbegf < 6) {
                s->spxbegf += 2;
            }
            else {
                s->spxbegf = s->spxbegf * 2 - 3;
            }
            if(s->spxendf < 3) {
                s->spxendf += 5;
            }
            else {
                s->spxendf = s->spxendf * 2 + 3;
            }
            GET_BITS(s->spxbndstrce, gbc, 1);

            assert(blk || s->spxbndstrce); // TODO default values..

            if(s->spxbndstrce)
            {
                for(bnd = s->spxbegf+1; bnd < s->spxendf; bnd++) {
                    GET_BITS(s->spxbndstrc[bnd], gbc, 1);
                }
            }
            // calculate number of spectral extension bands
            s->nspxbnds = 1;
            s->spxbndsztab[0] = 12;
            for (bnd = s->spxbegf+1; bnd < s->spxendf; bnd ++)
            {
                if (s->spxbndstrc[bnd] == 0)
                {
                    s->spxbndsztab[s->nspxbnds] = 12;
                    s->nspxbnds++;
                }
                else
                {
                    s->spxbndsztab[s->nspxbnds - 1] += 12;
                }
            }
            if(s->nspxbnds >= MAX_SPX_CODES){
                av_log(s, AV_LOG_ERROR, "s->nspxbnds >= MAX_SPX_CODES");
                return -1;
            }
        }
        else /* !spxinu */
        {
            for(ch = 1; ch <= s->nfchans; ch++)
            {
                s->chinspx[ch] = 0;
                s->firstspxcos[ch] = 1;
            }
        }
    }


    /* These fields for spectral extension coordinates */
    if(s->spxinu)
    {
        for(ch = 1; ch <= s->nfchans; ch++)
        {
            if(s->chinspx[ch])
            {
                if(s->firstspxcos[ch])
                {
                    s->spxcoe[ch] = 1;
                    s->firstspxcos[ch] = 0;
                }
                else /* !firstspxcos[ch] */ {
                    GET_BITS(s->spxcoe[ch], gbc, 1);
                }
                if(!blk && !s->spxcoe[ch]){
                    av_log(s, AV_LOG_ERROR, "no spectral extension coordinates in first block");
                    return -1;
                }

                if(s->spxcoe[ch])
                {
                    GET_BITS(s->spxblnd[ch], gbc, 5);
                    GET_BITS(s->mstrspxco[ch], gbc, 2);
                    /* nspxbnds determined from spxbegf, spxendf, and spxbndstrc[ ] */
                    if(s->nspxbnds >= MAX_SPX_CODES){
                        av_log(s, AV_LOG_ERROR, "s->nspxbnds >= MAX_SPX_CODES");
                        return -1;
                    }
                    for(bnd = 0; bnd < s->nspxbnds; bnd++)
                    {
                        GET_BITS(s->spxcoexp[ch][bnd], gbc, 4);
                        GET_BITS(s->spxcomant[ch][bnd], gbc, 2);
                    }
                }
            }
            else /* !chinspx[ch] */
            {
                s->firstspxcos[ch] = 1;
            }
        }
    }
    /* These fields for coupling strategy and enhanced coupling strategy information */
    if(s->cplstre[blk])
    {
        if (s->cplinu[blk])
        {
            GET_BITS(s->ecplinu, gbc, 1);
            assert(!s->ecplinu && "TODO");
            if (s->acmod == 0x2)
            {
                s->chincpl[0] = 1;
                s->chincpl[1] = 1;
            }
            else
            {
                for(ch = 1; ch <= s->nfchans; ch++) {
                    GET_BITS(s->chincpl[ch], gbc, 1);
                }
            }
            if (s->ecplinu == 0) /* standard coupling in use */
            {
                if(s->acmod == 0x2) /* if in 2/0 mode */ {
                    GET_BITS(s->phsflginu, gbc, 1);
                }
                GET_BITS(s->cplbegf, gbc, 4);
                if (s->spxinu == 0) /* if SPX not in use */
                {
                    GET_BITS(s->cplendf, gbc, 4);
                    s->cplendf = s->cplendf + 3;
                }
                else /* SPX in use */
                {
                    s->cplendf = s->spxbegf - 1;
                }

                // calc
                s->cplstrtmant = 37 + (12 * s->cplbegf);
                s->cplendmant = 37 + (12 * (s->cplendf + 3));

                GET_BITS(s->cplbndstrce, gbc, 1);
                if(s->cplbndstrce)
                {
                    for(bnd = s->cplbegf+1; bnd < s->cplendf; bnd++) {
                        GET_BITS(s->cplbndstrc[bnd], gbc, 1);
                    }
                }
                //TODO calc ncplsubnd ?
                s->ncplsubnd =  s->cplendf - s->spxbegf + 1;
                s->ncplbnd = s->ncplsubnd;
                for(bnd = 1; bnd < s->cplendf; bnd++){
                    s->ncplbnd -= s->cplbndstrc[bnd];
                }

            }
            else /* enhanced coupling in use */
            {
                av_log(s, AV_LOG_ERROR,  "enhanced couplin NOT IMPLEMENTED");
                return -1;
#if 0
                GET_BITS(s->ecplbegf, gbc, 4);
                if(s->ecplbegf < 3) {
                    s->ecpl_start_subbnd = s->ecplbegf * 2;
                }
                else if(s->ecplbegf < 13) {
                    s->ecpl_start_subbnd = s->ecplbegf + 2;
                }
                else {
                    s->ecpl_start_subbnd = s->ecplbegf * 2 - 10;
                }
                if (s->spxinu == 0) /* if SPX not in use */
                {
                    GET_BITS(s->ecplendf, gbc, 4);
                    s->ecpl_end_subbnd = s->ecplendf + 7;
                }
                else /* SPX in use */
                {
                    if(s->spxbegf < 6) {
                        s->ecpl_end_subbnd = s->spxbegf + 5;
                    }
                    else {
                        s->ecpl_end_subbnd = s->spxbegf * 2;
                    }
                }
                GET_BITS(s->ecplbndstrce, gbc, 1);
                if (s->ecplbndstrce)
                {
                    for(sbnd = FFMAX(9, s->ecpl_start_subbnd+1);
                            sbnd < s->ecpl_end_subbnd; sbnd++)
                    {
                        GET_BITS(s->ecplbndstrc[sbnd], gbc, 1);
                    }
                }
                //necplbnd = ecpl_end_subbnd - ecpl_start_subbnd;
                //necplbnd -= ecplbndstrc[ecpl_start_subbnd] + ... + ecplbndstrc[ecpl_end_subbnd -1]
                // CALC
                s->necplbnd = s->ecpl_end_subbnd - s->ecpl_start_subbnd;
                for(bnd=s->ecpl_start_subbnd; bnd<s->ecpl_end_subbnd; bnd++){
                    s->necplbnd -= s->ecplbndstrc[bnd];
                }
#endif

            } /* ecplinu[blk] */
        }
        else /* !cplinu[blk] */
        {
            for(ch = 1; ch <= s->nfchans; ch++)
            {
                s->chincpl[ch] = 0;
                s->firstcplcos[ch] = 1;
            }
            s->firstcplleak = 1;
            s->phsflginu = 0;
            s->ecplinu = 0;
        }
    } /* cplstre[blk] */
    /* These fields for coupling coordinates */
    if(s->cplinu[blk])
    {
        av_log(s, AV_LOG_INFO, "NOT TESTED");

        if(s->ecplinu == 0) /* standard coupling in use */
        {
            for(ch = 1; ch <= s->nfchans; ch++)
            {
                if(s->chincpl[ch])
                {
                    if (s->firstcplcos[ch])
                    {
                        s->cplcoe[ch] = 1;
                        s->firstcplcos[ch] = 0;
                    }
                    else /* !firstcplcos[ch] */ {
                        GET_BITS(s->cplcoe[ch], gbc, 1);
                    }
                    if(s->cplcoe[ch])
                    {
                        GET_BITS(s->mstrcplco[ch], gbc, 2);
                        /* ncplbnd derived from cplbegf, cplendf, and cplbndstrc */
                        for(bnd = 0; bnd < s->ncplbnd; bnd++)
                        {
                            GET_BITS(s->cplcoexp[ch][bnd], gbc, 4);
                            GET_BITS(s->cplcomant[ch][bnd], gbc, 4);
                        }
                    } /* cplcoe[ch] */
                }
                else /* ! chincpl[ch] */
                {
                    s->firstcplcos[ch] = 1;
                }
            } /* ch */
            if((s->acmod == 0x2) && s->phsflginu && (s->cplcoe[0] || s->cplcoe[1]))
            {
                for(bnd = 0; bnd < s->ncplbnd; bnd++) {
                    GET_BITS(s->phsflg[bnd], gbc, 1);
                }
            }
        }
        else /* enhanced coupling in use */
        {
            s->firstchincpl = -1;
            GET_BITS(s->ecplangleintrp, gbc, 1);
            for(ch = 1; ch <= s->nfchans; ch++)
            {
                if(s->chincpl[ch])
                {
                    if(s->firstchincpl == -1) {
                        s->firstchincpl = ch;
                    }
                    if(s->firstcplcos[ch])
                    {
                        s->ecplparam1e[ch] = 1;
                        if (ch > s->firstchincpl) {
                            s->ecplparam2e[ch] = 1;
                        }
                        else {
                            s->ecplparam2e[ch] = 0;
                        }
                        s->firstcplcos[ch] = 0;
                    }
                    else /* !firstcplcos[ch] */
                    {
                        GET_BITS(s->ecplparam1e[ch], gbc, 1);
                        if(ch > s->firstchincpl) {
                            GET_BITS(s->ecplparam2e[ch], gbc, 1);
                        }
                        else {
                            s->ecplparam2e[ch] = 0;
                        }
                    }
                    if(s->ecplparam1e[ch])
                    {
                        /* necplbnd derived from ecpl_start_subbnd, ecpl_end_subbnd, and ecplbndstrc */
                        for(bnd = 0; bnd < s->necplbnd; bnd++) {
                            GET_BITS(s->ecplamp[ch][bnd], gbc, 5);
                        }
                    }
                    if(s->ecplparam2e[ch])
                    {
                        /* necplbnd derived from ecpl_start_subbnd, ecpl_end_subbnd, and ecplbndstrc */
                        for(bnd = 0; bnd < s->necplbnd; bnd++)
                        {
                            GET_BITS(s->ecplangle[ch][bnd], gbc, 6);
                            GET_BITS(s->ecplchaos[ch][bnd], gbc, 3);
                        }
                    }
                    if(ch > s->firstchincpl) {
                        GET_BITS(s->ecpltrans[ch], gbc, 1);
                    }
                }
                else /* !chincpl[ch] */
                {
                    s->firstcplcos[ch] = 1;
                }
            } /* ch */
        } /* ecplinu[blk] */
    } /* cplinu[blk] */
    /* These fields for rematrixing operation in the 2/0 mode */
    if(s->acmod == 0x2) /* if in 2/0 mode */
    {
        if (blk == 0) {
            s->rematstr = 1;
        }
        else {
            GET_BITS(s->rematstr, gbc, 1);
        }
        if(s->rematstr)
        {
            /* nrematbnds determined from cplinu, ecplinu, spxinu, cplbegf, ecplbegf and spxbegf
             * TODO: how ? */
            av_log(s, AV_LOG_ERROR,  "NOT IMPLEMENTED");
            return -1;
            /*for(bnd = 0; bnd < s->nrematbnds; bnd++) {
                GET_BITS(s->rematflg[bnd], gbc, 1);
            }*/
        }
    }
    /* This field for channel bandwidth code */
    for(ch = 1; ch <= s->nfchans; ch++)
    {
        if(!blk && s->chexpstr[blk][ch]==EXP_REUSE){
            av_log(s, AV_LOG_ERROR,  "no channel exponent strategy in first block");
            return -1;
        }
        if(s->chexpstr[blk][ch] != EXP_REUSE)
        {
            if((!s->chincpl[ch]) && (!s->chinspx[ch])) {
                GET_BITS(s->chbwcod[ch], gbc, 6);
              if(s->chbwcod[ch] > 60)
                  return 1;
            }
        }
    }

    // calc
    for(ch = 1; ch<=s->nfchans; ch++){
        int grpsize = 3 << (s->chexpstr[blk][ch] - 1);

        if(s->chincpl[ch])
            s->endmant[ch] = s->cplstrtmant; /* channel is coupled */
        else
            s->endmant[ch] = ((s->chbwcod[ch] + 12) * 3) + 37; /* (ch is not coupled) */

        s->nchgrps[ch] = (s->endmant[ch] + grpsize - 4) / grpsize;
    }

    /* These fields for exponents */
    if(s->cplinu[blk]) /* exponents for the coupling channel */
    {
        if(s->cplexpstr[blk] != EXP_REUSE)
        {
            GET_BITS(s->cplabsexp, gbc, 4);
            /* ncplgrps derived from cplbegf, ecplbegf, cplendf, ecplendf, and cplexpstr */
            /* how... ? */
            av_log(s, AV_LOG_ERROR,  "NOT IMPLEMENTED");
            return -1;
            /*for(grp = 0; grp< s->ncplgrps; grp++) {
                GET_BITS(s->cplexps[grp], gbc, 7);
            }*/
        }
    }
    for(ch = 1; ch <= s->nfchans; ch++) /* exponents for full bandwidth channels */
    {
        if(!blk && !s->chexpstr[blk][ch] == EXP_REUSE){
            av_log(s, AV_LOG_ERROR,  "no channel exponent strategy in first block");
            return -1;
        }
        if(s->chexpstr[blk][ch] != EXP_REUSE)
        {
            GET_BITS(s->dexps[ch][0], gbc, 4);

            ff_ac3_decode_exponents(gbc, s->chexpstr[blk][ch], s->nchgrps[ch], s->dexps[ch][0],
                    s->dexps[ch] + 1);


            GET_BITS(s->gainrng[ch], gbc, 2);
        }
    }
    if(s->lfeon) /* exponents for the low frequency effects channel */
    {
        if(s->lfeexpstr[blk] != EXP_REUSE)
        {
            ch = s->lfe_channel;
            GET_BITS(s->dexps[ch][0], gbc, 4);
            s->nlfegrps = 2;
            ff_ac3_decode_exponents(gbc, s->lfeexpstr[blk], s->nlfegrps, s->dexps[ch][0],
                    s->dexps[ch] + 1);
        }
    }
    /* These fields for bit-allocation parametric information */
    if(s->bamode)
    {
        GET_BITS(s->baie, gbc, 1);
        assert(s->baie || blk);
        if(s->baie)
        {
            GET_BITS(s->sdcycod, gbc, 2);
            GET_BITS(s->fdcycod, gbc, 2);
            GET_BITS(s->sgaincod, gbc, 2);
            GET_BITS(s->dbpbcod, gbc, 2);
            GET_BITS(s->floorcod, gbc, 3);
        }
    }
    else
    {
        s->sdcycod = 0x2;
        s->fdcycod = 0x1;
        s->sgaincod = 0x1;
        s->dbpbcod = 0x2;
        s->floorcod = 0x7;
    }

    s->sdecay = ff_sdecaytab[s->sdcycod];   /* Table 7.6 */
    s->fdecay = ff_fdecaytab[s->fdcycod];   /* Table 7.7 */
    s->sgain = ff_sgaintab[s->sgaincod];    /* Table 7.8 */
    s->dbknee = ff_dbkneetab[s->dbpbcod];   /* Table 7.9 */
    s->floor = ff_floortab[s->floorcod];    /* Table 7.10 */

    if(s->snroffststr == 0x0)
    {
        if(s->cplinu[blk]) {
            s->cplfsnroffst = s->frmfsnroffst;
        }
        s->csnroffst = s->frmcsnroffst;
        for(ch = 1; ch <= s->nfchans; ch++) {
            s->fsnroffst[ch] = s->frmfsnroffst;
        }
        if(s->lfeon) {
            s->lfefsnroffst = s->frmfsnroffst;
        }
    }
    else
    {
        av_log(s, AV_LOG_INFO, "NOT TESTED");
        if(blk == 0) {
            s->snroffste = 1;
        }
        else {
            GET_BITS(s->snroffste, gbc, 1);
        }
        if(s->snroffste)
        {
            GET_BITS(s->csnroffst, gbc, 6);
            if(s->snroffststr == 0x1)
            {
                GET_BITS(s->blkfsnroffst, gbc, 4);
                s->cplfsnroffst = s->blkfsnroffst;
                for(ch = 1; ch <= s->nfchans; ch++) {
                    s->fsnroffst[ch] = s->blkfsnroffst;
                }
                s->lfefsnroffst = s->blkfsnroffst;
            }
            else if(s->snroffststr == 0x2)
            {
                if(s->cplinu[blk]) {
                    GET_BITS(s->cplfsnroffst, gbc, 4);
                }
                for(ch = 1; ch <= s->nfchans; ch++) {
                    GET_BITS(s->fsnroffst[ch], gbc, 4);
                }
                if(s->lfeon) {
                    GET_BITS(s->lfefsnroffst, gbc, 4);
                }
            }
        }
    }

    if(s->lfeon){
        s->fsnroffst[s->lfe_channel] = s->lfefsnroffst;
    }

    if(s->frmfgaincode) {
        GET_BITS(s->fgaincode, gbc, 1);
    }
    else {
        s->fgaincode = 0;
    }
    if(s->fgaincode)
    {
        if(s->cplinu[blk]) {
            GET_BITS(s->cplfgaincod, gbc, 3);
        }
        for(ch = 1; ch <= s->nfchans; ch++) {
            GET_BITS(s->fgaincod[ch], gbc, 3);
        }
        if(s->lfeon) {
            GET_BITS(s->lfefgaincod, gbc, 3);
        }
    }
    else
    {
        if(s->cplinu[blk]) {
            s->cplfgaincod = 0x4;
        }
        for(ch = 1; ch <= s->nfchans; ch++) {
            s->fgaincod[ch] = 0x4;
        }
        if(s->lfeon) {
            s->lfefgaincod = 0x4;
        }
    }
    if(s->lfeon){
        s->fgaincod[s->lfe_channel] = s->lfefgaincod;
    }
    if(s->strmtyp == 0x0)
    {
        GET_BITS(s->convsnroffste, gbc, 1);
        if(s->convsnroffste) {
            GET_BITS(s->convsnroffst, gbc, 10);
        }
    }
    if(s->cplinu[blk])
    {
        if (s->firstcplleak)
        {
            s->cplleake = 1;
            s->firstcplleak = 0;
        }
        else /* !firstcplleak */
        {
            GET_BITS(s->cplleake, gbc, 1);
        }
        if(s->cplleake)
        {
            GET_BITS(s->cplfleak, gbc, 3);
            GET_BITS(s->cplsleak, gbc, 3);
        }
    }
    /* Delta bit allocation information */
    if(s->dbaflde)
    {
        av_log(s, AV_LOG_INFO, "NOT TESTED");
        GET_BITS(s->deltbaie, gbc, 1);
        if(s->deltbaie)
        {
            if(s->cplinu[blk]) {
                GET_BITS(s->cpldeltbae, gbc, 2);
            }
            for(ch = 1; ch <= s->nfchans; ch++) {
                GET_BITS(s->deltbae[ch], gbc, 2);
            }
            if(s->cplinu[blk])
            {
                //TODO
                if(s->cpldeltbae==DBA_NEW)
                {
                    GET_BITS(s->cpldeltnseg, gbc, 3);
                    for(seg = 0; seg <= s->cpldeltnseg; seg++)
                    {
                        GET_BITS(s->cpldeltoffst[seg], gbc, 5);
                        GET_BITS(s->cpldeltlen[seg], gbc, 4);
                        GET_BITS(s->cpldeltba[seg], gbc, 3);
                    }
                }
            }
            for(ch = 1; ch <= s->nfchans; ch++)
            {
                if(s->deltbae[ch]==DBA_NEW)
                {
                    GET_BITS(s->deltnseg[ch], gbc, 3);
                    for(seg = 0; seg <= s->deltnseg[ch]; seg++)
                    {
                        GET_BITS(s->deltoffst[ch][seg], gbc, 5);
                        GET_BITS(s->deltlen[ch][seg], gbc, 4);
                        GET_BITS(s->deltba[ch][seg], gbc, 3);
                    }
                }
            }
        } /* if(s->deltbaie) */
    }/* if(s->dbaflde) */


    /* These fields for inclusion of unused dummy data */
    if(s->skipflde)
    {
        GET_BITS(s->skiple, gbc, 1);
        if(s->skiple)
        {
            GET_BITS(s->skipl, gbc, 9);
            // TODO skip :)
            for(i=0; i<s->skipl; i++)
                GET_BITS(s->skipfld, gbc, 8);
        }
    }

    /* run bit allocation */
    if(s->cplinu[blk]) {
        av_log(s, AV_LOG_ERROR,  "NOT IMPLEMENTED");
        return -1;
    }

    for(ch = 1; ch<=s->nfchans+s->lfeon; ch++) {
        int start=0, end=0;
        end = s->endmant[ch];

        ff_ac3_bit_alloc_calc_psd((int8_t *)s->dexps[ch], start, end,
                s->psd[ch], s->bndpsd[ch]);


        // TODO hmm... :)
        s->bit_alloc_params.fscod = s->fscod;
        s->bit_alloc_params.halfratecod = 0; // TODO
        s->bit_alloc_params.sdecay = s->sdecay;
        s->bit_alloc_params.fdecay = s->fdecay;
        s->bit_alloc_params.sgain = s->sgain;
        s->bit_alloc_params.dbknee = s->dbknee;
        s->bit_alloc_params.floor = s->floor;
        s->bit_alloc_params.cplfleak = s->cplfleak;
        s->bit_alloc_params.cplsleak = s->cplsleak;

        {
            int fgain = ff_fgaintab[s->fgaincod[ch]];
            if(ch == s->lfe_channel){
                s->deltbae[ch] = DBA_NONE;
            }
            ff_ac3_bit_alloc_calc_mask(&s->bit_alloc_params,
                    s->bndpsd[ch], start, end, fgain,
                    (ch == s->lfe_channel),
                    s->deltbae[ch], s->deltnseg[ch],
                    s->deltoffst[ch], s->deltlen[ch],
                    s->deltba[ch], s->mask[ch]);
        }

        {
            int snroffst = (((s->csnroffst - 15) << 4) + s->fsnroffst[ch]) << 2;
            //av_log(NULL, AV_LOG_INFO, "s->csnroffst=%i s->fsnroffst=%i snroffst = %i\n",
             //       s->csnroffst, s->fsnroffst[ch], snroffst);
            ff_ac3_bit_alloc_calc_bap(s->mask[ch], s->psd[ch], start, end,
                    snroffst, s->bit_alloc_params.floor,
                    s->bap[ch]);
        }


    }


    /* These fields for quantized mantissa values */

    got_cplchan = 0;
    ff_ac3_get_transform_coeffs(gbc, s->bap, s->dexps, s->nfchans+s->lfeon, s->chincpl, s->dithflag, s->transform_coeffs, s->strtmant, s->endmant, &s->dith_state);

    for(ch = 1; ch <= s->nfchans; ch++)
    {
        if(s->chahtinu[ch] == 0)
        {
            // ff_ac3_get_transform
        }
        else if(s->chahtinu[ch] == 1)
        {
            av_log(s, AV_LOG_ERROR,  "AHT NOT IMPLEMENTED");
            return -1;

#if 0
            GET_BITS(s->chgaqmod[ch], gbc, 2);
            if((s->chgaqmod[ch] > 0x0) && (s->chgaqmod[ch] < 0x3) )
            {
                for(n = 0; n < s->chgaqsections[ch]; n++) { // TODO chgaqsections ?
                    GET_BITS(s->chgaqgain[ch][n], gbc, 1);
                }
            }
            else if(s->chgaqmod[ch] == 0x3)
            {
                for(n = 0; n < s->chgaqsections[ch]; n++) { //TODO chgaqsections ?
                    GET_BITS(s->chgaqgain[ch][n], gbc, 5);
                }
            }
            for(bin = 0; bin < s->nchmant[ch]; bin++) // TODO nchmant ?
            {
                if(s->chgaqbin[ch][bin]) // TODO chgaqbin ?
                {
                    for(n = 0; n < 6; n++) {
                        GET_BITS(s->pre_chmant[n][ch][bin], gbc, (0-16)); // TODO 0-16 :]
                    }
                }
                else {
                    GET_BITS(s->pre_chmant[0][ch][bin], gbc, (0-9)); //  TODO 0-9 :]
                }
            }
            s->chahtinu[ch] = -1; /* AHT info for this frame has been read ? do not read again */
#endif
        }
        if(s->cplinu[blk] && s->chincpl[ch] && !got_cplchan)
        {
            av_log(s, AV_LOG_ERROR,  "NOT IMPLEMENTED");
            return -1;
#if 0
            if(s->cplahtinu == 0)
            {
                for(bin = 0; bin < s->ncplmant; bin++) { // TODO ncplmant ?
                    GET_BITS(s->cplmant[bin], gbc, (0-16)); // TODO 0-16 :]
                }
                got_cplchan = 1;
            }
            else if(s->cplahtinu == 1)
            {
                GET_BITS(s->cplgaqmod, gbc, 2);
                if((s->cplgaqmod > 0x0) && (s->cplgaqmod < 0x3) )
                {
                    for(n = 0; n < s->cplgaqsections; n++) { // TODO cplgaqsections?
                        GET_BITS(s->cplgaqgain[n], gbc, 1);
                    }
                }
                else if(s->cplgaqmod == 0x3)
                {
                    for(n = 0; n < s->cplgaqsections; n++) {
                        GET_BITS(s->cplgaqgain[n], gbc, 5);
                    }
                }
                for(bin = 0; bin < s->ncplmant; bin++) // TODO ncplmant?
                {
                    if(s->cplgaqbin[bin])
                    {
                        for(n = 0; n < 6; n++) {
                            GET_BITS(s->pre_cplmant[n][bin], gbc, (0-16)); // TODO 0-16 :]
                        }
                    }
                    else {
                        GET_BITS(s->pre_cplmant[0][bin], gbc, (0-9));
                    }
                }
                got_cplchan = 1;
                s->cplahtinu = -1; /* AHT info for this frame has been read ? do not read again */
            }
            else {
                got_cplchan = 1;
            }
#endif
        }
    }
    if(s->lfeon) /* mantissas of low frequency effects channel */
    {
        if(s->lfeahtinu == 0)
        {
            //ff_ac3_get_transform
        }
        else if(s->lfeahtinu == 1)
        {
            av_log(s, AV_LOG_ERROR,  "NOT IMPLEMENTED");
            return -1;
#if 0
            assert(0 && "TODO: AHT");
            GET_BITS(s->lfegaqmod, gbc, 2);
            if( (s->lfegaqmod > 0x0) && (s->lfegaqmod < 0x3) )
            {
                for(n = 0; n < s->lfegaqsections; n++) { // TODO  lfegaqsections?
                    GET_BITS(s->lfegaqgain[n], gbc, 1);
                }
            }
            else if(s->lfegaqmod == 0x3)
            {
                for(n = 0; n < s->lfegaqsections; n++) { // TODO
                    GET_BITS(s->lfegaqgain[n], gbc, 5);
                }
            }
            for(bin = 0; bin < s->nlfemant; bin++)
            {
                if(s->lfegaqbin[bin])
                {
                    for(n = 0; n < 6; n++) {
                        GET_BITS(s->pre_lfemant[n][bin], gbc, (0-16)); // TODO 0-16 :]
                    }
                }
                else {
                    GET_BITS(s->pre_lfemant[0][bin], gbc, (0-9)); // TODO 0-9 :]
                }
            }
            s->lfeahtinu = -1; /* AHT info for this frame has been read ? do not read again */
#endif
        }
    }
    return 0;
}

int eac3_parse_auxdata(GetBitContext *gbc, EAC3Context *s){
    // TODO
    return 0;
}
