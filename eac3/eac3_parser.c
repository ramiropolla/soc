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

#define DEBUG

#include "avcodec.h"
#include "bitstream.h"
#include "ac3.h"
#include "random.h"

#undef DEBUG

#ifdef DEBUG
#define GET_BITS(a, gbc, n) {a = get_bits(gbc, n); av_log(NULL, AV_LOG_INFO, "%s: %i\n", __STRING(a), a);}
#define GET_SBITS(a, gbc, n) {a = get_sbits(gbc, n); av_log(NULL, AV_LOG_INFO, "%s: %i\n", __STRING(a), a);}
#else
#define GET_BITS(a, gbc, n) a = get_bits(gbc, n)
#define GET_SBITS(a, gbc, n) a = get_sbits(gbc, n)
#endif //DEBUG

#include "eac3.h"
#include "ac3dec.h"

int ff_eac3_parse_syncinfo(GetBitContext *gbc, EAC3Context *s){
    GET_BITS(s->syncword, gbc, 16);
    return 0;
}

int ff_eac3_parse_bsi(GetBitContext *gbc, EAC3Context *s){
    int i, blk;

    GET_BITS(s->strmtyp, gbc, 2);
    GET_BITS(s->substreamid, gbc, 3);
    GET_BITS(s->frmsiz, gbc, 11);
    GET_BITS(s->fscod, gbc, 2);
    if(s->fscod == 0x3)
    {
        av_log(s->avctx, AV_LOG_ERROR, "NOT IMPLEMENTED");
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
    s->lfe_channel = s->ntchans+1;
    if(s->lfeon){
        s->strtmant[s->lfe_channel] = 0;
        s->endmant[s->lfe_channel] = 7;
        s->nchgrps[s->lfe_channel] = 2;

        s->ntchans ++ ;
    }

    GET_BITS(s->bsid, gbc, 5);
    if(s->bsid < 11 || s->bsid > 16){
        av_log(s->avctx, AV_LOG_ERROR, "bsid should be between 11 and 16");
        return -1;
    }

    GET_BITS(s->dialnorm[0], gbc, 5);
    if(get_bits1(gbc)) {
        GET_BITS(s->compr, gbc, 8);
    }else{
        //TODO default compr
    }
    if(s->acmod == 0x0) /* if 1+1 mode (dual mono, so some items need a second value) */
    {
        GET_BITS(s->dialnorm[1], gbc, 5);
        if(get_bits1(gbc)) {
            GET_BITS(s->compr2, gbc, 8);
        }else{
            //TODO default compr2
        }
    }
    if(s->strmtyp == 0x1) /* if dependent stream */
    {
        if(get_bits1(gbc)) {
            GET_BITS(s->chanmap, gbc, 16);
        }else{
            //TODO default channel map
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
            if(get_bits1(gbc)) {
                GET_BITS(s->pgmscl, gbc, 6);
            }else{
                //TODO program scale factor = 0dB
            }
            if(s->acmod == 0x0) /* if 1+1 mode (dual mono, so some items need a second value) */
            {
                if(get_bits1(gbc)) {
                    GET_BITS(s->pgmscl2, gbc, 6);
                }else{
                    //TODO program scale factor 2 = 0dB
                }
            }
            if(get_bits1(gbc)) {
                GET_BITS(s->extpgmscl, gbc, 6);
            }
            GET_BITS(s->mixdef, gbc, 2);
            if(s->mixdef == 0x1) /* mixing option 2 */ {
                skip_bits(gbc, 5);
            }
            else if(s->mixdef == 0x2) /* mixing option 3 */ {
                skip_bits(gbc, 12);
            }
            else if(s->mixdef == 0x3) /* mixing option 4 */
            {
                GET_BITS(s->mixdeflen, gbc, 5);
                skip_bits(gbc, 8*(s->mixdeflen+2));
            }
            if(s->acmod < 0x2) /* if mono or dual mono source */
            {
                if(get_bits1(gbc)) {
                    GET_BITS(s->paninfo, gbc, 14);
                }
                if(s->acmod == 0x0) /* if 1+1 mode (dual mono, so some items need a second value) */
                {
                    if(get_bits1(gbc)) {
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
        skip_bits1(gbc); //converter synchronization flag
    }
    if(s->strmtyp == 0x2) /* if bit stream converted from AC-3 */
    {
        if(s->numblkscod == 0x3 || get_bits1(gbc)) /* 6 blocks per frame */ {
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


int ff_eac3_parse_audfrm(GetBitContext *gbc, EAC3Context *s){
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
            for(ch = !s->cplinu[blk]; ch <= s->nfchans; ch++) {
                GET_BITS(s->chexpstr[blk][ch], gbc, 2); // TODO
            }
        }
    }
    else
    {
        /* cplexpstr[blk] and chexpstr[blk][ch] derived from table lookups. see Table E2.14 */
        if( (s->acmod > 0x1) && (s->ncplblks > 0) ) {
            GET_BITS(s->frmchexpstr[CPL_CH], gbc, 5);
            for(blk=0; blk<6; blk++){
                s->chexpstr[blk][CPL_CH] = ff_eac3_frm_expstr[s->frmchexpstr[0]][blk];
            }
        }
        for(ch = 1; ch <= s->nfchans; ch++) {
            GET_BITS(s->frmchexpstr[ch], gbc, 5);
            for(blk=0; blk<6; blk++){
                s->chexpstr[blk][ch] = ff_eac3_frm_expstr[s->frmchexpstr[ch]][blk];
            }
        }
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
        if(s->numblkscod == 0x3 || get_bits1(gbc)){
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
        av_log(s->avctx, AV_LOG_ERROR, "AHT NOT IMPLEMENTED");
        return -1;

#if 0
        /* AHT is only available in 6 block mode (numblkscod ==0x3) */

        s->ncplregs = 0;
        for(blk = 0; blk < 6; blk++){
            if(s->cplstre[blk]==1 || s->chexpstr[blk][CPL_CH] != EXP_REUSE)
                s->ncplregs++;
        }
        s->nchregs[0] = s->ncplregs;

        for(ch = 1; ch <= s->nfchans+s->lfeon; ch++){
            s->nchregs[ch] = 0;
            for(blk = 0; blk < 6; blk++){
                if(s->chexpstr[blk][ch] != EXP_REUSE)
                    s->nchregs[ch]++;
            }
        }

        /*
        s->nlferegs = 0;
        for(blk = 0; blk < 6; blk++){
            if(s->lfeexpstr[blk] != EXP_REUSE)
                s->nlferegs++;
        }
        */

        if( (s->ncplblks == 6) && (s->ncplregs ==1) ) {
            GET_BITS(s->cplahtinu, gbc, 1);
        }
        else {
            s->cplahtinu = 0;
        }
        s->chahtinu[0] = s->cplahtinu;

        for(ch = 1; ch <= s->nfchans; ch++)
        {
            // nchregs derived from chexpstr ? see Section E3.3.2
            if(s->nchregs[ch] == 1) {
                GET_BITS(s->chahtinu[ch], gbc, 1);
            }
            else {
                s->chahtinu[ch] = 0;
            }
        }

        if(s->lfeon)
        {
            if(s->nchregs[s->lfe_channel] == 1) {
                GET_BITS(s->lfeahtinu, gbc, 1);
            }
            else {
                s->lfeahtinu = 0;
            }
            s->chahtinu[s->lfe_channel] = s->lfeahtinu;
        }
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
        av_log(s->avctx, AV_LOG_INFO, "nblkstrtbits = %i\n", nblkstrtbits);
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

int ff_eac3_parse_audblk(GetBitContext *gbc, EAC3Context *s, const int blk){
    //int grp, sbnd, n, bin;
    int seg, bnd, ch, i;
    int got_cplchan;
    mant_groups m;

    m.b1ptr = m.b2ptr = m.b4ptr = 3;

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
    s->dithflag[CPL_CH] = s->dithflag[s->lfe_channel] = 0;
    /* These fields for dynamic range control */

    for(i = 0; i < (s->acmod?1:2); i++){
        GET_BITS(s->dynrnge[i], gbc, 1);
        if(s->dynrnge[i]){
            GET_BITS(s->dynrng[i], gbc, 8);
        }else{
            if(blk==0){
                s->dynrng[i] = 0;
            }
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
        av_log(s->avctx, AV_LOG_ERROR, "no spectral extension strategy in first block");
        return -1;
    }
    if(s->spxstre)
    {
        GET_BITS(s->spxinu, gbc, 1);
        if(s->spxinu)
        {
            if(s->acmod == 0x1)
            {
                s->chinspx[1] = 1;
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
                av_log(s->avctx, AV_LOG_ERROR, "s->nspxbnds >= MAX_SPX_CODES");
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
                    av_log(s->avctx, AV_LOG_ERROR, "no spectral extension coordinates in first block");
                    return -1;
                }

                if(s->spxcoe[ch])
                {
                    GET_BITS(s->spxblnd[ch], gbc, 5);
                    GET_BITS(s->mstrspxco[ch], gbc, 2);
                    /* nspxbnds determined from spxbegf, spxendf, and spxbndstrc[ ] */
                    if(s->nspxbnds >= MAX_SPX_CODES){
                        av_log(s->avctx, AV_LOG_ERROR, "s->nspxbnds >= MAX_SPX_CODES");
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
                s->chincpl[1] = 1;
                s->chincpl[2] = 1;
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

                av_log(s->avctx, AV_LOG_DEBUG, "cplbegf=%i cplendf=%i\n", s->cplbegf, s->cplendf);
                // calc
                s->strtmant[CPL_CH] = 37 + (12 * s->cplbegf);
                s->endmant[CPL_CH] = 37 + (12 * (s->cplendf));
                if(s->strtmant[CPL_CH] > s->endmant[CPL_CH]){
                    av_log(s->avctx, AV_LOG_ERROR, "cplstrtmant > cplendmant [blk=%i]\n", blk);
                    return -1;
                }

                GET_BITS(s->cplbndstrce, gbc, 1);
                if(s->cplbndstrce)
                {
                    for(bnd = s->cplbegf+1; bnd < s->cplendf; bnd++) {
                        GET_BITS(s->cplbndstrc[bnd], gbc, 1);
                    }
                }else if(!blk){
                    for(bnd = 0; bnd < 18; bnd++)
                        s->cplbndstrc[bnd] = ff_eac3_defcplbndstrc[bnd];
                }
                //TODO calc ncplsubnd ?
                s->ncplsubnd =  s->cplendf - s->cplbegf;
                s->ncplbnd = s->ncplsubnd;
                assert(!s->cplbndstrc[0]);
                for(bnd = s->cplbegf+1; bnd < s->cplendf; bnd++){
                    s->ncplbnd -= s->cplbndstrc[bnd];
                }
            }
            else /* enhanced coupling in use */
            {
                av_log(s->avctx, AV_LOG_ERROR,  "enhanced couplin NOT IMPLEMENTED");
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
//        av_log(s->avctx, AV_LOG_INFO, "NOT TESTED CPLINU\n");

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
                        int cplcoexp, cplcomant, mstrcplco;
                        GET_BITS(mstrcplco, gbc, 2);
                        mstrcplco = 3 * mstrcplco;
                        /* ncplbnd derived from cplbegf, cplendf, and cplbndstrc */
                        for(bnd = 0; bnd < s->ncplbnd; bnd++)
                        {
                            GET_BITS(cplcoexp, gbc, 4);
                            GET_BITS(cplcomant, gbc, 4);
                            if(cplcoexp==15)
                                s->cplco[ch][bnd] = cplcomant / 16.0f;
                            else
                                s->cplco[ch][bnd] = (cplcomant + 16.0f) / 32.0f;
                            s->cplco[ch][bnd] *=  ff_ac3_scale_factors[cplcoexp + mstrcplco];
                        }
                    } /* cplcoe[ch] */
                    else{
                        if(!blk){
                            av_log(s->avctx, AV_LOG_ERROR,  "no coupling coordinates in first block");
                            return -1;
                        }
                    }
                }
                else /* ! chincpl[ch] */
                {
                    s->firstcplcos[ch] = 1;
                }
            } /* ch */
            if((s->acmod == 0x2) && s->phsflginu && (s->cplcoe[1] || s->cplcoe[2]))
            {
                for(bnd = 0; bnd < s->ncplbnd; bnd++) {
                    GET_BITS(s->phsflg[bnd], gbc, 1);
                }
            }
        }
        else /* enhanced coupling in use */
        {
#if 0
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
#endif
        } /* ecplinu[blk] */
    } /* cplinu[blk] */
    /* These fields for rematrixing operation in the 2/0 mode */
    if(s->acmod == 0x2) /* if in 2/0 mode */
    {
        if (blk == 0 || get_bits1(gbc)) {
            s->rematstr = 1;
        }
        if(s->rematstr)
        {
            /* nrematbnds determined from cplinu, ecplinu, spxinu, cplbegf, ecplbegf and spxbegf
             * TODO XXX (code from AC-3) */
            s->nrematbnds = 4;
            if(s->cplinu && s->cplbegf <= 2)
                s->nrematbnds -= 1 + (s->cplbegf == 0);
            for(bnd = 0; bnd < s->nrematbnds; bnd++) {
                GET_BITS(s->rematflg[bnd], gbc, 1);
            }
        }
    }
    /* This field for channel bandwidth code */
    for(ch = 1; ch <= s->nfchans; ch++)
    {
        if(!blk && s->chexpstr[blk][ch]==EXP_REUSE){
            av_log(s->avctx, AV_LOG_ERROR,  "no channel exponent strategy in first block");
            return -1;
        }
        if(s->chexpstr[blk][ch] != EXP_REUSE)
        {
            if((!s->chincpl[ch]) && (!s->chinspx[ch])) {
                GET_BITS(s->chbwcod[ch], gbc, 6);
              if(s->chbwcod[ch] > 60){
                  av_log(s->avctx, AV_LOG_ERROR, "s->chbwcod[ch] > 60\n");
                  return -1;
              }
            }
        }
    }

    // calc
    for(ch = 1; ch<=s->nfchans; ch++){
        int grpsize = 3 << (s->chexpstr[blk][ch] - 1);

        s->strtmant[ch] = 0;
        if(s->chincpl[ch])
            s->endmant[ch] = s->strtmant[CPL_CH]; /* channel is coupled */
        else
            s->endmant[ch] = ((s->chbwcod[ch] + 12) * 3) + 37; /* (ch is not coupled) */

        s->nchgrps[ch] = (s->endmant[ch] + grpsize - 4) / grpsize;
    }

    /* These fields for exponents */
    if(s->cplinu[blk]) /* exponents for the coupling channel */
    {
        if(s->chexpstr[blk][CPL_CH] != EXP_REUSE)
        {
            /* ncplgrps derived from cplbegf, ecplbegf, cplendf, ecplendf, and cplexpstr */
            /* TODO add support for enhanced coupling */
            switch(s->chexpstr[blk][CPL_CH]){
                case EXP_D15:
                    s->nchgrps[CPL_CH] = (s->endmant[CPL_CH] - s->strtmant[CPL_CH])/3;
                    break;
                case EXP_D25:
                    s->nchgrps[CPL_CH] = (s->endmant[CPL_CH] - s->strtmant[CPL_CH])/6;
                    break;
                case EXP_D45:
                    s->nchgrps[CPL_CH] = (s->endmant[CPL_CH] - s->strtmant[CPL_CH])/12;
                    break;
            }
            GET_BITS(s->cplabsexp, gbc, 4);
/*            for(grp = 0; grp < s->nchgrps[CPL_CH]; grp++) {
                GET_BITS(s->cplexps[grp], gbc, 7);
            }*/
            ff_ac3_decode_exponents(gbc, s->chexpstr[blk][CPL_CH], s->nchgrps[CPL_CH],
                    s->cplabsexp<<1, s->dexps[CPL_CH] + s->strtmant[CPL_CH]);
        }
    }
    for(ch = 1; ch <= s->nfchans; ch++) /* exponents for full bandwidth channels */
    {
        if(!blk && s->chexpstr[blk][ch] == EXP_REUSE){
            av_log(s->avctx, AV_LOG_ERROR,  "no channel exponent strategy in first block");
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
        if(!blk && !s->baie){
            av_log(s->avctx, AV_LOG_ERROR, "no bit allocation information in first block\n");
            return -1;
        }
        if(s->baie)
        {
            s->bit_alloc_params.sdecay = ff_sdecaytab[get_bits(gbc, 2)];   /* Table 7.6 */
            s->bit_alloc_params.fdecay = ff_fdecaytab[get_bits(gbc, 2)];   /* Table 7.7 */
            s->bit_alloc_params.sgain = ff_sgaintab[get_bits(gbc, 2)];     /* Table 7.8 */
            s->bit_alloc_params.dbknee = ff_dbkneetab[get_bits(gbc, 2)];   /* Table 7.9 */
            s->bit_alloc_params.floor = ff_floortab[get_bits(gbc, 3)];     /* Table 7.10 */
        }
    }
    else
    {
        s->bit_alloc_params.sdecay = ff_sdecaytab[0x2];   /* Table 7.6 */
        s->bit_alloc_params.fdecay = ff_fdecaytab[0x1];   /* Table 7.7 */
        s->bit_alloc_params.sgain = ff_sgaintab[0x1];     /* Table 7.8 */
        s->bit_alloc_params.dbknee = ff_dbkneetab[0x2];   /* Table 7.9 */
        s->bit_alloc_params.floor = ff_floortab[0x7];     /* Table 7.10 */
    }

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
        av_log(s->avctx, AV_LOG_INFO, "NOT TESTED");
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
    s->fsnroffst[CPL_CH] = s->cplfsnroffst;

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
        for(ch = !s->cplinu[blk]; ch <= s->nfchans+s->lfeon; ch++)
            GET_BITS(s->fgaincod[ch], gbc, 3);
    }
    else
    {
        for(ch = !s->cplinu[blk]; ch <= s->nfchans+s->lfeon; ch++)
            s->fgaincod[ch] = 0x4;
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
            GET_BITS(s->bit_alloc_params.cplfleak, gbc, 3);
            GET_BITS(s->bit_alloc_params.cplsleak, gbc, 3);
        }
    }
    /* Delta bit allocation information */
    if(s->dbaflde)
    {
        av_log(s->avctx, AV_LOG_INFO, "NOT TESTED");
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
    else{
        if(!blk){
            for(ch=0; ch<=s->nfchans+s->lfe_channel; ch++){
                s->deltbae[ch] = DBA_NONE;
            }
        }
    }


    /* These fields for inclusion of unused dummy data */
    if(s->skipflde)
    {
        if(get_bits1(gbc)){
            int skipl = get_bits(gbc, 9);
            while(skipl--) skip_bits(gbc, 8);
        }
    }

    /* run bit allocation */
    /*if(s->cplinu[blk]) {
        av_log(s->avctx, AV_LOG_ERROR,  "NOT IMPLEMENTED (cplinu && run bit allocation)");
        return -1;
    }*/

    for(ch = !s->cplinu[blk]; ch<=s->nfchans+s->lfeon; ch++) {

        int start=0, end=0;
        start = s->strtmant[ch];
        end = s->endmant[ch];

        ff_ac3_bit_alloc_calc_psd((int8_t *)s->dexps[ch], start, end,
                s->psd[ch], s->bndpsd[ch]);


        // TODO hmm... :)
        s->bit_alloc_params.fscod = s->fscod;
        s->bit_alloc_params.halfratecod = 0; // TODO

        {
            int fgain = ff_fgaintab[s->fgaincod[ch]];
            ff_ac3_bit_alloc_calc_mask(&s->bit_alloc_params,
                    s->bndpsd[ch], start, end, fgain,
                    (ch == s->lfe_channel),
                    s->deltbae[ch], s->deltnseg[ch],
                    s->deltoffst[ch], s->deltlen[ch],
                    s->deltba[ch], s->mask[ch]);
        }

        {
            int snroffst = (((s->csnroffst - 15) << 4) + s->fsnroffst[ch]) << 2;
            //av_log(s->avctx, AV_LOG_INFO, "s->csnroffst=%i s->fsnroffst=%i snroffst = %i\n",
             //       s->csnroffst, s->fsnroffst[ch], snroffst);

            // TODO calculate hebap
            ff_ac3_bit_alloc_calc_bap(s->mask[ch], s->psd[ch], start, end,
                    snroffst, s->bit_alloc_params.floor,
                    s->bap[ch]);
        }


    }


    /* These fields for quantized mantissa values */

    got_cplchan = 0;
//    ff_ac3_get_transform_coeffs(gbc, s->bap, s->dexps, s->nfchans+s->lfeon, s->chincpl, s->dithflag, s->transform_coeffs, s->strtmant, s->endmant, &s->dith_state, s->ncplbnd, s->cplbndstrc, s->cplco);

    // TODO only for debug
    for(ch=0; ch<=s->nfchans+s->lfeon; ch++)
        memset(s->transform_coeffs[ch], 0, 256*sizeof(float));

    for(ch = 1; ch <= s->nfchans; ch++)
    {
        if(s->chahtinu[ch] == 0)
        {
            ff_ac3_get_transform_coeffs_ch(&m, gbc, s->dexps[ch], s->bap[ch], s->transform_coeffs[ch], s->strtmant[ch], s->endmant[ch], &s->dith_state);
            // ff_ac3_get_transform
            //memset(s->transform_coeffs[ch], 0, (256)*sizeof(float));
        }
        else if(s->chahtinu[ch] == 1)
        {
            int endbap, bin, n, m;
            int bg, g, bits, pre_chmant, remap;
            float mant;

            av_log(s->avctx, AV_LOG_INFO,  "AHT NOT TESTED");

            GET_BITS(s->chgaqmod[ch], gbc, 2);

            if (s->chgaqmod[ch] < 2){
                endbap = 12;
            }
            else{
                endbap = 17;
            }

            s->chactivegaqbins[ch] = 0;
            for(bin = 0; bin < s->endmant[ch]; bin++){
                if(s->hebap[ch][bin] > 7 && s->hebap[ch][bin] < endbap){
                    s->chgaqbin[ch][bin] = 1; /* Gain word is present */
                    s->chactivegaqbins[ch]++;
                }
                else if (s->hebap[ch][bin] >= endbap){
                    s->chgaqbin[ch][bin] = -1;/* Gain word not present */
                }else{
                    s->chgaqbin[ch][bin] = 0;
                }
            }


            switch(s->chgaqmod[ch]){
                    case EAC3_GAQ_NO: /* No GAQ gains present */
                        s->chgaqsections[ch] = 0;
                        break;
                    case EAC3_GAQ_12: /* GAQ gains 1 and 2 */
                    case EAC3_GAQ_14: /* GAQ gains 1 and 4 */
                        s->chgaqsections[ch] = s->chactivegaqbins[ch];
                        /* chactivegaqbins[ch] was computed earlier */
                        break;
                    case EAC3_GAQ_124: /* GAQ gains 1, 2, and 4 */
                        s->chgaqsections[ch] = s->chactivegaqbins[ch] / 3;
                        if (s->chactivegaqbins[ch] % 3) s->chgaqsections[ch]++;
                        break;
            }

            if((s->chgaqmod[ch] > 0x0) && (s->chgaqmod[ch] < 0x3) )
            {
                for(n = 0; n < s->chgaqsections[ch]; n++) { // TODO chgaqsections ?
                    GET_BITS(s->chgaqgain[ch][n], gbc, 1);
                }
            }
            else if(s->chgaqmod[ch] == 0x3)
            {
                int grpgain;
                for(n = 0; n < s->chgaqsections[ch]; n++) {
                    GET_BITS(grpgain, gbc, 5);
                    s->chgaqgain[ch][3*n]   = grpgain/9;
                    s->chgaqgain[ch][3*n+1] = (grpgain%9)/3;
                    s->chgaqgain[ch][3*n+2] = grpgain%3;
                }
            }

            // TODO test VQ and GAQ
            m=0;
            ///TODO calculate nchmant
            for(bin = 0; bin < s->nchmant[ch]; bin++)
            {
                if(s->chgaqbin[ch][bin]!=0)
                {
                    // GAQ (E3.3.4.2)
                    // XXX what about gaqmod = 0 ?
                    // difference between Gk=1 and gaqmod=0 ?
                    if(s->chgaqbin[ch][bin]>0){
                        // hebap in active range
                        // Gk = 1<<bg
                        bg = ff_gaq_gk[s->chgaqmod[ch]][s->chgaqgain[ch][m++]];
                    }else{
                        bg = 0;
                    }
                    bits = ff_bits_vs_hebap[s->hebap[ch][bin]];

                    for(n = 0; n < 6; n++){
                        // pre_chmant[n][ch][bin]
                        GET_SBITS(pre_chmant, gbc, bits-bg);
                        if(s->chgaqbin[ch][bin]>0 && bg && pre_chmant == -(1<<(bits-bg-1))){
                            // large mantissa
                            GET_SBITS(pre_chmant, gbc, bits - ((bg==1)?1:0));
                            mant = (float) pre_chmant / (1<<(bits - ((bg==1)?2:1)));
                            g = 0;
                            remap = 1;
                        }else{
                            // small mantissa
                            mant = (float) pre_chmant / (1<<(bits-bg-1));
                            g = bg;
                            remap = bg?0:1;
                        }

                        //TODO when remap needed ?
                        if(remap){
                            if(mant>=0){
                                mant = (float)
                                    ((ff_eac3_gaq_remap[s->hebap[ch][bin]-8][0][g][0] + 1.0f)
                                     * mant / (1<<g) +
                                     ff_eac3_gaq_remap[s->hebap[ch][bin]-8][0][g][1]) / 32768.0f;
                            }else{
                                mant = (float)
                                    ((ff_eac3_gaq_remap[s->hebap[ch][bin]-8][1][g][0] + 1.0f)
                                     * mant / (1<<g) +
                                     ff_eac3_gaq_remap[s->hebap[ch][bin]-8][1][g][1]) / 32768.0f;
                            }
                        }
                        s->pre_chmant[n][ch][bin] = mant;
                    }
                }
                else {
                    // hebap = 0 or VQ
                    if(s->hebap[ch][bin]){
                        GET_BITS(pre_chmant, gbc, ff_bits_vs_hebap[s->hebap[ch][bin]]);
                        for(n = 0; n < 6; n++){
                            s->pre_chmant[n][ch][bin] = ff_vq_hebap[s->hebap[ch][bin]][pre_chmant][n];
                        }
                    }else{
                        for(n = 0; n < 6; n++){
                            s->pre_chmant[n][ch][bin] = 0;
                        }
                    }
                }
            }
            s->chahtinu[ch] = -1; /* AHT info for this frame has been read ? do not read again */
        }
        if(s->chahtinu[ch] != 0){
            // TODO fast DCT
            int bin;
            float tmp;
            for(bin=0; bin<s->nchmant[ch]; bin++){
                tmp = 0;
                for(i=0; i<6; i++){
                    tmp += (i?sqrt(2):1) * s->pre_chmant[i][ch][bin] * cos(M_PI*i*(2*blk + 1)/12);
                }
                s->transform_coeffs[ch][bin] = tmp;
            }
        }
        if(s->cplinu[blk] && s->chincpl[ch] && !got_cplchan)
        {

            if(s->cplahtinu == 0)
            {
                /*
                for(bin = 0; bin < s->ncplmant; bin++) { // TODO ncplmant ?
                    GET_BITS(s->cplmant[bin], gbc, (0-16)); // TODO 0-16 :]
                }
                */
                ff_ac3_get_transform_coeffs_ch(&m, gbc, s->dexps[CPL_CH], s->bap[CPL_CH], s->transform_coeffs[CPL_CH], s->strtmant[CPL_CH], s->endmant[CPL_CH], &s->dith_state);
               // uncouple_channels(ctx);
                {
                    //TODO (form ac3)
                    int i, j, ch, bnd, subbnd;

                    subbnd = s->cplbegf+1;
                    i = s->strtmant[CPL_CH];
                    av_log(NULL, AV_LOG_DEBUG, "strtmant=%i endmant=%i\n", s->strtmant[CPL_CH], s->endmant[CPL_CH]);
                    av_log(NULL, AV_LOG_DEBUG, "ncplbnd=%i ncplsubbnd=%i\n", s->ncplbnd, s->ncplsubnd);
                    /*
                    for(bnd=0; bnd<256; bnd++){
                        av_log(NULL, AV_LOG_INFO, "%i: %f\n", bnd, transform_coeffs[CPL_CH][bnd]);
                    }*/
                    for(bnd=0; bnd<s->ncplbnd; bnd++) {
                        do {
                            for(j=0; j<12; j++) {
                                for(ch=1; ch<=s->nfchans; ch++) {// TODO lfe?
                                    if(s->chincpl[ch]) {
                                        s->transform_coeffs[ch][i] =
                                            s->transform_coeffs[CPL_CH][i] *
                                            s->cplco[ch][bnd] * 8.0f;
                                    }
                                }
                                av_log(NULL, AV_LOG_DEBUG, "%i ", i);
                                i++;
                            }
                            av_log(NULL, AV_LOG_DEBUG, "cplbndstrc[%i] = %i bnd=%i\n ", subbnd,
                                    s->cplbndstrc[subbnd], bnd);
                        } while(s->cplbndstrc[subbnd++] && subbnd<=s->cplendf);
                    }
                    av_log(NULL, AV_LOG_DEBUG, "\n");
                }
                got_cplchan = 1;
            }
            else if(s->cplahtinu == 1)
            {
                av_log(s->avctx, AV_LOG_ERROR,  "NOT IMPLEMENTED CPLINU && AHT");
                return -1;
#if 0
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
#endif
            }
            else {
                got_cplchan = 1;
            }
        }
    }

    if(s->lfeon) /* mantissas of low frequency effects channel */
    {
        if(s->lfeahtinu == 0)
        {
            //ff_ac3_get_transform
            ff_ac3_get_transform_coeffs_ch(&m, gbc, s->dexps[s->lfe_channel], s->bap[s->lfe_channel], s->transform_coeffs[s->lfe_channel], s->strtmant[s->lfe_channel], s->endmant[s->lfe_channel], &s->dith_state);
        }
        else if(s->lfeahtinu == 1)
        {
            av_log(s->avctx, AV_LOG_ERROR,  "NOT IMPLEMENTED");
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

int ff_eac3_parse_auxdata(GetBitContext *gbc, EAC3Context *s){
    // TODO
    return 0;
}
