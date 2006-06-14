/* AC3 Audio Decoder.
 *
 * Copyright (c) 2006 Kartikey Mahendra BHATT (bhattkm at gmail dot com)
 *
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
 */

#ifndef _AC3_DECODER_H
#define _AC3_DECODER_H

#include "avcodec.h"
#include "bitstream.h"
#include "dsputil.h"

/* Synchronization information. */
typedef struct {
    uint16_t sync_word;    //synchronization word = always 0x0b77
    uint16_t crc1;         //crc for the first 5/8 of the frame
    uint8_t  fscod;        //sampling rate code
    uint8_t  frmsizecod;   //frame size code

    /* Derived Attributes */
    int sampling_rate;    //sampling rate - 48, 44.1 or 32 kHz (value in Hz)
    int bit_rate;         //nominal bit rate (value in kbps)
    int frame_size;       //16 bit words before next sync word
} ac3_sync_info;



/* flags for the BSI. */
#define AC3_BSI_LFEON      0x00000001      //low frequency effects channel on
#define AC3_BSI_COMPRE     0x00000002      //compression exists
#define AC3_BSI_LANGCODE   0x00000004      //langcode exists
#define AC3_BSI_AUDPRODIE  0x00000008      //audio production information exists
#define AC3_BSI_COMPR2E    0x00000010      //compr2 exists
#define AC3_BSI_LANGCOD2E  0x00000020      //langcod2 exists
#define AC3_BSI_AUDPRODI2E 0x00000040      //audio production information 2 exists
#define AC3_BSI_COPYRIGHTB 0x00000080      //copyright
#define AC3_BSI_ORIGBS     0x00000100      //original bit stream
#define AC3_BSI_TIMECOD1E  0x00000200      //timecod1 exists
#define AC3_BSI_TIMECOD2E  0x00000400      //timecod2 exists
#define AC3_BSI_ADDBSIE    0x00000800      //additional bit stream information exists

/* Bit Stream Information. */
typedef struct {
    uint32_t flags;
    uint8_t  bsid;          //bit stream identification
    uint8_t  bsmod;         //bit stream mode - type of service
    uint8_t  acmod;         //audio coding mode - which channels are in use
    uint8_t  cmixlev;       //center mix level
    uint8_t  surmixlev;     //surround mix level
    uint8_t  dsurmod;       //dynamic surround encoded
    uint8_t  dialnorm;      //dialog normalization
    uint8_t  compr;         //compression gain word
    uint8_t  langcod;       //language code
    uint8_t  mixlevel;      //mixing level
    uint8_t  roomtyp;       //room type
    uint8_t  dialnorm2;     //dialogue normalization for 1+1 mode
    uint8_t  compr2;        //compression gain word for 1+1 mode
    uint8_t  langcod2;      //language code for 1+1 mode
    uint8_t  mixlevel2;     //mixing level for 1+1 mode
    uint8_t  roomtyp2;      //room type for 1+1 mode
    uint16_t timecod1;      //timecode 1
    uint16_t timecod2;      //timecode 2
    uint8_t  addbsil;       //additional bit stream information length

    /* Dervied Attributes */
    int      nfchans;      //number of full bandwidth channels - derived from acmod
} ac3_bsi;



/* #defs relevant to Audio Block. */
#define MAX_FBW_CHANNELS  5      //maximum full bandwidth channels
#define NUM_LFE_GROUPS    3      //number of LFE Groups
#define MAX_NUM_SEGS      8      //maximum number of segments per delta bit allocation
#define NUM_LFE_MANTS     7      //number of lfe mantissas
#define MAX_CPL_SUBNDS    18     //maximum number of coupling sub bands
#define MAX_CPL_BNDS      18     //maximum number of coupling bands
#define MAX_CPL_GRPS      253    //maximum number of coupling groups
#define MAX_CHNL_GRPS     88     //maximum number of channel groups
#define MAX_NUM_MANTISSAS 256    //maximum number of mantissas

/* flags for the Audio Block. */
#define AC3_AB_DYNRNGE   0x00000001    //dynamic range control exists
#define AC3_AB_DYNRNG2E  0x00000002    //dynamic range control 2 exists
#define AC3_AB_CPLSTRE   0x00000004    //coupling strategy exists
#define AC3_AB_CPLINU    0x00000008    //coupling in use
#define AC3_AB_PHSFLGINU 0x00000010    //phase flag in use
#define AC3_AB_REMATSTR  0x00000020    //rematrixing required
#define AC3_AB_LFEEXPSTR 0x00000100    //lfe exponent strategy
#define AC3_AB_BAIE      0x00000200    //bit allocation information exists
#define AC3_AB_SNROFFSTE 0x00000400    //SNR offset exists
#define AC3_AB_CPLLEAKE  0x00000800    //coupling leak initialization exists
#define AC3_AB_DELTBAIE  0x00001000    //delta bit allocation information exists
#define AC3_AB_SKIPLE    0x00002000    //skip length exists

/* Exponent strategies. */
#define AC3_EXPSTR_D15   0x01
#define AC3_EXPSTR_D25   0x02
#define AC3_EXPSTR_D45   0x03
#define AC3_EXPSTR_REUSE 0x00

/* Audio Block */
typedef struct {
    uint32_t flags;
    uint8_t  blksw;              //block switch flags for channels in use
    uint8_t  dithflag;           //dithering flags for channels in use
    uint8_t  dynrng;             //dynamic range word
    uint8_t  dynrng2;            //dynamic range word for 1+1 mode
    uint8_t  chincpl;            //channel in coupling flags for channels in use
    uint8_t  cplbegf;            //coupling begin frequency code
    uint8_t  cplendf;            //coupling end frequency code
    uint32_t cplbndstrc;         //coupling band structure
    uint8_t  cplcoe;             //coupling co-ordinates exists for the channel in use
    uint8_t  mstrcplco[5];       //master coupling co-ordinate for channels in use
    uint8_t  cplcoexp[5][18];    //coupling co-ordinate exponenets
    uint8_t  cplcomant[5][18];   //coupling co-ordinate mantissas
    uint32_t phsflg;             //phase flag per band
    uint8_t  rematflg;           //rematrixing flag
    uint8_t  cplexpstr;          //coupling exponent strategy
    uint8_t  chexpstr[5];        //channel exponent strategy
    uint8_t  lfeexpstr;          //lfe exponent strategy
    uint8_t  chbwcod[5];         //channel bandwdith code for channels in use
    uint8_t  cplabsexp;          //coupling absolute exponent
    uint8_t  cplexps[72];        //coupling exponents
    uint8_t  exps[5][88];        //channel exponents
    uint8_t  gainrng[5];         //gain range
    uint8_t  lfeexps[3];         //LFE exponents
    uint8_t  sdcycod;            //slow decay code
    uint8_t  fdcycod;            //fast decay code
    uint8_t  sgaincod;           //slow gain code
    uint8_t  dbpbcod;            //dB per bit code
    uint8_t  floorcod;           //masking floor code
    uint8_t  csnroffst;          //coarse SNR offset
    uint8_t  cplfsnroffst;       //coupling fine SNR offset
    uint8_t  cplfgaincod;        //coupling fast gain code
    uint8_t  fsnroffst[5];       //fine SNR offset for channels in use
    uint8_t  fgaincod[5];        //fast gain code for channels in use
    uint8_t  lfefsnroffst;       //lfe fine SNR offset
    uint8_t  lfefgaincod;        //lfe fast gain code
    uint8_t  cplfleak;           //coupling fast leak initialization value
    uint8_t  cplsleak;           //coupling slow leak initialization value
    uint8_t  cpldeltbae;         //coupling delta bit allocation exists
    uint8_t  deltbae[5];         //delta bit allocation exists for channels in use
    uint8_t  cpldeltnseg;        //coupling delta bit allocation number of segments
    uint8_t  cpldeltoffst[8];    //coupling delta offset
    uint8_t  cpldeltlen[8];      //coupling delta len
    uint8_t  cpldeltba[8];       //coupling delta bit allocation
    uint8_t  deltnseg[5];        //delta bit allocation number of segments per channel
    uint8_t  deltoffst[5][8];    //delta offset for channels in use
    uint8_t  deltlen[5][8];      //delta len for channels in use
    uint8_t  deltba[5][8];       //delta bit allocation
    uint16_t skipl;              //skip length

    /* Derived Attributes */
    int      ncplsubnd;          //number of active coupling sub bands = 3 + cplendf - cplbegf
    int      ncplbnd;            //derived from ncplsubnd and cplbndstrc
    int      ncplgrps;           //derived from ncplsubnd, cplexpstr
    int      nchgrps[5];         //derived from chexpstr, and cplbegf or chbwcod
    int      nchmant[5];         //derived from cplbegf or chbwcod
    int      ncplmant;           //derived from ncplsubnd = 12 * ncplsubnd

    uint8_t  cplstrtbnd;         //coupling start band for bit allocation
    uint8_t  cplstrtmant;        //coupling start mantissa
    uint8_t  cplendmant;         //coupling end mantissa
    uint8_t  endmant[5];         //channel end mantissas

    uint8_t  dcplexps[256];      //decoded coupling exponents
    uint8_t  dexps[5][256];      //decoded fbw channel exponents
    uint8_t  dlfeexps[256];      //decoded lfe exponents
    uint8_t  cplbap[256];        //coupling bit allocation parameters table
    uint8_t  bap[5][256];        //fbw channels bit allocation parameters table
    uint8_t  lfebap[256];        //lfe bit allocaiton parameters table

    float    cplcoeffs[256];     //temporary storage for coupling transform coefficients
    float    cplco[5][18];       //coupling co-ordinates
    float    *ab_samples;        //points to the samples for "this" audio block
} ac3_audio_block;



/* AC3 Context. */
typedef struct {
    ac3_sync_info   sync_info;
    ac3_bsi         bsi;
    ac3_audio_block audio_block;
    float           *samples;
    MDCTContext     mdct_ctx_256;
    MDCTContext     mdct_ctx_512;
    GetBitContext   gb;
} AC3DecodeContext;

#endif
