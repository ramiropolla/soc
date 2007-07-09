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
#ifndef EAC3_H
#define EAC3_H

// TODO clean...

#include "dsputil.h"
#include "avcodec.h"
#include "ac3.h"
#include "random.h"

#define DBA_NEW      0x01
#define DBA_NONE     0x02
#define DBA_RESERVED 0x03
#define DBA_REUSE    0x00

// TODO
#define MAX_CHANNELS AC3_MAX_CHANNELS

/** adjustments in dB gain */
#define LEVEL_MINUS_3DB         (0.7071067811865476f) /* sqrt(2)/2 */
#define LEVEL_MINUS_4POINT5DB   (0.5946035575013605f)
#define LEVEL_MINUS_6DB         (0.5000000000000000f)
#define LEVEL_PLUS_3DB          (1.4142135623730951f) /* sqrt(2) */
#define LEVEL_PLUS_6DB          (2.0000000000000000f)
#define LEVEL_ZERO              (0.0000000000000000f)

#define AC3_BLOCK_SIZE  256
#define MAX_BLOCKS 6
#define MAX_SPX_CODES 18
#define TODO_SIZE 1024


/** output configurations. */
#define AC3_OUTPUT_UNMODIFIED 0x01
#define AC3_OUTPUT_MONO       0x02
#define AC3_OUTPUT_STEREO     0x04
#define AC3_OUTPUT_DOLBY      0x08
#define AC3_OUTPUT_LFEON      0x10
/** AC-3 channel mode (audio coding mode) */

typedef enum {
    AC3_CHANNEL_MODE_DUALMONO = 0,
    AC3_CHANNEL_MODE_MONO,
    AC3_CHANNEL_MODE_STEREO,
    AC3_CHANNEL_MODE_3F,
    AC3_CHANNEL_MODE_2F_1R,
    AC3_CHANNEL_MODE_3F_1R,
    AC3_CHANNEL_MODE_2F_2R,
    AC3_CHANNEL_MODE_3F_2R
} AC3ChannelMode;

typedef struct EAC3Context{
    // TODO erase unused variables
    // TODO add documentation
//Syncinfo
    int syncword;

//BSI
    int strmtyp; // 2);              ///< Stream type
    int substreamid; // 3);          ///< Substream identification
    int frmsiz; // 11);              ///< Frame size
    int fscod; // 2);                ///< Sample rate code
    int fscod2; // 2);               ///< Sample rate code 2
    int numblkscod; // 2);           ///< Number of audio blocks
    int acmod; // 3);                ///< Audio coding mode
    int lfeon; // 1);                ///< Low frequency effect channel on
    int bsid; // 5);                 ///< Bit stream identification
    int dialnorm; // 5);             ///<Dialogue normalization
    int compre; // 1);               ///< Compression gain word exists
    int compr; // 8);                ///< Compression gain word
    int dialnorm2; // 5);            ///< Dialog normalization ch2
    int compr2e; // 1);              ///< Compression gain word ch2 exists
    int compr2; // 8);               ///< Compression gain word ch2
    int chanmape; // 1);             ///< Custom channel map exists
    int chanmap; // 16);             ///< Custom channel map
    int mixmdate; // 1);             ///< Mixing meta-data exists
    int dmixmod; // 2);              ///<
    int ltrtcmixlev; // 3);
    int lorocmixlev; // 3);
    int ltrtsurmixlev; // 3);
    int lorosurmixlev; // 3);
    int lfemixlevcode; // 1);        ///< lfe mix level code exists
    int lfemixlevcod; // 5);         ///< lfe mix level code
    int pgmscle; // 1);              ///<  Program scale factor exists
    int pgmscl; // 6);               ///< Program scale factor
    int pgmscl2e; // 1);             ///< Program scale factor #2 exists
    int pgmscl2; // 6);              ///< Program scale factor #2
    int extpgmscle; // 1);           ///< External program scale factor exists
    int extpgmscl; // 6);            ///< External program scale factor
    int mixdef; // 2);               ///< Mix control type
    int mixdata; // 12);             ///< Mixing parameter data
    int mixdeflen; // 5);            ///< Length of mixing parameter data field
//    int *mixdata; // 8*(mixdeflen+2));
    int paninfoe; // 1);             ///< Pan information exists
    int paninfo; // 14);             ///< Pan information
    int paninfo2e; // 1);            ///< Pan information 2 exists
    int paninfo2; // 14);            ///< Pan information 2
    int frmmixcfginfoe; // 1);       ///< Frame mixing configuration information exists
    int blkmixcfginfo0; // 5);
    int blkmixcfginfoe; // 1);       ///< Block mixing configuration information exists
    int blkmixcfginfoblk; // 5);     ///< Block mixing configuration information
    int infomdate; // 1);            ///< Informational meta-data exists
    int bsmod; // 3);                ///< Bit stream mode
    int copyrightb; // 1);           ///< Copyright bit
    int origbs; // 1); ///< Original bit stream
    int dsurmod; // 2);
    int dheadphonmod; // 2);
    int dsurexmod; // 2);
    int audprodie; // 1);
    int mixlevel; // 5);             ///< Mix level
    int roomtyp; // 2);              ///< Room type
    int adconvtyp; // 1);            ///< A/D converter type
    int audprodi2e; // 1);           ///< Audio production information exists ch2
    int mixlevel2; // 5);            ///< Mixing level ch2
    int roomtyp2; // 2);             ///< room type ch2
    int adconvtyp2; // 1);           ///< A/D converter type
    int sourcefscod; // 1);          ///< Source sample rate code
    int convsync; // 1);             ///< Converter synchronization flag
    int blkid; // 1);                ///< Block identification
    int frmsizecod; // 6);           ///< Frame size code
    int addbsie; // 1);              ///< Additional bit stream information exists
    int addbsil; // 6);              ///< Additional bit stream information length
    int addbsi[64];                  ///< Additional bit stream information

//Audfrm
    int expstre; // 1);              ///< Exponent strategy syntax enabled
    int ahte; // 1);                 ///< Adaptive hybrid transform enabled
    int snroffststr; // 2);          ///< SNR offset strategy
    int transproce; // 1);           ///< Transient pre-noise processing enabled
    int blkswe; // 1);               ///< Block switch syntax enabled
    int dithflage; // 1);            ///< Dither flag syntax enabled
    int bamode; // 1);               ///< Bit allocation model syntax enabled
    int frmfgaincode; // 1);         ///< Fast gain codes enabled
    int dbaflde; // 1);              ///< Delta bit allocation syntax enabled
    int skipflde; // 1);             ///< Skip Filed syntax enabled
    int spxattene; // 1);            ///< Spectral extension attenuation enabled
    int cplinu[MAX_BLOCKS]; // 1);   ///< Coupling in use
    int cplstre[MAX_BLOCKS]; // 1);  ///< Coupling strategy exists
    int cplexpstr[MAX_BLOCKS]; // 2);///< Coupling exponents strategy
    int chexpstr[MAX_BLOCKS][MAX_CHANNELS]; // 2); ///< Channel exponent strategy
    int frmcplexpstr; // 5);         ///< Frame based coupling exponent strategy
    int frmchexpstr[MAX_CHANNELS];   ///< frame based channel exponent strategy
    int lfeexpstr[MAX_BLOCKS];       ///< Lfe exponent strategy
    int convexpstre; // 1);          ///< Converter exponent strategy exists
    int convexpstr[MAX_CHANNELS];    ///< Converter channel exponent strategy
    int cplahtinu; // 1);            ///< Coupling channel AHT in use
    int chahtinu[MAX_CHANNELS];      ///< Channel AHT in use
    int lfeahtinu; // 1);            ///< Lfe channel AHT in use
    int frmcsnroffst; // 6);         ///< Frame coarse SNR offset
    int frmfsnroffst; // 4);         ///< Frame fine SNR offset
    int chintransproc[MAX_CHANNELS]; ///< Channel in transient pre-noise processing
    int transprocloc[MAX_CHANNELS];  ///< Transient location relative to start of frame
    int transproclen[MAX_CHANNELS];  ///< Transient processing length
    int chinspxatten[MAX_CHANNELS];  ///< Channel in spectral extension attenuation process
    int spxattencod[MAX_CHANNELS];   ///< spectral extension attenuation code
    int blkstrtinfoe; // 1);         ///< Block start information exists
    uint32_t blkstrtinfo;            ///< Block start information
    int ncplblks;

// EAC3Audblk
    int blksw[MAX_CHANNELS]; // 1);  ///< Block switch flag
    int dithflag[MAX_CHANNELS];      ///< Dither flag
    int dynrnge; // 1);              ///< Dynamic range gain word exists
    int dynrng; // 8);               ///< Dynamic range gain word
    int dynrng2e; // 1);             ///< Dynamic range gain word exists, ch2
    int dynrng2; // 8);              ///< Dynamic range gain word
    int spxstre; // 1);              ///< Spectral extension strategy exists
    int spxinu; // 1);               ///< spectral extension in use
    int chinspx[MAX_CHANNELS];       ///< Channel in spectral extension
    int spxstrtf; // 2);             ///< Spectral extension start copy frequency code
    int spxbegf; // 3);              ///< Spectral extension begin frequency code
    int spxendf; // 3);              ///< Spectral extension end frequency code
    int spxbndstrce; // 1);          ///< Spectral extension band structure exists
    int spxbndstrc[MAX_SPX_CODES];   ///< Spectral extension band structure
    int spxcoe[MAX_CHANNELS];        ///< Spectral extension coordinates exists
    int spxblnd[MAX_CHANNELS];       ///< Spectral extension blend
    int mstrspxco[MAX_CHANNELS];     ///< Master spectral extension coordinates
    int spxcoexp[MAX_CHANNELS][MAX_SPX_CODES];  ///< Spectral extension coordinate exponent
    int spxcomant[MAX_CHANNELS][MAX_SPX_CODES]; ///< Spectral extension coordinate mantissa
    int ecplinu; // 1);              ///< Enhanced coupling in use
    int chincpl[MAX_CHANNELS];       ///< Channel in coupling
    int phsflginu; // 1);            ///< Phase flag in use
    int cplbegf; // 4);              ///< Coupling begin frequency code
    int cplendf; // 4);              ///< Coupling end frequency code
    int cplbndstrce; // 1);          ///< Coupling band structure exists
    int cplbndstrc[19]; // 1);       ///< Coupling band structure
    int ecplbegf; // 4);             ///< Enhanced coupling begin frequency code
    int ecplendf; // 4);             ///< Enhanced coupling end frequency code
    int ecplbndstrce; // 1);         ///< Enhanced coupling band structure exists
    int ecplbndstrc[TODO_SIZE];      ///< Enhanced coupling band structure
    int cplcoe[MAX_CHANNELS];        ///< Coupling coordinates exists
    int mstrcplco[MAX_CHANNELS];     ///< Master coupling coordinates
    int cplcoexp[MAX_CHANNELS][TODO_SIZE]; ///< Coupling coordinate exponent
    int cplcomant[MAX_CHANNELS][TODO_SIZE]; ///< Coupling coordinate mantissa
    int phsflg[TODO_SIZE]; // 1);    ///< Phase flag
    int ecplangleintrp; // 1);       ///< Enhanced coupling angle interpolation flag
    int ecplparam1e[MAX_CHANNELS];   ///< Enhanced coupling parameters 1 exists
    int ecplparam2e[MAX_CHANNELS];   ///< Enhanced coupling parameters 2 exists
    int ecplamp[MAX_CHANNELS][TODO_SIZE];   ///< Enhanced coupling amplitude scaling
    int ecplangle[MAX_CHANNELS][TODO_SIZE]; ///< Enhanced coupling angle
    int ecplchaos[MAX_CHANNELS][TODO_SIZE]; ///< Enhanced coupling chaos
    int ecpltrans[MAX_CHANNELS];     ///< Enhanced coupling transient present
    int rematstr; // 1);             ///< Rematrixing strategy
    int rematflg[TODO_SIZE]; // 1);  ///< Rematrixing flag
    int chbwcod[MAX_CHANNELS];       ///< Rematrixing strategy
    int cplabsexp; // 4);            ///< Coupling absolute exponent
    int cplexps[TODO_SIZE]; // 7);   ///< Coupling exponent

    int gainrng[MAX_CHANNELS];  ///< Channel Gain range code
//    int lfeexps[3]; // 7); // 0...nlfegrps = const 0...2
    int baie; // 1);                 ///< Bit allocation information exists
    int sdcycod; // 2);              ///< Slow decay code
    int fdcycod; // 2);              ///< Fast decay code
    int sgaincod; // 2);             ///< Slow gain code
    int dbpbcod; // 2);              ///< dB per bit code
    int floorcod; // 3);             ///< Masking floor code
    int snroffste; // 1);            ///< SNR offset exists
    int csnroffst; // 6);            ///< Coarse SNR offset
    int blkfsnroffst; // 4);         ///< Block fine SNR offset
    int cplfsnroffst; // 4);         ///< Coupling fine SNR offset
    int fsnroffst[MAX_CHANNELS];     ///< Channel fine SNR offset
    int lfefsnroffst; // 4);         ///< Lfe fine SNR offset
    int fgaincode; // 1);            ///< Channel fast gain code enabled
    int cplfgaincod; // 3);          ///< Coupling fast gain code code
    int fgaincod[MAX_CHANNELS];      ///< Channel fast gain code
    int lfefgaincod; // 3);          ///< Lfe fast gain code
    int convsnroffste; // 1);        ///< Converter SNR offset exists
    int convsnroffst; // 10);        ///< Converter SNR offset
    int cplleake; // 1);             ///< Coupling leak initialization exists
    int cplfleak; // 3);             ///< Coupling fast leak initialization
    int cplsleak; // 3);             ///< Coupling slow leak initialization
    int deltbaie; // 1);             ///< Delta bit allocation information exists
    int cpldeltbae; // 2);           ///< Coupling delta bit allocation exists
    uint8_t deltbae[MAX_CHANNELS];   ///< Delta bit allocation exists
    int cpldeltnseg; // 3);          ///< Coupling delta bit allocation number of segments
    int cpldeltoffst[9]; // 5);      ///< Coupling bit allocation offset
    int cpldeltlen[9]; // 4);        ///< Coupling delta bit allocation length
    int cpldeltba[9]; // 3);         ///< Coupling delta bit allocation
    uint8_t deltnseg[MAX_CHANNELS];  ///< Channel delta bit allocation number of segments
    uint8_t deltoffst[MAX_CHANNELS][9]; ///< Channel delta bit allocation offset
    uint8_t deltlen[MAX_CHANNELS][9];   ///< Channel delta bit allocation length
    uint8_t deltba[MAX_CHANNELS][9];    ///< Channel delta bit allocation
    int skiple; // 1);               ///< Skip length exists
    int skipl; // 9);                ///< Skip length
    int skipfld; // TODO skipl * 8); ///< Skip field

    int got_cplchan;
//    int chmant[MAX_CHANNELS][TODO_SIZE];
    int chgaqmod[MAX_CHANNELS];                 ///< Channel gain adaptive quantization mode
    int chgaqgain[MAX_CHANNELS][TODO_SIZE];     ///< Channel gain adaptive quantization gain
    int pre_chmant[6][MAX_CHANNELS][TODO_SIZE]; ///< Pre channel mantissas
    int cplmant[TODO_SIZE]; // cplmant          ///< Coupling mantisass
    int cplgaqmod;                              ///< Coupling channel gain adaptive quantization mode
    int cplgaqgain[TODO_SIZE];                  ///< Coupling gain adaptive quantization gain
    int pre_cplmant[6][TODO_SIZE];              ///< Pre coupling channel_mantissas
    int lfegaqmod;                              ///<  Lfe channel gain adaptive quantization mode
    int lfegaqgain[TODO_SIZE];                  ///<  LFE channel gain adaptive quantization mode
    int pre_lfemant[6][TODO_SIZE];              ///< Pre lfe channel mantissas




    int firstspxcos[MAX_CHANNELS]; // TODO type ? ///< First spectral extension coordinates states
    int firstcplcos[MAX_CHANNELS]; // TODO type ? ///< First coupling coordinates states
    int firstcplleak; // TODO type ?              ///< First coupling leak state



    // TODO
    int lfegaqbin[TODO_SIZE];
    int lfegaqsections;
    int cplgaqbin[TODO_SIZE];
    int cplgaqsections;
    int chgaqbin[MAX_CHANNELS][TODO_SIZE]; // [][nchmant]
    int chgaqsections[MAX_CHANNELS];
    int nchmant[MAX_CHANNELS];         ///< Number of fbw channel mantissas
    int ncplsubnd;                     ///< Number of coupling sub-bands
    int ncplbnd;                       ///< Number of structured coupled bands
    int nrematbnds;
    int ncplgrps;                      ///< Number of coupled exponent groups
    int ncplmant;                      ///< Number of coupled mantissas

    int nlfemant;                      ///< Number of lfe mantissas
    int nchgrps[MAX_CHANNELS];         ///< Number of fbw channel exponent groups
    uint8_t dexps[MAX_CHANNELS][AC3_MAX_COEFS]; ///< Differential exponents

    int cplstrtmant;
    int cplendmant;
    int endmant[MAX_CHANNELS];
    int strtmant[MAX_CHANNELS];
    int firstchincpl;
    int ecpl_start_subbnd, ecpl_end_subbnd;
    int nlfegrps; ///< Number of lfe channel exponent groups

    int necplbnd;
    int nspxbnds;
    int spxbndsztab[MAX_SPX_CODES]; // max_spxendf + 1 = 18 ?
    int nfchans; ///< Number of fbw channels

    uint8_t bap[MAX_CHANNELS][AC3_MAX_COEFS]; ///< bit allocation pointers
    int16_t psd[MAX_CHANNELS][AC3_MAX_COEFS]; ///< scaled exponents
    int16_t bndpsd[MAX_CHANNELS][350];        ///< interpolated exponents FIXME in ac3dec [50] !?
    int16_t mask[MAX_CHANNELS][350];          ///< masking values


    DECLARE_ALIGNED_16(float, transform_coeffs[MAX_CHANNELS][AC3_MAX_COEFS]);

    /// TODO move to AC3BitAllocParameters
    int sdecay;
    int fdecay;
    int sgain;
    int dbknee;
    int floor;

    AC3BitAllocParameters bit_alloc_params;

    AVRandomState dith_state;        ///< for dither generation


    int ntchans; ///< Total of all channels
    int lfe_channel;

    GetBitContext *gbc;

    MDCTContext imdct_512;           ///< for 512 sample imdct transform
    MDCTContext imdct_256;           ///< for 256 sample imdct transform
    DSPContext  dsp;                 ///< for optimization

    DECLARE_ALIGNED_16(float, delay[MAX_CHANNELS][AC3_BLOCK_SIZE]);  ///< delay - added to the next block
    DECLARE_ALIGNED_16(float, window[AC3_BLOCK_SIZE]);               ///< window coefficients
    DECLARE_ALIGNED_16(float, tmp_output[AC3_BLOCK_SIZE * 24]);      ///< temp storage for output before windowing
    DECLARE_ALIGNED_16(float, tmp_imdct[AC3_BLOCK_SIZE * 24]);       ///< temp storage for imdct transform
    DECLARE_ALIGNED_16(float, output[MAX_CHANNELS][AC3_BLOCK_SIZE]); ///< output after imdct transform and windowing
    DECLARE_ALIGNED_16(int16_t, int_output[MAX_CHANNELS][AC3_BLOCK_SIZE]);///< final 16-bit integer output


    float add_bias;                  ///< offset for float_to_int16 conversion
    float mul_bias;                  ///< scaling for float_to_int16 conversion

    AC3ChannelMode  blkoutput;
}EAC3Context;


int eac3_parse_syncinfo(GetBitContext *gbc, EAC3Context *s);
int eac3_parse_bsi(GetBitContext *gbc, EAC3Context *s);
int eac3_parse_audfrm(GetBitContext *gbc, EAC3Context *s);
int eac3_parse_audblk(GetBitContext *gbc, EAC3Context *s, const int blk);
int eac3_parse_auxdata(GetBitContext *gbc, EAC3Context *s);

#endif
