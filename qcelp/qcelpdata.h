/*
 * QCELP decoder
 * Copyright (c) 2007 Reynaldo H. Verdejo Pinochet
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
 * @file qcelpdata.h
 * QCELP decoder
 */

typedef enum
{
    RATE_FULL   = 0,
    RATE_HALF   = 1,
    RATE_QUARTER= 2,
    RATE_OCTAVE = 3,
    I_F_Q,          /*!< insufficient frame quality */
    BLANK,
    RATE_UNKNOWN
} qcelp_packet_rate;

static const uint16_t qcelp_bits_per_rate[]={266,124,54,20};

typedef struct {
    uint8_t index;  /*!< index into the reference frame */
    uint8_t bitpos; /*!< bit position in the value's byte */
} QCELPBitmap;


/**
 * WARNING
 *
 * YOU WONT SEE ANY mention of a REFERENCE nor an UNIVERSAL frame
 * in the specs, this is just some internal way of handling the
 * reordering needed to unify the decoding process _inside_ this
 * code, nothing more.
 *
 *
 * UNIVERSAL FRAME
 * ---------------
 *
 * Format of QCELPFrame.data
 *
 *     QCELP_X0_POS
 *           |
 * CBSIGNs   0     1  2  3  4  5  6  7  8  9 10 11 12 13 14 15
 * CBGAINs  16    17 18 19 20 21 22 23 24 25 26 27 28 29 30 31
 * CINDEXs  32    33 34 35 36 37 38 39 40 41 43 43 44 45 46 47
 * PLAGs    48    49 50 51
 * PFRACs   52    53 54 55
 * PGAINs   56    57 58 59
 * LSPVs    60    61 62 63 64
 * RSVD     65
 * LSP      66    67 68 69 70 71 72 73 74 75
 * CBSEED   76
 *
 *
 * REFERENCE FRAME
 * ---------------
 *
 *
 * What follows are the reference frame slices. Each tuple will be mapped
 * to a QCELPBitmap showing the location of each bit in the input with respect
 * to a transmission code in the 'universal frame'.
 *
 * FIXME
 * it would be really nice if someone reviewed these numbers :)
 *---------------------------------------------------------------------------*/

#define QCELP_RATE_FULL_BITMAP \
{15,0},{47,0},{47,1},{47,2},{47,3},{47,4},{47,5},{47,6},\
{65,0},{65,1},{45,6},{30,0},{30,1},{30,2},{30,3},{14,0},\
{46,0},{46,1},{46,2},{46,3},{46,4},{46,5},{46,6},{31,0},\
{31,1},{31,2},{44,2},{44,3},{44,4},{44,5},{44,6},{29,0},\
{29,1},{29,2},{29,3},{13,0},{45,0},{45,1},{45,2},{45,3},\
{45,4},{45,5},{59,2},{51,0},{51,1},{51,2},{51,3},{51,4},\
{51,5},{51,6},{55,0},{28,0},{28,1},{28,2},{28,3},{12,0},\
{44,0},{44,1},{42,4},{42,5},{42,6},{27,0},{27,1},{27,2},\
{11,0},{43,0},{43,1},{43,2},{43,3},{43,4},{43,5},{43,6},\
{59,0},{59,1},{41,0},{41,1},{41,2},{41,3},{41,4},{41,5},\
{41,6},{26,0},{26,1},{26,2},{26,3},{10,0},{42,0},{42,1},\
{42,2},{42,3},{24,1},{24,2},{24,3},{ 8,0},{40,0},{40,1},\
{40,2},{40,3},{40,4},{40,5},{40,6},{24,0},{24,1},{24,2},\
{24,3},{ 9,0},{39,3},{39,4},{39,5},{39,6},{58,0},{58,1},\
{58,3},{50,0},{50,1},{50,2},{50,3},{50,4},{50,5},{50,6},\
{54,0},{24,0},{22,3},{ 6,0},{38,0},{38,1},{38,2},{38,3},\
{38,4},{38,5},{38,6},{23,0},{23,1},{23,2},{ 7,0},{39,0},\
{39,1},{39,2},{36,6},{21,0},{21,1},{21,2},{21,3},{ 5,0},\
{37,0},{37,1},{37,2},{37,3},{37,4},{37,5},{37,6},{22,0},\
{22,1},{22,2},{49,3},{49,4},{49,5},{49,6},{53,0},{20,0},\
{20,1},{20,2},{20,3},{ 4,0},{36,0},{36,1},{36,2},{36,3},\
{36,4},{36,5},{19,1},{19,2},{ 3,0},{35,0},{35,1},{35,2},\
{35,3},{35,4},{35,5},{35,6},{57,0},{57,1},{57,2},{49,0},\
{49,1},{49,2},{33,4},{33,5},{33,6},{18,0},{18,1},{18,2},\
{18,3},{ 2,0},{34,0},{34,1},{34,2},{34,3},{34,4},{34,5},\
{34,6},{19,0},{32,0},{32,1},{32,2},{32,3},{32,4},{32,5},\
{32,6},{17,0},{17,1},{17,2},{17,3},{ 1,0},{33,0},{33,1},\
{33,2},{33,3},{56,0},{56,1},{56,2},{48,0},{48,1},{48,2},\
{48,3},{48,4},{48,5},{48,6},{52,0},{16,0},{16,1},{16,2},\
{ 0,3},{ 0,0},{62,3},{62,4},{62,5},{62,6},{63,0},{63,1},\
{63,2},{63,3},{63,4},{63,5},{64,0},{64,1},{64,2},{64,3},\
{64,4},{64,5},{60,0},{60,1},{60,2},{60,3},{60,4},{60,5},\
{61,0},{61,1},{61,2},{61,3},{61,4},{61,5},{61,6},{62,0},\
{62,1},{62,2}

#define QCELP_RATE_HALF_BITMAP \
{19,0},{19,1},{19,2},{19,3},{ 3,0},{35,0},{35,1},{35,2},\
{35,3},{35,4},{35,5},{35,6},{34,2},{34,3},{34,4},{34,5},\
{34,6},{59,0},{59,1},{59,2},{51,0},{51,1},{51,2},{51,3},\
{51,4},{51,5},{51,6},{55,0},{58,2},{50,0},{50,1},{50,2},\
{50,3},{50,4},{50,5},{50,6},{54,0},{18,0},{18,1},{18,2},\
{18,3},{ 2,0},{34,0},{34,1},{49,6},{53,0},{17,0},{17,1},\
{17,2},{17,3},{ 1,0},{33,0},{33,1},{33,2},{33,3},{33,4},\
{33,5},{33,6},{58,0},{58,1},{32,0},{32,1},{32,2},{32,3},\
{32,4},{32,5},{32,6},{57,0},{57,1},{57,2},{49,0},{49,1},\
{49,2},{49,3},{49,4},{49,5},{56,0},{56,1},{56,2},{48,0},\
{48,1},{48,2},{48,3},{48,4},{48,5},{48,6},{52,1},{16,0},\
{16,1},{16,2},{16,3},{ 0,0},{62,3},{62,4},{62,5},{62,6},\
{63,0},{63,1},{63,2},{63,3},{63,4},{63,5},{64,0},{64,1},\
{64,2},{64,3},{64,4},{64,5},{60,0},{60,1},{60,2},{60,3},\
{60,4},{60,5},{61,0},{61,1},{61,2},{61,3},{61,4},{61,5},\
{61,6},{62,0},{62,1},{62,2}

#define QCELP_RATE_4THR_BITMAP \
{20,0},{20,1},{20,2},{20,3},{65,0},{65,1},{16,0},{16,1},\
{16,2},{16,3},{17,0},{17,1},{17,2},{17,3},{18,0},{18,1},\
{18,2},{18,3},{19,0},{19,1},{19,2},{19,3},{62,3},{62,4},\
{62,5},{62,6},{63,0},{63,1},{63,2},{63,3},{63,4},{63,5},\
{64,0},{64,1},{64,2},{64,3},{64,4},{64,5},{60,0},{60,1},\
{60,2},{60,3},{60,4},{60,5},{61,0},{61,1},{61,2},{61,3},\
{61,4},{61,5},{61,6},{62,0},{62,1},{62,2}

#define QCELP_RATE_8THR_BITMAP \
{76,3},{66,0},{67,0},{68,0},{76,2},{69,0},{70,0},{71,0},\
{76,1},{72,0},{73,0},{74,0},{76,0},{75,0},{16,1},{16,0},\
{65,3},{65,2},{65,1},{65,0}


/**
 * Position of the bitmapping data for each pkt type in
 * the big REFERENCE FRAME array
 */

#define QCELP_FULLPKT_REFERENCE_POS 0
#define QCELP_HALFPKT_REFERENCE_POS 266
#define QCELP_4THRPKT_REFERENCE_POS 390
#define QCELP_8THRPKT_REFERENCE_POS 444

static const QCELPBitmap QCELP_REFERENCE_FRAME[]={QCELP_RATE_FULL_BITMAP,
                                                  QCELP_RATE_HALF_BITMAP,
                                                  QCELP_RATE_4THR_BITMAP,
                                                  QCELP_RATE_8THR_BITMAP};

/**
 * Position of the transmission codes inside the universal frame.
 */

#define QCELP_CBSIGN0_POS 0
#define QCELP_CBGAIN0_POS 16
#define QCELP_CINDEX0_POS 32
#define QCELP_PLAG0_POS   48
#define QCELP_PFRAC0_POS  52
#define QCELP_PGAIN0_POS  56
#define QCELP_LSPV0_POS   60
#define QCELP_RSRVD_POS   65    /*!< on all but rate 1/2 packets */
#define QCELP_LSP0_POS    66    /*!< only in rate 1/8 packets    */
#define QCELP_CBSEED_POS  76    /*!< only in rate 1/8 packets    */

/* rest is currently unused */

static const int   qcelp_cumulative_gainloss[]={0,1,2,6};
static const float qcelp_cumulative_pitchsaturation[]={0.9,0.6,0.3,0.0};

static const float qcelp_fullrate_ccodebook[]=
{
    0.10,-0.65,-0.59, 0.12, 1.10, 0.34,-1.34, 1.57,
    1.04,-0.84,-0.34,-1.15, 0.23,-1.01, 0.03, 0.45,
   -1.01,-0.16,-0.59, 0.28,-0.45, 1.34,-0.67, 0.22,
    0.61,-0.29, 2.26,-0.26,-0.55,-1.79, 1.57,-0.51,
   -2.20,-0.93,-0.37, 0.60, 1.18, 0.74,-0.48,-0.95,
   -1.81, 1.11, 0.36,-0.52,-2.15, 0.78,-1.12, 0.39,
   -0.17,-0.47,-2.23, 0.19, 0.12,-0.98,-1.42, 1.30,
    0.54,-1.27, 0.21,-0.12, 0.39,-0.48, 0.12, 1.28,
    0.06,-1.67, 0.82,-1.02,-0.79, 0.55,-0.44, 0.48,
   -0.20,-0.53, 0.08,-0.61, 0.11,-0.70,-1.57,-1.68,
    0.20,-0.56,-0.74, 0.78, 0.33,-0.63,-1.73,-0.02,
   -0.75,-0.53,-1.46, 0.77, 0.66,-0.29, 0.09,-0.75,
    0.65, 1.19,-0.43, 0.76, 2.33, 0.98, 1.25,-1.56,
   -0.27, 0.78,-0.09, 1.70, 1.76, 1.43,-1.48,-0.07,
    0.27,-1.36, 0.05, 0.27, 0.18, 1.39, 2.04, 0.07,
   -1.84,-1.97, 0.52,-0.03, 0.78,-1.89, 0.08,-0.65
};

static const float qcelp_halfrate_ccodebook[]=
{
    0.0, -2.0,  0.0, -1.5,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0, -1.5, -1.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  2.5,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  2.0,  0.0,
    0.0,  1.5,  1.0,  0.0,  1.5,  2.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  1.5,  0.0,  0.0,
   -1.5,  1.5,  0.0,  0.0, -1.0,  0.0,  1.5,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -2.5,  0.0,
    0.0,  0.0,  0.0,  1.5,  0.0,  0.0,  0.0,  1.5,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  2.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  1.5,  3.0, -1.5, -2.0,  0.0, -1.5, -1.5,
    1.5, -1.5,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,
    0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0
};
