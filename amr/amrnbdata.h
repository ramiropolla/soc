/*
 * AMR narrowband data and definitions
 * Copyright (c) 2006-2007 Robert Swain
 * Copyright (c) 2009 Colin McQuillan
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
 * @file libavcodec/amrnbdata.h
 * AMR narrowband data and definitions
 */

#ifndef AVCODEC_AMRNBDATA_H
#define AVCODEC_AMRNBDATA_H

#include <stdint.h>

#include "libavutil/common.h"      /* offsetof */
#include "libavutil/mathematics.h" /* M_PI */

#define AMR_BLOCK_SIZE              160   ///< Samples per frame
#define AMR_SUBFRAME_SIZE            40   ///< Samples per subframe
#define AMR_SAMPLE_BOUND        32768.0   ///< Threshold for synthesis overflow

/**
 * Scale from constructed speech to [-1,1]
 *
 * AMR is designed to produce 16-bit PCM samples (3GPP TS 26.090 4.2) but
 * upscales by two (section 6.2.2).
 *
 * Fundamentally, this scale is determined by energy_mean through
 * the fixed vector contribution to the excitation vector.
 */
#define AMR_SAMPLE_SCALE  (2.0/32768.0)


/** Frame type (Table 1a in 3GPP TS 26.101) */
enum Mode {
    MODE_475 = 0,                         ///< 4.75 kbit/s
    MODE_515,                             ///< 5.15 kbit/s
    MODE_59,                              ///< 5.90 kbit/s
    MODE_67,                              ///< 6.70 kbit/s
    MODE_74,                              ///< 7.40 kbit/s
    MODE_795,                             ///< 7.95 kbit/s
    MODE_102,                             ///< 10.2 kbit/s
    MODE_122,                             ///< 12.2 kbit/s
    MODE_DTX,                             ///< silent frame
    N_MODES,                              ///< number of modes
    NO_DATA = 15                          ///< no transmission
};

#define LP_FILTER_ORDER 10        ///< linear predictive coding filter order

/**
 * Bit-order table type
 */
typedef struct AMROrder {
    unsigned index : 6;           ///< index in (uint16_t *)AMRNBFrame
    unsigned len   : 4;           ///< number of bits
    unsigned shift : 4;           ///< magnitude of shift to apply
    unsigned rev   : 1;           ///< reverse uint16_t and right shift instead of left shift
} AMROrder;


/**
 * AMRNB unpacked data subframe
 */
typedef struct {
    uint16_t p_lag;      ///< index to decode the pitch lag
    uint16_t p_gain;     ///< index to decode the pitch gain
    uint16_t fixed_gain; ///< index to decode the fixed gain factor, for MODE_122 and MODE_795
    uint16_t pulses[10]; ///< pulses: 10 for MODE_122, 7 for MODE_102, and index and sign for others
} AMRNBSubframe;

/**
 * AMRNB SID frame parameters
 */
typedef struct {
    uint16_t ref_vector; ///< index of reference vector
    uint16_t energy;     ///< index of logarithmic frame energy
} AMRNBSIDFrame;

/**
 * AMRNB unpacked data frame
 */
typedef struct {
    uint16_t lsf[5];           ///< lsf parameters: 5 parameters for MODE_122, only 3 for other modes
    AMRNBSubframe subframe[4]; ///< unpacked data for each subframe
    AMRNBSIDFrame sid;
} AMRNBFrame;


// The following order* tables are used to convert AMR frame parameters to and
// from a bitstream. See 3GPP TS 26.101 for more information.

#define AMR_BIT(field, l,s,r)                  {offsetof(AMRNBFrame, field) >> 1, l,s,r}
/** Specify an LSF parameter bit */
#define AMR_LSF(variable, l,s,r)               AMR_BIT(lsf[variable], l,s,r)
/** Specify a subframe-specific bit */
#define AMR_OF(frame_num, variable, l,s,r)     AMR_BIT(subframe[frame_num].variable, l,s,r)
/** Specify a pitch gain bit */
#define AMR_PGAIN(frame_num, l,s,r)            AMR_OF(frame_num, p_gain, l,s,r)
/** Specify a fixed gain bit */
#define AMR_FGAIN(frame_num, l,s,r)            AMR_OF(frame_num, fixed_gain, l,s,r)
/** Specify a pitch lag bit */
#define AMR_PLAG(frame_num, l,s,r)             AMR_OF(frame_num, p_lag, l,s,r)
/** Specify a pulse bit */
#define AMR_PULSE(frame_num, pulse_id, l,s,r)  AMR_OF(frame_num, pulses[pulse_id], l,s,r)
/** Specify an SID reference vector bit */
#define AMR_SVECTOR(l,s,r)                     AMR_BIT(sid.ref_vector, l,s,r)
/** Specify an SID energy index bit */
#define AMR_SENERGY(l,s,r)                     AMR_BIT(sid.energy,     l,s,r)

static AMROrder order_MODE_475[] = {
    AMR_LSF  (  0,8, 0,0), AMR_LSF  (  1,8, 0,0), AMR_PLAG (0  ,6, 2,0),
    AMR_PLAG (1  ,2, 2,0), AMR_PLAG (2  ,2, 2,0), AMR_PLAG (3  ,2, 2,0),
    AMR_PGAIN(0  ,4,12,1), AMR_PGAIN(2  ,4,12,1), AMR_LSF  (  2,2, 4,0),
    AMR_LSF  (  2,1, 2,0), AMR_LSF  (  2,1, 0,0), AMR_PGAIN(2  ,4, 8,1),
    AMR_PLAG (0  ,2, 0,0), AMR_PGAIN(0  ,4, 8,1), AMR_PULSE(0,1,2, 0,0),
    AMR_LSF  (  2,1, 6,0), AMR_LSF  (  2,1, 3,0), AMR_LSF  (  2,1, 1,0),
    AMR_PLAG (1  ,2, 0,0), AMR_PULSE(1,1,2, 0,0), AMR_PLAG (2  ,2, 0,0),
    AMR_PULSE(2,1,2, 0,0), AMR_PLAG (3  ,2, 0,0), AMR_PULSE(3,1,2, 0,0),
    AMR_PULSE(0,0,2, 4,0), AMR_PULSE(0,0,2, 1,0), AMR_PULSE(1,0,2, 4,0),
    AMR_PULSE(1,0,2, 1,0), AMR_PULSE(2,0,2, 4,0), AMR_PULSE(2,0,2, 1,0),
    AMR_PULSE(3,0,2, 4,0), AMR_PULSE(3,0,2, 1,0), AMR_PULSE(0,0,1, 3,0),
    AMR_PULSE(1,0,1, 3,0), AMR_PULSE(2,0,1, 3,0), AMR_PULSE(3,0,1, 3,0),
    AMR_PULSE(0,0,1, 0,0), AMR_PULSE(1,0,1, 0,0), AMR_PULSE(2,0,1, 0,0),
    AMR_PULSE(3,0,1, 0,0), AMR_PULSE(0,0,1, 6,0), AMR_PULSE(1,0,1, 6,0),
    AMR_PULSE(2,0,1, 6,0), AMR_PULSE(3,0,1, 6,0)
};

static AMROrder order_MODE_515[] = {
    AMR_LSF  (  0,8, 8,1), AMR_LSF  (  1,8, 8,1), AMR_PLAG (0  ,5, 3,0),
    AMR_PLAG (1  ,1, 3,0), AMR_PLAG (2  ,1, 3,0), AMR_PLAG (3  ,1, 3,0),
    AMR_PGAIN(0  ,3,13,1), AMR_PGAIN(1  ,3,13,1), AMR_PGAIN(2  ,3,13,1),
    AMR_PGAIN(3  ,3,13,1), AMR_PGAIN(0  ,1, 3,0), AMR_PGAIN(1  ,1, 3,0),
    AMR_PGAIN(2  ,1, 3,0), AMR_PGAIN(3  ,1, 3,0), AMR_PLAG (0  ,1, 2,0),
    AMR_PLAG (1  ,1, 2,0), AMR_PLAG (2  ,1, 2,0), AMR_PLAG (3  ,1, 2,0),
    AMR_LSF  (  2,1, 4,0), AMR_PGAIN(0  ,1, 4,0), AMR_PGAIN(1  ,1, 4,0),
    AMR_PGAIN(2  ,1, 4,0), AMR_PGAIN(3  ,1, 4,0), AMR_PLAG (0  ,1, 1,0),
    AMR_PLAG (1  ,1, 1,0), AMR_PLAG (2  ,1, 1,0), AMR_LSF  (  2,1, 5,0),
    AMR_LSF  (  2,1, 2,0), AMR_LSF  (  2,1, 0,0), AMR_PGAIN(0  ,1, 5,0),
    AMR_PGAIN(1  ,1, 5,0), AMR_PGAIN(2  ,1, 5,0), AMR_PGAIN(3  ,1, 5,0),
    AMR_LSF  (  2,1, 1,0), AMR_PLAG (0  ,1, 0,0), AMR_PLAG (1  ,1, 0,0),
    AMR_PLAG (2  ,1, 0,0), AMR_PLAG (3  ,1, 1,0), AMR_LSF  (  2,1, 3,0),
    AMR_LSF  (  2,1, 6,0), AMR_PLAG (3  ,1, 0,0), AMR_PULSE(0,1,2,14,1),
    AMR_PULSE(1,1,2,14,1), AMR_PULSE(2,1,1, 0,0), AMR_PULSE(0,0,1, 2,0),
    AMR_PULSE(1,0,1, 2,0), AMR_PULSE(2,0,1, 2,0), AMR_PULSE(3,0,1, 2,0),
    AMR_PULSE(2,1,1, 1,0), AMR_PULSE(3,1,2,14,1), AMR_PULSE(0,0,1, 1,0),
    AMR_PULSE(1,0,1, 1,0), AMR_PULSE(2,0,1, 1,0), AMR_PULSE(3,0,1, 1,0),
    AMR_PULSE(0,0,1, 5,0), AMR_PULSE(1,0,1, 5,0), AMR_PULSE(0,0,1, 4,0),
    AMR_PULSE(1,0,1, 4,0), AMR_PULSE(2,0,2, 4,0), AMR_PULSE(3,0,2, 4,0),
    AMR_PULSE(0,0,1, 6,0), AMR_PULSE(1,0,1, 6,0), AMR_PULSE(2,0,1, 6,0),
    AMR_PULSE(3,0,1, 6,0), AMR_PULSE(0,0,1, 0,0), AMR_PULSE(1,0,1, 0,0),
    AMR_PULSE(2,0,1, 0,0), AMR_PULSE(3,0,1, 0,0), AMR_PULSE(0,0,1, 3,0),
    AMR_PULSE(1,0,1, 3,0), AMR_PULSE(2,0,1, 3,0), AMR_PULSE(3,0,1, 3,0)
};

static AMROrder order_MODE_59[] = {
    AMR_LSF  (  0,2, 6,0), AMR_LSF  (  0,2, 2,0), AMR_LSF  (  0,1, 4,0),
    AMR_LSF  (  0,2, 0,0), AMR_LSF  (  0,1, 5,0), AMR_LSF  (  1,1, 3,0),
    AMR_LSF  (  1,1, 1,0), AMR_LSF  (  1,2, 7,0), AMR_LSF  (  1,2, 4,0),
    AMR_LSF  (  1,1, 2,0), AMR_LSF  (  1,1, 6,0), AMR_LSF  (  1,1, 0,0),
    AMR_PLAG (0  ,1, 5,0), AMR_PLAG (2  ,1, 5,0), AMR_PLAG (0  ,1, 4,0),
    AMR_PLAG (2  ,1, 4,0), AMR_PLAG (0  ,1, 6,0), AMR_PLAG (2  ,1, 6,0),
    AMR_PLAG (0  ,1, 7,0), AMR_PLAG (2  ,1, 7,0), AMR_PLAG (0  ,1, 3,0),
    AMR_PLAG (2  ,1, 3,0), AMR_PLAG (1  ,1, 3,0), AMR_PLAG (3  ,1, 3,0),
    AMR_PGAIN(0  ,1, 0,0), AMR_PGAIN(1  ,1, 0,0), AMR_PGAIN(2  ,1, 0,0),
    AMR_PGAIN(3  ,1, 0,0), AMR_PLAG (0  ,1, 2,0), AMR_PLAG (2  ,1, 2,0),
    AMR_PLAG (1  ,1, 2,0), AMR_PLAG (3  ,1, 2,0), AMR_PGAIN(0  ,1, 1,0),
    AMR_PGAIN(1  ,1, 1,0), AMR_PGAIN(2  ,1, 1,0), AMR_PGAIN(3  ,1, 1,0),
    AMR_PLAG (1  ,1, 1,0), AMR_PLAG (3  ,1, 1,0), AMR_PLAG (0  ,1, 1,0),
    AMR_PLAG (2  ,1, 1,0), AMR_PLAG (0  ,1, 0,0), AMR_PLAG (2  ,1, 0,0),
    AMR_PGAIN(0  ,1, 2,0), AMR_PGAIN(1  ,1, 2,0), AMR_PGAIN(2  ,1, 2,0),
    AMR_PGAIN(3  ,1, 2,0), AMR_PGAIN(0  ,1, 3,0), AMR_PGAIN(1  ,1, 3,0),
    AMR_PGAIN(2  ,1, 3,0), AMR_PGAIN(3  ,1, 3,0), AMR_PGAIN(0  ,1, 4,0),
    AMR_PGAIN(1  ,1, 4,0), AMR_PGAIN(2  ,1, 4,0), AMR_PGAIN(3  ,1, 4,0),
    AMR_LSF  (  2,1, 6,0), AMR_LSF  (  2,1, 4,0), AMR_LSF  (  2,2,12,1),
    AMR_LSF  (  2,2, 7,1), AMR_LSF  (  2,1, 5,0), AMR_LSF  (  2,1, 1,0),
    AMR_PULSE(3,1,1, 0,0), AMR_PULSE(0,1,1, 1,0), AMR_PULSE(2,1,1, 1,0),
    AMR_PULSE(3,1,1, 1,0), AMR_PULSE(1,1,2, 0,0), AMR_PULSE(0,1,1, 0,0),
    AMR_PULSE(2,1,1, 0,0), AMR_LSF  (  2,1, 0,0), AMR_PGAIN(0  ,1, 5,0),
    AMR_PGAIN(1  ,1, 5,0), AMR_PGAIN(2  ,1, 5,0), AMR_PGAIN(3  ,1, 5,0),
    AMR_PLAG (1  ,1, 0,0), AMR_PLAG (3  ,1, 0,0), AMR_PULSE(0,0,1, 2,0),
    AMR_PULSE(1,0,1, 2,0), AMR_PULSE(2,0,1, 2,0), AMR_PULSE(3,0,1, 2,0),
    AMR_PULSE(0,0,1, 3,0), AMR_PULSE(1,0,1, 3,0), AMR_PULSE(2,0,1, 3,0),
    AMR_PULSE(3,0,1, 3,0), AMR_PULSE(0,0,1, 6,0), AMR_PULSE(1,0,1, 6,0),
    AMR_PULSE(2,0,1, 6,0), AMR_PULSE(3,0,1, 6,0), AMR_PULSE(0,0,1, 7,0),
    AMR_PULSE(1,0,1, 7,0), AMR_PULSE(2,0,1, 7,0), AMR_PULSE(3,0,1, 7,0),
    AMR_PULSE(0,0,1, 8,0), AMR_PULSE(1,0,1, 8,0), AMR_PULSE(2,0,1, 8,0),
    AMR_PULSE(3,0,1, 8,0), AMR_PULSE(0,0,1, 0,0), AMR_PULSE(1,0,1, 0,0),
    AMR_PULSE(2,0,1, 0,0), AMR_PULSE(3,0,1, 0,0), AMR_PULSE(0,0,1, 1,0),
    AMR_PULSE(1,0,1, 1,0), AMR_PULSE(2,0,1, 1,0), AMR_PULSE(3,0,1, 1,0),
    AMR_PULSE(0,0,1, 4,0), AMR_PULSE(1,0,1, 4,0), AMR_PULSE(2,0,1, 4,0),
    AMR_PULSE(3,0,1, 4,0), AMR_PULSE(0,0,1, 5,0), AMR_PULSE(1,0,1, 5,0),
    AMR_PULSE(2,0,1, 5,0), AMR_PULSE(3,0,1, 5,0)
};

static AMROrder order_MODE_67[] = {
    AMR_LSF  (  0,2, 6,0), AMR_LSF  (  0,2,11,1), AMR_LSF  (  0,2, 1,0),
    AMR_LSF  (  1,1, 3,0), AMR_LSF  (  0,1, 0,0), AMR_LSF  (  0,1, 5,0),
    AMR_LSF  (  1,2, 7,0), AMR_LSF  (  1,1, 5,0), AMR_LSF  (  1,1, 1,0),
    AMR_LSF  (  1,1, 4,0), AMR_LSF  (  1,1, 2,0), AMR_LSF  (  1,1, 6,0),
    AMR_PLAG (0  ,1, 5,0), AMR_PLAG (2  ,1, 5,0), AMR_PLAG (0  ,1, 4,0),
    AMR_PLAG (2  ,1, 4,0), AMR_PLAG (0  ,1, 6,0), AMR_PLAG (2  ,1, 6,0),
    AMR_PLAG (0  ,1, 7,0), AMR_PLAG (2  ,1, 7,0), AMR_PLAG (0  ,1, 3,0),
    AMR_PLAG (2  ,1, 3,0), AMR_LSF  (  1,1, 0,0), AMR_PLAG (1  ,1, 3,0),
    AMR_PLAG (3  ,1, 3,0), AMR_PLAG (1  ,1, 2,0), AMR_PLAG (3  ,1, 2,0),
    AMR_PLAG (0  ,1, 2,0), AMR_PLAG (2  ,1, 2,0), AMR_PLAG (1  ,1, 1,0),
    AMR_PLAG (3  ,1, 1,0), AMR_PGAIN(0  ,1, 6,0), AMR_PGAIN(1  ,1, 6,0),
    AMR_PGAIN(2  ,1, 6,0), AMR_PGAIN(3  ,1, 6,0), AMR_PLAG (0  ,1, 1,0),
    AMR_PLAG (2  ,1, 1,0), AMR_PGAIN(0  ,1, 3,0), AMR_PGAIN(1  ,1, 3,0),
    AMR_PGAIN(2  ,1, 3,0), AMR_PGAIN(3  ,1, 3,0), AMR_PGAIN(0  ,1, 2,0),
    AMR_PGAIN(1  ,1, 2,0), AMR_PGAIN(2  ,1, 2,0), AMR_PGAIN(3  ,1, 2,0),
    AMR_PLAG (1  ,1, 0,0), AMR_PLAG (3  ,1, 0,0), AMR_PLAG (0  ,1, 0,0),
    AMR_PLAG (2  ,1, 0,0), AMR_LSF  (  2,1, 6,0), AMR_LSF  (  2,1, 2,0),
    AMR_PGAIN(0  ,1, 1,0), AMR_PGAIN(1  ,1, 1,0), AMR_PGAIN(2  ,1, 1,0),
    AMR_PGAIN(3  ,1, 1,0), AMR_LSF  (  2,2, 3,0), AMR_LSF  (  2,2, 7,1),
    AMR_LSF  (  2,1, 5,0), AMR_LSF  (  2,2, 0,0), AMR_PGAIN(0  ,1, 4,0),
    AMR_PGAIN(1  ,1, 4,0), AMR_PGAIN(2  ,1, 4,0), AMR_PGAIN(3  ,1, 4,0),
    AMR_PULSE(0,1,1, 0,0), AMR_PULSE(1,1,1, 0,0), AMR_PULSE(2,1,1, 0,0),
    AMR_PULSE(3,1,1, 0,0), AMR_PGAIN(0  ,1, 0,0), AMR_PGAIN(1  ,1, 0,0),
    AMR_PGAIN(2  ,1, 0,0), AMR_PGAIN(3  ,1, 0,0), AMR_PULSE(0,1,1, 1,0),
    AMR_PULSE(1,1,1, 1,0), AMR_PULSE(2,1,1, 1,0), AMR_PULSE(3,1,1, 1,0),
    AMR_PGAIN(3  ,1, 5,0), AMR_PGAIN(2  ,1, 5,0), AMR_PGAIN(1  ,1, 5,0),
    AMR_PGAIN(0  ,1, 5,0), AMR_PULSE(0,1,1, 2,0), AMR_PULSE(1,1,1, 2,0),
    AMR_PULSE(2,1,1, 2,0), AMR_PULSE(3,1,1, 2,0), AMR_PULSE(0,0,1, 2,0),
    AMR_PULSE(1,0,1, 2,0), AMR_PULSE(2,0,1, 2,0), AMR_PULSE(3,0,1, 2,0),
    AMR_PULSE(0,0,1, 5,0), AMR_PULSE(1,0,1, 5,0), AMR_PULSE(2,0,1, 5,0),
    AMR_PULSE(3,0,1, 5,0), AMR_PULSE(0,0,1, 6,0), AMR_PULSE(1,0,1, 6,0),
    AMR_PULSE(2,0,1, 6,0), AMR_PULSE(3,0,1, 6,0), AMR_PULSE(0,0,1, 9,0),
    AMR_PULSE(1,0,1, 9,0), AMR_PULSE(2,0,1, 9,0), AMR_PULSE(3,0,1, 9,0),
    AMR_PULSE(0,0,1,10,0), AMR_PULSE(1,0,1,10,0), AMR_PULSE(2,0,1,10,0),
    AMR_PULSE(3,0,1,10,0), AMR_PULSE(0,0,1, 0,0), AMR_PULSE(1,0,1, 0,0),
    AMR_PULSE(2,0,1, 0,0), AMR_PULSE(3,0,1, 0,0), AMR_PULSE(0,0,1, 1,0),
    AMR_PULSE(1,0,1, 1,0), AMR_PULSE(2,0,1, 1,0), AMR_PULSE(3,0,1, 1,0),
    AMR_PULSE(0,0,1, 3,0), AMR_PULSE(1,0,1, 3,0), AMR_PULSE(2,0,1, 3,0),
    AMR_PULSE(3,0,1, 3,0), AMR_PULSE(0,0,1, 4,0), AMR_PULSE(1,0,1, 4,0),
    AMR_PULSE(2,0,1, 4,0), AMR_PULSE(3,0,1, 4,0), AMR_PULSE(0,0,1, 7,0),
    AMR_PULSE(1,0,1, 7,0), AMR_PULSE(2,0,1, 7,0), AMR_PULSE(3,0,1, 7,0),
    AMR_PULSE(0,0,1, 8,0), AMR_PULSE(1,0,1, 8,0), AMR_PULSE(2,0,1, 8,0),
    AMR_PULSE(3,0,1, 8,0)
};

static AMROrder order_MODE_74[] = {
    AMR_LSF  (  0,8, 0,0), AMR_LSF  (  1,9, 0,0), AMR_PLAG (0  ,1, 7,0),
    AMR_PLAG (2  ,1, 7,0), AMR_PLAG (0  ,1, 6,0), AMR_PLAG (2  ,1, 6,0),
    AMR_PLAG (0  ,1, 5,0), AMR_PLAG (2  ,1, 5,0), AMR_PLAG (0  ,1, 4,0),
    AMR_PLAG (2  ,1, 4,0), AMR_PLAG (0  ,1, 3,0), AMR_PLAG (2  ,1, 3,0),
    AMR_PGAIN(0  ,1, 6,0), AMR_PGAIN(1  ,1, 6,0), AMR_PGAIN(2  ,1, 6,0),
    AMR_PGAIN(3  ,1, 6,0), AMR_PGAIN(0  ,1, 5,0), AMR_PGAIN(1  ,1, 5,0),
    AMR_PGAIN(2  ,1, 5,0), AMR_PGAIN(3  ,1, 5,0), AMR_PGAIN(0  ,1, 3,0),
    AMR_PGAIN(1  ,1, 3,0), AMR_PGAIN(2  ,1, 3,0), AMR_PGAIN(3  ,1, 3,0),
    AMR_PGAIN(0  ,1, 2,0), AMR_PGAIN(1  ,1, 2,0), AMR_PGAIN(2  ,1, 2,0),
    AMR_PGAIN(3  ,1, 2,0), AMR_PLAG (1  ,1, 4,0), AMR_PLAG (3  ,1, 4,0),
    AMR_PLAG (1  ,1, 3,0), AMR_PLAG (3  ,1, 3,0), AMR_LSF  (  2,3, 2,0),
    AMR_LSF  (  2,3, 6,0), AMR_PLAG (0  ,1, 2,0), AMR_PLAG (1  ,1, 2,0),
    AMR_PLAG (2  ,1, 2,0), AMR_PLAG (3  ,1, 2,0), AMR_PGAIN(0  ,1, 1,0),
    AMR_PGAIN(1  ,1, 1,0), AMR_PGAIN(2  ,1, 1,0), AMR_PGAIN(3  ,1, 1,0),
    AMR_LSF  (  2,1, 5,0), AMR_LSF  (  2,2, 0,0), AMR_PULSE(0,1,1, 0,0),
    AMR_PULSE(1,1,1, 0,0), AMR_PULSE(2,1,1, 0,0), AMR_PULSE(3,1,1, 0,0),
    AMR_PGAIN(0  ,1, 0,0), AMR_PGAIN(1  ,1, 0,0), AMR_PGAIN(2  ,1, 0,0),
    AMR_PGAIN(3  ,1, 0,0), AMR_PULSE(0,1,1, 1,0), AMR_PULSE(1,1,1, 1,0),
    AMR_PULSE(2,1,1, 1,0), AMR_PULSE(3,1,1, 1,0), AMR_PULSE(0,1,1, 2,0),
    AMR_PULSE(1,1,1, 2,0), AMR_PGAIN(0  ,1, 4,0), AMR_PGAIN(1  ,1, 4,0),
    AMR_PGAIN(2  ,1, 4,0), AMR_PGAIN(3  ,1, 4,0), AMR_PULSE(2,1,1, 2,0),
    AMR_PULSE(3,1,1, 2,0), AMR_PULSE(0,1,1, 3,0), AMR_PULSE(1,1,1, 3,0),
    AMR_PULSE(2,1,1, 3,0), AMR_PULSE(3,1,1, 3,0), AMR_PLAG (0  ,2, 0,0),
    AMR_PLAG (1  ,2, 0,0), AMR_PLAG (2  ,2, 0,0), AMR_PLAG (3  ,2, 0,0),
    AMR_PULSE(0,0,6, 0,0), AMR_PULSE(1,0,6, 0,0), AMR_PULSE(2,0,6, 0,0),
    AMR_PULSE(3,0,6, 0,0), AMR_PULSE(0,0,1,12,0), AMR_PULSE(1,0,1,12,0),
    AMR_PULSE(2,0,1,12,0), AMR_PULSE(3,0,1,12,0), AMR_PULSE(0,0,1,11,0),
    AMR_PULSE(1,0,1,11,0), AMR_PULSE(2,0,1,11,0), AMR_PULSE(3,0,1,11,0),
    AMR_PULSE(0,0,1,10,0), AMR_PULSE(1,0,1,10,0), AMR_PULSE(2,0,1,10,0),
    AMR_PULSE(3,0,1,10,0), AMR_PULSE(0,0,1, 9,0), AMR_PULSE(1,0,1, 9,0),
    AMR_PULSE(2,0,1, 9,0), AMR_PULSE(3,0,1, 9,0), AMR_PULSE(0,0,1, 8,0),
    AMR_PULSE(1,0,1, 8,0), AMR_PULSE(2,0,1, 8,0), AMR_PULSE(3,0,1, 8,0),
    AMR_PULSE(0,0,1, 7,0), AMR_PULSE(1,0,1, 7,0), AMR_PULSE(2,0,1, 7,0),
    AMR_PULSE(3,0,1, 7,0), AMR_PULSE(0,0,1, 6,0), AMR_PULSE(1,0,1, 6,0),
    AMR_PULSE(2,0,1, 6,0), AMR_PULSE(3,0,1, 6,0)
};

static AMROrder order_MODE_795[] = {
    AMR_LSF  (  0,7, 9,1), AMR_LSF  (  1,1, 3,0), AMR_LSF  (  1,1, 1,0),
    AMR_LSF  (  1,2, 7,0), AMR_LSF  (  1,2, 4,0), AMR_LSF  (  1,1, 2,0),
    AMR_LSF  (  1,1, 6,0), AMR_LSF  (  1,1, 0,0), AMR_LSF  (  2,1, 6,0),
    AMR_LSF  (  2,1, 4,0), AMR_LSF  (  2,2,12,1), AMR_LSF  (  2,2, 7,1),
    AMR_LSF  (  2,1, 5,0), AMR_FGAIN(0  ,1, 4,0), AMR_FGAIN(1  ,1, 4,0),
    AMR_FGAIN(2  ,1, 4,0), AMR_FGAIN(3  ,1, 4,0), AMR_FGAIN(0  ,1, 3,0),
    AMR_FGAIN(1  ,1, 3,0), AMR_FGAIN(2  ,1, 3,0), AMR_FGAIN(3  ,1, 3,0),
    AMR_FGAIN(0  ,1, 2,0), AMR_FGAIN(1  ,1, 2,0), AMR_FGAIN(2  ,1, 2,0),
    AMR_FGAIN(3  ,1, 2,0), AMR_PGAIN(0  ,1, 3,0), AMR_PGAIN(1  ,1, 3,0),
    AMR_PGAIN(2  ,1, 3,0), AMR_PGAIN(3  ,1, 3,0), AMR_PGAIN(0  ,1, 2,0),
    AMR_PGAIN(1  ,1, 2,0), AMR_PGAIN(2  ,1, 2,0), AMR_PGAIN(3  ,1, 2,0),
    AMR_PLAG (0  ,1, 7,0), AMR_PLAG (2  ,1, 7,0), AMR_PLAG (0  ,1, 6,0),
    AMR_PLAG (2  ,1, 6,0), AMR_PLAG (0  ,1, 5,0), AMR_PLAG (2  ,1, 5,0),
    AMR_PLAG (0  ,1, 4,0), AMR_PLAG (2  ,1, 4,0), AMR_PLAG (0  ,1, 3,0),
    AMR_PLAG (2  ,1, 3,0), AMR_PLAG (1  ,1, 5,0), AMR_PLAG (3  ,1, 5,0),
    AMR_PLAG (1  ,1, 4,0), AMR_PLAG (3  ,1, 4,0), AMR_PLAG (1  ,1, 3,0),
    AMR_PLAG (3  ,1, 3,0), AMR_FGAIN(0  ,1, 1,0), AMR_FGAIN(1  ,1, 1,0),
    AMR_FGAIN(2  ,1, 1,0), AMR_FGAIN(3  ,1, 1,0), AMR_PLAG (0  ,1, 2,0),
    AMR_PLAG (2  ,1, 2,0), AMR_PLAG (1  ,1, 2,0), AMR_PLAG (3  ,1, 2,0),
    AMR_LSF  (  0,2, 7,1), AMR_LSF  (  2,2, 0,0), AMR_PLAG (0  ,1, 1,0),
    AMR_PLAG (2  ,1, 1,0), AMR_PLAG (0  ,1, 0,0), AMR_PLAG (2  ,1, 0,0),
    AMR_PLAG (1  ,1, 1,0), AMR_PLAG (3  ,1, 1,0), AMR_PLAG (1  ,1, 0,0),
    AMR_PLAG (3  ,1, 0,0), AMR_PGAIN(0  ,1, 1,0), AMR_PGAIN(1  ,1, 1,0),
    AMR_PGAIN(2  ,1, 1,0), AMR_PGAIN(3  ,1, 1,0), AMR_FGAIN(0  ,1, 0,0),
    AMR_FGAIN(1  ,1, 0,0), AMR_FGAIN(2  ,1, 0,0), AMR_FGAIN(3  ,1, 0,0),
    AMR_PGAIN(0  ,1, 0,0), AMR_PGAIN(1  ,1, 0,0), AMR_PGAIN(2  ,1, 0,0),
    AMR_PGAIN(3  ,1, 0,0), AMR_PULSE(2,1,3,13,1), AMR_PULSE(0,0,1, 1,0),
    AMR_PULSE(1,0,1, 1,0), AMR_PULSE(2,0,1, 1,0), AMR_PULSE(3,0,1, 1,0),
    AMR_PULSE(0,0,1, 4,0), AMR_PULSE(1,0,1, 4,0), AMR_PULSE(2,0,1, 4,0),
    AMR_PULSE(3,0,1, 4,0), AMR_PULSE(0,0,1, 7,0), AMR_PULSE(1,0,1, 7,0),
    AMR_PULSE(2,0,1, 7,0), AMR_PULSE(3,0,1, 7,0), AMR_PULSE(0,0,1,11,0),
    AMR_PULSE(1,0,1,11,0), AMR_PULSE(2,0,1,11,0), AMR_PULSE(3,0,1,11,0),
    AMR_PULSE(2,1,1, 3,0), AMR_PULSE(3,1,4,12,1), AMR_PULSE(1,1,4,12,1),
    AMR_PULSE(0,1,4,12,1), AMR_PULSE(0,0,1, 0,0), AMR_PULSE(0,0,2,12,1),
    AMR_PULSE(0,0,1, 5,0), AMR_PULSE(0,0,1, 8,0), AMR_PULSE(0,0,1,12,0),
    AMR_PULSE(1,0,1, 0,0), AMR_PULSE(1,0,2,12,1), AMR_PULSE(1,0,1, 5,0),
    AMR_PULSE(1,0,1, 8,0), AMR_PULSE(1,0,1,12,0), AMR_PULSE(2,0,1, 0,0),
    AMR_PULSE(2,0,2,12,1), AMR_PULSE(2,0,1, 5,0), AMR_PULSE(2,0,1, 8,0),
    AMR_PULSE(2,0,1,12,0), AMR_PULSE(3,0,1, 0,0), AMR_PULSE(3,0,2,12,1),
    AMR_PULSE(3,0,1, 5,0), AMR_PULSE(3,0,1, 8,0), AMR_PULSE(3,0,1,12,0),
    AMR_PULSE(0,0,1, 6,0), AMR_PULSE(1,0,1, 6,0), AMR_PULSE(2,0,1, 6,0),
    AMR_PULSE(3,0,1, 6,0), AMR_PULSE(0,0,1,10,0), AMR_PULSE(1,0,1,10,0),
    AMR_PULSE(2,0,1,10,0), AMR_PULSE(3,0,1,10,0), AMR_PULSE(0,0,1, 9,0),
    AMR_PULSE(1,0,1, 9,0), AMR_PULSE(2,0,1, 9,0), AMR_PULSE(3,0,1, 9,0)
};

static AMROrder order_MODE_102[] = {
    AMR_LSF  (  0,8, 8,1), AMR_LSF  (  1,9, 7,1), AMR_PLAG (0  ,6, 2,0),
    AMR_PLAG (2  ,6, 2,0), AMR_PLAG (1  ,2, 3,0), AMR_PLAG (3  ,2, 3,0),
    AMR_PGAIN(0  ,1, 6,0), AMR_PGAIN(0  ,2, 2,0), AMR_PGAIN(1  ,1, 6,0),
    AMR_PGAIN(1  ,2, 2,0), AMR_PGAIN(2  ,1, 6,0), AMR_PGAIN(2  ,2, 2,0),
    AMR_PGAIN(3  ,1, 6,0), AMR_PGAIN(3  ,2, 2,0), AMR_PLAG (0  ,2, 0,0),
    AMR_PLAG (2  ,2, 0,0), AMR_PLAG (1  ,2, 1,0), AMR_PLAG (3  ,2, 1,0),
    AMR_PGAIN(0  ,1, 5,0), AMR_PGAIN(1  ,1, 5,0), AMR_PGAIN(2  ,1, 5,0),
    AMR_PGAIN(3  ,1, 5,0), AMR_LSF  (  2,1, 6,0), AMR_LSF  (  2,1, 2,0),
    AMR_LSF  (  2,2, 3,0), AMR_LSF  (  2,2, 7,1), AMR_LSF  (  2,1, 5,0),
    AMR_LSF  (  2,2, 0,0), AMR_PULSE(0,3,1, 0,0), AMR_PULSE(0,2,1, 0,0),
    AMR_PULSE(0,1,1, 0,0), AMR_PULSE(0,0,1, 0,0), AMR_PULSE(1,3,1, 0,0),
    AMR_PULSE(1,2,1, 0,0), AMR_PULSE(1,1,1, 0,0), AMR_PULSE(1,0,1, 0,0),
    AMR_PULSE(2,3,1, 0,0), AMR_PULSE(2,2,1, 0,0), AMR_PULSE(2,1,1, 0,0),
    AMR_PULSE(2,0,1, 0,0), AMR_PULSE(3,3,1, 0,0), AMR_PULSE(3,2,1, 0,0),
    AMR_PULSE(3,1,1, 0,0), AMR_PULSE(3,0,1, 0,0), AMR_PGAIN(0  ,1, 1,0),
    AMR_PGAIN(0  ,1, 4,0), AMR_PGAIN(0  ,1, 0,0), AMR_PGAIN(1  ,1, 1,0),
    AMR_PGAIN(1  ,1, 4,0), AMR_PGAIN(1  ,1, 0,0), AMR_PGAIN(2  ,1, 1,0),
    AMR_PGAIN(2  ,1, 4,0), AMR_PGAIN(2  ,1, 0,0), AMR_PGAIN(3  ,1, 1,0),
    AMR_PGAIN(3  ,1, 4,0), AMR_PGAIN(3  ,1, 0,0), AMR_PLAG (1  ,1, 0,0),
    AMR_PLAG (3  ,1, 0,0), AMR_PULSE(1,4,2, 8,0), AMR_PULSE(1,5,2, 6,1),
    AMR_PULSE(1,5,1, 7,0), AMR_PULSE(1,4,1, 7,0), AMR_PULSE(1,5,1, 5,0),
    AMR_PULSE(1,4,2, 9,1), AMR_PULSE(1,5,1, 6,0), AMR_PULSE(1,6,2, 5,0),
    AMR_PULSE(1,5,1, 4,0), AMR_PULSE(1,6,1, 3,0), AMR_PULSE(1,4,1, 4,0),
    AMR_PULSE(1,6,1, 4,0), AMR_PULSE(1,4,1, 3,0), AMR_PULSE(1,5,1, 3,0),
    AMR_PULSE(2,4,2, 8,0), AMR_PULSE(2,5,2, 6,1), AMR_PULSE(2,5,1, 7,0),
    AMR_PULSE(2,4,1, 7,0), AMR_PULSE(2,5,1, 5,0), AMR_PULSE(2,4,2, 9,1),
    AMR_PULSE(2,5,1, 6,0), AMR_PULSE(2,6,2, 5,0), AMR_PULSE(2,5,1, 4,0),
    AMR_PULSE(2,6,1, 3,0), AMR_PULSE(2,4,1, 4,0), AMR_PULSE(2,6,1, 4,0),
    AMR_PULSE(2,4,1, 3,0), AMR_PULSE(2,5,1, 3,0), AMR_PULSE(3,4,2, 8,0),
    AMR_PULSE(3,5,2, 6,1), AMR_PULSE(3,5,1, 7,0), AMR_PULSE(3,4,1, 7,0),
    AMR_PULSE(3,5,1, 5,0), AMR_PULSE(3,4,2, 9,1), AMR_PULSE(3,5,1, 6,0),
    AMR_PULSE(3,6,2, 5,0), AMR_PULSE(3,5,1, 4,0), AMR_PULSE(3,6,1, 3,0),
    AMR_PULSE(3,4,1, 4,0), AMR_PULSE(3,6,1, 4,0), AMR_PULSE(3,4,1, 3,0),
    AMR_PULSE(3,5,1, 3,0), AMR_PULSE(0,4,2, 8,0), AMR_PULSE(0,5,2, 6,1),
    AMR_PULSE(0,5,1, 7,0), AMR_PULSE(0,4,1, 7,0), AMR_PULSE(0,5,1, 5,0),
    AMR_PULSE(0,4,2, 9,1), AMR_PULSE(0,5,1, 6,0), AMR_PULSE(0,6,2, 5,0),
    AMR_PULSE(0,5,1, 4,0), AMR_PULSE(0,6,1, 3,0), AMR_PULSE(0,4,1, 4,0),
    AMR_PULSE(0,6,1, 4,0), AMR_PULSE(0,4,1, 3,0), AMR_PULSE(0,5,1, 3,0),
    AMR_PULSE(3,6,1, 2,0), AMR_PULSE(3,4,1, 0,0), AMR_PULSE(3,5,1, 0,0),
    AMR_PULSE(3,6,1, 0,0), AMR_PULSE(3,4,1, 2,0), AMR_PULSE(3,6,1, 1,0),
    AMR_PULSE(3,4,1, 1,0), AMR_PULSE(3,5,2, 1,0), AMR_PULSE(2,6,1, 2,0),
    AMR_PULSE(2,4,1, 0,0), AMR_PULSE(2,5,1, 0,0), AMR_PULSE(2,6,1, 0,0),
    AMR_PULSE(2,4,1, 2,0), AMR_PULSE(2,6,1, 1,0), AMR_PULSE(2,4,1, 1,0),
    AMR_PULSE(2,5,2, 1,0), AMR_PULSE(1,6,1, 2,0), AMR_PULSE(1,4,1, 0,0),
    AMR_PULSE(1,5,1, 0,0), AMR_PULSE(1,6,1, 0,0), AMR_PULSE(1,4,1, 2,0),
    AMR_PULSE(1,6,1, 1,0), AMR_PULSE(1,4,1, 1,0), AMR_PULSE(1,5,2, 1,0),
    AMR_PULSE(0,6,1, 2,0), AMR_PULSE(0,4,1, 0,0), AMR_PULSE(0,5,1, 0,0),
    AMR_PULSE(0,6,1, 0,0), AMR_PULSE(0,4,1, 2,0), AMR_PULSE(0,6,1, 1,0),
    AMR_PULSE(0,4,1, 1,0), AMR_PULSE(0,5,2, 1,0)
};

static AMROrder order_MODE_122[] = {
    AMR_LSF  (  0,7, 0,0), AMR_LSF  (  1,8, 0,0), AMR_LSF  (  2,1, 0,0),
    AMR_LSF  (  2,8, 1,0), AMR_LSF  (  3,5, 3,0), AMR_PLAG (0  ,1, 8,0),
    AMR_PLAG (2  ,1, 8,0), AMR_PLAG (0  ,1, 7,0), AMR_PLAG (2  ,1, 7,0),
    AMR_PLAG (0  ,1, 6,0), AMR_PLAG (2  ,1, 6,0), AMR_PLAG (0  ,1, 5,0),
    AMR_PLAG (2  ,1, 5,0), AMR_PLAG (0  ,1, 4,0), AMR_PLAG (2  ,1, 4,0),
    AMR_PLAG (0  ,1, 3,0), AMR_PLAG (2  ,1, 3,0), AMR_PLAG (0  ,1, 2,0),
    AMR_PLAG (2  ,1, 2,0), AMR_PLAG (0  ,1, 1,0), AMR_PLAG (2  ,1, 1,0),
    AMR_PLAG (0  ,1, 0,0), AMR_PLAG (2  ,1, 0,0), AMR_PGAIN(0  ,1, 3,0),
    AMR_PGAIN(1  ,1, 3,0), AMR_PGAIN(2  ,1, 3,0), AMR_PGAIN(3  ,1, 3,0),
    AMR_PGAIN(0  ,1, 2,0), AMR_PGAIN(1  ,1, 2,0), AMR_PGAIN(2  ,1, 2,0),
    AMR_PGAIN(3  ,1, 2,0), AMR_PGAIN(0  ,1, 1,0), AMR_PGAIN(1  ,1, 1,0),
    AMR_PGAIN(2  ,1, 1,0), AMR_PGAIN(3  ,1, 1,0), AMR_FGAIN(0  ,1, 4,0),
    AMR_FGAIN(1  ,1, 4,0), AMR_FGAIN(2  ,1, 4,0), AMR_FGAIN(3  ,1, 4,0),
    AMR_FGAIN(0  ,1, 3,0), AMR_FGAIN(1  ,1, 3,0), AMR_FGAIN(2  ,1, 3,0),
    AMR_FGAIN(3  ,1, 3,0), AMR_FGAIN(0  ,1, 2,0), AMR_FGAIN(1  ,1, 2,0),
    AMR_FGAIN(2  ,1, 2,0), AMR_FGAIN(3  ,1, 2,0), AMR_PLAG (1  ,1, 5,0),
    AMR_PLAG (3  ,1, 5,0), AMR_PLAG (1  ,1, 4,0), AMR_PLAG (3  ,1, 4,0),
    AMR_PLAG (1  ,1, 3,0), AMR_PLAG (3  ,1, 3,0), AMR_PLAG (1  ,1, 2,0),
    AMR_PLAG (3  ,1, 2,0), AMR_PLAG (1  ,1, 1,0), AMR_PLAG (3  ,1, 1,0),
    AMR_LSF  (  3,3, 0,0), AMR_LSF  (  4,4, 2,0), AMR_PGAIN(0  ,1, 0,0),
    AMR_PGAIN(1  ,1, 0,0), AMR_PGAIN(2  ,1, 0,0), AMR_PGAIN(3  ,1, 0,0),
    AMR_FGAIN(0  ,1, 1,0), AMR_FGAIN(1  ,1, 1,0), AMR_FGAIN(2  ,1, 1,0),
    AMR_FGAIN(3  ,1, 1,0), AMR_PULSE(0,0,1, 3,0), AMR_PULSE(1,0,1, 3,0),
    AMR_PULSE(2,0,1, 3,0), AMR_PULSE(3,0,1, 3,0), AMR_PULSE(0,1,1, 3,0),
    AMR_PULSE(1,1,1, 3,0), AMR_PULSE(2,1,1, 3,0), AMR_PULSE(3,1,1, 3,0),
    AMR_FGAIN(0  ,1, 0,0), AMR_FGAIN(1  ,1, 0,0), AMR_FGAIN(2  ,1, 0,0),
    AMR_FGAIN(3  ,1, 0,0), AMR_PULSE(0,2,1, 3,0), AMR_PULSE(1,2,1, 3,0),
    AMR_PULSE(2,2,1, 3,0), AMR_PULSE(3,2,1, 3,0), AMR_PULSE(0,3,1, 3,0),
    AMR_PULSE(1,3,1, 3,0), AMR_PULSE(2,3,1, 3,0), AMR_PULSE(3,3,1, 3,0),
    AMR_PULSE(0,4,1, 3,0), AMR_PULSE(1,4,1, 3,0), AMR_PULSE(2,4,1, 3,0),
    AMR_PULSE(3,4,1, 3,0), AMR_LSF  (  4,2, 0,0), AMR_PULSE(0,0,3,13,1),
    AMR_PULSE(0,1,3,13,1), AMR_PULSE(0,2,3,13,1), AMR_PULSE(0,3,3,13,1),
    AMR_PULSE(0,4,3,13,1), AMR_PULSE(1,0,3,13,1), AMR_PULSE(1,1,3,13,1),
    AMR_PULSE(1,2,3,13,1), AMR_PULSE(1,3,3,13,1), AMR_PULSE(1,4,3,13,1),
    AMR_PULSE(2,0,3,13,1), AMR_PULSE(2,1,3,13,1), AMR_PULSE(2,2,3,13,1),
    AMR_PULSE(2,3,3,13,1), AMR_PULSE(2,4,3,13,1), AMR_PULSE(3,0,3,13,1),
    AMR_PULSE(3,1,3,13,1), AMR_PULSE(3,2,3,13,1), AMR_PULSE(3,3,3,13,1),
    AMR_PULSE(3,4,3,13,1), AMR_PULSE(0,5,3,13,1), AMR_PULSE(0,6,3,13,1),
    AMR_PULSE(0,7,3,13,1), AMR_PULSE(0,8,3,13,1), AMR_PULSE(0,9,3,13,1),
    AMR_PULSE(1,5,3,13,1), AMR_PULSE(1,6,3,13,1), AMR_PULSE(1,7,3,13,1),
    AMR_PULSE(1,8,3,13,1), AMR_PULSE(1,9,3,13,1), AMR_PULSE(2,5,3,13,1),
    AMR_PULSE(2,6,3,13,1), AMR_PULSE(2,7,3,13,1), AMR_PULSE(2,8,3,13,1),
    AMR_PULSE(2,9,3,13,1), AMR_PULSE(3,5,3,13,1), AMR_PULSE(3,6,3,13,1),
    AMR_PULSE(3,7,3,13,1), AMR_PULSE(3,8,3,13,1), AMR_PULSE(3,9,3,13,1),
    AMR_PLAG (1  ,1, 0,0), AMR_PLAG (3  ,1, 0,0)
};

static AMROrder order_MODE_DTX[] = {
    AMR_SVECTOR(  3, 0,0), AMR_LSF    (0,8, 0,0), AMR_LSF    (1,9, 0,0),
    AMR_LSF    (2,9, 0,0), AMR_SENERGY(  6, 0,0)
};

/**
 * position of the bitmapping data for each packet type in
 * the AMRNBFrame
 */
static const AMROrder * const amr_unpacking_bitmaps_per_mode[9] = {
    order_MODE_475,
    order_MODE_515,
    order_MODE_59,
    order_MODE_67,
    order_MODE_74,
    order_MODE_795,
    order_MODE_102,
    order_MODE_122,
    order_MODE_DTX,
};

/** number of runs of consecutive bits for each mode */
static const uint8_t mode_runs[N_MODES] = {
    FF_ARRAY_ELEMS(order_MODE_475),
    FF_ARRAY_ELEMS(order_MODE_515),
    FF_ARRAY_ELEMS(order_MODE_59),
    FF_ARRAY_ELEMS(order_MODE_67),
    FF_ARRAY_ELEMS(order_MODE_74),
    FF_ARRAY_ELEMS(order_MODE_795),
    FF_ARRAY_ELEMS(order_MODE_102),
    FF_ARRAY_ELEMS(order_MODE_122),
    FF_ARRAY_ELEMS(order_MODE_DTX)
};

/** number of bits for each mode */
static const uint8_t mode_bits[N_MODES] = { 95, 103, 118, 134, 148, 159, 204, 244, 35 };

/**
 * Values for the lsp vector from the 4th subframe of the
 * previous subframe values.
 *
 * @note: Taken from Decoder_amr_reset in Q15 using val/1000
 */
static const int8_t lsp_sub4_init[LP_FILTER_ORDER] = { 30, 26, 21, 15, 8, 0, -8, -15, -21, -26 };

/**
 * Mean lsp values.
 *
 * @note: Taken from Decoder_amr_reset in Q15
 */
static const int16_t lsp_avg_init[LP_FILTER_ORDER] = {
    1384, 2077, 3420, 5108, 6742, 8122, 9863, 11092, 12714, 13701
};

// LSF tables

// These are stored as integers to save space. The values are taken from
// q_plsf_3.tab and q_plsf_5.tab in 3GPP TS 26.090.

static const int16_t lsf_3_3_MODE_515[128][4] = {
{  419,  163,  -30, -262},{ -455, -789,-1430, -721},{ 1006,  664,  269,   25},
{  619,  260,  183,   96},{ -968,-1358, -388,  135},{ -693,  835,  456,  154},
{ 1105,  703,  569,  363},{ 1625, 1326,  985,  748},{ -220,  219,   76, -208},
{-1455,-1662,   49,  149},{ -964, -172, -752, -336},{  625,  209, -250,  -66},
{-1017, -838,   -2,  317},{-2168,-1485, -138,  123},{-1876,-2099, -521,   85},
{ -967, -366, -695, -881},{ -921,-1011, -763, -949},{ -124, -256, -352, -660},
{  178,  463,  354,  304},{-1744, -591, -282,   79},{-2249,  175,  867,  499},
{ -138, -180, -181,  -21},{-2291,-1241, -460, -520},{ -771,  451,  -10, -308},
{  271,  -65,    4,  214},{ -279, -435,  -43, -348},{ -670,   35,  -65, -211},
{  806,  535,   85,  297},{   57,  239,  722,  493},{  225,  661,  840,  547},
{ -540, -376,   14,  349},{  469,  721,  331,  162},{ -544, -752,  -62,  -10},
{  398,  -88,  724,  701},{  -19, -533,  -94,  601},{  136,  -71, -681, -747},
{ -166, -344,  261,  -50},{  161,  -52,  485,  337},{-1675,   50,  190,  -93},
{-2282, -231, -194,  -82},{  -95, -595, -154,  128},{  894,  501,  588,  457},
{ -345,  206,  122,  110},{ -631, -227, -569,    3},{  408,  239,  397,  226},
{ -197,   -2,  128,  491},{ 1281,  904,  292,  215},{  538,  306,  259,  509},
{ -677,-1047,   13,  321},{ -679, -588, -358, -212},{ -558,  243,  646,  479},
{  486,  342,  634,  532},{  107,  802,  331,  136},{ -112, -398,-1031, -286},
{ -326, -705,  288,  272},{ 1299, 1144, 1178,  860},{ -423,  121, -385, -148},
{ -295, -302, -834, -819},{   16,  -24, -201, -476},{  555,   91, -245,  294},
{  -38, -379, -962,-1221},{-1191,-1518, -273, -395},{ -390,-1013, -645,  573},
{-1843,-1030,  505,  468},{  744,  947,  609,  493},{ -689,-1172, -628, -135},
{-1026,  195,  411,  196},{ 1582, 1147,  575,  337},{-1239, -777, -648, -142},
{  595,  825,  967,  735},{-1206, -970,  -81, -342},{ -745,   13,  -72,  375},
{  454,   19, 1407,  921},{-1647, -172,  861,  562},{  928, 1537, 1063,  740},
{-2472, -952,  264,   82},{ -502, -965,-1334,  123},{  867, 1236,  534,  171},
{-2320, -460,  780,  363},{-1190, -617,  252,  -61},{ -174,   34, 1011,  788},
{-2333,  247,  423,  153},{  -16, -355,  262,  449},{-1576,-1073, -544, -371},
{ -615, -305, 1051,  805},{  687,  528,    6, -182},{  935,  875, 1002,  809},
{  199,  257,  126,   76},{ -584,-1138,  599,  556},{-1105,-1391,-1591, -519},
{ -977,-1325,  108,  347},{ -722, -975,  365,  101},{ -145,  681,  249, -153},
{    0, -334, -570,  159},{  412,  285, -336, -617},{ -953, -966,  887,  689},
{-1251,   84, -185, -398},{ -592,  433, 1044,  653},{   85,  329,  -40,  361},
{ -433, -705,  466,  574},{ -154,  654,  592,  290},{ -167,   72,  349,  175},
{  674,  297,  977,  720},{ 1235, 1204,  757,  488},{ -400, -269,  538,  372},
{-1350,-1387,-1194,  -91},{ 1262,  876,  775,  700},{ -599,  -38, -430, -722},
{ 1976, 1630,  991,  608},{  111,  276, -226,  -96},{ -947, -388,  -11,   -7},
{ -303, -531, -839,  338},{ 1734, 1710, 1405, 1013},{ -516, -855, -645,  210},
{ -688, -416,  513,  230},{ -822, -637,-1146, -320},{ -952, -658, -694,  183},
{ -114, -623,  818,  674},{ -191, -204,  731,  635},{   51, 1221,  883,  576},
{ -954, -431,  826,  598},{ -342, -755, -900, -407},{-1126, -354, -206, -512},
{ -547, -810, -357, -620},{   66,  515,  -73, -410},{ -872, -945,-1444,-1227},
{  191,  -17, -544, -231},{-1540, -544, -901, -886}
};

static const int16_t lsf_3_1_MODE_795[512][3] = {
{ -890,-1550,-2541},{ -819, -970,  175},{ -826,-1234, -762},
{ -599,  -22,  634},{ -811, -987, -902},{ -323,  203,   26},
{ -383, -235, -781},{ -399, 1262,  906},{ -932,-1399,-1380},
{ -624,   93,   87},{ -414, -539, -691},{   37,  633,  510},
{ -387, -476,-1330},{  399,   66,  263},{ -407,  -49, -335},
{ -417, 1041, 1865},{ -779,-1089,-1440},{ -746, -858,  832},
{ -581, -759, -371},{ -673, -506, 2088},{ -560, -634,-1179},
{  271,  241,   14},{ -438, -244, -397},{  463, 1202, 1047},
{ -606, -797,-1438},{  -51, -323,  481},{ -224, -584, -527},
{  494,  881,  682},{ -433, -306,-1002},{  554,  659,  222},
{  171, -160, -353},{  681, 1798, 1565},{ -852,-1181,-1695},
{ -336, -666,  114},{ -581, -756, -744},{ -195,  375,  497},
{ -465, -804,-1098},{  154,  282, -131},{  -50, -191, -719},
{  323,  732, 1542},{ -722, -819,-1404},{  105, -250,  185},
{ -178, -502, -742},{  321,  510, 1111},{ -323, -567, -966},
{  127,  484,  338},{ -160,   52, -338},{  732, 1367, 1554},
{ -626, -802,-1696},{ -286, -586,  676},{ -695, -343, -370},
{ -490,  295, 1893},{ -630, -574,-1014},{  -80,  645,  -69},
{   -6, -318, -364},{  782, 1450, 1038},{ -313, -733,-1395},
{  120,   60,  477},{ -264, -585, -123},{  711, 1245,  633},
{  -91, -355,-1016},{  771,  758,  261},{  253,   81, -474},
{  930, 2215, 1720},{ -808,-1099,-1925},{ -560, -782,  169},
{ -804,-1074, -188},{ -626,  -55, 1405},{ -694, -716,-1194},
{ -660,  354,  329},{ -514,  -55, -543},{  366, 1033, 1182},
{ -658, -959,-1357},{  -55, -184,   93},{ -605, -286, -662},
{  404,  449,  827},{ -286, -350,-1263},{  628,  306,  227},
{  -16,  147, -623},{  186,  923, 2146},{ -674, -890,-1606},
{ -443, -228,  339},{ -369, -790, -409},{  231,   86, 1469},
{ -448, -581,-1061},{  594,  450, -177},{ -124, -170, -447},
{  671, 1159, 1404},{ -476, -667,-1511},{  -77, -138,  716},
{ -177, -372, -381},{  451,  934,  915},{ -250, -432, -822},
{  272,  828,  446},{   26,   19,  -31},{  698, 1692, 2168},
{ -646, -977,-1924},{ -179, -473,  268},{ -379, -745, -691},
{   11,  127, 1033},{ -488, -917, -825},{   61,  323,  135},
{  147, -145, -686},{  685,  786, 1682},{ -506, -848,-1297},
{   35,   90,  222},{  -23, -346, -670},{  455,  591, 1287},
{ -203, -593,-1086},{  652,  352,  437},{   39,   63, -457},
{  841, 1265, 2105},{ -520, -882,-1584},{ -328, -711, 1421},
{ -596, -342,  -70},{  209,  173, 1928},{ -423, -598, -921},
{  421,  605,  -38},{   -2, -245, -127},{  896, 1969, 1135},
{ -379, -518,-1579},{  173,  118,  753},{  -55, -381,  -52},
{  985, 1021,  753},{   -2, -291, -891},{  753,  992,  423},
{  264,  131, -196},{  895, 2274, 2543},{ -635,-1088,-2499},
{ -529, -982,  526},{ -764, -830, -548},{ -436,  316,  599},
{ -675, -940, -746},{  -57,  236,  -11},{ -201,  -81, -798},
{   16,  845, 1558},{ -737, -985,-1212},{ -468,   17,  290},
{ -279, -584, -700},{  183,  822,  705},{ -265, -492,-1187},
{  421,  152,  468},{ -390,  166, -268},{   39, 1550, 1868},
{ -635, -966,-1571},{ -453, -492,  910},{ -284,-1027,  -75},
{ -181, -133, 1852},{ -445, -624,-1174},{  420,  367,  -49},
{ -389, -212, -169},{  707, 1073, 1208},{ -539, -710,-1449},
{   83, -163,  484},{ -236, -543, -355},{  338, 1175,  814},
{ -246, -309, -958},{  606,  760,   60},{  166,   -8, -163},
{ -306, 1849, 2563},{ -747,-1025,-1783},{ -419, -446,  209},
{ -718, -566, -534},{ -506,  693,  857},{ -463, -697,-1082},
{  325,  431, -206},{  -15,   -8, -763},{  545,  919, 1518},
{ -611, -783,-1313},{  256,  -55,  208},{ -165, -348, -662},
{  321,  680,  930},{ -326, -429, -951},{  484,  446,  570},
{ -197,   72,  -73},{  909, 1455, 1741},{ -563, -737,-1974},
{ -124, -416,  718},{ -478, -404, -314},{  -16,  446, 1636},
{ -551, -537, -750},{  -58,  638,  214},{   55, -185, -271},
{ 1148, 1301, 1212},{ -483, -671,-1264},{  117,  285,  543},
{ -204, -391, -111},{  513, 1538,  854},{ -114, -190, -978},
{  877,  595,  464},{  260,  260, -311},{  748, 2283, 2216},
{ -517, -945,-2171},{ -326, -708,  378},{ -812, -691, -232},
{ -560,  687, 1409},{ -732, -690, -836},{ -359,  645,  386},
{ -265,   62, -678},{  145, 1644, 1208},{ -555, -988,-1233},
{  -78,   14,  114},{ -327, -358, -489},{  392,  677,  697},
{ -201, -236,-1140},{  693,  449,  178},{ -243,  256, -433},
{  611, 1385, 2456},{ -612, -901,-1464},{ -307,  -17,  499},
{ -315, -667, -254},{  256,  428, 1463},{ -486, -422,-1056},
{  655,  370,   18},{ -102, -185, -276},{  755, 1578, 1335},
{ -488, -603,-1418},{  182,  -93,  870},{  -73, -458, -348},
{  835,  862,  957},{ -282, -333, -746},{  547,  839,  428},
{  273,  -89,   13},{  940, 1708, 2576},{ -418,-1084,-1758},
{  -44, -358,  259},{ -497, -643, -560},{   99,  557,  961},
{ -421, -766, -917},{  295,  326,  184},{  175,   15, -626},
{  532,  878, 1981},{ -443, -768,-1275},{  221,  156,  268},
{   39, -363, -505},{  695,  772, 1140},{ -162, -459, -912},
{  709,  444,  658},{   25,  303, -312},{ 1268, 1410, 1715},
{ -297, -766,-1836},{ -263, -108, 1070},{ -406,  -13, -129},
{   57,  438, 2734},{ -374, -487, -835},{  304,  696,  164},
{  104, -235,    5},{ 1611, 1900, 1399},{ -229, -582,-1325},
{  405,  192,  817},{  -87, -438,  111},{ 1028, 1199,  993},
{   68, -175, -934},{ 1033, 1117,  451},{  478,  200, -248},
{ 2127, 2696, 2042},{ -835,-1323,-2131},{ -799, -692,  466},
{ -812,-1032, -469},{ -622,  288,  920},{ -701, -841,-1070},
{ -411,  512,    8},{ -390,  -91, -744},{  -30, 1043, 1161},
{ -822,-1148,-1156},{ -294,  -46,  110},{ -411, -374, -678},
{  214,  531,  668},{ -406, -420,-1194},{  487,  232,  303},
{ -318,   91, -472},{  123, 1232, 2445},{ -722, -952,-1495},
{ -738, -675, 1332},{ -543, -606, -211},{  -95,  -98, 1508},
{ -549, -514,-1193},{  473,  211,   73},{ -288, -112, -389},
{  537, 1332, 1258},{ -567, -755,-1545},{   71, -283,  632},
{ -170, -481, -493},{  681, 1002,  817},{ -356, -331, -877},
{  419,  706,  346},{  241,  -34, -326},{  377, 1950, 1883},
{ -727,-1075,-1625},{ -233, -543,  116},{ -524, -806, -585},
{  -73,  478,  729},{ -288, -925,-1143},{  173,  447,  -52},
{   68, -229, -606},{  449,  529, 1797},{ -591, -875,-1363},
{  183, -144,  324},{ -103, -452, -666},{  623,  488, 1176},
{ -238, -511,-1004},{  326,  552,  458},{  136,  108, -319},
{  626, 1343, 1883},{ -490, -646,-1730},{ -186, -449,  984},
{ -738,  -76, -170},{ -550,  755, 2560},{ -496, -510, -947},
{  210,  694,  -52},{   84, -322, -199},{ 1090, 1625, 1224},
{ -376, -603,-1396},{  343,   74,  632},{ -175, -502,  -32},
{  972, 1332,  734},{   52, -295,-1113},{ 1065,  918,  160},
{  393,  107, -397},{ 1214, 2649, 1741},{ -632,-1201,-1891},
{ -719, -277,  353},{ -651, -880, -122},{ -211,  209, 1338},
{ -562, -714,-1059},{ -208,  388,  159},{ -320,  -61, -551},
{  293, 1092, 1443},{ -648, -865,-1253},{  -49, -143,  305},
{ -401, -227, -585},{  561,  532,  927},{ -117, -443,-1188},
{  507,  436,  292},{  -79,  233, -458},{  671, 1025, 2396},
{ -633, -842,-1525},{ -308, -286,  640},{ -373, -621, -407},
{  418,  253, 1305},{ -315, -581,-1137},{  572,  685, -281},
{   61,  -68, -371},{  991, 1101, 1498},{ -493, -683,-1362},
{  -47,  164,  704},{ -256, -314, -268},{  631,  949, 1052},
{ -118, -348, -833},{   68, 1180,  568},{  152,  117,   34},
{ 1113, 1902, 2239},{ -601, -959,-1706},{ -143, -489,  480},
{ -332, -655, -574},{   54,  353, 1192},{ -462, -652, -796},
{  150,  549,  112},{  195, -111, -515},{  679, 1108, 1647},
{ -558, -749,-1217},{   -9,  272,  341},{  -53, -265, -535},
{  489,  843, 1298},{ -120, -482,-1032},{  632,  543,  408},
{  179,  306, -526},{ 1124, 1464, 2244},{ -417, -786,-1562},
{ -224, -384, 1364},{ -377, -459,  -25},{  385,  489, 2174},
{ -332, -651, -829},{  544,  553,   61},{   22, -113,  -89},
{ 1128, 1725, 1524},{ -216, -373,-1653},{  161,  316,  908},
{ -165, -222,  -67},{ 1362, 1175,  789},{   73, -252, -767},
{  738,  932,  616},{  362,  246, -126},{  787, 2654, 3027},
{ -691,-1106,-2190},{ -565, -588,  524},{ -590, -979, -490},
{ -263,  397,  982},{ -577, -837, -945},{  -22,  435,  -49},
{ -190, -118, -629},{  -88, 1240, 1513},{ -636,-1051,-1019},
{ -291,  189,  259},{ -257, -470, -629},{  145,  945,  894},
{ -326, -364,-1094},{  543,  260,  630},{ -202,  189, -209},
{  357, 1379, 2091},{ -569,-1075,-1449},{ -714, -239,  919},
{ -420, -705,  -84},{ -109, -114, 2407},{ -413, -529,-1177},
{  482,  368,  131},{ -186,  -72, -131},{  861, 1255, 1220},
{ -611, -658,-1341},{  227, -121,  631},{ -176, -489, -218},
{  745, 1175,  957},{ -321, -148, -936},{  671,  966,  216},
{  340,   -3, -143},{  469, 1848, 2437},{ -729, -961,-1683},
{ -213, -254,  321},{ -511, -438, -521},{ -126,  725,  903},
{ -340, -685,-1032},{  316,  480,   20},{   23,  -89, -551},
{  353, 1051, 1789},{ -544, -757,-1364},{  298,  -25,  436},
{ -100, -392, -519},{  467,  754, 1078},{ -210, -398,-1078},
{  620,  658,  630},{   33,  147, -178},{  921, 1687, 1921},
{ -325, -528,-1978},{    2, -285,  910},{ -371, -490, -230},
{    0,  597, 2010},{ -496, -395, -834},{   37,  945,  245},
{  181, -160, -144},{ 1481, 1373, 1357},{ -355, -601,-1270},
{  298,  322,  672},{ -193, -336,   77},{ 1089, 1533,  922},
{  177,  -39,-1125},{  996,  781,  536},{  456,  366, -432},
{ 1415, 2440, 2279},{ -466, -758,-2325},{ -303, -509,  387},
{ -727, -557,   66},{ -145,  643, 1248},{ -544, -676, -916},
{ -225,  862,  588},{ -152,   40, -533},{  423, 1423, 1558},
{ -572, -843,-1145},{ -128,   85,  461},{ -238, -257, -584},
{  605,  748,  861},{   24, -202,-1409},{  797,  487,  303},
{ -181,  364, -182},{  616, 1378, 2942},{ -494, -852,-1441},
{ -292,   61,  812},{  -84, -723, -182},{  555,  532, 1506},
{ -365, -493,-1057},{  822,  588,   11},{  -14,  -18, -230},
{ 1001, 1401, 1451},{ -474, -569,-1292},{  302,   62, 1062},
{  -70, -376, -222},{  982,  974, 1149},{ -196, -234, -795},
{  479, 1098,  499},{  362,   58,   70},{ 1147, 2069, 2857},
{ -487, -878,-1824},{   73, -288,  348},{ -358, -500, -508},
{  199,  721, 1242},{  -78, -697, -795},{  361,  536,  196},
{  374,  110, -735},{  847, 1051, 1896},{ -366, -713,-1182},
{  315,  320,  429},{   72, -215, -450},{  759,  886, 1363},
{  -30, -428, -834},{  861,  627,  796},{  118,  468, -279},
{ 1355, 1883, 1893},{ -188, -642,-1612},{   63, -175, 1198},
{ -418, -211,   51},{  414,  587, 2601},{ -234, -557, -858},
{  424,  889,  222},{  136, -101,   83},{ 1413, 2278, 1383},
{  -84, -445,-1389},{  414,  313, 1045},{   29, -343,   65},
{ 1552, 1647,  980},{  183,  -91, -829},{ 1273, 1413,  360},
{  553,  272, -107},{ 1587, 3149, 2603}
};

static const int16_t lsf_3_1[256][3] = {
{    6,   82, -131},{  154,  -56, -735},{  183,  -65, -265},
{    9, -210, -361},{  113,  718, 1817},{ 1010, 1214, 1573},
{  857, 1333, 2276},{  827, 1568, 1933},{  717, 1989, 2206},
{  838, 1172, 1823},{  721, 1000, 2154},{  286,  476, 1509},
{ -247, -531,  230},{  147,  -82,  569},{   26, -177, -944},
{  -27, -273,  692},{ -164, -264, -183},{  224,  790, 1039},
{  899,  946,  601},{  485,  771, 1150},{  524,  677,  903},
{ -140,  375,  778},{  410,  676,  429},{  301,  530, 1009},
{  719,  646,   38},{  226,  367,   40},{  145,  -45, -505},
{  290,  121, -121},{  302,  127,  166},{ -124, -383, -956},
{ -358, -455, -977},{  715,  878,  894},{  978,  923,  211},
{  477,  272,   64},{  188,  -78,   17},{ -143,  -65,   38},
{  643,  586,  621},{ -134, -426, -651},{  347,  545, 2820},
{ 1188, 2726, 2442},{  142,  -80, 1735},{  283,  130,  461},
{ -262, -399,-1145},{ -411,  155,  430},{  329,  375,  779},
{   53, -226, -139},{ -129, -236, 1682},{  285,  744, 1327},
{  738,  697, 1664},{  312,  409,  266},{  325,  720,  135},
{    1,  221,  453},{    8,  203,  145},{  299,  640,  760},
{   29,  468,  638},{  103,  429,  379},{  420,  954,  932},
{ 1326, 1210, 1258},{  704, 1012, 1152},{ -166, -444, -266},
{ -316, -130, -376},{  191, 1151, 1904},{ -240, -543,-1260},
{ -112,  268, 1207},{   70, 1062, 1583},{  278, 1360, 1574},
{ -258, -272, -768},{   19,  563, 2240},{   -3, -265,  135},
{ -295, -591, -388},{  140,  354, -206},{ -260, -504, -795},
{ -433, -718,-1319},{  109,  331,  962},{ -429,  -87,  652},
{ -296,  426, 1019},{ -239,  775,  851},{  489, 1334, 1073},
{ -334, -332,   25},{  543, 1206, 1807},{  326,   61,  727},
{  578,  849, 1405},{ -208, -277,  329},{ -152,   64,  669},
{ -434, -678, -727},{ -454,  -71,  251},{  605,  480,  254},
{ -482,   11,  996},{ -289,  395,  486},{  722, 1049, 1440},
{  -30, -316, -786},{ -106, -115, -619},{  861, 1474, 1412},
{ 1055, 1366, 1184},{  812, 1237,  925},{   42, -251, -576},
{  342,  141, -454},{ -168,  -80, 1359},{ -342, -656,-1763},
{  100,  821,  725},{  990,  747,  800},{  332,  440,  568},
{  663,  379,  852},{  112,  165, -369},{  597,  910,  282},
{   -8,  834, 1281},{ -352,  572,  695},{  462, 2246, 1806},
{  345,  190, 1374},{  416,  915, 2166},{  168,  -82,  280},
{ -516, -446,  840},{   47,  533,   44},{ -362, -711,-1143},
{   22,  193, 1472},{  -85,  233, 1813},{  -62,  579, 1504},
{  550,  944, 1749},{  723,  650, 1148},{  972,  884, 1395},
{ -425,  643,    0},{ 1000,  952, 1098},{  249, 1446,  672},
{ -334,  -87, 2172},{ -554, 1882, 2672},{  140, 1826, 1853},
{  920, 1749, 2590},{ 1076, 1933, 2038},{ -137, -443,-1555},
{ 1269, 1174,  468},{ -493, -122, 1521},{ -451, 1033, 1214},
{  482, 1695, 1118},{  815,  649,  384},{ -446, -692,  107},
{ -319, -605, -118},{ -207, -505,  525},{ -468,  -12, 2736},
{   75, 1934, 1305},{  880, 2358, 2267},{ 1285, 1575, 2004},
{  -48, -304,-1186},{ -435, -461, -251},{ -366, -404, -547},
{ -289, -605, -597},{ -538, -810, -165},{ -120,    3,  356},
{  639, 1241, 1502},{   96,  177,  750},{ -435, -585,-1174},
{ -356,  109,  -79},{ -485,  288, 2005},{    9, 1116,  731},
{  880, 2134,  946},{ -265, 1585, 1065},{ 1157, 1210,  843},
{ -498, -668,  431},{  374,  321, -229},{ 1440, 2101, 1381},
{  449,  461, 1155},{ -105,   39, -384},{ -263,  367,  182},
{ -371, -660,  773},{ -188, 1151,  971},{ 1333, 1632, 1435},
{  774, 1267, 1221},{ -482, -832,-1489},{ -237, -210,  860},
{  890, 1615, 1064},{  472, 1062, 1192},{  185, 1077,  989},
{ -568, -992,-1704},{ -449, -902,-2043},{ -142, -377, -458},
{ -210, -554,-1029},{  -11, 1133, 2265},{ -329, -675, -893},
{ -250,  657, 1187},{  519, 1510, 1779},{  520,  539, 1403},
{  527, 1421, 1302},{ -563, -871,-1248},{ -147, -463,  879},
{  -76, 2334, 2840},{  563, 2573, 2385},{  632, 1926, 2920},
{  719, 2023, 1840},{ -545, -723, 1108},{  129, -125,  884},
{ 1417, 1632,  925},{  -94, 1566, 1751},{ -341, 1533, 1551},
{  591,  395, -274},{  -76,  981, 2831},{  153, 2985, 1844},
{ 1032, 2565, 2749},{ 1508, 2832, 1879},{  791, 1199,  538},
{ -190, -453, 1489},{ -278, -548, 1158},{ -245, 1941, 2044},
{ 1024, 1560, 1650},{  512,  253,  466},{  -62, -323, 1151},
{ -473, -376,  507},{ -433, 1380, 2162},{  899, 1943, 1445},
{  134,  704,  440},{  460,  525,  -28},{ -450,  279, 1338},
{    0,  971,  252},{ -445, -627, -991},{ -348, -602,-1424},
{  398,  712, 1656},{ -107,  314, -178},{   93, 2226, 2238},
{  518,  849,  656},{ -462, -711, -447},{  174,  -34, 1191},
{ -119,   42, 1005},{ -372,  274,  758},{ 1036, 2352, 1838},
{  675, 1724, 1498},{  430, 1286, 2133},{ -129, -439,    0},
{ -373,  800, 2144},{    6, 1587, 2478},{  478,  596, 2128},
{ -428, -736, 1505},{  385,  178,  980},{  139,  449, 1225},
{ -526, -842, -982},{  145, 1554, 1242},{  623, 1448,  656},
{  349, 1016, 1482},{   31, -280,  415},{ -316,  724, 1641},
{  360, 1058,  556},{ -436, -358, 1201},{ -355, 1123, 1939},
{  401, 1584, 2248},{ -527,-1012,  355},{  233,  238, 2233},
{ -550, -897, -639},{ -365, -501, 1957},{  389, 1860, 1621},
{  162, 1132, 1264},{ -237, 1174, 1390},{ -640, -411,  116},
{ -228, 1694, 2298},{ 1639, 2186, 2267},{  562, 1273, 2658},
{  323,  338, 1774},{  578, 1107,  852},{   22,  594,  934},
{ -143,  718,  446}
};


static const int16_t lsf_3_2[512][3] = {
{   50,   71,   -9},{ -338, -698,-1407},{  102, -138, -820},
{ -310, -469,-1147},{  414,   67, -267},{ 1060,  814, 1441},
{ 1548, 1360, 1272},{ 1754, 1895, 1661},{ 2019, 2133, 1820},
{ 1808, 2318, 1845},{  644,  -93,  454},{  858,  329, -136},
{  489, -258, -128},{ -198, -745,  -41},{  -52, -265, -985},
{  346,  137,  479},{-1741, -748, -684},{-1163,-1725, -367},
{ -895,-1145, -784},{ -488, -946, -968},{  -85, -390, -725},
{  215, -340, -171},{ 1020,  916, 1969},{  564,  179,  746},
{  662,  977, 1734},{  887,  622,  914},{  939,  856, 1165},
{  309,  688,  803},{  917,  161,  570},{  118,  -20, -283},
{ -816,  -42,  204},{-1228, -325, -462},{ -963, -202, -143},
{ -988, -484, -361},{ -702, -978, -477},{ -302, -790,-1188},
{ -100, -786,-1088},{-1054, -947,-1684},{ -202, -843, -782},
{-1039,-1378, -901},{ -624, -110,  -85},{  356,  213,  -10},
{ -493,  364,  774},{  425,  822,  479},{  -83,  557,  520},
{ -992,-1560, -572},{ -603, -741,  -26},{ -502, -638, -903},
{  209,  306,  147},{ -316, -593, -596},{  -85, -211, -225},
{ -918, -529,  117},{  233, -439, -738},{ 1101,  751,  633},
{ 1457, 1716, 1511},{ 1765, 1457,  910},{ 1122, 1156,  849},
{ 1354,  868,  470},{ -871,-1150,-1796},{ -871, -861, -992},
{ -118,  155,  212},{-1051, -849, -606},{-1117,-1849,-2750},
{-1019,-1427,-1869},{  370, -184, -414},{  959,  493,  104},
{  958, 1039,  543},{  154,  653,  201},{ 1249,  507,  150},
{  663,  503,  230},{  623,  777,  675},{  659,   88, -110},
{  843,  244,  224},{  382,  541,  302},{  724,  433,  666},
{ 1166,  734,  341},{ -138,   20, -397},{-1183, -424,  -46},
{ -321, -352, -124},{ 1333, 1021, 1080},{  262,  366,  723},
{  922,  283, -551},{   31, -636, -611},{ -689, -697, -415},
{ -952, -779, -201},{-1329, -598, -359},{ -953,-1285,  166},
{  493,  305,  221},{  846,  703,  610},{  840,  936,  774},
{ -723,-1324,-1261},{ -357,-1025,-1388},{-1096,-1376, -365},
{-1416,-1881, -608},{-1798,-1727, -674},{ -545,-1173, -703},
{  678,  786,  148},{ -123,  696, 1288},{  644,  350,  -10},
{  414,  614,   15},{  137,  344, -211},{ -814,-1512, -819},
{ -391, -930, -588},{   47, -591, -898},{ -909,-1097, -163},
{-1272,-1167, -157},{-1464,-1525, -389},{-1274,-1188, -624},
{  671,  213,  454},{  124, -274, -525},{ -729, -496, -152},
{-1344,  122,  135},{-2905, -589, -394},{-1728,  441,  -50},
{ 1476,  904,  787},{  316,  236, -440},{ -347,  217,  413},
{ -911, -917,  121},{ -455, -932,  202},{  -92, -465, -375},
{  488,  390,  474},{  876,  729,  316},{-1815,-1312, -669},
{   87,  962,  432},{  563, -249,-1058},{  250,  285, 1105},
{ 1141,  427,  696},{-1038,-1664,-1582},{ -948,  346,  160},
{ -309, -272, -858},{  670,  624, 1250},{ -944, -408, -666},
{ -606, -320, -384},{ -492,  230,   65},{  334,  -50,  -16},
{  -16, -690,-1397},{ 1791, 1716, 1399},{ 2478, 2063, 1404},
{ 1245, 1471, 1426},{ -382,-1037,   -2},{  173, -398, 1145},
{ 1491, 2024, 1801},{  772, 1274, 1506},{ 1429, 1735, 2001},
{ 1079, 1218, 1273},{-1154,-1851,-1329},{ -808,-1133,-1096},
{ -451,-1033,-1722},{   65,  578,  -84},{-1476,-2434,-1778},
{ -765,-1366, -494},{ -218, -594, -931},{  337, -236,  562},
{ 2357, 2662, 1938},{ 1489, 1276,  874},{  189,  358,  374},
{-1519,-2281,-2346},{ -967,-1271,-2095},{ -628,-1188,-1542},
{ 1661, 1043,  546},{  565, 1061,  732},{  -64, -836, -434},
{ -436,  -96,  203},{ 1078, 1216, 1636},{  907, 1534,  986},
{  326,  965,  845},{  142,  -84,  197},{  470, 2379, 1570},
{ 1133,  470, 1214},{  395, 1376, 1200},{ 1125, 1042,  348},
{ -543,-1234, -376},{ -215, -181,  481},{-1947,-1621, -210},
{ -750,-1185,  390},{   29, -399,   27},{  820, 1236,  755},
{  695,  979,  409},{ -174, 1197, 1035},{  912, 1356, 1846},
{ -992,-1437,  484},{-1485,-1700,  208},{ -412, 1204, 1432},
{ -271,  896, 1144},{ -416, 1777, 1434},{-1696,-2644, -204},
{-1789,-1551, 1033},{-1656,-1559, 1303},{-1253,-1589, 1081},
{ -669,-1095,  -66},{ -682,  320, -345},{  659,  305, 1069},
{-1292, -804,  -19},{-1635,-1291,   29},{-1683, -497,   71},
{ -287,   -7, -100},{ -494, -962, -237},{  852, 1881, 1740},
{-1217,-1387,  227},{ -660,  302,  373},{   96, 1087, 1257},
{-1074,-1669,  160},{  485, 2076, 1798},{ -934, -220,  552},
{ -596, -612,  237},{  336, 1720,  879},{  643,  629,  434},
{ 1267,  522, 1633},{   15,  244, -441},{ 1475,  717,  184},
{ 1819, 1590, 1709},{  988,  261,  937},{ 2093, 2345, 1520},
{ 2139, 1858, 1606},{ -577, -579,-1203},{ -956,  135, -488},
{ -464,   51, -338},{ -629, -348, -723},{ 1146, 2073, 1442},
{ 2192, 1466,  911},{-1444,-1572,-2278},{ 1400,  710, 1297},
{ 1335,  633,  928},{ 1434, 2194, 2594},{ 2422, 2204, 1881},
{  982, 2242, 1854},{  380,  792, 1145},{  -63, -539,  414},
{ -252, -964, -314},{-1261, -683, -780},{ -831, -526,-1005},
{-1666,-1135, -424},{-1611, -452, -299},{ 1268, 1048,  642},
{ 1147,  853,  856},{ -675, -336,  139},{ 2268, 1343, 1418},
{   29,  768,  797},{-1224,  423,  564},{-1318,-1082,  245},
{-1302, -812,  573},{-1298,-1617,  646},{ -968,  834,  723},
{  993, 1652, 2027},{ -191, -817,  432},{  662,   60,  198},
{  626,  997, 1330},{ 1648, 1963, 1289},{-1597,  -93,  -45},
{-1088,   37,  -84},{ 1653, 2607, 2337},{ 1065, 2040, 2377},
{ 1139, 2326, 2118},{  859,  357, 1510},{  664, 1227, 1099},
{  479, 1360,  912},{ 1897, 1754, 2019},{ 1168, 1909, 1784},
{  399,   34,  256},{ -593, -304,-1053},{  547, 1694, 1407},
{  647,  -99, -341},{ 1492, 1647, 1190},{   38, -644, -212},
{  395,  846,  222},{ -704, -765, -716},{ -724,-1964,-2804},
{ -150,  291,  -82},{ 1233, 1459, 1007},{ -140, -155,  153},
{  439,  297, 1568},{-1529, -410, -636},{ 1536,  455, -237},
{-1328, -139, -260},{  531,  554,  868},{  269, 1264,  606},
{ -233,  883,  463},{  742,  600, -120},{  -73,  421,  212},
{ -439,  -58,  804},{-1286,-1241,  728},{  294, -490,   50},
{ -591, -905,-1254},{   42, -687,  147},{  -25,  273,  596},
{ -311, 1213,  601},{ -754,  849,  584},{  429,  607,  587},
{ -602, -166,  461},{ -796, -823,  777},{ 1380,  910, 1755},
{  119, 1417,  972},{ -219, -880,-1596},{-1049,-1010,  438},
{ -713,-1379,   78},{    0, -447,-1179},{-1136,-1319,-1573},
{ 2248, 1767, 1309},{  946, 1583, 1432},{ 1150,  482,  436},
{ -469,-1108,  618},{ -447, -966, 1088},{-1252,-1515, -114},
{-1104,-2008, -579},{  210,  613,  497},{-1975,-1437,  642},
{-1269, -856, 1011},{-1646,-1185, 1063},{-1555, -672, 1204},
{-1692,-1114,  623},{ -979,-1326,-1277},{  539, -147,  894},
{-1354, -897, -434},{  888,  475,  428},{  153, -384,  338},
{-1492, -511,  359},{ -974,-1115, -470},{  105, -550,  677},
{ -937,-1145,  877},{  380, -260,  210},{ 1685,  924, 1256},
{ 1775, 1190, 1095},{ 1419,  631,  533},{  627,  299, -347},
{ -411, -534,  647},{ -650,   29, -595},{ -378,-1367, 1563},
{ 1402, 1121, 1465},{ 1089, 1410,  648},{-2096,-1090,   -6},
{  311, -194, -869},{ -639, -831,  416},{-1162,-1224, 1349},
{-1247, -941, 1813},{-2193,-1987,  453},{ -619,-1367, -956},
{-1606,-1972,-1507},{-1175,-1057,-1104},{ -377,  601,  201},
{ 1876,  825,  374},{ -430,-1323,   29},{-1397,-1249,-1331},
{-1007,-1504,  960},{-1401,-2009,  197},{-1379,-1949, -236},
{-1077,  123,  422},{  615, 1269,  546},{ -306, 1526,  904},
{ 1194, 1788, 1177},{ -626, -884,-1526},{  199,  766, 1504},
{-1065,  862,  197},{-1034,-1773, -887},{ -800,  145,  599},
{-1134, -519,  626},{-1205,-1926,  500},{ -910,-1041,-1395},
{-1476,-1567, -969},{ -523,  842,   34},{ 1794,  646,  862},
{-1207,-1888,-1002},{  -78,   -9, -672},{ 1044,  759,   80},
{ -600, 1139, 1019},{   57, 2000, 1422},{ -833, 1414, 1121},
{-1202, 1630, 1260},{ -461, 1420, 1244},{ 1537,  975,  253},
{ -283,  324, -359},{  599, -195,  106},{  588,   62, -587},
{ -757,  645,  205},{   51, 1201,  758},{-1209,  673, -390},
{ -624, 1581,  941},{ -151, 1023,  735},{ 2820, 1301,  690},
{ -302,  524,  -99},{ -900,-1588,-1189},{ 1084,  251,  238},
{ 2014, 1792, 1010},{ 1245, 1633, 1741},{-1227,-1540,-1208},
{ -621,  456, -109},{   40,  -65,  788},{ -805, -699,-1350},
{ -583,  904,  832},{ -801,  532,  594},{ 1972, 1408, 1351},
{-1177,-1880,-2114},{ -773,  568,  948},{-1015, 1079, 1260},
{-1111,  482, -130},{ 1778, 1044,  780},{-1491,  245,  912},
{ -316,-1141, -917},{ -536,-1442,-2346},{ -785,-1546,-1988},
{-2003,  257,  909},{-1849, -633,-1209},{-1538,-1918,-1054},
{ 1606, 2239, 1576},{ -567,-1500,-1544},{-1279,  195, 1369},
{ -817,  293, 1219},{ -525,  630, 1197},{-1698,-2425,-1840},
{ -303,  731,  747},{-1169, -251,  269},{ -950,  -75, 1684},
{-1182, -453, 1005},{-1599,  585,  378},{-2075, -571, -427},
{ -529,-1159,-1171},{ -283, -205, -564},{ -796, 1246,  717},
{ 2277,  927,  539},{ -454,  559,  440},{ -717, 1460, 1615},
{-1030, 1052, 1610},{-1169, -138,  847},{  226,   39, -612},
{-1251, -106, -729},{ -651,  968, 1302},{ -714, -636, 1727},
{  353, 1069,  410},{ -798, -156, 1099},{ -574,  918,  446},
{-1310, 1012,  466},{ 1408, 1591,  765},{ 1429, 1380, 1757},
{ 1949, 1956, 2378},{ 1578, 2047, 2148},{  916,   98,   -7},
{ 1893, 1418, 2141},{  348, 1405, 1579},{  152, 1134, 1801},
{ -267,  154, 1395},{-1166,  469, 1054},{-1142, -405,-1073},
{-1341,-2264,-1581},{ -364,  869, 1706},{-1162,  549, 1550},
{-1225,-1932,-1666},{-1485,-1977,-2055},{-1727, -906,  -98},
{-1897,  233, 1492},{  892,  108, -331},{-1728,-1170,-1700},
{-1060, 1980, 1790},{-1070,-1741,-1909},{  -11, 1539, 1317},
{-1600,   94,  497},{  421,  443, -197},{-1578, -349, -994},
{ -599, -539, 1140},{ -965,-1419, -129},{-1341,  175, -447},
{ -375, 1311, 2055},{ -371, -650, -307},{-1073,  605,  365},
{-2057, -113,  430},{  652,  914,  967},{-1012,-1586,-2323},
{ 1505, 1248,  559},{  262, -486, -401},{-1727, 1342, 1546},
{   50,   56,  432},{ -330,  119, -604},{-1517,-1080, -810},
{  946, 1127, 1055},{-1400,-1703,-1712},{-1270, -704,-1317},
{  807, 1821, 1143},{ 2760, 1606, 2171},{ 1120,  409, -150},
{ -147,  404,  959},{ 2439, 1911, 2189},{ -906, -141, -866},
{ -904, -142, -458},{ -557, -708,-1679},{ -830,-1431,-1583},
{-1842,-1346,-1086},{-1604, -272,  915},{-1196,  772, 1056},
{ -638,-1234,-1897},{ -500,  -81, -822},{-1289,-1613, -735},
{ -117,  785,  168},{-1090, 1133,  922},{-1096, -746, 1384},
{  287, -547,-1063},{-1376,-2201,-1204},{-2176,-1570,-1757},
{-1511,-2241, -771},{-1737, 1099,  830},{-1588,  724, 1243},
{-1542,  693,  805},{-1690, -240, 1665},{-1700,   -4, -668},
{ 2149,  816, 1042},{ -818,-1841,   22},{ -764, -507,  449},
{-1151, -617,  289},{ -843,-1596, -240},{  498, -234, -657},
{ -752,  480, 1678},{ -319, -481,  193},{ -811,  171, -119},
{-2128, -202, -848},{ 1717, 1140, 1700}
};

static const int16_t lsf_3_3[512][4] = {
{   67,  -17,   66,  -12},{-1690, -581, -104, -272},{-1076,-1186,-1845, -376},
{-1140, -926, -420,  -58},{ -259, -656,-1134, -553},{ 1788, 1227,  455,  129},
{  462,  441, -240, -528},{  840,  514,  130,  -75},{ 1114,  623,  153,  216},
{ 1068,  564,   -6, -276},{ 1119,  727,  190,  -68},{  704,  306,  119, -264},
{  329,   61, -100,  156},{  364,  123,  183, -208},{ -171, -123,  220,  -65},
{ -306,  -62,  402,   17},{ -660, -938, -266,    0},{  385,  235,  276,  285},
{  320,  268, -336, -200},{ -724,   17,  -84,  381},{ -544,  429,  494,  519},
{ -117,  288,  304,  329},{  643,  157,  701,  508},{ 1200,  625,  796,  608},
{  998,  421,  492,  632},{ 1204,  780,  446,  132},{ 1257,  844,  547,  449},
{  829,  658,  541,  470},{ 1132, 1258,  918,  639},{  547,   51,  423,  279},
{    9,  392,   83,   94},{  542,  543,  229, -147},{ -198,  129,  194, -185},
{ -863,-1321, -302,   30},{ -597, -629,  -19,  114},{ -900,-1081,  466,  353},
{-1483,-1573,   15, -143},{-1708,-2059, -751,  196},{-1876,-2067, -642, -258},
{-2335,-1470, -450, -564},{ -584, -186, -872, -414},{-1805, -988,-1125,-1310},
{ -726,-1129,   28,  169},{-1039, -864, -718, -246},{  484,   36, -233,  -49},
{  265,   67,  289,  467},{  178,  543,  810,  540},{   84,  282,  672,  703},
{ -975, -777,  129,  287},{ -938, -227,  955,  595},{-1617, -289,  836,  649},
{-1847, -215, 1106,  718},{-2034,-1085,  650,  440},{-2101, -529,  907,  575},
{-2011, -336,  670,  204},{-2389, -692,  360,  137},{-2156,-2204,   -9,  280},
{ -266,  119,   39,  193},{   78,  -59, -120,  226},{ -975, -858, -781,-1095},
{ -619, -413, -451, -842},{-1216,-1321, -813, -883},{-1376,-1615, -394, -428},
{ -737,-1113, -549, -790},{ -880, -975, -967, -642},{ -985, -886,-1273,-1361},
{ -473, -804,-1401,-1407},{  160, -265, -919, -275},{ -248, -250, -718, -380},
{   97, -103, -375, -229},{ -415, -193, -135, -555},{  628,  361,  119,  216},
{  579,  364,  391,  209},{  634,  522, -154, -148},{  526,  389,  170,   33},
{  105,  267,   64,  380},{-1503,-1000,  -30, -369},{-1070,   58,  647,  223},
{-1520, -291,  621,  307},{-1531,  156,  762,  404},{-2029,  141,  734,  499},
{-1849, -650,  306,  512},{ -187, -104,  -59,  438},{  134, -230,  156, -186},
{  -61, -260,  -16,   10},{ -569,   -3, -421, -297},{-1725, -521, -346,  178},
{-1362,  -59,  -44,  157},{-2146, -461, -470, -349},{-2170,   -1, -369, -121},
{-1579, -373, -900,-1015},{-1117, -591, -613, -784},{ -561,  122,  -75, -449},
{   -4, -171, -123, -372},{  192,  168,  -76, -132},{  252, -107,  340,  210},
{  392,  509,  272,  181},{ -109,  145,  218,  119},{ -416, -263,  485,  265},
{ -181,   -8, -286,  226},{ -244, -218,   69, -290},{ -158,  191,   -1,  -64},
{ -592,  -90,  213,  -96},{  255,  435,  178,  -80},{ -369,  -18,  -33,  -80},
{  -42,  415,  140, -222},{ 1143,  651,  649,  329},{  767,  556,  249,  235},
{  948,  413,  442,  279},{  141,  339,  356,  557},{ -470, -170,   99,  237},
{ -569, -800,  352,  565},{  282,  473,  470,  332},{ -199, -690,-1284, -917},
{ -193, -426, -800,-1122},{  -26, -371, -490, -193},{  637,  595,  519,  330},
{  408, -115,   79,   12},{  477,   87, -103, -376},{ -666, -347, -277, -291},
{ -510, -481,  169,  297},{ -829, -738, -205, -171},{ -320, -540,  328,  283},
{ -859, -958,  442,   -2},{  556,  686,  130,   56},{ 1383, 1012,  755,  427},
{  612,  741,  628,  553},{ -339, -796,  134,  277},{ -633,-1085,   -2, -246},
{ -880,-1035,-1607,-1064},{ -994, -474,-1138, -488},{ -414, -795,   73, -206},
{   -8, -139,  439,  204},{ -176, -578,   23,  131},{ -269, -757, -191,  245},
{ -109, -338,  112,  316},{  120, -406, -118,  611},{ -180, -186, -645,  115},
{ -173,   34, -518, -489},{ -151,   61, -583, -844},{  220, -138, -681,-1020},
{  391,  -17, -598, -321},{  157, -295,  129,  155},{ -926, -875, -987,  285},
{  241,  -83, -125, -125},{  620,  597,  432,   92},{  393,   78,  409,   61},
{ -393, -739, -413, -748},{   83,   54,  361,   27},{-1084,  130, -337, -694},
{-1565,  297,  318,  -19},{-1873,   36,   51, -317},{-2323, -246,  231,  -84},
{-2306, -783,   40, -179},{-2233, -930, -474, -462},{ -754,  -86, -288, -626},
{-2411, -455,  -63,  171},{-1099,-1094,  -26, -143},{-1193, -455, -406, -381},
{ -605, -210,  -96,  -51},{ -580, -476, -276,  -15},{-1195, -634,-1203, -881},
{ -378, -221, -669, -952},{  594,  178, -403, -676},{  763,  327,  601,  290},
{  172,  300,  203,  157},{  -56, -336,  356,   24},{ -228, -296, -259,  -29},
{ -186,  263,  416,   14},{ -353,  373,  -12, -216},{  257,   96,  174,   57},
{-1526, -616, -954, -499},{ -497, -152, -333,  125},{  105,  200,  179,  -97},
{ -331, -224,  765,  697},{  760,  256,  301,   59},{  455,  -85,  204,  288},
{ -514,  240,  251, -109},{  256,  417,  -34, -413},{  101,  430,  384,  156},
{  -31,  -10,  206,  426},{  589,  145,  143,   71},{  808,  906,  333,  349},
{  986,  938,  589,  331},{ 1300,  824,  187,  509},{ 1062,  653,  379,  466},
{ 1462,  937,  401,  274},{  787,  861,  265,    2},{  609,  553,   28,  305},
{  926,  340,  106,  386},{  241, -267, -147,  225},{ -178, -534,  347,  502},
{ -643, -381,  397,   30},{ -651, -733, -435,  398},{ -407, -726, -484, -248},
{ -789, -914, -438, -476},{ -498, -390,   75, -295},{ -964, -590, -606,  150},
{ -121,  -49, -155,  -78},{  935,  550,  389,   38},{ -321,  127,  424,  315},
{ -285, -113,  283,  259},{  658,  203,  322,  486},{  903,  505,  748,  417},
{  611,  423,  555,  512},{  239,  -83, -578,  -19},{ -339, -731,  349,   13},
{ -934,-1399, -114, -360},{  107,  692,  182,   90},{-1243,-1538,-1551, -725},
{ -568, -903,-1363, -525},{ -517, -853, -861,-1004},{ -168, -690, -835,   63},
{ -137, -556, -547,  144},{ -286, -817,  485,  319},{ -147, -408,  526,  246},
{ -347, -434,  297,  -28},{ -290, -471,-1110,-1285},{ -460, -359, -988, -794},
{ 1347, 1299,  690,  523},{ 1216, 1068, 1094,  757},{  825, 1140,  752,  494},
{ 1252, 1365, 1195,  898},{  521, 1053,  532,  432},{ -334, -216, -313, -263},
{ -160,   52, -472, -155},{  127,  136, -380,   44},{  851,  410, -162, -489},
{  123, -255, -796, -667},{ 1090,  917,  789,  493},{ 1397, 1197,  558,  202},
{  -51, -118, -342, -701},{   83,  108,  -42, -441},{   61,   95,  287,  256},
{  -27,   89,  524,  531},{  351,  227,  592,  545},{  697,  155, -164,  307},
{  638,  274, -489,  -50},{  754,  240, -166, -124},{ -116, -579,-1212,  -63},
{  190, -295,-1040,-1296},{  147, -376, -177, -113},{  841, 1241, 1051,  668},
{    2,  293,  551,  304},{-1096, -953, -248,  376},{ -750, -965,   87,  516},
{ -275, -516,  689,  391},{ -379, -643,  876,  594},{ -390,-1013, -645,  573},
{ -107, -568, -689, -826},{-1025,  -27, -328, -203},{  861,  749,  548,  233},
{-1660,-1043,  451,  108},{ -660, -620,  430,  236},{   21, -396,-1158, -631},
{ 1372, 1298,  967,  577},{ 1125, 1125,  589,  454},{ -323, -865, -467,  153},
{ -468, -699, -804, -509},{ -392, -718, -204,  -35},{ -603,-1093, -567, -162},
{ -505,-1004, -102,  350},{  219,  224,  423,  252},{  395,  591,  608,  363},
{ -746,  -96,  373,  172},{  171,  295,  714,  339},{  233,   77,  107,  277},
{  157,  153, -499, -356},{ 1547, 1073,  576,  494},{ -292, -339, -504, -592},
{ -903,  -72, -619, -481},{-1594,-1117, -567, -254},{ -793, -507, -564, -291},
{ -492, -532,  502,  560},{ -382,  427,  600,  230},{ -227,  477,  251,   75},
{  285,  842,  813,  476},{-1310,-1333,  186,  377},{ -587, -917,  643,  381},
{-1186, -553,  411,   82},{-1127, -820, -174, -540},{ -604,  119,  543,  205},
{ -380,  657,  909,  567},{  112, -298, -374,  114},{ -857, -251,   56,  159},
{  401,  345,  -34, -140},{ -111, -607,   41,  614},{  355, -114,  -77,  474},
{  578,   56, 1450,  924},{ 1098, 1420,  741,  400},{  246,   22,  588,  313},
{ -121,  327,  831,  472},{-1138, -608,  856,  552},{-1241,-1072,  638,  600},
{ -358,  254, -333, -303},{ -646,  739,  358,   74},{ 1226, 1671, 1221,  849},
{ 2241, 1624,  983,  636},{ 1841, 1477,  749,  384},{  350,  263,   87,  128},
{-1902, -941, -144,  -64},{-1734, -255,  288,  -31},{-2644,-1238,  366,  235},
{-1643,-1092,-1344, -304},{ -541,-1075,-1116,  123},{-1178, -252, -816, -180},
{-1016,  533,  565,  233},{ -487, -430, -188,  334},{  867, 1236,  534,  171},
{-1590,-1607,  635,  630},{-2196,  310,  924,  412},{-2358, -328,  956,  529},
{-2639, -377,  630,  278},{-2602,  317,  799,  299},{-2406,  133,  340,   31},
{-2156,-1468,  131,  125},{-1184, -490, -139,   46},{ -744,  447,  891,  564},
{   67, -451,  646,  604},{ -553, -429, -876,  396},{  162,  -66, 1305,  915},
{  479,  579, 1088,  794},{  450,  278,  566,  324},{-1057, -154,  148, -177},
{-2545,  168, 1070,  592},{-2351,  -42,  819,  345},{-2344, -707,  721,  250},
{-2175,-1497, -309,  122},{  -78,  -73,  120,  173},{   -4,  262, -263, -261},
{ -431,  -64, -405, -732},{-2609,  116,  -83, -193},{-1525, -944, -477, -725},
{ -508,  307,  170,  172},{  832,  417,  832,  686},{ -225,  177,  894,  818},
{ -482, -389, 1279, 1039},{ -383,  201, -350,   40},{  730,  635,  226,  526},
{  503,  462,  338,  398},{  535,  714,   40, -282},{ 1482, 1471, 1085,  731},
{ 1561, 1072,  909,  693},{ 1419, 1282,  889,  879},{ 1153,  728, 1186,  840},
{ -226, 1130,  949,  689},{ -494, -986,-1556, -128},{ -568, -721, -713,  -26},
{  317,  524,   70,  135},{ -405, -865,-1766, -652},{ -174, -801,  885,  773},
{ -153,  -91, 1099,  751},{ -506,-1149,  853,  646},{  241,  782,  519,  539},
{ 1853, 1700, 1101,  684},{-1249,-1486, -464,  188},{ -893,-1409,-1312, -341},
{ -135,  438, -175,   18},{ 1111,  976,  319,  208},{-1430,-1768,   83,  458},
{ -530,-1000,  307,  129},{ -840,  -15,  -29, -356},{ -911, -924,-1147, -242},
{ -119, -528,  127, -133},{ -761, -765,  190,  -83},{ -315,  895,  522,  231},
{ -222,  102,  -63, -428},{  316,  699,  379,   70},{   25,  716,  314, -108},
{  507,  874,  566,  238},{  108,  941,  519,  195},{  425,  -60, -427,  257},
{  139, -103, -630,  446},{  334,  370,  412,   48},{ -172, -690, -283,  557},
{  187, -286,  158,  483},{  140,  270, -344, -631},{  924,  579, -116,  132},
{  142,  466,  -68,  -64},{  230, -145, -302, -542},{ -803, -912, 1018,  737},
{ -773, 1015,  630,  297},{-2596,   95,  445,  336},{-2122,  491,  510,  191},
{-1253,  161,   -2, -324},{-1450, -633, -712, -105},{ -842, -254, -411,  100},
{ -640, -290, 1010,  763},{ -650,  313, 1169,  730},{  140,  505, 1030,  766},
{  772,  287, 1067,  823},{  495,  749,  305,  323},{ -164,  462,   78,  399},
{ -342, -874,   69,  597},{  -16,  620,  621,  337},{ -138, -444, -265,  218},
{   84, -450,  953,  666},{ -222, -803,  541,  604},{ -921,-1376,  244,  116},
{ -841, -723,  630,  588},{  140,  663,  294,  368},{  935, 1046,  881,  759},
{ 1746, 1464,  916,  628},{  436,  963,  281,    1},{ -119,   74,  542,  213},
{    1, -567,  301,  241},{  260,  435,  222,  396},{  936,  957, 1108,  703},
{  510,  506,  808,  478},{  601,  694,  960,  620},{  972,  741,  980,  600},
{  834,  717,  767,  684},{  643,  972,  935,  638},{  501,  661,  720,  851},
{ -105, -632, -303, -117},{ -429,  130,  789,  442},{ -522, -188,  704,  373},
{ -759,   42,  814,  523},{ -531,-1137,  373,  578},{ -682,-1203, -455,  285},
{-1163,-1577,-1098,   44},{   81,  -82,  712,  363},{  477,  246,  954,  622},
{ 1604, 1622, 1277,  891},{ 1409,  859,  924,  892},{  774, 1041,  947, 1142},
{   40, -546,  -75,  288},{ -616, -106, -697,  -26},{ -169, -160, -891, -739},
{ -279, -384,-1029, -350},{ 1781, 1308, 1046,  816},{ 1580, 1533, 1472, 1178},
{ 1505, 1076, 1216,  899},{  890,  904,  564,  654},{  920,  692, 1021,  856},
{ -493,  132,  177,  505},{   71,  195,  -28,   97},{  456,  351, -164,   88},
{  439,  278,  -40,  350},{ 1395,  949,  234,  -95},{ -805, -472,   38, -163},
{  367,  -98,  489,  523},{ 1025, 1178, 1212,  906},{  319, 1314,  814,  461},
{ -123, -543, -804,  447},{ -748, -324, -897,-1127},{ -737, -501, -789, -713},
{  715,  777, 1239,  922},{ 1949, 1939, 1368,  865},{  730,  880,  758,  388},
{ -871,  454,   17, -251},{ -381, -810,-1583,  239},{ -521, -966, -792,  259},
{ -890,-1358, -770,  -73},{  166,  349, -212,  323},{ -840, -301,  473,  435},
{ -679, -464,  728,  351},{ -156, -199,  667,  432},{   29, -252,  415,  480},
{ -731, -379,  145,  559},{ -528, -631,-1158, -159},{  445,  273,  123,  639},
{  373, -126,  800,  568},{   84, -162,  720,  712},{ -830, -536, -185,  222},
{  408,  452,  501,  771},{ -897,-1355,  -67,  442},{ -792,-1406,  566,  602},
{  167, -326,  509,  330},{  -95, -626, -730, -344},{ 1668, 1217,  779,  455},
{ 1316,  828,  584,  719},{  404,  -31, 1013,  789},{   89,  107,  891,  549},
{  871, 1581,  917,  671},{  866, 1479, 1289,  854},{  391, 1068, 1122,  812},
{   78, -562,  345,  563},{  429, -103,  417,  787},{ -122, -437,  411,  788},
{ -913, -417,  602,  754},{ -226,  -16,  151,  760},{ -700,  118, -104,  -14},
{-1128,   48,  284,  393},{ -390, -419, -639, -116},{ -910,  306,  316,  -13},
{ 1207,  984,  821,  669},{-1195, -693,  140, -213},{ -884, -416, -199, -558},
{ -616,  245, -404, -664},{  262,   56, -617, -724},{  -85, -491, -320, -656},
{ -570, -831, -129, -528},{-1506,  -63, -367, -385},{ -358, -321,    4,   51},
{ -366, -214,  319,  511},{  146,  671,  -17, -291},{ -110,  464, -139, -496},
{ -202,  220, -312, -631},{ -660,  -73, -655, -820},{ -662, -653,-1288, -857},
{ -430, -953, -959, -264},{  -49, -468,  -72, -381},{ -350, -563, -193, -407},
{   55, -408, -803,   11},{ -309,  649,  188, -198},{ -512,  461,  -79, -458},
{-1318, -263, -134, -523},{-1657, -435, -495, -765},{   57, -347, -414,  434},
{-1141, -242, -664, -857},{   34,  -68, -707, -338}
};

static const int16_t lsf_5_1[128][4] = {
{ -451,-1065, -529,-1305},{ -450, -756, -497, -863},{ -384, -619, -413, -669},
{ -317, -538, -331, -556},{ -414, -508, -424, -378},{ -274, -324, -434, -614},
{ -226, -500, -232, -514},{ -263, -377, -298, -410},{ -151, -710, -174, -818},
{ -149, -412, -156, -429},{ -288, -462, -186, -203},{ -170, -302, -191, -321},
{ -131, -147, -297, -395},{ -228, -214, -245, -192},{  -67, -316,  -71, -327},
{ -104, -205,  -94, -183},{ -143,  -38, -193,  -95},{   16,  -76, -124, -248},
{   23, -237,   24, -244},{   18, -136,   44, -111},{  -33,  -24,  -25,    0},
{  149,   19,   23, -143},{  158, -169,  174, -181},{  133,  -55,  165,  -26},
{  111,   84,   98,   75},{   87,  183, -115,  -11},{   -8,  130,   11,  170},
{  254,   77,  205,   17},{  183,  112,  262,  194},{  202,  287,   95,  189},
{  -42, -105,  234,  179},{   39,  186,  163,  345},{  332,  199,  299,  161},
{  -54,  285,  -78,  281},{ -133,  141, -182,  111},{  249,  341,  271,  364},
{   93,  403,   75,  391},{   92,  510, -138,  220},{ -185,  -29,  -34,  361},
{ -115,  320,    3,  554},{   99,  286,  218,  591},{ -245,  406, -268,  453},
{    0,  580,   25,  606},{  275,  532,  148,  450},{  -73,  739, -285,  518},
{ -288,   94, -203,  674},{ -140,  -74,  205,  714},{ -114,  299,  176,  923},
{  182,  557,  240,  705},{  -16,  513,  485,  593},{  293,  384,  451,  617},
{  -38,   50,  563,  529},{  303,  209,  459,  363},{  433,  452,  450,  454},
{  367,  606,  477,  741},{  432,  353,  368,  267},{  361,  716,  273,  583},
{  453,  166,  510,  172},{  201,  629,  274,  191},{  568,  639,  302,  298},
{  634,  387,  643,  350},{  587,  560,  612,  565},{  600,  788,  487,  672},
{  512, 1015,  321,  333},{  357,  854, -125,  413},{  474,  712,   17, -151},
{  564,  285,  270, -241},{  971,  889,  489,  220},{  510,  896,  549,  924},
{  327,  825,  290,  911},{  540, 1108,  158,  805},{  199,  957,  511,  730},
{  100,  874,   13,  791},{  435,  632,  676,  972},{  249,  900,  467, 1218},
{  781, 1074,  585,  785},{  -23,  669,  267, 1043},{  619, 1084,  615, 1145},
{  622,  905,  916, 1049},{   80,  331,  584, 1075},{   89,  639,  988,  961},
{  770,  720,  798,  699},{  492,  447,  899,  627},{  271, 1188,  725, 1333},
{   87,  603,  832, 1603},{  616, 1127,  890, 1505},{ 1000, 1156,  866, 1009},
{  995,  827, 1149,  858},{  817, 1450,  773, 1320},{  500, 1389,  312, 1153},
{  -20, 1084,   64, 1283},{    2, 1172,  399, 1869},{  514, 1706,  502, 1636},
{  886, 1522,  416,  600},{ 1131, 1350, 1275, 1390},{  889, 1795,  914, 1766},
{  227, 1183, 1250, 1826},{  505, 1854,  919, 2353},{ -199,  431,  152, 1735},
{ -213,  -28,  392, 1334},{ -153,  -52,  978, 1151},{ -323, -400,  813, 1703},
{ -136,   84, 1449, 2015},{ -331, -143, -137, 1192},{ -256,  534, -157, 1031},
{ -307, -439,  542,  731},{ -329, -420,  -97,  616},{ -362, -168, -322,  366},
{ -247, -110, -211,   89},{ -196, -309,   20,   59},{ -364, -463, -286,   89},
{ -336,  175, -432,  141},{ -379, -190, -434, -196},{  -79,  150, -278, -227},
{ -280,  166, -555, -422},{ -155,  541, -366,   54},{  -29,  -83, -301, -774},
{  186,  628, -397, -264},{  242,  293, -197, -585},{  124,  410,   53, -133},
{   10,  340, -570,-1065},{   65, -446,   68, -493},{  383,  937, -357, -711},
{ -359, -250, -677,-1068},{  292,  -26,  363,    6},{  607, 1313, -127,  -10},
{ 1513, 1886,  713,  972},{ 1469, 2181, 1443, 2016}
};

static const int16_t lsf_5_2[256][4] = {
{-1631,-1600,-1796,-2290},{-1027,-1770,-1100,-2025},{-1277,-1388,-1367,-1534},
{ -947,-1461, -972,-1524},{ -999,-1222,-1020,-1172},{ -815, -987, -992,-1371},
{-1216,-1006,-1289,-1094},{ -744,-1268, -755,-1293},{ -862, -923, -905, -984},
{ -678,-1051, -685,-1050},{-1087, -985,-1062, -679},{ -989, -641,-1127, -976},
{ -762, -654, -890, -806},{ -833,-1091, -706, -629},{ -621, -806, -640, -812},
{ -775, -634, -779, -543},{ -996, -565,-1075, -580},{ -546, -611, -572, -619},
{ -760, -290, -879, -526},{ -823, -462, -795, -253},{ -553, -415, -589, -439},
{ -533, -340, -692, -935},{ -505, -772, -702,-1131},{ -263, -306, -971, -483},
{ -445,  -74, -555, -548},{ -614, -129, -693, -234},{ -396, -246, -475, -250},
{ -265, -404, -376, -514},{ -417, -510, -300, -313},{ -334, -664, -463, -814},
{ -386, -704, -337, -615},{ -234, -201, -233, -239},{ -167, -567, -203, -619},
{ -147, -415, -115, -352},{ -166, -750, -171, -761},{ -270, -879, -264, -903},
{ -367, -744,   43, -475},{   14, -653,   43, -670},{   11, -448,  -59, -521},
{ -126, -119, -155, -613},{  -42, -863,  -27, -931},{  136, -483,  183, -468},
{   55, -298,   55, -304},{  313, -609,  313, -720},{  322, -167,  100, -541},
{   -3, -119, -111, -187},{  233, -236,  260, -234},{   26, -165,  134,  -45},
{  -40, -549,  360, -203},{  378, -388,  450, -383},{  275,   20,  182, -103},
{  246, -111,  431,   37},{  462, -146,  487, -157},{ -284,  -59,  503, -184},
{   24,   53,   -3,   54},{  122,  259,  333,   66},{  484,  104,  436,   68},
{  195,  116,  190,  206},{  269,   -9,  482,  352},{  382,  285,  399,  277},
{  452,  256,   69,  186},{   13,  297,  -13,  259},{  -95,   30,   56,  394},
{  196,  425,  205,  456},{  281,  577,   15,  191},{  375,  290,  407,  576},
{  -56,  227,  544,  405},{    0,  549,  -92,  528},{ -229,  351, -245,  338},
{ -362,  435,  167,  527},{  -75,  302,   91,  824},{  129,  599,  496,  679},
{  186,  749,  153,  737},{ -281,  600, -348,  615},{ -236,  769,   41,  881},
{   38,  890, -220,  841},{ -357,  883, -393,  903},{ -634,  474, -444,  850},
{ -175,  678, -493,  242},{ -519,  785, -714,  582},{ -541,  366, -543,  434},
{ -597,  500, -765,  222},{ -702,  917, -743,  962},{ -869,  501, -899,  548},
{ -379,  200, -435,  157},{ -819,  214, -861,  157},{ -614,   40, -632,   94},
{ -883,  -54, -741,  516},{ -501,  298, -614, -171},{ -870, -161, -865,  -23},
{ -818,   93,-1015, -267},{ -662, -359, -549,    2},{ -442, -121, -377,    0},
{ -227,   33, -414, -126},{ -129,  212, -934,   34},{-1082, -282,-1119, -268},
{ -710, -825, -420, -191},{-1076, -928, -917,  -93},{ -628, -358,   97,    7},
{ -206, -393, -101,   24},{ -203,   38, -168,   83},{ -599, -423, -279,  426},
{ -700,  118,  -75,  206},{ -981, -673, -680,  417},{ -367,   37, -279,  474},
{ -129, -318,  319,  296},{ -626,  -39,  343,  602},{ -696,  -39, -303,  940},
{  104,  233, -380,  137},{  -36,  269,  -75, -214},{  120,   43, -529, -477},
{  459,  164, -202, -229},{  -49, -167,  609,  792},{   98, -220,  915,  148},
{  293,  283,  869,   91},{  575,  394,  326,  -78},{  717,   67,  365, -323},
{  616,  -36,  731,   27},{  619,  238,  632,  273},{  448,   99,  801,  476},
{  869,  273,  685,   64},{  789,   72, 1021,  217},{  793,  459,  734,  360},
{  646,  480,  360,  322},{  429,  464,  638,  430},{  756,  363, 1000,  404},
{  683,  528,  602,  615},{  655,  413,  946,  687},{  937,  602,  904,  604},
{  555,  737,  786,  662},{  467,  654,  362,  589},{  929,  710,  498,  478},
{  415,  420,  693,  883},{  813,  683,  781,  925},{  913,  939,  726,  732},
{  491,  853,  531,  948},{  734,  963,  315,  808},{  761,  755, 1144,  760},
{  655, 1076,  826, 1057},{ 1091,  838, 1003,  808},{ 1047, 1133,  659, 1101},
{  992, 1050, 1074, 1075},{  971,  694, 1226, 1054},{  571,  841,  884, 1404},
{ 1379, 1096, 1080,  861},{ 1231,  735, 1284,  760},{ 1272,  991, 1367, 1053},
{ 1257,  700, 1050,  534},{  988,  453, 1264,  599},{ 1140,  679, 1621,  815},
{ 1384,  521, 1317,  393},{ 1564,  805, 1448,  686},{ 1068,  648,  875,  307},
{ 1083,  361, 1047,  317},{ 1417,  964,  675,  571},{ 1152,   79, 1114,  -47},
{ 1530,  311, 1721,  314},{ 1166,  689,  514,  -94},{  349,  282, 1412,  328},
{ 1025,  487,  -65,   57},{  805,  970,   36,   62},{  769, -263,  791, -346},
{  637,  699, -137,  620},{  534,  541, -735,  194},{  711,  300, -268, -863},
{  926,  769, -708, -428},{  506,  174, -892, -630},{  435,  547,-1435, -258},
{  621,  471,-1018,-1368},{ -393,  521, -920, -686},{  -25,   20, -982,-1156},
{  340,    9,-1558,-1135},{ -352,   48,-1579, -402},{ -887,    6,-1156, -888},
{ -548, -352,-1643,-1168},{ -159,  610,-2024, -963},{ -225,  193,-1656,-1960},
{ -245, -493, -964,-1680},{ -936, -635,-1299,-1744},{-1388, -604,-1540, -835},
{-1397, -135,-1588, -290},{-1670, -712,-2011,-1632},{-1663,  -27,-2258, -811},
{-1157,  184,-1265,  189},{-1367,  586,-2011,  201},{ -790,  712,-1210,    3},
{-1033,  808,-1251,  830},{ -111,  635,-1636,  447},{ -463, -949, -445, -928},
{ -504,-1162, -501,-1211},{  144, -351, -372,-1052},{ -283,-1059, -279,-1123},
{ -575,-1438, -587,-1614},{ -935, -984,  229,  690},{ -921, -719, -403, 1362},
{ -685, -465,  874,  397},{ -509,  -46,  317, 1334},{ -485,  456,  813,  439},
{ -411,  339,  898, 1067},{ -425,   46, 1441,  497},{ -909, -800, 1465, 1046},
{ -254, -321, 1430, 1165},{   68,  350, 1034,  666},{  370,   11, 1311,  790},
{  143,  232, 1041, 1562},{ -114,  663, 1616, 1078},{  454,  579, 1275, 1040},
{  -76,  909,  752, 1067},{  153,  512,  348, 1214},{  614,  385, 1843,  808},
{  269, 1034,  203, 1086},{  652, 1017, 1783, 1130},{  429, 1327,  387, 1384},
{  -49, 1183,  -72, 1215},{ -416, 1001,  544, 1749},{ -352, 1223, -502, 1199},
{ -589,  569, -227, 1630},{ -142, 1578, -230, 1715},{ -714, 1288, -838, 1398},
{ 1131, 1357, -208, 1232},{  437,  965, -929,  818},{  811, 1410,  859, 1507},
{  164, 1212, 1387, 1793},{  484, 1874,  456, 2063},{  996, 1170, 1326, 1402},
{ 1316, 1360, 1135, 1262},{ 1234, 1618, 1361, 1768},{ 1421, 1227, 1584, 1347},
{  854,  672, 1685, 1566},{ 1139, 1270, 2016, 1825},{ 1773, 1581, 1532, 1460},
{ 1487,  946, 1659, 1021},{ 1744, 1212, 1392,  977},{ 1772, 1161, 1826, 1164},
{ 1718, 1429, 1973, 1591},{ 1185,  864, 2132, 1061},{ 1799,  814, 1838,  757},
{ 2104, 1315, 2054, 1258},{ 2113,  915, 2331,  930},{ 1467, 1147, 2590, 1439},
{ 2245, 1744, 2090, 1620},{ 2358, 1454, 2666, 1506},{ 1876, 1837, 2070, 1975},
{ 1739, 1577,  682, 1289},{ 1584, 2045, 1454, 2098},{ 2498, 2004, 2711, 2066},
{  726, 1588, 2756, 2336},{  228,  847, 2456, 1659},{   36,  301, 1942, 1957},
{ -446,  -96, 2154, 1396},{ 1533, 1101,   14,  608},{ -923, -732, 1383, 1982},
{ 1345,  952, -680,  321},{ 1281, 1268,-1594,  365},{  941,  946,-1737, -822},
{ 2374, 2787, 1821, 2788}
};

static const int16_t lsf_5_3[256][4] = {
{-1812,-2275,-1879,-2537},{-1640,-1848,-1695,-2004},{-1220,-1912,-1221,-2106},
{-1559,-1588,-1573,-1556},{-1195,-1615,-1224,-1727},{-1359,-1151,-1616,-1948},
{-1274,-1391,-1305,-1403},{-1607,-1179,-1676,-1311},{-1443,-1478,-1367, -898},
{-1256,-1059,-1331,-1134},{ -982,-1133,-1149,-1504},{-1080,-1308,-1020,-1183},
{ -980,-1486, -967,-1495},{ -988, -922,-1047,-1077},{ -838,-1179, -858,-1222},
{-1131,-1041,-1064, -767},{ -872,-1157, -701, -880},{ -706, -906, -774,-1016},
{ -578,-1080, -801,-1478},{ -591,-1111, -592,-1146},{ -713,-1388, -640,-1376},
{ -597,-1059, -416, -903},{ -686, -832, -661, -708},{ -444, -868, -490, -921},
{ -374, -776, -619,-1170},{ -585, -549, -769, -795},{ -435, -659, -530, -741},
{ -498, -837, -357, -597},{ -279, -871, -243, -887},{ -282, -665, -280, -667},
{ -165, -560, -394, -903},{ -362, -410, -448, -583},{ -409, -574, -313, -357},
{ -637, -548, -570, -436},{ -896, -504, -382, -757},{  -58, -481, -165, -618},
{ -191, -374, -234, -382},{ -222, -683,  -25, -480},{ -418, -359, -730, -353},
{ -324, -157, -432, -322},{ -394, -303, -284, -104},{ -601, -289, -556, -196},
{ -588, -150, -659, -608},{ -473,  -24,  -68, -448},{ -474,   -8, -506,  -45},
{ -748, -184, -844, -252},{ -901,  -91, -584,  -97},{ -652,  138, -764, -131},
{ -678,  -12, -670,  165},{ -259,   -3, -840, -107},{ -909,   37, -992,   44},
{ -854, -415, -839,   13},{-1001, -271,-1026, -309},{ -798, -478, -832, -488},
{ -943,  168,-1112, -387},{-1185, -101,-1183,  -40},{ -941, -316,-1030, -770},
{-1044, -625,-1081, -538},{-1224, -299,-1312, -436},{-1197, -663,-1167, -161},
{-1216, -690,-1237, -831},{-1432, -720,-1403, -493},{ -898, -740, -922, -801},
{-1102, -402,-1579, -964},{-1061, -638,-1269,-1438},{-1499, -934,-1502, -895},
{-1598, -564,-1723, -717},{ -606, -597,-1166,-1085},{-1369, -468,-1946,-1493},
{-1838, -953,-1932, -931},{-1499, -188,-1635, -421},{-1457, -338,-1448,  -22},
{-1942, -422,-2006, -249},{ -496, -114,-1910, -755},{-1289,  174,-1451, -109},
{ -482, -257,-1221, -508},{-1617,  151,-1694,  208},{ -654,  107,-1651,   29},
{-1141,  279,-1215,  306},{-1228, -506, -730, -175},{-1236, -101, -969,  551},
{ -870,  278, -823,  315},{ -563,  376,-1051,  228},{ -507,  280, -599,  281},
{ -758,  253, -305,  379},{ -755, -134, -611,  660},{ -824,  536, -817,  646},
{ -413,   49, -341,  177},{ -453,  526, -482,  589},{  -71,  339, -657,  264},
{ -244,  295, -237,  315},{ -387,  569, -506,   -9},{ -377,   14, -160,  661},
{ -216,   40, -308,  -46},{   95,  214, -242,  167},{  -86,  192,  -56,   27},
{  -76,   31,   36,  309},{ -106, -182, -113,   74},{ -441,  -22,   23,  139},
{   81,  -11,   44,   15},{  -87, -137, -118, -207},{ -158,  -58,  272,  -92},
{ -156, -441,    8, -136},{  128, -221,  101, -218},{   40, -197,  -76, -456},
{    9, -445,   33, -423},{  226,   60,   73, -222},{  156, -399,  280, -318},
{  245, -341,  166, -499},{  339, -190,  327, -219},{  325, -137,  -89, -596},
{  100, -627,  144, -677},{  487,   28,  252, -391},{  214,  -41,  282,  -28},
{   99, -286,  331,   49},{  459, -388,  565, -369},{  436,   28,  336,   -9},
{  397, -167,  618,   34},{  596,  -17,  561, -140},{  299,   79,  522,  125},
{  203,    2,  244,  288},{  255,  211,  175,   82},{  596,  187,  517,  108},
{  381,  255,  365,  297},{  497,  352,  327,  -82},{   25,  210,  371,  245},
{  261,    3,  545,  449},{  140,  294,   44,  295},{  212,  347,  244,  494},
{  331,  528,  201,  307},{  349,  411,  613,  284},{  614,  413,  464,  322},
{  624,  397,   97,  200},{ -160,  384,  149,  362},{  495,  525,  269,  585},
{   33,  491, -121,  433},{  427,  611,  498,  516},{  171,  443,  497,  666},
{  440,  275,  566,  575},{  146,  639,  155,  670},{  -33,  173,  212,  696},
{ -166,  601, -191,  695},{ -489,  503,  175,  742},{  214,  476,  372, 1083},
{  578,  530,  586,  777},{  425,  874,  315,  841},{  374,  848, -165,  565},
{   35,  991,  -39, 1062},{  329,  712,  786,  840},{  645,  795,  661,  676},
{  571,  918,  632, 1079},{  673,  817,  318,  388},{  874, 1012,  564,  848},
{  880,  620,  557,  479},{  671,  453,  692,  468},{  840,  642,  844,  645},
{  506,  428,  897,  567},{  837,  387,  962,  499},{  691,  561,  939,  926},
{  783,  296,  790,  268},{ 1028,  530,  874,  329},{  548,  143,  675,  291},
{  503,   66, 1041,  359},{  786,   97,  805,   33},{  837,  470,  511,   49},
{ 1092,  327, 1174,  323},{    3,  242,  872,  474},{  689,  429, 1329,  678},
{ 1042,  620, 1109,  664},{  321,  193,  889,  950},{ 1153,  874,  893,  635},
{  877,  862,  948,  913},{ 1293,  665, 1320,  639},{  997,  793, 1402, 1030},
{ 1176, 1012, 1110,  959},{ 1410,  925, 1403,  915},{  543,  862, 1116, 1222},
{  835, 1190,  835, 1190},{  959, 1148, 1147, 1376},{ 1300, 1193, 1415, 1231},
{ 1335, 1341,  746, 1092},{ 1711, 1283, 1389, 1073},{ 1334, 1566, 1153, 1475},
{ 1645, 1137, 1825, 1220},{ 1056, 1382, 1521, 1730},{ 1632, 1545, 1620, 1542},
{  855, 1596,  865, 1667},{  693,  885, 1716, 1519},{ 1167, 1296, 2209, 1760},
{ 1952, 1493, 2020, 1482},{ 1534, 1866, 1694, 2008},{ 1566,  748, 1761,  825},
{  294, 1392, 1084, 2058},{  621, 1315,  365, 1287},{  198, 1028,  488, 1408},
{  249,  403, 1014, 1561},{  324,  363, 1645, 1044},{  193,  367, 2034, 1859},
{ -251,  579,  750,  994},{ -243,   30, 1325,  879},{  -28, -169,  624,  917},
{ -453,  159,  186, 1370},{ -614,    6,  537,  392},{  -94, -291,  781,  229},
{ -128, -298,  245,  491},{ -701, -648,  972,  789},{ -501, -640,  178,  255},
{ -365, -390, -255,  317},{ -958, -294, -191,  228},{ -775, -447,  157, -237},
{ -657, -720, -407,   92},{ -117, -611,  334, -230},{ -679,-1084, -144, -317},
{ -901, -861, -738, -360},{  -85, -727,  -90, -787},{  100,  -22, -391, -263},
{  -56,  -73, -337, -754},{    5, -189, -706, -624},{   89, -344, -135,-1113},
{ -353, -237, -684,-1135},{ -275,-1102, -269,-1203},{  152,  145, -722,-1232},
{   49,   80,-1248, -776},{ -248,  391, -732, -547},{  469,  218, -255, -864},
{   69,  366, -166, -485},{ -688,  191,-1212,-1196},{ -170, -169,-1308,-1631},
{  321,  470,-1419,-1243},{  -64,  272,-1361, -248},{  492,  565, -721, -609},
{  195,  485, -573, -133},{  427,  202, -171, -118},{  199,  575,    2,  -31},
{  694,  755,-1366,  -39},{  552,  557, -489,  271},{  680,  537,   13, -453},
{  855,  954, -133,  -52},{  -81,  738,-1169,  637},{ 1055, 1059,  -95,  676},
{ 1259, 1081,  489,  305},{ -449,  954, -534,  996},{ -969,  866,-1058, 1059},
{-1294,  618,-1416,  617},{ -458, 1366, -159, 1821},{ -774, -528,  -14, 1110},
{-1202, -901, -772,  433},{-1256,-1255,-1011, -302},{ -602, -585, -759,-1618},
{ -760,-1549, -840,-1921},{ -816, -539,-1769,-2235},{ -227,  -36,-2034,-1831},
{-2107,-1126,-2471,-1816},{-1470,  252,-2701, -415},{ -571, -467, 1509, 1554},
{ 2180, 1975, 2326, 2020}
};

static const int16_t lsf_5_4[256][4] = {
{-1857,-1681,-1857,-1755},{-2056,-1150,-2134,-1654},{-1619,-1099,-1704,-1131},
{-1345,-1608,-1359,-1638},{-1338,-1293,-1325,-1265},{-1664,-1649,-1487, -851},
{-1346,-1832,-1413,-2188},{-1282, -681,-1785,-1649},{ -966,-1082,-1183,-1676},
{-1054,-1073,-1142,-1158},{-1207, -744,-1274, -997},{ -934,-1383, -927,-1416},
{-1010,-1305, -783, -955},{-1049, -900, -993, -817},{ -737, -823, -972,-1189},
{ -738,-1094, -738,-1154},{ -784, -801, -810, -786},{ -892, -520,-1000, -818},
{ -644, -965, -577, -882},{ -541, -694, -671, -917},{ -595, -642, -646, -615},
{ -956, -621, -925, -515},{ -727, -483, -815, -485},{ -840, -578, -440, -713},
{ -578, -325, -657, -670},{ -386, -570, -441, -666},{ -514, -787, -392, -529},
{ -522, -453, -487, -423},{ -616, -585, -617, -157},{ -662, -268, -680, -348},
{ -322, -323, -632, -444},{ -304, -430, -332, -458},{ -277, -468, -659, -793},
{ -319, -636, -227, -554},{ -373, -347, -334, -210},{ -456, -192, -530, -242},
{ -216, -198, -366, -370},{ -338, -161, -409, -748},{ -107, -380, -294, -643},
{ -223, -665, -234, -741},{ -141, -496, -130, -510},{ -139, -327, -172, -305},
{ -306, -580, -164, -263},{ -262, -172,  -67, -402},{   31, -366,  -10, -436},
{  -86, -527,   71, -377},{  -22, -609,  -12, -678},{  -67, -319,   63, -191},
{   35, -181,  -39, -242},{  126, -167, -140, -544},{  155, -297,  174, -297},
{   38,   -8,  117, -380},{  197, -452,  240, -522},{  223, -103,  110, -187},
{   87, -155,  169,  -47},{  157,   26,  -83, -100},{  128,   80,  209,  -62},
{    6,    7,   22,    5},{  318,  -20,  248,  -45},{ -200,  -63,  156,  -69},
{  250, -183,  369, -126},{ -113,  -76, -142, -122},{  -64, -254,  -31,   35},
{ -177,  -71,   -7,  171},{   93,   27,  108,  212},{ -330, -209, -123,  -70},
{ -279,   95,  -96,   20},{ -188,  -61, -314,   87},{ -300,  -78, -354, -134},
{   11,  122, -140,  122},{ -275,  152, -293,  140},{  -82,  138, -321, -111},
{ -480, -156, -359,   76},{ -254,  -40, -635,  -96},{ -522,   79, -507,    8},
{ -268,  303, -539,   68},{ -446,   61, -522,  306},{  111,  189, -435,  122},
{ -379,  166, -571, -398},{ -632,  -74, -747,  -95},{ -455,  194, -952,   83},
{ -798,  192, -755,  192},{ -781, -162, -619,  234},{ -663, -297, -488, -109},
{ -964, -132, -838,  -68},{ -843,   58,-1112,  -86},{ -805, -299, -944, -253},
{ -778,  -50, -965, -549},{ -352,  -98, -992, -343},{-1117, -315,-1117, -307},
{-1155, -374, -637, -230},{-1166,  -43,-1299, -100},{ -925, -393,-1274, -600},
{ -689, -130,-1479, -312},{-1321, -254,-1464, -442},{-1292, -613,-1261, -503},
{-1501, -368,-1322,   26},{-1432,  -66,-1743, -161},{-1644, -467,-1760, -548},
{-1393, -568,-1556, -871},{-1495,-1034,-1387, -571},{-1917, -528,-1783, -123},
{-1897, -231,-2054, -323},{-2052, -906,-1976, -567},{-1917, -620,-2047, -989},
{-1077, -370,-2031, -704},{-2355, -749,-2740,-1089},{-1909,  159,-2012,  248},
{ -626, -123,-2339, -962},{ -669, -408,-1379,-1174},{ -452, -364,-1044, -735},
{ -132,  183,-1620, -752},{ -547, -307, -777,-1261},{  -98,   41, -880,-1091},
{ -257,   97,-1602,-1833},{   31,  -26, -644, -561},{ -180, -546, -385,-1095},
{ -410, -802, -414, -827},{ -457, -970, -490,-1109},{ -215, -916, -144, -937},
{ -493,-1269, -517,-1507},{  181,  101, -332, -889},{ -836, -937, -559, -429},
{ -629, -547, -183, -337},{ -545,  -82, -250, -286},{    5, -132, -348, -252},
{ -293, -472, -158,  100},{  -29,  197, -236, -424},{ -861, -213, -140,   -7},
{ -427, -443,  187,  -97},{ -684, -736, -293,  258},{ -368, -152, -150,  392},
{ -609,  175, -142,  299},{ -138,  152, -119,  329},{ -486,  -52,  293,  198},
{ -183,  117,  175,  331},{  -58, -274,  231,  300},{ -288,  330, -305,  372},
{ -111,  409,   -9,  423},{   83,  256,   67,  367},{  -19,  248,   91,  113},
{  -35,  406, -191,  154},{  238,  296,    5,  197},{  141,  221,  313,  198},
{  211,  421,  244,  334},{   88,  426, -243,  454},{  202,  552,   -5,  403},
{  291,  185,  219,  301},{  251,  138,  128,   69},{  197,  288, -140,  -61},
{  188,  361,  197,  598},{  442,  273,  290,  143},{  472,  482,  157,  370},
{  415,  321,  372,  385},{  402,  552,  155,   24},{  550,  263,  -11,   21},
{  360,  227,  147, -254},{  424,   97,  366,  -13},{  375,  141,  449,  232},
{  396,  507,  474,  272},{  701,  324,  362,  -47},{  587,  148,  543,   69},
{  400,  -51,  561,   59},{  220,  -10,  352,  147},{  206,  211,  653,  185},
{  563,  297,  565,  284},{  594,  121,  766,  192},{  398,  118,  642,  434},
{  233,  264,  481,  467},{  129, -165,  699,  239},{   90,   26,  342,  474},
{  -55,   27,  388,   94},{ -172,    0,  725,  379},{  -60,  337,  370,  465},
{   95,  319,  806,  595},{   78,  260,  497,  851},{  210,  560,  458,  574},
{ -464,  202,  497,  625},{ -202,  152,   48,  712},{  -20,  566,  100,  715},
{  455,  468,  411,  605},{  319,  646,  195,  615},{  401,  538,  680,  739},
{  201,  667,  434,  954},{  454,  425,  646,  491},{  606,  681,  416,  508},
{  497,  822,  426,  815},{  660,  647,  628,  716},{  697,  466,  618,  457},
{  685,  460,  365,  309},{  721,  567,  836,  601},{  609,  300,  825,  459},
{  943,  687,  681,  533},{  915,  598,  591,  243},{  876,  451,  874,  420},
{  786,  317,  732,  220},{  922,  317, 1108,  367},{  531,  466, 1028,  649},
{ 1053,  615, 1034,  553},{  829,  602, 1021,  799},{  927,  803,  878,  763},
{  799,  496, 1373,  773},{  585,  770,  803,  930},{ 1099,  793, 1222,  862},
{ 1209,  895, 1025,  727},{  772,  845, 1172, 1115},{  867, 1021,  830, 1013},
{  841,  910,  506,  703},{ 1239, 1077,  620,  819},{ 1196, 1083, 1155, 1081},
{ 1142,  907, 1547, 1121},{ 1309,  648, 1343,  612},{ 1484,  988, 1479,  937},
{  985, 1328,  955, 1341},{  429,  910,  841, 1338},{  564, 1179,  412, 1156},
{ 1427, 1320, 1434, 1330},{  640,  760, 1726, 1410},{  190,  555, 1073, 1005},
{  426,  257,  839,  980},{  235,  231, 1520, 1167},{  109,  293, 1014, 1569},
{  305,  142, 1148,  539},{ -291, -108, 1213,  972},{   22, -216,  667,  828},
{ -482,  438,  453, 1431},{ -581, -422,  789,  387},{ -358, -454,  174,  780},
{  -36, -372,  390, -134},{ -629,  160, -306,  751},{-1258, -331,  177,  522},
{ -248,  574, -251,  639},{ -531,  407, -596,  394},{ -419,  789, -617,  801},
{ -986,  399, -857,  727},{   -7,  518, -703,  310},{-1143,  -24,-1002,  287},
{ -960,  363,-1299,  312},{-1534,  245,-1557,  305},{   28,  153, -859, -175},
{  -33,  332,-1398, -154},{  212,  410, -593, -197},{-1092, -704, -904,  -65},
{  282,  367, -918, -686},{  345,   93, -258, -357},{  696,  644, -693,  -28},
{  448,  493, -273,  193},{  527,  546, -243, -513},{  384, -136,  273, -353},
{  512, -142,  537, -198},{  941,  750,   83,  248},{  578,  861,  -56,  592},
{  842,   44,  892,   24},{   33,  890,  -16,  982},{  831, 1398, 1535, 1898},
{ 1716, 1376, 1948, 1465}
};

static const int16_t lsf_5_5[64][4] = {
{-1002, -929,-1096,-1203},{ -641, -931, -604, -961},{ -779, -673, -835, -788},
{ -416, -664, -458, -766},{ -652, -521, -662, -495},{-1023, -509,-1023, -428},
{ -444, -552, -368, -449},{ -479, -211,-1054, -903},{ -316, -249, -569, -591},
{ -569, -275, -541, -191},{ -716, -188, -842, -264},{ -333, -248, -318, -228},
{ -275,    1, -567, -228},{ -115, -221, -238, -374},{ -197, -507, -222, -579},
{ -258, -432,  -61, -244},{ -345,    2, -338,   39},{ -215, -169,  -58,    0},
{  -56,   -6, -203, -131},{    1, -186,   -5, -211},{    6, -380,   11, -418},
{ -116,  131, -134,  113},{   89,   -4,   71,   -2},{  -19, -192,  262,   24},
{  189,  151, -133, -109},{  186, -153,  166, -219},{   37,  139,  193,  171},
{  337,  124,  158,  -61},{  141,  226,  -13,  190},{  231,   34,  354,  109},
{  316,  201,  244,  164},{  330,  -85,  390,  -84},{  254,  327,  257,  335},
{  491,  147,  476,  105},{   54,   77,  437,  370},{  421,  314,  449,  342},
{  329,  126,  673,  292},{  571,  388,  243,  193},{  653,  320,  621,  280},
{  194,  380,  517,  581},{   45,  323,  111,  422},{  489,  395,  734,  534},
{  622,  546,  486,  502},{  318,  572,  189,  550},{  385,  422, -157,  153},
{ -125,  382, -197,  386},{ -263,  334,  228,  697},{ -188,    1,   51,  297},
{ -507,  213, -376,  397},{  -24,  255, -547,   89},{ -502,  -94,  387,  179},
{ -620,   68, -684,  112},{ -642, -350, -260,  172},{ -438, -324,  264,  648},
{ -964,   -4,-1121,    7},{ -134,  134,-1133, -306},{  143,   96, -420, -497},
{-1221, -350,-1527, -685},{ -161,   72,  873,  691},{  732,  283,  921,  353},
{  334,  475, 1095,  821},{  864,  524,  843,  497},{  714,  711,  788,  750},
{ 1076,  714, 1204,  753}
};

static const float lsf_3_mean[LP_FILTER_ORDER] = {
 377.441,  554.688,  922.363, 1339.84 , 1702.15 ,
2046.39 , 2452.88 , 2741.46 , 3116.70 , 3348.14 ,
};

static const float lsf_5_mean[LP_FILTER_ORDER] = {
 337.891,  507.080,  834.961, 1247.07 , 1646.00 ,
1982.91 , 2407.96 , 2708.01 , 3104.00 , 3344.97 ,
};

/** Prediction factor table for modes other than 12.2kbit/s */
static const float pred_fac[LP_FILTER_ORDER] = {
0.291626, 0.328644, 0.383636, 0.405640, 0.438873,
0.355560, 0.323120, 0.298065, 0.262238, 0.197876,
};

/** Prediction factor for 12.2kbit/s mode */
#define PRED_FAC_MODE_122              0.65

#define LSF_R_FAC          (8000.0/32768.0) ///< LSF residual tables to Hertz
#define MIN_LSF_SPACING             50.0488 ///< Ensures stability of LPC filter
#define PITCH_LAG_MAX                   143 ///< Upper bound on decoded lag search
#define PITCH_LAG_MIN                    20 ///< Lower bound on decoded lag search
#define PITCH_LAG_MIN_MODE_122           18 ///< Lower bound on decoded lag search in 12.2kbit/s mode

/** b60 hamming windowed sinc function coefficients */
static const float b60[61] = {
 0.898529  ,  0.865051  ,  0.769257  ,  0.624054  ,  0.448639  ,  0.265289   ,
 0.0959167 , -0.0412598 , -0.134338  , -0.178986  , -0.178528  , -0.142609   ,
-0.0849304 , -0.0205078 ,  0.0369568 ,  0.0773926 ,  0.0955200 ,  0.0912781  ,
 0.0689392 ,  0.0357056 ,  0.        , -0.0305481 , -0.0504150 , -0.0570068  ,
-0.0508423 , -0.0350037 , -0.0141602 ,  0.00665283,  0.0230713 ,  0.0323486  ,
 0.0335388 ,  0.0275879 ,  0.0167847 ,  0.00411987, -0.00747681, -0.0156860  ,
-0.0193481 , -0.0183716 , -0.0137634 , -0.00704956,  0.        ,  0.00582886 ,
 0.00939941,  0.0103760 ,  0.00903320,  0.00604248,  0.00238037, -0.00109863 ,
-0.00366211, -0.00497437, -0.00503540, -0.00402832, -0.00241089, -0.000579834,
 0.00103760,  0.00222778,  0.00277710,  0.00271606,  0.00213623,  0.00115967 ,
 0.
};


// fixed tables

/** In 12.2kbit/s mode, positions are divided into TRACKS classes. */
#define TRACKS          5
/** In 10.2kbit/s mode, positions are divided into TRACKS_MODE_102 classes. */
#define TRACKS_MODE_102 4

/**
 * number of pulses per mode
 */
static const uint8_t pulses_nb_per_mode[] = {2, 2, 2, 3, 4, 4, 8, 10};

/** track start positions for algebraic code book routines */
static const uint8_t track_position[16] = { 0, 2, 0, 3, 0, 2, 0, 3, 1, 3, 2, 4, 1, 4, 1, 4 };

/** 3-bit Gray code to binary lookup table */
static const uint8_t gray_decode[8] = { 0, 1, 3, 2, 5, 6, 4, 7 };


// gain tables

/** Initial energy in dB. Also used for bad frames (unimplemented). */
#define MIN_ENERGY -14.0

/** scalar quantized pitch gain table for 7.95 and 12.2 kbps modes */
static const uint16_t qua_gain_pit[16] = {
     0,  3277,  6556,  8192,  9830, 11469, 12288, 13107,
 13926, 14746, 15565, 16384, 17203, 18022, 18842, 19661
};

/** scalar quantized fixed gain table for 7.95 and 12.2 kbps modes */
static const uint16_t qua_gain_code[32] = {
   159,   206,   268,   349,   419,   482,   554,   637,
   733,   842,   969,  1114,  1281,  1473,  1694,  1948,
  2241,  2577,  2963,  3408,  3919,  4507,  5183,  5960,
  6855,  7883,  9065, 10425, 12510, 16263, 21142, 27485
};

/** desired mean innovation energy, indexed by active mode */
static const float energy_mean[8] = { 33.0, 33.0, 33.0, 28.75, 30.0, 36.0, 33.0, 36.0 };

/** 4-tap moving average prediction coefficients in reverse order */
static const float energy_pred_fac[4] = { 0.19, 0.34, 0.58, 0.68 };

/** gain table for 4.75 kbps mode
 *
 * first index has even/odd indexes for subframes 0,2/1,3
 * second index is {pitch_gain, fixed_gain_factor} */
static const uint16_t gains_MODE_475[512][2] = {
{  812,  128}, {  542,  140}, { 2873, 1135}, { 2266, 3402}, { 2067,  563},
{12677,  647}, { 4132, 1798}, { 5601, 5285}, { 7689,  374}, { 3735,  441},
{10912, 2638}, {11807, 2494}, {20490,  797}, { 5218,  675}, { 6724, 8354},
{ 5282, 1696}, { 1488,  428}, { 5882,  452}, { 5332, 4072}, { 3583, 1268},
{ 2469,  901}, {15894, 1005}, {14982, 3271}, {10331, 4858}, { 3635, 2021},
{ 2596,  835}, {12360, 4892}, {12206, 1704}, {13432, 1604}, { 9118, 2341},
{ 3968, 1538}, { 5479, 9936}, { 3795,  417}, { 1359,  414}, { 3640, 1569},
{ 7995, 3541}, {11405,  645}, { 8552,  635}, { 4056, 1377}, {16608, 6124},
{11420,  700}, { 2007,  607}, {12415, 1578}, {11119, 4654}, {13680, 1708},
{11990, 1229}, { 7996, 7297}, {13231, 5715}, { 2428, 1159}, { 2073, 1941},
{ 6218, 6121}, { 3546, 1804}, { 8925, 1802}, { 8679, 1580}, {13935, 3576},
{13313, 6237}, { 6142, 1130}, { 5994, 1734}, {14141, 4662}, {11271, 3321},
{12226, 1551}, {13931, 3015}, { 5081,10464}, { 9444, 6706}, { 1689,  683},
{ 1436, 1306}, { 7212, 3933}, { 4082, 2713}, { 7793,  704}, {15070,  802},
{ 6299, 5212}, { 4337, 5357}, { 6676,  541}, { 6062,  626}, {13651, 3700},
{11498, 2408}, {16156,  716}, {12177,  751}, { 8065,11489}, { 6314, 2256},
{ 4466,  496}, { 7293,  523}, {10213, 3833}, { 8394, 3037}, { 8403,  966},
{14228, 1880}, { 8703, 5409}, {16395, 4863}, { 7420, 1979}, { 6089, 1230},
{ 9371, 4398}, {14558, 3363}, {13559, 2873}, {13163, 1465}, { 5534, 1678},
{13138,14771}, { 7338,  600}, { 1318,  548}, { 4252, 3539}, {10044, 2364},
{10587,  622}, {13088,  669}, {14126, 3526}, { 5039, 9784}, {15338,  619},
{ 3115,  590}, {16442, 3013}, {15542, 4168}, {15537, 1611}, {15405, 1228},
{16023, 9299}, { 7534, 4976}, { 1990, 1213}, {11447, 1157}, {12512, 5519},
{ 9475, 2644}, { 7716, 2034}, {13280, 2239}, {16011, 5093}, { 8066, 6761},
{10083, 1413}, { 5002, 2347}, {12523, 5975}, {15126, 2899}, {18264, 2289},
{15827, 2527}, {16265,10254}, {14651,11319}, { 1797,  337}, { 3115,  397},
{ 3510, 2928}, { 4592, 2670}, { 7519,  628}, {11415,  656}, { 5946, 2435},
{ 6544, 7367}, { 8238,  829}, { 4000,  863}, {10032, 2492}, {16057, 3551},
{18204, 1054}, { 6103, 1454}, { 5884, 7900}, {18752, 3468}, { 1864,  544},
{ 9198,  683}, {11623, 4160}, { 4594, 1644}, { 3158, 1157}, {15953, 2560},
{12349, 3733}, {17420, 5260}, { 6106, 2004}, { 2917, 1742}, {16467, 5257},
{16787, 1680}, {17205, 1759}, { 4773, 3231}, { 7386, 6035}, {14342,10012},
{ 4035,  442}, { 4194,  458}, { 9214, 2242}, { 7427, 4217}, {12860,  801},
{11186,  825}, {12648, 2084}, {12956, 6554}, { 9505,  996}, { 6629,  985},
{10537, 2502}, {15289, 5006}, {12602, 2055}, {15484, 1653}, {16194, 6921},
{14231, 5790}, { 2626,  828}, { 5615, 1686}, {13663, 5778}, { 3668, 1554},
{11313, 2633}, { 9770, 1459}, {14003, 4733}, {15897, 6291}, { 6278, 1870},
{ 7910, 2285}, {16978, 4571}, {16576, 3849}, {15248, 2311}, {16023, 3244},
{14459,17808}, {11847, 2763}, { 1981, 1407}, { 1400,  876}, { 4335, 3547},
{ 4391, 4210}, { 5405,  680}, {17461,  781}, { 6501, 5118}, { 8091, 7677},
{ 7355,  794}, { 8333, 1182}, {15041, 3160}, {14928, 3039}, {20421,  880},
{14545,  852}, {12337,14708}, { 6904, 1920}, { 4225,  933}, { 8218, 1087},
{10659, 4084}, {10082, 4533}, { 2735,  840}, {20657, 1081}, {16711, 5966},
{15873, 4578}, {10871, 2574}, { 3773, 1166}, {14519, 4044}, {20699, 2627},
{15219, 2734}, {15274, 2186}, { 6257, 3226}, {13125,19480}, { 7196,  930},
{ 2462, 1618}, { 4515, 3092}, {13852, 4277}, {10460,  833}, {17339,  810},
{16891, 2289}, {15546, 8217}, {13603, 1684}, { 3197, 1834}, {15948, 2820},
{15812, 5327}, {17006, 2438}, {16788, 1326}, {15671, 8156}, {11726, 8556},
{ 3762, 2053}, { 9563, 1317}, {13561, 6790}, {12227, 1936}, { 8180, 3550},
{13287, 1778}, {16299, 6599}, {16291, 7758}, { 8521, 2551}, { 7225, 2645},
{18269, 7489}, {16885, 2248}, {17882, 2884}, {17265, 3328}, { 9417,20162},
{11042, 8320}, { 1286,  620}, { 1431,  583}, { 5993, 2289}, { 3978, 3626},
{ 5144,  752}, {13409,  830}, { 5553, 2860}, {11764, 5908}, {10737,  560},
{ 5446,  564}, {13321, 3008}, {11946, 3683}, {19887,  798}, { 9825,  728},
{13663, 8748}, { 7391, 3053}, { 2515,  778}, { 6050,  833}, { 6469, 5074},
{ 8305, 2463}, { 6141, 1865}, {15308, 1262}, {14408, 4547}, {13663, 4515},
{ 3137, 2983}, { 2479, 1259}, {15088, 4647}, {15382, 2607}, {14492, 2392},
{12462, 2537}, { 7539, 2949}, {12909,12060}, { 5468,  684}, { 3141,  722},
{ 5081, 1274}, {12732, 4200}, {15302,  681}, { 7819,  592}, { 6534, 2021},
{16478, 8737}, {13364,  882}, { 5397,  899}, {14656, 2178}, {14741, 4227},
{14270, 1298}, {13929, 2029}, {15477, 7482}, {15815, 4572}, { 2521, 2013},
{ 5062, 1804}, { 5159, 6582}, { 7130, 3597}, {10920, 1611}, {11729, 1708},
{16903, 3455}, {16268, 6640}, { 9306, 1007}, { 9369, 2106}, {19182, 5037},
{12441, 4269}, {15919, 1332}, {15357, 3512}, {11898,14141}, {16101, 6854},
{ 2010,  737}, { 3779,  861}, {11454, 2880}, { 3564, 3540}, { 9057, 1241},
{12391,  896}, { 8546, 4629}, {11561, 5776}, { 8129,  589}, { 8218,  588},
{18728, 3755}, {12973, 3149}, {15729,  758}, {16634,  754}, {15222,11138},
{15871, 2208}, { 4673,  610}, {10218,  678}, {15257, 4146}, { 5729, 3327},
{ 8377, 1670}, {19862, 2321}, {15450, 5511}, {14054, 5481}, { 5728, 2888},
{ 7580, 1346}, {14384, 5325}, {16236, 3950}, {15118, 3744}, {15306, 1435},
{14597, 4070}, {12301,15696}, { 7617, 1699}, { 2170,  884}, { 4459, 4567},
{18094, 3306}, {12742,  815}, {14926,  907}, {15016, 4281}, {15518, 8368},
{17994, 1087}, { 2358,  865}, {16281, 3787}, {15679, 4596}, {16356, 1534},
{16584, 2210}, {16833, 9697}, {15929, 4513}, { 3277, 1085}, { 9643, 2187},
{11973, 6068}, { 9199, 4462}, { 8955, 1629}, {10289, 3062}, {16481, 5155},
{15466, 7066}, {13678, 2543}, { 5273, 2277}, {16746, 6213}, {16655, 3408},
{20304, 3363}, {18688, 1985}, {14172,12867}, {15154,15703}, { 4473, 1020},
{ 1681,  886}, { 4311, 4301}, { 8952, 3657}, { 5893, 1147}, {11647, 1452},
{15886, 2227}, { 4582, 6644}, { 6929, 1205}, { 6220,  799}, {12415, 3409},
{15968, 3877}, {19859, 2109}, { 9689, 2141}, {14742, 8830}, {14480, 2599},
{ 1817, 1238}, { 7771,  813}, {19079, 4410}, { 5554, 2064}, { 3687, 2844},
{17435, 2256}, {16697, 4486}, {16199, 5388}, { 8028, 2763}, { 3405, 2119},
{17426, 5477}, {13698, 2786}, {19879, 2720}, { 9098, 3880}, {18172, 4833},
{17336,12207}, { 5116,  996}, { 4935,  988}, { 9888, 3081}, { 6014, 5371},
{15881, 1667}, { 8405, 1183}, {15087, 2366}, {19777, 7002}, {11963, 1562},
{ 7279, 1128}, {16859, 1532}, {15762, 5381}, {14708, 2065}, {20105, 2155},
{17158, 8245}, {17911, 6318}, { 5467, 1504}, { 4100, 2574}, {17421, 6810},
{ 5673, 2888}, {16636, 3382}, { 8975, 1831}, {20159, 4737}, {19550, 7294},
{ 6658, 2781}, {11472, 3321}, {19397, 5054}, {18878, 4722}, {16439, 2373},
{20430, 4386}, {11353,26526}, {11593, 3068}, { 2866, 1566}, { 5108, 1070},
{ 9614, 4915}, { 4939, 3536}, { 7541,  878}, {20717,  851}, { 6938, 4395},
{16799, 7733}, {10137, 1019}, { 9845,  964}, {15494, 3955}, {15459, 3430},
{18863,  982}, {20120,  963}, {16876,12887}, {14334, 4200}, { 6599, 1220},
{ 9222,  814}, {16942, 5134}, { 5661, 4898}, { 5488, 1798}, {20258, 3962},
{17005, 6178}, {17929, 5929}, { 9365, 3420}, { 7474, 1971}, {19537, 5177},
{19003, 3006}, {16454, 3788}, {16070, 2367}, { 8664, 2743}, { 9445,26358},
{10856, 1287}, { 3555, 1009}, { 5606, 3622}, {19453, 5512}, {12453,  797},
{20634,  911}, {15427, 3066}, {17037,10275}, {18883, 2633}, { 3913, 1268},
{19519, 3371}, {18052, 5230}, {19291, 1678}, {19508, 3172}, {18072,10754},
{16625, 6845}, { 3134, 2298}, {10869, 2437}, {15580, 6913}, {12597, 3381},
{11116, 3297}, {16762, 2424}, {18853, 6715}, {17171, 9887}, {12743, 2605},
{ 8937, 3140}, {19033, 7764}, {18347, 3880}, {20475, 3682}, {19602, 3380},
{13044,19373}, {10526,23124}
};

/** gain table for 6.70, 7.40 and 10.2 kbps modes
 *
 * second index is {pitch_gain, fixed_gain_factor} */
static const uint16_t gains_high[128][2] = {
{  577,  662}, {  806, 1836}, { 3109, 1052}, { 4181, 1387}, { 2373, 1425},
{ 3248, 1985}, { 1827, 2320}, {  941, 3314}, { 2351, 2977}, { 3616, 2420},
{ 3451, 3096}, { 2955, 4301}, { 1848, 4500}, { 3884, 5416}, { 1187, 7210},
{ 3083, 9000}, { 7384,  883}, { 5962, 1506}, { 5155, 2134}, { 7944, 2009},
{ 6507, 2250}, { 7670, 2752}, { 5952, 3016}, { 4898, 3764}, { 6989, 3588},
{ 8174, 3978}, { 6064, 4404}, { 7709, 5087}, { 5523, 6021}, { 7769, 7126},
{ 6060, 7938}, { 5594,11487}, {10581, 1356}, { 9049, 1597}, { 9794, 2035},
{ 8946, 2415}, {10296, 2584}, { 9407, 2734}, { 8700, 3218}, { 9757, 3395},
{10177, 3892}, { 9170, 4528}, {10152, 5004}, { 9114, 5735}, {10500, 6266},
{10110, 7631}, { 8844, 8727}, { 8956,12496}, {12924,  976}, {11435, 1755},
{12138, 2328}, {11388, 2368}, {10700, 3064}, {12332, 2861}, {11722, 3327},
{11270, 3700}, {10861, 4413}, {12082, 4533}, {11283, 5205}, {11960, 6305},
{11167, 7534}, {12128, 8329}, {10969,10777}, {10300,17376}, {13899, 1681},
{12580, 2045}, {13265, 2439}, {14033, 2989}, {13452, 3098}, {12396, 3658},
{13510, 3780}, {12880, 4272}, {13533, 4861}, {12667, 5457}, {13854, 6106},
{13031, 6483}, {13557, 7721}, {12957, 9311}, {13714,11551}, {12591,15206},
{15113, 1540}, {15072, 2333}, {14527, 2511}, {14692, 3199}, {15382, 3560},
{14133, 3960}, {15102, 4236}, {14332, 4824}, {14846, 5451}, {15306, 6083},
{14329, 6888}, {15060, 7689}, {14406, 9426}, {15387, 9741}, {14824,14271},
{13600,24939}, {16396, 1969}, {16817, 2832}, {15713, 2843}, {16104, 3336},
{16384, 3963}, {16940, 4579}, {15711, 4599}, {16222, 5448}, {16832, 6382},
{15745, 7141}, {16326, 7469}, {16611, 8624}, {17028,10418}, {15905,11817},
{16878,14690}, {16515,20870}, {18142, 2083}, {19401, 3178}, {17508, 3426},
{20054, 4027}, {18069, 4249}, {18952, 5066}, {17711, 5402}, {19835, 6192},
{17950, 7014}, {21318, 7877}, {17910, 9289}, {19144, 9290}, {20517,11381},
{18075,14485}, {19999,17882}, {18842,32764}
};

/** gain table for 5.15 and 5.90 kbps modes
 *
 * second index is {pitch_gain, fixed_gain_factor} */
static const uint16_t gains_low[64][2] = {
{10813,28753}, {20480, 2785}, {18841, 6594}, { 6225, 7413}, {17203,10444},
{21626, 1269}, {21135, 4423}, {11304, 1556}, {19005,12820}, {17367, 2498},
{17858, 4833}, { 9994, 2498}, {17530, 7864}, {14254, 1884}, {15892, 3153},
{ 6717, 1802}, {18186,20193}, {18022, 3031}, {16711, 5857}, { 8847, 4014},
{15892, 8970}, {18022, 1392}, {16711, 4096}, { 8192,  655}, {15237,13926},
{14254, 3112}, {14090, 4669}, { 5406, 2703}, {13434, 6553}, {12451,  901},
{12451, 2662}, { 3768,  655}, {14745,23511}, {19169, 2457}, {20152, 5079},
{ 6881, 4096}, {20480, 8560}, {19660,  737}, {19005, 4259}, { 7864, 2088},
{11468,12288}, {15892, 1474}, {15728, 4628}, { 9175, 1433}, {16056, 7004},
{14827,  737}, {15073, 2252}, { 5079, 1228}, {13271,17326}, {16547, 2334},
{15073, 5816}, { 3932, 3686}, {14254, 8601}, {16875,  778}, {15073, 3809},
{ 6062,  614}, { 9338, 9256}, {13271, 1761}, {13271, 3522}, { 2457, 1966},
{11468, 5529}, {10485,  737}, {11632, 3194}, { 1474,  778}
};


// pre-processing tables

/** Maximum sharpening factor
 *
 * The specification says 0.8, which should be 13107, but the reference C code
 * uses 13017 instead. (Amusingly the same applies to SHARP_MAX in g729dec.c.)
 */
#define SHARP_MAX 0.79449462890625

/** impulse response filter tables converted to float from Q15 int32_t
 * used for anti-sparseness processing */
static const float ir_filter_strong_MODE_795[AMR_SUBFRAME_SIZE] = {
 0.817169,  0.024445,  0.076447, -0.020844, -0.042175,  0.017761,  0.018433,
-0.038879,  0.107147, -0.179871,  0.138367, -0.015228, -0.059204,  0.091888,
-0.154358,  0.171326, -0.060730, -0.032379, -0.044525,  0.135559, -0.021362,
-0.162811,  0.140656,  0.013794, -0.017975, -0.102295,  0.090118,  0.038666,
-0.036987, -0.079041,  0.052826,  0.112000, -0.136566, -0.029755,  0.134003,
-0.077423,  0.028961, -0.041595, -0.029877,  0.174988,
};

static const float ir_filter_strong[AMR_SUBFRAME_SIZE] = {
 0.448303,  0.351501,  0.038696, -0.084259, -0.173065,  0.229309, -0.001068,
-0.085663, -0.092773,  0.147186,  0.090088, -0.257080,  0.115509,  0.044403,
 0.066498, -0.263580,  0.245697, -0.064178, -0.044373,  0.023712,  0.033813,
-0.072784,  0.068787, -0.011078, -0.020569, -0.064178,  0.184509, -0.173370,
 0.032715,  0.095306, -0.154358,  0.162109, -0.071075, -0.113770,  0.211304,
-0.118683,  0.020599, -0.054169,  0.000885,  0.309601,
};

static const float ir_filter_medium[AMR_SUBFRAME_SIZE] = {
 0.923889,  0.116913, -0.123169,  0.090698, -0.031982, -0.030579,  0.075592,
-0.092865,  0.085907, -0.068085,  0.053497, -0.049164,  0.052307, -0.054169,
 0.047089, -0.030762,  0.013092, -0.005157,  0.014404, -0.038574,  0.066406,
-0.082581,  0.076996, -0.049469,  0.010498,  0.025208, -0.046661,  0.052612,
-0.050568,  0.051910, -0.062958,  0.080688, -0.093384,  0.088409, -0.060364,
 0.016998,  0.023804, -0.041779,  0.025696,  0.019989,
};

static const float *ir_filters_lookup[2]          = { ir_filter_strong,          ir_filter_medium };
static const float *ir_filters_lookup_MODE_795[2] = { ir_filter_strong_MODE_795, ir_filter_medium };


// postfilter tables

/** Powers of formant factor gamma_n=0.7 for modes 12.2 and 10.2kbit/s */
static const float formant_high_n[10] = {
0.700000, 0.490000, 0.343000, 0.240100, 0.168070,
0.117649, 0.082354, 0.057648, 0.040354, 0.028248
};
/** Powers of formant factor gamma_d=0.75 for modes 12.2 and 10.2kbit/s */
static const float formant_high_d[10] = {
0.750000, 0.562500, 0.421875, 0.316406, 0.237305,
0.177979, 0.133484, 0.100113, 0.075085, 0.056314
};
/** Powers of formant factor gamma_n=0.55 for modes less than 10.2kbit/s */
static const float formant_low_n[10] = {
0.550000, 0.302500, 0.166375, 0.091506, 0.050328,
0.027681, 0.015224, 0.008373, 0.004605, 0.002533
};
/** Powers of formant factor gamma_d=0.7 for modes less than 10.2kbit/s */
static const float *formant_low_d = formant_high_n;

/** Number of impulse response coefficients used for tilt factor */
#define AMR_TILT_RESPONSE 22
/** Tilt factor = 1st reflection coefficient * gamma_t */
#define AMR_TILT_GAMMA_T 0.8
/** Adaptive gain control factor used in post-filter */
#define AMR_AGC_ALPHA 0.9


#endif /* AVCODEC_AMRNBDATA_H */

