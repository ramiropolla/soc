/*
 * ALS header file for common data
 * Copyright (c) 2009 Thilo Borgmann <thilo.borgmann _at_ googlemail.com>
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

#ifndef AVCODEC_ALS_DATA_H
#define AVCODEC_ALS_DATA_H

/**
 * @file libavcodec/als_data.h
 * MPEG-4 ALS common data tables
 * @author Thilo Borgmann <thilo.borgmann _at_ googlemail.com>
 */


#include <stdint.h>

/** Rice parameters and corresponding index offsets for decoding the
 *  indices of scaled PARCOR values. The table choosen is set globally
 *  by the encoder and stored in ALSSpecificConfig.
 */
static const int8_t parcor_rice_table[3][20][2] = {
    { {-52, 4}, {-29, 5}, {-31, 4}, { 19, 4}, {-16, 4},
      { 12, 3}, { -7, 3}, {  9, 3}, { -5, 3}, {  6, 3},
      { -4, 3}, {  3, 3}, { -3, 2}, {  3, 2}, { -2, 2},
      {  3, 2}, { -1, 2}, {  2, 2}, { -1, 2}, {  2, 2} },
    { {-58, 3}, {-42, 4}, {-46, 4}, { 37, 5}, {-36, 4},
      { 29, 4}, {-29, 4}, { 25, 4}, {-23, 4}, { 20, 4},
      {-17, 4}, { 16, 4}, {-12, 4}, { 12, 3}, {-10, 4},
      {  7, 3}, { -4, 4}, {  3, 3}, { -1, 3}, {  1, 3} },
    { {-59, 3}, {-45, 5}, {-50, 4}, { 38, 4}, {-39, 4},
      { 32, 4}, {-30, 4}, { 25, 3}, {-23, 3}, { 20, 3},
      {-20, 3}, { 16, 3}, {-13, 3}, { 10, 3}, { -7, 3},
      {  3, 3}, {  0, 3}, { -1, 3}, {  2, 3}, { -1, 2} }
};


/** Scaled PARCOR values used for the first two PARCOR coefficients.
 *  To be indexed by the Rice coded indices.
 *  Generated by: parcor_scaled_values[i] = 32 + ((i * (i+1)) << 7) - (1 << 20)
 *  Actual values are divided by 32 in order to be stored in 16 bits.
 */
static const int16_t parcor_scaled_values[] = {
    -1048544 / 32, -1048288 / 32, -1047776 / 32, -1047008 / 32,
    -1045984 / 32, -1044704 / 32, -1043168 / 32, -1041376 / 32,
    -1039328 / 32, -1037024 / 32, -1034464 / 32, -1031648 / 32,
    -1028576 / 32, -1025248 / 32, -1021664 / 32, -1017824 / 32,
    -1013728 / 32, -1009376 / 32, -1004768 / 32,  -999904 / 32,
     -994784 / 32,  -989408 / 32,  -983776 / 32,  -977888 / 32,
     -971744 / 32,  -965344 / 32,  -958688 / 32,  -951776 / 32,
     -944608 / 32,  -937184 / 32,  -929504 / 32,  -921568 / 32,
     -913376 / 32,  -904928 / 32,  -896224 / 32,  -887264 / 32,
     -878048 / 32,  -868576 / 32,  -858848 / 32,  -848864 / 32,
     -838624 / 32,  -828128 / 32,  -817376 / 32,  -806368 / 32,
     -795104 / 32,  -783584 / 32,  -771808 / 32,  -759776 / 32,
     -747488 / 32,  -734944 / 32,  -722144 / 32,  -709088 / 32,
     -695776 / 32,  -682208 / 32,  -668384 / 32,  -654304 / 32,
     -639968 / 32,  -625376 / 32,  -610528 / 32,  -595424 / 32,
     -580064 / 32,  -564448 / 32,  -548576 / 32,  -532448 / 32,
     -516064 / 32,  -499424 / 32,  -482528 / 32,  -465376 / 32,
     -447968 / 32,  -430304 / 32,  -412384 / 32,  -394208 / 32,
     -375776 / 32,  -357088 / 32,  -338144 / 32,  -318944 / 32,
     -299488 / 32,  -279776 / 32,  -259808 / 32,  -239584 / 32,
     -219104 / 32,  -198368 / 32,  -177376 / 32,  -156128 / 32,
     -134624 / 32,  -112864 / 32,   -90848 / 32,   -68576 / 32,
      -46048 / 32,   -23264 / 32,     -224 / 32,    23072 / 32,
       46624 / 32,    70432 / 32,    94496 / 32,   118816 / 32,
      143392 / 32,   168224 / 32,   193312 / 32,   218656 / 32,
      244256 / 32,   270112 / 32,   296224 / 32,   322592 / 32,
      349216 / 32,   376096 / 32,   403232 / 32,   430624 / 32,
      458272 / 32,   486176 / 32,   514336 / 32,   542752 / 32,
      571424 / 32,   600352 / 32,   629536 / 32,   658976 / 32,
      688672 / 32,   718624 / 32,   748832 / 32,   779296 / 32,
      810016 / 32,   840992 / 32,   872224 / 32,   903712 / 32,
      935456 / 32,   967456 / 32,   999712 / 32,  1032224 / 32
};


/** Gain values of p(0) for long-term prediction.
 *  To be indexed by the Rice coded indices.
 */
static const int ltp_gain [] = {
    0, 8, 16, 24, 32, 40, 48, 56, 64, 70, 76, 82, 88, 92, 96, 100
};

#endif /* AVCODEC_ALS_DATA_H */
