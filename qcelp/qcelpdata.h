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

#ifndef AVCODEC_QCELPDATA_H
#define AVCODEC_QCELPDATA_H

/**
 * @file qcelpdata.h
 * Data tables for the QCELP decoder
 * @author Reynaldo H. Verdejo Pinochet
 */

typedef enum
{
    RATE_UNKNOWN = -2,
    I_F_Q,             /*!< insufficient frame quality */
    BLANK,
    RATE_OCTAVE,
    RATE_QUARTER,
    RATE_HALF,
    RATE_FULL
} qcelp_packet_rate;

static const uint16_t qcelp_bits_per_rate[]={0, 20, 54, 124, 266};
static const float    qcelp_hammsinc_table[]={-0.006822, 0.041249,-0.143459,
                                               0.588863, 0.588863,-0.143459,
                                               0.041249,-0.006822};

typedef struct
{
    uint8_t index;  /*!< index into the reference frame */
    uint8_t bitpos; /*!< bit position in the value's byte */
} QCELPBitmap;


/**
 * WARNING
 *
 * YOU WILL NOT SEE ANY mention of a REFERENCE nor an UNIVERSAL frame
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
 * CINDEXs  32    33 34 35 36 37 38 39 40 41 42 43 44 45 46 47
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
 *---------------------------------------------------------------------------*/


#define QCELP_RATE_FULL_BITMAP \
{62,2},{62,1},{62,0},{61,6},{61,5},{61,4},{61,3},{61,2},\
{61,1},{61,0},{60,5},{60,4},{60,3},{60,2},{60,1},{60,0},\
{64,5},{64,4},{64,3},{64,2},{64,1},{64,0},{63,5},{63,4},\
{63,3},{63,2},{63,1},{63,0},{62,6},{62,5},{62,4},{62,3},\
{ 0,0},{16,3},{16,2},{16,1},{16,0},{52,0},{48,6},{48,5},\
{48,4},{48,3},{48,2},{48,1},{48,0},{56,2},{56,1},{56,0},\
{33,3},{33,2},{33,1},{33,0},{ 1,0},{17,3},{17,2},{17,1},\
{17,0},{32,6},{32,5},{32,4},{32,3},{32,2},{32,1},{32,0},\
{19,0},{34,6},{34,5},{34,4},{34,3},{34,2},{34,1},{34,0},\
{ 2,0},{18,3},{18,2},{18,1},{18,0},{33,6},{33,5},{33,4},\
{49,2},{49,1},{49,0},{57,2},{57,1},{57,0},{35,6},{35,5},\
{35,4},{35,3},{35,2},{35,1},{35,0},{ 3,0},{19,2},{19,1},\
{36,5},{36,4},{36,3},{36,2},{36,1},{36,0},{ 4,0},{20,3},\
{20,2},{20,1},{20,0},{53,0},{49,6},{49,5},{49,4},{49,3},\
{22,2},{22,1},{22,0},{37,6},{37,5},{37,4},{37,3},{37,2},\
{37,1},{37,0},{ 5,0},{21,3},{21,2},{21,1},{21,0},{36,6},\
{39,2},{39,1},{39,0},{ 7,0},{23,2},{23,1},{23,0},{38,6},\
{38,5},{38,4},{38,3},{38,2},{38,1},{38,0},{ 6,0},{22,3},\
{24,0},{54,0},{50,6},{50,5},{50,4},{50,3},{50,2},{50,1},\
{50,0},{58,2},{58,1},{58,0},{39,6},{39,5},{39,4},{39,3},\
{ 9,0},{25,3},{25,2},{25,1},{25,0},{40,6},{40,5},{40,4},\
{40,3},{40,2},{40,1},{40,0},{ 8,0},{24,3},{24,2},{24,1},\
{42,3},{42,2},{42,1},{42,0},{10,0},{26,3},{26,2},{26,1},\
{26,0},{41,6},{41,5},{41,4},{41,3},{41,2},{41,1},{41,0},\
{59,1},{59,0},{43,6},{43,5},{43,4},{43,3},{43,2},{43,1},\
{43,0},{11,0},{27,2},{27,1},{27,0},{42,6},{42,5},{42,4},\
{44,1},{44,0},{12,0},{28,3},{28,2},{28,1},{28,0},{55,0},\
{51,6},{51,5},{51,4},{51,3},{51,2},{51,1},{51,0},{59,2},\
{45,5},{45,4},{45,3},{45,2},{45,1},{45,0},{13,0},{29,3},\
{29,2},{29,1},{29,0},{44,6},{44,5},{44,4},{44,3},{44,2},\
{31,2},{31,1},{31,0},{46,6},{46,5},{46,4},{46,3},{46,2},\
{46,1},{46,0},{14,0},{30,3},{30,2},{30,1},{30,0},{45,6},\
{65,1},{65,0},{47,6},{47,5},{47,4},{47,3},{47,2},{47,1},\
{47,0},{15,0}

#define QCELP_RATE_HALF_BITMAP \
{62,2},{62,1},{62,0},{61,6},{61,5},{61,4},{61,3},{61,2},\
{61,1},{61,0},{60,5},{60,4},{60,3},{60,2},{60,1},{60,0},\
{64,5},{64,4},{64,3},{64,2},{64,1},{64,0},{63,5},{63,4},\
{63,3},{63,2},{63,1},{63,0},{62,6},{62,5},{62,4},{62,3},\
{ 0,0},{16,3},{16,2},{16,1},{16,0},{52,1},{48,6},{48,5},\
{48,4},{48,3},{48,2},{48,1},{48,0},{56,2},{56,1},{56,0},\
{49,5},{49,4},{49,3},{49,2},{49,1},{49,0},{57,2},{57,1},\
{57,0},{32,6},{32,5},{32,4},{32,3},{32,2},{32,1},{32,0},\
{58,1},{58,0},{33,6},{33,5},{33,4},{33,3},{33,2},{33,1},\
{33,0},{ 1,0},{17,3},{17,2},{17,1},{17,0},{53,0},{49,6},\
{34,1},{34,0},{ 2,0},{18,3},{18,2},{18,1},{18,0},{54,0},\
{50,6},{50,5},{50,4},{50,3},{50,2},{50,1},{50,0},{58,2},\
{55,0},{51,6},{51,5},{51,4},{51,3},{51,2},{51,1},{51,0},\
{59,2},{59,1},{59,0},{34,6},{34,5},{34,4},{34,3},{34,2},\
{35,6},{35,5},{35,4},{35,3},{35,2},{35,1},{35,0},{ 3,0},\
{19,3},{19,2},{19,1},{19,0}

#define QCELP_RATE_4THR_BITMAP \
{62,2},{62,1},{62,0},{61,6},{61,5},{61,4},{61,3},{61,2},\
{61,1},{61,0},{60,5},{60,4},{60,3},{60,2},{60,1},{60,0},\
{64,5},{64,4},{64,3},{64,2},{64,1},{64,0},{63,5},{63,4},\
{63,3},{63,2},{63,1},{63,0},{62,6},{62,5},{62,4},{62,3},\
{19,3},{19,2},{19,1},{19,0},{18,3},{18,2},{18,1},{18,0},\
{17,3},{17,2},{17,1},{17,0},{16,3},{16,2},{16,1},{16,0},\
{65,1},{65,0},{20,3},{20,2},{20,1},{20,0}

#define QCELP_RATE_8THR_BITMAP \
{65,0},{65,1},{65,2},{65,3},{16,0},{16,1},{75,0},{76,0},\
{74,0},{73,0},{72,0},{76,1},{71,0},{70,0},{69,0},{76,2},\
{68,0},{67,0},{66,0},{76,3}

/**
 * position of the bitmapping data for each pkt type in
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
 * position of the transmission codes inside the universal frame
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

typedef struct
{
    uint16_t x;
    uint16_t y;
} qcelp_vector;

/**
 * LSP Vector quantization tables in x*10000 form
 *
 * TIA/EIA/IS-733 tables 2.4.3.2.6.3-1 through 2.4.3.2.6.3-5
 */

static const qcelp_vector qcelp_lspvq1[]={
{ 327, 118},{ 919, 111},{ 427, 440},{1327, 185},
{ 469,  50},{1272,  91},{ 892,  59},{1771, 193},
{ 222, 158},{1100, 127},{ 827,  55},{ 978, 791},
{ 665,  47},{ 700,1401},{ 670, 859},{1913,1048},
{ 471, 215},{1046, 125},{ 645, 298},{1599, 160},
{ 593,  39},{1187, 462},{ 749, 341},{1520, 511},
{ 290, 792},{ 909, 362},{ 753,  81},{1111,1058},
{ 519, 253},{ 828, 839},{ 685, 541},{1421,1258},
{ 386, 130},{ 962, 119},{ 542, 387},{1431, 185},
{ 526,  51},{1175, 260},{ 831, 167},{1728, 510},
{ 273, 437},{1172, 113},{ 771, 144},{1122, 751},
{ 619, 119},{ 492,1276},{ 658, 695},{1882, 615},
{ 415, 200},{1018,  88},{ 681, 339},{1436, 325},
{ 555, 122},{1042, 485},{ 826, 345},{1374, 743},
{ 383,1018},{1005, 358},{ 704,  86},{1301, 586},
{ 597, 241},{ 832, 621},{ 555, 573},{1504, 839}};

static const qcelp_vector qcelp_lspvq2[]={
{ 255, 293},{ 904, 219},{ 151,1211},{1447, 498},
{ 470, 253},{1559, 177},{1547, 994},{2394, 242},
{  91, 813},{ 857, 590},{ 934,1326},{1889, 282},
{ 813, 472},{1057,1494},{ 450,3315},{2163,1895},
{ 538, 532},{1399, 218},{ 146,1552},{1755, 626},
{ 822, 202},{1299, 663},{ 706,1732},{2656, 401},
{ 418, 745},{ 762,1038},{ 583,1748},{1746,1285},
{ 527,1169},{1314, 830},{ 556,2116},{1073,2321},
{ 297, 570},{ 981, 403},{ 468,1103},{1740, 243},
{ 725, 179},{1255, 474},{1374,1362},{1922, 912},
{ 285, 947},{ 930, 700},{ 593,1372},{1909, 576},
{ 588, 916},{1110,1116},{ 224,2719},{1633,2220},
{ 402, 520},{1061, 448},{ 402,1352},{1499, 775},
{ 664, 589},{1081, 727},{ 801,2206},{2165,1157},
{ 566, 802},{ 911,1116},{ 306,1703},{1792, 836},
{ 655, 999},{1061,1038},{ 298,2089},{1110,1753},
{ 361, 311},{ 970, 239},{ 265,1231},{1495, 573},
{ 566, 262},{1569, 293},{1341,1144},{2271, 544},
{ 214, 877},{ 847, 719},{ 794,1384},{2067, 274},
{ 703, 688},{1099,1306},{ 391,2947},{2024,1670},
{ 471, 525},{1245, 290},{ 264,1557},{1568, 807},
{ 718, 399},{1193, 685},{ 883,1594},{2729, 764},
{ 500, 754},{ 809,1108},{ 541,1648},{1523,1385},
{ 614,1196},{1209, 847},{ 345,2242},{1442,1747},
{ 199, 560},{1092, 194},{ 349,1253},{1653, 507},
{ 625, 354},{1376, 431},{1187,1465},{2164, 872},
{ 360, 974},{1008, 698},{ 704,1346},{2114, 452},
{ 720, 816},{1240,1089},{ 439,2475},{1498,2040},
{ 336, 718},{1213, 187},{ 451,1450},{1368, 885},
{ 592, 578},{1131, 531},{ 861,1855},{1764,1500},
{ 444, 970},{ 935, 903},{ 424,1687},{1633,1102},
{ 793, 897},{1060, 897},{ 185,2011},{1205,1855}};

static const qcelp_vector qcelp_lspvq3[]={
{ 225, 283},{1296, 355},{ 543, 343},{2073, 274},
{ 204,1099},{1562, 523},{1388, 161},{2784, 274},
{ 112, 849},{1870, 175},{1189, 160},{1490,1088},
{ 969,1115},{ 659,3322},{1158,1073},{3183,1363},
{ 517, 223},{1740, 223},{ 704, 387},{2637, 234},
{ 692,1005},{1287,1610},{ 952, 532},{2393, 646},
{ 490, 552},{1619, 657},{ 845, 670},{1784,2280},
{ 191,1775},{ 272,2868},{ 942, 952},{2628,1479},
{ 278, 579},{1565, 218},{ 814, 180},{2379, 187},
{ 276,1444},{1199,1223},{1200, 349},{3009, 307},
{ 312, 844},{1898, 306},{ 863, 470},{1685,1241},
{ 513,1727},{ 711,2233},{1085, 864},{3398, 527},
{ 414, 440},{1356, 612},{ 964, 147},{2173, 738},
{ 465,1292},{ 877,1749},{1104, 689},{2105,1311},
{ 580, 864},{1895, 752},{ 652, 609},{1485,1699},
{ 514,1400},{ 386,2131},{ 933, 798},{2473, 986},
{ 334, 360},{1375, 398},{ 621, 276},{2183, 280},
{ 311,1114},{1382, 807},{1284, 175},{2605, 636},
{ 230, 816},{1739, 408},{1074, 176},{1619,1120},
{ 784,1371},{ 448,3050},{1189, 880},{3039,1165},
{ 424, 241},{1672, 186},{ 815, 333},{2432, 324},
{ 584,1029},{1137,1546},{1015, 585},{2198, 995},
{ 574, 581},{1746, 647},{ 733, 740},{1938,1737},
{ 347,1710},{ 373,2429},{ 787,1061},{2439,1438},
{ 185, 536},{1489, 178},{ 703, 216},{2178, 487},
{ 154,1421},{1414, 994},{1103, 352},{3072, 473},
{ 408, 819},{2055, 168},{ 998, 354},{1917,1140},
{ 665,1799},{ 993,2213},{1234, 631},{3003, 762},
{ 373, 620},{1518, 425},{ 913, 300},{1966, 836},
{ 402,1185},{ 948,1385},{1121, 555},{1802,1509},
{ 474, 886},{1888, 610},{ 739, 585},{1231,2379},
{ 661,1335},{ 205,2211},{ 823, 822},{2480,1179}};

static const qcelp_vector qcelp_lspvq4[]={
{ 348, 311},{ 812,1145},{ 552, 461},{1826, 263},
{ 601, 675},{1730, 172},{1523, 193},{2449, 277},
{ 334, 668},{ 805,1441},{1319, 207},{1684, 910},
{ 582,1318},{1403,1098},{ 979, 832},{2700,1359},
{ 624, 228},{1292, 979},{ 800, 195},{2226, 285},
{ 730, 862},{1537, 601},{1115, 509},{2720, 354},
{ 218,1167},{1212,1538},{1074, 247},{1674,1710},
{ 322,2142},{1263, 777},{ 981, 556},{2119,1710},
{ 193, 596},{1035, 957},{ 694, 397},{1997, 253},
{ 743, 603},{1584, 321},{1346, 346},{2221, 708},
{ 451, 732},{1040,1415},{1184, 230},{1853, 919},
{ 310,1661},{1625, 706},{ 856, 843},{2902, 702},
{ 467, 348},{1108,1048},{ 859, 306},{1964, 463},
{ 560,1013},{1425, 533},{1142, 634},{2391, 879},
{ 397,1084},{1345,1700},{ 976, 248},{1887,1189},
{ 644,2087},{1262, 603},{ 877, 550},{2203,1307}};

static const qcelp_vector qcelp_lspvq5[]={
{ 360, 222},{ 820,1097},{ 601, 319},{1656, 198},
{ 604, 513},{1552, 141},{1391, 155},{2474, 261},
{ 269, 785},{1463, 646},{1123, 191},{2015, 223},
{ 785, 844},{1202,1011},{ 980, 807},{3014, 793},
{ 570, 180},{1135,1382},{ 778, 256},{1901, 179},
{ 807, 622},{1461, 458},{1231, 178},{2028, 821},
{ 387, 927},{1496,1004},{ 888, 392},{2246, 341},
{ 295,1462},{1156, 694},{1022, 473},{2226,1364},
{ 210, 478},{1029,1020},{ 722, 181},{1730, 251},
{ 730, 488},{1465, 293},{1303, 326},{2595, 387},
{ 458, 584},{1569, 742},{1029, 173},{1910, 495},
{ 605,1159},{1268, 719},{ 973, 646},{2872, 428},
{ 443, 334},{ 835,1465},{ 912, 138},{1716, 442},
{ 620, 778},{1316, 450},{1186, 335},{1446,1665},
{ 486,1050},{1675,1019},{ 880, 278},{2214, 202},
{ 539,1564},{1142, 533},{ 984, 391},{2130,1089}};

/**
 * table for computing Ga (decoded linear codebook gain magnitude)
 *
 * TIA/EIA/IS-733 2.4.6.2.1-3
 */

static const float qcelp_g12ga[]={
   1.000,   1.125,   1.250,   1.375,   1.625,  1.750,  2.000,  2.250,
   2.500,   2.875,   3.125,   3.500,   4.000,  4.500,  5.000,  5.625,
   6.250,   7.125,   8.000,   8.875,  10.000, 11.250, 12.625, 14.125,
  15.875,  17.750,  20.000,  22.375,  25.125, 28.125, 31.625, 35.500,
  39.750,  44.625,  50.125,  56.250,  63.125, 70.750, 79.375, 89.125,
 100.000, 112.250, 125.875, 141.250, 158.500,177.875,199.500,223.875,
 251.250, 281.875, 316.250, 354.875, 398.125,446.625,501.125,563.375,
 631.000, 708.000, 794.375, 891.250,1000.000};

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

#define QCELP_SQRT1887 1.373681186

static const double qcelp_rnd_fir_coefs[]=
{
            0,-1.344519e-1, 1.735384e-2,-6.905826e-2,
  2.434368e-2,-8.210701e-2, 3.041388e-2,-9.251384e-2,
  3.501983e-2,-9.918777e-2, 3.749518e-2, 8.985137e-1,
  3.749518e-2,-9.918777e-2, 3.501983e-2,-9.251384e-2,
  3.041388e-2,-8.210701e-2, 2.434368e-2,-6.905826e-2,
  1.735384e-2,-1.344519e-1
}; /*!< Start reading from [1]. */

#endif /* AVCODEC_QCELPDATA_H */
