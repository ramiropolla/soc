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

typedef struct {
    float x;
    float y;
} qcelp_fvector;

static const qcelp_fvector qcelp_lspvq1[]={
{.0327,.0118},{.0919,.0111},{.0427,.0440},{.1327,.0185},
{.0469,.0050},{.1272,.0091},{.0892,.0059},{.1771,.0193},
{.0222,.0158},{.1100,.0127},{.0827,.0055},{.0978,.0791},
{.0665,.0047},{.0700,.1401},{.0670,.0859},{.1913,.1048},
{.0471,.0215},{.1046,.0125},{.0645,.0298},{.1599,.0160},
{.0593,.0039},{.1187,.0462},{.0749,.0341},{.1520,.0511},
{.0290,.0792},{.0909,.0362},{.0753,.0081},{.1111,.1058},
{.0519,.0253},{.0828,.0839},{.0685,.0541},{.1421,.1258},
{.0386,.0130},{.0962,.0119},{.0542,.0387},{.1431,.0185},
{.0526,.0051},{.1175,.0260},{.0831,.0167},{.1728,.0510},
{.0273,.0437},{.1172,.0113},{.0771,.0144},{.1122,.0751},
{.0619,.0119},{.0492,.1276},{.0658,.0695},{.1882,.0615},
{.0415,.0200},{.1018,.0088},{.0681,.0339},{.1436,.0325},
{.0555,.0122},{.1042,.0485},{.0826,.0345},{.1374,.0743},
{.0383,.1018},{.1005,.0358},{.0704,.0086},{.1301,.0586},
{.0597,.0241},{.0832,.0621},{.0555,.0573},{.1504,.0839}};

static const qcelp_fvector qcelp_lspvq2[]={
{.0255,.0293},{.0904,.0219},{.0151,.1211},{.1447,.0498},
{.0470,.0253},{.1559,.0177},{.1547,.0994},{.2394,.0242},
{.0091,.0813},{.0857,.0590},{.0934,.1326},{.1889,.0282},
{.0813,.0472},{.1057,.1494},{.0450,.3315},{.2163,.1895},
{.0538,.0532},{.1399,.0218},{.0146,.1552},{.1755,.0626},
{.0822,.0202},{.1299,.0663},{.0706,.1732},{.2656,.0401},
{.0418,.0745},{.0762,.1038},{.0583,.1748},{.1746,.1285},
{.0527,.1169},{.1314,.0830},{.0556,.2116},{.1073,.2321},
{.0297,.0570},{.0981,.0403},{.0468,.1103},{.1740,.0243},
{.0725,.0179},{.1255,.0474},{.1374,.1362},{.1922,.0912},
{.0285,.0947},{.0930,.0700},{.0593,.1372},{.1909,.0576},
{.0588,.0916},{.1110,.1116},{.0224,.2719},{.1633,.2220},
{.0402,.0520},{.1061,.0448},{.0402,.1352},{.1499,.0775},
{.0664,.0589},{.1081,.0727},{.0801,.2206},{.2165,.1157},
{.0566,.0802},{.0911,.1116},{.0306,.1703},{.1792,.0836},
{.0655,.0999},{.1061,.1038},{.0298,.2089},{.1110,.1753},
{.0361,.0311},{.0970,.0239},{.0265,.1231},{.1495,.0573},
{.0566,.0262},{.1569,.0293},{.1341,.1144},{.2271,.0544},
{.0214,.0877},{.0847,.0719},{.0794,.1384},{.2067,.0274},
{.0703,.0688},{.1099,.1306},{.0391,.2947},{.2024,.1670},
{.0471,.0525},{.1245,.0290},{.0264,.1557},{.1568,.0807},
{.0718,.0399},{.1193,.0685},{.0883,.1594},{.2729,.0764},
{.0500,.0754},{.0809,.1108},{.0541,.1648},{.1523,.1385},
{.0614,.1196},{.1209,.0847},{.0345,.2242},{.1442,.1747},
{.0199,.0560},{.1092,.0194},{.0349,.1253},{.1653,.0507},
{.0625,.0354},{.1376,.0431},{.1187,.1465},{.2164,.0872},
{.0360,.0974},{.1008,.0698},{.0704,.1346},{.2114,.0452},
{.0720,.0816},{.1240,.1089},{.0439,.2475},{.1498,.2040},
{.0336,.0718},{.1213,.0187},{.0451,.1450},{.1368,.0885},
{.0592,.0578},{.1131,.0531},{.0861,.1855},{.1764,.1500},
{.0444,.0970},{.0935,.0903},{.0424,.1687},{.1633,.1102},
{.0793,.0897},{.1060,.0897},{.0185,.2011},{.1205,.1855}};

static const qcelp_fvector qcelp_lspvq3[]={
{.0225,.0283},{.1296,.0355},{.0543,.0343},{.2073,.0274},
{.0204,.1099},{.1562,.0523},{.1388,.0161},{.2784,.0274},
{.0112,.0849},{.1870,.0175},{.1189,.0160},{.1490,.1088},
{.0969,.1115},{.0659,.3322},{.1158,.1073},{.3183,.1363},
{.0517,.0223},{.1740,.0223},{.0704,.0387},{.2637,.0234},
{.0692,.1005},{.1287,.1610},{.0952,.0532},{.2393,.0646},
{.0490,.0552},{.1619,.0657},{.0845,.0670},{.1784,.2280},
{.0191,.1775},{.0272,.2868},{.0942,.0952},{.2628,.1479},
{.0278,.0579},{.1565,.0218},{.0814,.0180},{.2379,.0187},
{.0276,.1444},{.1199,.1223},{.1200,.0349},{.3009,.0307},
{.0312,.0844},{.1898,.0306},{.0863,.0470},{.1685,.1241},
{.0513,.1727},{.0711,.2233},{.1085,.0864},{.3398,.0527},
{.0414,.0440},{.1356,.0612},{.0964,.0147},{.2173,.0738},
{.0465,.1292},{.0877,.1749},{.1104,.0689},{.2105,.1311},
{.0580,.0864},{.1895,.0752},{.0652,.0609},{.1485,.1699},
{.0514,.1400},{.0386,.2131},{.0933,.0798},{.2473,.0986},
{.0334,.0360},{.1375,.0398},{.0621,.0276},{.2183,.0280},
{.0311,.1114},{.1382,.0807},{.1284,.0175},{.2605,.0636},
{.0230,.0816},{.1739,.0408},{.1074,.0176},{.1619,.1120},
{.0784,.1371},{.0448,.3050},{.1189,.0880},{.3039,.1165},
{.0424,.0241},{.1672,.0186},{.0815,.0333},{.2432,.0324},
{.0584,.1029},{.1137,.1546},{.1015,.0585},{.2198,.0995},
{.0574,.0581},{.1746,.0647},{.0733,.0740},{.1938,.1737},
{.0347,.1710},{.0373,.2429},{.0787,.1061},{.2439,.1438},
{.0185,.0536},{.1489,.0178},{.0703,.0216},{.2178,.0487},
{.0154,.1421},{.1414,.0994},{.1103,.0352},{.3072,.0473},
{.0408,.0819},{.2055,.0168},{.0998,.0354},{.1917,.1140},
{.0665,.1799},{.0993,.2213},{.1234,.0631},{.3003,.0762},
{.0373,.0620},{.1518,.0425},{.0913,.0300},{.1966,.0836},
{.0402,.1185},{.0948,.1385},{.1121,.0555},{.1802,.1509},
{.0474,.0886},{.1888,.0610},{.0739,.0585},{.1231,.2379},
{.0661,.1335},{.0205,.2211},{.0823,.0822},{.2480,.1179}};

static const qcelp_fvector qcelp_lspvq4[]={
{.0348,.0311},{.0812,.1145},{.0552,.0461},{.1826,.0263},
{.0601,.0675},{.1730,.0172},{.1523,.0193},{.2449,.0277},
{.0334,.0668},{.0805,.1441},{.1319,.0207},{.1684,.0910},
{.0582,.1318},{.1403,.1098},{.0979,.0832},{.2700,.1359},
{.0624,.0228},{.1292,.0979},{.0800,.0195},{.2226,.0285},
{.0730,.0862},{.1537,.0601},{.1115,.0509},{.2720,.0354},
{.0218,.1167},{.1212,.1538},{.1074,.0247},{.1674,.1710},
{.0322,.2142},{.1263,.0777},{.0981,.0556},{.2119,.1710},
{.0193,.0596},{.1035,.0957},{.0694,.0397},{.1997,.0253},
{.0743,.0603},{.1584,.0321},{.1346,.0346},{.2221,.0708},
{.0451,.0732},{.1040,.1415},{.1184,.0230},{.1853,.0919},
{.0310,.1661},{.1625,.0706},{.0856,.0843},{.2902,.0702},
{.0467,.0348},{.1108,.1048},{.0859,.0306},{.1964,.0463},
{.0560,.1013},{.1425,.0533},{.1142,.0634},{.2391,.0879},
{.0397,.1084},{.1345,.1700},{.0976,.0248},{.1887,.1189},
{.0644,.2087},{.1262,.0603},{.0877,.0550},{.2203,.1307}};

static const qcelp_fvector qcelp_lspvq5[]={
{.0360,.0222},{.0820,.1097},{.0601,.0319},{.1656,.0198},
{.0604,.0513},{.1552,.0141},{.1391,.0155},{.2474,.0261},
{.0269,.0785},{.1463,.0646},{.1123,.0191},{.2015,.0223},
{.0785,.0844},{.1202,.1011},{.0980,.0807},{.3014,.0793},
{.0570,.0180},{.1135,.1382},{.0778,.0256},{.1901,.0179},
{.0807,.0622},{.1461,.0458},{.1231,.0178},{.2028,.0821},
{.0387,.0927},{.1496,.1004},{.0888,.0392},{.2246,.0341},
{.0295,.1462},{.1156,.0694},{.1022,.0473},{.2226,.1364},
{.0210,.0478},{.1029,.1020},{.0722,.0181},{.1730,.0251},
{.0730,.0488},{.1465,.0293},{.1303,.0326},{.2595,.0387},
{.0458,.0584},{.1569,.0742},{.1029,.0173},{.1910,.0495},
{.0605,.1159},{.1268,.0719},{.0973,.0646},{.2872,.0428},
{.0443,.0334},{.0835,.1465},{.0912,.0138},{.1716,.0442},
{.0620,.0778},{.1316,.0450},{.1186,.0335},{.1446,.1665},
{.0486,.1050},{.1675,.1019},{.0880,.0278},{.2214,.0202},
{.0539,.1564},{.1142,.0533},{.0984,.0391},{.2130,.1089}};

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
