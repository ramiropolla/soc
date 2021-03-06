/* -*-  indent-tabs-mode:nil; c-basic-offset:4;  -*- */
/*
 * Copyright (C) 2007 Marco Gerards <marco@gnu.org>
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
 * @file
 * Arithmetic decoder for Dirac
 * @author Marco Gerards <marco@gnu.org>
 */

#include "dirac_arith.h"

static uint16_t arith_lookup[256] = {
    0,    2,    5,    8,    11,   15,   20,   24,
    29,   35,   41,   47,   53,   60,   67,   74,
    82,   89,   97,   106,  114,  123,  132,  141,
    150,  160,  170,  180,  190,  201,  211,  222,
    233,  244,  256,  267,  279,  291,  303,  315,
    327,  340,  353,  366,  379,  392,  405,  419,
    433,  447,  461,  475,  489,  504,  518,  533,
    548,  563,  578,  593,  609,  624,  640,  656,
    672,  688,  705,  721,  738,  754,  771,  788,
    805,  822,  840,  857,  875,  892,  910,  928,
    946,  964,  983,  1001, 1020, 1038, 1057, 1076,
    1095, 1114, 1133, 1153, 1172, 1192, 1211, 1231,
    1251, 1271, 1291, 1311, 1332, 1352, 1373, 1393,
    1414, 1435, 1456, 1477, 1498, 1520, 1541, 1562,
    1584, 1606, 1628, 1649, 1671, 1694, 1716, 1738,
    1760, 1783, 1806, 1828, 1851, 1874, 1897, 1920,
    1935, 1942, 1949, 1955, 1961, 1968, 1974, 1980,
    1985, 1991, 1996, 2001, 2006, 2011, 2016, 2021,
    2025, 2029, 2033, 2037, 2040, 2044, 2047, 2050,
    2053, 2056, 2058, 2061, 2063, 2065, 2066, 2068,
    2069, 2070, 2071, 2072, 2072, 2072, 2072, 2072,
    2072, 2071, 2070, 2069, 2068, 2066, 2065, 2063,
    2060, 2058, 2055, 2052, 2049, 2045, 2042, 2038,
    2033, 2029, 2024, 2019, 2013, 2008, 2002, 1996,
    1989, 1982, 1975, 1968, 1960, 1952, 1943, 1934,
    1925, 1916, 1906, 1896, 1885, 1874, 1863, 1851,
    1839, 1827, 1814, 1800, 1786, 1772, 1757, 1742,
    1727, 1710, 1694, 1676, 1659, 1640, 1622, 1602,
    1582, 1561, 1540, 1518, 1495, 1471, 1447, 1422,
    1396, 1369, 1341, 1312, 1282, 1251, 1219, 1186,
    1151, 1114, 1077, 1037, 995,  952,  906,  857,
    805, 750,   690,  625,  553,  471,  376,  255
};

static uint8_t next_ctx[ARITH_CONTEXT_COUNT] = {
    ARITH_CONTEXT_ZP_F2,    // ZPZN_F1
    ARITH_CONTEXT_ZP_F2,    // ZPNN_F1
    ARITH_CONTEXT_NP_F2,    // NPZN_F1
    ARITH_CONTEXT_NP_F2,    // NPNN_F1
    ARITH_CONTEXT_ZP_F3,    // ZP_F2
    ARITH_CONTEXT_ZP_F4,    // ZP_F3
    ARITH_CONTEXT_ZP_F5,    // ZP_F4
    ARITH_CONTEXT_ZP_F6,    // ZP_F5
    ARITH_CONTEXT_ZP_F6,    // ZP_F6
    ARITH_CONTEXT_NP_F3,    // NP_F2
    ARITH_CONTEXT_NP_F4,    // NP_F3
    ARITH_CONTEXT_NP_F5,    // NP_F4
    ARITH_CONTEXT_NP_F6,    // NP_F5
    ARITH_CONTEXT_NP_F6,    // NP_F6
    0,
    0,
    0,
    0,
    0,
    ARITH_CONTEXT_Q_OFFSET_FOLLOW,
    0,
    0,

    ARITH_CONTEXT_SB_F2,    // SB_F1
    ARITH_CONTEXT_SB_F2,    // SB_F2
    0,
    0,
    0,
    0,
    ARITH_CONTEXT_VECTOR_F2,    // VECTOR_F1
    ARITH_CONTEXT_VECTOR_F3,    // VECTOR_F2
    ARITH_CONTEXT_VECTOR_F4,    // VECTOR_F3
    ARITH_CONTEXT_VECTOR_F5,    // VECTOR_F4
    ARITH_CONTEXT_VECTOR_F5,    // VECTOR_F5
    0,
    0,
    ARITH_CONTEXT_DC_F2,        // DC_F2
    ARITH_CONTEXT_DC_F2,        // DC_F2
};


static void dirac_arith_init_common(dirac_arith_state *arith)
{
    int i;

    arith->range = 0xFFFF;

    /* Initialize contexts. */
    for (i = 0; i < ARITH_CONTEXT_COUNT; i++) {
        arith->contexts[i] = 0x8000;
    }
}

/**
 * Initialize arithmetic decoding
 *
 * @param arith state used for further calls to the arithmetic decoder
 * @param gb GetBitContext to read from
 * @param length amount of bytes to decode
 */
void dirac_init_arith_decoder(dirac_arith_state *arith,
                              GetBitContext *gb, int length)
{
    align_get_bits(gb);
    arith->pb        = NULL;
    arith->counter   = 16;

    arith->bytestream_start =
    arith->bytestream       = gb->buffer + get_bits_count(gb)/8;
    arith->bytestream_end   = arith->bytestream_start + length;
    skip_bits_long(gb, length*8);

    arith->low  = (*arith->bytestream++) << 24;
    arith->low |= (*arith->bytestream++) << 16;
    arith->low |= (*arith->bytestream++) << 8;
    arith->low |= (*arith->bytestream++);

    dirac_arith_init_common(arith);
}

void dirac_init_arith_encoder(dirac_arith_state *arith, PutBitContext *pb)
{
    arith->pb        = pb;
    arith->carry     = 0;
    arith->low       = 0;

    dirac_arith_init_common(arith);
}

static inline void renorm_arith_decoder(dirac_arith_state *arith)
{
    while (arith->range <= 0x4000) {
        arith->low   <<= 1;
        arith->range <<= 1;

        if (!--arith->counter) {
            if (arith->bytestream < arith->bytestream_end)
                arith->low |= *arith->bytestream << 8;
            else
                arith->low |= 0xff00;
            arith->bytestream++;

            if (arith->bytestream < arith->bytestream_end)
                arith->low |= *arith->bytestream;
            else
                arith->low |= 0xff;
            arith->bytestream++;

            arith->counter = 16;
        }
    }
}

/**
 * Read a single bit using the arithmetic decoder
 *
 * @param arith state of arithmetic decoder
 * @param context the context of the bit to read
 * @return the bit read
 */
int dirac_get_arith_bit(dirac_arith_state *arith, int context)
{
    unsigned int prob_zero = arith->contexts[context];
    unsigned int range_times_prob;
    unsigned int ret;

    assert(!arith->pb);

    range_times_prob  = (arith->range * prob_zero) >> 16;
    if ((arith->low >> 16) >= range_times_prob) {
        ret = 1;
        arith->low   -= range_times_prob << 16;
        arith->range -= range_times_prob;
        arith->contexts[context] -= arith_lookup[arith->contexts[context] >> 8];
    } else {
        ret = 0;
        arith->range  = range_times_prob;
        arith->contexts[context] += arith_lookup[255 - (arith->contexts[context] >> 8)];
    }

    renorm_arith_decoder(arith);

    return ret;
}

/**
 * Write a single bit using the arithmetic coder
 *
 * @param arith state of arithmetic coder
 * @param context the context of the bit to write
 * @param bit the bit to write
 */
void dirac_put_arith_bit(dirac_arith_state *arith, int context, int bit)
{
    PutBitContext *pb = arith->pb;
    unsigned int prob_zero = arith->contexts[context];

    assert(!arith->gb);

    if (bit == 0) {
        arith->range  = (arith->range * prob_zero) >> 16;
        arith->contexts[context] += arith_lookup[255 - (arith->contexts[context] >> 8)];
    } else {
        arith->low   += (arith->range * prob_zero) >> 16;
        arith->range -= (arith->range * prob_zero) >> 16;
        arith->contexts[context] -= arith_lookup[arith->contexts[context] >> 8];
    }

    /* Renormalisation and output. */
    while (arith->range <= 0x4000) {
        if (((arith->low + arith->range - 1)^arith->low) >= 0x8000) {
            arith->low ^= 0x4000;
            arith->carry++;
        } else {
            put_bits(pb, 1, (arith->low >> 15)&1);
            while (arith->carry > 0) {
                put_bits(pb, 1, !((arith->low >> 15)&1));
                arith->carry--;
            }
        }
        arith->low   <<= 1;
        arith->range <<= 1;
        arith->low    &= 0xFFFF;
    }
}

/**
 * Read an unsigned int using the arithmetic decoder
 * @param arith state of arithmetic decoder
 * @param context_set the collection of contexts to read the unsigned int
 * @return value read by arithmetic decoder
 */
unsigned int dirac_get_arith_uint(dirac_arith_state *arith,
                                  int follow_ctx, int data_ctx)
{
    int ret = 1;

    while (dirac_get_arith_bit (arith, follow_ctx) == 0) {
        ret <<= 1;
        if (dirac_get_arith_bit (arith, data_ctx))
            ret++;
        follow_ctx = next_ctx[follow_ctx];
    }
    ret--;
    return ret;
}

/**
 * Write an unsigned int using the arithmetic coder
 *
 * @param arith        state of arithmetic coder
 * @param context_set  the collection of contexts used to write the unsigned int
 * @param i            value to write
 */
void dirac_put_arith_uint(dirac_arith_state *arith,
                          int follow_ctx, int data_ctx,
                          unsigned int i)
{
    int log = av_log2(++i);
    while (log) {
        log--;
        dirac_put_arith_bit(arith, follow_ctx, 0);
        dirac_put_arith_bit(arith, data_ctx, (i >> log)&1);
        follow_ctx = next_ctx[follow_ctx];
    }
    dirac_put_arith_bit(arith, follow_ctx, 1);
}

/**
 * Read a signed int using the arithmetic decoder
 * @param arith state of arithmetic decoder
 * @param context_set the collection of contexts to read the signed int
 * @return value read by arithmetic decoder
 */
int dirac_get_arith_int(dirac_arith_state *arith, int follow_ctx,
                        int data_ctx, int sign_ctx)
{
    int ret = dirac_get_arith_uint(arith, follow_ctx, data_ctx);
    if (ret != 0 && dirac_get_arith_bit(arith, sign_ctx))
        ret = -ret;
    return ret;
}

/**
 * Write a signed int using the arithmetic coder
 *
 * @param arith        state of arithmetic coder
 * @param context_set  the collection of contexts used to write the signed int
 * @param i            value to write
 */
void dirac_put_arith_int(dirac_arith_state *arith,
                         int follow_ctx, int data_ctx, int sign_ctx, int i)
{
    dirac_put_arith_uint(arith, follow_ctx, data_ctx, FFABS(i));
    if (i)
        dirac_put_arith_bit(arith, sign_ctx, i < 0);
}


/**
 * Flush the arithmetic decoder, consume all bytes up to the
 * initialized length.
 */
void dirac_get_arith_terminate(dirac_arith_state *arith)
{
    assert(!arith->pb);
}


/**
 * Flush the arithmetic coder.
 */
void dirac_put_arith_terminate(dirac_arith_state *arith)
{
    int i;
    int rem;
    assert(!arith->gb);

    /* Output remaining resolved msbs. */
    while (((arith->low + arith->range - 1)^arith->low) < 0x8000) {
        put_bits(arith->pb, 1, (arith->low >> 15)&1);
        while (arith->carry > 0) {
            put_bits(arith->pb, 1, !((arith->low >> 15)&1));
            arith->carry--;
        }
        arith->low   <<= 1;
        arith->low    &= 0xFFFF;
        arith->range <<= 1;
    }

    /* Resolve remaining straddle conditions. */
    while ((arith->low & 0x4000)
           && !((arith->low + arith->range - 1) & 0x4000)) {
        arith->carry++;
        arith->low    ^= 0x4000;
        arith->low   <<= 1;
        arith->range <<= 1;
        arith->low    &= 0xFFFF;
    }

    /* Discharge carry bits. */
    put_bits(arith->pb, 1, (arith->low >> 14)&1);
    for (i = 0; i <= arith->carry; i++)
        put_bits(arith->pb, 1, !((arith->low >> 14)&1));

    /* Add padding bits. */
    rem = 8 - (put_bits_count(arith->pb) % 8);
    for (i = 0; i < rem; i++)
        put_bits(arith->pb, 1, 0);
}

#if 0
void dirac_arith_test(void)
{
    struct dirac_arith_state arith;
    char in[] = "**** Test arithmetic coding and decoding ****";
    char out[100];
    PutBitContext pb;
    GetBitContext gb;
    char buf[2000];
    int i, c;
    int length;

    /* Code the string. */
    init_put_bits(&pb, buf, sizeof(buf)*8);
    dirac_init_arith_encoder(&arith, &pb);
    for (c = 0; c < sizeof(in); c++) {
        for (i = 0; i < 8; i++) {
            int bit = (in[c] >> (7 - i)) & 1;
            dirac_put_arith_bit(&arith, i, bit);
        }
    }

    dirac_put_arith_uint(&arith, ARITH_CONTEXT_ZPZN_F1,
                                 ARITH_CONTEXT_COEFF_DATA, 50);
    dirac_put_arith_uint(&arith, ARITH_CONTEXT_ZPZN_F1,
                                 ARITH_CONTEXT_COEFF_DATA, 100000);
    dirac_put_arith_uint(&arith, ARITH_CONTEXT_ZPZN_F1,
                                 ARITH_CONTEXT_COEFF_DATA, 0);
    dirac_put_arith_uint(&arith, ARITH_CONTEXT_ZPZN_F1,
                                 ARITH_CONTEXT_COEFF_DATA, 123);
    dirac_put_arith_uint(&arith, ARITH_CONTEXT_ZPZN_F1,
                                 ARITH_CONTEXT_COEFF_DATA, 4321);

    dirac_put_arith_int(&arith, ARITH_CONTEXT_ZPZN_F1,
                                ARITH_CONTEXT_COEFF_DATA,
                                ARITH_CONTEXT_SIGN_ZERO, -100);
    dirac_put_arith_int(&arith, ARITH_CONTEXT_ZPZN_F1,
                                ARITH_CONTEXT_COEFF_DATA,
                                ARITH_CONTEXT_SIGN_ZERO, -12345);
    dirac_put_arith_int(&arith, ARITH_CONTEXT_ZPZN_F1,
                                ARITH_CONTEXT_COEFF_DATA,
                                ARITH_CONTEXT_SIGN_ZERO, 0);
    dirac_put_arith_int(&arith, ARITH_CONTEXT_ZPZN_F1,
                                ARITH_CONTEXT_COEFF_DATA,
                                ARITH_CONTEXT_SIGN_ZERO, 1234);
    dirac_put_arith_int(&arith, ARITH_CONTEXT_ZPZN_F1,
                                ARITH_CONTEXT_COEFF_DATA,
                                ARITH_CONTEXT_SIGN_ZERO, -1);

    dirac_put_arith_terminate(&arith);
    flush_put_bits(&pb);
    length = put_bits_count(&pb);

    /* Now decode the string. */
    init_get_bits(&gb, buf, sizeof(buf)*8);
    dirac_init_arith_decoder(&arith, &gb, length);
    for(c = 0; 1; c++) {
        out[c] = 0;
        for (i = 0; i < 8; i++) {
            int bit = dirac_get_arith_bit(&arith, i);
            out[c] = (out[c] << 1) | bit;
        }
        if (out[c] == 0)
            break;
    }

    for (i = 0; i < 5; i++)
        dprintf(0, "UINT: %d\n", dirac_get_arith_uint(&arith,
                            ARITH_CONTEXT_ZPZN_F1, ARITH_CONTEXT_COEFF_DATA));
    for (i = 0; i < 5; i++)
        dprintf(0, "INT: %d\n", dirac_get_arith_int(&arith,
                              ARITH_CONTEXT_ZPZN_F1, ARITH_CONTEXT_COEFF_DATA,
                              ARITH_CONTEXT_SIGN_ZERO));

    dirac_get_arith_terminate(&arith);

    dprintf(0, "Encoder input : `%s'\n", in);
    dprintf(0, "Decoder output: `%s'\n", out);
    dprintf(0, "size decoder: %d, size encoded: %d\n",
            sizeof(in) * 8, length);
}
#endif
