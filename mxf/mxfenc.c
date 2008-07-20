/*
 * MXF muxer
 * Copyright (c) 2008 GUCAS, Zhentan Feng<spyfeng at gmail dot com>
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

/*
 * References
 * SMPTE 336M KLV Data Encoding Protocol Using Key-Length-Value
 * SMPTE 377M MXF File Format Specifications
 * SMPTE 379M MXF Generic Container
 * SMPTE 381M Mapping MPEG Streams into the MXF Generic Container
 * SMPTE RP210: SMPTE Metadata Dictionary
 * SMPTE RP224: Registry of SMPTE Universal Labels
 */

#define DEBUG

#include "libavutil/random.h"
#include "avformat.h"
#include "libavcodec/bytestream.h"

typedef uint8_t UID[16];
typedef uint8_t UMID[32];

enum MXFMetadataSetType {
    MaterialPackage,
    SourcePackage,
};

typedef struct {
    UID key;
    offset_t offset;
    uint64_t length;
} KLVPacket;

typedef struct {
    UID uid;
    unsigned matching_len;
    enum CodecID id;
} MXFCodecUL;

typedef struct {
    int local_tag;
    UID uid;
} MXFLocalTagPair;

typedef struct {
    UID uid;
    enum CodecType type;
} MXFDataDefinitionUL;

typedef struct {
    UID uid;
    enum CodecID type;
} MXFEssenceElementKey;

typedef struct {
    UID *identification;
    UID *content_storage;
    UID *package;
    UID *track;
    UID **sequence;
    UID **structural_component;
} MXFReferenceContext;

typedef struct MXFContext {
    UMID top_src_package_uid;
    int64_t header_byte_count;
    int64_t header_start;
    int64_t header_byte_count_offset;
    int64_t header_footer_partition_offset;
    AVRandomState random_state;
    MXFReferenceContext *reference;
    char *track_number_sign;
    UID *track_essence_element_key;
    int type_num;
    const MXFCodecUL *video_container_ul;
    const MXFCodecUL *audio_container_ul;
} MXFContext;

static const uint8_t umid_base[] = {0x06, 0x0a, 0x2b, 0x34, 0x01, 0x01, 0x01, 0x01, 0x01, 0x01, 0x0f, 0x00, 0x13, 0x00, 0x00, 0x00};

/* complete key */
static const uint8_t op1a_ul[]            = { 0x06, 0x0e, 0x2b, 0x34, 0x04, 0x01, 0x01, 0x01, 0x0d, 0x01, 0x02, 0x01, 0x01, 0x01, 0x01, 0x00 };
static const uint8_t header_partition_key[]            = { 0x06, 0x0e, 0x2b, 0x34, 0x02, 0x05, 0x01, 0x01, 0x0d, 0x01, 0x02, 0x01, 0x01, 0x02, 0x04, 0x00 }; // ClosedComplete
static const uint8_t footer_partition_key[] = {0x06, 0x0e, 0x2b, 0x34, 0x02, 0x05, 0x01, 0x01, 0x0d, 0x01, 0x02, 0x01, 0x01, 0x04, 0x04, 0x00}; // ClosedComplete
static const uint8_t primer_pack_key[] = { 0x06,0x0E,0x2B,0x34,0x02,0x05,0x01,0x01,0x0d,0x01,0x02,0x01,0x01,0x05,0x01,0x00 };

static const MXFEssenceElementKey mxf_essence_element_key[] = {
    { { 0x06,0x0e,0x2b,0x34,0x01,0x02,0x01,0x01,0x0d,0x01,0x03,0x01,0x15,0x01,0x05,0x00}, CODEC_ID_MPEG2VIDEO},
    { { 0x06,0x0e,0x2b,0x34,0x01,0x02,0x01,0x01,0x0d,0x01,0x03,0x01,0x16,0x01,0x01,0x00}, CODEC_ID_PCM_S16LE},
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00}, CODEC_ID_NONE},
};

/* partial key for header metadata */
static const uint8_t header_metadata_key[] = {0x06,0x0E,0x2B,0x34,0x02,0x53,0x01,0x01,0x0d,0x01,0x01,0x01,0x01};

/* SMPTE RP224 http://www.smpte-ra.org/mdd/index.html */
static const MXFDataDefinitionUL mxf_data_definition_uls[] = {
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x01,0x03,0x02,0x02,0x01,0x00,0x00,0x00 }, CODEC_TYPE_VIDEO },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x01,0x03,0x02,0x02,0x02,0x00,0x00,0x00 }, CODEC_TYPE_AUDIO },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x05,0x01,0x03,0x02,0x02,0x02,0x02,0x00,0x00 }, CODEC_TYPE_AUDIO },
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  CODEC_TYPE_DATA },
};


static const MXFCodecUL mxf_picture_essence_container_uls[] = {
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x02,0x0D,0x01,0x03,0x01,0x02,0x04,0x60,0x01 }, 14, CODEC_ID_MPEG2VIDEO }, /* MPEG-ES Frame wrapped */
//    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x03,0x01,0x02,0x02,0x41,0x01 }, 14,    CODEC_ID_DVVIDEO }, /* DV 625 25mbps */
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  0,       CODEC_ID_NONE },
};

static const MXFCodecUL mxf_sound_essence_container_uls[] = {
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x03,0x01,0x02,0x06,0x01,0x00 }, 14, CODEC_ID_PCM_S16LE }, /* BWF Frame wrapped */
//    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x02,0x0D,0x01,0x03,0x01,0x02,0x04,0x40,0x01 }, 14,       CODEC_ID_MP2 }, /* MPEG-ES Frame wrapped, 0x40 ??? stream id */
//    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x03,0x01,0x02,0x01,0x01,0x01 }, 14, CODEC_ID_PCM_S16LE }, /* D-10 Mapping 50Mbps PAL Extended Template */
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  0,      CODEC_ID_NONE },
};

/* SMPTE RP210 http://www.smpte-ra.org/mdd/index.html */
static const MXFLocalTagPair mxf_local_tag_batch[] = {
    // preface set
    { 0x3C0A, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x01,0x01,0x01,0x15,0x02,0x00,0x00,0x00,0x00}}, /* Instance UID */
    { 0x3B02, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x07,0x02,0x01,0x10,0x02,0x04,0x00,0x00}}, /* Last Modified Date */
    { 0x3B05, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x03,0x01,0x02,0x01,0x05,0x00,0x00,0x00}}, /* Version */
    { 0x3B06, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x06,0x04,0x00,0x00}}, /* Identifications reference */
    { 0x3B03, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x02,0x01,0x00,0x00}}, /* Content Storage reference */
    { 0x3B09, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x05,0x01,0x02,0x02,0x03,0x00,0x00,0x00,0x00}}, /* Operational Pattern UL */
    { 0x3B0A, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x05,0x01,0x02,0x02,0x10,0x02,0x01,0x00,0x00}}, /* Essence Containers UL batch */
    { 0x3B0B, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x05,0x01,0x02,0x02,0x10,0x02,0x02,0x00,0x00}}, /* DM Schemes UL batch */
    // Identification
    { 0x3C09, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x05,0x20,0x07,0x01,0x01,0x00,0x00,0x00}}, /* This Generation UID */
    { 0x3C01, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x05,0x20,0x07,0x01,0x02,0x01,0x00,0x00}}, /* Company Name */
    { 0x3C02, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x05,0x20,0x07,0x01,0x03,0x01,0x00,0x00}}, /* Product Name */
    { 0x3C04, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x05,0x20,0x07,0x01,0x04,0x00,0x00,0x00}}, /* Version String */
    { 0x3C05, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x05,0x20,0x07,0x01,0x07,0x00,0x00,0x00}}, /* Product ID */
    { 0x3C06, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x07,0x02,0x01,0x10,0x02,0x03,0x00,0x00}}, /* Modification Date */
    // Content Storage
    { 0x1901, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x05,0x01,0x00,0x00}}, /* Package strong reference batch */
    // Essence Container Data
    { 0x2701, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x06,0x01,0x00,0x00,0x00}}, /* Linked Package UID */
    { 0x3F07, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x04,0x01,0x03,0x04,0x04,0x00,0x00,0x00,0x00}}, /* BodySID */
    // Package
    { 0x4401, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x01,0x01,0x01,0x15,0x10,0x00,0x00,0x00,0x00}}, /* Package UID */
    { 0x4405, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x07,0x02,0x01,0x10,0x01,0x03,0x00,0x00}}, /* Package Creation Date */
    { 0x4404, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x07,0x02,0x01,0x10,0x02,0x05,0x00,0x00}}, /* Package Modified Date */
    { 0x4403, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x06,0x05,0x00,0x00}}, /* Tracks Strong reference array */
    // Track
    { 0x4801, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x01,0x07,0x01,0x01,0x00,0x00,0x00,0x00}}, /* Track ID */
    { 0x4804, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x01,0x03,0x00,0x00}}, /* Track Numberr */
    { 0x4B01, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x05,0x30,0x04,0x05,0x00,0x00,0x00,0x00}}, /* Edit Rate */
    { 0x4B02, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x07,0x02,0x01,0x03,0x01,0x03,0x00,0x00}}, /* Origin */
    { 0x4803, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x02,0x04,0x00,0x00}}, /* Sequence reference */
    // Sequence
    { 0x0201, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x04,0x07,0x01,0x00,0x00,0x00,0x00,0x00}}, /* Data Definition UL */
    { 0x0202, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x07,0x02,0x02,0x01,0x01,0x03,0x00,0x00}}, /* Duration */
    { 0x1001, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x06,0x09,0x00,0x00}}, /* Structural Components reference array */
    // Source Clip
    { 0x1201, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x05,0x07,0x02,0x01,0x03,0x01,0x0A,0x00,0x00}}, /* Start position */
    { 0x1101, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x03,0x01,0x00,0x00,0x00}}, /* SourcePackageID */
    { 0x1102, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x03,0x02,0x00,0x00,0x00}}, /* SourceTrackID */
};

#define PRINT_KEY(pc, s, x) dprintf(pc, "%s %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X\n", s, \
                             (x)[0], (x)[1], (x)[2], (x)[3], (x)[4], (x)[5], (x)[6], (x)[7], (x)[8], (x)[9], (x)[10], (x)[11], (x)[12], (x)[13], (x)[14], (x)[15])
static void mxf_generate_uuid(AVFormatContext *s, UID uuid)
{
    MXFContext *mxf = s->priv_data;
    int rand_num, i;

    for (i = 0; i < 16; i++) {
        rand_num = av_random(&mxf->random_state);
        rand_num = rand_num & 0x00ff;

        // the 7th byte is version according to ISO 11578
        if (i == 6) {
            rand_num &= 0x0f;
            rand_num |= 0x40;
        }

        // the 8th byte is variant for current use according to ISO 11578
        if (i == 8) {
            rand_num &= 0x3f;
            rand_num |= 0x80;
        }
        uuid[i] = rand_num;
    }
}

static void mxf_generate_umid(AVFormatContext *s, UMID umid)
{
    memcpy(umid, umid_base, 16);
    mxf_generate_uuid(s, umid + 16);
}

static int mxf_generate_reference(AVFormatContext *s, UID **refs, int ref_count)
{
    int i;
    UID *p;
    *refs = av_mallocz(ref_count * sizeof(UID));
    if (!refs)
        return -1;
    p = *refs;
    for (i = 0; i < ref_count; i++) {
        mxf_generate_uuid(s, *p);
        p += 16;
    }
    p = 0;
    return 0;
}

static int klv_encode_ber_length(ByteIOContext *pb, uint64_t len)
{
    // Determine the best BER size
    int size = 0, i;
    if (len < 128) {
        //short form
        size = 1;
        put_byte(pb, len);
        return size;
    }

    while (len) {
        len >>= 8;
        size ++;
    }

    // long form
    put_byte(pb, 0x80 + size);
    i = size;
    while(i) {
        put_byte(pb, len & 0xff);
        len >>= 8;
        i--;
    }
    return size;
}

static int mxf_write_primer_pack(AVFormatContext *s)
{
    ByteIOContext *pb = s->pb;
    const MXFLocalTagPair *local_tag_batch;
    int i,local_tag_number = 0;

    local_tag_number = sizeof(mxf_local_tag_batch) / sizeof(MXFLocalTagPair);

    put_buffer(pb, primer_pack_key, 16);
    klv_encode_ber_length(pb, local_tag_number * 18 + 8);

    put_be32(pb, local_tag_number); // local_tag num
    put_be32(pb, 18); // item size, always 18 according to the specs

    for (local_tag_batch = mxf_local_tag_batch; i < local_tag_number; local_tag_batch++, i++) {
        put_be16(pb, local_tag_batch->local_tag);
        put_buffer(pb, local_tag_batch->uid, 16);
    }
    return 0;
}

static int mxf_write_local_tag(ByteIOContext *pb, int value_size, int tag)
{
    put_be16(pb, tag);
    put_be16(pb, value_size);
    return 0;
}

static int mxf_write_reference(ByteIOContext *pb, int ref_count, UID *value)
{
    put_be32(pb, ref_count);
    put_be32(pb, 16);
    put_buffer(pb, *value, 16 * ref_count);
    return 0;
}

static int utf8len(const uint8_t *b){
    int len=0;
    int val;
    while(*b){
        GET_UTF8(val, *b++, return -1;)
        len++;
    }
    return len;
}

static void mxf_free(AVFormatContext *s)
{
    MXFContext *mxf = s->priv_data;
    int i;

    av_freep(&mxf->reference->identification);
    av_freep(&mxf->reference->content_storage);
    av_freep(&mxf->reference->package);
    av_freep(&mxf->reference->track);
    for (i = 0; i < s->nb_streams; i++) {
        av_freep(&mxf->reference->sequence[i]);
        av_freep(&mxf->reference->structural_component[i]);
    }
    av_freep(&mxf->reference->sequence);
    av_freep(&mxf->reference->structural_component);
    av_freep(&mxf->reference);
    av_freep(&mxf->track_essence_element_key);
    av_freep(&mxf->track_number_sign);
}

static const MXFDataDefinitionUL *mxf_get_data_definition_ul(const MXFDataDefinitionUL *uls, enum CodecType type)
{
    while (uls->type != CODEC_TYPE_DATA) {
        if (type == uls->type)
            break;
        uls ++;
    }
    return uls;
}

static int mxf_write_preface(AVFormatContext *s, KLVPacket *klv)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = mxf->reference;
    UID uid;
    ByteIOContext *pb = s->pb;

    klv->key[13] = 0x01;
    klv->key[14] = 0x2f;
    klv->key[15] = 0x00;

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 146);

    // write preface set uid
    mxf_generate_uuid(s, uid);
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, uid, 16);

#ifdef DEBUG
    PRINT_KEY(s, "preface uid", uid);
#endif

    // write create date as unknown
    mxf_write_local_tag(pb, 8, 0x3B02);
    put_buffer(pb, "0", 8);

    // write version
    mxf_write_local_tag(pb, 2, 0x3B05);
    put_be16(pb, 1);

    // write identification_refs
    if (mxf_generate_reference(s, &refs->identification, 1) < 0)
        return -1;
    mxf_write_local_tag(pb, 16 + 8, 0x3B06);
    mxf_write_reference(pb, 1, refs->identification);

    // write content_storage_refs
    if (mxf_generate_reference(s, &refs->content_storage, 1) < 0)
        return -1;
    mxf_write_local_tag(pb, 16, 0x3B03);
    put_buffer(pb, *refs->content_storage, 16);

    mxf_write_local_tag(pb, 16, 0x3B09);
    put_buffer(pb, op1a_ul, 16);

    // write essence_container_refs
    mxf_write_local_tag(pb, 8 + 16 * mxf->type_num, 0x3B0A);
    put_be32(pb,mxf->type_num);
    put_be32(pb,16);
    if (mxf->video_container_ul != 0)
        put_buffer(pb, mxf->video_container_ul->uid, 16);
    if (mxf->audio_container_ul != 0)
        put_buffer(pb, mxf->audio_container_ul->uid, 16);

    // write dm_scheme_refs
    mxf_write_local_tag(pb, 8, 0x3B0B);
    put_be64(pb, 0);
    return 0;
}

static int mxf_write_identification(AVFormatContext *s, KLVPacket *klv)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = mxf->reference;
    ByteIOContext *pb = s->pb;
    UID uid;
    int length, company_name_len, product_name_len, version_string_len;

    klv->key[13] = 0x01;
    klv->key[14] = 0x30;
    klv->key[15] = 0x00;

    put_buffer(pb, klv->key, 16);

    company_name_len = utf8len("FFmpeg") + 1;
    product_name_len = utf8len("OP1a Muxer") + 1;
    version_string_len = utf8len("version 0.0.1") + 1;
    length = 84 + company_name_len + product_name_len + version_string_len;

    klv_encode_ber_length(pb, length);

    // write uid
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->identification, 16);

    // write generation uid
    mxf_generate_uuid(s, uid);
    mxf_write_local_tag(pb, 16, 0x3C09);
    put_buffer(pb, uid, 16);

    mxf_write_local_tag(pb, company_name_len, 0x3C01);
    put_buffer(pb, "FFmpeg", company_name_len);

    mxf_write_local_tag(pb, product_name_len, 0x3C02);
    put_buffer(pb, "OP1a Muxer", product_name_len);

    mxf_write_local_tag(pb, version_string_len, 0x3C04);
    put_buffer(pb, "version 0.0.1", version_string_len);

    // write product uid
    mxf_generate_uuid(s, uid);
    mxf_write_local_tag(pb, 16, 0x3C05);
    put_buffer(pb, uid, 16);

    // write modified date
    mxf_write_local_tag(pb, 8, 0x3C06);
    put_buffer(pb, "0", 8);
    return 0;
}

static int mxf_write_content_storage(AVFormatContext *s, KLVPacket *klv)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = mxf->reference;
    ByteIOContext *pb = s->pb;

    klv->key[13] = 0x01;
    klv->key[14] = 0x18;
    klv->key[15] = 0x00;

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 64);

    // write uid
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->content_storage, 16);

    // write package reference
    if (mxf_generate_reference(s, &refs->package, 2) < 0)
        return -1;
    mxf_write_local_tag(pb, 16 * 2 + 8, 0x1901);
    mxf_write_reference(pb, 2, refs->package);
    return 0;
}

static int mxf_write_package(AVFormatContext *s, KLVPacket *klv, enum MXFMetadataSetType type)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = mxf->reference;
    ByteIOContext *pb = s->pb;
    UMID umid;
    UID *ref;

    klv->key[13] = 0x01;
    klv->key[14] = type == MaterialPackage ? 0x36 : 0x37;
    klv->key[15] = 0x00;

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 92 + 16 * s->nb_streams);

    // write uid
    ref = &refs->package[type == SourcePackage];
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *ref, 16);

    // write package umid
    mxf_write_local_tag(pb, 32, 0x4401);
    if (type == MaterialPackage) {
        mxf_generate_umid(s, umid);
        put_buffer(pb, umid, 32);
    } else {
        put_buffer(pb, mxf->top_src_package_uid, 32);
    }

    // write create date
    mxf_write_local_tag(pb, 8, 0x4405);
    put_buffer(pb, "0", 8);

    // write modified date
    mxf_write_local_tag(pb, 8, 0x4404);
    put_buffer(pb, "0", 8);

    // write track refs
    if (mxf_generate_reference(s, &refs->track, s->nb_streams) < 0)
        return -1;
    mxf_write_local_tag(pb, s->nb_streams * 16 + 8, 0x4403);
    mxf_write_reference(pb, s->nb_streams, refs->track);

    // every track have 1 sequence and 1 structural componet, malloc memory for the refs pointer
    refs->sequence = av_mallocz(s->nb_streams * sizeof(*refs->sequence));
    if (!refs->sequence)
        return -1;
    refs->structural_component = av_mallocz(s->nb_streams * sizeof(*refs->structural_component));
    if (refs->structural_component)
        return -1;

    // malloc memory for track number sign
    if (type == SourcePackage) {
        mxf->track_number_sign = av_mallocz(sizeof(mxf_essence_element_key)/sizeof(MXFEssenceElementKey));
        if (!mxf->track_number_sign)
            return -1;
    }

    // malloc memory for essence element key of each track
    mxf->track_essence_element_key = av_mallocz(s->nb_streams * sizeof(UID));
    if (!mxf->track_essence_element_key)
        return -1;
    return 0;
}

static int mxf_write_track(AVFormatContext *s, KLVPacket *klv, int stream_index, enum MXFMetadataSetType type)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFEssenceElementKey *element;
    int i = 0;

    klv->key[13] = 0x01;
    klv->key[14] = 0x3b;
    klv->key[15] = 0x00;

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 80);

    st = s->streams[stream_index];

    // set pts information
    if (st->codec->codec_type == CODEC_TYPE_VIDEO) {
        av_set_pts_info(st, 64, 1, st->codec->time_base.den);
    } else if (st->codec->codec_type == CODEC_TYPE_AUDIO) {
        av_set_pts_info(st, 64, 1, st->codec->sample_rate);
    }

    // write track uid
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, refs->track[stream_index], 16);

    // write track id
    mxf_write_local_tag(pb, 4, 0x4801);
    put_be32(pb, stream_index + 1);

    if (type != MaterialPackage) {
        for (element = mxf_essence_element_key; element->type != CODEC_ID_NONE; element++) {
            if (st->codec->codec_id== element->type) {
                // write track number
                mxf_write_local_tag(pb, 4, 0x4804);
                put_buffer(pb, element->uid + 12, 3);
                put_byte(pb, element->uid[15] + mxf->track_number_sign[i]);
                mxf->track_number_sign[i] ++;

                // set essence_element key
                memcpy(mxf->track_essence_element_key[stream_index], element->uid, 16);
                break;
            }
            i++;
        }
    } else {
        put_buffer(pb, "0", 4); // track number of material package is 0
    }

    mxf_write_local_tag(pb, 8, 0x4B01);
    put_be32(pb, st->time_base.num);
    put_be32(pb, st->time_base.den);

    // write origin
    mxf_write_local_tag(pb, 8, 0x4B02);
    put_be64(pb, 0);

    // write sequence refs
    if (mxf_generate_reference(s, &refs->sequence[stream_index], 1) < 0)
        return -1;
    mxf_write_local_tag(pb, 16, 0x4803);
    put_buffer(pb, *refs->sequence[stream_index], 16);
    return 0;
}

static int mxf_write_sequence(AVFormatContext *s, KLVPacket *klv, int stream_index)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFDataDefinitionUL * data_def_ul;

    klv->key[13] = 0x01;
    klv->key[14] = 0x0f;
    klv->key[15] = 0x00;

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 80);

    st = s->streams[stream_index];

    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->sequence[stream_index], 16);

    // find data define uls
    data_def_ul = mxf_get_data_definition_ul(mxf_data_definition_uls, st->codec->codec_type);
    mxf_write_local_tag(pb, 16, 0x0201);
    put_buffer(pb, data_def_ul->uid, 16);

    mxf_write_local_tag(pb, 8, 0x0202);
    put_be32(pb, st->duration);

    // write structural component
    if (mxf_generate_reference(s, &refs->structural_component[stream_index], 1) < 0)
        return -1;
    mxf_write_local_tag(pb, 16 + 8, 0x1001);
    mxf_write_reference(pb, 1, refs->structural_component[stream_index]);
    return 0;
}

static int mxf_write_structural_component(AVFormatContext *s, KLVPacket *klv, int stream_index, enum MXFMetadataSetType type)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFDataDefinitionUL * data_def_ul;

    klv->key[13] = 0x01;
    klv->key[14] = 0x11;
    klv->key[15] = 0x00;
    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 90);

    st = s->streams[stream_index];

    // write uid
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->structural_component[stream_index], 16);

    data_def_ul = mxf_get_data_definition_ul(mxf_data_definition_uls, st->codec->codec_type);
    mxf_write_local_tag(pb, 16, 0x0201);
    put_buffer(pb, data_def_ul->uid, 16);

    // write start_position
    mxf_write_local_tag(pb, 8, 0x1201);
    put_be64(pb, 0);

    // write duration
    mxf_write_local_tag(pb, 8, 0x0202);
    put_be64(pb, st->duration);

    if (type == SourcePackage) {
        // write source package uid, end of the reference
        mxf_write_local_tag(pb, 32, 0x1101);
        put_buffer(pb, "0", 32);

        // write source track id
        mxf_write_local_tag(pb, 4, 0x1102);
        put_be64(pb, 0);
    } else {
        mxf_write_local_tag(pb, 32, 0x1101);
        put_buffer(pb, mxf->top_src_package_uid, 32);

        mxf_write_local_tag(pb, 4, 0x1102);
        put_be64(pb, stream_index + 1);
    }
    return 0;
}

static int mxf_build_structural_metadata(AVFormatContext *s, KLVPacket* klv, enum MXFMetadataSetType type)
{
    int i;

    if (mxf_write_package(s, klv, type) < 0)
        return -1;

    for (i = 0;i < s->nb_streams; i++) {
        if (mxf_write_track(s, klv, i, type) < 0)
            return -1;

        if (mxf_write_sequence(s, klv, i) < 0)
            return -1;

        if (mxf_write_structural_component(s, klv, i, type) < 0)
            return -1;
    }
    return 0;
}

static int mxf_write_header_metadata_sets(AVFormatContext *s)
{
    KLVPacket klv;
    memcpy(klv.key, header_metadata_key, 13);
    if (mxf_write_preface(s, &klv) < 0)
        return -1;

    if (mxf_write_identification(s,&klv) < 0)
        return -1;

    if (mxf_write_content_storage(s, &klv) < 0)
        return -1;

    if (mxf_build_structural_metadata(s, &klv, MaterialPackage) < 0)
        return -1;

    if (mxf_build_structural_metadata(s, &klv, SourcePackage) < 0)
        return -1;
    return 0;
}

static const MXFCodecUL *mxf_get_essence_container_ul(const MXFCodecUL *uls, enum CodecID type)
{
    while (uls->id != CODEC_ID_NONE) {
        if (uls->id == type)
            break;
        uls++;
    }
    return uls;
}

static void mxf_set_essence_number(AVFormatContext *s)
{
    MXFContext *mxf = s->priv_data;
    AVStream *st;
    int i, video_type = 0, audio_type = 0;

    for (i = 0; i < s->nb_streams; i++) {
        st = s->streams[i];
        if (!video_type && st->codec->codec_type == CODEC_TYPE_VIDEO) {
            mxf->video_container_ul = mxf_get_essence_container_ul(mxf_picture_essence_container_uls, st->codec->codec_id);
            video_type++;
        }
        if (!audio_type && st->codec->codec_type == CODEC_TYPE_AUDIO) {
            mxf->audio_container_ul = mxf_get_essence_container_ul(mxf_sound_essence_container_uls, st->codec->codec_id);
            audio_type++;
        }
        if (video_type && audio_type)
            break;
    }
    mxf->type_num = video_type + audio_type;
}

static void mxf_write_partition(AVFormatContext *s, int64_t this_partition, int bodysid, const uint8_t *key)
{
    MXFContext *mxf = s->priv_data;
    ByteIOContext *pb = s->pb;

    // write klv
    put_buffer(pb, key, 16);
    klv_encode_ber_length(pb, 88 + 16 * mxf->type_num);

    // write partition value
    put_be16(pb, 1); // majorVersion
    put_be16(pb, 2); // minorVersion
    put_be32(pb, 1); // kagSize

    put_be64(pb, this_partition); // thisPartition
    put_be64(pb, 0); // previousPartition

    // set offset
    if (!this_partition)
        mxf->header_footer_partition_offset = url_ftell(pb);
    put_be64(pb, this_partition); // footerPartition,update later

    // set offset
    if (!this_partition)
        mxf->header_byte_count_offset = url_ftell(pb);
    put_be64(pb, 0); // headerByteCount, update later

    // no indexTable
    put_be64(pb, 0); // indexByteCount
    put_be32(pb, 0); // indexSID
    put_be64(pb, 0); // bodyOffset

    put_be32(pb, bodysid); // bodySID
    put_buffer(pb, op1a_ul, 16); // operational pattern

    // essence container
    put_be32(pb,mxf->type_num);
    put_be32(pb,16);
    if (mxf->video_container_ul != 0)
        put_buffer(pb, mxf->video_container_ul->uid, 16);
    if (mxf->audio_container_ul != 0)
        put_buffer(pb, mxf->audio_container_ul->uid, 16);
}

static int mux_write_header(AVFormatContext *s)
{
    MXFContext *mxf = s->priv_data;
    ByteIOContext *pb = s->pb;
    int64_t header_metadata_start;

    av_init_random(0xbeefdead, &mxf->random_state);
    // intial MXFReferenceContext
    mxf->reference = av_mallocz(sizeof(MXFReferenceContext));
    if (!mxf->reference)
        goto fail;

    // mark the header start position, for some fields update later
    mxf->header_start = url_ftell(pb);

    // calculate the numner of essence container type
    mxf_set_essence_number(s);
    mxf_write_partition(s, 0, 0, header_partition_key);

    // generate Source Package Set UMID for op1a
    // will be used by material_package->source_track->sequence->structual_component->source_package_id
    mxf_generate_umid(s, mxf->top_src_package_uid);

    // mark the start of the headermetadata and calculate metadata size
    header_metadata_start = url_ftell(s->pb);
    mxf_write_primer_pack(s);
    if (mxf_write_header_metadata_sets(s) < 0)
        goto fail;
    mxf->header_byte_count = url_ftell(s->pb) - header_metadata_start;

    put_flush_packet(pb);
    return 0;
fail:
    mxf_free(s);
    return -1;
}

static int mux_write_packet(AVFormatContext *s, AVPacket *pkt)
{
    MXFContext *mxf = s->priv_data;
    ByteIOContext *pb = s->pb;

    put_buffer(pb, mxf->track_essence_element_key[pkt->stream_index], 16); // write key
    klv_encode_ber_length(pb, pkt->size); // write length
    put_buffer(pb, pkt->data, pkt->size); // write value

    put_flush_packet(pb);
    return 0;
}

static int mxf_update_header_partition(AVFormatContext *s, int64_t footer_partition_offset)
{
    MXFContext *mxf = s->priv_data;
    ByteIOContext *pb = s->pb;

    url_fseek(pb, mxf->header_byte_count_offset, SEEK_SET);
    put_be64(pb, mxf->header_byte_count);
    put_flush_packet(pb);

    url_fseek(pb, mxf->header_footer_partition_offset, SEEK_SET);
    put_be64(pb, footer_partition_offset);
    put_flush_packet(pb);
    return 0;
}


static int mux_write_footer(AVFormatContext *s)
{
    MXFContext *mxf = s->priv_data;
    ByteIOContext *pb = s->pb;

    int64_t this_partition = url_ftell(pb) - mxf->header_start;
    mxf_write_partition(s, this_partition, 1, footer_partition_key);

    put_flush_packet(pb);

    mxf_update_header_partition(s, this_partition);
    mxf_free(s);
    return 0;
}

AVOutputFormat mxf_muxer = {
    "mxf",
    NULL_IF_CONFIG_SMALL("Material eXchange Format"),
    NULL,
    "mxf",
    sizeof(MXFContext),
    CODEC_ID_PCM_S16LE,
    CODEC_ID_MPEG2VIDEO,
    mux_write_header,
    mux_write_packet,
    mux_write_footer,
};


