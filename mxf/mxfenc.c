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
    UID **package;
    UID **track;
    UID **sequence;
    UID **structural_component;
    UID *mul_desc;
    UID **sub_desc;
} MXFReferenceContext;

typedef struct MXFContext {
    UMID top_src_package_uid;
    int64_t header_byte_count;
    int64_t header_byte_count_offset;
    int64_t header_footer_partition_offset;
    unsigned int random_state;
    MXFReferenceContext reference;
    UID *track_essence_element_key;
    int essence_container_count;
    UID *essence_container_uls;
} MXFContext;

typedef struct {
    const UID key;
    int (*write)();
    enum CodecType type;
} MXFDescriptorWriteTableEntry;

static const uint8_t umid_base[]            = { 0x06,0x0A,0x2B,0x34,0x01,0x01,0x01,0x01,0x01,0x01,0x0F,0x00,0x13,0x00,0x00,0x00 };

/* complete key */
static const uint8_t op1a_ul[]              = { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x02,0x01,0x01,0x01,0x01,0x00 };
static const uint8_t header_partition_key[] = { 0x06,0x0E,0x2B,0x34,0x02,0x05,0x01,0x01,0x0D,0x01,0x02,0x01,0x01,0x02,0x04,0x00 }; // ClosedComplete
static const uint8_t footer_partition_key[] = { 0x06,0x0E,0x2B,0x34,0x02,0x05,0x01,0x01,0x0D,0x01,0x02,0x01,0x01,0x04,0x04,0x00 }; // ClosedComplete
static const uint8_t primer_pack_key[]      = { 0x06,0x0E,0x2B,0x34,0x02,0x05,0x01,0x01,0x0D,0x01,0x02,0x01,0x01,0x05,0x01,0x00 };

/* partial key for header metadata */
static const uint8_t header_metadata_key[]  = { 0x06,0x0E,0x2B,0x34,0x02,0x53,0x01,0x01,0x0D,0x01,0x01,0x01,0x01 };

static const MXFEssenceElementKey mxf_essence_element_key[] = {
    { { 0x06,0x0E,0x2B,0x34,0x01,0x02,0x01,0x01,0x0D,0x01,0x03,0x01,0x15,0x01,0x05,0x00 }, CODEC_ID_MPEG2VIDEO},
    { { 0x06,0x0E,0x2B,0x34,0x01,0x02,0x01,0x01,0x0D,0x01,0x03,0x01,0x16,0x01,0x01,0x00 }, CODEC_ID_PCM_S16LE},
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 }, CODEC_ID_NONE},
};

static track_number_sign[sizeof(mxf_essence_element_key)/sizeof(MXFEssenceElementKey)] = { 0 };

/* SMPTE RP224 http://www.smpte-ra.org/mdd/index.html */
static const MXFDataDefinitionUL mxf_data_definition_uls[] = {
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x01,0x03,0x02,0x02,0x01,0x00,0x00,0x00 }, CODEC_TYPE_VIDEO },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x01,0x03,0x02,0x02,0x02,0x00,0x00,0x00 }, CODEC_TYPE_AUDIO },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x05,0x01,0x03,0x02,0x02,0x02,0x02,0x00,0x00 }, CODEC_TYPE_AUDIO },
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  CODEC_TYPE_DATA },
};

static const MXFCodecUL mxf_codec_uls[] = {
    /* PictureEssenceCoding */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x01,0x11,0x00 }, 14, CODEC_ID_MPEG2VIDEO }, /* MP@ML Long GoP */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x01,0x02,0x02,0x01,0x02,0x01,0x01 }, 14, CODEC_ID_MPEG2VIDEO }, /* D-10 50Mbps PAL */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x03,0x03,0x00 }, 14, CODEC_ID_MPEG2VIDEO }, /* MP@HL Long GoP */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x04,0x02,0x00 }, 14, CODEC_ID_MPEG2VIDEO }, /* 422P@HL I-Frame */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x04,0x01,0x02,0x02,0x01,0x20,0x02,0x03 }, 14,      CODEC_ID_MPEG4 }, /* XDCAM proxy_pal030926.mxf */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x01,0x02,0x02,0x02,0x01,0x02,0x00 }, 13,    CODEC_ID_DVVIDEO }, /* DV25 IEC PAL */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x07,0x04,0x01,0x02,0x02,0x03,0x01,0x01,0x00 }, 14,   CODEC_ID_JPEG2000 }, /* JPEG2000 Codestream */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x01,0x02,0x01,0x7F,0x00,0x00,0x00 }, 13,   CODEC_ID_RAWVIDEO }, /* Uncompressed */
    /* SoundEssenceCompression */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x01,0x00,0x00,0x00,0x00 }, 13,  CODEC_ID_PCM_S16LE }, /* Uncompressed */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x01,0x7F,0x00,0x00,0x00 }, 13,  CODEC_ID_PCM_S16LE },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x07,0x04,0x02,0x02,0x01,0x7E,0x00,0x00,0x00 }, 13,  CODEC_ID_PCM_S16BE }, /* From Omneon MXF file */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x04,0x04,0x02,0x02,0x02,0x03,0x01,0x01,0x00 }, 15,   CODEC_ID_PCM_ALAW }, /* XDCAM Proxy C0023S01.mxf */
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x02,0x03,0x02,0x01,0x00 }, 15,        CODEC_ID_AC3 },
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x02,0x03,0x02,0x05,0x00 }, 15,        CODEC_ID_MP2 }, /* MP2 or MP3 */
  //{ { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x04,0x02,0x02,0x02,0x03,0x02,0x1C,0x00 }, 15,    CODEC_ID_DOLBY_E }, /* Dolby-E */
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 },  0,       CODEC_ID_NONE },
};

static const uint8_t multiple_desc_ul[] = { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x03,0x0D,0x01,0x03,0x01,0x02,0x7F,0x01,0x00 };

static const MXFCodecUL mxf_essence_container_uls[] = {
    // picture essence container
    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x02,0x0D,0x01,0x03,0x01,0x02,0x04,0x60,0x01 }, 14, CODEC_ID_MPEG2VIDEO }, /* MPEG-ES Frame wrapped */
//    { { 0x06,0x0E,0x2B,0x34,0x04,0x01,0x01,0x01,0x0D,0x01,0x03,0x01,0x02,0x02,0x41,0x01 }, 14,    CODEC_ID_DVVIDEO }, /* DV 625 25mbps */
    // audio essence conatiner uls
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
    { 0x4701, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x02,0x03,0x00,0x00}}, /* Descriptor */
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
    // file descriptor
    { 0x3F01, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x04,0x06,0x01,0x01,0x04,0x06,0x0B,0x00,0x00}}, /* sub descriptor uid*/
    { 0x3006, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x05,0x06,0x01,0x01,0x03,0x05,0x00,0x00,0x00}}, /* Linked Track ID */
    { 0x3004, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x06,0x01,0x01,0x04,0x01,0x02,0x00,0x00}}, /* essence container ul */
    // generic picture eseence descriptor
    { 0x3203, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x01,0x04,0x01,0x05,0x02,0x02,0x00,0x00,0x00}}, /* stored width */
    { 0x3202, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x01,0x04,0x01,0x05,0x02,0x01,0x00,0x00,0x00}}, /* stored heigth */
    { 0x320E, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x01,0x04,0x01,0x01,0x01,0x01,0x00,0x00,0x00}}, /* aspect ratio*/
    { 0x3201, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x04,0x01,0x06,0x01,0x00,0x00,0x00,0x00}}, /* picture essence coding*/
    // generic picture sound essence descriptor
    { 0x3D03, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x05,0x04,0x02,0x03,0x01,0x01,0x01,0x00,0x00}}, /* audio sampling rate */
    { 0x3D07, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x05,0x04,0x02,0x01,0x01,0x04,0x00,0x00,0x00}}, /* channel count */
    { 0x3D01, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x04,0x04,0x02,0x03,0x03,0x04,0x00,0x00,0x00}}, /* quantization bits */
    { 0x3D06, {0x06,0x0E,0x2B,0x34,0x01,0x01,0x01,0x02,0x04,0x01,0x06,0x01,0x00,0x00,0x00,0x00}}, /* sound essence compression */
};

#define PRINT_KEY(pc, s, x) dprintf(pc, "%s %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X %02X\n", s, \
                             (x)[0], (x)[1], (x)[2], (x)[3], (x)[4], (x)[5], (x)[6], (x)[7], (x)[8], (x)[9], (x)[10], (x)[11], (x)[12], (x)[13], (x)[14], (x)[15])
static void mxf_generate_uuid(AVFormatContext *s, UID uuid)
{
    MXFContext *mxf = s->priv_data;
    int i;

    for (i = 0; i < 16; i++) {
        mxf->random_state= mxf->random_state*1664525+10139042;
        uuid[i]= mxf->random_state>>24;
    }
    // the 7th byte is version according to ISO 11578
    uuid[6] &= 0x0f;
    uuid[6] |= 0x40;

    // the 8th byte is variant for current use according to ISO 11578
    uuid[8] &= 0x3f;
    uuid[8] |= 0x80;
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
    if (!*refs)
        return AVERROR(ENOMEM);
    p = *refs;
    for (i = 0; i < ref_count; i++) {
        mxf_generate_uuid(s, *p);
        p ++;
    }
    return 0;
}

static int klv_encode_ber_length(ByteIOContext *pb, uint64_t len)
{
    // Determine the best BER size
    int size = 0;
    uint64_t tmp = len;
    if (len < 128) {
        //short form
        put_byte(pb, len);
        return 1;
    }

    while (tmp) {
        tmp >>= 8;
        size ++;
    }

    // long form
    put_byte(pb, 0x80 + size);
    while(size) {
        size --;
        put_byte(pb, len >> 8 * size & 0xff);
    }
    return 0;
}

static const MXFCodecUL *mxf_get_essence_container_ul(enum CodecID type)
{
    const MXFCodecUL *uls = mxf_essence_container_uls;
    while (uls->id != CODEC_ID_NONE) {
        if (uls->id == type)
            break;
        uls++;
    }
    return uls;
}

static int mxf_write_primer_pack(AVFormatContext *s)
{
    ByteIOContext *pb = s->pb;
    const MXFLocalTagPair *local_tag_batch;
    int local_tag_number, i = 0;

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

static void mxf_write_local_tag(ByteIOContext *pb, int value_size, int tag)
{
    put_be16(pb, tag);
    put_be16(pb, value_size);
}

static void mxf_write_reference(ByteIOContext *pb, int ref_count, UID value)
{
    put_be32(pb, ref_count);
    put_be32(pb, 16);
    put_buffer(pb, value, sizeof(UID) * ref_count);
}

static void mxf_free(AVFormatContext *s)
{
    MXFContext *mxf = s->priv_data;
    int i;

    av_freep(&mxf->reference.identification);
    av_freep(&mxf->reference.content_storage);
    av_freep(&mxf->reference.package);
    av_freep(&mxf->reference.track);
    for (i = 0; i < s->nb_streams; i++) {
        av_freep(&mxf->reference.sequence[i]);
        av_freep(&mxf->reference.structural_component[i]);
    }
    av_freep(&mxf->reference.sequence);
    av_freep(&mxf->reference.structural_component);
    av_freep(&mxf->reference.track);
    av_freep(&mxf->reference.sub_desc);
    av_freep(&mxf->reference.mul_desc);
    av_freep(&mxf->reference.package);
    av_freep(&mxf->track_essence_element_key);
    av_freep(&mxf->essence_container_uls);
}

static const MXFDataDefinitionUL *mxf_get_data_definition_ul(enum CodecType type)
{
    const MXFDataDefinitionUL *uls = mxf_data_definition_uls;
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
    MXFReferenceContext *refs = &mxf->reference;
    UID uid;
    ByteIOContext *pb = s->pb;

    AV_WB24(klv->key + 13, 0x012f00);

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 146);

    // write preface set uid
    mxf_generate_uuid(s, uid);
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, uid, 16);

#ifdef DEBUG
    PRINT_KEY(s, "preface key", klv->key);
    PRINT_KEY(s, "preface uid", uid);
#endif

    // write create date as unknown
    mxf_write_local_tag(pb, 8, 0x3B02);
    put_be64(pb, 0);

    // write version
    mxf_write_local_tag(pb, 2, 0x3B05);
    put_be16(pb, 1);

    // write identification_refs
    if (mxf_generate_reference(s, &refs->identification, 1) < 0)
        return -1;
    mxf_write_local_tag(pb, 16 + 8, 0x3B06);
    mxf_write_reference(pb, 1, *refs->identification);

    // write content_storage_refs
    if (mxf_generate_reference(s, &refs->content_storage, 1) < 0)
        return -1;
    mxf_write_local_tag(pb, 16, 0x3B03);
    put_buffer(pb, *refs->content_storage, 16);

    mxf_write_local_tag(pb, 16, 0x3B09);
    put_buffer(pb, op1a_ul, 16);

    // write essence_container_refs
    mxf_write_local_tag(pb, 8 + 16 * mxf->essence_container_count, 0x3B0A);
    mxf_write_reference(pb, mxf->essence_container_count, *mxf->essence_container_uls);

    // write dm_scheme_refs
    mxf_write_local_tag(pb, 8, 0x3B0B);
    put_be64(pb, 0);
    return 0;
}

static int mxf_write_identification(AVFormatContext *s, KLVPacket *klv)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;
    UID uid;
    int length, company_name_len, product_name_len, version_string_len;

    AV_WB24(klv->key + 13, 0x013000);

    put_buffer(pb, klv->key, 16);

    company_name_len = sizeof("FFmpeg");
    product_name_len = sizeof("OP1a Muxer");

    length = 80 + company_name_len + product_name_len;
    if (!(s->streams[0]->codec->flags & CODEC_FLAG_BITEXACT)) {
        version_string_len = sizeof(LIBAVFORMAT_IDENT);
        length += 4 + version_string_len;
    }
    klv_encode_ber_length(pb, length);

    // write uid
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->identification, 16);
#ifdef DEBUG
    PRINT_KEY(s, "identification key", klv->key);
    PRINT_KEY(s, "identification uid", *refs->identification);
#endif
    // write generation uid
    mxf_generate_uuid(s, uid);
    mxf_write_local_tag(pb, 16, 0x3C09);
    put_buffer(pb, uid, 16);

    mxf_write_local_tag(pb, company_name_len, 0x3C01);
    put_buffer(pb, "FFmpeg", company_name_len);

    mxf_write_local_tag(pb, product_name_len, 0x3C02);
    put_buffer(pb, "OP1a Muxer", product_name_len);

    if (!(s->streams[0]->codec->flags & CODEC_FLAG_BITEXACT)) {
        mxf_write_local_tag(pb, version_string_len, 0x3C04);
        put_buffer(pb, LIBAVFORMAT_IDENT, version_string_len);
    }

    // write product uid
    mxf_generate_uuid(s, uid);
    mxf_write_local_tag(pb, 16, 0x3C05);
    put_buffer(pb, uid, 16);

    // write modified date
    mxf_write_local_tag(pb, 8, 0x3C06);
    put_be64(pb, 0);
    return 0;
}

static int mxf_write_content_storage(AVFormatContext *s, KLVPacket *klv)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;

    AV_WB24(klv->key + 13, 0x011800);

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 64);

    // write uid
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->content_storage, 16);
#ifdef DEBUG
    PRINT_KEY(s, "content storage key", klv->key);
    PRINT_KEY(s, "content storage uid", *refs->content_storage);
#endif
    // write package reference
    refs->package= av_mallocz(s->nb_streams * sizeof(*refs->package));
    if (!refs->package)
        return AVERROR(ENOMEM);
    if (mxf_generate_reference(s, refs->package, 2) < 0)
        return -1;
    mxf_write_local_tag(pb, 16 * 2 + 8, 0x1901);
    mxf_write_reference(pb, 2, **refs->package);
    return 0;
}

static int mxf_write_package(AVFormatContext *s, KLVPacket *klv, enum MXFMetadataSetType type)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;
    UMID umid;
    int i;

    klv->key[13] = 0x01;
    klv->key[14] = type == MaterialPackage ? 0x36 : 0x37;
    klv->key[15] = 0x00;

    put_buffer(pb, klv->key, 16);
    if (type == MaterialPackage)
        klv_encode_ber_length(pb, 92 + 16 * s->nb_streams);
    else
        klv_encode_ber_length(pb, 112 + 16 * s->nb_streams); // 20 bytes length for descriptor reference

    // write uid
    i = type == MaterialPackage ? 0 : 1;
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, (*refs->package)[i], 16);
#ifdef DEBUG
    av_log(s,AV_LOG_DEBUG, "package type:%d\n", type);
    PRINT_KEY(s, "package", klv->key);
    PRINT_KEY(s, "package uid", (*refs->package)[i]);
    PRINT_KEY(s, "package umid first part", umid);
    PRINT_KEY(s, "package umid second part", umid + 16);
#endif

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
    put_be64(pb, 0);

    // write modified date
    mxf_write_local_tag(pb, 8, 0x4404);
    put_be64(pb, 0);

    // write track refs
    refs->track = av_mallocz(s->nb_streams * sizeof(*refs->track));
    if (!refs->track)
        return AVERROR(ENOMEM);
    if (mxf_generate_reference(s, refs->track, s->nb_streams) < 0)
        return -1;
    mxf_write_local_tag(pb, s->nb_streams * 16 + 8, 0x4403);
    mxf_write_reference(pb, s->nb_streams, **refs->track);

    // every track have 1 sequence and 1 structural componet, malloc memory for the refs pointer
    refs->sequence = av_mallocz(s->nb_streams * sizeof(*refs->sequence));
    if (!refs->sequence)
        return AVERROR(ENOMEM);
    refs->structural_component = av_mallocz(s->nb_streams * sizeof(*refs->structural_component));
    if (!refs->structural_component)
        return AVERROR(ENOMEM);

    // malloc memory for track number sign
    if (type == SourcePackage) {
        // write multiple descriptor reference
        if (mxf_generate_reference(s, &refs->mul_desc, 1) < 0)
            return -1;
        mxf_write_local_tag(pb, 16, 0x4701);
        put_buffer(pb, *refs->mul_desc, 16);
    }

    // malloc memory for essence element key of each track
    mxf->track_essence_element_key = av_mallocz(s->nb_streams * sizeof(UID));
    if (!mxf->track_essence_element_key)
        return AVERROR(ENOMEM);
    return 0;
}

static int mxf_write_track(AVFormatContext *s, KLVPacket *klv, int stream_index, enum MXFMetadataSetType type)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFEssenceElementKey *element;
    int i = 0;

    AV_WB24(klv->key + 13, 0x013b00);

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
    put_buffer(pb, (*refs->track)[stream_index], 16);
#ifdef DEBUG
    PRINT_KEY(s, "track key", klv->key);
    PRINT_KEY(s, "track uid", (*refs->track)[stream_index]);
#endif
    // write track id
    mxf_write_local_tag(pb, 4, 0x4801);
    put_be32(pb, stream_index);

    mxf_write_local_tag(pb, 4, 0x4804);
    if (type != MaterialPackage) {
        for (element = mxf_essence_element_key; element->type != CODEC_ID_NONE; element++) {
            if (st->codec->codec_id== element->type) {
                // write track number
                put_buffer(pb, element->uid + 12, 3);
                put_byte(pb, element->uid[15] + track_number_sign[i]);
                track_number_sign[i] ++;

                // set essence_element key
                memcpy(mxf->track_essence_element_key[stream_index], element->uid, 16);
                break;
            }
            i++;
        }
    } else {
        put_be32(pb, 0); // track number of material package is 0
    }

    mxf_write_local_tag(pb, 8, 0x4B01);
    put_be32(pb, st->time_base.den);
    put_be32(pb, st->time_base.num);

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
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFDataDefinitionUL * data_def_ul;

    AV_WB24(klv->key + 13, 0x010f00);

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 80);

    st = s->streams[stream_index];

    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->sequence[stream_index], 16);

#ifdef DEBUG
    PRINT_KEY(s, "sequence key", klv->key);
    PRINT_KEY(s, "sequence uid", *refs->sequence[stream_index]);
#endif
    // find data define uls
    data_def_ul = mxf_get_data_definition_ul(st->codec->codec_type);
    mxf_write_local_tag(pb, 16, 0x0201);
    put_buffer(pb, data_def_ul->uid, 16);

    mxf_write_local_tag(pb, 8, 0x0202);
    put_be64(pb, st->duration);

    // write structural component
    if (mxf_generate_reference(s, &refs->structural_component[stream_index], 1) < 0)
        return -1;
    mxf_write_local_tag(pb, 16 + 8, 0x1001);
    mxf_write_reference(pb, 1, *refs->structural_component[stream_index]);
    return 0;
}

static int mxf_write_structural_component(AVFormatContext *s, KLVPacket *klv, int stream_index, enum MXFMetadataSetType type)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFDataDefinitionUL * data_def_ul;
    int i;

    AV_WB24(klv->key + 13, 0x011100);

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 108);

    st = s->streams[stream_index];

    // write uid
    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->structural_component[stream_index], 16);

#ifdef DEBUG
    PRINT_KEY(s, "structural component key", klv->key);
    PRINT_KEY(s, "structural component uid", *refs->structural_component[stream_index]);
#endif
    data_def_ul = mxf_get_data_definition_ul(st->codec->codec_type);
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
        for (i = 0; i < 4; i++) {
            put_be64(pb, 0);
        }

        // write source track id
        mxf_write_local_tag(pb, 4, 0x1102);
        put_be32(pb, 0);
    } else {
        mxf_write_local_tag(pb, 32, 0x1101);
        put_buffer(pb, mxf->top_src_package_uid, 32);

        mxf_write_local_tag(pb, 4, 0x1102);
        put_be32(pb, stream_index);
    }
    return 0;
}

static int mxf_write_multi_descriptor(AVFormatContext *s, KLVPacket *klv)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;

    AV_WB24(klv->key + 13, 0x014400);

    put_buffer(pb, klv->key, 16);
    klv_encode_ber_length(pb, 64 + 16 * s->nb_streams);

    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, *refs->mul_desc, 16);
#ifdef DEBUG
    PRINT_KEY(s, "multi_desc uid", *refs->mul_desc);
#endif
    // write sample rate
    // SMPTE377M D.1 says this field is necessary,
    // but mxf.c actually do not read the field,so we set 0 as default.
    mxf_write_local_tag(pb, 8, 0x3001);
    put_be64(pb, 0);

    // write essence container ul
    mxf_write_local_tag(pb, 16, 0x3004);
    put_buffer(pb, multiple_desc_ul, 16);

    // write sub descriptor refs
    refs->sub_desc= av_mallocz(s->nb_streams * sizeof(*refs->sub_desc));
    if (!refs->sub_desc)
        return AVERROR(ENOMEM);
    if (mxf_generate_reference(s, refs->sub_desc, s->nb_streams) < 0)
        return -1;
    mxf_write_local_tag(pb, s->nb_streams * 16 + 8, 0x3F01);
    mxf_write_reference(pb, s->nb_streams, **refs->sub_desc);
    return 0;
}

static int mxf_write_mpeg_video_desc(AVFormatContext *s, const MXFDescriptorWriteTableEntry *desc_tbl, int stream_index)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFCodecUL *codec_ul = NULL;

    st = s->streams[stream_index];

    put_buffer(pb, desc_tbl->key, 16);
    klv_encode_ber_length(pb, 96);

    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, (*refs->sub_desc)[stream_index], 16);

    mxf_write_local_tag(pb, 4, 0x3006);
    put_be32(pb, stream_index);
#ifdef DEBUG
    PRINT_KEY(s, "mpeg2video uid", (*refs->track)[stream_index]);
    av_log(s, AV_LOG_DEBUG, "linked track ID:%d\n", stream_index);
#endif

    codec_ul = mxf_get_essence_container_ul(st->codec->codec_id);
    mxf_write_local_tag(pb, 16, 0x3004);
    put_buffer(pb, codec_ul->uid, 16);

    mxf_write_local_tag(pb, 4, 0x3203);
    put_be32(pb, st->codec->width);

    mxf_write_local_tag(pb, 4, 0x3202);
    put_be32(pb, st->codec->height);

    mxf_write_local_tag(pb, 8, 0x320E);
    put_be32(pb, st->time_base.den);
    put_be32(pb, st->time_base.num);

    // tmp write, will modified later
    mxf_write_local_tag(pb, 16, 0x3201);
    put_buffer(pb, mxf_codec_uls->uid, 16);
    return 0;
}

static int mxf_write_wav_desc(AVFormatContext *s, const MXFDescriptorWriteTableEntry *desc_tbl, int stream_index)
{
    MXFContext *mxf = s->priv_data;
    MXFReferenceContext *refs = &mxf->reference;
    ByteIOContext *pb = s->pb;
    AVStream *st;
    const MXFCodecUL *codec_ul = NULL;

    st = s->streams[stream_index];

    put_buffer(pb, desc_tbl->key, 16);
    klv_encode_ber_length(pb, 96);

    mxf_write_local_tag(pb, 16, 0x3C0A);
    put_buffer(pb, (*refs->sub_desc)[stream_index], 16);

    mxf_write_local_tag(pb, 4, 0x3006);
    put_be32(pb, stream_index);
#ifdef DEBUG
    PRINT_KEY(s, "wav desc uid", (*refs->track)[stream_index]);
#endif
    codec_ul = mxf_get_essence_container_ul(st->codec->codec_id);
    mxf_write_local_tag(pb, 16, 0x3004);
    put_buffer(pb, codec_ul->uid, 16);

    // write audio sampling rate
    mxf_write_local_tag(pb, 8, 0x3D03);
    put_be32(pb, st->codec->sample_rate);
    put_be32(pb, 1);

    mxf_write_local_tag(pb, 4, 0x3D07);
    put_be32(pb, st->codec->channels);

    mxf_write_local_tag(pb, 4, 0x3D01);
    put_be32(pb, st->codec->bits_per_sample);

    // tmp write, will modified later
    mxf_write_local_tag(pb, 16, 0x3201);
    put_buffer(pb, (mxf_codec_uls + 8) ->uid, 16);
    return 0;
}

static const MXFDescriptorWriteTableEntry mxf_descriptor_read_table[] = {
    { { 0x06,0x0E,0x2B,0x34,0x02,0x53,0x01,0x01,0x0d,0x01,0x01,0x01,0x01,0x01,0x51,0x00 }, mxf_write_mpeg_video_desc, CODEC_ID_MPEG2VIDEO},
    { { 0x06,0x0E,0x2B,0x34,0x02,0x53,0x01,0x01,0x0d,0x01,0x01,0x01,0x01,0x01,0x48,0x00 }, mxf_write_wav_desc, CODEC_ID_PCM_S16LE},
    { { 0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00,0x00 }, NULL, CODEC_ID_NONE},
};

static int mxf_build_structural_metadata(AVFormatContext *s, KLVPacket* klv, enum MXFMetadataSetType type)
{
    int i;
    const MXFDescriptorWriteTableEntry *desc = NULL;

    if (mxf_write_package(s, klv, type) < 0)
        return -1;
    if (type == SourcePackage) {
        if (mxf_write_multi_descriptor(s, klv) < 0)
            return -1;
    }

    for (i = 0;i < s->nb_streams; i++) {
        if (mxf_write_track(s, klv, i, type) < 0)
            return -1;

        if (mxf_write_sequence(s, klv, i) < 0)
            return -1;

        if (mxf_write_structural_component(s, klv, i, type) < 0)
            return -1;

        if (type == SourcePackage) {
            for (desc = mxf_descriptor_read_table; desc->write; desc++) {
                if (s->streams[i]->codec->codec_id == desc->type) {
                    if (desc->write(s, desc, i) < 0) {
                        av_log(s, AV_LOG_ERROR, "error writing descriptor\n");
                        return -1;
                    }
                    break;
                }
            }
        }
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

static int mxf_add_essence_container_ul(MXFContext *mxf, const MXFCodecUL *codec_ul)
{
    int i;
    mxf->essence_container_uls = av_realloc(mxf->essence_container_uls, (mxf->essence_container_count + 1) * 16);
    if (!mxf->essence_container_uls)
        return -1;
    memcpy(mxf->essence_container_uls[mxf->essence_container_count], codec_ul->uid, 16);
    mxf->essence_container_count++;
    return mxf->essence_container_count;
}

static int mxf_build_essence_container_refs(AVFormatContext *s)
{
    MXFContext *mxf = s->priv_data;
    AVStream *st;
    int i;
    const MXFCodecUL *codec_ul = NULL;

    for (codec_ul = mxf_essence_container_uls; codec_ul->id; codec_ul++) {
        for (i = 0; i < s->nb_streams; i++) {
            st = s->streams[i];
            if (st->codec->codec_id == codec_ul->id) {
                if (mxf_add_essence_container_ul(mxf, codec_ul) < 0 )
                    return -1;
                break;
            }
        }
    }
    return 0;
}

static void mxf_write_partition(AVFormatContext *s, int64_t this_partition, int bodysid, const uint8_t *key)
{
    MXFContext *mxf = s->priv_data;
    ByteIOContext *pb = s->pb;
#ifdef DEBUG
    int i;
#endif
    // write klv
    put_buffer(pb, key, 16);
    klv_encode_ber_length(pb, 88 + 16 * mxf->essence_container_count);

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
    mxf_write_reference(pb, mxf->essence_container_count, *mxf->essence_container_uls);
#ifdef DEBUG
    av_log(s,AV_LOG_DEBUG, "essence container count:%d\n", mxf->essence_container_count);
    for (i = 0; i < mxf->essence_container_count; i++)
        PRINT_KEY(s, "essence container ul:\n", mxf->essence_container_uls[i]);

#endif
}

static int mux_write_header(AVFormatContext *s)
{
    MXFContext *mxf = s->priv_data;
    ByteIOContext *pb = s->pb;
    int64_t header_metadata_start;

    // calculate the numner of essence container type
    mxf_build_essence_container_refs(s);
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

    int64_t this_partition = url_ftell(pb);
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


