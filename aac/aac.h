
#define MAX_CHANNELS 64
#define MAX_TAGID 16

// Audio Object Types
enum {
    AOT_NULL = 0x0,
    AOT_AAC_MAIN,
    AOT_AAC_LC,
    AOT_AAC_SSR,
    AOT_AAC_LTP,
    AOT_SBR,
    AOT_AAC_SCALABLE,
    AOT_TWINVQ,
    AOT_CELP,
    AOT_HVXC
};

// IDs for raw_data_block
enum {
    ID_SCE = 0x0,
    ID_CPE,
    ID_CCE,
    ID_LFE,
    ID_DSE,
    ID_PCE,
    ID_FIL,
    ID_END
};

// window sequences
enum {
    ONLY_LONG_SEQUENCE = 0,
    LONG_START_SEQUENCE,
    EIGHT_SHORT_SEQUENCE,
    LONG_STOP_SEQUENCE
};

// special codebooks
#define ZERO_HCB 0
#define FIRST_PAIR_HCB 5
#define ESC_HCB 11
#define NOISE_HCB 13
#define INTENSITY_HCB2 14
#define INTENSITY_HCB 15
#define ESC_FLAG 16

//tns
#define TNS_MAX_ORDER 20

// dithering
#define N 624
#define M 397
#define MATRIX_A    0x9908b0df
#define UPPER_MASK  0x80000000
#define LOWER_MASK  0x7fffffff

typedef struct {
    uint32_t mt[N];
    int      mti;
} dither_state;

typedef struct {
    int present;

    int num_channels;

    int num_front;
    int front_cpe;
    int front_tag[MAX_TAGID];

    int num_side;
    int side_cpe;
    int side_tag[MAX_TAGID];

    int num_back;
    int back_cpe;
    int back_tag[MAX_TAGID];

    int num_lfe;
    int lfe_tag[MAX_TAGID];

    int num_assoc_data;
    int assoc_data_tag[MAX_TAGID];

    int num_cc;
    int cc_ind_sw;
    int cc_tag[MAX_TAGID];

    int mono_mixdown;
    int stereo_mixdown;
    int matrix_mixdown;
    int pseudo_surround;
} program_config_struct;

typedef struct {
    int intensity_present;
    int noise_present;

    int max_sfb;
    int window_sequence;
    int window_shape;
    int window_shape_prev;
    int predictor;
    int num_window_groups;
    uint8_t grouping;
    uint8_t group_len[8];
    // calculated
    const uint16_t *swb_offset;
    int num_swb;
    int num_windows;
    int tns_max_bands;
} ics_struct;

typedef struct {
    int present;
    int n_filt[8];
    int length[8][4];
    int direction[8][4];
    int order[8][4];
    float *tmp2_map;
    int coef[8][4][TNS_MAX_ORDER];
} tns_struct;

typedef struct {
    int present;
    int mask[8][64];
} ms_struct;

typedef struct {
    int present;
    int num_pulse;
    int start;
    int offset[4];
    int amp[4];
} pulse_struct;

typedef struct {
    int ind_sw;
    int domain;

    int num_coupled;
    int is_cpe[9];
    int tag_select[9];
    int l[9];
    int r[9];

    float gain[18][8][64];
} coupling_struct;

// individual channel element
typedef struct {
    int global_gain;
    ics_struct ics;
    tns_struct tns;
    int cb[8][64];   // codebooks
    float sf[8][64];
    DECLARE_ALIGNED_16(float, coeffs[1024]);
    DECLARE_ALIGNED_16(float, saved[1024]);
    DECLARE_ALIGNED_16(float, ret[1024]);
} sce_struct;

// channel element
typedef struct {
    ms_struct ms;
    sce_struct *ch[2];
} cpe_struct;

typedef struct {
    coupling_struct coup;
    sce_struct *ch;
} cc_struct;

typedef struct {
    // objects
    AVCodecContext * avccontext;
    GetBitContext gb;
    VLC mainvlc;
    VLC books[11];
    dither_state dither;

    // main config
    int audioObjectType;
    int ext_audioObjectType;
    int sbr_present;
    int sampling_index;
    int ext_sampling_index;
    int sample_rate;
    int ext_sample_rate;
    int channels;
    int frame_length;

    // decoder param
    program_config_struct pcs;
    sce_struct * che_sce[MAX_TAGID];
    cpe_struct * che_cpe[MAX_TAGID];
    sce_struct * che_lfe[MAX_TAGID];
    cc_struct * che_cc[MAX_TAGID];

    DECLARE_ALIGNED_16(float, buf_mdct[2048]);
    int is_saved;

    //cashes
    const uint16_t *swb_offset_1024;
    const uint16_t *swb_offset_128;
    int num_swb_1024;
    int num_swb_128;
    int tns_max_bands_1024;
    int tns_max_bands_128;

    // tables
    DECLARE_ALIGNED_16(float, kbd_long_1024[1024]);
    DECLARE_ALIGNED_16(float, kbd_short_128[128]);
    DECLARE_ALIGNED_16(float, sine_long_1024[1024]);
    DECLARE_ALIGNED_16(float, sine_short_128[128]);
    DECLARE_ALIGNED_16(float, pow2sf_tab[256]);
    DECLARE_ALIGNED_16(float, intensity_tab[256]);
    DECLARE_ALIGNED_16(float, ivquant_tab[256]);

    MDCTContext mdct;
    MDCTContext mdct_small;
    int * vq[11];

    // statistics
    int num_frame;
} aac_context_t;

// sampling table
static const int sampling_table[] = { 96000, 88200, 64000, 48000, 44100, 32000, 24000, 22050, 16000, 12000, 11025, 8000, 7350 };

// scalefactor bands
static const uint16_t swb_offset_1024_96[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56,
    64, 72, 80, 88, 96, 108, 120, 132, 144, 156, 172, 188, 212, 240,
    276, 320, 384, 448, 512, 576, 640, 704, 768, 832, 896, 960, 1024
};

static const uint16_t swb_offset_128_96[] = {
    0, 4, 8, 12, 16, 20, 24, 32, 40, 48, 64, 92, 128
};

static const uint16_t swb_offset_1024_64[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56,
    64, 72, 80, 88, 100, 112, 124, 140, 156, 172, 192, 216, 240, 268,
    304, 344, 384, 424, 464, 504, 544, 584, 624, 664, 704, 744, 784, 824,
    864, 904, 944, 984, 1024
};

static const uint16_t swb_offset_128_64[] = {
    0, 4, 8, 12, 16, 20, 24, 32, 40, 48, 64, 92, 128
};


static const uint16_t swb_offset_1024_48[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 56, 64, 72,
    80, 88, 96, 108, 120, 132, 144, 160, 176, 196, 216, 240, 264, 292,
    320, 352, 384, 416, 448, 480, 512, 544, 576, 608, 640, 672, 704, 736,
    768, 800, 832, 864, 896, 928, 1024
};

static const uint16_t swb_offset_128_48[] = {
    0, 4, 8, 12, 16, 20, 28, 36, 44, 56, 68, 80, 96, 112, 128
};

static const uint16_t swb_offset_1024_32[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 48, 56, 64, 72,
    80, 88, 96, 108, 120, 132, 144, 160, 176, 196, 216, 240, 264, 292,
    320, 352, 384, 416, 448, 480, 512, 544, 576, 608, 640, 672, 704, 736,
    768, 800, 832, 864, 896, 928, 960, 992, 1024
};


static const uint16_t swb_offset_1024_24[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 52, 60, 68,
    76, 84, 92, 100, 108, 116, 124, 136, 148, 160, 172, 188, 204, 220,
    240, 260, 284, 308, 336, 364, 396, 432, 468, 508, 552, 600, 652, 704,
    768, 832, 896, 960, 1024
};

static const uint16_t swb_offset_128_24[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 36, 44, 52, 64, 76, 92, 108, 128
};

static const uint16_t swb_offset_1024_16[] = {
    0, 8, 16, 24, 32, 40, 48, 56, 64, 72, 80, 88, 100, 112, 124,
    136, 148, 160, 172, 184, 196, 212, 228, 244, 260, 280, 300, 320, 344,
    368, 396, 424, 456, 492, 532, 572, 616, 664, 716, 772, 832, 896, 960, 1024
};

static const uint16_t swb_offset_128_16[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 32, 40, 48, 60, 72, 88, 108, 128
};

static const uint16_t swb_offset_1024_8[] = {
    0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 120, 132, 144, 156, 172,
    188, 204, 220, 236, 252, 268, 288, 308, 328, 348, 372, 396, 420, 448,
    476, 508, 544, 580, 620, 664, 712, 764, 820, 880, 944, 1024
};

static const uint16_t swb_offset_128_8[] = {
    0, 4, 8, 12, 16, 20, 24, 28, 36, 44, 52, 60, 72, 88, 108, 128
};

static const uint16_t *swb_offset_1024[] = {
    swb_offset_1024_96, swb_offset_1024_96, swb_offset_1024_64,
    swb_offset_1024_48, swb_offset_1024_48, swb_offset_1024_32,
    swb_offset_1024_24, swb_offset_1024_24, swb_offset_1024_16,
    swb_offset_1024_16, swb_offset_1024_16, swb_offset_1024_8
};

static const uint8_t num_swb_1024[] = {
    41, 41, 47, 49, 49, 51, 47, 47, 43, 43, 43, 40
};

static const uint16_t *swb_offset_128[] = {
    swb_offset_128_96, swb_offset_128_96, swb_offset_128_64,
    swb_offset_128_48, swb_offset_128_48, swb_offset_128_48,
    swb_offset_128_24, swb_offset_128_24, swb_offset_128_16,
    swb_offset_128_16, swb_offset_128_16, swb_offset_128_8
};

static const uint8_t num_swb_128[] = {
    12, 12, 12, 14, 14, 14, 15, 15, 15, 15, 15, 15
};

// TNS tables
static const uint8_t tns_max_bands_1024[] = {
    31, 31, 34, 40, 42, 51, 46, 46, 42, 42, 42, 39
};

static const uint8_t tns_max_bands_128[] = {
    9, 9, 10, 14, 14, 14, 14, 14, 14, 14, 14, 14
};

static float tns_tmp2_map_1_3[TNS_MAX_ORDER] = {
    0.00000000,    0.43388373, -0.64278758, -0.34202015,
    0.97492790, 0.78183150, -0.64278758, -0.34202015,
    -0.43388373, -0.78183150, -0.64278758, -0.34202015,
    -0.78183150, -0.43388373, -0.64278758, -0.34202015,
    0.78183150, 0.97492790, -0.64278758, -0.34202015
};

static float tns_tmp2_map_0_3[TNS_MAX_ORDER] = {
    0.00000000, 0.43388373, 0.78183150, 0.97492790,
    -0.98480773, -0.86602539, -0.64278758, -0.34202015,
    -0.43388373, -0.78183150, -0.97492790, -0.97492790,
    -0.98480773, -0.86602539, -0.64278758, -0.34202015,
    0.78183150, 0.97492790, 0.97492790, 0.78183150
};

static float tns_tmp2_map_1_4[TNS_MAX_ORDER] = {
    0.00000000, 0.20791170, 0.40673664, 0.58778524,
    -0.67369562, -0.52643216, -0.36124167, -0.18374951,
    0.99452192, 0.95105648, 0.86602539, 0.74314481,
    -0.67369562, -0.52643216, -0.36124167, -0.18374951,
    -0.20791176, -0.40673670, -0.58778530, -0.74314487
};

static float tns_tmp2_map_0_4[TNS_MAX_ORDER] = {
    0.00000000, 0.20791170, 0.40673664, 0.58778524,
    0.74314481, 0.86602539, 0.95105654, 0.99452192,
    -0.99573416, -0.96182561, -0.89516330, -0.79801720,
    -0.67369562, -0.52643216, -0.36124167, -0.18374951,
    -0.20791176, -0.40673670, -0.58778530, -0.74314487
};

static float *tns_tmp2_map[4] = {
    tns_tmp2_map_0_3,
    tns_tmp2_map_0_4,
    tns_tmp2_map_1_3,
    tns_tmp2_map_1_4
};

// Huffman tables
static const unsigned int aac_scalefactor_huffman_table[][2] = {
    /* codeword, code length */
    { 0x3FFE8, 18 },
    { 0x3FFE6, 18 },
    { 0x3FFE7, 18 },
    { 0x3FFE5, 18 },
    { 0x7FFF5, 19 },
    { 0x7FFF1, 19 },
    { 0x7FFED, 19 },
    { 0x7FFF6, 19 },
    { 0x7FFEE, 19 },
    { 0x7FFEF, 19 },
    { 0x7FFF0, 19 },
    { 0x7FFFC, 19 },
    { 0x7FFFD, 19 },
    { 0x7FFFF, 19 },
    { 0x7FFFE, 19 },
    { 0x7FFF7, 19 },
    { 0x7FFF8, 19 },
    { 0x7FFFB, 19 },
    { 0x7FFF9, 19 },
    { 0x3FFE4, 18 },
    { 0x7FFFA, 19 },
    { 0x3FFE3, 18 },
    { 0x1FFEF, 17 },
    { 0x1FFF0, 17 },
    {  0xFFF5, 16 },
    { 0x1FFEE, 17 },
    {  0xFFF2, 16 },
    {  0xFFF3, 16 },
    {  0xFFF4, 16 },
    {  0xFFF1, 16 },
    {  0x7FF6, 15 },
    {  0x7FF7, 15 },
    {  0x3FF9, 14 },
    {  0x3FF5, 14 },
    {  0x3FF7, 14 },
    {  0x3FF3, 14 },
    {  0x3FF6, 14 },
    {  0x3FF2, 14 },
    {  0x1FF7, 13 },
    {  0x1FF5, 13 },
    {   0xFF9, 12 },
    {   0xFF7, 12 },
    {   0xFF6, 12 },
    {   0x7F9, 11 },
    {   0xFF4, 12 },
    {   0x7F8, 11 },
    {   0x3F9, 10 },
    {   0x3F7, 10 },
    {   0x3F5, 10 },
    {   0x1F8,  9 },
    {   0x1F7,  9 },
    {    0xFA,  8 },
    {    0xF8,  8 },
    {    0xF6,  8 },
    {    0x79,  7 },
    {    0x3A,  6 },
    {    0x38,  6 },
    {    0x1A,  5 },
    {     0xB,  4 },
    {     0x4,  3 },
    {     0x0,  1 },
    {     0xA,  4 },
    {     0xC,  4 },
    {    0x1B,  5 },
    {    0x39,  6 },
    {    0x3B,  6 },
    {    0x78,  7 },
    {    0x7A,  7 },
    {    0xF7,  8 },
    {    0xF9,  8 },
    {   0x1F6,  9 },
    {   0x1F9,  9 },
    {   0x3F4, 10 },
    {   0x3F6, 10 },
    {   0x3F8, 10 },
    {   0x7F5, 11 },
    {   0x7F4, 11 },
    {   0x7F6, 11 },
    {   0x7F7, 11 },
    {   0xFF5, 12 },
    {   0xFF8, 12 },
    {  0x1FF4, 13 },
    {  0x1FF6, 13 },
    {  0x1FF8, 13 },
    {  0x3FF8, 14 },
    {  0x3FF4, 14 },
    {  0xFFF0, 16 },
    {  0x7FF4, 15 },
    {  0xFFF6, 16 },
    {  0x7FF5, 15 },
    { 0x3FFE2, 18 },
    { 0x7FFD9, 19 },
    { 0x7FFDA, 19 },
    { 0x7FFDB, 19 },
    { 0x7FFDC, 19 },
    { 0x7FFDD, 19 },
    { 0x7FFDE, 19 },
    { 0x7FFD8, 19 },
    { 0x7FFD2, 19 },
    { 0x7FFD3, 19 },
    { 0x7FFD4, 19 },
    { 0x7FFD5, 19 },
    { 0x7FFD6, 19 },
    { 0x7FFF2, 19 },
    { 0x7FFDF, 19 },
    { 0x7FFE7, 19 },
    { 0x7FFE8, 19 },
    { 0x7FFE9, 19 },
    { 0x7FFEA, 19 },
    { 0x7FFEB, 19 },
    { 0x7FFE6, 19 },
    { 0x7FFE0, 19 },
    { 0x7FFE1, 19 },
    { 0x7FFE2, 19 },
    { 0x7FFE3, 19 },
    { 0x7FFE4, 19 },
    { 0x7FFE5, 19 },
    { 0x7FFD7, 19 },
    { 0x7FFEC, 19 },
    { 0x7FFF4, 19 },
    { 0x7FFF3, 19 },
};

static const uint16_t codebook1[][2] = {
    {   0x7F8, 11 },
    {   0x1F1,  9 }, {   0x7FD, 11 }, {   0x3F5, 10 }, {    0x68,  7 },
    {   0x3F0, 10 }, {   0x7F7, 11 }, {   0x1EC,  9 }, {   0x7F5, 11 },
    {   0x3F1, 10 }, {    0x72,  7 }, {   0x3F4, 10 }, {    0x74,  7 },
    {    0x11,  5 }, {    0x76,  7 }, {   0x1EB,  9 }, {    0x6C,  7 },
    {   0x3F6, 10 }, {   0x7FC, 11 }, {   0x1E1,  9 }, {   0x7F1, 11 },
    {   0x1F0,  9 }, {    0x61,  7 }, {   0x1F6,  9 }, {   0x7F2, 11 },
    {   0x1EA,  9 }, {   0x7FB, 11 }, {   0x1F2,  9 }, {    0x69,  7 },
    {   0x1ED,  9 }, {    0x77,  7 }, {    0x17,  5 }, {    0x6F,  7 },
    {   0x1E6,  9 }, {    0x64,  7 }, {   0x1E5,  9 }, {    0x67,  7 },
    {    0x15,  5 }, {    0x62,  7 }, {    0x12,  5 }, {     0x0,  1 },
    {    0x14,  5 }, {    0x65,  7 }, {    0x16,  5 }, {    0x6D,  7 },
    {   0x1E9,  9 }, {    0x63,  7 }, {   0x1E4,  9 }, {    0x6B,  7 },
    {    0x13,  5 }, {    0x71,  7 }, {   0x1E3,  9 }, {    0x70,  7 },
    {   0x1F3,  9 }, {   0x7FE, 11 }, {   0x1E7,  9 }, {   0x7F3, 11 },
    {   0x1EF,  9 }, {    0x60,  7 }, {   0x1EE,  9 }, {   0x7F0, 11 },
    {   0x1E2,  9 }, {   0x7FA, 11 }, {   0x3F3, 10 }, {    0x6A,  7 },
    {   0x1E8,  9 }, {    0x75,  7 }, {    0x10,  5 }, {    0x73,  7 },
    {   0x1F4,  9 }, {    0x6E,  7 }, {   0x3F7, 10 }, {   0x7F6, 11 },
    {   0x1E0,  9 }, {   0x7F9, 11 }, {   0x3F2, 10 }, {    0x66,  7 },
    {   0x1F5,  9 }, {   0x7FF, 11 }, {   0x1F7,  9 }, {   0x7F4, 11 },
};

static const uint16_t codebook2[][2] = {
    {   0x1F3,  9 },
    {    0x6F,  7 }, {   0x1FD,  9 }, {    0xEB,  8 }, {    0x23,  6 },
    {    0xEA,  8 }, {   0x1F7,  9 }, {    0xE8,  8 }, {   0x1FA,  9 },
    {    0xF2,  8 }, {    0x2D,  6 }, {    0x70,  7 }, {    0x20,  6 },
    {     0x6,  5 }, {    0x2B,  6 }, {    0x6E,  7 }, {    0x28,  6 },
    {    0xE9,  8 }, {   0x1F9,  9 }, {    0x66,  7 }, {    0xF8,  8 },
    {    0xE7,  8 }, {    0x1B,  6 }, {    0xF1,  8 }, {   0x1F4,  9 },
    {    0x6B,  7 }, {   0x1F5,  9 }, {    0xEC,  8 }, {    0x2A,  6 },
    {    0x6C,  7 }, {    0x2C,  6 }, {     0xA,  5 }, {    0x27,  6 },
    {    0x67,  7 }, {    0x1A,  6 }, {    0xF5,  8 }, {    0x24,  6 },
    {     0x8,  5 }, {    0x1F,  6 }, {     0x9,  5 }, {     0x0,  3 },
    {     0x7,  5 }, {    0x1D,  6 }, {     0xB,  5 }, {    0x30,  6 },
    {    0xEF,  8 }, {    0x1C,  6 }, {    0x64,  7 }, {    0x1E,  6 },
    {     0xC,  5 }, {    0x29,  6 }, {    0xF3,  8 }, {    0x2F,  6 },
    {    0xF0,  8 }, {   0x1FC,  9 }, {    0x71,  7 }, {   0x1F2,  9 },
    {    0xF4,  8 }, {    0x21,  6 }, {    0xE6,  8 }, {    0xF7,  8 },
    {    0x68,  7 }, {   0x1F8,  9 }, {    0xEE,  8 }, {    0x22,  6 },
    {    0x65,  7 }, {    0x31,  6 }, {     0x2,  4 }, {    0x26,  6 },
    {    0xED,  8 }, {    0x25,  6 }, {    0x6A,  7 }, {   0x1FB,  9 },
    {    0x72,  7 }, {   0x1FE,  9 }, {    0x69,  7 }, {    0x2E,  6 },
    {    0xF6,  8 }, {   0x1FF,  9 }, {    0x6D,  7 }, {   0x1F6,  9 },
};

static const uint16_t codebook3[][2] = {
    {     0x0,  1 },
    {     0x9,  4 }, {    0xEF,  8 }, {     0xB,  4 }, {    0x19,  5 },
    {    0xF0,  8 }, {   0x1EB,  9 }, {   0x1E6,  9 }, {   0x3F2, 10 },
    {     0xA,  4 }, {    0x35,  6 }, {   0x1EF,  9 }, {    0x34,  6 },
    {    0x37,  6 }, {   0x1E9,  9 }, {   0x1ED,  9 }, {   0x1E7,  9 },
    {   0x3F3, 10 }, {   0x1EE,  9 }, {   0x3ED, 10 }, {  0x1FFA, 13 },
    {   0x1EC,  9 }, {   0x1F2,  9 }, {   0x7F9, 11 }, {   0x7F8, 11 },
    {   0x3F8, 10 }, {   0xFF8, 12 }, {     0x8,  4 }, {    0x38,  6 },
    {   0x3F6, 10 }, {    0x36,  6 }, {    0x75,  7 }, {   0x3F1, 10 },
    {   0x3EB, 10 }, {   0x3EC, 10 }, {   0xFF4, 12 }, {    0x18,  5 },
    {    0x76,  7 }, {   0x7F4, 11 }, {    0x39,  6 }, {    0x74,  7 },
    {   0x3EF, 10 }, {   0x1F3,  9 }, {   0x1F4,  9 }, {   0x7F6, 11 },
    {   0x1E8,  9 }, {   0x3EA, 10 }, {  0x1FFC, 13 }, {    0xF2,  8 },
    {   0x1F1,  9 }, {   0xFFB, 12 }, {   0x3F5, 10 }, {   0x7F3, 11 },
    {   0xFFC, 12 }, {    0xEE,  8 }, {   0x3F7, 10 }, {  0x7FFE, 15 },
    {   0x1F0,  9 }, {   0x7F5, 11 }, {  0x7FFD, 15 }, {  0x1FFB, 13 },
    {  0x3FFA, 14 }, {  0xFFFF, 16 }, {    0xF1,  8 }, {   0x3F0, 10 },
    {  0x3FFC, 14 }, {   0x1EA,  9 }, {   0x3EE, 10 }, {  0x3FFB, 14 },
    {   0xFF6, 12 }, {   0xFFA, 12 }, {  0x7FFC, 15 }, {   0x7F2, 11 },
    {   0xFF5, 12 }, {  0xFFFE, 16 }, {   0x3F4, 10 }, {   0x7F7, 11 },
    {  0x7FFB, 15 }, {   0xFF7, 12 }, {   0xFF9, 12 }, {  0x7FFA, 15 },
};

static const uint16_t codebook4[][2] = {
    {     0x7,  4 },
    {    0x16,  5 }, {    0xF6,  8 }, {    0x18,  5 }, {     0x8,  4 },
    {    0xEF,  8 }, {   0x1EF,  9 }, {    0xF3,  8 }, {   0x7F8, 11 },
    {    0x19,  5 }, {    0x17,  5 }, {    0xED,  8 }, {    0x15,  5 },
    {     0x1,  4 }, {    0xE2,  8 }, {    0xF0,  8 }, {    0x70,  7 },
    {   0x3F0, 10 }, {   0x1EE,  9 }, {    0xF1,  8 }, {   0x7FA, 11 },
    {    0xEE,  8 }, {    0xE4,  8 }, {   0x3F2, 10 }, {   0x7F6, 11 },
    {   0x3EF, 10 }, {   0x7FD, 11 }, {     0x5,  4 }, {    0x14,  5 },
    {    0xF2,  8 }, {     0x9,  4 }, {     0x4,  4 }, {    0xE5,  8 },
    {    0xF4,  8 }, {    0xE8,  8 }, {   0x3F4, 10 }, {     0x6,  4 },
    {     0x2,  4 }, {    0xE7,  8 }, {     0x3,  4 }, {     0x0,  4 },
    {    0x6B,  7 }, {    0xE3,  8 }, {    0x69,  7 }, {   0x1F3,  9 },
    {    0xEB,  8 }, {    0xE6,  8 }, {   0x3F6, 10 }, {    0x6E,  7 },
    {    0x6A,  7 }, {   0x1F4,  9 }, {   0x3EC, 10 }, {   0x1F0,  9 },
    {   0x3F9, 10 }, {    0xF5,  8 }, {    0xEC,  8 }, {   0x7FB, 11 },
    {    0xEA,  8 }, {    0x6F,  7 }, {   0x3F7, 10 }, {   0x7F9, 11 },
    {   0x3F3, 10 }, {   0xFFF, 12 }, {    0xE9,  8 }, {    0x6D,  7 },
    {   0x3F8, 10 }, {    0x6C,  7 }, {    0x68,  7 }, {   0x1F5,  9 },
    {   0x3EE, 10 }, {   0x1F2,  9 }, {   0x7F4, 11 }, {   0x7F7, 11 },
    {   0x3F1, 10 }, {   0xFFE, 12 }, {   0x3ED, 10 }, {   0x1F1,  9 },
    {   0x7F5, 11 }, {   0x7FE, 11 }, {   0x3F5, 10 }, {   0x7FC, 11 },
};

static const uint16_t codebook5[][2] = {
    {  0x1FFF, 13 },
    {   0xFF7, 12 }, {   0x7F4, 11 }, {   0x7E8, 11 }, {   0x3F1, 10 },
    {   0x7EE, 11 }, {   0x7F9, 11 }, {   0xFF8, 12 }, {  0x1FFD, 13 },
    {   0xFFD, 12 }, {   0x7F1, 11 }, {   0x3E8, 10 }, {   0x1E8,  9 },
    {    0xF0,  8 }, {   0x1EC,  9 }, {   0x3EE, 10 }, {   0x7F2, 11 },
    {   0xFFA, 12 }, {   0xFF4, 12 }, {   0x3EF, 10 }, {   0x1F2,  9 },
    {    0xE8,  8 }, {    0x70,  7 }, {    0xEC,  8 }, {   0x1F0,  9 },
    {   0x3EA, 10 }, {   0x7F3, 11 }, {   0x7EB, 11 }, {   0x1EB,  9 },
    {    0xEA,  8 }, {    0x1A,  5 }, {     0x8,  4 }, {    0x19,  5 },
    {    0xEE,  8 }, {   0x1EF,  9 }, {   0x7ED, 11 }, {   0x3F0, 10 },
    {    0xF2,  8 }, {    0x73,  7 }, {     0xB,  4 }, {     0x0,  1 },
    {     0xA,  4 }, {    0x71,  7 }, {    0xF3,  8 }, {   0x7E9, 11 },
    {   0x7EF, 11 }, {   0x1EE,  9 }, {    0xEF,  8 }, {    0x18,  5 },
    {     0x9,  4 }, {    0x1B,  5 }, {    0xEB,  8 }, {   0x1E9,  9 },
    {   0x7EC, 11 }, {   0x7F6, 11 }, {   0x3EB, 10 }, {   0x1F3,  9 },
    {    0xED,  8 }, {    0x72,  7 }, {    0xE9,  8 }, {   0x1F1,  9 },
    {   0x3ED, 10 }, {   0x7F7, 11 }, {   0xFF6, 12 }, {   0x7F0, 11 },
    {   0x3E9, 10 }, {   0x1ED,  9 }, {    0xF1,  8 }, {   0x1EA,  9 },
    {   0x3EC, 10 }, {   0x7F8, 11 }, {   0xFF9, 12 }, {  0x1FFC, 13 },
    {   0xFFC, 12 }, {   0xFF5, 12 }, {   0x7EA, 11 }, {   0x3F3, 10 },
    {   0x3F2, 10 }, {   0x7F5, 11 }, {   0xFFB, 12 }, {  0x1FFE, 13 },
};

static const uint16_t codebook6[][2] = {
    {   0x7FE, 11 },
    {   0x3FD, 10 }, {   0x1F1,  9 }, {   0x1EB,  9 }, {   0x1F4,  9 },
    {   0x1EA,  9 }, {   0x1F0,  9 }, {   0x3FC, 10 }, {   0x7FD, 11 },
    {   0x3F6, 10 }, {   0x1E5,  9 }, {    0xEA,  8 }, {    0x6C,  7 },
    {    0x71,  7 }, {    0x68,  7 }, {    0xF0,  8 }, {   0x1E6,  9 },
    {   0x3F7, 10 }, {   0x1F3,  9 }, {    0xEF,  8 }, {    0x32,  6 },
    {    0x27,  6 }, {    0x28,  6 }, {    0x26,  6 }, {    0x31,  6 },
    {    0xEB,  8 }, {   0x1F7,  9 }, {   0x1E8,  9 }, {    0x6F,  7 },
    {    0x2E,  6 }, {     0x8,  4 }, {     0x4,  4 }, {     0x6,  4 },
    {    0x29,  6 }, {    0x6B,  7 }, {   0x1EE,  9 }, {   0x1EF,  9 },
    {    0x72,  7 }, {    0x2D,  6 }, {     0x2,  4 }, {     0x0,  4 },
    {     0x3,  4 }, {    0x2F,  6 }, {    0x73,  7 }, {   0x1FA,  9 },
    {   0x1E7,  9 }, {    0x6E,  7 }, {    0x2B,  6 }, {     0x7,  4 },
    {     0x1,  4 }, {     0x5,  4 }, {    0x2C,  6 }, {    0x6D,  7 },
    {   0x1EC,  9 }, {   0x1F9,  9 }, {    0xEE,  8 }, {    0x30,  6 },
    {    0x24,  6 }, {    0x2A,  6 }, {    0x25,  6 }, {    0x33,  6 },
    {    0xEC,  8 }, {   0x1F2,  9 }, {   0x3F8, 10 }, {   0x1E4,  9 },
    {    0xED,  8 }, {    0x6A,  7 }, {    0x70,  7 }, {    0x69,  7 },
    {    0x74,  7 }, {    0xF1,  8 }, {   0x3FA, 10 }, {   0x7FF, 11 },
    {   0x3F9, 10 }, {   0x1F6,  9 }, {   0x1ED,  9 }, {   0x1F8,  9 },
    {   0x1E9,  9 }, {   0x1F5,  9 }, {   0x3FB, 10 }, {   0x7FC, 11 },
};

static const uint16_t codebook7[][2] = {
    {     0x0,  1 },
    {     0x5,  3 }, {    0x37,  6 }, {    0x74,  7 }, {    0xF2,  8 },
    {   0x1EB,  9 }, {   0x3ED, 10 }, {   0x7F7, 11 }, {     0x4,  3 },
    {     0xC,  4 }, {    0x35,  6 }, {    0x71,  7 }, {    0xEC,  8 },
    {    0xEE,  8 }, {   0x1EE,  9 }, {   0x1F5,  9 }, {    0x36,  6 },
    {    0x34,  6 }, {    0x72,  7 }, {    0xEA,  8 }, {    0xF1,  8 },
    {   0x1E9,  9 }, {   0x1F3,  9 }, {   0x3F5, 10 }, {    0x73,  7 },
    {    0x70,  7 }, {    0xEB,  8 }, {    0xF0,  8 }, {   0x1F1,  9 },
    {   0x1F0,  9 }, {   0x3EC, 10 }, {   0x3FA, 10 }, {    0xF3,  8 },
    {    0xED,  8 }, {   0x1E8,  9 }, {   0x1EF,  9 }, {   0x3EF, 10 },
    {   0x3F1, 10 }, {   0x3F9, 10 }, {   0x7FB, 11 }, {   0x1ED,  9 },
    {    0xEF,  8 }, {   0x1EA,  9 }, {   0x1F2,  9 }, {   0x3F3, 10 },
    {   0x3F8, 10 }, {   0x7F9, 11 }, {   0x7FC, 11 }, {   0x3EE, 10 },
    {   0x1EC,  9 }, {   0x1F4,  9 }, {   0x3F4, 10 }, {   0x3F7, 10 },
    {   0x7F8, 11 }, {   0xFFD, 12 }, {   0xFFE, 12 }, {   0x7F6, 11 },
    {   0x3F0, 10 }, {   0x3F2, 10 }, {   0x3F6, 10 }, {   0x7FA, 11 },
    {   0x7FD, 11 }, {   0xFFC, 12 }, {   0xFFF, 12 },
};

static const uint16_t codebook8[][2] = {
    {     0xE,  5 },
    {     0x5,  4 }, {    0x10,  5 }, {    0x30,  6 }, {    0x6F,  7 },
    {    0xF1,  8 }, {   0x1FA,  9 }, {   0x3FE, 10 }, {     0x3,  4 },
    {     0x0,  3 }, {     0x4,  4 }, {    0x12,  5 }, {    0x2C,  6 },
    {    0x6A,  7 }, {    0x75,  7 }, {    0xF8,  8 }, {     0xF,  5 },
    {     0x2,  4 }, {     0x6,  4 }, {    0x14,  5 }, {    0x2E,  6 },
    {    0x69,  7 }, {    0x72,  7 }, {    0xF5,  8 }, {    0x2F,  6 },
    {    0x11,  5 }, {    0x13,  5 }, {    0x2A,  6 }, {    0x32,  6 },
    {    0x6C,  7 }, {    0xEC,  8 }, {    0xFA,  8 }, {    0x71,  7 },
    {    0x2B,  6 }, {    0x2D,  6 }, {    0x31,  6 }, {    0x6D,  7 },
    {    0x70,  7 }, {    0xF2,  8 }, {   0x1F9,  9 }, {    0xEF,  8 },
    {    0x68,  7 }, {    0x33,  6 }, {    0x6B,  7 }, {    0x6E,  7 },
    {    0xEE,  8 }, {    0xF9,  8 }, {   0x3FC, 10 }, {   0x1F8,  9 },
    {    0x74,  7 }, {    0x73,  7 }, {    0xED,  8 }, {    0xF0,  8 },
    {    0xF6,  8 }, {   0x1F6,  9 }, {   0x1FD,  9 }, {   0x3FD, 10 },
    {    0xF3,  8 }, {    0xF4,  8 }, {    0xF7,  8 }, {   0x1F7,  9 },
    {   0x1FB,  9 }, {   0x1FC,  9 }, {   0x3FF, 10 },
};

static const uint16_t codebook9[][2] = {
    {     0x0,  1 },
    {     0x5,  3 }, {    0x37,  6 }, {    0xE7,  8 }, {   0x1DE,  9 },
    {   0x3CE, 10 }, {   0x3D9, 10 }, {   0x7C8, 11 }, {   0x7CD, 11 },
    {   0xFC8, 12 }, {   0xFDD, 12 }, {  0x1FE4, 13 }, {  0x1FEC, 13 },
    {     0x4,  3 }, {     0xC,  4 }, {    0x35,  6 }, {    0x72,  7 },
    {    0xEA,  8 }, {    0xED,  8 }, {   0x1E2,  9 }, {   0x3D1, 10 },
    {   0x3D3, 10 }, {   0x3E0, 10 }, {   0x7D8, 11 }, {   0xFCF, 12 },
    {   0xFD5, 12 }, {    0x36,  6 }, {    0x34,  6 }, {    0x71,  7 },
    {    0xE8,  8 }, {    0xEC,  8 }, {   0x1E1,  9 }, {   0x3CF, 10 },
    {   0x3DD, 10 }, {   0x3DB, 10 }, {   0x7D0, 11 }, {   0xFC7, 12 },
    {   0xFD4, 12 }, {   0xFE4, 12 }, {    0xE6,  8 }, {    0x70,  7 },
    {    0xE9,  8 }, {   0x1DD,  9 }, {   0x1E3,  9 }, {   0x3D2, 10 },
    {   0x3DC, 10 }, {   0x7CC, 11 }, {   0x7CA, 11 }, {   0x7DE, 11 },
    {   0xFD8, 12 }, {   0xFEA, 12 }, {  0x1FDB, 13 }, {   0x1DF,  9 },
    {    0xEB,  8 }, {   0x1DC,  9 }, {   0x1E6,  9 }, {   0x3D5, 10 },
    {   0x3DE, 10 }, {   0x7CB, 11 }, {   0x7DD, 11 }, {   0x7DC, 11 },
    {   0xFCD, 12 }, {   0xFE2, 12 }, {   0xFE7, 12 }, {  0x1FE1, 13 },
    {   0x3D0, 10 }, {   0x1E0,  9 }, {   0x1E4,  9 }, {   0x3D6, 10 },
    {   0x7C5, 11 }, {   0x7D1, 11 }, {   0x7DB, 11 }, {   0xFD2, 12 },
    {   0x7E0, 11 }, {   0xFD9, 12 }, {   0xFEB, 12 }, {  0x1FE3, 13 },
    {  0x1FE9, 13 }, {   0x7C4, 11 }, {   0x1E5,  9 }, {   0x3D7, 10 },
    {   0x7C6, 11 }, {   0x7CF, 11 }, {   0x7DA, 11 }, {   0xFCB, 12 },
    {   0xFDA, 12 }, {   0xFE3, 12 }, {   0xFE9, 12 }, {  0x1FE6, 13 },
    {  0x1FF3, 13 }, {  0x1FF7, 13 }, {   0x7D3, 11 }, {   0x3D8, 10 },
    {   0x3E1, 10 }, {   0x7D4, 11 }, {   0x7D9, 11 }, {   0xFD3, 12 },
    {   0xFDE, 12 }, {  0x1FDD, 13 }, {  0x1FD9, 13 }, {  0x1FE2, 13 },
    {  0x1FEA, 13 }, {  0x1FF1, 13 }, {  0x1FF6, 13 }, {   0x7D2, 11 },
    {   0x3D4, 10 }, {   0x3DA, 10 }, {   0x7C7, 11 }, {   0x7D7, 11 },
    {   0x7E2, 11 }, {   0xFCE, 12 }, {   0xFDB, 12 }, {  0x1FD8, 13 },
    {  0x1FEE, 13 }, {  0x3FF0, 14 }, {  0x1FF4, 13 }, {  0x3FF2, 14 },
    {   0x7E1, 11 }, {   0x3DF, 10 }, {   0x7C9, 11 }, {   0x7D6, 11 },
    {   0xFCA, 12 }, {   0xFD0, 12 }, {   0xFE5, 12 }, {   0xFE6, 12 },
    {  0x1FEB, 13 }, {  0x1FEF, 13 }, {  0x3FF3, 14 }, {  0x3FF4, 14 },
    {  0x3FF5, 14 }, {   0xFE0, 12 }, {   0x7CE, 11 }, {   0x7D5, 11 },
    {   0xFC6, 12 }, {   0xFD1, 12 }, {   0xFE1, 12 }, {  0x1FE0, 13 },
    {  0x1FE8, 13 }, {  0x1FF0, 13 }, {  0x3FF1, 14 }, {  0x3FF8, 14 },
    {  0x3FF6, 14 }, {  0x7FFC, 15 }, {   0xFE8, 12 }, {   0x7DF, 11 },
    {   0xFC9, 12 }, {   0xFD7, 12 }, {   0xFDC, 12 }, {  0x1FDC, 13 },
    {  0x1FDF, 13 }, {  0x1FED, 13 }, {  0x1FF5, 13 }, {  0x3FF9, 14 },
    {  0x3FFB, 14 }, {  0x7FFD, 15 }, {  0x7FFE, 15 }, {  0x1FE7, 13 },
    {   0xFCC, 12 }, {   0xFD6, 12 }, {   0xFDF, 12 }, {  0x1FDE, 13 },
    {  0x1FDA, 13 }, {  0x1FE5, 13 }, {  0x1FF2, 13 }, {  0x3FFA, 14 },
    {  0x3FF7, 14 }, {  0x3FFC, 14 }, {  0x3FFD, 14 }, {  0x7FFF, 15 },
};

static const uint16_t codebook10[][2] = {
    {    0x22,  6 },
    {     0x8,  5 }, {    0x1D,  6 }, {    0x26,  6 }, {    0x5F,  7 },
    {    0xD3,  8 }, {   0x1CF,  9 }, {   0x3D0, 10 }, {   0x3D7, 10 },
    {   0x3ED, 10 }, {   0x7F0, 11 }, {   0x7F6, 11 }, {   0xFFD, 12 },
    {     0x7,  5 }, {     0x0,  4 }, {     0x1,  4 }, {     0x9,  5 },
    {    0x20,  6 }, {    0x54,  7 }, {    0x60,  7 }, {    0xD5,  8 },
    {    0xDC,  8 }, {   0x1D4,  9 }, {   0x3CD, 10 }, {   0x3DE, 10 },
    {   0x7E7, 11 }, {    0x1C,  6 }, {     0x2,  4 }, {     0x6,  5 },
    {     0xC,  5 }, {    0x1E,  6 }, {    0x28,  6 }, {    0x5B,  7 },
    {    0xCD,  8 }, {    0xD9,  8 }, {   0x1CE,  9 }, {   0x1DC,  9 },
    {   0x3D9, 10 }, {   0x3F1, 10 }, {    0x25,  6 }, {     0xB,  5 },
    {     0xA,  5 }, {     0xD,  5 }, {    0x24,  6 }, {    0x57,  7 },
    {    0x61,  7 }, {    0xCC,  8 }, {    0xDD,  8 }, {   0x1CC,  9 },
    {   0x1DE,  9 }, {   0x3D3, 10 }, {   0x3E7, 10 }, {    0x5D,  7 },
    {    0x21,  6 }, {    0x1F,  6 }, {    0x23,  6 }, {    0x27,  6 },
    {    0x59,  7 }, {    0x64,  7 }, {    0xD8,  8 }, {    0xDF,  8 },
    {   0x1D2,  9 }, {   0x1E2,  9 }, {   0x3DD, 10 }, {   0x3EE, 10 },
    {    0xD1,  8 }, {    0x55,  7 }, {    0x29,  6 }, {    0x56,  7 },
    {    0x58,  7 }, {    0x62,  7 }, {    0xCE,  8 }, {    0xE0,  8 },
    {    0xE2,  8 }, {   0x1DA,  9 }, {   0x3D4, 10 }, {   0x3E3, 10 },
    {   0x7EB, 11 }, {   0x1C9,  9 }, {    0x5E,  7 }, {    0x5A,  7 },
    {    0x5C,  7 }, {    0x63,  7 }, {    0xCA,  8 }, {    0xDA,  8 },
    {   0x1C7,  9 }, {   0x1CA,  9 }, {   0x1E0,  9 }, {   0x3DB, 10 },
    {   0x3E8, 10 }, {   0x7EC, 11 }, {   0x1E3,  9 }, {    0xD2,  8 },
    {    0xCB,  8 }, {    0xD0,  8 }, {    0xD7,  8 }, {    0xDB,  8 },
    {   0x1C6,  9 }, {   0x1D5,  9 }, {   0x1D8,  9 }, {   0x3CA, 10 },
    {   0x3DA, 10 }, {   0x7EA, 11 }, {   0x7F1, 11 }, {   0x1E1,  9 },
    {    0xD4,  8 }, {    0xCF,  8 }, {    0xD6,  8 }, {    0xDE,  8 },
    {    0xE1,  8 }, {   0x1D0,  9 }, {   0x1D6,  9 }, {   0x3D1, 10 },
    {   0x3D5, 10 }, {   0x3F2, 10 }, {   0x7EE, 11 }, {   0x7FB, 11 },
    {   0x3E9, 10 }, {   0x1CD,  9 }, {   0x1C8,  9 }, {   0x1CB,  9 },
    {   0x1D1,  9 }, {   0x1D7,  9 }, {   0x1DF,  9 }, {   0x3CF, 10 },
    {   0x3E0, 10 }, {   0x3EF, 10 }, {   0x7E6, 11 }, {   0x7F8, 11 },
    {   0xFFA, 12 }, {   0x3EB, 10 }, {   0x1DD,  9 }, {   0x1D3,  9 },
    {   0x1D9,  9 }, {   0x1DB,  9 }, {   0x3D2, 10 }, {   0x3CC, 10 },
    {   0x3DC, 10 }, {   0x3EA, 10 }, {   0x7ED, 11 }, {   0x7F3, 11 },
    {   0x7F9, 11 }, {   0xFF9, 12 }, {   0x7F2, 11 }, {   0x3CE, 10 },
    {   0x1E4,  9 }, {   0x3CB, 10 }, {   0x3D8, 10 }, {   0x3D6, 10 },
    {   0x3E2, 10 }, {   0x3E5, 10 }, {   0x7E8, 11 }, {   0x7F4, 11 },
    {   0x7F5, 11 }, {   0x7F7, 11 }, {   0xFFB, 12 }, {   0x7FA, 11 },
    {   0x3EC, 10 }, {   0x3DF, 10 }, {   0x3E1, 10 }, {   0x3E4, 10 },
    {   0x3E6, 10 }, {   0x3F0, 10 }, {   0x7E9, 11 }, {   0x7EF, 11 },
    {   0xFF8, 12 }, {   0xFFE, 12 }, {   0xFFC, 12 }, {   0xFFF, 12 },
};

static const uint16_t codebook11[][2] = {
    {     0x0,  4 },
    {     0x6,  5 }, {    0x19,  6 }, {    0x3D,  7 }, {    0x9C,  8 },
    {    0xC6,  8 }, {   0x1A7,  9 }, {   0x390, 10 }, {   0x3C2, 10 },
    {   0x3DF, 10 }, {   0x7E6, 11 }, {   0x7F3, 11 }, {   0xFFB, 12 },
    {   0x7EC, 11 }, {   0xFFA, 12 }, {   0xFFE, 12 }, {   0x38E, 10 },
    {     0x5,  5 }, {     0x1,  4 }, {     0x8,  5 }, {    0x14,  6 },
    {    0x37,  7 }, {    0x42,  7 }, {    0x92,  8 }, {    0xAF,  8 },
    {   0x191,  9 }, {   0x1A5,  9 }, {   0x1B5,  9 }, {   0x39E, 10 },
    {   0x3C0, 10 }, {   0x3A2, 10 }, {   0x3CD, 10 }, {   0x7D6, 11 },
    {    0xAE,  8 }, {    0x17,  6 }, {     0x7,  5 }, {     0x9,  5 },
    {    0x18,  6 }, {    0x39,  7 }, {    0x40,  7 }, {    0x8E,  8 },
    {    0xA3,  8 }, {    0xB8,  8 }, {   0x199,  9 }, {   0x1AC,  9 },
    {   0x1C1,  9 }, {   0x3B1, 10 }, {   0x396, 10 }, {   0x3BE, 10 },
    {   0x3CA, 10 }, {    0x9D,  8 }, {    0x3C,  7 }, {    0x15,  6 },
    {    0x16,  6 }, {    0x1A,  6 }, {    0x3B,  7 }, {    0x44,  7 },
    {    0x91,  8 }, {    0xA5,  8 }, {    0xBE,  8 }, {   0x196,  9 },
    {   0x1AE,  9 }, {   0x1B9,  9 }, {   0x3A1, 10 }, {   0x391, 10 },
    {   0x3A5, 10 }, {   0x3D5, 10 }, {    0x94,  8 }, {    0x9A,  8 },
    {    0x36,  7 }, {    0x38,  7 }, {    0x3A,  7 }, {    0x41,  7 },
    {    0x8C,  8 }, {    0x9B,  8 }, {    0xB0,  8 }, {    0xC3,  8 },
    {   0x19E,  9 }, {   0x1AB,  9 }, {   0x1BC,  9 }, {   0x39F, 10 },
    {   0x38F, 10 }, {   0x3A9, 10 }, {   0x3CF, 10 }, {    0x93,  8 },
    {    0xBF,  8 }, {    0x3E,  7 }, {    0x3F,  7 }, {    0x43,  7 },
    {    0x45,  7 }, {    0x9E,  8 }, {    0xA7,  8 }, {    0xB9,  8 },
    {   0x194,  9 }, {   0x1A2,  9 }, {   0x1BA,  9 }, {   0x1C3,  9 },
    {   0x3A6, 10 }, {   0x3A7, 10 }, {   0x3BB, 10 }, {   0x3D4, 10 },
    {    0x9F,  8 }, {   0x1A0,  9 }, {    0x8F,  8 }, {    0x8D,  8 },
    {    0x90,  8 }, {    0x98,  8 }, {    0xA6,  8 }, {    0xB6,  8 },
    {    0xC4,  8 }, {   0x19F,  9 }, {   0x1AF,  9 }, {   0x1BF,  9 },
    {   0x399, 10 }, {   0x3BF, 10 }, {   0x3B4, 10 }, {   0x3C9, 10 },
    {   0x3E7, 10 }, {    0xA8,  8 }, {   0x1B6,  9 }, {    0xAB,  8 },
    {    0xA4,  8 }, {    0xAA,  8 }, {    0xB2,  8 }, {    0xC2,  8 },
    {    0xC5,  8 }, {   0x198,  9 }, {   0x1A4,  9 }, {   0x1B8,  9 },
    {   0x38C, 10 }, {   0x3A4, 10 }, {   0x3C4, 10 }, {   0x3C6, 10 },
    {   0x3DD, 10 }, {   0x3E8, 10 }, {    0xAD,  8 }, {   0x3AF, 10 },
    {   0x192,  9 }, {    0xBD,  8 }, {    0xBC,  8 }, {   0x18E,  9 },
    {   0x197,  9 }, {   0x19A,  9 }, {   0x1A3,  9 }, {   0x1B1,  9 },
    {   0x38D, 10 }, {   0x398, 10 }, {   0x3B7, 10 }, {   0x3D3, 10 },
    {   0x3D1, 10 }, {   0x3DB, 10 }, {   0x7DD, 11 }, {    0xB4,  8 },
    {   0x3DE, 10 }, {   0x1A9,  9 }, {   0x19B,  9 }, {   0x19C,  9 },
    {   0x1A1,  9 }, {   0x1AA,  9 }, {   0x1AD,  9 }, {   0x1B3,  9 },
    {   0x38B, 10 }, {   0x3B2, 10 }, {   0x3B8, 10 }, {   0x3CE, 10 },
    {   0x3E1, 10 }, {   0x3E0, 10 }, {   0x7D2, 11 }, {   0x7E5, 11 },
    {    0xB7,  8 }, {   0x7E3, 11 }, {   0x1BB,  9 }, {   0x1A8,  9 },
    {   0x1A6,  9 }, {   0x1B0,  9 }, {   0x1B2,  9 }, {   0x1B7,  9 },
    {   0x39B, 10 }, {   0x39A, 10 }, {   0x3BA, 10 }, {   0x3B5, 10 },
    {   0x3D6, 10 }, {   0x7D7, 11 }, {   0x3E4, 10 }, {   0x7D8, 11 },
    {   0x7EA, 11 }, {    0xBA,  8 }, {   0x7E8, 11 }, {   0x3A0, 10 },
    {   0x1BD,  9 }, {   0x1B4,  9 }, {   0x38A, 10 }, {   0x1C4,  9 },
    {   0x392, 10 }, {   0x3AA, 10 }, {   0x3B0, 10 }, {   0x3BC, 10 },
    {   0x3D7, 10 }, {   0x7D4, 11 }, {   0x7DC, 11 }, {   0x7DB, 11 },
    {   0x7D5, 11 }, {   0x7F0, 11 }, {    0xC1,  8 }, {   0x7FB, 11 },
    {   0x3C8, 10 }, {   0x3A3, 10 }, {   0x395, 10 }, {   0x39D, 10 },
    {   0x3AC, 10 }, {   0x3AE, 10 }, {   0x3C5, 10 }, {   0x3D8, 10 },
    {   0x3E2, 10 }, {   0x3E6, 10 }, {   0x7E4, 11 }, {   0x7E7, 11 },
    {   0x7E0, 11 }, {   0x7E9, 11 }, {   0x7F7, 11 }, {   0x190,  9 },
    {   0x7F2, 11 }, {   0x393, 10 }, {   0x1BE,  9 }, {   0x1C0,  9 },
    {   0x394, 10 }, {   0x397, 10 }, {   0x3AD, 10 }, {   0x3C3, 10 },
    {   0x3C1, 10 }, {   0x3D2, 10 }, {   0x7DA, 11 }, {   0x7D9, 11 },
    {   0x7DF, 11 }, {   0x7EB, 11 }, {   0x7F4, 11 }, {   0x7FA, 11 },
    {   0x195,  9 }, {   0x7F8, 11 }, {   0x3BD, 10 }, {   0x39C, 10 },
    {   0x3AB, 10 }, {   0x3A8, 10 }, {   0x3B3, 10 }, {   0x3B9, 10 },
    {   0x3D0, 10 }, {   0x3E3, 10 }, {   0x3E5, 10 }, {   0x7E2, 11 },
    {   0x7DE, 11 }, {   0x7ED, 11 }, {   0x7F1, 11 }, {   0x7F9, 11 },
    {   0x7FC, 11 }, {   0x193,  9 }, {   0xFFD, 12 }, {   0x3DC, 10 },
    {   0x3B6, 10 }, {   0x3C7, 10 }, {   0x3CC, 10 }, {   0x3CB, 10 },
    {   0x3D9, 10 }, {   0x3DA, 10 }, {   0x7D3, 11 }, {   0x7E1, 11 },
    {   0x7EE, 11 }, {   0x7EF, 11 }, {   0x7F5, 11 }, {   0x7F6, 11 },
    {   0xFFC, 12 }, {   0xFFF, 12 }, {   0x19D,  9 }, {   0x1C2,  9 },
    {    0xB5,  8 }, {    0xA1,  8 }, {    0x96,  8 }, {    0x97,  8 },
    {    0x95,  8 }, {    0x99,  8 }, {    0xA0,  8 }, {    0xA2,  8 },
    {    0xAC,  8 }, {    0xA9,  8 }, {    0xB1,  8 }, {    0xB3,  8 },
    {    0xBB,  8 }, {    0xC0,  8 }, {   0x18F,  9 }, {     0x4,  5 },
};
