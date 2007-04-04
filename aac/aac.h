
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
    AOT_HVXC,
    AOT_TTSI = 12,
    AOT_MAINSYNTH,
    AOT_WAVESYNTH,
    AOT_MIDI,
    AOT_SAFX,
    AOT_ER_AAC_LC,
    AOT_ER_AAC_LTP = 19,
    AOT_ER_AAC_SCALABLE,
    AOT_ER_TWINVQ,
    AOT_ER_BSAC,
    AOT_ER_AAC_LD,
    AOT_ER_CELP,
    AOT_ER_HVXC,
    AOT_ER_HILN,
    AOT_ER_PARAM,
    AOT_SSC
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

// IDs for extension_payload
enum {
    EXT_FILL = 0x0,
    EXT_FILL_DATA,
    EXT_DATA_ELEMENT,
    EXT_DYNAMIC_RANGE = 0xb,
    EXT_SBR_DATA = 0xd,
    EXT_SBR_DATA_CRC = 0xe
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

//ltp
#define MAX_LTP_LONG_SFB 40

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
    int generated;

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

enum {
    MIXMODE_DEFAULT = 0,
    MIXMODE_1TO1,
    MIXMODE_2TO1,
    MIXMODE_1TO2,
    MIXMODE_2TO2,
    MIXMODE_MATRIX1,
    MIXMODE_MATRIX2,
    MIXMODE_UNKNOWN
};

typedef struct {
    int mode;
    int c_tag;
    int lr_tag;
    int sur_tag;
    float sce_gain[MAX_TAGID];
    float cpe_gain[MAX_TAGID][2];
    float lfe_gain[MAX_TAGID];
} mix_config_struct;

typedef struct {
    int present;
    int lag;
    float coef;
    int used[MAX_LTP_LONG_SFB];
} ltp_struct;

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
    // ltp
    ltp_struct ltp;
    ltp_struct ltp2;
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
    int max_band;
    int adjust_num[4][8];
    int alev[4][8][8];
    int aloc[4][8][8];
    float buf[4][24];
} ssr_struct;

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
    float mixing_gain;
    ics_struct ics;
    tns_struct tns;
    int cb[8][64];   // codebooks
    float sf[8][64];
    DECLARE_ALIGNED_16(float, coeffs[1024]);
    DECLARE_ALIGNED_16(float, saved[1024]);
    DECLARE_ALIGNED_16(float, ret[1024]);
    int16_t *ltp_state;
    ssr_struct *ssr;
} sce_struct;

// channel element
typedef struct {
    int common_window;
    ms_struct ms;
    sce_struct ch[2];
} cpe_struct;

typedef struct {
    coupling_struct coup;
    sce_struct ch;
} cc_struct;

typedef struct {
    float q[4][4];
    float t0[4][12];
    float t1[4][12];
} ssr_context;

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
    mix_config_struct mix;
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
    MDCTContext *mdct_ltp;
    DSPContext dsp;
    int * vq[11];
    ssr_context * ssrctx;

    // statistics
    int num_frame;
} aac_context_t;

