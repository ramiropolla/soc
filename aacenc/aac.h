/* this file is a mere collection of things that were borrowed from
 * GSoC AAC decoder and should be put into common aac.h
 * while the merge is not done, declarations should reside here
 */
#ifndef ERSATZ_AAC_H
#define ERSATZ_AAC_H

#define MAX_SWB_SIZE  51

DECLARE_ALIGNED_16(static float, kbd_long_1024[1024]);
DECLARE_ALIGNED_16(static float, kbd_short_128[128]);
DECLARE_ALIGNED_16(static float, sine_long_1024[1024]);
DECLARE_ALIGNED_16(static float, sine_short_128[128]);

/**
 * window sequences
 */
enum WindowSequence {
    ONLY_LONG_SEQUENCE,
    LONG_START_SEQUENCE,
    EIGHT_SHORT_SEQUENCE,
    LONG_STOP_SEQUENCE,
};

/**
 * IDs for raw_data_block
 */
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

/**
 * special codebooks
 */
enum Codebook {
    ZERO_HCB       = 0,
    FIRST_PAIR_HCB = 5,
    ESC_HCB        = 11,
    NOISE_HCB      = 13,
    INTENSITY_HCB2 = 14,
    INTENSITY_HCB  = 15,
    ESC_FLAG       = 16,
};

/**
 * pulse tool
 */
typedef struct {
    int present;
    int num_pulse;
    int start;
    int offset[4];
    int amp[4];
} Pulse;

#define MAX_TAGID 16

/**
 * Program configuration - describes how channels are arranged. Either read from
 * stream (ID_PCE) or created based on a default fixed channel arrangement.
 */
typedef struct {
    int che_type[4][MAX_TAGID]; ///< channel element type with the first index as the first 4 raw_data_block IDs
    int mono_mixdown;           ///< The SCE tag to use if user requests mono   output, -1 if not available.
    int stereo_mixdown;         ///< The CPE tag to use if user requests stereo output, -1 if not available.
    int matrix_mixdown;         ///< The CPE tag to use if user requests matrixed stereo output, -1 if not available.
    int mixdown_coeff_index;    ///< 0-3
    int pseudo_surround;        ///< Mix surround channels out of phase.
} ProgramConfig;


/**
 * Individual Channel Stream
 */
typedef struct {
    int intensity_present;
    uint8_t max_sfb;            ///< number of scalefactor bands per group
    enum WindowSequence window_sequence;
    enum WindowSequence window_sequence_prev;
    uint8_t use_kb_window[2];   ///< If set, use Kaiser-Bessel window, otherwise use a sinus window.
    int num_window_groups;
    uint8_t grouping;
    uint8_t group_len[8];
    const uint8_t *swb_sizes;
    int num_swb;
    int num_windows;
    int tns_max_bands;
} IndividualChannelStream;

#define TNS_MAX_ORDER 20
/**
 * Temporal Noise Shaping
 */
typedef struct {
    int present;
    int n_filt[8];
    int length[8][4];
    int direction[8][4];
    int order[8][4];
    int coef_res[8];
    int coef_compress[8][4];
    int coef_len[8][4];
    const float *tmp2_map[8][4];
    int coef[8][4][TNS_MAX_ORDER];
} TemporalNoiseShaping;

/**
 * M/S joint channel coding
 */
typedef struct {
    int present;
    uint8_t mask[8][64];
} MidSideStereo;

/**
 * Single Channel Element
 * Used for both SCE and LFE elements
 */
typedef struct {
    int gain;                                 /**< Channel gain (not used by AAC bitstream).
                                               *   Note that this is applied before joint stereo decoding.
                                               *   Thus, when used inside CPE elements, both channels must have equal gain.
                                               */
    IndividualChannelStream ics;
    TemporalNoiseShaping tns;
    Pulse pulse;
    int zeroes[8][64];
    int sf_idx[8][64];
    enum Codebook cb[8][64];                  ///< codebooks
    int cb_run_end[8][64];                    ///< codebook run end points
    float sf[8][64];                          ///< scalefactors
    DECLARE_ALIGNED_16(float, coeffs[1024]);  ///< coefficients for IMDCT
    DECLARE_ALIGNED_16(float, saved[1024]);   ///< overlap
    DECLARE_ALIGNED_16(float, ret[1024]);     ///< PCM output
    DECLARE_ALIGNED_16(int,   icoefs[1024]);  ///< integer coefficients for coding
} SingleChannelElement;

/**
 * channel element - generic struct for SCE/CPE/CCE/LFE
 */
typedef struct {
    // CPE specific
    int common_window;     ///< Set if channels share a common 'IndividualChannelStream' in bitstream.
    MidSideStereo ms;
    // shared
    SingleChannelElement ch[2];
    // CCE specific
//    ChannelCoupling coup;
} ChannelElement;

//my stuff

#define SCALE_ONE_POS   140    ///< scalefactor index that corresponds to scale=1.0
#define SCALE_MAX_POS   255    ///< scalefactor index maximum value
#define SCALE_MAX_DIFF   60    ///< maximum scalefactor difference allowed by standard
#define SCALE_DIFF_ZERO  60    ///< codebook index corresponding to zero scalefactor indices difference


#endif

