#ifndef KISS_IMDCT_H
#define KISS_IMDCT_H

#define kiss_fft_scalar float

typedef struct kiss_imdct_s kiss_imdct_t;

kiss_imdct_t * kiss_imdct_init(int n);
void kiss_imdct(kiss_imdct_t * s, kiss_fft_scalar * inbuffer, kiss_fft_scalar * outbuffer);
void kiss_imdct_end(kiss_imdct_t * s);

#endif
