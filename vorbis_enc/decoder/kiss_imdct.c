/* Copyright (c) 2003-2004, Mark Borgerding
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * * Neither the author nor the names of any contributors may be used to
 *   endorse or promote products derived from this software without
 *   specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <inttypes.h>
#include <malloc.h>

#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif

#ifdef USE_FFTW
#include <fftw3.h>
#endif

#ifdef USE_SIMD
# include <xmmintrin.h>
# define kiss_fft_scalar __m128
#else
//# define KISS_FFT_MALLOC malloc
#endif
# define KISS_FFT_MALLOC(nbytes) memalign(16,nbytes)

#ifdef FIXED_POINT
#include <sys/types.h>
# if (FIXED_POINT == 32)
#  define kiss_fft_scalar int32_t
# else
#  define kiss_fft_scalar int16_t
# endif
#else
# ifndef kiss_fft_scalar
#   define kiss_fft_scalar float
# endif
#endif

typedef struct {
    kiss_fft_scalar r;
    kiss_fft_scalar i;
} kiss_fft_cpx;

#define MAXFACTORS 32
// e.g. an fft of length 128 has 4 factors as far as kissfft is concerned 4*4*4*2

typedef struct kiss_fft_state {
    int nfft;
    int inverse;
    int factors[2*MAXFACTORS];
    kiss_fft_cpx twiddles[1];
} * kiss_fft_cfg;

/* Explanation of macros dealing with complex math:

   C_MUL(m,a,b)         : m = a*b
   C_FIXDIV( c , div )  : if a fixed point impl., c /= div. noop otherwise
   C_SUB( res, a,b)     : res = a - b
   C_SUBFROM( res , a)  : res -= a
   C_ADDTO( res , a)    : res += a
*/

#ifdef FIXED_POINT
#if (FIXED_POINT==32)
# define FRACBITS 14
# define SAMPPROD int64_t
# define SAMP_MAX (1<<FRACBITS)
#else
# define FRACBITS 15
# define SAMPPROD int32_t
# define SAMP_MAX 32767
#endif

#define SAMP_MIN -SAMP_MAX

#if defined(CHECK_OVERFLOW)
#  define CHECK_OVERFLOW_OP(a,op,b)  \
	if ( (SAMPPROD)(a) op (SAMPPROD)(b) > SAMP_MAX || (SAMPPROD)(a) op (SAMPPROD)(b) < SAMP_MIN ) { \
		fprintf(stderr,"WARNING:overflow @ " __FILE__ "(%d): (%d " #op" %d) = %ld\n",__LINE__,(a),(b),(SAMPPROD)(a) op (SAMPPROD)(b) );  }
#endif


#   define smul(a,b) ( (SAMPPROD)(a)*(b) )
#   define sround( x )  (kiss_fft_scalar)( ( (x) + (1<<(FRACBITS-1)) ) >> FRACBITS )

#   define S_MUL(a,b) sround( smul(a,b) )

#   define C_MUL(m,a,b) \
      do{ (m).r = sround( smul((a).r,(b).r) - smul((a).i,(b).i) ); \
          (m).i = sround( smul((a).r,(b).i) + smul((a).i,(b).r) ); }while(0)

#   define DIVSCALAR(x,k) \
	(x) = sround( smul(  x, SAMP_MAX/k ) )

#   define C_FIXDIV(c,div) \
	do {    DIVSCALAR( (c).r , div);  \
		DIVSCALAR( (c).i  , div); }while (0)

#   define C_MULBYSCALAR( c, s ) \
    do{ (c).r =  sround( smul( (c).r , s ) ) ;\
        (c).i =  sround( smul( (c).i , s ) ) ; }while(0)

#else  /* not FIXED_POINT*/

#   define S_MUL(a,b) ( (a)*(b) )
#define C_MUL(m,a,b) \
    do{ (m).r = (a).r*(b).r - (a).i*(b).i;\
        (m).i = (a).r*(b).i + (a).i*(b).r; }while(0)
#   define C_FIXDIV(c,div) /* NOOP */
#   define C_MULBYSCALAR( c, s ) \
    do{ (c).r *= (s);\
        (c).i *= (s); }while(0)
#endif

#ifndef CHECK_OVERFLOW_OP
#  define CHECK_OVERFLOW_OP(a,op,b) /* noop */
#endif

#define  C_ADD( res, a,b)\
    do { \
	    CHECK_OVERFLOW_OP((a).r,+,(b).r)\
	    CHECK_OVERFLOW_OP((a).i,+,(b).i)\
	    (res).r=(a).r+(b).r;  (res).i=(a).i+(b).i; \
    }while(0)
#define  C_SUB( res, a,b)\
    do { \
	    CHECK_OVERFLOW_OP((a).r,-,(b).r)\
	    CHECK_OVERFLOW_OP((a).i,-,(b).i)\
	    (res).r=(a).r-(b).r;  (res).i=(a).i-(b).i; \
    }while(0)
#define C_ADDTO( res , a)\
    do { \
	    CHECK_OVERFLOW_OP((res).r,+,(a).r)\
	    CHECK_OVERFLOW_OP((res).i,+,(a).i)\
	    (res).r += (a).r;  (res).i += (a).i;\
    }while(0)

#define C_SUBFROM( res , a)\
    do {\
	    CHECK_OVERFLOW_OP((res).r,-,(a).r)\
	    CHECK_OVERFLOW_OP((res).i,-,(a).i)\
	    (res).r -= (a).r;  (res).i -= (a).i; \
    }while(0)


#ifdef FIXED_POINT
#  define KISS_FFT_COS(phase)  floor(.5+SAMP_MAX * cos (phase))
#  define KISS_FFT_SIN(phase)  floor(.5+SAMP_MAX * sin (phase))
#  define HALF_OF(x) ((x)>>1)
#elif defined(USE_SIMD)
#  define KISS_FFT_COS(phase) _mm_set1_ps( cos(phase) )
#  define KISS_FFT_SIN(phase) _mm_set1_ps( sin(phase) )
#  define HALF_OF(x) ((x)*_mm_set1_ps(.5))
#else
#  define KISS_FFT_COS(phase) (kiss_fft_scalar) cos(phase)
#  define KISS_FFT_SIN(phase) (kiss_fft_scalar) sin(phase)
#  define HALF_OF(x) ((x)*.5)
#endif

#define kf_cexp(x,phase) do{ (x)->r = KISS_FFT_COS(phase); (x)->i = KISS_FFT_SIN(phase); }while(0)

#if !defined(USE_FFTW) && !defined(USE_FFMPEG)

static kiss_fft_cpx *scratchbuf=NULL;
static size_t nscratchbuf=0;
static kiss_fft_cpx *tmpbuf=NULL;
static size_t ntmpbuf=0;

#define CHECKBUF(buf,nbuf,n) \
    do { \
        if ( nbuf < (size_t)(n) ) {\
            free(buf); \
            buf = (kiss_fft_cpx*)KISS_FFT_MALLOC(sizeof(kiss_fft_cpx)*(n)); \
            nbuf = (size_t)(n); \
        } \
   }while(0)


static void kf_bfly2(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx * Fout2;
    kiss_fft_cpx * tw1 = st->twiddles;
    kiss_fft_cpx t;
    Fout2 = Fout + m;
    do{
        //C_FIXDIV(*Fout,2); C_FIXDIV(*Fout2,2);

        C_MUL (t,  *Fout2 , *tw1);
        tw1 += fstride;
        C_SUB( *Fout2 ,  *Fout , t );
        C_ADDTO( *Fout ,  t );
        ++Fout2;
        ++Fout;
    }while (--m);
}

static void kf_bfly4(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        const size_t m
        )
{
    kiss_fft_cpx *tw1,*tw2,*tw3;
    kiss_fft_cpx scratch[6];
    size_t k=m;
    const size_t m2=2*m;
    const size_t m3=3*m;

    tw3 = tw2 = tw1 = st->twiddles;

    do {
        //C_FIXDIV(*Fout,2); C_FIXDIV(Fout[m],2); C_FIXDIV(Fout[m2],2); C_FIXDIV(Fout[m3],2);

        C_MUL(scratch[0],Fout[m] , *tw1 );
        C_MUL(scratch[1],Fout[m2] , *tw2 );
        C_MUL(scratch[2],Fout[m3] , *tw3 );

        C_SUB( scratch[5] , *Fout, scratch[1] );
        C_ADDTO(*Fout, scratch[1]);
        C_ADD( scratch[3] , scratch[0] , scratch[2] );
        C_SUB( scratch[4] , scratch[0] , scratch[2] );
        C_SUB( Fout[m2], *Fout, scratch[3] );
        tw1 += fstride;
        tw2 += fstride*2;
        tw3 += fstride*3;
        C_ADDTO( *Fout , scratch[3] );

        if(st->inverse) {
            Fout[m].r = scratch[5].r - scratch[4].i;
            Fout[m].i = scratch[5].i + scratch[4].r;
            Fout[m3].r = scratch[5].r + scratch[4].i;
            Fout[m3].i = scratch[5].i - scratch[4].r;
        }else{
            Fout[m].r = scratch[5].r + scratch[4].i;
            Fout[m].i = scratch[5].i - scratch[4].r;
            Fout[m3].r = scratch[5].r - scratch[4].i;
            Fout[m3].i = scratch[5].i + scratch[4].r;
        }
        ++Fout;
    }while(--k);
}

static void kf_bfly3(
         kiss_fft_cpx * Fout,
         const size_t fstride,
         const kiss_fft_cfg st,
         size_t m
         )
{
     size_t k=m;
     const size_t m2 = 2*m;
     kiss_fft_cpx *tw1,*tw2;
     kiss_fft_cpx scratch[5];
     kiss_fft_cpx epi3;
     epi3 = st->twiddles[fstride*m];

     tw1=tw2=st->twiddles;

     do{
         C_FIXDIV(*Fout,3); C_FIXDIV(Fout[m],3); C_FIXDIV(Fout[m2],3);

         C_MUL(scratch[1],Fout[m] , *tw1);
         C_MUL(scratch[2],Fout[m2] , *tw2);

         C_ADD(scratch[3],scratch[1],scratch[2]);
         C_SUB(scratch[0],scratch[1],scratch[2]);
         tw1 += fstride;
         tw2 += fstride*2;

         Fout[m].r = Fout->r - HALF_OF(scratch[3].r);
         Fout[m].i = Fout->i - HALF_OF(scratch[3].i);

         C_MULBYSCALAR( scratch[0] , epi3.i );

         C_ADDTO(*Fout,scratch[3]);

         Fout[m2].r = Fout[m].r + scratch[0].i;
         Fout[m2].i = Fout[m].i - scratch[0].r;

         Fout[m].r -= scratch[0].i;
         Fout[m].i += scratch[0].r;

         ++Fout;
     }while(--k);
}

static void kf_bfly5(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m
        )
{
    kiss_fft_cpx *Fout0,*Fout1,*Fout2,*Fout3,*Fout4;
    int u;
    kiss_fft_cpx scratch[13];
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx *tw;
    kiss_fft_cpx ya,yb;
    ya = twiddles[fstride*m];
    yb = twiddles[fstride*2*m];

    Fout0=Fout;
    Fout1=Fout0+m;
    Fout2=Fout0+2*m;
    Fout3=Fout0+3*m;
    Fout4=Fout0+4*m;

    tw=st->twiddles;
    for ( u=0; u<m; ++u ) {
        C_FIXDIV( *Fout0,5); C_FIXDIV( *Fout1,5); C_FIXDIV( *Fout2,5); C_FIXDIV( *Fout3,5); C_FIXDIV( *Fout4,5);
        scratch[0] = *Fout0;

        C_MUL(scratch[1] ,*Fout1, tw[u*fstride]);
        C_MUL(scratch[2] ,*Fout2, tw[2*u*fstride]);
        C_MUL(scratch[3] ,*Fout3, tw[3*u*fstride]);
        C_MUL(scratch[4] ,*Fout4, tw[4*u*fstride]);

        C_ADD( scratch[7],scratch[1],scratch[4]);
        C_SUB( scratch[10],scratch[1],scratch[4]);
        C_ADD( scratch[8],scratch[2],scratch[3]);
        C_SUB( scratch[9],scratch[2],scratch[3]);

        Fout0->r += scratch[7].r + scratch[8].r;
        Fout0->i += scratch[7].i + scratch[8].i;

        scratch[5].r = scratch[0].r + S_MUL(scratch[7].r,ya.r) + S_MUL(scratch[8].r,yb.r);
        scratch[5].i = scratch[0].i + S_MUL(scratch[7].i,ya.r) + S_MUL(scratch[8].i,yb.r);

        scratch[6].r =  S_MUL(scratch[10].i,ya.i) + S_MUL(scratch[9].i,yb.i);
        scratch[6].i = -S_MUL(scratch[10].r,ya.i) - S_MUL(scratch[9].r,yb.i);

        C_SUB(*Fout1,scratch[5],scratch[6]);
        C_ADD(*Fout4,scratch[5],scratch[6]);

        scratch[11].r = scratch[0].r + S_MUL(scratch[7].r,yb.r) + S_MUL(scratch[8].r,ya.r);
        scratch[11].i = scratch[0].i + S_MUL(scratch[7].i,yb.r) + S_MUL(scratch[8].i,ya.r);
        scratch[12].r = - S_MUL(scratch[10].i,yb.i) + S_MUL(scratch[9].i,ya.i);
        scratch[12].i = S_MUL(scratch[10].r,yb.i) - S_MUL(scratch[9].r,ya.i);

        C_ADD(*Fout2,scratch[11],scratch[12]);
        C_SUB(*Fout3,scratch[11],scratch[12]);

        ++Fout0;++Fout1;++Fout2;++Fout3;++Fout4;
    }
}

/* perform the butterfly for one stage of a mixed radix FFT */
static void kf_bfly_generic(
        kiss_fft_cpx * Fout,
        const size_t fstride,
        const kiss_fft_cfg st,
        int m,
        int p
        )
{
    int u,k,q1,q;
    kiss_fft_cpx * twiddles = st->twiddles;
    kiss_fft_cpx t;
    int Norig = st->nfft;

    CHECKBUF(scratchbuf,nscratchbuf,p);

    for ( u=0; u<m; ++u ) {
        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            scratchbuf[q1] = Fout[ k  ];
            C_FIXDIV(scratchbuf[q1],p);
            k += m;
        }

        k=u;
        for ( q1=0 ; q1<p ; ++q1 ) {
            int twidx=0;
            Fout[ k ] = scratchbuf[0];
            for (q=1;q<p;++q ) {
                twidx += fstride * k;
                if (twidx>=Norig) twidx-=Norig;
                C_MUL(t,scratchbuf[q] , twiddles[twidx] );
                C_ADDTO( Fout[ k ] ,t);
            }
            k += m;
        }
    }
}

static
void kf_work(
        kiss_fft_cpx * Fout,
        const kiss_fft_cpx * f,
        const size_t fstride,
        int in_stride,
        int * factors,
        const kiss_fft_cfg st
        )
{
    kiss_fft_cpx * Fout_beg=Fout;
    const int p=*factors++; /* the radix  */
    const int m=*factors++; /* stage's fft length/p */
    const kiss_fft_cpx * Fout_end = Fout + p*m;

    if (m==1) {
        do{
            *Fout = *f;
            f += fstride*in_stride;
        }while(++Fout != Fout_end );
    }else{
        do{
            kf_work( Fout , f, fstride*p, in_stride, factors,st);
            f += fstride*in_stride;
        }while( (Fout += m) != Fout_end );
    }

    Fout=Fout_beg;

    switch (p) {
        case 2: kf_bfly2(Fout,fstride,st,m); break;
        case 3: kf_bfly3(Fout,fstride,st,m); break;
        case 4: kf_bfly4(Fout,fstride,st,m); break;
        case 5: kf_bfly5(Fout,fstride,st,m); break;
        default: kf_bfly_generic(Fout,fstride,st,m,p); break;
    }
}

/*  facbuf is populated by p1,m1,p2,m2, ...
    where
    p[i] * m[i] = m[i-1]
    m0 = n                  */
static
void kf_factor(int n,int * facbuf)
{
    int p=4;
    double floor_sqrt;
    floor_sqrt = floor( sqrt((double)n) );

    /*factor out powers of 4, powers of 2, then any remaining primes */
    do {
        while (n % p) {
            switch (p) {
                case 4: p = 2; break;
                case 2: p = 3; break;
                default: p += 2; break;
            }
            if (p > floor_sqrt)
                p = n;          /* no more factors, skip to end */
        }
        n /= p;
        *facbuf++ = p;
        *facbuf++ = n;
    } while (n > 1);
}

/*
 *
 * User-callable function to allocate all necessary storage space for the fft.
 *
 * The return value is a contiguous block of memory, allocated with malloc.  As such,
 * It can be freed with free(), rather than a kiss_fft-specific function.
 * */
kiss_fft_cfg kiss_fft_alloc(int nfft,int inverse_fft,void * mem,size_t * lenmem )
{
    kiss_fft_cfg st=NULL;
    size_t memneeded = sizeof(struct kiss_fft_state)
        + sizeof(kiss_fft_cpx)*(nfft-1); /* twiddle factors*/

    if ( lenmem==NULL ) {
        st = ( kiss_fft_cfg)KISS_FFT_MALLOC( memneeded );
    }else{
        if (*lenmem >= memneeded)
            st = (kiss_fft_cfg)mem;
        *lenmem = memneeded;
    }
    if (st) {
        int i;
        st->nfft=nfft;
        st->inverse = inverse_fft;

        for (i=0;i<nfft;++i) {
            const double pi=3.14159265358979323846264338327;
            double phase = ( -2*pi /nfft ) * i;
            if (st->inverse)
                phase *= -1;
            kf_cexp(st->twiddles+i, phase );
        }

        kf_factor(nfft,st->factors);
    }
    return st;
}

static void kiss_fft_stride(kiss_fft_cfg st,const kiss_fft_cpx *fin,kiss_fft_cpx *fout,int in_stride)
{
    if (fin == fout) {
        CHECKBUF(tmpbuf,ntmpbuf,st->nfft);
        kf_work(tmpbuf,fin,1,in_stride, st->factors,st);
        memcpy(fout,tmpbuf,sizeof(kiss_fft_cpx)*st->nfft);
    }else{
        kf_work( fout, fin, 1,in_stride, st->factors,st );
    }
}

static void kiss_fft(kiss_fft_cfg cfg,const kiss_fft_cpx *fin,kiss_fft_cpx *fout) {
    kiss_fft_stride(cfg,fin,fout,1);
}


/* not really necessary to call, but if someone is doing in-place ffts, they may want to free the
   buffers from CHECKBUF
 */
static void kiss_fft_cleanup(void) {
    free(scratchbuf);
    scratchbuf = NULL;
    nscratchbuf=0;
    free(tmpbuf);
    tmpbuf=NULL;
    ntmpbuf=0;
}

#endif

#ifdef USE_FFMPEG
typedef float FFTSample;

typedef struct FFTComplex {
    FFTSample re, im;
} FFTComplex;

typedef struct FFTContext {
    int nbits;
    int inverse;
    uint16_t *revtab;
    FFTComplex *exptab;
    FFTComplex *exptab1; /* only used by SSE code */
    void (*fft_calc)(struct FFTContext *s, FFTComplex *z);
} FFTContext;

/**
 * The size of the FFT is 2^nbits. If inverse is TRUE, inverse FFT is
 * done
 */
void ff_fft_calc_c(FFTContext *s, FFTComplex *z);

#if defined(HAVE_MMX)
#include <xmmintrin.h>

static const float p1p1p1m1[4] __attribute__((aligned(16))) =
    { 1.0, 1.0, 1.0, -1.0 };

static const float p1p1m1p1[4] __attribute__((aligned(16))) =
    { 1.0, 1.0, -1.0, 1.0 };

static const float p1p1m1m1[4] __attribute__((aligned(16))) =
    { 1.0, 1.0, -1.0, -1.0 };

void ff_fft_calc_sse(FFTContext *s, FFTComplex *z)
{
    int ln = s->nbits;
    int	j, np, np2;
    int	nblocks, nloops;
    register FFTComplex *p, *q;
    FFTComplex *cptr, *cptr1;
    int k;

    np = 1 << ln;

    {
        __m128 *r, a, b, a1, c1, c2;

        r = (__m128 *)&z[0];
        c1 = *(__m128 *)p1p1m1m1;
        c2 = *(__m128 *)p1p1p1m1;
        if (s->inverse)
            c2 = *(__m128 *)p1p1m1p1;
        else
            c2 = *(__m128 *)p1p1p1m1;

        j = (np >> 2);
        do {
            a = r[0];
            b = _mm_shuffle_ps(a, a, _MM_SHUFFLE(1, 0, 3, 2));
            a = _mm_mul_ps(a, c1);
            /* do the pass 0 butterfly */
            a = _mm_add_ps(a, b);

            a1 = r[1];
            b = _mm_shuffle_ps(a1, a1, _MM_SHUFFLE(1, 0, 3, 2));
            a1 = _mm_mul_ps(a1, c1);
            /* do the pass 0 butterfly */
            b = _mm_add_ps(a1, b);

            /* multiply third by -i */
            b = _mm_shuffle_ps(b, b, _MM_SHUFFLE(2, 3, 1, 0));
            b = _mm_mul_ps(b, c2);

            /* do the pass 1 butterfly */
            r[0] = _mm_add_ps(a, b);
            r[1] = _mm_sub_ps(a, b);
            r += 2;
        } while (--j != 0);
    }
    /* pass 2 .. ln-1 */

    nblocks = np >> 3;
    nloops = 1 << 2;
    np2 = np >> 1;

    cptr1 = s->exptab1;
    do {
        p = z;
        q = z + nloops;
        j = nblocks;
        do {
            cptr = cptr1;
            k = nloops >> 1;
            do {
                __m128 a, b, c, t1, t2;

                a = *(__m128 *)p;
                b = *(__m128 *)q;

                /* complex mul */
                c = *(__m128 *)cptr;
                /*  cre*re cim*re */
                t1 = _mm_mul_ps(c,
                                _mm_shuffle_ps(b, b, _MM_SHUFFLE(2, 2, 0, 0)));
                c = *(__m128 *)(cptr + 2);
                /*  -cim*im cre*im */
                t2 = _mm_mul_ps(c,
                                _mm_shuffle_ps(b, b, _MM_SHUFFLE(3, 3, 1, 1)));
                b = _mm_add_ps(t1, t2);

                /* butterfly */
                *(__m128 *)p = _mm_add_ps(a, b);
                *(__m128 *)q = _mm_sub_ps(a, b);

                p += 2;
                q += 2;
                cptr += 4;
            } while (--k);

            p += nloops;
            q += nloops;
        } while (--j);
        cptr1 += nloops * 2;
        nblocks = nblocks >> 1;
        nloops = nloops << 1;
    } while (nblocks != 0);
}
#endif

int ff_fft_init(FFTContext *s, int nbits, int inverse)
{
    int i, j, m, n;
    float alpha, c1, s1, s2;

    s->nbits = nbits;
    n = 1 << nbits;

    s->exptab = KISS_FFT_MALLOC((n / 2) * sizeof(FFTComplex));
    s->revtab = KISS_FFT_MALLOC(n * sizeof(uint16_t));
    s->inverse = inverse;

    s2 = inverse ? 1.0 : -1.0;

    for(i=0;i<(n/2);i++) {
        alpha = 2 * M_PI * (float)i / (float)n;
        c1 = cos(alpha);
        s1 = sin(alpha) * s2;
        s->exptab[i].re = c1;
        s->exptab[i].im = s1;
    }
    s->fft_calc = ff_fft_calc_c;
    s->exptab1 = NULL;

    /* compute constant table for HAVE_SSE version */
#if defined(HAVE_MMX)
    {
            int np, nblocks, np2, l;
            FFTComplex *q;

            np = 1 << nbits;
            nblocks = np >> 3;
            np2 = np >> 1;
            s->exptab1 = KISS_FFT_MALLOC(np * 2 * sizeof(FFTComplex));
            q = s->exptab1;
            do {
                for(l = 0; l < np2; l += 2 * nblocks) {
                    *q++ = s->exptab[l];
                    *q++ = s->exptab[l + nblocks];

                    q->re = -s->exptab[l].im;
                    q->im = s->exptab[l].re;
                    q++;
                    q->re = -s->exptab[l + nblocks].im;
                    q->im = s->exptab[l + nblocks].re;
                    q++;
                }
                nblocks = nblocks >> 1;
            } while (nblocks != 0);
            free(s->exptab); s->exptab = NULL;
            s->fft_calc = ff_fft_calc_sse;
    }
#endif

    /* compute bit reverse table */

    for(i=0;i<n;i++) {
        m=0;
        for(j=0;j<nbits;j++) {
            m |= ((i >> j) & 1) << (nbits-j-1);
        }
        s->revtab[i]=m;
    }
    return 0;
}

/* butter fly op */
#define BF(pre, pim, qre, qim, pre1, pim1, qre1, qim1) \
{\
  FFTSample ax, ay, bx, by;\
  bx=pre1;\
  by=pim1;\
  ax=qre1;\
  ay=qim1;\
  pre = (bx + ax);\
  pim = (by + ay);\
  qre = (bx - ax);\
  qim = (by - ay);\
}

#define MUL16(a,b) ((a) * (b))

#define CMUL(pre, pim, are, aim, bre, bim) \
{\
   pre = (MUL16(are, bre) - MUL16(aim, bim));\
   pim = (MUL16(are, bim) + MUL16(bre, aim));\
}

/**
 * Do a complex FFT with the parameters defined in ff_fft_init(). The
 * input data must be permuted before with s->revtab table. No
 * 1.0/sqrt(n) normalization is done.
 */
void ff_fft_calc_c(FFTContext *s, FFTComplex *z)
{
    int ln = s->nbits;
    int	j, np, np2;
    int	nblocks, nloops;
    register FFTComplex *p, *q;
    FFTComplex *exptab = s->exptab;
    int l;
    FFTSample tmp_re, tmp_im;

    np = 1 << ln;

    /* pass 0 */

    p=&z[0];
    j=(np >> 1);
    do {
        BF(p[0].re, p[0].im, p[1].re, p[1].im,
           p[0].re, p[0].im, p[1].re, p[1].im);
        p+=2;
    } while (--j != 0);

    /* pass 1 */


    p=&z[0];
    j=np >> 2;
    if (s->inverse) {
        do {
            BF(p[0].re, p[0].im, p[2].re, p[2].im,
               p[0].re, p[0].im, p[2].re, p[2].im);
            BF(p[1].re, p[1].im, p[3].re, p[3].im,
               p[1].re, p[1].im, -p[3].im, p[3].re);
            p+=4;
        } while (--j != 0);
    } else {
        do {
            BF(p[0].re, p[0].im, p[2].re, p[2].im,
               p[0].re, p[0].im, p[2].re, p[2].im);
            BF(p[1].re, p[1].im, p[3].re, p[3].im,
               p[1].re, p[1].im, p[3].im, -p[3].re);
            p+=4;
        } while (--j != 0);
    }
    /* pass 2 .. ln-1 */

    nblocks = np >> 3;
    nloops = 1 << 2;
    np2 = np >> 1;
    do {
        p = z;
        q = z + nloops;
        for (j = 0; j < nblocks; ++j) {
            BF(p->re, p->im, q->re, q->im,
               p->re, p->im, q->re, q->im);

            p++;
            q++;
            for(l = nblocks; l < np2; l += nblocks) {
                CMUL(tmp_re, tmp_im, exptab[l].re, exptab[l].im, q->re, q->im);
                BF(p->re, p->im, q->re, q->im,
                   p->re, p->im, tmp_re, tmp_im);
                p++;
                q++;
            }

            p += nloops;
            q += nloops;
        }
        nblocks = nblocks >> 1;
        nloops = nloops << 1;
    } while (nblocks != 0);
}

/**
 * Do the permutation needed BEFORE calling ff_fft_calc()
 */
void ff_fft_permute(FFTContext *s, FFTComplex *z)
{
    int j, k, np;
    FFTComplex tmp;
    const uint16_t *revtab = s->revtab;

    /* reverse */
    np = 1 << s->nbits;
    for(j=0;j<np;j++) {
        k = revtab[j];
        if (k < j) {
            tmp = z[k];
            z[k] = z[j];
            z[j] = tmp;
        }
    }
}

void ff_fft_end(FFTContext *s)
{
    free(s->revtab);
    free(s->exptab);
    free(s->exptab1);
}

#endif

#define PRECISION 1

typedef struct kiss_mdct_s {
	int mdct_size;
	//Rotation tables
	kiss_fft_cpx * mdct_cs;
	//tmp buffers
	kiss_fft_cpx * mdct_tmp_buf_in;
	kiss_fft_cpx * mdct_tmp_buf_out;
#ifdef USE_FFTW
	fftwf_plan plan;
#else
#ifdef USE_FFMPEG
	FFTContext fc;
#else
	kiss_fft_cfg mdct_fft_cfg;
#endif
#endif
} kiss_mdct_t;

kiss_mdct_t * kiss_imdct_init(int n) {
	//The code is based on the algorithm in liba52 specs page 100
	//N/4 fft with pre and post rotation steps
	kiss_mdct_t * s;
	int i;
	int fft_size = n/4;
	float angle;

	s = malloc(sizeof(kiss_mdct_t));
	s->mdct_size = n;

	s->mdct_cs = KISS_FFT_MALLOC(fft_size * sizeof(kiss_fft_cpx));
	s->mdct_tmp_buf_in = KISS_FFT_MALLOC(fft_size * sizeof(kiss_fft_cpx));
	s->mdct_tmp_buf_out = KISS_FFT_MALLOC(fft_size * sizeof(kiss_fft_cpx));

	for(i = 0; i < fft_size; i++) {
		angle = 2 * M_PI * (i + 0.125) / s->mdct_size;
		s->mdct_cs[i].r = -cosf(angle);
		s->mdct_cs[i].i = -sinf(angle);
	}

#ifdef USE_FFTW
	s->plan = fftwf_plan_dft_1d(fft_size, (void*)s->mdct_tmp_buf_in, (void*)s->mdct_tmp_buf_out, FFTW_BACKWARD, FFTW_MEASURE);
#else
#ifdef USE_FFMPEG
	for (i = 1; fft_size >> i; i++); i--;
	ff_fft_init(&s->fc, i, 1);
#else
	s->mdct_fft_cfg = kiss_fft_alloc(fft_size, 1, 0, 0);
#endif
#endif
	return s;
}

void kiss_imdct(kiss_mdct_t * s, kiss_fft_scalar * inbuffer, kiss_fft_scalar * outbuffer) {
	int k, size, size2, size4, size8;
	kiss_fft_cpx tmp;

	size = s->mdct_size;
	size2 = size >> 1;
	size4 = size >> 2;
	size8 = size >> 3;

	//FIXME the inbuffer reordering might be inefficient
	//C_MUL for this might not be the best thing to use
	//pre rotation
	for(k = 0; k < size4; k++) {
		tmp.r = inbuffer[size2-2*k-1];
		tmp.i = inbuffer[2*k];
		C_MUL(s->mdct_tmp_buf_in[k], tmp, s->mdct_cs[k]);
	}

#define TMPIN mdct_tmp_buf_in
#define TMPOUT mdct_tmp_buf_out

#ifdef USE_FFTW
	fftwf_execute(s->plan);
#else
#ifdef USE_FFMPEG
	ff_fft_permute(&s->fc, (void*)s->mdct_tmp_buf_in);
	s->fc.fft_calc(&s->fc, (void*)s->mdct_tmp_buf_in);
#undef TMPIN
#undef TMPOUT
#define TMPIN mdct_tmp_buf_out
#define TMPOUT mdct_tmp_buf_in
#else
	kiss_fft(s->mdct_fft_cfg, s->mdct_tmp_buf_in, s->mdct_tmp_buf_out);
#endif
#endif

	//reuse the buf_in buffer as the output buffer
	//post rotation
	for(k = 0; k < size4; k++){
		C_MUL(s->TMPIN[k], s->mdct_cs[k], s->TMPOUT[k]);
	}

	//S_MUL should be used here with the window
	//deinterleave, windowing should be done here also
	//maybe do it in a better order ?
	for(k = 0; k < size8; k++) {
		outbuffer[k*2]           = -s->TMPIN[size8+k].i;
		outbuffer[k*2+1]         =  s->TMPIN[size8-k-1].r;
		outbuffer[size4+k*2]     = -s->TMPIN[k].r;
		outbuffer[size4+k*2+1]   =  s->TMPIN[size4-k-1].i;
		outbuffer[size2+k*2]     = -s->TMPIN[size8+k].r;
		outbuffer[size2+k*2+1]   =  s->TMPIN[size8-k-1].i;
		outbuffer[3*size4+k*2]   =  s->TMPIN[k].i;
		outbuffer[3*size4+k*2+1] = -s->TMPIN[size4-k-1].r;
	}
}

void kiss_imdct_end(kiss_mdct_t * s) {
	if (!s) return;
	free(s->mdct_tmp_buf_out);
	free(s->mdct_tmp_buf_in);
	free(s->mdct_cs);
#ifdef USE_FFTW
	fftwf_destroy_plan(s->plan);
#else
#ifdef USE_FFMPEG
	ff_fft_end(&s->fc);
#else
	free(s->mdct_fft_cfg);
	kiss_fft_cleanup();
#endif
#endif
	free(s);
}
