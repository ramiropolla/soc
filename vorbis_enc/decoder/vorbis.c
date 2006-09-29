#define NEW_BOOK

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "vorbis.h"
#include "kiss_imdct.h"

#define MIN(a,b) ((a) > (b) ? (b) : (a))
#define MAX(a,b) ((a) < (b) ? (b) : (a))
#define ABS(a) ((a) > 0 ? (a) : -(a))

//#define debug_msg(...) fprintf(stderr, __VA_ARGS__)
//#define SLOW_DEBUG(...) __VA_ARGS__
#define debug_msg(...) ((void)0)
#define SLOW_DEBUG(...) /* nothing */
#define DEBUG_RATE 44100

#define malloc(x) (debug_msg("%d at %s:%d, %s()\n", (int)(x), __FILE__, __LINE__, __PRETTY_FUNCTION__), malloc(x))

#ifndef M_PI
# define M_PI           3.14159265358979323846
#endif

#define SMALL_BITS 6
#define SMALL_SIZE (1<<SMALL_BITS)

static inline long long read_time(void) {
        long long l;
        asm volatile("rdtsc\n\t":"=A"(l));
        return l;
}

#define START_TIMER \
uint64_t tend;\
uint64_t tstart= read_time();\

#define STOP_TIMER(id) \
tend= read_time();\
{\
  static uint64_t tsum=0;\
  static int tcount=0;\
  static int tskip_count=0;\
  if(tcount<2 || tend - tstart < 8*tsum/tcount){\
      tsum+= tend - tstart;\
      tcount++;\
  }else\
      tskip_count++;\
  if(256*256*256*64%(tcount+tskip_count)==0){\
      fprintf(stderr, "%"PRIu64" dezicycles in %s, %d runs, %d skips\n", tsum*10/tcount, id, tcount, tskip_count);\
  }\
}

typedef struct bit_packer_s {
	int index;
	int total;
	uint8_t * buf;
} bit_packer_t;

typedef struct cb_entry_s {
	int len;
	int real_num;
	uint32_t codeword;
} cb_entry_t;

typedef struct codebook_s {
	int nentries;
	int ndimentions;
	float * dimentions;
	cb_entry_t * entries;
	int max_len;
	int lookup;
	uint32_t small_table[SMALL_SIZE];
} codebook_t;

typedef struct floor0_context_s {
	int order;
	int rate;
	int bark_size;
	int amplitude_bits;
	int amplitude_offset;
	int nbooks;
	int * books;
	float * lsp;
	float * map0;
	int * map0_count;
	float * map1;
	int * map1_count;
} floor0_context_t;

typedef struct floor1_class_s {
	int dim;
	int subclass;
	int masterbook;
	int * books;
} floor1_class_t;

typedef struct floor1_list_s {
	int x;
	int y;
	int sort;
	int low;
	int high;
	int flag;
} floor1_list_t;

typedef struct floor1_context_s {
	int partitions;
	int multiplier;
	int rangebits;
	int * partition_to_class;
	floor1_list_t * list;
	int nclasses;
	floor1_class_t * classes;
	int values;
} floor1_context_t;

typedef struct residue_context_s {
	int type;
	int begin;
	int end;
	int partition_size;
	int classifications;
	int classbook;
	int (*books)[8];
} residue_context_t;

typedef struct floor_context_s {
	int type;
	union {
		floor0_context_t f0;
		floor1_context_t f1;
	} u;
} floor_context_t;

typedef struct mapping_context_s {
	int submaps;
	int * mux;
	int coupling_steps;
	int * magnitude;
	int * angle;
	int * floor;
	int * residue;
} mapping_context_t;

typedef struct mode_context_s {
	int blockflag;
	int mapping;
} mode_context_t;

typedef struct vorbis_context_s {
	int channels;
	int sample_rate;
	int blocksize0;
	int blocksize1;

	const float * win0;
	const float * win1;

	int ncodebooks;
	codebook_t * codebooks;
	int nfloors;
	floor_context_t * floors;
	int nresidues;
	residue_context_t * residues;
	int nmappings;
	mapping_context_t * mappings;
	int nmodes;
	mode_context_t * modes;

	// all of these are per channel
	float * floor_buf;
	float * residue_buf; // multiplied in place
	float * saved_buf;
	float * buf; // mdct output

	int saved;
	kiss_imdct_t * imdct0;
	kiss_imdct_t * imdct1;
} vorbis_context_t;

static float float32_unpack(uint32_t x) {
	float mant = x & 0x001fffff;
	int exp = (x & 0x7fe00000) >> 21;
	if (x & 0x80000000) mant *= -1;
	return ldexpf(mant, exp - 20 - 768);
}

static int ilog(int a) {
	int i;
	for (i = 0; (a >> i) > 0; i++);
	return i;
}

static void init_bp(bit_packer_t * bp, uint8_t * buf, int bytes) {
	bp->buf = buf;
	bp->total = bytes*8;
	bp->index = 0;
}

static inline uint32_t bswap_32(uint32_t x) {
	asm ("bswap %0" : "=r" (x) : "0" (x));
	return x;
}

static inline uint64_t bswap_64(uint64_t x) {
	union {
		uint64_t ll;
		struct { uint32_t l, h; } l;
	} r;
	r.l.h = bswap_32 (x);
	r.l.l = bswap_32 (x>>32);
	return r.ll;
}

#define bswap_32(x) (x)
#define bswap_64(x) (x)

static inline int get_bits(bit_packer_t * bp, int bits, uint32_t * res) {
	int pos = bp->index;
	bp->index += bits;
	if (bp->index > bp->total) return 1;
	if (!bits) { if (res) *res = 0; return 0; }
	if (res) *res = (bswap_64(*(uint64_t*)(bp->buf + (pos >> 3))) >> (pos & 7)) & ((1ULL << bits) - 1); // >> (64-bits)
	//fprintf(stderr, "%3d %2d 0x%016llX\n", *res, bits, bswap_64(*(uint64_t*)(bp->buf + (pos >> 3))));}
	return 0;
}

#define ERROR(x, n) do{ if (x) { err = n; goto err_out; } }while(0)
#define CHECK(x) do{ if ((err = (x))) goto err_out; }while(0)
#define GET_B(bp, i, x) do{ uint32_t _num; CHECK(get_bits(bp, i, &_num)); (x) = _num; }while(0)

static int cb_lookup_vals(int lookup, int dimentions, int entries) {
	int tmp, i;
	if (lookup == 1) {
		for (tmp = 0; ; tmp++) {
			int n = 1;
			for (i = 0; i < dimentions; i++) n *= tmp;
			if (n > entries) break;
		}
		tmp--;
	} else tmp = dimentions * entries;
	return tmp;
}

static void uninit_codebook(codebook_t * cb) {
	if (cb->nentries == -1) return;

	free(cb->dimentions);
	free(cb->entries);
	cb->entries = NULL;
	cb->nentries = -1;
}
static uint32_t FLIP_SMALL(uint32_t a) {
	int i;
	uint32_t tmp = 0;
	for (i = 0; i < SMALL_BITS; i++) tmp |= ((a >> i) & 1) << (SMALL_BITS - i - 1);
	return tmp;
}
static int read_codebook_header(vorbis_context_t * vc, codebook_t * cb, bit_packer_t * bp) {
	int tmp;
	int h[33] = { 1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
	int i;
	int err;
	int tot_entries;

	cb->entries = NULL;
	cb->dimentions = NULL;

	GET_B(bp, 24, tmp);
	ERROR(tmp != 0x564342, 5);
	GET_B(bp, 16, cb->ndimentions);
	GET_B(bp, 24, cb->nentries);
	ERROR(cb->nentries == 0, 5); // empty codebook?
	tot_entries = cb->nentries;

	GET_B(bp, 1, tmp);
	if (tmp) { // ordered
		int len;
		cb->entries = malloc(sizeof(cb_entry_t) * cb->nentries);
		for (i = 0; i < cb->nentries; i++) cb->entries[i].real_num = i;
		i = 0;
		GET_B(bp, 5, len);
		len++;
		while (i < cb->nentries) {
			int j;
			int num;
			GET_B(bp, ilog(cb->nentries - i), num);
			for (j = i; j < i+num; j++) cb->entries[j].len = len;
			i += num;
			len++;
		}
		ERROR(i > cb->nentries, 6);
	} else { // not ordered
		GET_B(bp, 1, tmp);
		if (tmp) { // sparse
			int e[tot_entries];
			for (i = 0; i < tot_entries; i++) {
				GET_B(bp, 1, e[i]);
				if (!e[i]) cb->nentries--;
				else {
					GET_B(bp, 5, e[i]);
					e[i]++;
				}
			}
			cb->entries = malloc(sizeof(cb_entry_t) * cb->nentries);
			tmp = 0;
			for (i = 0; i < tot_entries; i++) {
				if (e[i]) {
					cb->entries[tmp].real_num = i;
					cb->entries[tmp].len = e[i];
					tmp++;
				}
			}
		} else { // not sparse
			cb->entries = malloc(sizeof(cb_entry_t) * cb->nentries);
			for (i = 0; i < cb->nentries; i++) {
				cb->entries[i].real_num = i;
				GET_B(bp, 5, cb->entries[i].len);
				cb->entries[i].len++;
			}
		}
	}

	// huffman
	cb->max_len = 0;
	for (i = 0; i < cb->nentries; i++) {
		cb_entry_t * e = &cb->entries[i];
		int j = 0;
		cb->max_len = MAX(cb->max_len, e->len);
		if (h[0]) h[0] = 0;
		else for (j = e->len; !h[j]; j--) ERROR(!j, 9);
		e->codeword = h[j];
		h[j] = 0;
		for (j++; j <= e->len; j++) h[j] = e->codeword | (1 << (j - 1));
	}
	for (i = 0; i < 33; i++) ERROR(h[i], 9);

#ifdef NEW_BOOK
	for (i = 0; i < cb->nentries; i++) {
		cb_entry_t * e = &cb->entries[i];
		int j, tmp = 0;

		for (j = 0; j < e->len; j++) tmp |= ((e->codeword >> j) & 1) << (e->len - j - 1);

		e->codeword = tmp << (32 - e->len);
	}
#endif

	for (i = 0; i < cb->nentries - 1; i++) {
		int j;
		for (j = i + 1; j < cb->nentries; j++) {
#ifdef NEW_BOOK
			if (cb->entries[i].codeword > cb->entries[j].codeword) {
#else
			if (cb->entries[i].len > cb->entries[j].len) {
#endif
				cb_entry_t tmp = cb->entries[i];
				cb->entries[i] = cb->entries[j];
				cb->entries[j] = tmp;
			}
		}
	}

	for (i = 0; i < SMALL_SIZE; i++) cb->small_table[i] = 0;
	for (i = 0; i < cb->nentries; i++) {
		cb_entry_t * e = &cb->entries[i];
		int j, tmp = 0;
#ifdef NEW_BOOK
		tmp = e->codeword >> (32 - e->len);
#else
		for (j = 0; j < e->len; j++) tmp |= ((e->codeword >> j) & 1) << (e->len - j - 1);
#endif
		if (e->len <= SMALL_BITS) {
			tmp <<= SMALL_BITS - e->len;
			for (j = 0; j < (1<<(SMALL_BITS - e->len)); j++)
				cb->small_table[FLIP_SMALL(tmp + j)] = i;
		} else {
			int cur, hi = 0, lo = cb->nentries;
			tmp >>= e->len - SMALL_BITS;
			cur = cb->small_table[FLIP_SMALL(tmp)];
			if (cur & 0x80000000U) {
				hi = cb->nentries - ((cur>>15) & 0x7fff);
				lo = cur & 0x7fff;
			}
			hi = cb->nentries - MAX(hi, i) - 1;
			lo = MIN(lo, i);
			cb->small_table[FLIP_SMALL(tmp)] = 0x80000000U | ((hi&0x7fff) << 15) | (lo&0x7fff);
		}
	}

	GET_B(bp, 4, cb->lookup);
	ERROR(cb->lookup > 2, 7);
	if (cb->lookup) {
		int vals = cb_lookup_vals(cb->lookup, cb->ndimentions, tot_entries);
		int multilicands[vals];
		int bits;
		int seq_p;
		uint32_t tmp2;
		float min, delta;
		float * d = malloc(sizeof(float) * cb->ndimentions * cb->nentries);

		cb->dimentions = d;

		GET_B(bp, 32, tmp2); min = float32_unpack(tmp2);
		GET_B(bp, 32, tmp2); delta = float32_unpack(tmp2);
		GET_B(bp, 4, bits);
		bits++;
		GET_B(bp, 1, seq_p);
		for (i = 0; i < vals; i++) GET_B(bp, bits, multilicands[i]);

		debug_msg("codebook %d lookup %d seq_p %d\n", cb - vc->codebooks, cb->lookup, seq_p);
		debug_msg("min %lf delta %lf\n", min, delta);
		for (i = 0; i < vals; i++) debug_msg("%d, ", multilicands[i]);
		debug_msg("\n");

		for (i = 0; i < cb->nentries; i++) {
			cb_entry_t * e = &cb->entries[i];
			float last = 0;
			int j;
			int div = 1;
			for (j = 0; j < cb->ndimentions; j++) {
				int off;
				int pos = e->real_num;
				if (cb->lookup == 1) off = (pos / div) % vals; // lookup type 1
				else off = pos * cb->ndimentions + j; // lookup type 2

				d[i * cb->ndimentions + j] = last + min + multilicands[off] * delta;
				if (seq_p) last = d[i * cb->ndimentions + j];
				div *= vals;
			}
		}
	}
#if 0
	debug_msg("codebook %d %p\n", cb - vc->codebooks, cb);
	for (i = 0; i < cb->nentries; i++) {
		cb_entry_t * e = &cb->entries[i];
		int j;
		//if (!cb->lookup) continue;
		debug_msg("%4d %2d   ", e->real_num, e->len);
#ifdef NEW_BOOK
		for (j = 0; j < e->len; j++) debug_msg(((e->codeword >> (31 - j))&1) ? "1" : "0");
#else
		for (j = 0; j < e->len; j++) debug_msg(((e->codeword >> j)&1) ? "1" : "0");
#endif
		for (j = e->len; j < 21; j++) debug_msg(" ");
		if (cb->lookup) {
			float * d = cb->dimentions;
			debug_msg(" { ");
			for (j = 0; j < cb->ndimentions; j++) debug_msg("%+8.3f, ", d[i * cb->ndimentions + j]);
			debug_msg("}");
		}
		debug_msg("\n");
	}
#endif
	return 0;
err_out:
	uninit_codebook(cb);
	return err;
}

static uint8_t flip_bits8(uint8_t x) {
	x = ((x>> 4)&0x0f) | ((x<< 4)&0xf0);
	x = ((x>> 2)&0x33) | ((x<< 2)&0xcc);
	x = ((x>> 1)&0x55) | ((x<< 1)&0xaa);
	return x;
}

static inline uint32_t FLIP_N(uint32_t a, int n) {
	int i;
	uint32_t tmp = 0;
	for (i = 0; i < n; i++) tmp |= ((a >> i) & 1) << (n - i - 1);
	return tmp;
}
static inline int peek_bits(bit_packer_t * bp, int bits, uint32_t * res) {
	int pos = bp->index;
	if (pos + bits > bp->total) return 1;
	*res = FLIP_N((bswap_64(*(uint64_t*)(bp->buf + (pos >> 3))) >> (pos & 7)) & ((1ULL << bits) - 1), bits);
	return 0;
}
static inline void skip_bits(bit_packer_t * bp, int bits) {
	bp->index += bits;
}

static inline int get_book(bit_packer_t * bp, codebook_t * cb, int * res) {
#ifndef NEW_BOOK
	uint32_t tmp;
	int i;
	int j = 0;
	static int miss, hit;

	if (bp->total - bp->index >= SMALL_BITS) {
		peek_bits(bp, SMALL_BITS, &tmp);
		tmp = cb->small_table[tmp & (SMALL_SIZE-1)];
		if (tmp & 0x80000000U) {
			miss++;
			if (256*256*256*64 % miss == 0) debug_msg("%d %d\n", miss, hit);
		} else {
			hit++;
			skip_bits(bp, cb->entries[tmp].len);
			*res = tmp;
			return 0;
		}
	}

	tmp = 0;

	for (i = 1; i <= cb->max_len; i++) {
		uint32_t a;
		if (get_bits(bp, 1, &a)) return 1;
		//if (--bp->left < 0) return 1;
		tmp |= a << (i - 1);
		//if (++bp->pos == 8) { bp->buf_ptr++; bp->pos = 0; }

		for (; j < cb->nentries && cb->entries[j].len == i; j++) {
			if (cb->entries[j].codeword != tmp) continue;
			*res = j;
			return 0;
		}
	}
	return 1;
#else
	uint32_t tmp;
	int read;
	int lo = 0, hi = cb->nentries;
	//static int miss, hit;

	if (bp->total - bp->index >= SMALL_BITS) {
		tmp = bswap_32(*(uint32_t*)(bp->buf + (bp->index >> 3))) >> (bp->index & 7);
		//tmp = *bp->buf_ptr >> bp->pos;
		//if ((8 - bp->pos) < SMALL_BITS) tmp |= bp->buf_ptr[1] << (8 - bp->pos);
		//printf("%d ", tmp & (SMALL_SIZE-1));
		//peek_bits(bp, SMALL_BITS, &tmp);
		//printf("%d\n", tmp & (SMALL_SIZE-1));
		tmp = cb->small_table[tmp & (SMALL_SIZE-1)];
		if (tmp & 0x80000000U) {
			//miss++;
			//if (256*256*256*64 % miss == 0) debug_msg("%d %d\n", miss, hit);
			hi = cb->nentries - ((tmp>>15) & 0x7fff);
			lo = tmp & 0x7fff;
		} else {
			//hit++;
			skip_bits(bp, cb->entries[tmp].len);
			*res = tmp;
			return 0;
		}
	}

	read = MIN(bp->total - bp->index, cb->max_len);
	peek_bits(bp, read, &tmp);

	tmp <<= 32 - read;
	while (hi - lo > 1) {
		long p = (hi - lo) >> 1;
		long test = cb->entries[lo+p].codeword > tmp;
		lo += p & (test - 1);
		hi -= p & (-test);
	}
	if (cb->entries[lo].len > read) return 1;
	skip_bits(bp, cb->entries[lo].len);
	*res = lo;
	return 0;
#endif
}

#define GET_SCALAR(bp, i, x) do{ \
	int _num; \
	CHECK(get_book(bp, &vc->codebooks[i], &_num)); \
	(x) = vc->codebooks[i].entries[_num].real_num; \
}while(0)
#define GET_VECTOR(bp, i, x) do{ \
	float * d = vc->codebooks[i].dimentions; \
	CHECK(get_book(bp, &vc->codebooks[i], &_num)); \
	(x) = &d[_num * vc->codebooks[i].ndimentions]; \
}while(0)

static void uninit_floor(floor_context_t * fcc) {
	if (fcc->type == -1) return;

	if (fcc->type == 1) {
		floor1_context_t * fc = &fcc->u.f1;
		int i;
		if (fc->classes) for (i = 0; i < fc->nclasses; i++) free(fc->classes[i].books);
		free(fc->partition_to_class);
		free(fc->classes);
		free(fc->list);
	} else {
		floor0_context_t * fc = &fcc->u.f0;
		free(fc->books);
		free(fc->lsp);
		free(fc->map0);
		free(fc->map0_count);
		free(fc->map1);
		free(fc->map1_count);
	}
	fcc->type = -1;
}

static void fill_floor0_map(vorbis_context_t * vc, floor_context_t * fcc, int blocktype) {
	floor0_context_t * fc = &fcc->u.f0;
	int * map_count = blocktype ? fc->map1_count : fc->map0_count;

	if (map_count) return;
	else {
		int n = (blocktype ? vc->blocksize1 : vc->blocksize0)/2;
		int tmp_map[n];
		int i;
		int count = 1;
		float * map;
		float wstep = M_PI / fc->bark_size;

		for (i = 0; i < n; i++) {
#define BARK(x) (13.1f*atan(0.00074f*(x))+2.24f*atan(1.85e-8f*(x)*(x))+1e-4f*(x))
			tmp_map[i] = floor(BARK((fc->rate * i) / (2. * n)) * fc->bark_size / BARK(0.5 * fc->rate));
			tmp_map[i] = MIN(tmp_map[i], fc->bark_size - 1);
			if (i && tmp_map[i] != tmp_map[i-1]) count++;
		}

		map = malloc(count * sizeof(float));
		map_count = malloc((count+1) * sizeof(int));
		if (blocktype) {
			fc->map1 = map;
			fc->map1_count = map_count;
		} else {
			fc->map0 = map;
			fc->map0_count = map_count;
		}
		count = 0;
		for (i = 0; i < n; i++) {
			if (i && tmp_map[i] != tmp_map[i-1]) {
				*map++ = 2.f * cos(wstep * tmp_map[i-1]);
				*map_count++ = count;
				count = 1;
			} else count++;
		}
		*map++ = 2.f * cosf(wstep * tmp_map[i-1]);
		*map_count++ = count;
		*map_count++ = 0;
	}
}

static int read_floor0_header(vorbis_context_t * vc, floor_context_t * fcc, bit_packer_t * bp) {
	floor0_context_t * fc = &fcc->u.f0;
	int i;
	int err;
	int max_dim = 0;

	fc->books = NULL;
	fc->lsp = NULL;
	fc->map0 = NULL;
	fc->map0_count = NULL;
	fc->map1 = NULL;
	fc->map1_count = NULL;

	GET_B(bp, 8, fc->order);
	GET_B(bp, 16, fc->rate);
	GET_B(bp, 16, fc->bark_size);
	GET_B(bp, 6, fc->amplitude_bits);
	GET_B(bp, 8, fc->amplitude_offset);
	GET_B(bp, 4, fc->nbooks);
	fc->nbooks++;

	fc->books = malloc(fc->nbooks * sizeof(int));
	for (i = 0; i < fc->nbooks; i++) {
		GET_B(bp, 8, fc->books[i]);
		ERROR(fc->books[i] >= vc->ncodebooks, 21);
		max_dim = MAX(max_dim, vc->codebooks[fc->books[i]].ndimentions);
	}

	fc->lsp = malloc((fc->order + max_dim - 1) * sizeof(float));
	fill_floor0_map(vc, fcc, 0);
	fill_floor0_map(vc, fcc, 1);

	return 0;
err_out:
	uninit_floor(fcc);
	return err;
}

static int floor0_decode(vorbis_context_t * vc, floor_context_t * fcc, bit_packer_t * bp, float * buf, int blocktype) {
	floor0_context_t * fc = &fcc->u.f0;
	unsigned int n = (blocktype ? vc->blocksize1 : vc->blocksize0)/2;
	float * map = blocktype ? fc->map1 : fc->map0;
	int * map_count = blocktype ? fc->map1_count : fc->map0_count;
	int i;
	int err = 0;
	int amp;
	float last = 0.;
	int book, dim;
	float amplitude;

	GET_B(bp, fc->amplitude_bits, amp);
	if (!amp) { // no decode
		memset(buf, 0, sizeof(float)*n);
		return -1;
	}
	amplitude = (double)amp * fc->amplitude_offset / ((1<<fc->amplitude_bits) - 1);

	GET_B(bp, ilog(fc->nbooks), book);
	ERROR(book >= fc->nbooks, 23);
	book = fc->books[book];
	dim = vc->codebooks[book].ndimentions;
	i = 0;
	while (i < fc->order) {
		int _num;
		int j;
		float * tmp;
		GET_VECTOR(bp, book, tmp);
		for (j = 0; j < dim; j++, i++) fc->lsp[i] = 2.f * cos(tmp[j] + last);
		last = tmp[j-1] + last;
	}
	while (*map_count) {
		float p = 0.5f, q = 0.5f;
		float cos_w = *map++;

		for (i = 0; i < fc->order-1; i += 2) {
			q *= (fc->lsp[i]   - cos_w);
			p *= (fc->lsp[i+1] - cos_w);
		}
		if (i == fc->order) { // even
			p *= p * (2.f - cos_w);
			q *= q * (2.f + cos_w);
		} else { // odd
			q *= (cos_w - fc->lsp[i]);

			p *= p * (4.f - cos_w * cos_w);
			q *= q;
		}
		q = exp((amplitude / sqrt(p + q) - fc->amplitude_offset) * 0.11512925f);
		i = *map_count++;
		while (i--) *buf++ = q;
	}
err_out:
	return err;
}

static int read_floor1_header(vorbis_context_t * vc, floor_context_t * fcc, bit_packer_t * bp) {
	floor1_context_t * fc = &fcc->u.f1;
	int tmp;
	int i;
	int err;

	fc->partition_to_class = NULL;
	fc->classes = NULL;
	fc->list = NULL;

	GET_B(bp, 5, fc->partitions);
	fc->partition_to_class = malloc(sizeof(int) * fc->partitions);
	tmp = -1;
	for (i = 0; i < fc->partitions; i++) {
		GET_B(bp, 4, fc->partition_to_class[i]);
		tmp = MAX(tmp, fc->partition_to_class[i]);
		debug_msg("[%d] partition-class %d\n", i, fc->partition_to_class[i]);
	}
	fc->nclasses = tmp + 1;
	fc->classes = malloc(sizeof(floor1_class_t) * fc->nclasses);
	for (i = 0; i < fc->nclasses; i++) fc->classes[i].books = NULL;

	for (i = 0; i < fc->nclasses; i++) {
		int j, books;
		GET_B(bp, 3, fc->classes[i].dim);
		fc->classes[i].dim++;
		debug_msg("[%d] dim %d\n", i, fc->classes[i].dim);
		GET_B(bp, 2, fc->classes[i].subclass);
		if (fc->classes[i].subclass) {
			GET_B(bp, 8, fc->classes[i].masterbook);
			ERROR(fc->classes[i].masterbook >= vc->ncodebooks, 11);
			debug_msg("masterbook[%d] %d\n", i, fc->classes[i].masterbook);
		}

		books = (1 << fc->classes[i].subclass);
		fc->classes[i].books = malloc(sizeof(int) * books);
		for (j = 0; j < books; j++) {
			GET_B(bp, 8, fc->classes[i].books[j]);
			fc->classes[i].books[j]--;
			ERROR(fc->classes[i].books[j] >= vc->ncodebooks, 11);
			debug_msg("classbook[%d][%d] %d\n", i, j, fc->classes[i].books[j]);
		}
	}
	GET_B(bp, 2, fc->multiplier);
	fc->multiplier++;
	GET_B(bp, 4, fc->rangebits);
	debug_msg("multiplier %d, rangebits %d\n", fc->multiplier, fc->rangebits);

	fc->values = 2;
	for (i = 0; i < fc->partitions; i++)
		fc->values += fc->classes[fc->partition_to_class[i]].dim;

	fc->list = malloc(sizeof(floor1_list_t) * fc->values);
	fc->list[0].x = 0;
	fc->list[1].x = 1 << fc->rangebits;
	fc->list[0].sort = 0;
	fc->list[1].sort = 1;

	for (i = 2; i < fc->values; i++) {
		int j;
		GET_B(bp, fc->rangebits, fc->list[i].x);
		debug_msg("[%d] X %d\n", i, fc->list[i].x);
		fc->list[i].low = 0;
		fc->list[i].high = 1;
		fc->list[i].sort = i;
		for (j = 2; j < i; j++) {
			tmp = fc->list[j].x;
			if (tmp < fc->list[i].x && tmp > fc->list[fc->list[i].low].x)
				fc->list[i].low = j;
			if (tmp > fc->list[i].x && tmp < fc->list[fc->list[i].high].x)
				fc->list[i].high = j;
		}
	}
	for (i = 0; i < fc->values - 1; i++) {
		int j;
		for (j = i + 1; j < fc->values; j++) {
			if (fc->list[fc->list[i].sort].x > fc->list[fc->list[j].sort].x) {
				tmp = fc->list[i].sort;
				fc->list[i].sort = fc->list[j].sort;
				fc->list[j].sort = tmp;
			}
		}
	}
	return 0;
err_out:
	uninit_floor(fcc);
	return err;
}

static int render_point(int x0, int y0, int x1, int y1, int x) {
	return y0 +  (x - x0) * (y1 - y0) / (x1 - x0);
}

static void render_line(int x0, int y0, int x1, int y1, float * buf, int n) {
	int dy = y1 - y0;
	int adx = x1 - x0;
	int ady = MAX(dy, -dy);
	int base = dy / adx;
	int x = x0;
	int y = y0;
	int err = 0;
	int sy;
	if (dy < 0) sy = base - 1;
	else sy = base + 1;
	ady = ady - MAX(base, -base) * adx;
	if (x >= n) return;
	buf[x] = floor1_inverse[y];
	for (x = x0 + 1; x < x1; x++) {
		if (x >= n) return;
		err += ady;
		if (err >= adx) {
			err -= adx;
			y += sy;
		} else {
			y += base;
		}
		buf[x] = floor1_inverse[y];
	}
}

static int floor1_decode(vorbis_context_t * vc, floor_context_t * fcc, bit_packer_t * bp, float * buf, int n) {
	floor1_context_t * fc = &fcc->u.f1;
	int err = 0;
	int tmp;
	int i;
	int range = 255 / fc->multiplier + 1;
	int lx, ly;

	GET_B(bp, 1, tmp); // nozero
	if (!tmp) { // no decode
		memset(buf, 0, sizeof(float)*n);
		return -1;
	}

	GET_B(bp, ilog(range - 1), fc->list[0].y);
	GET_B(bp, ilog(range - 1), fc->list[1].y);

	tmp = 2;
	for (i = 0; i < fc->partitions; i++) {
		floor1_class_t * c = &fc->classes[fc->partition_to_class[i]];
		int cval = 0;
		int j;
		int csub = (1 << c->subclass) - 1;
		if (csub) GET_SCALAR(bp, c->masterbook, cval);
		for (j = 0; j < c->dim; j++) {
			int book = c->books[cval & csub];
			cval >>= c->subclass;
			if (book != -1) GET_SCALAR(bp, book, fc->list[tmp].y);
			else fc->list[tmp].y = 0;
			tmp++;
		}
	}
	fc->list[0].flag = 1;
	fc->list[1].flag = 1;
	for (i = 2; i < fc->values; i++) {
		int predicted = render_point(fc->list[fc->list[i].low].x,
					     fc->list[fc->list[i].low].y,
					     fc->list[fc->list[i].high].x,
					     fc->list[fc->list[i].high].y,
					     fc->list[i].x);
		int val = fc->list[i].y;
		int highroom = range -  predicted;
		int lowroom = predicted;
		int room = MIN(highroom, lowroom) * 2;
		if (val) {
			fc->list[fc->list[i].low].flag = 1;
			fc->list[fc->list[i].high].flag = 1;
			fc->list[i].flag = 1;
			if (val >= room) {
				if (highroom > lowroom) fc->list[i].y = val - lowroom + predicted;
				else fc->list[i].y = predicted - val + highroom - 1;
			} else {
				if (val & 1) val = -val;
				fc->list[i].y = predicted + (val >> 1);
			}
		} else {
			fc->list[i].flag = 0;
			fc->list[i].y = predicted;
		}
	}
	lx = 0;
	ly = fc->list[0].y  * fc->multiplier; // sorted 0 is still 0
	for (i = 1; i < fc->values; i++) {
		int pos = fc->list[i].sort;
		if (fc->list[pos].flag) {
			render_line(lx, ly, fc->list[pos].x, fc->list[pos].y * fc->multiplier, buf, n);
			lx = fc->list[pos].x;
			ly = fc->list[pos].y * fc->multiplier;
		}
		if (lx >= n) break;
	}
	if (lx < n) render_line(lx, ly, n, ly, buf, n);
err_out:
	return err;
}

static void uninit_residue(residue_context_t * rc) {
	if (rc->type == -1) return;

	free(rc->books);
	rc->books = NULL;
	rc->type = -1;
}

static int read_residue_header(vorbis_context_t * vc, residue_context_t * rc, bit_packer_t * bp) {
	int err;
	int i;
	int fast;

	rc->books = NULL;

	GET_B(bp, 16, rc->type);
	ERROR(rc->type > 2, 12);

	GET_B(bp, 24, rc->begin);
	GET_B(bp, 24, rc->end);
	GET_B(bp, 24, rc->partition_size); rc->partition_size++;
	GET_B(bp, 6, rc->classifications); rc->classifications++;
	GET_B(bp, 8, rc->classbook);
	ERROR(rc->classbook >= vc->ncodebooks, 13);

	debug_msg("begin %d end %d size %d classifs %d classbook %d\n", rc->begin, rc->end, rc->partition_size, rc->classifications, rc->classbook);

	{
	int cascade[rc->classifications];
	rc->books = malloc(sizeof(int[8]) * rc->classifications);

	for (i = 0; i < rc->classifications; i++) {
		int tmp;
		GET_B(bp, 3, cascade[i]);
		GET_B(bp, 1, tmp);
		if (tmp) {
			GET_B(bp, 5, tmp);
			cascade[i] |= tmp << 3;
		}
	}
	fast = rc->type == 2 && !(rc->begin & 1);
	for (i = 0; i < rc->classifications; i++) {
		int bit;
		debug_msg("[%d] { ", i);
		for (bit = 0; bit < 8; bit++) {
			if (cascade[i] & (1 << bit)) {
				GET_B(bp, 8, rc->books[i][bit]);
				ERROR(rc->books[i][bit] >= vc->ncodebooks, 13);
				fast &= !(vc->codebooks[rc->books[i][bit]].ndimentions & 1);
			} else {
				rc->books[i][bit] = -1;
			}
			debug_msg("%d, ", rc->books[i][bit]);
		}
		debug_msg("}\n");
	}
	}
	if (fast) rc->type = 3;
	return 0;
err_out:
	uninit_residue(rc);
	return err;
}

static int residue_decode(vorbis_context_t * vc, residue_context_t * rc, bit_packer_t * bp, float * residue_buf, int n, int real_ch, int * do_not_decode) {
	float * buf = residue_buf;
	int classwords = vc->codebooks[rc->classbook].ndimentions;
	int psize = rc->partition_size;
	int partitions = (rc->end - rc->begin) / psize;
	int err = 0;
	int i; // partition
	int j; // channel
	int k; // vector
	int l; // position in vector
	int dim, step; // dim of small book, partition/dim
	int pass;
	int ch = rc->type >= 2 ? 1 : real_ch;
	int classes[ch][partitions];
	int book;
	int p;
	float * tmp;
	int _num;
	float * buf1, * buf2;
	int * cs = classes[0];
	for (i = 0; i < real_ch; i++) for (p = 0; p < n; p++) buf[i * vc->blocksize1/2 + p] = 0.;
	if (rc->type >= 2) {
		for (i = 1; i < real_ch; i++) do_not_decode[0] &= do_not_decode[i];
		if (do_not_decode[0]) return 0;
	}
	if (rc->type == 3 && real_ch == 2) {
		for (pass = 0; pass < 8; pass++) {
			buf1 = buf                    + rc->begin;
			buf2 = buf + vc->blocksize1/2 + rc->begin;
			for (p = 0; p < partitions; ) {
				if (pass == 0) {
					GET_SCALAR(bp, rc->classbook, k);
					for (i = classwords; i--; ) {
						cs[p + i] = k % rc->classifications;
						k /= rc->classifications;
					}
				}
				for (i = 0; i < classwords && p < partitions; i++, p++) {
					book = rc->books[cs[p]][pass];
					if (book == -1) {
						buf1 += psize/2;
						buf2 += psize/2;
						continue;
					}
					dim = vc->codebooks[book].ndimentions;
					for (k = 0; k < psize; k += dim) {
						GET_VECTOR(bp, book, tmp);
						for (l = 0; l < dim; l += 2) {
							*buf1++ += *tmp++;
							*buf2++ += *tmp++;
						}
					}
				}
			}
		}
		return err;
	}
	for (pass = 0; pass < 8; pass++) {
		p = 0;
		while (p < partitions) {
			if (pass == 0) for (j = 0; j < ch; j++) {
				if (do_not_decode[j]) continue;
				GET_SCALAR(bp, rc->classbook, k);
				for (i = classwords; i--; ) {
					classes[j][p + i] = k % rc->classifications;
					k /= rc->classifications;
				}
			}
			for (i = 0; i < classwords && p < partitions; i++, p++) {
				for (j = 0; j < ch; j++) {
					if (do_not_decode[j]) continue;
					book = rc->books[classes[j][p]][pass];
					if (book == -1) continue;
					dim = vc->codebooks[book].ndimentions;
					step = psize / dim;
					buf2 = buf + j * vc->blocksize1/2 + rc->begin + p*psize;
					if (rc->type == 0) for (k = 0; k < step; k++) {
						GET_VECTOR(bp, book, tmp);
						for (l = k; l < psize; l += step)
							buf2[l] += tmp[l/step];
					} else if (rc->type == 1) for (k = 0; k < psize; k += dim) {
						GET_VECTOR(bp, book, tmp);
						for (l = 0; l < dim; l++)
							buf2[l+k] += tmp[l];
					} else {
						int s = rc->begin + p * psize;
						for (k = 0; k < psize; k += dim) {
							GET_VECTOR(bp, book, tmp);
							for (l = s + k; l < s + k + dim; l++)
								buf[(l % real_ch) * vc->blocksize1/2 + l / real_ch] += tmp[l-k-s];
						}
					}
				}
			}
		}
	}
err_out:
	return err;
}

static void priv_vorbis_uninit(vorbis_context_t * vc) {
	int i;
	if (vc->codebooks) for (i = 0; i < vc->ncodebooks; i++) uninit_codebook(&vc->codebooks[i]);
	if (vc->floors) for (i = 0; i < vc->nfloors; i++) uninit_floor(&vc->floors[i]);
	if (vc->residues) for (i = 0; i < vc->nresidues; i++) uninit_residue(&vc->residues[i]);
	if (vc->mappings) for (i = 0; i < vc->nmappings; i++) {
		free(vc->mappings[i].magnitude);
		free(vc->mappings[i].angle);
		free(vc->mappings[i].mux);
		free(vc->mappings[i].floor);
		free(vc->mappings[i].residue);
	}
	free(vc->codebooks);   vc->codebooks = NULL;
	free(vc->floors);      vc->floors = NULL;
	free(vc->residues);    vc->residues = NULL;
	free(vc->mappings);    vc->mappings = NULL;
	free(vc->modes);       vc->modes = NULL;
	free(vc->floor_buf);   vc->floor_buf = NULL;
	free(vc->residue_buf); vc->residue_buf = NULL;
	free(vc->saved_buf);   vc->saved_buf = NULL;
	free(vc->buf);         vc->buf = NULL;
	kiss_imdct_end(vc->imdct0); vc->imdct0 = NULL;
	kiss_imdct_end(vc->imdct1); vc->imdct1 = NULL;
}

void vorbis_uninit(vorbis_context_t * vc) {
	priv_vorbis_uninit(vc);
	free(vc);
}

int vorbis_read_headers(vorbis_context_t * vc, uint8_t * buf, int len) {
	bit_packer_t bp;
	int err;
	int comment_size;
	int i;

	ERROR(len < 2, 1);
	ERROR(buf[0] != 2, 2);
	ERROR(buf[1] != 30, 3);

	for (i = 2; buf[i] == 255; i++) if (i > len) return 1;
	comment_size = buf[i] + (i - 2)*255;
	i++;

	buf += i; len -= i;
	ERROR(len < 30 + comment_size, 1);

	init_bp(&bp, buf, 30);
	GET_B(&bp, 8, i); // magic
	ERROR(i != 1, 4);
	for (i = 0; "vorbis"[i]; i++) { int tmp; GET_B(&bp, 8, tmp); ERROR(tmp != "vorbis"[i], 4); }

	GET_B(&bp, 32, i); // version
	ERROR(i != 0, 4);
	GET_B(&bp, 8, vc->channels);
	ERROR(vc->channels == 0, 4);
	GET_B(&bp, 32, vc->sample_rate);
	ERROR(vc->sample_rate == 0, 4);
	GET_B(&bp, 32, i); // bitrate
	GET_B(&bp, 32, i); // bitrate
	GET_B(&bp, 32, i); // bitrate
	GET_B(&bp, 4, vc->blocksize0);
	GET_B(&bp, 4, vc->blocksize1);
	ERROR(vc->blocksize0 > vc->blocksize1, 4);
	ERROR(vc->blocksize0 > 13 || vc->blocksize0 < 6, 4);
	ERROR(vc->blocksize1 > 13 || vc->blocksize1 < 6, 4);
	{
		const float * vwin[8] = { vwin64, vwin128, vwin256, vwin512, vwin1024, vwin2048,
					vwin4096, vwin8192 };
		vc->win0 = vwin[vc->blocksize0 - 6];
		vc->win1 = vwin[vc->blocksize1 - 6];
	}

	debug_msg("==== blocksize 0: %d 1: %d\n", vc->blocksize0, vc->blocksize1);

	vc->blocksize0 = 1 << vc->blocksize0;
	vc->blocksize1 = 1 << vc->blocksize1;
	GET_B(&bp, 1, i); // framing
	ERROR(!i, 4);

	init_bp(&bp, buf + 30 + comment_size, len - 30 - comment_size);

	GET_B(&bp, 8, i);
	ERROR(i != 5, 4);
	for (i = 0; "vorbis"[i]; i++) { int tmp; GET_B(&bp, 8, tmp); ERROR(tmp != "vorbis"[i], 4); }

	// codebook
	GET_B(&bp, 8, vc->ncodebooks);
	vc->ncodebooks++;
	vc->codebooks = malloc(sizeof(codebook_t) * vc->ncodebooks);
	for (i = 0; i < vc->ncodebooks; i++) vc->codebooks[i].nentries = -1;
	for (i = 0; i < vc->ncodebooks; i++) {
		CHECK(read_codebook_header(vc, &vc->codebooks[i], &bp));
	}

	// time domain
	GET_B(&bp, 6, i);
	for (i++; i > 0; i--) {
		int tmp;
		GET_B(&bp, 16, tmp);
		ERROR(tmp != 0, 8);
	}

	// floors
	GET_B(&bp, 6, vc->nfloors);
	vc->nfloors++;
	vc->floors = malloc(sizeof(floor_context_t) * vc->nfloors);
	for (i = 0; i < vc->nfloors; i++) vc->floors[i].type = -1;
	for (i = 0; i < vc->nfloors; i++) {
		GET_B(&bp, 16, vc->floors[i].type);
		ERROR(vc->floors[i].type > 1, 10);
		if (vc->floors[i].type == 0) CHECK(read_floor0_header(vc, &vc->floors[i], &bp));
		else CHECK(read_floor1_header(vc, &vc->floors[i], &bp));
	}

	// residues
	GET_B(&bp, 6, vc->nresidues);
	vc->nresidues++;
	vc->residues = malloc(sizeof(residue_context_t) * vc->nresidues);
	for (i = 0; i < vc->nresidues; i++) vc->residues[i].type = -1;
	for (i = 0; i < vc->nresidues; i++) {
		CHECK(read_residue_header(vc, &vc->residues[i], &bp));
	}

	// mappings
	GET_B(&bp, 6, vc->nmappings);
	vc->nmappings++;
	vc->mappings = malloc(sizeof(mapping_context_t) * vc->nmappings);
	for (i = 0; i < vc->nmappings; i++) {
		vc->mappings[i].magnitude = NULL;
		vc->mappings[i].angle = NULL;
		vc->mappings[i].mux = NULL;
		vc->mappings[i].floor = NULL;
		vc->mappings[i].residue = NULL;
	}
	for (i = 0; i < vc->nmappings; i++) {
		mapping_context_t * mc = &vc->mappings[i];
		int tmp;
		int j;

		GET_B(&bp, 16, tmp);
		ERROR(tmp, 13); // bad mapping type

		GET_B(&bp, 1, mc->submaps);
		if (mc->submaps) GET_B(&bp, 4, mc->submaps);
		mc->submaps++;

		mc->coupling_steps = 0;
		GET_B(&bp, 1, tmp);
		if (tmp) { // square polar
			GET_B(&bp, 8, mc->coupling_steps);
			mc->coupling_steps++;
			mc->magnitude = malloc(sizeof(int) * mc->coupling_steps);
			mc->angle = malloc(sizeof(int) * mc->coupling_steps);
			for (j = 0; j < mc->coupling_steps; j++) {
				GET_B(&bp, ilog(vc->channels - 1), mc->magnitude[j]);
				GET_B(&bp, ilog(vc->channels - 1), mc->angle[j]);
				ERROR(mc->magnitude[j] >= vc->channels, 14);
				ERROR(mc->angle[j] >= vc->channels, 14);
			}
		}
		GET_B(&bp, 2, tmp); // reserved
		ERROR(tmp, 15);
		mc->mux = malloc(sizeof(int) * vc->channels);
		for (j = 0; j < vc->channels; j++) {
			mc->mux[j] = 0;
			if (mc->submaps > 1) GET_B(&bp, 4, mc->mux[j]);
			ERROR(mc->mux[j] >= mc->submaps, 15);
		}
		mc->floor = malloc(sizeof(int) * mc->submaps);
		mc->residue = malloc(sizeof(int) * mc->submaps);
		for (j = 0; j < mc->submaps; j++) {
			GET_B(&bp, 8, tmp); // ignored time configuration
			GET_B(&bp, 8, mc->floor[j]);
			ERROR(mc->floor[j] >= vc->nfloors, 16);
			GET_B(&bp, 8, mc->residue[j]);
			ERROR(mc->residue[j] >= vc->nresidues, 17);
		}
	}

	// modes
	GET_B(&bp, 6, vc->nmodes);
	vc->nmodes++;
	vc->modes = malloc(sizeof(mode_context_t) * vc->nmodes);
	for (i = 0; i < vc->nmodes; i++) {
		int tmp;
		GET_B(&bp, 1, vc->modes[i].blockflag);
		GET_B(&bp, 16, tmp); // window type
		ERROR(tmp, 18);
		GET_B(&bp, 16, tmp); // transform type
		ERROR(tmp, 18);
		GET_B(&bp, 8, vc->modes[i].mapping);
	}
	GET_B(&bp, 1, i); // framing
	ERROR(!i, 19);

	vc->floor_buf = malloc(sizeof(float[vc->blocksize1/2]) * vc->channels);
	vc->residue_buf = malloc(sizeof(float[vc->blocksize1/2]) * vc->channels);
	vc->saved_buf = malloc(sizeof(float[vc->blocksize1/2]) * vc->channels);
	vc->buf = malloc(sizeof(float[vc->blocksize1]) * vc->channels);

	vc->imdct0 = kiss_imdct_init(vc->blocksize0);
	vc->imdct1 = kiss_imdct_init(vc->blocksize1);

	return 0;
err_out:
	priv_vorbis_uninit(vc);
	return err;
}

vorbis_context_t * vorbis_init() {
	vorbis_context_t * vc = malloc(sizeof(vorbis_context_t));
	vc->codebooks = NULL;
	vc->floors = NULL;
	vc->residues = NULL;
	vc->mappings = NULL;
	vc->modes = NULL;
	vc->floor_buf = NULL;
	vc->residue_buf = NULL;
	vc->saved_buf = NULL;
	vc->buf = NULL;
	vc->imdct0 = NULL;
	vc->imdct1 = NULL;
	vc->saved = 0;
	return vc;
}

static int window(vorbis_context_t * vc, int bs, int prevblock, int nextblock, int16_t(*out)[2]) {
	float * floor_buf = vc->floor_buf;
	float * saved_buf = vc->saved_buf;
	float * buf_tmp = vc->buf;
	float * floor;
	float * saved;
	float * buf;
	int i, j;
	int start, end;
	int begin_win = prevblock ? 0 : (bs - vc->blocksize0) / 4;
	int begin_middle_win = prevblock ? bs/2 : (bs + vc->blocksize0) / 4;
	int end_middle_win = bs / 2 + (nextblock ? 0 : (bs - vc->blocksize0) / 4);
	int end_win = bs / 2 + (nextblock ? bs/2 : (bs + vc->blocksize0) / 4);
	const float * win;
	SLOW_DEBUG(static int top = 0; static int dcount = 0;)

	// imdct
	for (i = 0; i < vc->channels; i++) {
		floor = floor_buf + i * vc->blocksize1/2;
		saved = saved_buf + i * vc->blocksize1/2;
		buf = buf_tmp + i * vc->blocksize1;

		kiss_imdct(bs == vc->blocksize1 ? vc->imdct1 : vc->imdct0, floor, buf);

#define BIAS 385
		// overlap_add with prev
		if (vc->saved) {
			start = begin_win;
			end = begin_middle_win;
			win = prevblock ? vc->win1 : vc->win0;
			saved -= start;
			win -= start;
			for (j = start; j < end; j++) {
				buf[j] = saved[j] + buf[j] * win[j] + BIAS;
			}
			saved += start;
			win += start;
		}
		start = begin_middle_win;
		end = end_middle_win;
		for (j = start; j < end; j++) buf[j] += BIAS;

		// prepare save for next
		start = end_middle_win;
		end = end_win;
		win = nextblock ? vc->win1 : vc->win0;
		saved -= start;
		for (j = start; j < end; j++) {
			saved[j] = buf[j] * win[end - j - 1];
		}
		saved += start;

		if (vc->saved) start = begin_win;
		else start = begin_middle_win;
		end = end_middle_win;

		out -= start;
		for (j = start; j < end; j++) {
			int tmp = ((int32_t*)buf)[j];
			if (tmp & 0xf0000) {
				if (tmp > 0x43c0ffff) tmp = 0xFFFF;
				else                  tmp = 0;
			}
			out[j][i] = tmp - 0x8000;
			//out[j][i] = buf[j]*32768;
			SLOW_DEBUG(top = MAX(top,ABS(out[j][i]));)
		}
		out += start;
	}
	SLOW_DEBUG(
		if ((dcount += bs) > DEBUG_RATE) {
			if (top > 30000) debug_msg("\033[01;40;31msample: %4d\033[00;m\n", top);
			else debug_msg("sample: %4d\033[00;m\n", top);
			dcount -= DEBUG_RATE;
			top = 0;
		}
	)

	if (vc->saved) start = begin_win;
	else start = begin_middle_win;
	end = end_middle_win;

	if (i == 1) {
		// single channel
		for (j = 0; j < end - start; j++) out[j][1] = 0;
	}

	vc->saved = 1;

	return end - start;
}

int vorbis_decode(vorbis_context_t * vc, uint8_t * in, int len, int16_t (*out)[2]) {
	bit_packer_t bp;
	mode_context_t * mode;
	mapping_context_t * mapping;
	float * residue_buf = vc->residue_buf;
	float * floor_buf = vc->floor_buf;
	int no_residue[vc->channels];
	int chan_to_res[vc->channels];
	int err = 0;
	int i;
	int bs;
	int tmp;
	int nextblock = 0, prevblock = 0;
	SLOW_DEBUG(static float top1 = 0, top2 = 0, top3 = 0; static int dcount = 0;)

	init_bp(&bp, in, len);

	GET_B(&bp, 1, tmp);
	ERROR(tmp, 2);
	GET_B(&bp, ilog(vc->nmodes - 1), tmp);
	ERROR(tmp >= vc->nmodes, 3);

	mode = &vc->modes[tmp];
	mapping = &vc->mappings[mode->mapping];

	if (mode->blockflag) {
		GET_B(&bp, 1, prevblock);
		GET_B(&bp, 1, nextblock);
	}

	bs = mode->blockflag ? vc->blocksize1 : vc->blocksize0;

 	for (i = 0; i < vc->channels; i++) {
		floor_context_t * floor = &vc->floors[mapping->floor[mapping->mux[i]]];
		no_residue[i] = 0;
		if (floor->type == 0) err = floor0_decode(vc, floor, &bp, floor_buf + i * vc->blocksize1/2, mode->blockflag);
		else err = floor1_decode(vc, floor, &bp, floor_buf + i * vc->blocksize1/2, bs / 2);
		if (err == -1) { err = 0; no_residue[i] = 1; }
		if (err) goto err_out;
	}

	tmp = 0;
	for (i = 0; i < mapping->submaps; i++) {
		residue_context_t * residue = &vc->residues[mapping->residue[i]];
		int ch = 0;
		int j;
		int do_not_decode[vc->channels];
		for (j = 0; j < vc->channels; j++) {
			if (mapping->mux[j] == i) {
				do_not_decode[ch] = no_residue[j];
				chan_to_res[j] = tmp + ch;
				ch++;
			}
		}
		if ((err = residue_decode(vc, residue, &bp, residue_buf + tmp * vc->blocksize1/2, bs / 2, ch, do_not_decode))) {
			if (err == 1) err = 0;
			else goto err_out;
		}
		tmp += ch;
	}

	// channel coupling
	for (i = mapping->coupling_steps - 1; i >= 0; i--) {
		float * mag = residue_buf + chan_to_res[mapping->magnitude[i]] * vc->blocksize1/2;
		float * ang = residue_buf + chan_to_res[mapping->angle[i]] * vc->blocksize1/2;
		int j;
		float m, a;
		if (no_residue[mapping->angle[i]]) continue;
		for (j = 0; j < bs/2; j++) {
			m = mag[j];
			a = ang[j];

			if (a <= 0) {
				ang[j] = m;
				if (m <= 0) mag[j] = m - a;
				else mag[j] = m + a;
			} else {
				if (m > 0) ang[j] = m - a;
				else ang[j] = m + a;
			}
		}
	}
	// dot product
	for (i = 0; i < vc->channels; i++) {
		float * floor = floor_buf + i * vc->blocksize1/2;
		float * residue = residue_buf + chan_to_res[i] * vc->blocksize1/2;
		int j;
		if (no_residue[i]) continue;
		for (j = 0; j < bs/2; j++) {
			SLOW_DEBUG(top2 = MAX(top2,ABS(floor[j]));)
			SLOW_DEBUG(top3 = MAX(top3,ABS(residue[j]));)
			floor[j] *= residue[j];
			SLOW_DEBUG(top1 = MAX(top1,ABS(floor[j]));)
		}
	}
	SLOW_DEBUG(
		if ((dcount += bs) > DEBUG_RATE) {
			debug_msg("\033[01;40;37mmax floor: %f residue: %6.2f coeff: %f ", top2, top3, top1);
			dcount -= DEBUG_RATE;
			top2 = top3 = top1 = 0;
		}
	)

	tmp = window(vc, bs, prevblock, nextblock, out);

	fwrite(out, 1, tmp*2*2, stdout);
err_out:
	return err;
}
