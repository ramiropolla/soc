SOURCES = kiss_imdct.c vorbis.c main.c libnut/demuxer.c
CFLAGS  = -lm -Ilibnut -Wall -D_LARGEFILE_SOURCE -ffast-math
CFLAGS += -g

CFLAGS += -fomit-frame-pointer
#CFLAGS += -pg

CFLAGS += -O4
#CFLAGS += -Os
CFLAGS += -finline-functions -finline-limit=1000 -frename-registers

CFLAGS += -march=athlon-xp
#CFLAGS += -march=i686

CC = gcc
#CC = /usr/bin/gcc-3.3

nutvorbis-ffmpeg-mmx: ${SOURCES}
	${CC} -DUSE_FFMPEG -DHAVE_MMX -msse ${CFLAGS} $^ -o $@
nutvorbis-kiss: ${SOURCES}
	${CC} ${CFLAGS} $^ -o $@
nutvorbis-fftw: ${SOURCES}
	${CC} -DUSE_FFTW -lfftw3f ${CFLAGS} $^ -o $@
nutvorbis-ffmpeg: ${SOURCES}
	${CC} -DUSE_FFMPEG ${CFLAGS} $^ -o $@

all: nutvorbis-kiss nutvorbis-fftw nutvorbis-ffmpeg nutvorbis-ffmpeg-mmx
