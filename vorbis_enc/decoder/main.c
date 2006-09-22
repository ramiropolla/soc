#include <sys/ioctl.h>
#include <linux/soundcard.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>
#include <nut.h>

typedef struct vorbis_context_s vorbis_context_t;
int vorbis_read_headers(vorbis_context_t * vc, uint8_t * buf, int len);
int vorbis_decode(vorbis_context_t * vc, uint8_t * in, int len, uint16_t * out);
vorbis_context_t * vorbis_init();
void vorbis_uninit(vorbis_context_t * vc);

int main(int argc, char * argv []) {
	nut_packet_t pd;
	FILE * file = argv[1] ? fopen(argv[1], "r") : stdin;
	nut_demuxer_opts_t dopts = { { file, NULL, NULL, NULL }, 1 };
	nut_context_t * nut = nut_demuxer_init(&dopts);
	vorbis_context_t * vc = vorbis_init();
	nut_stream_header_t * s = NULL;
	int err = 0;
	int i = 0;

	if ((err = nut_read_headers(nut, &s))) return err;

	{
	uint8_t buf[s->codec_specific_len+7];
	memcpy(buf, s->codec_specific, s->codec_specific_len);
	memset(buf+s->codec_specific_len, 0, 7);
	if ((err = vorbis_read_headers(vc, buf, s->codec_specific_len))) return err;
	}

	free(s);

	i = AFMT_S16_LE; ioctl(1, SNDCTL_DSP_SETFMT, &i);
	i = 2;           ioctl(1, SNDCTL_DSP_CHANNELS, &i);
	i = 44100;       ioctl(1, SNDCTL_DSP_SPEED, &i);

	fcntl(0, F_SETFL, O_NONBLOCK);

	if (argc > 2) nut_seek(nut, strtod(argv[2], NULL), 0, NULL);

	while (!nut_read_next_packet(nut, &pd)) {
		int len = pd.len;
		uint8_t buf[len+7];
		uint16_t out[2048*2];
		nut_read_frame(nut, &pd.len, buf);
		memset(buf+len, 0, 7);

		if ((err = vorbis_decode(vc, buf, len, out))) return err;
		fflush(stdout);
		while (read(0, buf, 1) > 0) {
			double n = 0;
			switch (*buf) {
				case 'a': n = -10; break;
				case 'd': n =  10; break;
				case 'w': n =  60; break;
				case 's': n = -60; break;
			}
			if (n) nut_seek(nut, n, 1 | ((n>0) ? 2: 0), NULL);
		}
	}
	nut_demuxer_uninit(nut);
	vorbis_uninit(vc);
	fclose(file);
	return err;
}
