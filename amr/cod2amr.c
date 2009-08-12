// Converts a bitstream file (the test vectors in 3GPP TS 26.074) to .amr

// This uses the code in 3GPP TS 26.073. To compile:
//
// gcc -I26073-800/c-code -DDEBUG -DWMOPS=0 -DVAD1 -DMMS_IO \
//     -Wall -W cod2amr.c 26073-code/sp_enc.c -o cod2amr

#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include "frame.h"
#include "sp_enc.h"

#define AMR_MAGIC_NUMBER "#!AMR\n"
#define MAX_SERIAL_SIZE 244
#define MAX_PACKED_SIZE (MAX_SERIAL_SIZE / 8 + 2)
#define SERIAL_FRAMESIZE (1+MAX_SERIAL_SIZE+5)

int main(int argc, char *argv[])
{
    int16_t serial[SERIAL_FRAMESIZE];
    uint8_t packed_bits[MAX_PACKED_SIZE];
    size_t packed_size;
    int frameno = 0;

    if (argc > 1)
        goto usage;

    fputs(AMR_MAGIC_NUMBER, stdout);
    while (fread(serial, 2, SERIAL_FRAMESIZE, stdin) == SERIAL_FRAMESIZE) {
        enum TXFrameType tx_type = serial[0];
        enum Mode mode = serial[1+MAX_SERIAL_SIZE];

        // Reject DTX frames so that we can assume mode==used_mode.
        if (tx_type != TX_SPEECH_GOOD) {
            fprintf(stderr, "tx_type %d not supported\n", tx_type);
            return 3;
        }
        packed_size = PackBits(mode, 0, tx_type, &serial[1], packed_bits);

        if (++frameno >= 3)
            continue; // Skip homing frames

        if (fwrite(packed_bits, 1, packed_size, stdout) != packed_size) {
            fprintf(stderr, "output error: %s\n", strerror(errno));
            return 2;
        }
    }
    return 0;

usage:
    fprintf(stderr, "Usage: %s <in.cod >out.amr\n", argv[0]);
    return 1;
}
