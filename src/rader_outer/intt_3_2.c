#include "rader_outer/intt_3_2.h"

#define Q 4591
#define Qprime (-15631)
#define OMEGA3SQlo (310)
#define OMEGA3SQhi (2213)
#define Qbar (29235)

static int16_t __attribute__((aligned (16))) iradix3_args[8] = {Q, Qprime, OMEGA3SQlo, OMEGA3SQhi, Qbar};

void __asm_radix32(int16_t*, int16_t*);

void intt_3_2(int16_t ntt[1632]) {
  __asm_radix32(ntt, iradix3_args);
}
