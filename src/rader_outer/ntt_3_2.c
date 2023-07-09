#include "rader_outer/ntt_3_2.h"

#define Q 4591
#define Qprime (-15631)
#define OMEGA3lo (-311)
#define OMEGA3hi (-2220)
#define Qbar (29235)

static int16_t __attribute__((aligned (16))) radix3_args[8] = {Q, Qprime, OMEGA3lo, OMEGA3hi, Qbar};

void __asm_radix32(int16_t*, int16_t*);

void ntt_3_2(int16_t ntt[1632]) {
  __asm_radix32(ntt, radix3_args);
}
