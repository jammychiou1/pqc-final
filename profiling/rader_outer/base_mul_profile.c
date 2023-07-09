#include <stdlib.h>

#include "rader_outer/base_mul.h"

int16_t ntt_a[1632];
int16_t ntt_b[1632];
int16_t ntt_c[1632];

int main() {
  for (int i = 0; i < 1632; i++) {
    ntt_a[i] = rand() % 4591 - 2295;
    ntt_b[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    ntt_a[0] = rand() % 4591 - 2295;
    ntt_b[0] = rand() % 4591 - 2295;
    base_mul(ntt_a, ntt_b, ntt_c);
  }
}
