#include <cstdlib>

#include "neon/base_mul.h"

int16_t ntt_a[9][2][10][8];
int16_t ntt_b[9][2][10][8];
int16_t ntt_c[9][2][10][8];

int main() {
  for (int t = 0; t < 1000000; t++) {
    ntt_a[0][0][0][0] = rand() % 4591 - 2295;
    ntt_b[0][0][0][0] = rand() % 4591 - 2295;
    base_mul(ntt_a, ntt_b, ntt_c);
  }
}
