#include <cstdlib>

#include "neon/ntt_9.h"

int16_t ntt[9][2][10][8];
int16_t a[800];
int main() {
  for (int t = 0; t < 1000000; t++) {
    a[0] = rand() % 4591 - 2295;
    ntt_9(ntt, a);
  }
}
