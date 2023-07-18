#include <cstdlib>

#include "neon/ntt_10.h"

int16_t ntt[9][2][10][8];
int main() {
  for (int t = 0; t < 1000000; t++) {
    ntt[0][0][0][0] = rand() % 4591 - 2295;
    ntt_10(ntt);
  }
}
