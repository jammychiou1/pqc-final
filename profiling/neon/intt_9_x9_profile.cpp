#include <cstdlib>

#include "neon/intt_9_x9.h"

int16_t ntt[9][2][10][8];
int16_t poly[1400];
int main() {
  for (int t = 0; t < 1000000; t++) {
    ntt[0][0][0][0] = rand() % 4591 - 2295;
    intt_9_x9(ntt, poly);
  }
}
