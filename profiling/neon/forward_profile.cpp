#include <cstdlib>

#include "neon/ntt_10.h"
#include "neon/ntt_9.h"

int16_t ntt[10][9][16];
int16_t a[800];
int main() {
  for (int i = 0; i < 761; i++) {
    a[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    a[0] = rand() % 4591 - 2295;
    ntt_9(ntt, a);
    ntt_10(ntt);
  }
}
