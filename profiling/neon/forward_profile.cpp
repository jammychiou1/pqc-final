#include <cstdlib>

#include "neon/ntt_10.h"
#include "neon/ntt_9.h"

int16_t ntt[10][9][16];
int main() {
  int16_t a[768];
  for (int i = 0; i < 768; i++) {
    if (i < 761) {
      a[i] = rand() % 4591 - 2295;
    }
    else {
      a[i] = 0;
    }
  }
  for (int t = 0; t < 1000000; t++) {
    a[0] = rand() % 4591 - 2295;
    ntt_10(ntt, a);
    ntt_9(ntt);
  }
}
