#include <cstdlib>

#include "neon/mult_low.h"

int16_t low_a[160];
int16_t low_b[160];
int16_t low_c[160];

int main() {
  for (int i = 0; i < 160; i++) {
    low_a[i] = rand() % 4591 - 2295;
    low_b[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    // low_a[0] = rand() % 4591 - 2295;
    // low_b[0] = rand() % 4591 - 2295;
    mult_low(low_a, low_b, low_c);
  }
}
