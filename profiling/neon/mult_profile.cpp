#include <cstring>
#include <cstdlib>

#include "neon/neon.h"

int16_t a[800];
int16_t b[800];
int16_t c[768];

int main() {
  for (int i = 0; i < 761; i++) {
    a[i] = rand() % 4591 - 2295;
    b[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    a[0] = rand() % 4591 - 2295;
    b[0] = rand() % 4591 - 2295;
    mult(a, b, c);
  }
}
