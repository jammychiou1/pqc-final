#include <cstring>
#include <cstdlib>

#include "mult.h"

int main() {
  int16_t a[768], b[768], c[768];
  std::memset(&a[761], 0, sizeof(int16_t) * 7);
  std::memset(&b[761], 0, sizeof(int16_t) * 7);
  for (int i = 0; i < 761; i++) {
    a[i] = rand() % 4591 - 2295;
    b[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 100000; t++) {
    a[0] = rand() % 4591 - 2295;
    b[0] = rand() % 4591 - 2295;
    mult(a, b, c);
  }
}
