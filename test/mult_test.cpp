#include <iostream>
#include <cstring>

#include "mult.h"
#include "pretty.h"

int main() {
  int16_t a[768], b[768], c[P];
  scan_poly(a);
  scan_poly(b);
  std::memset(&a[761], 0, sizeof(int16_t) * 7);
  std::memset(&b[761], 0, sizeof(int16_t) * 7);
  mult(a, b, c);
  print_poly(c);
}
