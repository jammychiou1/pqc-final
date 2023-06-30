#include <cstring>

#include "utils/pretty.h"

void mult(int16_t in1_poly[], int16_t in2_poly[], int16_t out_poly[]);

int main() {
  int16_t a[768], b[768], c[768];
  scan_poly(a);
  scan_poly(b);
  std::memset(&a[761], 0, sizeof(int16_t) * 7);
  std::memset(&b[761], 0, sizeof(int16_t) * 7);
  mult(a, b, c);
  print_poly(c);
}
