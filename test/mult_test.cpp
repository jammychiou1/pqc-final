#include "neon/neon.h"
#include "utils/pretty.h"

int16_t a[800];
int16_t b[800];
int16_t c[768];

int main() {
  scan_poly(a);
  scan_poly(b);
  mult(a, b, c);
  print_poly(c);
}
