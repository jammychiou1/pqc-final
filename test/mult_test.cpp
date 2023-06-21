#include <iostream>

#include "mult.h"
#include "pretty.h"

int main() {
  int16_t a[P], b[P], c[P];
  scan_poly(a);
  scan_poly(b);
  mult(a, b, c);
  print_poly(c);
}
