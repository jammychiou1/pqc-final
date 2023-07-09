#include <stdlib.h>

#include "rader_outer/polymul.h"

int16_t poly_a[1632];
int16_t poly_b[1632];
int16_t poly_c[1632];

int main() {
  for (int i = 0; i < 761; i++) {
    poly_a[i] = rand() % 4591 - 2295;
    poly_b[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    poly_a[0] = rand() % 4591 - 2295;
    poly_b[0] = rand() % 4591 - 2295;
    polymul(poly_a, poly_b, poly_c);
  }
}
