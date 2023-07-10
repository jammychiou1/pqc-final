#include <stdio.h>

#include "rader_outer/polymul.h"

int16_t poly_a[1632];
int16_t poly_b[1632];
int16_t poly_c[1632];

int main() {
  for (int i = 0; i < 761; i++) {
    scanf("%hd\n", &poly_a[i]);
  }
  for (int i = 0; i < 761; i++) {
    scanf("%hd\n", &poly_b[i]);
  }
  polymul(poly_c, poly_a, poly_b);
  for (int i = 0; i < 761; i++) {
    printf("%hd%c", poly_c[i], " \n"[i == 760]);
  }
}
