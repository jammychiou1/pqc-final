#include <iostream>

#include "utils/pretty.h"

void mult(int16_t in1_poly[], int16_t in2_poly[], int16_t out_poly[]);

int main() {
  int16_t a[768] = {}, b[768] = {}, c[768];
  int n;
  std::cin >> n;
  for (int i = 0; i <= n; i++) {
    std::cin >> a[i];
  }
  for (int i = 0; i <= n; i++) {
    std::cin >> b[i];
  }
  mult(a, b, c);
  for (int i = 0; i <= 2 * n; i++) {
    std::cout << c[i] << " \n"[i == 2 * n];
  }
}
