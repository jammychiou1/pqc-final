#include <iostream>

#include "neon/neon.h"

int16_t a[800];
int16_t b[800];
int16_t c[768];

int main() {
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
