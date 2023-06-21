#include <iostream>

#include "mult.h"

int main() {
  int16_t a[P] = {}, b[P] = {}, c[P];
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
