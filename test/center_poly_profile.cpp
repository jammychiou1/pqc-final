#include <cstdlib>

void center_poly(int16_t poly[]);

int16_t out_poly[768];
int main() {
  for (int i = 0; i < 768; i++) {
    out_poly[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    out_poly[0] = rand() % 4591 - 2295;
    center_poly(out_poly);
  }
}
