#include <cstdlib>

void ntt_10(int16_t ntt[10][9][16], int16_t poly[1440]);

int16_t ntt[10][9][16];
int main() {
  int16_t a[768];
  for (int i = 0; i < 768; i++) {
    if (i < 761) {
      a[i] = rand() % 4591 - 2295;
    }
    else {
      a[i] = 0;
    }
  }
  for (int t = 0; t < 1000000; t++) {
    a[0] = rand() % 4591 - 2295;
    ntt_10(ntt, a);
  }
}
