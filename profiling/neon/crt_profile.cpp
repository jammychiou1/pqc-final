#include <cstdlib>

void crt(int16_t poly[], int16_t main_poly[], int16_t low[]);

int16_t out_main[1448];
int16_t out_low[81];
int16_t out_poly[761];
int main() {
  for (int i = 0; i < 1440; i++) {
    out_main[i] = rand() % 4591 - 2295;
  }
  for (int i = 0; i < 81; i++) {
    out_low[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    out_main[0] = rand() % 4591 - 2295;
    out_low[0] = rand() % 4591 - 2295;
    crt(out_poly, out_poly, out_low);
  }
}
