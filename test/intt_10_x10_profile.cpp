#include <cstdlib>

void intt_10_x10(int16_t ntt[10][9][16], int16_t poly[1440]);

int16_t ntt[10][9][16];
int16_t poly[1440];
int main() {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      for (int k = 0; k < 16; k++) {
        ntt[i][j][k] = rand() % 4591 - 2295;
      }
    }
  }
  for (int t = 0; t < 1000000; t++) {
    ntt[0][0][0] = rand() % 4591 - 2295;
    intt_10_x10(ntt, poly);
  }
}
