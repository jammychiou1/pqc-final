#include <cstdlib>

void base_mul(int16_t in1_ntt[10][9][16], int16_t in2_ntt[10][9][16], int16_t out_ntt[10][9][16]);

int16_t ntt_a[10][9][16];
int16_t ntt_b[10][9][16];
int16_t ntt_c[10][9][16];

int main() {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      for (int k = 0; k < 16; k++) {
        ntt_a[i][j][k] = rand() % 4591 - 2295;
        ntt_b[i][j][k] = rand() % 4591 - 2295;
      }
    }
  }
  for (int t = 0; t < 1000000; t++) {
    ntt_a[0][0][0] = rand() % 4591 - 2295;
    ntt_b[0][0][0] = rand() % 4591 - 2295;
    base_mul(ntt_a, ntt_b, ntt_c);
  }
}
