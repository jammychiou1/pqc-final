// #include <arm_neon.h>
#include <iostream>
#include <cstring>
#include <cassert>
#include <array>

#include "mult.h"
#include "pretty.h"

template <size_t SZ>
constexpr std::array<int32_t, SZ> gen_pows(int32_t base) {
  std::array<int32_t, SZ> pows = {1};
  for (int i = 1; i < SZ; i++) {
    pows[i] = center_lift(pows[i - 1] * base);
  }
  return pows;
}

constexpr int ORD = 4590;
constexpr int32_t W_4590 = 11;

constexpr std::array<int32_t, ORD> ws = gen_pows<ORD>(W_4590);

constexpr int32_t W_10 = ws[ORD / 10];
constexpr int32_t W_5 = ws[ORD / 5];
constexpr int32_t W_9 = ws[ORD / 9];
constexpr int32_t W_3 = ws[ORD / 3];

constexpr std::array<int32_t, 5> w_5s = gen_pows<5>(W_5);
constexpr std::array<int32_t, 10> w_10s = gen_pows<10>(W_10);
constexpr std::array<int32_t, 9> w_9s = gen_pows<9>(W_9);
constexpr std::array<int32_t, 3> w_3s = gen_pows<3>(W_3);

void ntt_10(int16_t ntt[10][9][16]) {

  for (int j = 0; j < 9; j++) {
    for (int k0 = 0; k0 < 16; k0 += 8) {

      int16_t tmp1[5][2][8] = {};
      int16_t tmp2[5][2][8] = {};

      for (int i = j + 5; i < j + 11; i++) {
        for (int dk = 0; dk < 8; dk++) {
          tmp1[i % 5][i % 2][dk] = ntt[i % 10][j][k0 + dk];
        }
      }

      // TODO skip zero
      for (int i2 = 0; i2 < 2; i2++) {
        for (int i1 = 0; i1 < 5; i1++) {
          for (int ii1 = 0; ii1 < 5; ii1++) {
            for (int dk = 0; dk < 8; dk++) {
              tmp2[i1][i2][dk] += center_lift(w_5s[i1 * ii1 % 5] * tmp1[ii1][i2][dk]);
              assert(-3 * CENTER_MAG <= tmp2[i1][i2][dk] && tmp2[i1][i2][dk] <= 3 * CENTER_MAG);
            }
            // |tmp2| ~ 1.5Q
          }
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        for (int dk = 0; dk < 8; dk++) {
          tmp1[i1][0][dk] = tmp2[i1][0][dk] + tmp2[i1][1][dk];
          tmp1[i1][1][dk] = tmp2[i1][0][dk] - tmp2[i1][1][dk];
          assert(-6 * CENTER_MAG <= tmp1[i1][0][dk] && tmp1[i1][0][dk] <= 6 * CENTER_MAG);
          assert(-6 * CENTER_MAG <= tmp1[i1][1][dk] && tmp1[i1][1][dk] <= 6 * CENTER_MAG);
        }
        // |tmp1| ~ 3Q
      }

      for (int i = 0; i < 10; i++) {
        for (int dk = 0; dk < 8; dk++) {
          ntt[7 * i % 10][j][k0 + dk] = tmp1[i % 5][i % 2][dk];
        }
        // |ntt| ~ 3Q
      }

    }
  }

}

void intt_10_x10(int16_t ntt[10][9][16]) {

  for (int j = 0; j < 9; j++) {
    for (int k0 = 0; k0 < 16; k0 += 8) {

      int16_t tmp1[5][2][8] = {};
      int16_t tmp2[5][2][8] = {};

      for (int i = 0; i < 10; i++) {
        for (int dk = 0; dk < 8; dk++) {
          tmp1[i % 5][i % 2][dk] = ntt[(10 - i) % 10][j][k0 + dk];
        }
      }

      for (int i2 = 0; i2 < 2; i2++) {
        for (int i1 = 0; i1 < 5; i1++) {
          for (int ii1 = 0; ii1 < 5; ii1++) {
            for (int dk = 0; dk < 8; dk++) {
              tmp2[i1][i2][dk] += center_lift(w_5s[i1 * ii1 % 5] * tmp1[ii1][i2][dk]);
              assert(-5 * CENTER_MAG <= tmp2[i1][i2][dk] && tmp2[i1][i2][dk] <= 5 * CENTER_MAG);
            }
            // |tmp2| ~ 2.5Q
          }
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        for (int dk = 0; dk < 8; dk++) {
          tmp1[i1][0][dk] = tmp2[i1][0][dk] + tmp2[i1][1][dk];
          tmp1[i1][1][dk] = tmp2[i1][0][dk] - tmp2[i1][1][dk];
          assert(-10 * CENTER_MAG <= tmp1[i1][0][dk] && tmp1[i1][0][dk] <= 10 * CENTER_MAG);
          assert(-10 * CENTER_MAG <= tmp1[i1][1][dk] && tmp1[i1][1][dk] <= 10 * CENTER_MAG);
        }
        // |tmp1| ~ 5Q
      }

      for (int i = 0; i < 10; i++) {
        for (int dk = 0; dk < 8; dk++) {
          ntt[7 * i % 10][j][k0 + dk] = tmp1[i % 5][i % 2][dk];
        }
        // |ntt| ~ 3Q
      }

    }
  }

}

void ntt_9(int16_t ntt[10][9][16]) {

  for (int i = 0; i < 10; i++) {
    for (int k0 = 0; k0 < 16; k0 += 8) {

      int16_t tmp1[3][3][8] = {};
      int16_t tmp2[3][3][8] = {};

      for (int j = 0; j < 9; j++) {
        for (int dk = 0; dk < 8; dk++) {
          tmp1[j / 3][j % 3][dk] = ntt[i][j][k0 + dk];
        }
      }

      for (int j2 = 0; j2 < 3; j2++) {
        for (int j1 = 0; j1 < 3; j1++) {
          for (int jj1 = 0; jj1 < 3; jj1++) {
            for (int dk = 0; dk < 8; dk++) {
              tmp2[j1][j2][dk] += center_lift(w_3s[j1 * jj1 % 3] * tmp1[jj1][j2][dk]);
              assert(-3 * CENTER_MAG <= tmp2[j1][j2][dk] && tmp2[j1][j2][dk] <= 3 * CENTER_MAG);
            }
            // |tmp2| ~ 1.5Q
          }
        }
      }

      std::memset(tmp1, 0, sizeof(tmp1));

      for (int j1 = 0; j1 < 3; j1++) {
        for (int j2 = 0; j2 < 3; j2++) {
          for (int jj2 = 0; jj2 < 3; jj2++) {
            for (int dk = 0; dk < 8; dk++) {
              tmp1[j1][j2][dk] += center_lift(w_9s[(j1 + 3 * j2) * jj2 % 9] * tmp2[j1][jj2][dk]);
              assert(-3 * CENTER_MAG <= tmp1[j1][j2][dk] && tmp1[j1][j2][dk] <= 3 * CENTER_MAG);
            }
            // |tmp1| ~ 1.5Q
          }
        }
      }

      for (int j = 0; j < 9; j++) {
        for (int dk = 0; dk < 8; dk++) {
          ntt[i][j][k0 + dk] = tmp1[j % 3][j / 3][dk];
        }
        // |ntt| ~ 1.5Q
      }

    }
  }

}

void intt_9_x9(int16_t ntt[10][9][16]) {

  for (int i = 0; i < 10; i++) {
    for (int k0 = 0; k0 < 16; k0 += 8) {

      int16_t tmp1[3][3][8] = {};
      int16_t tmp2[3][3][8] = {};

      for (int j = 0; j < 9; j++) {
        for (int dk = 0; dk < 8; dk++) {
          tmp1[j / 3][j % 3][dk] = ntt[i][(9 - j) % 9][k0 + dk];
        }
      }

      for (int j2 = 0; j2 < 3; j2++) {
        for (int j1 = 0; j1 < 3; j1++) {
          for (int jj1 = 0; jj1 < 3; jj1++) {
            for (int dk = 0; dk < 8; dk++) {
              tmp2[j1][j2][dk] += center_lift(w_3s[j1 * jj1 % 3] * tmp1[jj1][j2][dk]);
              assert(-3 * CENTER_MAG <= tmp2[j1][j2][dk] && tmp2[j1][j2][dk] <= 3 * CENTER_MAG);
            }
            // |tmp2| ~ 1.5Q
          }
        }
      }

      std::memset(tmp1, 0, sizeof(tmp1));

      for (int j1 = 0; j1 < 3; j1++) {
        for (int j2 = 0; j2 < 3; j2++) {
          for (int jj2 = 0; jj2 < 3; jj2++) {
            for (int dk = 0; dk < 8; dk++) {
              tmp1[j1][j2][dk] += center_lift(w_9s[(j1 + 3 * j2) * jj2 % 9] * tmp2[j1][jj2][dk]);
              assert(-3 * CENTER_MAG <= tmp1[j1][j2][dk] && tmp1[j1][j2][dk] <= 3 * CENTER_MAG);
            }
            // |tmp1| ~ 1.5Q
          }
        }
      }

      for (int j = 0; j < 9; j++) {
        for (int dk = 0; dk < 8; dk++) {
          ntt[i][j][k0 + dk] = tmp1[j % 3][j / 3][dk];
        }
        // |ntt| ~ 1.5Q
      }

    }
  }

}


void base_mul(int16_t in1_ntt[10][9][16], int16_t in2_ntt[10][9][16], int16_t out_ntt[10][9][16]) {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      // TODO outer product
      int64_t tmp[10][9][16] = {};
      int64_t twiddle = center_lift(int32_t(1) * w_10s[i] * w_9s[j]);
      for (int k1 = 0; k1 < 16; k1++) {
        for (int k2 = 0; k2 < 16; k2++) {
          if (k1 + k2 >= 16) {
            tmp[i][j][k1 + k2 - 16] += twiddle * in1_ntt[i][j][k1] * in2_ntt[i][j][k2];
          }
          else {
            tmp[i][j][k1 + k2] += int64_t(1) * in1_ntt[i][j][k1] * in2_ntt[i][j][k2];
          }
        }
      }
      for (int k = 0; k < 16; k++) {
        out_ntt[i][j][k] = center_lift(tmp[i][j][k]);
      }
    }
  }
}

void forward(int16_t in_poly[P], int16_t out_ntt[10][9][16]) {
  for (int i = 0; i < 6 * 9 * 16; i++) {
    int t = i >> 4;
    out_ntt[t % 10][t % 9][i % 16] = i < P ? in_poly[i] : 0;
  }
  ntt_10(out_ntt);
  // print_ntt(out_ntt);
  ntt_9(out_ntt);
  // print_ntt(out_ntt);
}

constexpr int64_t INV90 = -51;
void backward(int16_t in_ntt[10][9][16], int16_t out_main[1440]) {
  intt_9_x9(in_ntt);
  intt_10_x10(in_ntt);
  for (int i = 0; i < 1440; i++) {
    int t = i >> 4;
    out_main[i] = center_lift(in_ntt[t % 10][t % 9][i % 16] * INV90);
  }
}

void mult_low(int16_t in1_low[81], int16_t in2_low[81], int16_t out_low[81]) {
  std::memset(out_low, 0, sizeof(int16_t[81]));
  for (int i = 0; i < 81; i++) {
    for (int j = 0; i + j < 81; j++) {
      out_low[i + j] = center_lift(out_low[i + j] + int64_t(1) * in1_low[i] * in2_low[j]);
    }
  }
}

void mult(int16_t in1_poly[P], int16_t in2_poly[P], int16_t out_poly[P]) {
  int16_t in1_ntt[10][9][16];
  int16_t in2_ntt[10][9][16];
  int16_t out_ntt[10][9][16];
  int16_t out_low[81];
  int16_t out_main[1440];

  forward(in1_poly, in1_ntt);
  forward(in2_poly, in2_ntt);
  mult_low(in1_poly, in2_poly, out_low);
  base_mul(in1_ntt, in2_ntt, out_ntt);
  backward(out_ntt, out_main);
  // for (int i = 0; i < 1440; i++) {
  //   std::cout << out_main[i] << " \n"[i == 1439];
  // }
  for (int i = 0; i < 761; i++) {
    out_poly[i] = out_main[(i + 760) % 1440] + out_main[(i + 761) % 1440];
  }
  out_poly[0] -= out_main[760];
  out_poly[760] -= out_main[81];
  for (int i = 81; i < 761; i++) {
    out_poly[i] += out_main[i];
  }
  for (int i = 680; i < 760; i++) {
    out_poly[i] -= out_low[i - 680] + out_low[i - 679];
  }
  out_poly[679] -= out_low[0];
  out_poly[760] -= out_low[80];
  for (int i = 0; i < 80; i++) {
    out_poly[i] += out_low[i];
  }
  out_poly[80] += out_low[80];
  for (int i = 0; i < 761; i++) {
    out_poly[i] = center_lift(out_poly[i]);
  }
}
