#include <arm_neon.h>
#include <iostream>
#include <cstring>
#include <cassert>
#include <array>
#include <cmath>

#include "mult.h"
#include "pretty.h"

template <size_t SZ>
constexpr std::array<int16_t, SZ> gen_pows(int16_t base) {
  std::array<int16_t, SZ> pows = {1};
  for (int i = 1; i < SZ; i++) {
    pows[i] = center_lift(int32_t(pows[i - 1]) * base);
  }
  return pows;
}

constexpr int32_t BARRET_Q = 467759;
constexpr int32_t BARRET_xW_9_Q = -751221681;
constexpr int32_t BARRET_xW_9_p2_Q = -423790064;
constexpr int32_t BARRET_xW_9_p4_Q = -445774759;
constexpr int32_t BARRET_TWIST_Q[3][3] = {
  {BARRET_Q, BARRET_Q, BARRET_Q},
  {BARRET_Q, BARRET_xW_9_Q, BARRET_xW_9_p2_Q},
  {BARRET_Q, BARRET_xW_9_p2_Q, BARRET_xW_9_p4_Q}
};

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;

constexpr std::array<int16_t, ORD> ws = gen_pows<ORD>(W_4590);

constexpr int16_t W_10 = ws[ORD / 10];
constexpr int16_t W_5 = ws[ORD / 5];
constexpr int16_t W_9 = ws[ORD / 9];
constexpr int16_t W_3 = ws[ORD / 3];

constexpr std::array<int16_t, 17> W_5S = gen_pows<17>(W_5);
constexpr std::array<int16_t, 10> W_10S = gen_pows<10>(W_10);
constexpr std::array<int16_t, 17> W_9S = gen_pows<17>(W_9);
constexpr std::array<int16_t, 5> W_3S = gen_pows<5>(W_3);

constexpr int16_t gen_bar(int16_t coef) {
  return std::round(double(coef) * (1 << 15) / Q);
}

template <size_t SZ>
constexpr std::array<int16_t, SZ> gen_bars(std::array<int16_t, SZ> arr) {
  std::array<int16_t, SZ> bars = {};
  for (int i = 0; i < SZ; i++) {
    bars[i] = gen_bar(arr[i]);
  }
  return bars;
}

constexpr int16_t ONE_BAR = gen_bar(1);
constexpr std::array<int16_t, 17> W_5_BARS = gen_bars<17>(W_5S);
constexpr std::array<int16_t, 10> W_10_BARS = gen_bars<10>(W_10S);
constexpr std::array<int16_t, 17> W_9_BARS = gen_bars<17>(W_9S);
constexpr std::array<int16_t, 5> W_3_BARS = gen_bars<5>(W_3S);

void debug_int16x8(int16x8_t val) {
  int16_t tmp[8];
  vst1q_s16(&tmp[0], val);
  for (int i = 0; i < 8; i++) {
    std::cerr << tmp[i] << " \n"[i == 7];
  }
}

void debug_int32x4(int32x4_t val) {
  int32_t tmp[4];
  vst1q_s32(&tmp[0], val);
  for (int i = 0; i < 4; i++) {
    std::cerr << tmp[i] << " \n"[i == 3];
  }
}

void barret_reduce(int16x8_t &v) {
  int16x8_t esti = vqrdmulhq_n_s16(v, ONE_BAR);
  v = vmlsq_n_s16(v, esti, Q);
}

int16x8_t barret_mul(int16x8_t v, int16_t coef, int16_t bar) {
  int16x8_t esti = vqrdmulhq_n_s16(v, bar);
  int16x8_t res = vmulq_n_s16(v, coef);
  res = vmlsq_n_s16(res, esti, Q);
  return res;
}

void barret_mla(int16x8_t &vd, int16x8_t v1, int16_t coef, int16_t bar) {
  int16x8_t esti = vqrdmulhq_n_s16(v1, bar);
  vd = vmlaq_n_s16(vd, v1, coef);
  vd = vmlsq_n_s16(vd, esti, Q);
}

void ntt_10(int16_t ntt[10][9][16], int16_t poly[1440]) {

  for (int j = 0; j < 9; j++) {
    {
      int16x8_t tmp_front[5][2];
      for (int i1 = 0; i1 < 5; i1++) {
        for (int i2 = 0; i2 < 2; i2++) {
          tmp_front[i1][i2] = vdupq_n_s16(0);
        }
      }

      for (int i = j, t = 0, idx = j; t < (j <= 2 ? 6 : 5); t++, i = (i + 9) % 10, idx += 9) {
        int ii1 = i % 5;
        int i2 = i % 2;
        int16x8_t fi_front = vld1q_s16(&poly[idx * 16]);
        if (i2 == 0) {
          tmp_front[0][0] = vaddq_s16(tmp_front[0][0], fi_front);
          barret_mla(tmp_front[1][0], fi_front, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla(tmp_front[2][0], fi_front, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla(tmp_front[3][0], fi_front, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla(tmp_front[4][0], fi_front, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
        }
        else {
          tmp_front[0][1] = vaddq_s16(tmp_front[0][1], fi_front);
          barret_mla(tmp_front[1][1], fi_front, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla(tmp_front[2][1], fi_front, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla(tmp_front[3][1], fi_front, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla(tmp_front[4][1], fi_front, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        vst1q_s16(&ntt[2 * i1][j][0], vaddq_s16(tmp_front[i1][0], tmp_front[i1][1]));
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][0], vsubq_s16(tmp_front[i1][0], tmp_front[i1][1]));
      }
    }

    {
      int16x8_t tmp_back[5][2];
      for (int i1 = 0; i1 < 5; i1++) {
        for (int i2 = 0; i2 < 2; i2++) {
          tmp_back[i1][i2] = vdupq_n_s16(0);
        }
      }

      for (int i = j, t = 0, idx = j; t < (j <= 2 ? 6 : 5); t++, i = (i + 9) % 10, idx += 9) {
        int ii1 = i % 5;
        int i2 = i % 2;
        int16x8_t fi_back = vld1q_s16(&poly[idx * 16 + 8]);
        if (i2 == 0) {
          tmp_back[0][0] = vaddq_s16(tmp_back[0][0], fi_back);
          barret_mla(tmp_back[1][0], fi_back, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla(tmp_back[2][0], fi_back, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla(tmp_back[3][0], fi_back, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla(tmp_back[4][0], fi_back, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
        }
        else {
          tmp_back[0][1] = vaddq_s16(tmp_back[0][1], fi_back);
          barret_mla(tmp_back[1][1], fi_back, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla(tmp_back[2][1], fi_back, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla(tmp_back[3][1], fi_back, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla(tmp_back[4][1], fi_back, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        vst1q_s16(&ntt[2 * i1][j][8], vaddq_s16(tmp_back[i1][0], tmp_back[i1][1]));
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][8], vsubq_s16(tmp_back[i1][0], tmp_back[i1][1]));
      }
    }

    // std::cerr << "ntt_10, j = " << j << '\n';
    // for (int i = 0; i < 10; i++) {
    //   std::cerr << "  i = " << i << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }
  }

}

void intt_10_x10(int16_t ntt[10][9][16]) {

  for (int j = 0; j < 9; j++) {
    {
      int16x8_t tmp_front[5][2];

      for (int i2 = 0; i2 < 2; i2++) {
        int16x8_t f0_front = vld1q_s16(&ntt[5 * i2][j][0]);
        tmp_front[0][i2] = f0_front;
        tmp_front[1][i2] = f0_front;
        tmp_front[2][i2] = f0_front;
        tmp_front[3][i2] = f0_front;
        tmp_front[4][i2] = f0_front;

        int16x8_t f1_front = vld1q_s16(&ntt[(5 * i2 + 4) % 10][j][0]);
        tmp_front[0][i2] = vaddq_s16(tmp_front[0][i2], f1_front);
        barret_mla(tmp_front[1][i2], f1_front, W_5S[1], W_5_BARS[1]);
        barret_mla(tmp_front[2][i2], f1_front, W_5S[2], W_5_BARS[2]);
        barret_mla(tmp_front[3][i2], f1_front, W_5S[3], W_5_BARS[3]);
        barret_mla(tmp_front[4][i2], f1_front, W_5S[4], W_5_BARS[4]);

        int16x8_t f2_front = vld1q_s16(&ntt[(5 * i2 + 8) % 10][j][0]);
        tmp_front[0][i2] = vaddq_s16(tmp_front[0][i2], f2_front);
        barret_mla(tmp_front[1][i2], f2_front, W_5S[2], W_5_BARS[2]);
        barret_mla(tmp_front[2][i2], f2_front, W_5S[4], W_5_BARS[4]);
        barret_mla(tmp_front[3][i2], f2_front, W_5S[1], W_5_BARS[1]);
        barret_mla(tmp_front[4][i2], f2_front, W_5S[3], W_5_BARS[3]);

        int16x8_t f3_front = vld1q_s16(&ntt[(5 * i2 + 2) % 10][j][0]);
        tmp_front[0][i2] = vaddq_s16(tmp_front[0][i2], f3_front);
        barret_mla(tmp_front[1][i2], f3_front, W_5S[3], W_5_BARS[3]);
        barret_mla(tmp_front[2][i2], f3_front, W_5S[1], W_5_BARS[1]);
        barret_mla(tmp_front[3][i2], f3_front, W_5S[4], W_5_BARS[4]);
        barret_mla(tmp_front[4][i2], f3_front, W_5S[2], W_5_BARS[2]);

        int16x8_t f4_front = vld1q_s16(&ntt[(5 * i2 + 6) % 10][j][0]);
        tmp_front[0][i2] = vaddq_s16(tmp_front[0][i2], f4_front);
        barret_mla(tmp_front[1][i2], f4_front, W_5S[4], W_5_BARS[4]);
        barret_mla(tmp_front[2][i2], f4_front, W_5S[3], W_5_BARS[3]);
        barret_mla(tmp_front[3][i2], f4_front, W_5S[2], W_5_BARS[2]);
        barret_mla(tmp_front[4][i2], f4_front, W_5S[1], W_5_BARS[1]);
      }

      for (int i1 = 0; i1 < 5; i1++) {
        vst1q_s16(&ntt[2 * i1][j][0], vaddq_s16(tmp_front[i1][0], tmp_front[i1][1]));
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][0], vsubq_s16(tmp_front[i1][0], tmp_front[i1][1]));
      }
    }

    {
      int16x8_t tmp_back[5][2];

      for (int i2 = 0; i2 < 2; i2++) {
        int16x8_t f0_back = vld1q_s16(&ntt[5 * i2][j][8]);
        tmp_back[0][i2] = f0_back;
        tmp_back[1][i2] = f0_back;
        tmp_back[2][i2] = f0_back;
        tmp_back[3][i2] = f0_back;
        tmp_back[4][i2] = f0_back;

        int16x8_t f1_back = vld1q_s16(&ntt[(5 * i2 + 4) % 10][j][8]);
        tmp_back[0][i2] = vaddq_s16(tmp_back[0][i2], f1_back);
        barret_mla(tmp_back[1][i2], f1_back, W_5S[1], W_5_BARS[1]);
        barret_mla(tmp_back[2][i2], f1_back, W_5S[2], W_5_BARS[2]);
        barret_mla(tmp_back[3][i2], f1_back, W_5S[3], W_5_BARS[3]);
        barret_mla(tmp_back[4][i2], f1_back, W_5S[4], W_5_BARS[4]);

        int16x8_t f2_back = vld1q_s16(&ntt[(5 * i2 + 8) % 10][j][8]);
        tmp_back[0][i2] = vaddq_s16(tmp_back[0][i2], f2_back);
        barret_mla(tmp_back[1][i2], f2_back, W_5S[2], W_5_BARS[2]);
        barret_mla(tmp_back[2][i2], f2_back, W_5S[4], W_5_BARS[4]);
        barret_mla(tmp_back[3][i2], f2_back, W_5S[1], W_5_BARS[1]);
        barret_mla(tmp_back[4][i2], f2_back, W_5S[3], W_5_BARS[3]);

        int16x8_t f3_back = vld1q_s16(&ntt[(5 * i2 + 2) % 10][j][8]);
        tmp_back[0][i2] = vaddq_s16(tmp_back[0][i2], f3_back);
        barret_mla(tmp_back[1][i2], f3_back, W_5S[3], W_5_BARS[3]);
        barret_mla(tmp_back[2][i2], f3_back, W_5S[1], W_5_BARS[1]);
        barret_mla(tmp_back[3][i2], f3_back, W_5S[4], W_5_BARS[4]);
        barret_mla(tmp_back[4][i2], f3_back, W_5S[2], W_5_BARS[2]);

        int16x8_t f4_back = vld1q_s16(&ntt[(5 * i2 + 6) % 10][j][8]);
        tmp_back[0][i2] = vaddq_s16(tmp_back[0][i2], f4_back);
        barret_mla(tmp_back[1][i2], f4_back, W_5S[4], W_5_BARS[4]);
        barret_mla(tmp_back[2][i2], f4_back, W_5S[3], W_5_BARS[3]);
        barret_mla(tmp_back[3][i2], f4_back, W_5S[2], W_5_BARS[2]);
        barret_mla(tmp_back[4][i2], f4_back, W_5S[1], W_5_BARS[1]);
      }

      for (int i1 = 0; i1 < 5; i1++) {
        vst1q_s16(&ntt[2 * i1][j][8], vaddq_s16(tmp_back[i1][0], tmp_back[i1][1]));
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][8], vsubq_s16(tmp_back[i1][0], tmp_back[i1][1]));
        // std::cerr << i1 << ' ' << j << ": put to " << (((2 * i1) * 81 + 10 * j) % 90 * 16) << " and " << ((2 * i1 + 5) * 81 + 10 * j) % 90 * 16 << '\n';
      }
    }

    // std::cerr << "intt_10_x10, j = " << j << '\n';
    // for (int i = 0; i < 10; i++) {
    //   std::cerr << "  i = " << i << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

  }

}

void ntt_9(int16_t ntt[10][9][16]) {

  for (int i = 0; i < 10; i++) {
    // std::cerr << "ntt_9 input, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

    int16x8_t tmp_front[3][3];
    int16x8_t tmp_back[3][3];

    for (int j2 = 0; j2 < 3; j2++) {
      int16x8_t f0_front = vld1q_s16(&ntt[i][j2][0]);
      int16x8_t f0_back = vld1q_s16(&ntt[i][j2][8]);
      int16x8_t f1_front = vld1q_s16(&ntt[i][3 + j2][0]);
      int16x8_t f1_back = vld1q_s16(&ntt[i][3 + j2][8]);
      int16x8_t f2_front = vld1q_s16(&ntt[i][6 + j2][0]);
      int16x8_t f2_back = vld1q_s16(&ntt[i][6 + j2][8]);

      barret_reduce(f0_front);
      barret_reduce(f0_back);
      barret_reduce(f1_front);
      barret_reduce(f1_back);
      barret_reduce(f2_front);
      barret_reduce(f2_back);

      int16x8_t f1_f2_front = vaddq_s16(f1_front, f2_front);
      int16x8_t f1_f2_back = vaddq_s16(f1_back, f2_back);
      barret_mla(f1_front, f2_front, W_3S[1], W_3_BARS[1]);
      barret_mla(f1_back, f2_back, W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_front = barret_mul(f1_front, W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_back = barret_mul(f1_back, W_3S[1], W_3_BARS[1]);

      tmp_front[0][j2] = vaddq_s16(f0_front, f1_f2_front);
      tmp_back[0][j2] = vaddq_s16(f0_back, f1_f2_back);
      tmp_front[1][j2] = vaddq_s16(f0_front, w3f1_w32f2_front);
      tmp_back[1][j2] = vaddq_s16(f0_back, w3f1_w32f2_back);
      tmp_front[2][j2] = vsubq_s16(vsubq_s16(f0_front, f1_f2_front), w3f1_w32f2_front);
      tmp_back[2][j2] = vsubq_s16(vsubq_s16(f0_back, f1_f2_back), w3f1_w32f2_back);
    }

    tmp_front[1][1] = barret_mul(tmp_front[1][1], W_9S[1], W_9_BARS[1]);
    tmp_back[1][1] = barret_mul(tmp_back[1][1], W_9S[1], W_9_BARS[1]);
    tmp_front[2][1] = barret_mul(tmp_front[2][1], W_9S[2], W_9_BARS[2]);
    tmp_back[2][1] = barret_mul(tmp_back[2][1], W_9S[2], W_9_BARS[2]);

    tmp_front[1][2] = barret_mul(tmp_front[1][2], W_9S[2], W_9_BARS[2]);
    tmp_back[1][2] = barret_mul(tmp_back[1][2], W_9S[2], W_9_BARS[2]);
    tmp_front[2][2] = barret_mul(tmp_front[2][2], W_9S[4], W_9_BARS[4]);
    tmp_back[2][2] = barret_mul(tmp_back[2][2], W_9S[4], W_9_BARS[4]);

    for (int j1 = 0; j1 < 3; j1++) {
      int16x8_t f1_f2_front = vaddq_s16(tmp_front[j1][1], tmp_front[j1][2]);
      int16x8_t f1_f2_back = vaddq_s16(tmp_back[j1][1], tmp_back[j1][2]);
      barret_mla(tmp_front[j1][1], tmp_front[j1][2], W_3S[1], W_3_BARS[1]);
      barret_mla(tmp_back[j1][1], tmp_back[j1][2], W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_front = barret_mul(tmp_front[j1][1], W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_back = barret_mul(tmp_back[j1][1], W_3S[1], W_3_BARS[1]);

      tmp_front[j1][1] = vaddq_s16(tmp_front[j1][0], w3f1_w32f2_front);
      tmp_back[j1][1] = vaddq_s16(tmp_back[j1][0], w3f1_w32f2_back);
      tmp_front[j1][2] = vsubq_s16(vsubq_s16(tmp_front[j1][0], f1_f2_front), w3f1_w32f2_front);
      tmp_back[j1][2] = vsubq_s16(vsubq_s16(tmp_back[j1][0], f1_f2_back), w3f1_w32f2_back);
      tmp_front[j1][0] = vaddq_s16(tmp_front[j1][0], f1_f2_front);
      tmp_back[j1][0] = vaddq_s16(tmp_back[j1][0], f1_f2_back);

      vst1q_s16(&ntt[i][j1][0], tmp_front[j1][0]);
      vst1q_s16(&ntt[i][j1][8], tmp_back[j1][0]);
      vst1q_s16(&ntt[i][j1 + 3][0], tmp_front[j1][1]);
      vst1q_s16(&ntt[i][j1 + 3][8], tmp_back[j1][1]);
      vst1q_s16(&ntt[i][j1 + 6][0], tmp_front[j1][2]);
      vst1q_s16(&ntt[i][j1 + 6][8], tmp_back[j1][2]);
    }

    // std::cerr << "ntt_9 output, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }
  }

}

void intt_9_x9(int16_t ntt[10][9][16], int16_t poly[1440]) {

  for (int i = 0; i < 10; i++) {
    // std::cerr << "intt_9 input, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

    int16x8_t tmp_front[3][3];
    int16x8_t tmp_back[3][3];

    for (int j2 = 0; j2 < 3; j2++) {
      int16x8_t f0_front = vld1q_s16(&ntt[i][(9 - j2) % 9][0]);
      int16x8_t f0_back = vld1q_s16(&ntt[i][(9 - j2) % 9][8]);
      int16x8_t f1_front = vld1q_s16(&ntt[i][6 - j2][0]);
      int16x8_t f1_back = vld1q_s16(&ntt[i][6 - j2][8]);
      int16x8_t f2_front = vld1q_s16(&ntt[i][3 - j2][0]);
      int16x8_t f2_back = vld1q_s16(&ntt[i][3 - j2][8]);

      barret_reduce(f0_front);
      barret_reduce(f0_back);
      barret_reduce(f1_front);
      barret_reduce(f1_back);
      barret_reduce(f2_front);
      barret_reduce(f2_back);

      int16x8_t f1_f2_front = vaddq_s16(f1_front, f2_front);
      int16x8_t f1_f2_back = vaddq_s16(f1_back, f2_back);
      barret_mla(f1_front, f2_front, W_3S[1], W_3_BARS[1]);
      barret_mla(f1_back, f2_back, W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_front = barret_mul(f1_front, W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_back = barret_mul(f1_back, W_3S[1], W_3_BARS[1]);

      tmp_front[0][j2] = vaddq_s16(f0_front, f1_f2_front);
      tmp_back[0][j2] = vaddq_s16(f0_back, f1_f2_back);
      tmp_front[1][j2] = vaddq_s16(f0_front, w3f1_w32f2_front);
      tmp_back[1][j2] = vaddq_s16(f0_back, w3f1_w32f2_back);
      tmp_front[2][j2] = vsubq_s16(vsubq_s16(f0_front, f1_f2_front), w3f1_w32f2_front);
      tmp_back[2][j2] = vsubq_s16(vsubq_s16(f0_back, f1_f2_back), w3f1_w32f2_back);
    }

    tmp_front[1][1] = barret_mul(tmp_front[1][1], W_9S[1], W_9_BARS[1]);
    tmp_back[1][1] = barret_mul(tmp_back[1][1], W_9S[1], W_9_BARS[1]);
    tmp_front[2][1] = barret_mul(tmp_front[2][1], W_9S[2], W_9_BARS[2]);
    tmp_back[2][1] = barret_mul(tmp_back[2][1], W_9S[2], W_9_BARS[2]);

    tmp_front[1][2] = barret_mul(tmp_front[1][2], W_9S[2], W_9_BARS[2]);
    tmp_back[1][2] = barret_mul(tmp_back[1][2], W_9S[2], W_9_BARS[2]);
    tmp_front[2][2] = barret_mul(tmp_front[2][2], W_9S[4], W_9_BARS[4]);
    tmp_back[2][2] = barret_mul(tmp_back[2][2], W_9S[4], W_9_BARS[4]);

    for (int j1 = 0; j1 < 3; j1++) {
      int16x8_t f1_f2_front = vaddq_s16(tmp_front[j1][1], tmp_front[j1][2]);
      int16x8_t f1_f2_back = vaddq_s16(tmp_back[j1][1], tmp_back[j1][2]);
      barret_mla(tmp_front[j1][1], tmp_front[j1][2], W_3S[1], W_3_BARS[1]);
      barret_mla(tmp_back[j1][1], tmp_back[j1][2], W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_front = barret_mul(tmp_front[j1][1], W_3S[1], W_3_BARS[1]);
      int16x8_t w3f1_w32f2_back = barret_mul(tmp_back[j1][1], W_3S[1], W_3_BARS[1]);

      tmp_front[j1][1] = vaddq_s16(tmp_front[j1][0], w3f1_w32f2_front);
      tmp_back[j1][1] = vaddq_s16(tmp_back[j1][0], w3f1_w32f2_back);
      tmp_front[j1][2] = vsubq_s16(vsubq_s16(tmp_front[j1][0], f1_f2_front), w3f1_w32f2_front);
      tmp_back[j1][2] = vsubq_s16(vsubq_s16(tmp_back[j1][0], f1_f2_back), w3f1_w32f2_back);
      tmp_front[j1][0] = vaddq_s16(tmp_front[j1][0], f1_f2_front);
      tmp_back[j1][0] = vaddq_s16(tmp_back[j1][0], f1_f2_back);

      vst1q_s16(&poly[(81 * i + 10 * j1) % 90 * 16], tmp_front[j1][0]);
      vst1q_s16(&poly[(81 * i + 10 * j1) % 90 * 16 + 8], tmp_back[j1][0]);
      vst1q_s16(&poly[(81 * i + 10 * (j1 + 3)) % 90 * 16], tmp_front[j1][1]);
      vst1q_s16(&poly[(81 * i + 10 * (j1 + 3)) % 90 * 16 + 8], tmp_back[j1][1]);
      vst1q_s16(&poly[(81 * i + 10 * (j1 + 6)) % 90 * 16], tmp_front[j1][2]);
      vst1q_s16(&poly[(81 * i + 10 * (j1 + 6)) % 90 * 16 + 8], tmp_back[j1][2]);
    }

    // std::cerr << "intt_9 output, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }
  }

}

void base_mul(int16_t in1_ntt[10][9][16], int16_t in2_ntt[10][9][16], int16_t out_ntt[10][9][16]) {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      // TODO outer product
      int64_t tmp[10][9][16] = {};
      int64_t twiddle = center_lift(int32_t(1) * W_10S[i] * W_9S[j]);
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

void forward(int16_t in_poly[], int16_t out_ntt[10][9][16]) {
  ntt_10(out_ntt, in_poly);
  // print_ntt(out_ntt);
  ntt_9(out_ntt);
  // print_ntt(out_ntt);
}


constexpr int16_t INV90 = -51;
constexpr int16_t INV90_BAR = gen_bar(INV90);
void div90_main(int16_t main_poly[1440]) {
  int16x8_t zeros = vdupq_n_s16(0);
  for (int i = 0; i < 1440; i += 8) {
    int16x8_t chunk = vld1q_s16(&main_poly[i]);
    chunk = barret_mul(chunk, INV90, INV90_BAR);
    vst1q_s16(&main_poly[i], chunk);
  }
}

void backward(int16_t in_ntt[10][9][16], int16_t out_main[1440]) {
  intt_10_x10(in_ntt);
  intt_9_x9(in_ntt, out_main);
  div90_main(out_main);
  // for (int i = 0; i < 1440; i++) {
  //   int t = i >> 4;
  //   out_main[i] = center_lift(in_ntt[t % 10][t % 9][i % 16] * INV90);
  // }
}

void mult_low(int16_t in1_low[81], int16_t in2_low[81], int16_t out_low[81]) {
  std::memset(out_low, 0, sizeof(int16_t[81]));
  for (int i = 0; i < 81; i++) {
    for (int j = 0; i + j < 81; j++) {
      out_low[i + j] = center_lift(out_low[i + j] + int64_t(1) * in1_low[i] * in2_low[j]);
    }
  }
}

// main_poly length must >= 1448
// poly length must >= 768
void crt(int16_t poly[], int16_t main_poly[], int16_t low[]) {
  main_poly[1440] = main_poly[0];
  for (int i = 0; i < 680; i += 8) {
    int16x8_t a = vld1q_s16(&main_poly[i + 760]);
    int16x8_t b = vld1q_s16(&main_poly[i + 761]);
    int16x8_t c;
    c = vaddq_s16(a, b);
    vst1q_s16(&poly[i], c);
  }
  for (int i = 680; i < 760; i += 8) {
    int16x8_t a = vld1q_s16(&main_poly[i - 680]);
    int16x8_t b = vld1q_s16(&main_poly[i - 679]);
    int16x8_t c;
    c = vaddq_s16(a, b);
    vst1q_s16(&poly[i], c);
  }
  poly[760] = main_poly[80];
  poly[0] -= main_poly[760];
  // for (int i = 81; i < 761; i++) {
  //   poly[i] += main_poly[i];
  // }
  for (int i = 81; i < 761; i += 8) {
    int16x8_t a = vld1q_s16(&main_poly[i]);
    int16x8_t c = vld1q_s16(&poly[i]);
    c = vaddq_s16(c, a);
    vst1q_s16(&poly[i], c);
  }
  for (int i = 680; i < 760; i += 8) {
    int16x8_t a = vld1q_s16(&low[i - 680]);
    int16x8_t b = vld1q_s16(&low[i - 679]);
    a = vaddq_s16(a, b);
    int16x8_t c = vld1q_s16(&poly[i]);
    c = vsubq_s16(c, a);
    vst1q_s16(&poly[i], c);
  }
  poly[679] -= low[0];
  poly[760] -= low[80];
  // for (int i = 0; i < 80; i++) {
  //   poly[i] += low[i];
  // }
  for (int i = 0; i < 80; i += 8) {
    int16x8_t a = vld1q_s16(&low[i]);
    int16x8_t c = vld1q_s16(&poly[i]);
    c = vaddq_s16(c, a);
    vst1q_s16(&poly[i], c);
  }
  poly[80] += low[80];
}

// poly length must >= 768
void center_poly(int16_t poly[]) {
  int16x8_t zeros = vdupq_n_s16(0);
  for (int i = 0; i < 768; i += 8) {
    int16x8_t chunk = vld1q_s16(&poly[i]);
    // debug_int16x8(chunk);
    int32x4_t wide_low = vmovl_s16(vget_low_s16(chunk));
    int32x4_t wide_high = vmovl_high_s16(chunk);
    // debug_int32x4(wide_low);
    // debug_int32x4(wide_high);

    int32x4_t esti;

    esti = vqrdmulhq_n_s32(wide_low, BARRET_Q);
    wide_low = vmlsq_n_s32(wide_low, esti, Q);

    esti = vqrdmulhq_n_s32(wide_high, BARRET_Q);
    wide_high = vmlsq_n_s32(wide_high, esti, Q);

    chunk = vuzp1q_s16(vreinterpretq_s16_s32(wide_low), vreinterpretq_s16_s32(wide_high));
    vst1q_s16(&poly[i], chunk);
  }
}

// inx_poly length must >= 768
// inx_poly[761 : 768] filled with zero
// out_poly length must >= 768
void mult(int16_t in1_poly[], int16_t in2_poly[], int16_t out_poly[]) {
  int16_t in1_ntt[10][9][16];
  int16_t in2_ntt[10][9][16];
  int16_t out_ntt[10][9][16];
  int16_t out_low[81];
  int16_t out_main[1448];

  forward(in1_poly, in1_ntt);
  forward(in2_poly, in2_ntt);
  mult_low(in1_poly, in2_poly, out_low);
  base_mul(in1_ntt, in2_ntt, out_ntt);
  backward(out_ntt, out_main);
  // for (int i = 0; i < 1440; i++) {
  //   std::cout << out_main[i] << " \n"[i == 1439];
  // }
  crt(out_poly, out_main, out_low);
  center_poly(out_poly);
  // for (int i = 0; i < 761; i++) {
  //   out_poly[i] = center_lift(out_poly[i]);
  // }
}
