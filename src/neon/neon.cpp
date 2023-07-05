#include <arm_neon.h>
#include <cstring>
#include <iostream>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"
#include "utils/debug.h"

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;

constexpr int16_t W_10 = gen_pow<int16_t, Q>(W_4590, ORD / 10);
constexpr int16_t W_5 = gen_pow<int16_t, Q>(W_4590, ORD / 5);
constexpr int16_t W_9 = gen_pow<int16_t, Q>(W_4590, ORD / 9);
constexpr int16_t W_3 = gen_pow<int16_t, Q>(W_4590, ORD / 3);

constexpr std::array<int16_t, 17> W_5S = gen_pows<int16_t, 17, Q>(W_5);
constexpr std::array<int16_t, 10> W_10S = gen_pows<int16_t, 10, Q>(W_10);
constexpr std::array<int16_t, 17> W_9S = gen_pows<int16_t, 17, Q>(W_9);
constexpr std::array<int16_t, 5> W_3S = gen_pows<int16_t, 5, Q>(W_3);

constexpr std::array<int16_t, 17> W_5_BARS = gen_bars<int16_t, 17, Q>(W_5S);

void ntt_10(int16_t ntt[10][9][16], int16_t poly[1440]) {

  for (int j = 0; j < 9; j++) {
    // std::cerr << "ntt_10 input, j = " << j << '\n';
    // for (int i = 0; i < 10; i++) {
    //   std::cerr << "  i = " << i << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

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
          barret_mla<Q>(tmp_front[1][0], fi_front, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla<Q>(tmp_front[2][0], fi_front, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla<Q>(tmp_front[3][0], fi_front, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla<Q>(tmp_front[4][0], fi_front, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
        }
        else {
          tmp_front[0][1] = vaddq_s16(tmp_front[0][1], fi_front);
          barret_mla<Q>(tmp_front[1][1], fi_front, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla<Q>(tmp_front[2][1], fi_front, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla<Q>(tmp_front[3][1], fi_front, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla<Q>(tmp_front[4][1], fi_front, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
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
          barret_mla<Q>(tmp_back[1][0], fi_back, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla<Q>(tmp_back[2][0], fi_back, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla<Q>(tmp_back[3][0], fi_back, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla<Q>(tmp_back[4][0], fi_back, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
        }
        else {
          tmp_back[0][1] = vaddq_s16(tmp_back[0][1], fi_back);
          barret_mla<Q>(tmp_back[1][1], fi_back, W_5S[ii1], W_5_BARS[ii1]);
          barret_mla<Q>(tmp_back[2][1], fi_back, W_5S[2 * ii1], W_5_BARS[2 * ii1]);
          barret_mla<Q>(tmp_back[3][1], fi_back, W_5S[3 * ii1], W_5_BARS[3 * ii1]);
          barret_mla<Q>(tmp_back[4][1], fi_back, W_5S[4 * ii1], W_5_BARS[4 * ii1]);
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        vst1q_s16(&ntt[2 * i1][j][8], vaddq_s16(tmp_back[i1][0], tmp_back[i1][1]));
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][8], vsubq_s16(tmp_back[i1][0], tmp_back[i1][1]));
      }
    }

    // std::cerr << "ntt_10 output, j = " << j << '\n';
    // for (int i = 0; i < 10; i++) {
    //   std::cerr << "  i = " << i << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }
  }

}

constexpr int16_t W5_W5_4 = 502;
constexpr int16_t W5_mW5_4 = 459;
constexpr int16_t W5_2_W5_3 = -503;
constexpr int16_t W5_2_mW5_3 = 868;
constexpr int16_t W5_W5_2_mW5_3_mW_4 = 1327;

constexpr int16_t INV2 = -2295;

void intt_10_x10(int16_t ntt[10][9][16]) {

  for (int j = 0; j < 9; j++) {

    // std::cerr << "intt_10_x10 input, j = " << j << '\n';
    // for (int i = 0; i < 10; i++) {
    //   std::cerr << "  i = " << i << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

    {
      int16x8_t tmp_front[5][2];

      for (int i2 = 0; i2 < 2; i2++) {
        int16x8_t f0_front = vld1q_s16(&ntt[5 * i2][j][0]);
        int16x8_t f1_front = vld1q_s16(&ntt[(5 * i2 + 4) % 10][j][0]);
        int16x8_t f2_front = vld1q_s16(&ntt[(5 * i2 + 8) % 10][j][0]);
        int16x8_t f3_front = vld1q_s16(&ntt[(5 * i2 + 2) % 10][j][0]);
        int16x8_t f4_front = vld1q_s16(&ntt[(5 * i2 + 6) % 10][j][0]);

        int16x8_t f1_f4_front = vaddq_s16(f1_front, f4_front);
        int16x8_t f1_mf4_front = vsubq_s16(f1_front, f4_front);

        int16x8_t f3_f2_front = vaddq_s16(f3_front, f2_front);
        int16x8_t f3_mf2_front = vsubq_s16(f3_front, f2_front);

        int16x8_t f1_f2_f3_f4_front = vaddq_s16(f1_f4_front, f3_f2_front);
        int16x8_t f1_mf2_f3_mf4_front = vaddq_s16(f1_mf4_front, f3_mf2_front);

        tmp_front[0][i2] = vaddq_s16(f0_front, f1_f2_f3_f4_front);

        int16x8_t neg_c0 = barret_mul_const<Q, -W5_W5_4>(f1_f4_front);
        barret_mla_const<Q, -W5_2_W5_3>(neg_c0, f3_f2_front);
        int16x8_t neg_c1 = vsubq_s16(f1_f2_f3_f4_front, neg_c0);

        f1_mf4_front = barret_mul_const<Q, W5_mW5_4>(f1_mf4_front);
        f3_mf2_front = barret_mul_const<Q, W5_2_mW5_3>(f3_mf2_front);
        int16x8_t neg_n0 = vsubq_s16(f3_mf2_front, f1_mf4_front);
        int16x8_t neg_n1 = vaddq_s16(f3_mf2_front, f1_mf4_front);
        barret_mla_const<Q, -W5_W5_2_mW5_3_mW_4>(neg_n1, f1_mf2_f3_mf4_front);

        tmp_front[1][i2] = vaddq_s16(neg_c0, neg_n0);
        tmp_front[2][i2] = vaddq_s16(neg_c1, neg_n1);
        tmp_front[4][i2] = vsubq_s16(neg_c0, neg_n0);
        tmp_front[3][i2] = vsubq_s16(neg_c1, neg_n1);

        tmp_front[1][i2] = barret_mul_const<Q, -INV2>(tmp_front[1][i2]);
        tmp_front[2][i2] = barret_mul_const<Q, -INV2>(tmp_front[2][i2]);
        tmp_front[3][i2] = barret_mul_const<Q, -INV2>(tmp_front[3][i2]);
        tmp_front[4][i2] = barret_mul_const<Q, -INV2>(tmp_front[4][i2]);

        tmp_front[1][i2] = vaddq_s16(tmp_front[1][i2], f0_front);
        tmp_front[2][i2] = vaddq_s16(tmp_front[2][i2], f0_front);
        tmp_front[3][i2] = vaddq_s16(tmp_front[3][i2], f0_front);
        tmp_front[4][i2] = vaddq_s16(tmp_front[4][i2], f0_front);
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
        int16x8_t f1_back = vld1q_s16(&ntt[(5 * i2 + 4) % 10][j][8]);
        int16x8_t f2_back = vld1q_s16(&ntt[(5 * i2 + 8) % 10][j][8]);
        int16x8_t f3_back = vld1q_s16(&ntt[(5 * i2 + 2) % 10][j][8]);
        int16x8_t f4_back = vld1q_s16(&ntt[(5 * i2 + 6) % 10][j][8]);

        int16x8_t f1_f4_back = vaddq_s16(f1_back, f4_back);
        int16x8_t f1_mf4_back = vsubq_s16(f1_back, f4_back);

        int16x8_t f3_f2_back = vaddq_s16(f3_back, f2_back);
        int16x8_t f3_mf2_back = vsubq_s16(f3_back, f2_back);

        int16x8_t f1_f2_f3_f4_back = vaddq_s16(f1_f4_back, f3_f2_back);
        int16x8_t f1_mf2_f3_mf4_back = vaddq_s16(f1_mf4_back, f3_mf2_back);

        tmp_back[0][i2] = vaddq_s16(f0_back, f1_f2_f3_f4_back);

        int16x8_t neg_c0 = barret_mul_const<Q, -W5_W5_4>(f1_f4_back);
        barret_mla_const<Q, -W5_2_W5_3>(neg_c0, f3_f2_back);
        int16x8_t neg_c1 = vsubq_s16(f1_f2_f3_f4_back, neg_c0);

        f1_mf4_back = barret_mul_const<Q, W5_mW5_4>(f1_mf4_back);
        f3_mf2_back = barret_mul_const<Q, W5_2_mW5_3>(f3_mf2_back);
        int16x8_t neg_n0 = vsubq_s16(f3_mf2_back, f1_mf4_back);
        int16x8_t neg_n1 = vaddq_s16(f3_mf2_back, f1_mf4_back);
        barret_mla_const<Q, -W5_W5_2_mW5_3_mW_4>(neg_n1, f1_mf2_f3_mf4_back);

        tmp_back[1][i2] = vaddq_s16(neg_c0, neg_n0);
        tmp_back[2][i2] = vaddq_s16(neg_c1, neg_n1);
        tmp_back[4][i2] = vsubq_s16(neg_c0, neg_n0);
        tmp_back[3][i2] = vsubq_s16(neg_c1, neg_n1);

        tmp_back[1][i2] = barret_mul_const<Q, -INV2>(tmp_back[1][i2]);
        tmp_back[2][i2] = barret_mul_const<Q, -INV2>(tmp_back[2][i2]);
        tmp_back[3][i2] = barret_mul_const<Q, -INV2>(tmp_back[3][i2]);
        tmp_back[4][i2] = barret_mul_const<Q, -INV2>(tmp_back[4][i2]);

        tmp_back[1][i2] = vaddq_s16(tmp_back[1][i2], f0_back);
        tmp_back[2][i2] = vaddq_s16(tmp_back[2][i2], f0_back);
        tmp_back[3][i2] = vaddq_s16(tmp_back[3][i2], f0_back);
        tmp_back[4][i2] = vaddq_s16(tmp_back[4][i2], f0_back);
      }

      for (int i1 = 0; i1 < 5; i1++) {
        vst1q_s16(&ntt[2 * i1][j][8], vaddq_s16(tmp_back[i1][0], tmp_back[i1][1]));
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][8], vsubq_s16(tmp_back[i1][0], tmp_back[i1][1]));
      }
    }

    // std::cerr << "intt_10_x10 output, j = " << j << '\n';
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

      barret_reduce<Q>(f0_front);
      barret_reduce<Q>(f0_back);
      barret_reduce<Q>(f1_front);
      barret_reduce<Q>(f1_back);
      barret_reduce<Q>(f2_front);
      barret_reduce<Q>(f2_back);

      int16x8_t f1_f2_front = vaddq_s16(f1_front, f2_front);
      int16x8_t f1_f2_back = vaddq_s16(f1_back, f2_back);
      barret_mla_const<Q, W_3S[1]>(f1_front, f2_front);
      barret_mla_const<Q, W_3S[1]>(f1_back, f2_back);
      int16x8_t w3f1_w32f2_front = barret_mul_const<Q, W_3S[1]>(f1_front);
      int16x8_t w3f1_w32f2_back = barret_mul_const<Q, W_3S[1]>(f1_back);

      tmp_front[0][j2] = vaddq_s16(f0_front, f1_f2_front);
      tmp_back[0][j2] = vaddq_s16(f0_back, f1_f2_back);
      tmp_front[1][j2] = vaddq_s16(f0_front, w3f1_w32f2_front);
      tmp_back[1][j2] = vaddq_s16(f0_back, w3f1_w32f2_back);
      tmp_front[2][j2] = vsubq_s16(vsubq_s16(f0_front, f1_f2_front), w3f1_w32f2_front);
      tmp_back[2][j2] = vsubq_s16(vsubq_s16(f0_back, f1_f2_back), w3f1_w32f2_back);

    }

    // std::cerr << "before twist:\n";
    // debug_int16x8(tmp_front[0][0]);
    // debug_int16x8(tmp_front[1][0]);
    // debug_int16x8(tmp_front[2][0]);
    // debug_int16x8(tmp_front[0][1]);
    // debug_int16x8(tmp_front[1][1]);
    // debug_int16x8(tmp_front[2][1]);
    // debug_int16x8(tmp_front[0][2]);
    // debug_int16x8(tmp_front[1][2]);
    // debug_int16x8(tmp_front[2][2]);

    tmp_front[1][1] = barret_mul_const<Q, W_9S[1]>(tmp_front[1][1]);
    tmp_back[1][1] = barret_mul_const<Q, W_9S[1]>(tmp_back[1][1]);
    tmp_front[2][1] = barret_mul_const<Q, W_9S[2]>(tmp_front[2][1]);
    tmp_back[2][1] = barret_mul_const<Q, W_9S[2]>(tmp_back[2][1]);

    tmp_front[1][2] = barret_mul_const<Q, W_9S[2]>(tmp_front[1][2]);
    tmp_back[1][2] = barret_mul_const<Q, W_9S[2]>(tmp_back[1][2]);
    tmp_front[2][2] = barret_mul_const<Q, W_9S[4]>(tmp_front[2][2]);
    tmp_back[2][2] = barret_mul_const<Q, W_9S[4]>(tmp_back[2][2]);

    // std::cerr << "after twist:\n";
    // debug_int16x8(tmp_front[0][0]);
    // debug_int16x8(tmp_front[1][0]);
    // debug_int16x8(tmp_front[2][0]);
    // debug_int16x8(tmp_front[0][1]);
    // debug_int16x8(tmp_front[1][1]);
    // debug_int16x8(tmp_front[2][1]);
    // debug_int16x8(tmp_front[0][2]);
    // debug_int16x8(tmp_front[1][2]);
    // debug_int16x8(tmp_front[2][2]);

    for (int j1 = 0; j1 < 3; j1++) {
      int16x8_t f1_f2_front = vaddq_s16(tmp_front[j1][1], tmp_front[j1][2]);
      int16x8_t f1_f2_back = vaddq_s16(tmp_back[j1][1], tmp_back[j1][2]);
      barret_mla_const<Q, W_3S[1]>(tmp_front[j1][1], tmp_front[j1][2]);
      barret_mla_const<Q, W_3S[1]>(tmp_back[j1][1], tmp_back[j1][2]);
      int16x8_t w3f1_w32f2_front = barret_mul_const<Q, W_3S[1]>(tmp_front[j1][1]);
      int16x8_t w3f1_w32f2_back = barret_mul_const<Q, W_3S[1]>(tmp_back[j1][1]);

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

      barret_reduce<Q>(f0_front);
      barret_reduce<Q>(f0_back);
      barret_reduce<Q>(f1_front);
      barret_reduce<Q>(f1_back);
      barret_reduce<Q>(f2_front);
      barret_reduce<Q>(f2_back);

      int16x8_t f1_f2_front = vaddq_s16(f1_front, f2_front);
      int16x8_t f1_f2_back = vaddq_s16(f1_back, f2_back);
      barret_mla_const<Q, W_3S[1]>(f1_front, f2_front);
      barret_mla_const<Q, W_3S[1]>(f1_back, f2_back);
      int16x8_t w3f1_w32f2_front = barret_mul_const<Q, W_3S[1]>(f1_front);
      int16x8_t w3f1_w32f2_back = barret_mul_const<Q, W_3S[1]>(f1_back);

      tmp_front[0][j2] = vaddq_s16(f0_front, f1_f2_front);
      tmp_back[0][j2] = vaddq_s16(f0_back, f1_f2_back);
      tmp_front[1][j2] = vaddq_s16(f0_front, w3f1_w32f2_front);
      tmp_back[1][j2] = vaddq_s16(f0_back, w3f1_w32f2_back);
      tmp_front[2][j2] = vsubq_s16(vsubq_s16(f0_front, f1_f2_front), w3f1_w32f2_front);
      tmp_back[2][j2] = vsubq_s16(vsubq_s16(f0_back, f1_f2_back), w3f1_w32f2_back);
    }

    tmp_front[1][1] = barret_mul_const<Q, W_9S[1]>(tmp_front[1][1]);
    tmp_back[1][1] = barret_mul_const<Q, W_9S[1]>(tmp_back[1][1]);
    tmp_front[2][1] = barret_mul_const<Q, W_9S[2]>(tmp_front[2][1]);
    tmp_back[2][1] = barret_mul_const<Q, W_9S[2]>(tmp_back[2][1]);

    tmp_front[1][2] = barret_mul_const<Q, W_9S[2]>(tmp_front[1][2]);
    tmp_back[1][2] = barret_mul_const<Q, W_9S[2]>(tmp_back[1][2]);
    tmp_front[2][2] = barret_mul_const<Q, W_9S[4]>(tmp_front[2][2]);
    tmp_back[2][2] = barret_mul_const<Q, W_9S[4]>(tmp_back[2][2]);

    for (int j1 = 0; j1 < 3; j1++) {
      int16x8_t f1_f2_front = vaddq_s16(tmp_front[j1][1], tmp_front[j1][2]);
      int16x8_t f1_f2_back = vaddq_s16(tmp_back[j1][1], tmp_back[j1][2]);
      barret_mla_const<Q, W_3S[1]>(tmp_front[j1][1], tmp_front[j1][2]);
      barret_mla_const<Q, W_3S[1]>(tmp_back[j1][1], tmp_back[j1][2]);
      int16x8_t w3f1_w32f2_front = barret_mul_const<Q, W_3S[1]>(tmp_front[j1][1]);
      int16x8_t w3f1_w32f2_back = barret_mul_const<Q, W_3S[1]>(tmp_back[j1][1]);

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

constexpr std::array<std::array<int16_t, 9>, 10> base_mul_twiddles = [] {
  std::array<std::array<int16_t, 9>, 10> res = {};
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      res[i][j] = center_lift<int64_t, Q>(int64_t(W_10S[i]) * W_9S[j]);
    }
  }
  return res;
} ();

constexpr std::array<std::array<int16_t, 9>, 10> base_mul_twiddle_bars = [] {
  std::array<std::array<int16_t, 9>, 10> res = {};
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      res[i][j] = gen_bar<int16_t, Q>(base_mul_twiddles[i][j]);
    }
  }
  return res;
} ();

template<int LANE>
void base_mul_col(int32x4_t &res_front_low, int32x4_t &res_front_high, int32x4_t &res_back_low, int32x4_t &res_back_high, int16x8_t col_front, int16x8_t col_back, int16x8_t weights) {
  res_front_low = vmlal_laneq_s16(res_front_low, vget_low_s16(col_front), weights, LANE);
  res_front_high = vmlal_high_laneq_s16(res_front_high, col_front, weights, LANE);
  res_back_low = vmlal_laneq_s16(res_back_low, vget_low_s16(col_back), weights, LANE);
  res_back_high = vmlal_high_laneq_s16(res_back_high, col_back, weights, LANE);
}

void base_mul(int16_t in1_ntt[10][9][16], int16_t in2_ntt[10][9][16], int16_t out_ntt[10][9][16]) {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      // std::cerr << "base_mul input, i = " << i << ", j = " << j << ", twiddle = " << base_mul_twiddles[i][j] << '\n';
      // std::cerr << "    ";
      // for (int k = 0; k < 16; k++) {
      //   std::cerr << in1_ntt[i][j][k]  << " \n"[k == 15];
      // }
      // std::cerr << "    ";
      // for (int k = 0; k < 16; k++) {
      //   std::cerr << in2_ntt[i][j][k]  << " \n"[k == 15];
      // }

      int32x4_t res_front_low = {};
      int32x4_t res_front_high = {};
      int32x4_t res_back_low = {};
      int32x4_t res_back_high = {};
      int16_t twiddle = base_mul_twiddles[i][j];
      int16_t twiddle_bar = base_mul_twiddle_bars[i][j];
      int16x8_t in1_front = vld1q_s16(&in1_ntt[i][j][0]);
      int16x8_t in1_back = vld1q_s16(&in1_ntt[i][j][8]);
      int16x8_t in2_front = vld1q_s16(&in2_ntt[i][j][0]);
      int16x8_t in2_back = vld1q_s16(&in2_ntt[i][j][8]);
      int16x8_t twiddle_in2_front = barret_mul<Q>(in2_front, twiddle, twiddle_bar);
      int16x8_t twiddle_in2_back = barret_mul<Q>(in2_back, twiddle, twiddle_bar);

      base_mul_col<0>(res_front_low, res_front_high, res_back_low, res_back_high,
          in2_front, in2_back, in1_front);
      base_mul_col<1>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_back, in2_front, 7), vextq_s16(in2_front, in2_back, 7), in1_front);
      base_mul_col<2>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_back, in2_front, 6), vextq_s16(in2_front, in2_back, 6), in1_front);
      base_mul_col<3>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_back, in2_front, 5), vextq_s16(in2_front, in2_back, 5), in1_front);
      base_mul_col<4>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_back, in2_front, 4), vextq_s16(in2_front, in2_back, 4), in1_front);
      base_mul_col<5>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_back, in2_front, 3), vextq_s16(in2_front, in2_back, 3), in1_front);
      base_mul_col<6>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_back, in2_front, 2), vextq_s16(in2_front, in2_back, 2), in1_front);
      base_mul_col<7>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_back, in2_front, 1), vextq_s16(in2_front, in2_back, 1), in1_front);

      base_mul_col<0>(res_front_low, res_front_high, res_back_low, res_back_high,
          twiddle_in2_back, in2_front, in1_back);
      base_mul_col<1>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_front, twiddle_in2_back, 7), vextq_s16(twiddle_in2_back, in2_front, 7), in1_back);
      base_mul_col<2>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_front, twiddle_in2_back, 6), vextq_s16(twiddle_in2_back, in2_front, 6), in1_back);
      base_mul_col<3>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_front, twiddle_in2_back, 5), vextq_s16(twiddle_in2_back, in2_front, 5), in1_back);
      base_mul_col<4>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_front, twiddle_in2_back, 4), vextq_s16(twiddle_in2_back, in2_front, 4), in1_back);
      base_mul_col<5>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_front, twiddle_in2_back, 3), vextq_s16(twiddle_in2_back, in2_front, 3), in1_back);
      base_mul_col<6>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_front, twiddle_in2_back, 2), vextq_s16(twiddle_in2_back, in2_front, 2), in1_back);
      base_mul_col<7>(res_front_low, res_front_high, res_back_low, res_back_high,
          vextq_s16(twiddle_in2_front, twiddle_in2_back, 1), vextq_s16(twiddle_in2_back, in2_front, 1), in1_back);

      // debug_int32x4(res_front_low);
      // debug_int32x4(res_front_high);
      // debug_int32x4(res_back_low);
      // debug_int32x4(res_back_high);

      int32x4_t esti;
      esti = vqrdmulhq_n_s32(res_front_low, 467759);
      res_front_low = vmlsq_n_s32(res_front_low, esti, Q);
      esti = vqrdmulhq_n_s32(res_front_high, 467759);
      res_front_high = vmlsq_n_s32(res_front_high, esti, Q);
      esti = vqrdmulhq_n_s32(res_back_low, 467759);
      res_back_low = vmlsq_n_s32(res_back_low, esti, Q);
      esti = vqrdmulhq_n_s32(res_back_high, 467759);
      res_back_high = vmlsq_n_s32(res_back_high, esti, Q);

      int16x8_t res_front = vuzp1q_s16(vreinterpretq_s16_s32(res_front_low), vreinterpretq_s16_s32(res_front_high));
      int16x8_t res_back = vuzp1q_s16(vreinterpretq_s16_s32(res_back_low), vreinterpretq_s16_s32(res_back_high));

      vst1q_s16(&out_ntt[i][j][0], res_front);
      vst1q_s16(&out_ntt[i][j][8], res_back);

      // std::cerr << "base_mul output, i = " << i << ", j = " << j << ", twiddle = " << base_mul_twiddles[i][j] << '\n';
      // std::cerr << "    ";
      // for (int k = 0; k < 16; k++) {
      //   std::cerr << out_ntt[i][j][k]  << " \n"[k == 15];
      // }
    }
  }
}

void forward(int16_t in_poly[], int16_t out_ntt[10][9][16]) {
  ntt_10(out_ntt, in_poly);
  ntt_9(out_ntt);
}

constexpr int16_t INV90 = -51;
constexpr int16_t INV90_BAR = gen_bar<int16_t, Q>(INV90);
void div90_main(int16_t main_poly[1440]) {
  int16x8_t zeros = vdupq_n_s16(0);
  for (int i = 0; i < 1440; i += 8) {
    int16x8_t chunk = vld1q_s16(&main_poly[i]);
    chunk = barret_mul_const<Q, INV90>(chunk);
    vst1q_s16(&main_poly[i], chunk);
  }
}

void backward(int16_t in_ntt[10][9][16], int16_t out_main[1440]) {
  intt_10_x10(in_ntt);
  intt_9_x9(in_ntt, out_main);
  div90_main(out_main);
}

void mult_low(int16_t in1_low[81], int16_t in2_low[81], int16_t out_low[81]) {
  std::memset(out_low, 0, sizeof(int16_t[81]));
  for (int i = 0; i < 81; i++) {
    for (int j = 0; i + j < 81; j++) {
      out_low[i + j] = center_lift<int64_t, Q>(out_low[i + j] + int64_t(1) * in1_low[i] * in2_low[j]);
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
  int16x8_t qs = vdupq_n_s16(Q);
  int16x8_t half_qs = vdupq_n_s16((Q - 1) / 2);
  int16x8_t neg_half_qs = vdupq_n_s16(-(Q - 1) / 2);
  for (int i = 0; i < 768; i += 8 * 4) {
    int16x8x4_t chunks = vld1q_s16_x4(&poly[i]);

    barret_reduce<Q>(chunks.val[0]);
    barret_reduce<Q>(chunks.val[1]);
    barret_reduce<Q>(chunks.val[2]);
    barret_reduce<Q>(chunks.val[3]);

    chunks.val[0] = vsubq_s16(chunks.val[0], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[0], half_qs)), qs));
    chunks.val[1] = vsubq_s16(chunks.val[1], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[1], half_qs)), qs));
    chunks.val[2] = vsubq_s16(chunks.val[2], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[2], half_qs)), qs));
    chunks.val[3] = vsubq_s16(chunks.val[3], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[3], half_qs)), qs));

    chunks.val[0] = vaddq_s16(chunks.val[0], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[0], neg_half_qs)), qs));
    chunks.val[1] = vaddq_s16(chunks.val[1], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[1], neg_half_qs)), qs));
    chunks.val[2] = vaddq_s16(chunks.val[2], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[2], neg_half_qs)), qs));
    chunks.val[3] = vaddq_s16(chunks.val[3], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[3], neg_half_qs)), qs));

    vst1q_s16_x4(&poly[i], chunks);
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
  //   std::cerr << out_main[i] << " \n"[i == 1439];
  // }
  crt(out_poly, out_main, out_low);
  center_poly(out_poly);
}
