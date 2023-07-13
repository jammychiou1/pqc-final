#include "neon/ntt_10.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>
#include <iostream>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"
#include "utils/debug.h"

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;

constexpr int16_t W_5 = gen_pow<int16_t, Q>(W_4590, ORD / 5);
constexpr std::array<int16_t, 17> W_5S = gen_pows<int16_t, 17, Q>(W_5);

constexpr int16_t W5_W5_4 = 502;
constexpr int16_t W5_mW5_4 = 459;
constexpr int16_t W5_2_W5_3 = -503;
constexpr int16_t W5_2_mW5_3 = 868;
constexpr int16_t W5_W5_2_mW5_3_mW5_4 = 1327;

constexpr int16_t INV2 = -2295;

void ntt_10(int16_t ntt[10][9][16]) {

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

      for (int i2 = 0; i2 < 2; i2++) {
        int16x8_t f0_front = vld1q_s16(&ntt[5 * i2][j][0]);
        int16x8_t f1_front = vld1q_s16(&ntt[(5 * i2 + 6) % 10][j][0]);
        int16x8_t f2_front = vld1q_s16(&ntt[(5 * i2 + 2) % 10][j][0]);
        int16x8_t f3_front = vld1q_s16(&ntt[(5 * i2 + 8) % 10][j][0]);
        int16x8_t f4_front = vld1q_s16(&ntt[(5 * i2 + 4) % 10][j][0]);

        int16x8_t f1_f4_front = vaddq_s16(f1_front, f4_front);
        int16x8_t f1_mf4_front = vsubq_s16(f1_front, f4_front);

        int16x8_t f3_f2_front = vaddq_s16(f3_front, f2_front);
        int16x8_t f3_mf2_front = vsubq_s16(f3_front, f2_front);

        int16x8_t f1_f2_f3_f4_front = vaddq_s16(f1_f4_front, f3_f2_front);
        int16x8_t f1_mf2_f3_mf4_front = vaddq_s16(f1_mf4_front, f3_mf2_front);

        tmp_front[0][i2] = vaddq_s16(f0_front, f1_f2_f3_f4_front);
        tmp_front[0][i2] = barret_mul_const<Q, -2>(tmp_front[0][i2]);

        int16x8_t neg_c0 = barret_mul_const<Q, -W5_W5_4>(f1_f4_front);
        barret_mla_const<Q, -W5_2_W5_3>(neg_c0, f3_f2_front);
        int16x8_t neg_c1 = vsubq_s16(f1_f2_f3_f4_front, neg_c0);
        barret_reduce<Q>(neg_c1);

        f1_mf4_front = barret_mul_const<Q, W5_mW5_4>(f1_mf4_front);
        f3_mf2_front = barret_mul_const<Q, W5_2_mW5_3>(f3_mf2_front);
        int16x8_t neg_n0 = vsubq_s16(f3_mf2_front, f1_mf4_front);
        int16x8_t neg_n1 = vaddq_s16(f3_mf2_front, f1_mf4_front);
        barret_mla_const<Q, -W5_W5_2_mW5_3_mW5_4>(neg_n1, f1_mf2_f3_mf4_front);

        tmp_front[1][i2] = vaddq_s16(neg_c0, neg_n0);
        tmp_front[2][i2] = vaddq_s16(neg_c1, neg_n1);
        tmp_front[4][i2] = vsubq_s16(neg_c0, neg_n0);
        tmp_front[3][i2] = vsubq_s16(neg_c1, neg_n1);

        int16x8_t df0_front = vshlq_n_s16(f0_front, 1);

        tmp_front[1][i2] = vsubq_s16(tmp_front[1][i2], df0_front);
        tmp_front[2][i2] = vsubq_s16(tmp_front[2][i2], df0_front);
        tmp_front[3][i2] = vsubq_s16(tmp_front[3][i2], df0_front);
        tmp_front[4][i2] = vsubq_s16(tmp_front[4][i2], df0_front);

        // tmp_front[1][i2] = barret_mul_const<Q, -INV2>(tmp_front[1][i2]);
        // tmp_front[2][i2] = barret_mul_const<Q, -INV2>(tmp_front[2][i2]);
        // tmp_front[3][i2] = barret_mul_const<Q, -INV2>(tmp_front[3][i2]);
        // tmp_front[4][i2] = barret_mul_const<Q, -INV2>(tmp_front[4][i2]);

        // barret_reduce<Q>(tmp_front[1][i2]);
        // barret_reduce<Q>(tmp_front[2][i2]);
        // barret_reduce<Q>(tmp_front[3][i2]);
        // barret_reduce<Q>(tmp_front[4][i2]);
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
        int16x8_t f1_back = vld1q_s16(&ntt[(5 * i2 + 6) % 10][j][8]);
        int16x8_t f2_back = vld1q_s16(&ntt[(5 * i2 + 2) % 10][j][8]);
        int16x8_t f3_back = vld1q_s16(&ntt[(5 * i2 + 8) % 10][j][8]);
        int16x8_t f4_back = vld1q_s16(&ntt[(5 * i2 + 4) % 10][j][8]);

        int16x8_t f1_f4_back = vaddq_s16(f1_back, f4_back);
        int16x8_t f1_mf4_back = vsubq_s16(f1_back, f4_back);

        int16x8_t f3_f2_back = vaddq_s16(f3_back, f2_back);
        int16x8_t f3_mf2_back = vsubq_s16(f3_back, f2_back);

        int16x8_t f1_f2_f3_f4_back = vaddq_s16(f1_f4_back, f3_f2_back);
        int16x8_t f1_mf2_f3_mf4_back = vaddq_s16(f1_mf4_back, f3_mf2_back);

        tmp_back[0][i2] = vaddq_s16(f0_back, f1_f2_f3_f4_back);
        tmp_back[0][i2] = barret_mul_const<Q, -2>(tmp_back[0][i2]);

        int16x8_t neg_c0 = barret_mul_const<Q, -W5_W5_4>(f1_f4_back);
        barret_mla_const<Q, -W5_2_W5_3>(neg_c0, f3_f2_back);
        int16x8_t neg_c1 = vsubq_s16(f1_f2_f3_f4_back, neg_c0);
        barret_reduce<Q>(neg_c1);

        f1_mf4_back = barret_mul_const<Q, W5_mW5_4>(f1_mf4_back);
        f3_mf2_back = barret_mul_const<Q, W5_2_mW5_3>(f3_mf2_back);
        int16x8_t neg_n0 = vsubq_s16(f3_mf2_back, f1_mf4_back);
        int16x8_t neg_n1 = vaddq_s16(f3_mf2_back, f1_mf4_back);
        barret_mla_const<Q, -W5_W5_2_mW5_3_mW5_4>(neg_n1, f1_mf2_f3_mf4_back);

        tmp_back[1][i2] = vaddq_s16(neg_c0, neg_n0);
        tmp_back[2][i2] = vaddq_s16(neg_c1, neg_n1);
        tmp_back[4][i2] = vsubq_s16(neg_c0, neg_n0);
        tmp_back[3][i2] = vsubq_s16(neg_c1, neg_n1);

        int16x8_t df0_back = vshlq_n_s16(f0_back, 1);

        tmp_back[1][i2] = vsubq_s16(tmp_back[1][i2], df0_back);
        tmp_back[2][i2] = vsubq_s16(tmp_back[2][i2], df0_back);
        tmp_back[3][i2] = vsubq_s16(tmp_back[3][i2], df0_back);
        tmp_back[4][i2] = vsubq_s16(tmp_back[4][i2], df0_back);

        // tmp_back[1][i2] = barret_mul_const<Q, -INV2>(tmp_back[1][i2]);
        // tmp_back[2][i2] = barret_mul_const<Q, -INV2>(tmp_back[2][i2]);
        // tmp_back[3][i2] = barret_mul_const<Q, -INV2>(tmp_back[3][i2]);
        // tmp_back[4][i2] = barret_mul_const<Q, -INV2>(tmp_back[4][i2]);

        // barret_reduce<Q>(tmp_back[1][i2]);
        // barret_reduce<Q>(tmp_back[2][i2]);
        // barret_reduce<Q>(tmp_back[3][i2]);
        // barret_reduce<Q>(tmp_back[4][i2]);
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
