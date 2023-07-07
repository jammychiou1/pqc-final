#include "neon/intt_9_x9.h"

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
constexpr int16_t W_9 = gen_pow<int16_t, Q>(W_4590, ORD / 9);
constexpr int16_t W_3 = gen_pow<int16_t, Q>(W_4590, ORD / 3);

constexpr std::array<int16_t, 17> W_9S = gen_pows<int16_t, 17, Q>(W_9);
constexpr std::array<int16_t, 5> W_3S = gen_pows<int16_t, 5, Q>(W_3);

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
