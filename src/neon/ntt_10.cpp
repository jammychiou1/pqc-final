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
constexpr std::array<int16_t, 17> W_5_BARS = gen_bars<int16_t, 17, Q>(W_5S);

void ntt_10(int16_t ntt[10][9][16], const int16_t poly[1440]) {

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
