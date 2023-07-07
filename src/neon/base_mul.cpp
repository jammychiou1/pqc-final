#include "neon/base_mul.h"

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

constexpr int16_t W_10 = gen_pow<int16_t, Q>(W_4590, ORD / 10);
constexpr int16_t W_9 = gen_pow<int16_t, Q>(W_4590, ORD / 9);

constexpr std::array<int16_t, 10> W_10S = gen_pows<int16_t, 10, Q>(W_10);
constexpr std::array<int16_t, 17> W_9S = gen_pows<int16_t, 17, Q>(W_9);

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
