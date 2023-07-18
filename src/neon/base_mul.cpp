#include "neon/base_mul.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"

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
void base_mul_col_8(int32x4_t &res_low, int32x4_t &res_high, int16x8_t col, int16x8_t weights) {
  res_low = vmlal_laneq_s16(res_low, vget_low_s16(col), weights, LANE);
  res_high = vmlal_high_laneq_s16(res_high, col, weights, LANE);
}

inline void base_mul_8(int32x4_t &res_low, int32x4_t &res_high, int16x8_t m1, int16x8_t m2, int16x8_t v) {
  base_mul_col_8<0>(res_low, res_high, m2, v);
  base_mul_col_8<1>(res_low, res_high, vextq_s16(m1, m2, 7), v);
  base_mul_col_8<2>(res_low, res_high, vextq_s16(m1, m2, 6), v);
  base_mul_col_8<3>(res_low, res_high, vextq_s16(m1, m2, 5), v);
  base_mul_col_8<4>(res_low, res_high, vextq_s16(m1, m2, 4), v);
  base_mul_col_8<5>(res_low, res_high, vextq_s16(m1, m2, 3), v);
  base_mul_col_8<6>(res_low, res_high, vextq_s16(m1, m2, 2), v);
  base_mul_col_8<7>(res_low, res_high, vextq_s16(m1, m2, 1), v);
}

void base_mul(int16_t in1_ntt[9][2][10][8], int16_t in2_ntt[9][2][10][8], int16_t out_ntt[9][2][10][8]) {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      int16_t twiddle = base_mul_twiddles[i][j];
      int16_t twiddle_bar = base_mul_twiddle_bars[i][j];

      int16x8_t a_fr = vld1q_s16(&in1_ntt[j][0][i][0]);
      int16x8_t a_bk = vld1q_s16(&in1_ntt[j][1][i][0]);

      int16x8_t a_fb = a_fr;
      barret_reduce<Q>(a_fb);
      a_fb = vaddq_s16(a_fb, a_bk);

      int16x8_t b_fr = vld1q_s16(&in2_ntt[j][0][i][0]);
      int16x8_t b_bk = vld1q_s16(&in2_ntt[j][1][i][0]);

      barret_reduce<Q>(b_fr);

      int16x8_t c12 = b_fr;
      int16x8_t c11 = barret_mul<Q>(b_bk, twiddle, twiddle_bar);
      int16x8_t c22 = vsubq_s16(b_bk, b_fr);
      int16x8_t c21 = vsubq_s16(b_fr, c11);
      int16x8_t c02 = vnegq_s16(c21);
      int16x8_t c01 = barret_mul<Q>(c22, -twiddle, -twiddle_bar);

      int32x4_t res_fr_low = {};
      int32x4_t res_fr_high = {};

      base_mul_8(res_fr_low, res_fr_high, c11, c12, a_fb);

      int32x4_t res_bk_low = res_fr_low;
      int32x4_t res_bk_high = res_fr_high;

      base_mul_8(res_fr_low, res_fr_high, c01, c02, a_bk);
      base_mul_8(res_bk_low, res_bk_high, c21, c22, a_fr);

      int32x4_t esti;
      esti = vqrdmulhq_n_s32(res_fr_low, 467759);
      res_fr_low = vmlsq_n_s32(res_fr_low, esti, Q);
      esti = vqrdmulhq_n_s32(res_fr_high, 467759);
      res_fr_high = vmlsq_n_s32(res_fr_high, esti, Q);
      esti = vqrdmulhq_n_s32(res_bk_low, 467759);
      res_bk_low = vmlsq_n_s32(res_bk_low, esti, Q);
      esti = vqrdmulhq_n_s32(res_bk_high, 467759);
      res_bk_high = vmlsq_n_s32(res_bk_high, esti, Q);

      int16x8_t res_fr = vuzp1q_s16(vreinterpretq_s16_s32(res_fr_low), vreinterpretq_s16_s32(res_fr_high));
      int16x8_t res_bk = vuzp1q_s16(vreinterpretq_s16_s32(res_bk_low), vreinterpretq_s16_s32(res_bk_high));

      vst1q_s16(&out_ntt[j][0][i][0], res_fr);
      vst1q_s16(&out_ntt[j][1][i][0], res_bk);
    }
  }
}
