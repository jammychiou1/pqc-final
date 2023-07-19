#include "neon/base_mul.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"

// #include <iostream>
// #include "utils/debug.h"

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;

constexpr int16_t W_10 = gen_pow<int16_t, Q>(W_4590, ORD / 10);
constexpr int16_t W_9 = gen_pow<int16_t, Q>(W_4590, ORD / 9);

constexpr std::array<int16_t, 10> W_10S = gen_pows<int16_t, 10, Q>(W_10);
constexpr std::array<int16_t, 17> W_9S = gen_pows<int16_t, 17, Q>(W_9);

constexpr static std::array<std::array<int32_t, 9>, 10> TWIDDLES = [] {
  std::array<std::array<int32_t, 9>, 10> res = {};
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      int16_t coef = center_lift<int64_t, Q>(int64_t(W_10S[i]) * W_9S[j]);
      int16_t bar = gen_bar<int16_t, Q>(coef);
      // note: `a << b` for signed and negative `a` is undefined until c++20
      res[i][j] = uint16_t(bar) << 16;
      res[i][j] |= uint16_t(coef);
    }
  }
  return res;
} ();

template<int LANE>
void mla_v_sc_8(int32x4_t &res_low, int32x4_t &res_high, int16x8_t vec, int16x8_t scl) {
  res_low = vmlal_laneq_s16(res_low, vget_low_s16(vec), scl, LANE);
  res_high = vmlal_high_laneq_s16(res_high, vec, scl, LANE);
}

inline void mla_tm_v_8(int32x4_t &res_low, int32x4_t &res_high, int16x8_t tm1, int16x8_t tm2, int16x8_t vec) {
  mla_v_sc_8<0>(res_low, res_high, tm2, vec);
  mla_v_sc_8<1>(res_low, res_high, vextq_s16(tm1, tm2, 7), vec);
  mla_v_sc_8<2>(res_low, res_high, vextq_s16(tm1, tm2, 6), vec);
  mla_v_sc_8<3>(res_low, res_high, vextq_s16(tm1, tm2, 5), vec);
  mla_v_sc_8<4>(res_low, res_high, vextq_s16(tm1, tm2, 4), vec);
  mla_v_sc_8<5>(res_low, res_high, vextq_s16(tm1, tm2, 3), vec);
  mla_v_sc_8<6>(res_low, res_high, vextq_s16(tm1, tm2, 2), vec);
  mla_v_sc_8<7>(res_low, res_high, vextq_s16(tm1, tm2, 1), vec);
}

template<int LANE>
void mls_v_sc_8(int32x4_t &res_low, int32x4_t &res_high, int16x8_t vec, int16x8_t scl) {
  res_low = vmlsl_laneq_s16(res_low, vget_low_s16(vec), scl, LANE);
  res_high = vmlsl_high_laneq_s16(res_high, vec, scl, LANE);
}

inline void mls_tm_v_8(int32x4_t &res_low, int32x4_t &res_high, int16x8_t tm1, int16x8_t tm2, int16x8_t vec) {
  mls_v_sc_8<0>(res_low, res_high, tm2, vec);
  mls_v_sc_8<1>(res_low, res_high, vextq_s16(tm1, tm2, 7), vec);
  mls_v_sc_8<2>(res_low, res_high, vextq_s16(tm1, tm2, 6), vec);
  mls_v_sc_8<3>(res_low, res_high, vextq_s16(tm1, tm2, 5), vec);
  mls_v_sc_8<4>(res_low, res_high, vextq_s16(tm1, tm2, 4), vec);
  mls_v_sc_8<5>(res_low, res_high, vextq_s16(tm1, tm2, 3), vec);
  mls_v_sc_8<6>(res_low, res_high, vextq_s16(tm1, tm2, 2), vec);
  mls_v_sc_8<7>(res_low, res_high, vextq_s16(tm1, tm2, 1), vec);
}

constexpr static int32x4_t COEF = [] {
  int16_t one_bar = gen_bar<int16_t, Q>(1);
  // note: `a << b` for signed and negative `a` is undefined until c++20
  int32_t comb = uint16_t(one_bar) << 16;
  comb |= uint16_t(Q);
  int32x4_t coefs = {0, comb, Q, 467759};
  return coefs;
} ();

inline void reduce_comb(int16x8_t &vd, int32x4_t comb_vec) {
  int16x4_t low = vget_low_s16(vreinterpretq_s16_s32(comb_vec));
  int16x8_t esti = vqrdmulhq_lane_s16(vd, low, 3);
  vd = vmlsq_lane_s16(vd, esti, low, 2);
}

inline int16x8_t twist_comb(int16x8_t v, int32x4_t comb_vec) {
  int16x4_t low = vget_low_s16(vreinterpretq_s16_s32(comb_vec));
  int16x8_t esti = vqrdmulhq_lane_s16(v, low, 1);
  int16x8_t res = vmulq_lane_s16(v, low, 0);
  res = vmlsq_lane_s16(res, esti, low, 2);
  return res;
}

void base_mul(int16_t in1_ntt[9][2][10][8], int16_t in2_ntt[9][2][10][8], int16_t out_ntt[9][2][10][8]) {
  int32x4_t coefs = COEF;
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      int32_t comb = TWIDDLES[i][j];
      coefs = vsetq_lane_s32(comb, coefs, 0);

      // std::cerr << "i = " << i << ", j = " << j << '\n';
      // int16_t ref_coef = center_lift<int64_t, Q>(int64_t(W_10S[i]) * W_9S[j]);
      // int16_t ref_bar = gen_bar<int16_t, Q>(ref_coef);
      // std::cerr << "  ref_coef = " << ref_coef << ", ref_bar = " << ref_bar << '\n';
      // debug_int32x4(coefs);
      // debug_int16x8(vreinterpretq_s16_s32(coefs));

      int16x8_t a_fr = vld1q_s16(&in1_ntt[j][0][i][0]);
      int16x8_t a_bk = vld1q_s16(&in1_ntt[j][1][i][0]);

      int16x8_t a_fb = a_fr;
      reduce_comb(a_fb, coefs);
      a_fb = vaddq_s16(a_fb, a_bk);

      int16x8_t b_fr = vld1q_s16(&in2_ntt[j][0][i][0]);
      int16x8_t b_bk = vld1q_s16(&in2_ntt[j][1][i][0]);

      reduce_comb(b_fr, coefs);

      int16x8_t c12 = b_fr;
      int16x8_t c11 = twist_comb(b_bk, coefs);
      int16x8_t c22 = vsubq_s16(b_bk, b_fr);
      int16x8_t c21 = vsubq_s16(b_fr, c11);
      int16x8_t c02 = c21;
      int16x8_t c01 = twist_comb(c22, coefs);

      int32x4_t res_fr_low = {};
      int32x4_t res_fr_high = {};

      mla_tm_v_8(res_fr_low, res_fr_high, c11, c12, a_fb);

      int32x4_t res_bk_low = res_fr_low;
      int32x4_t res_bk_high = res_fr_high;

      mls_tm_v_8(res_fr_low, res_fr_high, c01, c02, a_bk);
      mla_tm_v_8(res_bk_low, res_bk_high, c21, c22, a_fr);

      int32x4_t esti;
      esti = vqrdmulhq_laneq_s32(res_fr_low, coefs, 3);
      res_fr_low = vmlsq_laneq_s32(res_fr_low, esti, coefs, 2);
      esti = vqrdmulhq_laneq_s32(res_fr_high, coefs, 3);
      res_fr_high = vmlsq_laneq_s32(res_fr_high, esti, coefs, 2);
      esti = vqrdmulhq_laneq_s32(res_bk_low, coefs, 3);
      res_bk_low = vmlsq_laneq_s32(res_bk_low, esti, coefs, 2);
      esti = vqrdmulhq_laneq_s32(res_bk_high, coefs, 3);
      res_bk_high = vmlsq_laneq_s32(res_bk_high, esti, coefs, 2);

      int16x8_t res_fr = vuzp1q_s16(vreinterpretq_s16_s32(res_fr_low), vreinterpretq_s16_s32(res_fr_high));
      int16x8_t res_bk = vuzp1q_s16(vreinterpretq_s16_s32(res_bk_low), vreinterpretq_s16_s32(res_bk_high));

      vst1q_s16(&out_ntt[j][0][i][0], res_fr);
      vst1q_s16(&out_ntt[j][1][i][0], res_bk);
    }
  }
}
