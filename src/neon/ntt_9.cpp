#include "neon/ntt_9.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"
#include "arith_tmpl/neon_arith_opaque.h"
#include "arith_tmpl/reduce_tbl.h"

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;
constexpr int16_t W_9 = gen_pow<int16_t, Q>(W_4590, ORD / 9);
constexpr int16_t W_3 = gen_pow<int16_t, Q>(W_4590, ORD / 3);

constexpr static std::array<int16_t, 9> W_9S = gen_pows<int16_t, 9, Q>(W_9);
constexpr static std::array<int16_t, 3> W_3S = gen_pows<int16_t, 3, Q>(W_3);

constexpr static std::array<int16_t, 9> W_9_BARS = gen_bars<int16_t, 9, Q>(W_9S);

// constexpr static std::array<int16_t, 8> COEFS_ARR = {
//   W_3S[1], W_3S[2], W_9S[1], W_9S[2], W_9S[4]
// };
// constexpr static std::array<int16_t, 8> COEFS_MOD_ARR = back_mod<Q>(COEFS_ARR);
// constexpr static int16x8_t COEFS_MOD = arr_to_x8_t(COEFS_MOD_ARR);
constexpr static std::array<int16_t, 8> COEFS = {
  W_3S[1], W_3S[2], W_9S[1], W_9S[2], W_9S[4]
};
constexpr static std::array<int16_t, 8> COEFS_MOD = back_mod<Q>(COEFS);
constexpr static std::array<int16_t, 8> BARS = gen_bars<int16_t, 8, Q>(COEFS);
constexpr static std::array<int16_t, 8> BARS_RED = back_red<Q>(BARS);

constexpr static
std::pair<std::array<std::array<int16_t, 8>, 10>, std::array<std::array<int16_t, 8>, 10>>
TWISTS = [] {
  std::pair<std::array<std::array<int16_t, 8>, 10>, std::array<std::array<int16_t, 8>, 10>> res = {};
  for (int i = 0; i < 10; i++) {
    for (int j = 1; j < 9; j++) {
      res.first[i][j - 1] = W_9S[i * j % 9];
      res.second[i][j - 1] = W_9_BARS[i * j % 9];
    }
  }
  return res;
} ();

inline void btrfly3(
    int16x8_t x0, int16x8_t x1, int16x8_t x2,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2) {

  int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
  int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

  int16x8_t a02 = vaddq_s16(x0, x2);
  int16x8_t s02 = vsubq_s16(x0, x2);
  int16x8_t s12 = vsubq_s16(x1, x2);

  h0 = vaddq_s16(a02, x1);
  h1 = barret_mul_laneq_opaque<Q, 0>(s12, coefs_mod, bars_red, coefs_mod);
  h1 = vaddq_s16(h1, s02);
  h2 = barret_mul_laneq_opaque<Q, 1>(s12, coefs_mod, bars_red, coefs_mod);
  h2 = vaddq_s16(h2, s02);
}

inline void five_nonezero(
    int16x8_t x0, int16x8_t x1, int16x8_t x2, int16x8_t x3, int16x8_t x4,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2, int16x8_t &h3, int16x8_t &h4,
    int16x8_t &h5, int16x8_t &h6, int16x8_t &h7, int16x8_t &h8) {

  int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
  int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

  int16x8_t a0 = vaddq_s16(x0, x3);
  int16x8_t a1 = barret_mul_laneq_opaque<Q, 0>(x3, coefs_mod, bars_red, coefs_mod);
  int16x8_t a2 = vsubq_s16(x0, x3);
  a2 = vsubq_s16(a2, a1);
  a1 = vaddq_s16(a1, x0);

  int16x8_t b0 = vaddq_s16(x1, x4);
  int16x8_t b1 = barret_mul_laneq_opaque<Q, 0>(x4, coefs_mod, bars_red, coefs_mod);
  int16x8_t b2 = vsubq_s16(x1, x4);
  b2 = vsubq_s16(b2, b1);
  b1 = vaddq_s16(b1, x1);

  b1 = barret_mul_laneq_opaque<Q, 2>(b1, coefs_mod, bars_red, coefs_mod);
  b2 = barret_mul_laneq_opaque<Q, 3>(b2, coefs_mod, bars_red, coefs_mod);

  int16x8_t c0 = x2;
  int16x8_t c1 = barret_mul_laneq_opaque<Q, 3>(x2, coefs_mod, bars_red, coefs_mod);
  int16x8_t c2 = barret_mul_laneq_opaque<Q, 4>(x2, coefs_mod, bars_red, coefs_mod);

  btrfly3(a0, b0, c0, h0, h3, h6);
  btrfly3(a1, b1, c1, h1, h4, h7);
  btrfly3(a2, b2, c2, h2, h5, h8);
}

void ntt_9(int16_t ntt[9][2][10][8], const int16_t poly[800]) {
  for (int i = 0; i < 10; i++) {
    int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
    int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

    int16x8_t twist_coefs = vld1q_s16(&TWISTS.first[i][0]);
    int16x8_t twist_bars = vld1q_s16(&TWISTS.second[i][0]);

    int16x8x2_t x0 = vld1q_s16_x2(&poly[i * 16]);
    int16x8x2_t x1 = vld1q_s16_x2(&poly[(10 + i) * 16]);
    // int16x8x2_t x2 = vld1q_s16_x2(&poly[(20 + i) * 16]);
    int16x8x2_t x3 = vld1q_s16_x2(&poly[(30 + i) * 16]);
    int16x8x2_t x4 = vld1q_s16_x2(&poly[(40 + i) * 16]);

    {
      int16x8_t x0_fr = x0.val[0];
      int16x8_t x1_fr = x1.val[0];
      // int16x8_t x2_fr = x2.val[0];
      int16x8_t x2_fr = vld1q_s16(&poly[(20 + i) * 16]);
      int16x8_t x3_fr = x3.val[0];
      int16x8_t x4_fr = x4.val[0];

      int16x8_t h0_fr, h1_fr, h2_fr, h3_fr, h4_fr, h5_fr, h6_fr, h7_fr, h8_fr;

      five_nonezero(x0_fr, x1_fr, x2_fr, x3_fr, x4_fr,
          h0_fr, h1_fr, h2_fr, h3_fr, h4_fr,
          h5_fr, h6_fr, h7_fr, h8_fr);

      barret_reduce_laneq_opaque<Q>(h0_fr, bars_red, coefs_mod);
      vst1q_s16(&ntt[0][0][i][0], h0_fr);

      h1_fr = barret_mul_laneq_opaque<Q, 0>(h1_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[1][0][i][0], h1_fr);

      h2_fr = barret_mul_laneq_opaque<Q, 1>(h2_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[2][0][i][0], h2_fr);

      h3_fr = barret_mul_laneq_opaque<Q, 2>(h3_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[3][0][i][0], h3_fr);

      h4_fr = barret_mul_laneq_opaque<Q, 3>(h4_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[4][0][i][0], h4_fr);

      h5_fr = barret_mul_laneq_opaque<Q, 4>(h5_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[5][0][i][0], h5_fr);

      h6_fr = barret_mul_laneq_opaque<Q, 5>(h6_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[6][0][i][0], h6_fr);

      h7_fr = barret_mul_laneq_opaque<Q, 6>(h7_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[7][0][i][0], h7_fr);

      h8_fr = barret_mul_laneq_opaque<Q, 7>(h8_fr, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[8][0][i][0], h8_fr);
    }

    {
      int16x8_t x0_bk = x0.val[1];
      int16x8_t x1_bk = x1.val[1];
      // int16x8_t x2_bk = x2.val[1];
      int16x8_t x2_bk = vld1q_s16(&poly[(20 + i) * 16 + 8]);
      int16x8_t x3_bk = x3.val[1];
      int16x8_t x4_bk = x4.val[1];

      int16x8_t h0_bk, h1_bk, h2_bk, h3_bk, h4_bk, h5_bk, h6_bk, h7_bk, h8_bk;

      five_nonezero(x0_bk, x1_bk, x2_bk, x3_bk, x4_bk,
          h0_bk, h1_bk, h2_bk, h3_bk, h4_bk,
          h5_bk, h6_bk, h7_bk, h8_bk);

      barret_reduce_laneq_opaque<Q>(h0_bk, bars_red, coefs_mod);
      vst1q_s16(&ntt[0][1][i][0], h0_bk);

      h1_bk = barret_mul_laneq_opaque<Q, 0>(h1_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[1][1][i][0], h1_bk);

      h2_bk = barret_mul_laneq_opaque<Q, 1>(h2_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[2][1][i][0], h2_bk);

      h3_bk = barret_mul_laneq_opaque<Q, 2>(h3_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[3][1][i][0], h3_bk);

      h4_bk = barret_mul_laneq_opaque<Q, 3>(h4_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[4][1][i][0], h4_bk);

      h5_bk = barret_mul_laneq_opaque<Q, 4>(h5_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[5][1][i][0], h5_bk);

      h6_bk = barret_mul_laneq_opaque<Q, 5>(h6_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[6][1][i][0], h6_bk);

      h7_bk = barret_mul_laneq_opaque<Q, 6>(h7_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[7][1][i][0], h7_bk);

      h8_bk = barret_mul_laneq_opaque<Q, 7>(h8_bk, twist_coefs, twist_bars, coefs_mod);
      vst1q_s16(&ntt[8][1][i][0], h8_bk);
    }
  }
}
