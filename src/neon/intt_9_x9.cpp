#include "neon/intt_9_x9.h"

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

constexpr std::array<int16_t, 17> W_9S = gen_pows<int16_t, 17, Q>(W_9);
constexpr std::array<int16_t, 5> W_3S = gen_pows<int16_t, 5, Q>(W_3);

constexpr static std::array<int16_t, 8> COEFS = {
  W_3S[1], W_3S[2], W_9S[1], W_9S[2], W_9S[4]
};
constexpr static std::array<int16_t, 8> COEFS_MOD = back_mod<Q>(COEFS);
constexpr static std::array<int16_t, 8> BARS = gen_bars<int16_t, 8, Q>(COEFS);
constexpr static std::array<int16_t, 8> BARS_RED = back_red<Q>(BARS);


inline void btrfly3(
    int16x8_t x0, int16x8_t x1, int16x8_t x2,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2) {

  int16x8_t a02 = vaddq_s16(x0, x2);
  int16x8_t s02 = vsubq_s16(x0, x2);
  int16x8_t s12 = vsubq_s16(x1, x2);

  int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
  int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

  h0 = vaddq_s16(a02, x1);
  h1 = barret_mul_laneq_opaque<Q, 0>(s12, coefs_mod, bars_red, coefs_mod);
  h1 = vaddq_s16(h1, s02);
  h2 = barret_mul_laneq_opaque<Q, 1>(s12, coefs_mod, bars_red, coefs_mod);
  h2 = vaddq_s16(h2, s02);
}

inline void btrfly9(
    int16x8_t x0, int16x8_t x1, int16x8_t x2, int16x8_t x3, int16x8_t x4,
    int16x8_t x5, int16x8_t x6, int16x8_t x7, int16x8_t x8,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2, int16x8_t &h3, int16x8_t &h4,
    int16x8_t &h5, int16x8_t &h6, int16x8_t &h7, int16x8_t &h8) {

  int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
  int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

  int16x8_t a0, a1, a2;
  btrfly3(x0, x3, x6, a0, a1, a2);

  int16x8_t b0, b1, b2;
  btrfly3(x1, x4, x7, b0, b1, b2);

  int16x8_t c0, c1, c2;
  btrfly3(x2, x5, x8, c0, c1, c2);

  barret_crude_reduce_laneq_opaque<Q>(a0, bars_red, coefs_mod);
  barret_crude_reduce_laneq_opaque<Q>(a1, bars_red, coefs_mod);
  barret_crude_reduce_laneq_opaque<Q>(a2, bars_red, coefs_mod);
  barret_crude_reduce_laneq_opaque<Q>(b0, bars_red, coefs_mod);
  b1 = barret_mul_laneq_opaque<Q, 2>(b1, coefs_mod, bars_red, coefs_mod);
  b2 = barret_mul_laneq_opaque<Q, 3>(b2, coefs_mod, bars_red, coefs_mod);
  barret_crude_reduce_laneq_opaque<Q>(c0, bars_red, coefs_mod);
  c1 = barret_mul_laneq_opaque<Q, 3>(c1, coefs_mod, bars_red, coefs_mod);
  c2 = barret_mul_laneq_opaque<Q, 4>(c2, coefs_mod, bars_red, coefs_mod);

  btrfly3(a0, b0, c0, h0, h3, h6);
  btrfly3(a1, b1, c1, h1, h4, h7);
  btrfly3(a2, b2, c2, h2, h5, h8);
}


void intt_9_x9(int16_t ntt[9][2][10][8], int16_t poly[1440]) {

  int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
  int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

  for (int i = 0; i < 10; i++) {
    {
      int16x8_t x0_fr = vld1q_s16(&ntt[0][0][i][0]);
      int16x8_t x1_fr = vld1q_s16(&ntt[8][0][i][0]);
      int16x8_t x2_fr = vld1q_s16(&ntt[7][0][i][0]);
      int16x8_t x3_fr = vld1q_s16(&ntt[6][0][i][0]);
      int16x8_t x4_fr = vld1q_s16(&ntt[5][0][i][0]);
      int16x8_t x5_fr = vld1q_s16(&ntt[4][0][i][0]);
      int16x8_t x6_fr = vld1q_s16(&ntt[3][0][i][0]);
      int16x8_t x7_fr = vld1q_s16(&ntt[2][0][i][0]);
      int16x8_t x8_fr = vld1q_s16(&ntt[1][0][i][0]);

      barret_reduce_laneq_opaque<Q>(x0_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x1_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x2_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x3_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x4_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x5_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x6_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x7_fr, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x8_fr, bars_red, coefs_mod);

      int16x8_t h0_fr, h1_fr, h2_fr, h3_fr, h4_fr, h5_fr, h6_fr, h7_fr, h8_fr;

      btrfly9(x0_fr, x1_fr, x2_fr, x3_fr, x4_fr,
          x5_fr, x6_fr, x7_fr, x8_fr,
          h0_fr, h1_fr, h2_fr, h3_fr, h4_fr,
          h5_fr, h6_fr, h7_fr, h8_fr);

      vst1q_s16(&poly[81 * i % 90 * 16], h0_fr);
      vst1q_s16(&poly[(81 * i + 10) % 90 * 16], h1_fr);
      vst1q_s16(&poly[(81 * i + 20) % 90 * 16], h2_fr);
      vst1q_s16(&poly[(81 * i + 30) % 90 * 16], h3_fr);
      vst1q_s16(&poly[(81 * i + 40) % 90 * 16], h4_fr);
      vst1q_s16(&poly[(81 * i + 50) % 90 * 16], h5_fr);
      vst1q_s16(&poly[(81 * i + 60) % 90 * 16], h6_fr);
      vst1q_s16(&poly[(81 * i + 70) % 90 * 16], h7_fr);
      vst1q_s16(&poly[(81 * i + 80) % 90 * 16], h8_fr);
    }

    {
      int16x8_t x0_bk = vld1q_s16(&ntt[0][1][i][0]);
      int16x8_t x1_bk = vld1q_s16(&ntt[8][1][i][0]);
      int16x8_t x2_bk = vld1q_s16(&ntt[7][1][i][0]);
      int16x8_t x3_bk = vld1q_s16(&ntt[6][1][i][0]);
      int16x8_t x4_bk = vld1q_s16(&ntt[5][1][i][0]);
      int16x8_t x5_bk = vld1q_s16(&ntt[4][1][i][0]);
      int16x8_t x6_bk = vld1q_s16(&ntt[3][1][i][0]);
      int16x8_t x7_bk = vld1q_s16(&ntt[2][1][i][0]);
      int16x8_t x8_bk = vld1q_s16(&ntt[1][1][i][0]);

      barret_reduce_laneq_opaque<Q>(x0_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x1_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x2_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x3_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x4_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x5_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x6_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x7_bk, bars_red, coefs_mod);
      barret_reduce_laneq_opaque<Q>(x8_bk, bars_red, coefs_mod);

      int16x8_t h0_bk, h1_bk, h2_bk, h3_bk, h4_bk, h5_bk, h6_bk, h7_bk, h8_bk;

      btrfly9(x0_bk, x1_bk, x2_bk, x3_bk, x4_bk,
          x5_bk, x6_bk, x7_bk, x8_bk,
          h0_bk, h1_bk, h2_bk, h3_bk, h4_bk,
          h5_bk, h6_bk, h7_bk, h8_bk);

      vst1q_s16(&poly[81 * i % 90 * 16 + 8], h0_bk);
      vst1q_s16(&poly[(81 * i + 10) % 90 * 16 + 8], h1_bk);
      vst1q_s16(&poly[(81 * i + 20) % 90 * 16 + 8], h2_bk);
      vst1q_s16(&poly[(81 * i + 30) % 90 * 16 + 8], h3_bk);
      vst1q_s16(&poly[(81 * i + 40) % 90 * 16 + 8], h4_bk);
      vst1q_s16(&poly[(81 * i + 50) % 90 * 16 + 8], h5_bk);
      vst1q_s16(&poly[(81 * i + 60) % 90 * 16 + 8], h6_bk);
      vst1q_s16(&poly[(81 * i + 70) % 90 * 16 + 8], h7_bk);
      vst1q_s16(&poly[(81 * i + 80) % 90 * 16 + 8], h8_bk);
    }
  }

}
