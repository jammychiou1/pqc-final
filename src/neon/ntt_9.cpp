#include "neon/ntt_9.h"

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

constexpr std::array<int16_t, 9> W_9S = gen_pows<int16_t, 9, Q>(W_9);
constexpr std::array<int16_t, 3> W_3S = gen_pows<int16_t, 3, Q>(W_3);

constexpr std::array<int16_t, 9> W_9_BARS = gen_bars<int16_t, 9, Q>(W_9S);

constexpr std::array<std::array<std::pair<int16_t, int16_t>, 9>, 10> TWISTS = [] {
  std::array<std::array<std::pair<int16_t, int16_t>, 9>, 10> res = {};
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      res[i][j].first = W_9S[i * j % 9];
      res[i][j].second = W_9_BARS[i * j % 9];
    }
  }
  return res;
} ();

inline void btrfly3(
    int16x8_t x0, int16x8_t x1, int16x8_t x2,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2) {

  int16x8_t a02 = vaddq_s16(x0, x2);
  int16x8_t s02 = vsubq_s16(x0, x2);
  int16x8_t s12 = vsubq_s16(x1, x2);

  h0 = vaddq_s16(a02, x1);
  h1 = barret_mul_const<Q, W_3S[1]>(s12);
  h1 = vaddq_s16(h1, s02);
  h2 = barret_mul_const<Q, W_3S[2]>(s12);
  h2 = vaddq_s16(h2, s02);
}

inline void five_nonezero(
    int16x8_t x0, int16x8_t x1, int16x8_t x2, int16x8_t x3, int16x8_t x4,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2, int16x8_t &h3, int16x8_t &h4,
    int16x8_t &h5, int16x8_t &h6, int16x8_t &h7, int16x8_t &h8) {

  int16x8_t a0 = vaddq_s16(x0, x3);
  int16x8_t a1 = barret_mul_const<Q, W_3S[1]>(x3);
  int16x8_t a2 = vsubq_s16(x0, x3);
  a2 = vsubq_s16(a2, a1);
  a1 = vaddq_s16(a1, x0);

  int16x8_t b0 = vaddq_s16(x1, x4);
  int16x8_t b1 = barret_mul_const<Q, W_3S[1]>(x4);
  int16x8_t b2 = vsubq_s16(x1, x4);
  b2 = vsubq_s16(b2, b1);
  b1 = vaddq_s16(b1, x1);

  b1 = barret_mul_const<Q, W_9S[1]>(b1);
  b2 = barret_mul_const<Q, W_9S[2]>(b2);

  int16x8_t c0 = x2;
  int16x8_t c1 = barret_mul_const<Q, W_9S[2]>(x2);
  int16x8_t c2 = barret_mul_const<Q, W_9S[4]>(x2);

  btrfly3(a0, b0, c0, h0, h3, h6);
  btrfly3(a1, b1, c1, h1, h4, h7);
  btrfly3(a2, b2, c2, h2, h5, h8);
}

void ntt_9(int16_t ntt[10][9][16], const int16_t poly[800]) {

  for (int i = 0; i < 10; i++) {
    // std::cerr << "ntt_9 input, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

    {
      int16x8_t x0_fr = vld1q_s16(&poly[i * 16]);
      int16x8_t x1_fr = vld1q_s16(&poly[(10 + i) * 16]);
      int16x8_t x2_fr = vld1q_s16(&poly[(20 + i) * 16]);
      int16x8_t x3_fr = vld1q_s16(&poly[(30 + i) * 16]);
      int16x8_t x4_fr = vld1q_s16(&poly[(40 + i) * 16]);

      // debug_int16x8(x0_fr);
      // debug_int16x8(x1_fr);
      // debug_int16x8(x2_fr);
      // debug_int16x8(x3_fr);
      // debug_int16x8(x4_fr);

      int16x8_t h0_fr, h1_fr, h2_fr, h3_fr, h4_fr, h5_fr, h6_fr, h7_fr, h8_fr;
      five_nonezero(x0_fr, x1_fr, x2_fr, x3_fr, x4_fr,
          h0_fr, h1_fr, h2_fr, h3_fr, h4_fr,
          h5_fr, h6_fr, h7_fr, h8_fr);

      // debug_int16x8(h0_fr);
      // debug_int16x8(h1_fr);
      // debug_int16x8(h2_fr);
      // debug_int16x8(h3_fr);
      // debug_int16x8(h4_fr);
      // debug_int16x8(h5_fr);
      // debug_int16x8(h6_fr);
      // debug_int16x8(h7_fr);
      // debug_int16x8(h8_fr);

      h0_fr = barret_mul<Q>(h0_fr, TWISTS[i][0].first, TWISTS[i][0].second);
      h1_fr = barret_mul<Q>(h1_fr, TWISTS[i][1].first, TWISTS[i][1].second);
      h2_fr = barret_mul<Q>(h2_fr, TWISTS[i][2].first, TWISTS[i][2].second);
      h3_fr = barret_mul<Q>(h3_fr, TWISTS[i][3].first, TWISTS[i][3].second);
      h4_fr = barret_mul<Q>(h4_fr, TWISTS[i][4].first, TWISTS[i][4].second);
      h5_fr = barret_mul<Q>(h5_fr, TWISTS[i][5].first, TWISTS[i][5].second);
      h6_fr = barret_mul<Q>(h6_fr, TWISTS[i][6].first, TWISTS[i][6].second);
      h7_fr = barret_mul<Q>(h7_fr, TWISTS[i][7].first, TWISTS[i][7].second);
      h8_fr = barret_mul<Q>(h8_fr, TWISTS[i][8].first, TWISTS[i][8].second);

      // debug_int16x8(h0_fr);
      // debug_int16x8(h1_fr);
      // debug_int16x8(h2_fr);
      // debug_int16x8(h3_fr);
      // debug_int16x8(h4_fr);
      // debug_int16x8(h5_fr);
      // debug_int16x8(h6_fr);
      // debug_int16x8(h7_fr);
      // debug_int16x8(h8_fr);

      vst1q_s16(&ntt[i][0][0], h0_fr);
      vst1q_s16(&ntt[i][1][0], h1_fr);
      vst1q_s16(&ntt[i][2][0], h2_fr);
      vst1q_s16(&ntt[i][3][0], h3_fr);
      vst1q_s16(&ntt[i][4][0], h4_fr);
      vst1q_s16(&ntt[i][5][0], h5_fr);
      vst1q_s16(&ntt[i][6][0], h6_fr);
      vst1q_s16(&ntt[i][7][0], h7_fr);
      vst1q_s16(&ntt[i][8][0], h8_fr);
    }

    {
      int16x8_t x0_bk = vld1q_s16(&poly[i * 16 + 8]);
      int16x8_t x1_bk = vld1q_s16(&poly[(10 + i) * 16 + 8]);
      int16x8_t x2_bk = vld1q_s16(&poly[(20 + i) * 16 + 8]);
      int16x8_t x3_bk = vld1q_s16(&poly[(30 + i) * 16 + 8]);
      int16x8_t x4_bk = vld1q_s16(&poly[(40 + i) * 16 + 8]);

      // debug_int16x8(x0_bk);
      // debug_int16x8(x1_bk);
      // debug_int16x8(x2_bk);
      // debug_int16x8(x3_bk);
      // debug_int16x8(x4_bk);

      int16x8_t h0_bk, h1_bk, h2_bk, h3_bk, h4_bk, h5_bk, h6_bk, h7_bk, h8_bk;
      five_nonezero(x0_bk, x1_bk, x2_bk, x3_bk, x4_bk,
          h0_bk, h1_bk, h2_bk, h3_bk, h4_bk,
          h5_bk, h6_bk, h7_bk, h8_bk);

      h0_bk = barret_mul<Q>(h0_bk, TWISTS[i][0].first, TWISTS[i][0].second);
      h1_bk = barret_mul<Q>(h1_bk, TWISTS[i][1].first, TWISTS[i][1].second);
      h2_bk = barret_mul<Q>(h2_bk, TWISTS[i][2].first, TWISTS[i][2].second);
      h3_bk = barret_mul<Q>(h3_bk, TWISTS[i][3].first, TWISTS[i][3].second);
      h4_bk = barret_mul<Q>(h4_bk, TWISTS[i][4].first, TWISTS[i][4].second);
      h5_bk = barret_mul<Q>(h5_bk, TWISTS[i][5].first, TWISTS[i][5].second);
      h6_bk = barret_mul<Q>(h6_bk, TWISTS[i][6].first, TWISTS[i][6].second);
      h7_bk = barret_mul<Q>(h7_bk, TWISTS[i][7].first, TWISTS[i][7].second);
      h8_bk = barret_mul<Q>(h8_bk, TWISTS[i][8].first, TWISTS[i][8].second);

      vst1q_s16(&ntt[i][0][8], h0_bk);
      vst1q_s16(&ntt[i][1][8], h1_bk);
      vst1q_s16(&ntt[i][2][8], h2_bk);
      vst1q_s16(&ntt[i][3][8], h3_bk);
      vst1q_s16(&ntt[i][4][8], h4_bk);
      vst1q_s16(&ntt[i][5][8], h5_bk);
      vst1q_s16(&ntt[i][6][8], h6_bk);
      vst1q_s16(&ntt[i][7][8], h7_bk);
      vst1q_s16(&ntt[i][8][8], h8_bk);
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
