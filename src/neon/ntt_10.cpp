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
constexpr int16_t W_10 = gen_pow<int16_t, Q>(W_4590, ORD / 10);

constexpr std::array<int16_t, 17> W_5S = gen_pows<int16_t, 17, Q>(W_5);
constexpr std::array<int16_t, 17> W_5_BARS = gen_bars<int16_t, 17, Q>(W_5S);

constexpr std::array<int16_t, 10> W_10S = gen_pows<int16_t, 10, Q>(W_10);
constexpr std::array<int16_t, 10> W_10_BARS = gen_bars<int16_t, 10, Q>(W_10S);

constexpr int16_t K1 = -502; // -(W_5 + W_5^4)
constexpr int16_t K2 = 503; // -(W_5^2 + W_5^3)
constexpr int16_t K3 = -459; // -(W_5 - W_5^4)
constexpr int16_t K4 = -868; // -(W_5^2 - W_5^3)
constexpr int16_t K5 = 1327; // W_5 + W_5^2 - W_5^3 - W_5^4
constexpr int16_t INV2 = -2295;

inline void six_nonzero(
    int16x8_t x0, int16x8_t x1, int16x8_t x2, int16x8_t x3, int16x8_t x4, int16x8_t x5,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2, int16x8_t &h3, int16x8_t &h4,
    int16x8_t &h5, int16x8_t &h6, int16x8_t &h7, int16x8_t &h8, int16x8_t &h9) {

  int16x8_t a05 = vaddq_s16(x0, x5);
  int16x8_t s05 = vsubq_s16(x0, x5);

  // a part

  int16x8_t nc0a;
  int16x8_t nc1a;
  int16x8_t nn0s;
  int16x8_t nn1s;

  {
    int16x8_t a41 = vaddq_s16(x4, x1);
    int16x8_t a23 = vaddq_s16(x2, x3);

    int16x8_t aa = vaddq_s16(a41, a23);
    h0 = vaddq_s16(aa, a05);

    nc0a = barret_mul_const<Q, K1>(a41);
    barret_mla_const<Q, K2>(nc0a, a23);

    nc1a = vsubq_s16(aa, nc0a);

    int16x8_t b1a = barret_mul_const<Q, K3>(a41);
    int16x8_t b2a = barret_mul_const<Q, K4>(a23);

    nn0s = vsubq_s16(b2a, b1a);
    nn1s = vaddq_s16(b2a, b1a);
    barret_mla_const<Q, K5>(nn1s, aa);
  }

  // s part

  int16x8_t nc0s;
  int16x8_t nc1s;
  int16x8_t nn0a;
  int16x8_t nn1a;

  {
    int16x8_t s41 = vsubq_s16(x4, x1);
    int16x8_t s23 = vsubq_s16(x2, x3);

    int16x8_t as = vaddq_s16(s41, s23);
    h5 = vaddq_s16(as, s05);

    nc0s = barret_mul_const<Q, K1>(s41);
    barret_mla_const<Q, K2>(nc0s, s23);

    nc1s = vsubq_s16(as, nc0s);

    int16x8_t b1s = barret_mul_const<Q, K3>(s41);
    int16x8_t b2s = barret_mul_const<Q, K4>(s23);

    nn0a = vsubq_s16(b2s, b1s);
    nn1a = vaddq_s16(b2s, b1s);
    barret_mla_const<Q, K5>(nn1a, as);
  }

  h2 = vaddq_s16(nc0a, nn0a);
  h8 = vsubq_s16(nc0a, nn0a);
  h4 = vaddq_s16(nc1a, nn1a);
  h6 = vsubq_s16(nc1a, nn1a);
  h2 = barret_mul_const<Q, -INV2>(h2);
  h8 = barret_mul_const<Q, -INV2>(h8);
  h4 = barret_mul_const<Q, -INV2>(h4);
  h6 = barret_mul_const<Q, -INV2>(h6);
  h2 = vaddq_s16(h2, a05);
  h8 = vaddq_s16(h8, a05);
  h4 = vaddq_s16(h4, a05);
  h6 = vaddq_s16(h6, a05);

  h7 = vaddq_s16(nc0s, nn0s);
  h3 = vsubq_s16(nc0s, nn0s);
  h9 = vaddq_s16(nc1s, nn1s);
  h1 = vsubq_s16(nc1s, nn1s);
  h7 = barret_mul_const<Q, -INV2>(h7);
  h3 = barret_mul_const<Q, -INV2>(h3);
  h9 = barret_mul_const<Q, -INV2>(h9);
  h1 = barret_mul_const<Q, -INV2>(h1);
  h7 = vaddq_s16(h7, s05);
  h3 = vaddq_s16(h3, s05);
  h9 = vaddq_s16(h9, s05);
  h1 = vaddq_s16(h1, s05);
}

inline void five_nonzero(
    int16x8_t x0, int16x8_t x1, int16x8_t x2, int16x8_t x3, int16x8_t x4,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2, int16x8_t &h3, int16x8_t &h4,
    int16x8_t &h5, int16x8_t &h6, int16x8_t &h7, int16x8_t &h8, int16x8_t &h9) {

  // a part

  int16x8_t nc0a;
  int16x8_t nc1a;
  int16x8_t nn0s;
  int16x8_t nn1s;

  {
    int16x8_t a41 = vaddq_s16(x4, x1);
    int16x8_t a23 = vaddq_s16(x2, x3);

    int16x8_t aa = vaddq_s16(a41, a23);
    h0 = vaddq_s16(aa, x0);

    nc0a = barret_mul_const<Q, K1>(a41);
    barret_mla_const<Q, K2>(nc0a, a23);

    nc1a = vsubq_s16(aa, nc0a);

    int16x8_t b1a = barret_mul_const<Q, K3>(a41);
    int16x8_t b2a = barret_mul_const<Q, K4>(a23);

    nn0s = vsubq_s16(b2a, b1a);
    nn1s = vaddq_s16(b2a, b1a);
    barret_mla_const<Q, K5>(nn1s, aa);
  }

  // s part

  int16x8_t nc0s;
  int16x8_t nc1s;
  int16x8_t nn0a;
  int16x8_t nn1a;

  {
    int16x8_t s41 = vsubq_s16(x4, x1);
    int16x8_t s23 = vsubq_s16(x2, x3);

    int16x8_t as = vaddq_s16(s41, s23);
    h5 = vaddq_s16(as, x0);

    nc0s = barret_mul_const<Q, K1>(s41);
    barret_mla_const<Q, K2>(nc0s, s23);

    nc1s = vsubq_s16(as, nc0s);

    int16x8_t b1s = barret_mul_const<Q, K3>(s41);
    int16x8_t b2s = barret_mul_const<Q, K4>(s23);

    nn0a = vsubq_s16(b2s, b1s);
    nn1a = vaddq_s16(b2s, b1s);
    barret_mla_const<Q, K5>(nn1a, as);
  }

  h2 = vaddq_s16(nc0a, nn0a);
  h8 = vsubq_s16(nc0a, nn0a);
  h4 = vaddq_s16(nc1a, nn1a);
  h6 = vsubq_s16(nc1a, nn1a);
  h2 = barret_mul_const<Q, -INV2>(h2);
  h8 = barret_mul_const<Q, -INV2>(h8);
  h4 = barret_mul_const<Q, -INV2>(h4);
  h6 = barret_mul_const<Q, -INV2>(h6);
  h2 = vaddq_s16(h2, x0);
  h8 = vaddq_s16(h8, x0);
  h4 = vaddq_s16(h4, x0);
  h6 = vaddq_s16(h6, x0);

  h7 = vaddq_s16(nc0s, nn0s);
  h3 = vsubq_s16(nc0s, nn0s);
  h9 = vaddq_s16(nc1s, nn1s);
  h1 = vsubq_s16(nc1s, nn1s);
  h7 = barret_mul_const<Q, -INV2>(h7);
  h3 = barret_mul_const<Q, -INV2>(h3);
  h9 = barret_mul_const<Q, -INV2>(h9);
  h1 = barret_mul_const<Q, -INV2>(h1);
  h7 = vaddq_s16(h7, x0);
  h3 = vaddq_s16(h3, x0);
  h9 = vaddq_s16(h9, x0);
  h1 = vaddq_s16(h1, x0);
}

void ntt_10(int16_t ntt[10][9][16], const int16_t poly[1440]) {

  for (int j = 0; j < 9; j++) {
    if (j <= 2) {
      {
        int16x8_t x0_fr = vld1q_s16(&poly[(45 + j) * 16]);
        int16x8_t x1_fr = vld1q_s16(&poly[(36 + j) * 16]);
        int16x8_t x2_fr = vld1q_s16(&poly[(27 + j) * 16]);
        int16x8_t x3_fr = vld1q_s16(&poly[(18 + j) * 16]);
        int16x8_t x4_fr = vld1q_s16(&poly[(9 + j) * 16]);
        int16x8_t x5_fr = vld1q_s16(&poly[j * 16]);

        int16x8_t h0_fr, h1_fr, h2_fr, h3_fr, h4_fr;
        int16x8_t h5_fr, h6_fr, h7_fr, h8_fr, h9_fr;
        six_nonzero(x0_fr, x1_fr, x2_fr, x3_fr, x4_fr, x5_fr,
            h0_fr, h1_fr, h2_fr, h3_fr, h4_fr,
            h5_fr, h6_fr, h7_fr, h8_fr, h9_fr);

        h0_fr = barret_mul<Q>(h0_fr, W_10S[0], W_10_BARS[0]);
        h1_fr = barret_mul<Q>(h1_fr, W_10S[(j + 5) % 10], W_10_BARS[(j + 5) % 10]);
        h2_fr = barret_mul<Q>(h2_fr, W_10S[2 * (j + 5) % 10], W_10_BARS[2 * (j + 5) % 10]);
        h3_fr = barret_mul<Q>(h3_fr, W_10S[3 * (j + 5) % 10], W_10_BARS[3 * (j + 5) % 10]);
        h4_fr = barret_mul<Q>(h4_fr, W_10S[4 * (j + 5) % 10], W_10_BARS[4 * (j + 5) % 10]);
        h5_fr = barret_mul<Q>(h5_fr, W_10S[5 * (j + 5) % 10], W_10_BARS[5 * (j + 5) % 10]);
        h6_fr = barret_mul<Q>(h6_fr, W_10S[6 * (j + 5) % 10], W_10_BARS[6 * (j + 5) % 10]);
        h7_fr = barret_mul<Q>(h7_fr, W_10S[7 * (j + 5) % 10], W_10_BARS[7 * (j + 5) % 10]);
        h8_fr = barret_mul<Q>(h8_fr, W_10S[8 * (j + 5) % 10], W_10_BARS[8 * (j + 5) % 10]);
        h9_fr = barret_mul<Q>(h9_fr, W_10S[9 * (j + 5) % 10], W_10_BARS[9 * (j + 5) % 10]);

        vst1q_s16(&ntt[0][j][0], h0_fr);
        vst1q_s16(&ntt[1][j][0], h1_fr);
        vst1q_s16(&ntt[2][j][0], h2_fr);
        vst1q_s16(&ntt[3][j][0], h3_fr);
        vst1q_s16(&ntt[4][j][0], h4_fr);
        vst1q_s16(&ntt[5][j][0], h5_fr);
        vst1q_s16(&ntt[6][j][0], h6_fr);
        vst1q_s16(&ntt[7][j][0], h7_fr);
        vst1q_s16(&ntt[8][j][0], h8_fr);
        vst1q_s16(&ntt[9][j][0], h9_fr);
      }

      {
        int16x8_t x0_bk = vld1q_s16(&poly[(45 + j) * 16 + 8]);
        int16x8_t x1_bk = vld1q_s16(&poly[(36 + j) * 16 + 8]);
        int16x8_t x2_bk = vld1q_s16(&poly[(27 + j) * 16 + 8]);
        int16x8_t x3_bk = vld1q_s16(&poly[(18 + j) * 16 + 8]);
        int16x8_t x4_bk = vld1q_s16(&poly[(9 + j) * 16 + 8]);
        int16x8_t x5_bk = vld1q_s16(&poly[j * 16 + 8]);

        int16x8_t h0_bk, h1_bk, h2_bk, h3_bk, h4_bk;
        int16x8_t h5_bk, h6_bk, h7_bk, h8_bk, h9_bk;
        six_nonzero(x0_bk, x1_bk, x2_bk, x3_bk, x4_bk, x5_bk,
            h0_bk, h1_bk, h2_bk, h3_bk, h4_bk,
            h5_bk, h6_bk, h7_bk, h8_bk, h9_bk);

        h0_bk = barret_mul<Q>(h0_bk, W_10S[0], W_10_BARS[0]);
        h1_bk = barret_mul<Q>(h1_bk, W_10S[(j + 5) % 10], W_10_BARS[(j + 5) % 10]);
        h2_bk = barret_mul<Q>(h2_bk, W_10S[2 * (j + 5) % 10], W_10_BARS[2 * (j + 5) % 10]);
        h3_bk = barret_mul<Q>(h3_bk, W_10S[3 * (j + 5) % 10], W_10_BARS[3 * (j + 5) % 10]);
        h4_bk = barret_mul<Q>(h4_bk, W_10S[4 * (j + 5) % 10], W_10_BARS[4 * (j + 5) % 10]);
        h5_bk = barret_mul<Q>(h5_bk, W_10S[5 * (j + 5) % 10], W_10_BARS[5 * (j + 5) % 10]);
        h6_bk = barret_mul<Q>(h6_bk, W_10S[6 * (j + 5) % 10], W_10_BARS[6 * (j + 5) % 10]);
        h7_bk = barret_mul<Q>(h7_bk, W_10S[7 * (j + 5) % 10], W_10_BARS[7 * (j + 5) % 10]);
        h8_bk = barret_mul<Q>(h8_bk, W_10S[8 * (j + 5) % 10], W_10_BARS[8 * (j + 5) % 10]);
        h9_bk = barret_mul<Q>(h9_bk, W_10S[9 * (j + 5) % 10], W_10_BARS[9 * (j + 5) % 10]);

        vst1q_s16(&ntt[0][j][8], h0_bk);
        vst1q_s16(&ntt[1][j][8], h1_bk);
        vst1q_s16(&ntt[2][j][8], h2_bk);
        vst1q_s16(&ntt[3][j][8], h3_bk);
        vst1q_s16(&ntt[4][j][8], h4_bk);
        vst1q_s16(&ntt[5][j][8], h5_bk);
        vst1q_s16(&ntt[6][j][8], h6_bk);
        vst1q_s16(&ntt[7][j][8], h7_bk);
        vst1q_s16(&ntt[8][j][8], h8_bk);
        vst1q_s16(&ntt[9][j][8], h9_bk);
      }
    }
    else {
      {
        int16x8_t x0_fr = vld1q_s16(&poly[(36 + j) * 16]);
        int16x8_t x1_fr = vld1q_s16(&poly[(27 + j) * 16]);
        int16x8_t x2_fr = vld1q_s16(&poly[(18 + j) * 16]);
        int16x8_t x3_fr = vld1q_s16(&poly[(9 + j) * 16]);
        int16x8_t x4_fr = vld1q_s16(&poly[j * 16]);

        int16x8_t h0_fr, h1_fr, h2_fr, h3_fr, h4_fr;
        int16x8_t h5_fr, h6_fr, h7_fr, h8_fr, h9_fr;
        five_nonzero(x0_fr, x1_fr, x2_fr, x3_fr, x4_fr,
            h0_fr, h1_fr, h2_fr, h3_fr, h4_fr,
            h5_fr, h6_fr, h7_fr, h8_fr, h9_fr);

        h0_fr = barret_mul<Q>(h0_fr, W_10S[0], W_10_BARS[0]);
        h1_fr = barret_mul<Q>(h1_fr, W_10S[(j + 6) % 10], W_10_BARS[(j + 6) % 10]);
        h2_fr = barret_mul<Q>(h2_fr, W_10S[2 * (j + 6) % 10], W_10_BARS[2 * (j + 6) % 10]);
        h3_fr = barret_mul<Q>(h3_fr, W_10S[3 * (j + 6) % 10], W_10_BARS[3 * (j + 6) % 10]);
        h4_fr = barret_mul<Q>(h4_fr, W_10S[4 * (j + 6) % 10], W_10_BARS[4 * (j + 6) % 10]);
        h5_fr = barret_mul<Q>(h5_fr, W_10S[5 * (j + 6) % 10], W_10_BARS[5 * (j + 6) % 10]);
        h6_fr = barret_mul<Q>(h6_fr, W_10S[6 * (j + 6) % 10], W_10_BARS[6 * (j + 6) % 10]);
        h7_fr = barret_mul<Q>(h7_fr, W_10S[7 * (j + 6) % 10], W_10_BARS[7 * (j + 6) % 10]);
        h8_fr = barret_mul<Q>(h8_fr, W_10S[8 * (j + 6) % 10], W_10_BARS[8 * (j + 6) % 10]);
        h9_fr = barret_mul<Q>(h9_fr, W_10S[9 * (j + 6) % 10], W_10_BARS[9 * (j + 6) % 10]);

        vst1q_s16(&ntt[0][j][0], h0_fr);
        vst1q_s16(&ntt[1][j][0], h1_fr);
        vst1q_s16(&ntt[2][j][0], h2_fr);
        vst1q_s16(&ntt[3][j][0], h3_fr);
        vst1q_s16(&ntt[4][j][0], h4_fr);
        vst1q_s16(&ntt[5][j][0], h5_fr);
        vst1q_s16(&ntt[6][j][0], h6_fr);
        vst1q_s16(&ntt[7][j][0], h7_fr);
        vst1q_s16(&ntt[8][j][0], h8_fr);
        vst1q_s16(&ntt[9][j][0], h9_fr);
      }

      {
        int16x8_t x0_bk = vld1q_s16(&poly[(36 + j) * 16 + 8]);
        int16x8_t x1_bk = vld1q_s16(&poly[(27 + j) * 16 + 8]);
        int16x8_t x2_bk = vld1q_s16(&poly[(18 + j) * 16 + 8]);
        int16x8_t x3_bk = vld1q_s16(&poly[(9 + j) * 16 + 8]);
        int16x8_t x4_bk = vld1q_s16(&poly[j * 16 + 8]);

        int16x8_t h0_bk, h1_bk, h2_bk, h3_bk, h4_bk;
        int16x8_t h5_bk, h6_bk, h7_bk, h8_bk, h9_bk;
        five_nonzero(x0_bk, x1_bk, x2_bk, x3_bk, x4_bk,
            h0_bk, h1_bk, h2_bk, h3_bk, h4_bk,
            h5_bk, h6_bk, h7_bk, h8_bk, h9_bk);

        h0_bk = barret_mul<Q>(h0_bk, W_10S[0], W_10_BARS[0]);
        h1_bk = barret_mul<Q>(h1_bk, W_10S[(j + 6) % 10], W_10_BARS[(j + 6) % 10]);
        h2_bk = barret_mul<Q>(h2_bk, W_10S[2 * (j + 6) % 10], W_10_BARS[2 * (j + 6) % 10]);
        h3_bk = barret_mul<Q>(h3_bk, W_10S[3 * (j + 6) % 10], W_10_BARS[3 * (j + 6) % 10]);
        h4_bk = barret_mul<Q>(h4_bk, W_10S[4 * (j + 6) % 10], W_10_BARS[4 * (j + 6) % 10]);
        h5_bk = barret_mul<Q>(h5_bk, W_10S[5 * (j + 6) % 10], W_10_BARS[5 * (j + 6) % 10]);
        h6_bk = barret_mul<Q>(h6_bk, W_10S[6 * (j + 6) % 10], W_10_BARS[6 * (j + 6) % 10]);
        h7_bk = barret_mul<Q>(h7_bk, W_10S[7 * (j + 6) % 10], W_10_BARS[7 * (j + 6) % 10]);
        h8_bk = barret_mul<Q>(h8_bk, W_10S[8 * (j + 6) % 10], W_10_BARS[8 * (j + 6) % 10]);
        h9_bk = barret_mul<Q>(h9_bk, W_10S[9 * (j + 6) % 10], W_10_BARS[9 * (j + 6) % 10]);

        vst1q_s16(&ntt[0][j][8], h0_bk);
        vst1q_s16(&ntt[1][j][8], h1_bk);
        vst1q_s16(&ntt[2][j][8], h2_bk);
        vst1q_s16(&ntt[3][j][8], h3_bk);
        vst1q_s16(&ntt[4][j][8], h4_bk);
        vst1q_s16(&ntt[5][j][8], h5_bk);
        vst1q_s16(&ntt[6][j][8], h6_bk);
        vst1q_s16(&ntt[7][j][8], h7_bk);
        vst1q_s16(&ntt[8][j][8], h8_bk);
        vst1q_s16(&ntt[9][j][8], h9_bk);
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
