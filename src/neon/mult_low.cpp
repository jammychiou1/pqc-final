#include "neon/mult_low.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>
#include <iostream>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"
#include "utils/debug.h"

#include "arith_tmpl/arith.h"

constexpr int16_t K1 = -502; // -(W_5 + W_5^4)
constexpr int16_t K2 = 503; // -(W_5^2 + W_5^3)
constexpr int16_t K3 = -459; // -(W_5 - W_5^4)
constexpr int16_t K4 = -868; // -(W_5^2 - W_5^3)
constexpr int16_t K5 = 1327; // W_5 + W_5^2 - W_5^3 - W_5^4
constexpr int16_t INV2 = -2295;

void low_ntt_10(int16_t ntt[10][16], int16_t low[80]) {
  {
    int16x8_t x0_fr = vld1q_s16(&low[0]);
    int16x8_t x1_fr = vld1q_s16(&low[16]);
    int16x8_t x2_fr = vld1q_s16(&low[32]);
    int16x8_t x3_fr = vld1q_s16(&low[48]);
    int16x8_t x4_fr = vld1q_s16(&low[64]);

    // a part

    int16x8_t h0_fr;
    int16x8_t nc0a_fr;
    int16x8_t nc1a_fr;
    int16x8_t nn0s_fr;
    int16x8_t nn1s_fr;

    {
      int16x8_t a41_fr = vaddq_s16(x4_fr, x1_fr);
      int16x8_t a23_fr = vaddq_s16(x2_fr, x3_fr);

      int16x8_t aa_fr = vaddq_s16(a41_fr, a23_fr);
      h0_fr = vaddq_s16(aa_fr, x0_fr);

      nc0a_fr = barret_mul_const<Q, K1>(a41_fr);
      barret_mla_const<Q, K2>(nc0a_fr, a23_fr);

      nc1a_fr = vsubq_s16(aa_fr, nc0a_fr);

      int16x8_t b1a_fr = barret_mul_const<Q, K3>(a41_fr);
      int16x8_t b2a_fr = barret_mul_const<Q, K4>(a23_fr);

      nn0s_fr = vsubq_s16(b2a_fr, b1a_fr);
      nn1s_fr = vaddq_s16(b2a_fr, b1a_fr);
      barret_mla_const<Q, K5>(nn1s_fr, aa_fr);
    }

    // s part

    int16x8_t h5_fr;
    int16x8_t nc0s_fr;
    int16x8_t nc1s_fr;
    int16x8_t nn0a_fr;
    int16x8_t nn1a_fr;

    {
      int16x8_t s41_fr = vsubq_s16(x4_fr, x1_fr);
      int16x8_t s23_fr = vsubq_s16(x2_fr, x3_fr);

      int16x8_t as_fr = vaddq_s16(s41_fr, s23_fr);
      h5_fr = vaddq_s16(as_fr, x0_fr);

      nc0s_fr = barret_mul_const<Q, K1>(s41_fr);
      barret_mla_const<Q, K2>(nc0s_fr, s23_fr);

      nc1s_fr = vsubq_s16(as_fr, nc0s_fr);

      int16x8_t b1s_fr = barret_mul_const<Q, K3>(s41_fr);
      int16x8_t b2s_fr = barret_mul_const<Q, K4>(s23_fr);

      nn0a_fr = vsubq_s16(b2s_fr, b1s_fr);
      nn1a_fr = vaddq_s16(b2s_fr, b1s_fr);
      barret_mla_const<Q, K5>(nn1a_fr, as_fr);
    }

    int16x8_t h2_fr = vaddq_s16(nc0a_fr, nn0a_fr);
    int16x8_t h8_fr = vsubq_s16(nc0a_fr, nn0a_fr);
    int16x8_t h4_fr = vaddq_s16(nc1a_fr, nn1a_fr);
    int16x8_t h6_fr = vsubq_s16(nc1a_fr, nn1a_fr);
    h2_fr = barret_mul_const<Q, -INV2>(h2_fr);
    h8_fr = barret_mul_const<Q, -INV2>(h8_fr);
    h4_fr = barret_mul_const<Q, -INV2>(h4_fr);
    h6_fr = barret_mul_const<Q, -INV2>(h6_fr);
    h2_fr = vaddq_s16(h2_fr, x0_fr);
    h8_fr = vaddq_s16(h8_fr, x0_fr);
    h4_fr = vaddq_s16(h4_fr, x0_fr);
    h6_fr = vaddq_s16(h6_fr, x0_fr);

    int16x8_t h7_fr = vaddq_s16(nc0s_fr, nn0s_fr);
    int16x8_t h3_fr = vsubq_s16(nc0s_fr, nn0s_fr);
    int16x8_t h9_fr = vaddq_s16(nc1s_fr, nn1s_fr);
    int16x8_t h1_fr = vsubq_s16(nc1s_fr, nn1s_fr);
    h7_fr = barret_mul_const<Q, -INV2>(h7_fr);
    h3_fr = barret_mul_const<Q, -INV2>(h3_fr);
    h9_fr = barret_mul_const<Q, -INV2>(h9_fr);
    h1_fr = barret_mul_const<Q, -INV2>(h1_fr);
    h7_fr = vaddq_s16(h7_fr, x0_fr);
    h3_fr = vaddq_s16(h3_fr, x0_fr);
    h9_fr = vaddq_s16(h9_fr, x0_fr);
    h1_fr = vaddq_s16(h1_fr, x0_fr);

    vst1q_s16(&ntt[0][0], h0_fr);
    vst1q_s16(&ntt[1][0], h1_fr);
    vst1q_s16(&ntt[2][0], h2_fr);
    vst1q_s16(&ntt[3][0], h3_fr);
    vst1q_s16(&ntt[4][0], h4_fr);
    vst1q_s16(&ntt[5][0], h5_fr);
    vst1q_s16(&ntt[6][0], h6_fr);
    vst1q_s16(&ntt[7][0], h7_fr);
    vst1q_s16(&ntt[8][0], h8_fr);
    vst1q_s16(&ntt[9][0], h9_fr);
  }

  {
    int16x8_t x0_bk = vld1q_s16(&low[8]);
    int16x8_t x1_bk = vld1q_s16(&low[24]);
    int16x8_t x2_bk = vld1q_s16(&low[40]);
    int16x8_t x3_bk = vld1q_s16(&low[56]);
    int16x8_t x4_bk = vld1q_s16(&low[72]);

    // a part

    int16x8_t h0_bk;
    int16x8_t nc0a_bk;
    int16x8_t nc1a_bk;
    int16x8_t nn0s_bk;
    int16x8_t nn1s_bk;

    {
      int16x8_t a41_bk = vaddq_s16(x4_bk, x1_bk);
      int16x8_t a23_bk = vaddq_s16(x2_bk, x3_bk);

      int16x8_t aa_bk = vaddq_s16(a41_bk, a23_bk);
      h0_bk = vaddq_s16(aa_bk, x0_bk);

      nc0a_bk = barret_mul_const<Q, K1>(a41_bk);
      barret_mla_const<Q, K2>(nc0a_bk, a23_bk);

      nc1a_bk = vsubq_s16(aa_bk, nc0a_bk);

      int16x8_t b1a_bk = barret_mul_const<Q, K3>(a41_bk);
      int16x8_t b2a_bk = barret_mul_const<Q, K4>(a23_bk);

      nn0s_bk = vsubq_s16(b2a_bk, b1a_bk);
      nn1s_bk = vaddq_s16(b2a_bk, b1a_bk);
      barret_mla_const<Q, K5>(nn1s_bk, aa_bk);
    }

    // s part

    int16x8_t h5_bk;
    int16x8_t nc0s_bk;
    int16x8_t nc1s_bk;
    int16x8_t nn0a_bk;
    int16x8_t nn1a_bk;

    {
      int16x8_t s41_bk = vsubq_s16(x4_bk, x1_bk);
      int16x8_t s23_bk = vsubq_s16(x2_bk, x3_bk);

      int16x8_t as_bk = vaddq_s16(s41_bk, s23_bk);
      h5_bk = vaddq_s16(as_bk, x0_bk);

      nc0s_bk = barret_mul_const<Q, K1>(s41_bk);
      barret_mla_const<Q, K2>(nc0s_bk, s23_bk);

      nc1s_bk = vsubq_s16(as_bk, nc0s_bk);

      int16x8_t b1s_bk = barret_mul_const<Q, K3>(s41_bk);
      int16x8_t b2s_bk = barret_mul_const<Q, K4>(s23_bk);

      nn0a_bk = vsubq_s16(b2s_bk, b1s_bk);
      nn1a_bk = vaddq_s16(b2s_bk, b1s_bk);
      barret_mla_const<Q, K5>(nn1a_bk, as_bk);
    }

    int16x8_t h2_bk = vaddq_s16(nc0a_bk, nn0a_bk);
    int16x8_t h8_bk = vsubq_s16(nc0a_bk, nn0a_bk);
    int16x8_t h4_bk = vaddq_s16(nc1a_bk, nn1a_bk);
    int16x8_t h6_bk = vsubq_s16(nc1a_bk, nn1a_bk);
    h2_bk = barret_mul_const<Q, -INV2>(h2_bk);
    h8_bk = barret_mul_const<Q, -INV2>(h8_bk);
    h4_bk = barret_mul_const<Q, -INV2>(h4_bk);
    h6_bk = barret_mul_const<Q, -INV2>(h6_bk);
    h2_bk = vaddq_s16(h2_bk, x0_bk);
    h8_bk = vaddq_s16(h8_bk, x0_bk);
    h4_bk = vaddq_s16(h4_bk, x0_bk);
    h6_bk = vaddq_s16(h6_bk, x0_bk);

    int16x8_t h7_bk = vaddq_s16(nc0s_bk, nn0s_bk);
    int16x8_t h3_bk = vsubq_s16(nc0s_bk, nn0s_bk);
    int16x8_t h9_bk = vaddq_s16(nc1s_bk, nn1s_bk);
    int16x8_t h1_bk = vsubq_s16(nc1s_bk, nn1s_bk);
    h7_bk = barret_mul_const<Q, -INV2>(h7_bk);
    h3_bk = barret_mul_const<Q, -INV2>(h3_bk);
    h9_bk = barret_mul_const<Q, -INV2>(h9_bk);
    h1_bk = barret_mul_const<Q, -INV2>(h1_bk);
    h7_bk = vaddq_s16(h7_bk, x0_bk);
    h3_bk = vaddq_s16(h3_bk, x0_bk);
    h9_bk = vaddq_s16(h9_bk, x0_bk);
    h1_bk = vaddq_s16(h1_bk, x0_bk);

    vst1q_s16(&ntt[0][8], h0_bk);
    vst1q_s16(&ntt[1][8], h1_bk);
    vst1q_s16(&ntt[2][8], h2_bk);
    vst1q_s16(&ntt[3][8], h3_bk);
    vst1q_s16(&ntt[4][8], h4_bk);
    vst1q_s16(&ntt[5][8], h5_bk);
    vst1q_s16(&ntt[6][8], h6_bk);
    vst1q_s16(&ntt[7][8], h7_bk);
    vst1q_s16(&ntt[8][8], h8_bk);
    vst1q_s16(&ntt[9][8], h9_bk);
  }

  // std::cerr << "low_ntt_10 output\n";
  // for (int i = 0; i < 10; i++) {
  //   std::cerr << "  i = " << i << '\n';
  //   std::cerr << "    ";
  //   for (int k = 0; k < 16; k++) {
  //     std::cerr << ntt[i][k]  << " \n"[k == 15];
  //   }
  // }
}

int16_t in1_low_ntt[10][16];
int16_t in2_low_ntt[10][16];
int16_t out_low_ntt[10][16];

void mult_low(int16_t in1_low[81], int16_t in2_low[81], int16_t out_low[81]) {
  // for (int i = 0; i < 81; i++) {
  //   for (int j = 0; i + j < 81; j++) {
  //     out_low[i + j] = center_lift<int64_t, Q>(out_low[i + j] + int64_t(1) * in1_low[i] * in2_low[j]);
  //   }
  // }

  low_ntt_10(in1_low_ntt, in1_low);
  low_ntt_10(in2_low_ntt, in2_low);
}
