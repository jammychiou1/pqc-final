#include "neon/mult_low.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"

#include "arith_tmpl/arith.h"

constexpr int16_t K1 = -502; // -(W_5 + W_5^4)
constexpr int16_t K2 = 503; // -(W_5^2 + W_5^3)
constexpr int16_t K3 = -459; // -(W_5 - W_5^4)
constexpr int16_t K4 = -868; // -(W_5^2 - W_5^3)
constexpr int16_t K5 = 1327; // W_5 + W_5^2 - W_5^3 - W_5^4
constexpr int16_t INV2 = -2295;
constexpr int16_t INV10 = -459;

void low_ntt_10(int16_t ntt[10][16], const int16_t low[96]) {
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
      h0_fr = barret_mul_const<Q, -2>(h0_fr);

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
      h5_fr = barret_mul_const<Q, -2>(h5_fr);

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
    // h2_fr = barret_mul_const<Q, -INV2>(h2_fr);
    // h8_fr = barret_mul_const<Q, -INV2>(h8_fr);
    // h4_fr = barret_mul_const<Q, -INV2>(h4_fr);
    // h6_fr = barret_mul_const<Q, -INV2>(h6_fr);
    int16x8_t dx0_fr = vshlq_n_s16(x0_fr, 1);
    h2_fr = vsubq_s16(h2_fr, dx0_fr);
    h8_fr = vsubq_s16(h8_fr, dx0_fr);
    h4_fr = vsubq_s16(h4_fr, dx0_fr);
    h6_fr = vsubq_s16(h6_fr, dx0_fr);

    int16x8_t h7_fr = vaddq_s16(nc0s_fr, nn0s_fr);
    int16x8_t h3_fr = vsubq_s16(nc0s_fr, nn0s_fr);
    int16x8_t h9_fr = vaddq_s16(nc1s_fr, nn1s_fr);
    int16x8_t h1_fr = vsubq_s16(nc1s_fr, nn1s_fr);
    // h7_fr = barret_mul_const<Q, -INV2>(h7_fr);
    // h3_fr = barret_mul_const<Q, -INV2>(h3_fr);
    // h9_fr = barret_mul_const<Q, -INV2>(h9_fr);
    // h1_fr = barret_mul_const<Q, -INV2>(h1_fr);
    h7_fr = vsubq_s16(h7_fr, dx0_fr);
    h3_fr = vsubq_s16(h3_fr, dx0_fr);
    h9_fr = vsubq_s16(h9_fr, dx0_fr);
    h1_fr = vsubq_s16(h1_fr, dx0_fr);

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
      h0_bk = barret_mul_const<Q, -2>(h0_bk);

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
      h5_bk = barret_mul_const<Q, -2>(h5_bk);

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
    // h2_bk = barret_mul_const<Q, -INV2>(h2_bk);
    // h8_bk = barret_mul_const<Q, -INV2>(h8_bk);
    // h4_bk = barret_mul_const<Q, -INV2>(h4_bk);
    // h6_bk = barret_mul_const<Q, -INV2>(h6_bk);
    int16x8_t dx0_bk = vshlq_n_s16(x0_bk, 1);
    h2_bk = vsubq_s16(h2_bk, dx0_bk);
    h8_bk = vsubq_s16(h8_bk, dx0_bk);
    h4_bk = vsubq_s16(h4_bk, dx0_bk);
    h6_bk = vsubq_s16(h6_bk, dx0_bk);

    int16x8_t h7_bk = vaddq_s16(nc0s_bk, nn0s_bk);
    int16x8_t h3_bk = vsubq_s16(nc0s_bk, nn0s_bk);
    int16x8_t h9_bk = vaddq_s16(nc1s_bk, nn1s_bk);
    int16x8_t h1_bk = vsubq_s16(nc1s_bk, nn1s_bk);
    // h7_bk = barret_mul_const<Q, -INV2>(h7_bk);
    // h3_bk = barret_mul_const<Q, -INV2>(h3_bk);
    // h9_bk = barret_mul_const<Q, -INV2>(h9_bk);
    // h1_bk = barret_mul_const<Q, -INV2>(h1_bk);
    h7_bk = vsubq_s16(h7_bk, dx0_bk);
    h3_bk = vsubq_s16(h3_bk, dx0_bk);
    h9_bk = vsubq_s16(h9_bk, dx0_bk);
    h1_bk = vsubq_s16(h1_bk, dx0_bk);

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
}

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;
constexpr int16_t W_10 = gen_pow<int16_t, Q>(W_4590, ORD / 10);

constexpr std::array<int16_t, 10> W_10S = gen_pows<int16_t, 10, Q>(W_10);
constexpr std::array<int16_t, 10> W_10_BARS = gen_bars<int16_t, 10, Q>(W_10S);

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

void low_base_mul(int16_t in1_ntt[10][16], int16_t in2_ntt[10][16], int16_t out_ntt[10][16]) {
  for (int i = 0; i < 10; i++) {
    int16_t twiddle = W_10S[i];
    int16_t twiddle_bar = W_10_BARS[i];

    int16x8_t a_fr = vld1q_s16(&in1_ntt[i][0]);
    int16x8_t a_bk = vld1q_s16(&in1_ntt[i][8]);

    int16x8_t a_fb = a_fr;
    barret_reduce<Q>(a_fb);
    a_fb = vaddq_s16(a_fb, a_bk);

    int16x8_t b_fr = vld1q_s16(&in2_ntt[i][0]);
    int16x8_t b_bk = vld1q_s16(&in2_ntt[i][8]);

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

    vst1q_s16(&out_ntt[i][0], res_fr);
    vst1q_s16(&out_ntt[i][8], res_bk);
  }
}

void low_intt_10(int16_t ntt[10][16], int16_t low[96]) {

  {
    int16x8_t h0_fr = vld1q_s16(&ntt[0][0]);
    int16x8_t h1_fr = vld1q_s16(&ntt[9][0]);
    int16x8_t h2_fr = vld1q_s16(&ntt[8][0]);
    int16x8_t h3_fr = vld1q_s16(&ntt[7][0]);
    int16x8_t h4_fr = vld1q_s16(&ntt[6][0]);
    int16x8_t h5_fr = vld1q_s16(&ntt[5][0]);
    int16x8_t h6_fr = vld1q_s16(&ntt[4][0]);
    int16x8_t h7_fr = vld1q_s16(&ntt[3][0]);
    int16x8_t h8_fr = vld1q_s16(&ntt[2][0]);
    int16x8_t h9_fr = vld1q_s16(&ntt[1][0]);

    int16x8_t nc0a_fr = vaddq_s16(h2_fr, h8_fr);
    int16x8_t nn0a_fr = vsubq_s16(h2_fr, h8_fr);
    int16x8_t nc1a_fr = vaddq_s16(h4_fr, h6_fr);
    int16x8_t nn1a_fr = vsubq_s16(h4_fr, h6_fr);
    int16x8_t nc0s_fr = vaddq_s16(h7_fr, h3_fr);
    int16x8_t nn0s_fr = vsubq_s16(h7_fr, h3_fr);
    int16x8_t nc1s_fr = vaddq_s16(h9_fr, h1_fr);
    int16x8_t nn1s_fr = vsubq_s16(h9_fr, h1_fr);

    int16x8_t a05_fr = vaddq_s16(nc0a_fr, nc1a_fr);
    int16x8_t s05_fr = vaddq_s16(nc0s_fr, nc1s_fr);
    a05_fr = vaddq_s16(a05_fr, h0_fr);
    s05_fr = vaddq_s16(s05_fr, h5_fr);

    nc0a_fr = vsubq_s16(nc0a_fr, nc1a_fr);
    nc0s_fr = vsubq_s16(nc0s_fr, nc1s_fr);

    int16x8_t x0_fr = vaddq_s16(a05_fr, s05_fr);
    int16x8_t x5_fr = vsubq_s16(a05_fr, s05_fr);

    // nc0a_fr = barret_mul_const<Q, -INV2>(nc0a_fr);
    // nn0a_fr = barret_mul_const<Q, -INV2>(nn0a_fr);
    // nc1a_fr = barret_mul_const<Q, -INV2>(nc1a_fr);
    // nn1a_fr = barret_mul_const<Q, -INV2>(nn1a_fr);
    // nc0s_fr = barret_mul_const<Q, -INV2>(nc0s_fr);
    // nn0s_fr = barret_mul_const<Q, -INV2>(nn0s_fr);
    // nc1s_fr = barret_mul_const<Q, -INV2>(nc1s_fr);
    // nn1s_fr = barret_mul_const<Q, -INV2>(nn1s_fr);

    // a part

    int16x8_t a41_fr;
    int16x8_t a23_fr;
    {
      int16x8_t dh0_fr = barret_mul_const<Q, 2>(h0_fr);
      int16x8_t aa_fr = vsubq_s16(nc1a_fr, dh0_fr);
      barret_mla_const<Q, K5>(aa_fr, nn1s_fr);

      int16x8_t b1a_fr = vsubq_s16(nn1s_fr, nn0s_fr);
      int16x8_t b2a_fr = vaddq_s16(nn1s_fr, nn0s_fr);

      a41_fr = barret_mul_const<Q, K3>(b1a_fr);
      barret_mla_const<Q, K1>(a41_fr, nc0a_fr);
      a41_fr = vaddq_s16(a41_fr, aa_fr);

      a23_fr = barret_mul_const<Q, K4>(b2a_fr);
      barret_mla_const<Q, K2>(a23_fr, nc0a_fr);
      a23_fr = vaddq_s16(a23_fr, aa_fr);
    }

    // s part

    int16x8_t s41_fr;
    int16x8_t s23_fr;
    {
      int16x8_t dh5_fr = barret_mul_const<Q, 2>(h5_fr);
      int16x8_t as_fr = vsubq_s16(nc1s_fr, dh5_fr);
      barret_mla_const<Q, K5>(as_fr, nn1a_fr);

      int16x8_t b1s_fr = vsubq_s16(nn1a_fr, nn0a_fr);
      int16x8_t b2s_fr = vaddq_s16(nn1a_fr, nn0a_fr);

      s41_fr = barret_mul_const<Q, K3>(b1s_fr);
      barret_mla_const<Q, K1>(s41_fr, nc0s_fr);
      s41_fr = vaddq_s16(s41_fr, as_fr);

      s23_fr = barret_mul_const<Q, K4>(b2s_fr);
      barret_mla_const<Q, K2>(s23_fr, nc0s_fr);
      s23_fr = vaddq_s16(s23_fr, as_fr);
    }

    int16x8_t x4_fr = vaddq_s16(a41_fr, s41_fr);
    int16x8_t x1_fr = vsubq_s16(a41_fr, s41_fr);
    int16x8_t x2_fr = vaddq_s16(a23_fr, s23_fr);
    int16x8_t x3_fr = vsubq_s16(a23_fr, s23_fr);

    // x0_fr = barret_mul_const<Q, INV10>(x0_fr);
    // x1_fr = barret_mul_const<Q, INV10>(x1_fr);
    // x2_fr = barret_mul_const<Q, INV10>(x2_fr);
    // x3_fr = barret_mul_const<Q, INV10>(x3_fr);
    // x4_fr = barret_mul_const<Q, INV10>(x4_fr);
    // x5_fr = barret_mul_const<Q, INV10>(x5_fr);

    x0_fr = barret_mul_const<Q, -18>(x0_fr);
    x1_fr = barret_mul_const<Q, 9>(x1_fr);
    x2_fr = barret_mul_const<Q, 9>(x2_fr);
    x3_fr = barret_mul_const<Q, 9>(x3_fr);
    x4_fr = barret_mul_const<Q, 9>(x4_fr);
    x5_fr = barret_mul_const<Q, -18>(x5_fr);

    vst1q_s16(&low[0], x0_fr);
    vst1q_s16(&low[16], x1_fr);
    vst1q_s16(&low[32], x2_fr);
    vst1q_s16(&low[48], x3_fr);
    vst1q_s16(&low[64], x4_fr);
    vst1q_s16(&low[80], x5_fr);
  }

  {
    int16x8_t h0_bk = vld1q_s16(&ntt[0][8]);
    int16x8_t h1_bk = vld1q_s16(&ntt[9][8]);
    int16x8_t h2_bk = vld1q_s16(&ntt[8][8]);
    int16x8_t h3_bk = vld1q_s16(&ntt[7][8]);
    int16x8_t h4_bk = vld1q_s16(&ntt[6][8]);
    int16x8_t h5_bk = vld1q_s16(&ntt[5][8]);
    int16x8_t h6_bk = vld1q_s16(&ntt[4][8]);
    int16x8_t h7_bk = vld1q_s16(&ntt[3][8]);
    int16x8_t h8_bk = vld1q_s16(&ntt[2][8]);
    int16x8_t h9_bk = vld1q_s16(&ntt[1][8]);

    int16x8_t nc0a_bk = vaddq_s16(h2_bk, h8_bk);
    int16x8_t nn0a_bk = vsubq_s16(h2_bk, h8_bk);
    int16x8_t nc1a_bk = vaddq_s16(h4_bk, h6_bk);
    int16x8_t nn1a_bk = vsubq_s16(h4_bk, h6_bk);
    int16x8_t nc0s_bk = vaddq_s16(h7_bk, h3_bk);
    int16x8_t nn0s_bk = vsubq_s16(h7_bk, h3_bk);
    int16x8_t nc1s_bk = vaddq_s16(h9_bk, h1_bk);
    int16x8_t nn1s_bk = vsubq_s16(h9_bk, h1_bk);

    int16x8_t a05_bk = vaddq_s16(nc0a_bk, nc1a_bk);
    int16x8_t s05_bk = vaddq_s16(nc0s_bk, nc1s_bk);
    a05_bk = vaddq_s16(a05_bk, h0_bk);
    s05_bk = vaddq_s16(s05_bk, h5_bk);

    nc0a_bk = vsubq_s16(nc0a_bk, nc1a_bk);
    nc0s_bk = vsubq_s16(nc0s_bk, nc1s_bk);

    int16x8_t x0_bk = vaddq_s16(a05_bk, s05_bk);
    // int16x8_t x5_bk = vsubq_s16(a05_bk, s05_bk);

    // nc0a_bk = barret_mul_const<Q, -INV2>(nc0a_bk);
    // nn0a_bk = barret_mul_const<Q, -INV2>(nn0a_bk);
    // nc1a_bk = barret_mul_const<Q, -INV2>(nc1a_bk);
    // nn1a_bk = barret_mul_const<Q, -INV2>(nn1a_bk);
    // nc0s_bk = barret_mul_const<Q, -INV2>(nc0s_bk);
    // nn0s_bk = barret_mul_const<Q, -INV2>(nn0s_bk);
    // nc1s_bk = barret_mul_const<Q, -INV2>(nc1s_bk);
    // nn1s_bk = barret_mul_const<Q, -INV2>(nn1s_bk);

    // a part

    int16x8_t a41_bk;
    int16x8_t a23_bk;
    {
      int16x8_t dh0_bk = barret_mul_const<Q, 2>(h0_bk);
      int16x8_t aa_bk = vsubq_s16(nc1a_bk, dh0_bk);
      barret_mla_const<Q, K5>(aa_bk, nn1s_bk);

      int16x8_t b1a_bk = vsubq_s16(nn1s_bk, nn0s_bk);
      int16x8_t b2a_bk = vaddq_s16(nn1s_bk, nn0s_bk);

      a41_bk = barret_mul_const<Q, K3>(b1a_bk);
      barret_mla_const<Q, K1>(a41_bk, nc0a_bk);
      a41_bk = vaddq_s16(a41_bk, aa_bk);

      a23_bk = barret_mul_const<Q, K4>(b2a_bk);
      barret_mla_const<Q, K2>(a23_bk, nc0a_bk);
      a23_bk = vaddq_s16(a23_bk, aa_bk);
    }

    // s part

    int16x8_t s41_bk;
    int16x8_t s23_bk;
    {
      int16x8_t dh5_bk = barret_mul_const<Q, 2>(h5_bk);
      int16x8_t as_bk = vsubq_s16(nc1s_bk, dh5_bk);
      barret_mla_const<Q, K5>(as_bk, nn1a_bk);

      int16x8_t b1s_bk = vsubq_s16(nn1a_bk, nn0a_bk);
      int16x8_t b2s_bk = vaddq_s16(nn1a_bk, nn0a_bk);

      s41_bk = barret_mul_const<Q, K3>(b1s_bk);
      barret_mla_const<Q, K1>(s41_bk, nc0s_bk);
      s41_bk = vaddq_s16(s41_bk, as_bk);

      s23_bk = barret_mul_const<Q, K4>(b2s_bk);
      barret_mla_const<Q, K2>(s23_bk, nc0s_bk);
      s23_bk = vaddq_s16(s23_bk, as_bk);
    }

    int16x8_t x4_bk = vaddq_s16(a41_bk, s41_bk);
    int16x8_t x1_bk = vsubq_s16(a41_bk, s41_bk);
    int16x8_t x2_bk = vaddq_s16(a23_bk, s23_bk);
    int16x8_t x3_bk = vsubq_s16(a23_bk, s23_bk);

    // x0_bk = barret_mul_const<Q, INV10>(x0_bk);
    // x1_bk = barret_mul_const<Q, INV10>(x1_bk);
    // x2_bk = barret_mul_const<Q, INV10>(x2_bk);
    // x3_bk = barret_mul_const<Q, INV10>(x3_bk);
    // x4_bk = barret_mul_const<Q, INV10>(x4_bk);
    // x5_bk = barret_mul_const<Q, INV10>(x5_bk);

    x0_bk = barret_mul_const<Q, -18>(x0_bk);
    x1_bk = barret_mul_const<Q, 9>(x1_bk);
    x2_bk = barret_mul_const<Q, 9>(x2_bk);
    x3_bk = barret_mul_const<Q, 9>(x3_bk);
    x4_bk = barret_mul_const<Q, 9>(x4_bk);
    // x5_bk = barret_mul_const<Q, -18>(x5_bk);

    vst1q_s16(&low[8], x0_bk);
    vst1q_s16(&low[24], x1_bk);
    vst1q_s16(&low[40], x2_bk);
    vst1q_s16(&low[56], x3_bk);
    vst1q_s16(&low[72], x4_bk);
    // vst1q_s16(&low[88], x5_bk);
  }

}

void mult_low(const int16_t in1_low[96], const int16_t in2_low[96], int16_t out_low[96]) {
  static int16_t in1_low_ntt[10][16];
  static int16_t in2_low_ntt[10][16];
  static int16_t out_low_ntt[10][16];

  low_ntt_10(in1_low_ntt, in1_low);
  low_ntt_10(in2_low_ntt, in2_low);
  low_base_mul(in1_low_ntt, in2_low_ntt, out_low_ntt);
  low_intt_10(out_low_ntt, out_low);
  out_low[80] = center_lift<int64_t, Q>(out_low[80] - 720 * (int64_t(in1_low[0]) * in2_low[80] + int64_t(in1_low[80]) * in2_low[0]));
}
