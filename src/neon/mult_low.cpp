#include "neon/mult_low.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"
#include "arith_tmpl/neon_arith_opaque.h"

constexpr int16_t K1 = -502; // -(W_5 + W_5^4)
constexpr int16_t K2 = 503; // -(W_5^2 + W_5^3)
constexpr int16_t K3 = -459; // -(W_5 - W_5^4)
constexpr int16_t K4 = -868; // -(W_5^2 - W_5^3)
constexpr int16_t K5 = 1327; // W_5 + W_5^2 - W_5^3 - W_5^4

constexpr static std::array<int16_t, 8> COEFS = {
  -502, // -(W_5 + W_5^4)
  503, // -(W_5^2 + W_5^3)
  -459, // -(W_5 - W_5^4)
  -868, // -(W_5^2 - W_5^3)
  1327, // W_5 + W_5^2 - W_5^3 - W_5^4
  -2,
  2
};

constexpr static std::array<int16_t, 8> COEFS_MOD = back_mod<Q>(COEFS);
constexpr static std::array<int16_t, 8> BARS = gen_bars<int16_t, 8, Q>(COEFS);
constexpr static std::array<int16_t, 8> BARS_RED = back_red<Q>(BARS);

inline void five_nonezero(
    int16x8_t x0, int16x8_t x1, int16x8_t x2, int16x8_t x3, int16x8_t x4,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2, int16x8_t &h3, int16x8_t &h4,
    int16x8_t &h5, int16x8_t &h6, int16x8_t &h7, int16x8_t &h8, int16x8_t &h9) {

  int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
  int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

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
    h0 = barret_mul_n2_laneq_opaque<Q, 5>(h0, coefs_mod, bars_red, coefs_mod);

    nc0a = barret_mul_laneq_opaque<Q, 0>(a41, coefs_mod, bars_red, coefs_mod);
    barret_mla_laneq_opaque<Q, 1>(nc0a, a23, coefs_mod, bars_red, coefs_mod);

    nc1a = vsubq_s16(aa, nc0a);

    int16x8_t b1a = barret_mul_laneq_opaque<Q, 2>(a41, coefs_mod, bars_red, coefs_mod);
    int16x8_t b2a = barret_mul_laneq_opaque<Q, 3>(a23, coefs_mod, bars_red, coefs_mod);

    nn0s = vsubq_s16(b2a, b1a);
    nn1s = vaddq_s16(b2a, b1a);
    barret_mla_laneq_opaque<Q, 4>(nn1s, aa, coefs_mod, bars_red, coefs_mod);
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
    h5 = barret_mul_n2_laneq_opaque<Q, 5>(h5, coefs_mod, bars_red, coefs_mod);

    nc0s = barret_mul_laneq_opaque<Q, 0>(s41, coefs_mod, bars_red, coefs_mod);
    barret_mla_laneq_opaque<Q, 1>(nc0s, s23, coefs_mod, bars_red, coefs_mod);

    nc1s = vsubq_s16(as, nc0s);

    int16x8_t b1s = barret_mul_laneq_opaque<Q, 2>(s41, coefs_mod, bars_red, coefs_mod);
    int16x8_t b2s = barret_mul_laneq_opaque<Q, 3>(s23, coefs_mod, bars_red, coefs_mod);

    nn0a = vsubq_s16(b2s, b1s);
    nn1a = vaddq_s16(b2s, b1s);
    barret_mla_laneq_opaque<Q, 4>(nn1a, as, coefs_mod, bars_red, coefs_mod);
  }

  int16x8_t dx0 = vshlq_n_s16(x0, 1);

  h2 = vaddq_s16(nc0a, nn0a);
  h8 = vsubq_s16(nc0a, nn0a);
  h4 = vaddq_s16(nc1a, nn1a);
  h6 = vsubq_s16(nc1a, nn1a);
  h2 = vsubq_s16(h2, dx0);
  h8 = vsubq_s16(h8, dx0);
  h4 = vsubq_s16(h4, dx0);
  h6 = vsubq_s16(h6, dx0);

  h7 = vaddq_s16(nc0s, nn0s);
  h3 = vsubq_s16(nc0s, nn0s);
  h9 = vaddq_s16(nc1s, nn1s);
  h1 = vsubq_s16(nc1s, nn1s);
  h7 = vsubq_s16(h7, dx0);
  h3 = vsubq_s16(h3, dx0);
  h9 = vsubq_s16(h9, dx0);
  h1 = vsubq_s16(h1, dx0);
}

void low_ntt_10(int16_t ntt[10][16], const int16_t low[96]) {

  {
    int16x8_t x0_fr = vld1q_s16(&low[0]);
    int16x8_t x1_fr = vld1q_s16(&low[16]);
    int16x8_t x2_fr = vld1q_s16(&low[32]);
    int16x8_t x3_fr = vld1q_s16(&low[48]);
    int16x8_t x4_fr = vld1q_s16(&low[64]);

    int16x8_t h0_fr, h1_fr, h2_fr, h3_fr, h4_fr;
    int16x8_t h5_fr, h6_fr, h7_fr, h8_fr, h9_fr;

    five_nonezero(x0_fr, x1_fr, x2_fr, x3_fr, x4_fr,
        h0_fr, h1_fr, h2_fr, h3_fr, h4_fr, h5_fr, h6_fr, h7_fr, h8_fr, h9_fr);

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

    int16x8_t h0_bk, h1_bk, h2_bk, h3_bk, h4_bk;
    int16x8_t h5_bk, h6_bk, h7_bk, h8_bk, h9_bk;

    five_nonezero(x0_bk, x1_bk, x2_bk, x3_bk, x4_bk,
        h0_bk, h1_bk, h2_bk, h3_bk, h4_bk, h5_bk, h6_bk, h7_bk, h8_bk, h9_bk);

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
    a_fb = vaddq_s16(a_fb, a_bk);

    int16x8_t b_fr = vld1q_s16(&in2_ntt[i][0]);
    int16x8_t b_bk = vld1q_s16(&in2_ntt[i][8]);

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

inline void six_nonezero_transpose(
    int16x8_t h0, int16x8_t h1, int16x8_t h2, int16x8_t h3, int16x8_t h4,
    int16x8_t h5, int16x8_t h6, int16x8_t h7, int16x8_t h8, int16x8_t h9,
    int16x8_t &x0, int16x8_t &x1, int16x8_t &x2, int16x8_t &x3, int16x8_t &x4, int16x8_t &x5) {

  int16x8_t coefs_mod = vld1q_s16(&COEFS_MOD[0]);
  int16x8_t bars_red = vld1q_s16(&BARS_RED[0]);

  int16x8_t nc0a = vaddq_s16(h2, h8);
  int16x8_t nn0a = vsubq_s16(h2, h8);
  int16x8_t nc1a = vaddq_s16(h4, h6);
  int16x8_t nn1a = vsubq_s16(h4, h6);
  int16x8_t nc0s = vaddq_s16(h7, h3);
  int16x8_t nn0s = vsubq_s16(h7, h3);
  int16x8_t nc1s = vaddq_s16(h9, h1);
  int16x8_t nn1s = vsubq_s16(h9, h1);

  int16x8_t a05 = vaddq_s16(nc0a, nc1a);
  int16x8_t s05 = vaddq_s16(nc0s, nc1s);
  a05 = vaddq_s16(a05, h0);
  s05 = vaddq_s16(s05, h5);

  nc0a = vsubq_s16(nc0a, nc1a);
  nc0s = vsubq_s16(nc0s, nc1s);

  x0 = vaddq_s16(a05, s05);
  x5 = vsubq_s16(a05, s05);

  // a part

  int16x8_t a41;
  int16x8_t a23;
  {
    int16x8_t dh0 = barret_mul_2_laneq_opaque<Q, 6>(h0, coefs_mod, bars_red, coefs_mod);
    int16x8_t aa = vsubq_s16(nc1a, dh0);
    barret_mla_laneq_opaque<Q, 4>(aa, nn1s, coefs_mod, bars_red, coefs_mod);

    int16x8_t b1a = vsubq_s16(nn1s, nn0s);
    int16x8_t b2a = vaddq_s16(nn1s, nn0s);

    a41 = barret_mul_laneq_opaque<Q, 2>(b1a, coefs_mod, bars_red, coefs_mod);
    barret_mla_laneq_opaque<Q, 0>(a41, nc0a, coefs_mod, bars_red, coefs_mod);
    a41 = vaddq_s16(a41, aa);

    a23 = barret_mul_laneq_opaque<Q, 3>(b2a, coefs_mod, bars_red, coefs_mod);
    barret_mla_laneq_opaque<Q, 1>(a23, nc0a, coefs_mod, bars_red, coefs_mod);
    a23 = vaddq_s16(a23, aa);
  }

  // s part

  int16x8_t s41;
  int16x8_t s23;
  {
    int16x8_t dh5 = barret_mul_2_laneq_opaque<Q, 6>(h5, coefs_mod, bars_red, coefs_mod);
    int16x8_t as = vsubq_s16(nc1s, dh5);
    barret_mla_laneq_opaque<Q, 4>(as, nn1a, coefs_mod, bars_red, coefs_mod);

    int16x8_t b1s = vsubq_s16(nn1a, nn0a);
    int16x8_t b2s = vaddq_s16(nn1a, nn0a);

    s41 = barret_mul_laneq_opaque<Q, 2>(b1s, coefs_mod, bars_red, coefs_mod);
    barret_mla_laneq_opaque<Q, 0>(s41, nc0s, coefs_mod, bars_red, coefs_mod);
    s41 = vaddq_s16(s41, as);

    s23 = barret_mul_laneq_opaque<Q, 3>(b2s, coefs_mod, bars_red, coefs_mod);
    barret_mla_laneq_opaque<Q, 1>(s23, nc0s, coefs_mod, bars_red, coefs_mod);
    s23 = vaddq_s16(s23, as);
  }

  x4 = vaddq_s16(a41, s41);
  x1 = vsubq_s16(a41, s41);
  x2 = vaddq_s16(a23, s23);
  x3 = vsubq_s16(a23, s23);
}

constexpr static std::array<int16_t, 8> COEF_BAR_9_N19 = {
  9, -18, gen_bar<int16_t, Q>(9), gen_bar<int16_t, Q>(-18), Q
};

void low_intt_10(int16_t ntt[10][16], int16_t low[96]) {

  int16x8_t cb_9_n19 = vld1q_s16(&COEF_BAR_9_N19[0]);

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

    int16x8_t x0_fr, x1_fr, x2_fr, x3_fr, x4_fr, x5_fr;
    six_nonezero_transpose(
        h0_fr, h1_fr, h2_fr, h3_fr, h4_fr,
        h5_fr, h6_fr, h7_fr, h8_fr, h9_fr,
        x0_fr, x1_fr, x2_fr, x3_fr, x4_fr, x5_fr);

    x0_fr = barret_mul_laneq_mix_opaque<Q, 1, 3, 4>(x0_fr, cb_9_n19);
    x1_fr = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x1_fr, cb_9_n19);
    x2_fr = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x2_fr, cb_9_n19);
    x3_fr = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x3_fr, cb_9_n19);
    x4_fr = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x4_fr, cb_9_n19);
    x5_fr = barret_mul_laneq_mix_opaque<Q, 1, 3, 4>(x5_fr, cb_9_n19);

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

    int16x8_t x0_bk, x1_bk, x2_bk, x3_bk, x4_bk, x5_bk;
    six_nonezero_transpose(
        h0_bk, h1_bk, h2_bk, h3_bk, h4_bk,
        h5_bk, h6_bk, h7_bk, h8_bk, h9_bk,
        x0_bk, x1_bk, x2_bk, x3_bk, x4_bk, x5_bk);

    x0_bk = barret_mul_laneq_mix_opaque<Q, 1, 3, 4>(x0_bk, cb_9_n19);
    x1_bk = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x1_bk, cb_9_n19);
    x2_bk = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x2_bk, cb_9_n19);
    x3_bk = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x3_bk, cb_9_n19);
    x4_bk = barret_mul_laneq_mix_opaque<Q, 0, 2, 4>(x4_bk, cb_9_n19);
    // x5_bk = barret_mul_laneq_mix_opaque<Q, 1, 3, 4>(x5_bk, cb_9_n19);

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
  int64_t tmp = (out_low[80] - 720 * (int64_t(in1_low[0]) * in2_low[80] + int64_t(in1_low[80]) * in2_low[0]));
  int32_t esti = (tmp * 935519 + (1ll << 31)) >> 32;
  out_low[80] = tmp - 4591 * esti;
}
