#include "neon/ntt_10.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;

constexpr int16_t W_5 = gen_pow<int16_t, Q>(W_4590, ORD / 5);
constexpr std::array<int16_t, 17> W_5S = gen_pows<int16_t, 17, Q>(W_5);

constexpr int16_t K1 = -502; // -(W_5 + W_5^4)
constexpr int16_t K2 = 459; // W_5 - W_5^4
constexpr int16_t K3 = 503; // -(W_5^2 + W_5^3)
constexpr int16_t K4 = 868; // W_5^2 - W_5^3
constexpr int16_t K5 = -1327; // W_5 + W_5^2 - W_5^3 - W_5^4

inline void btrfly5_xn2(int16x8_t x0, int16x8_t x1, int16x8_t x2, int16x8_t x3, int16x8_t x4,
    int16x8_t &h0, int16x8_t &h1, int16x8_t &h2, int16x8_t &h3, int16x8_t &h4) {

  int16x8_t a14 = vaddq_s16(x1, x4);
  int16x8_t s14 = vsubq_s16(x1, x4);

  int16x8_t a32 = vaddq_s16(x3, x2);
  int16x8_t s32 = vsubq_s16(x3, x2);

  int16x8_t aa = vaddq_s16(a14, a32);
  int16x8_t as = vaddq_s16(s14, s32);

  h0 = vaddq_s16(x0, aa);
  h0 = barret_mul_const<Q, -2>(h0);

  int16x8_t nc0 = barret_mul_const<Q, K1>(a14);
  barret_mla_const<Q, K3>(nc0, a32);
  int16x8_t nc1 = vsubq_s16(aa, nc0);
  barret_reduce<Q>(nc1);

  s14 = barret_mul_const<Q, K2>(s14);
  s32 = barret_mul_const<Q, K4>(s32);
  int16x8_t nn0 = vsubq_s16(s32, s14);
  int16x8_t nn1 = vaddq_s16(s32, s14);
  barret_mla_const<Q, K5>(nn1, as);

  h1 = vaddq_s16(nc0, nn0);
  h2 = vaddq_s16(nc1, nn1);
  h4 = vsubq_s16(nc0, nn0);
  h3 = vsubq_s16(nc1, nn1);

  int16x8_t dx0 = vshlq_n_s16(x0, 1);

  h1 = vsubq_s16(h1, dx0);
  h2 = vsubq_s16(h2, dx0);
  h3 = vsubq_s16(h3, dx0);
  h4 = vsubq_s16(h4, dx0);
}

void ntt_10(int16_t ntt[9][2][10][8]) {
  int16_t (*flat)[10][8] = &ntt[0][0];

  for (int t = 0; t < 18; t++) {
    int16x8x4_t x0123 = vld1q_s16_x4(&flat[t][0][0]);
    int16x8x4_t x4567 = vld1q_s16_x4(&flat[t][4][0]);
    int16x8x2_t x89 = vld1q_s16_x2(&flat[t][8][0]);

    int16x8_t h00, h01, h02, h03, h04, h10, h11, h12, h13, h14;
    btrfly5_xn2(x0123.val[0], x4567.val[2], x0123.val[2], x89.val[0], x4567.val[0],
        h00, h01, h02, h03, h04);
    btrfly5_xn2(x4567.val[1], x0123.val[1], x4567.val[3], x0123.val[3], x89.val[1],
        h10, h11, h12, h13, h14);

    int16x8x4_t h0123;
    int16x8x4_t h4567;
    int16x8x2_t h89;

    h0123.val[0] = vaddq_s16(h00, h10);
    h0123.val[2] = vaddq_s16(h01, h11);
    h4567.val[0] = vaddq_s16(h02, h12);
    h4567.val[2] = vaddq_s16(h03, h13);
    h89.val[0] = vaddq_s16(h04, h14);

    h4567.val[1] = vsubq_s16(h00, h10);
    h4567.val[3] = vsubq_s16(h01, h11);
    h89.val[1] = vsubq_s16(h02, h12);
    h0123.val[1] = vsubq_s16(h03, h13);
    h0123.val[3] = vsubq_s16(h04, h14);

    vst1q_s16_x4(&flat[t][0][0], h0123);
    vst1q_s16_x4(&flat[t][4][0], h4567);
    vst1q_s16_x2(&flat[t][8][0], h89);
  }
}
