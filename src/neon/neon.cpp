#include <arm_neon.h>
#include <cstring>
#include <iostream>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"
#include "utils/debug.h"

#include "neon/ntt_10.h"
#include "neon/intt_10_x10.h"
#include "neon/ntt_9.h"
#include "neon/intt_9_x9.h"
#include "neon/base_mul.h"
#include "neon/mult_low.h"

void forward(int16_t in_poly[], int16_t out_ntt[10][9][16]) {
  ntt_10(out_ntt, in_poly);
  ntt_9(out_ntt);
}

constexpr int16_t INV90 = -51;
constexpr int16_t INV90_BAR = gen_bar<int16_t, Q>(INV90);
void div90_main(int16_t main_poly[1440]) {
  int16x8_t zeros = vdupq_n_s16(0);
  for (int i = 0; i < 1440; i += 8) {
    int16x8_t chunk = vld1q_s16(&main_poly[i]);
    chunk = barret_mul_const<Q, INV90>(chunk);
    vst1q_s16(&main_poly[i], chunk);
  }
}

void backward(int16_t in_ntt[10][9][16], int16_t out_main[1440]) {
  intt_10_x10(in_ntt);
  intt_9_x9(in_ntt, out_main);
  div90_main(out_main);
}

// main_poly length must >= 1448
// poly length must >= 768
void crt(int16_t poly[], int16_t main_poly[], int16_t low[]) {
  main_poly[1440] = main_poly[0];
  for (int i = 0; i < 680; i += 8) {
    int16x8_t a = vld1q_s16(&main_poly[i + 760]);
    int16x8_t b = vld1q_s16(&main_poly[i + 761]);
    int16x8_t c;
    c = vaddq_s16(a, b);
    vst1q_s16(&poly[i], c);
  }
  for (int i = 680; i < 760; i += 8) {
    int16x8_t a = vld1q_s16(&main_poly[i - 680]);
    int16x8_t b = vld1q_s16(&main_poly[i - 679]);
    int16x8_t c;
    c = vaddq_s16(a, b);
    vst1q_s16(&poly[i], c);
  }
  poly[760] = main_poly[80];
  poly[0] -= main_poly[760];
  for (int i = 81; i < 761; i += 8) {
    int16x8_t a = vld1q_s16(&main_poly[i]);
    int16x8_t c = vld1q_s16(&poly[i]);
    c = vaddq_s16(c, a);
    vst1q_s16(&poly[i], c);
  }
  for (int i = 680; i < 760; i += 8) {
    int16x8_t a = vld1q_s16(&low[i - 680]);
    int16x8_t b = vld1q_s16(&low[i - 679]);
    a = vaddq_s16(a, b);
    int16x8_t c = vld1q_s16(&poly[i]);
    c = vsubq_s16(c, a);
    vst1q_s16(&poly[i], c);
  }
  poly[679] -= low[0];
  poly[760] -= low[80];
  for (int i = 0; i < 80; i += 8) {
    int16x8_t a = vld1q_s16(&low[i]);
    int16x8_t c = vld1q_s16(&poly[i]);
    c = vaddq_s16(c, a);
    vst1q_s16(&poly[i], c);
  }
  poly[80] += low[80];
}

// poly length must >= 768
void center_poly(int16_t poly[]) {
  int16x8_t qs = vdupq_n_s16(Q);
  int16x8_t half_qs = vdupq_n_s16((Q - 1) / 2);
  int16x8_t neg_half_qs = vdupq_n_s16(-(Q - 1) / 2);
  for (int i = 0; i < 768; i += 8 * 4) {
    int16x8x4_t chunks = vld1q_s16_x4(&poly[i]);

    barret_reduce<Q>(chunks.val[0]);
    barret_reduce<Q>(chunks.val[1]);
    barret_reduce<Q>(chunks.val[2]);
    barret_reduce<Q>(chunks.val[3]);

    chunks.val[0] = vsubq_s16(chunks.val[0], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[0], half_qs)), qs));
    chunks.val[1] = vsubq_s16(chunks.val[1], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[1], half_qs)), qs));
    chunks.val[2] = vsubq_s16(chunks.val[2], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[2], half_qs)), qs));
    chunks.val[3] = vsubq_s16(chunks.val[3], vandq_s16(vreinterpretq_s16_u16(vcgtq_s16(chunks.val[3], half_qs)), qs));

    chunks.val[0] = vaddq_s16(chunks.val[0], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[0], neg_half_qs)), qs));
    chunks.val[1] = vaddq_s16(chunks.val[1], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[1], neg_half_qs)), qs));
    chunks.val[2] = vaddq_s16(chunks.val[2], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[2], neg_half_qs)), qs));
    chunks.val[3] = vaddq_s16(chunks.val[3], vandq_s16(vreinterpretq_s16_u16(vcltq_s16(chunks.val[3], neg_half_qs)), qs));

    vst1q_s16_x4(&poly[i], chunks);
  }
}

// inx_poly length must >= 768
// inx_poly[761 : 768] filled with zero
// out_poly length must >= 768
void mult(int16_t in1_poly[], int16_t in2_poly[], int16_t out_poly[]) {
  int16_t in1_ntt[10][9][16];
  int16_t in2_ntt[10][9][16];
  int16_t out_ntt[10][9][16];
  int16_t out_low[81];
  int16_t out_main[1448];

  forward(in1_poly, in1_ntt);
  forward(in2_poly, in2_ntt);
  mult_low(in1_poly, in2_poly, out_low);
  base_mul(in1_ntt, in2_ntt, out_ntt);
  backward(out_ntt, out_main);
  // for (int i = 0; i < 1440; i++) {
  //   std::cerr << out_main[i] << " \n"[i == 1439];
  // }
  crt(out_poly, out_main, out_low);
  center_poly(out_poly);
}
