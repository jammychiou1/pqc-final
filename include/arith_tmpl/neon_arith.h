#ifndef NEON_ARITH_H
#define NEON_ARITH_H

#include <arm_neon.h>
#include <cstdint>

template <int16_t MOD>
void barret_reduce(int16x8_t &v) {
  constexpr int16_t ONEBAR = gen_bar<int16_t, MOD>(1);
  int16x8_t esti = vqrdmulhq_n_s16(v, ONEBAR);
  v = vmlsq_n_s16(v, esti, MOD);
}

template <int16_t MOD>
int16x8_t barret_mul(int16x8_t v, int16_t coef, int16_t bar) {
  int16x8_t esti = vqrdmulhq_n_s16(v, bar);
  int16x8_t res = vmulq_n_s16(v, coef);
  res = vmlsq_n_s16(res, esti, MOD);
  return res;
}

template <int16_t MOD>
void barret_mla(int16x8_t &vd, int16x8_t v1, int16_t coef, int16_t bar) {
  int16x8_t esti = vqrdmulhq_n_s16(v1, bar);
  vd = vmlaq_n_s16(vd, v1, coef);
  vd = vmlsq_n_s16(vd, esti, MOD);
}

template <int16_t MOD>
void barret_mls(int16x8_t &vd, int16x8_t v1, int16_t coef, int16_t bar) {
  int16x8_t esti = vqrdmulhq_n_s16(v1, bar);
  vd = vmlsq_n_s16(vd, v1, coef);
  vd = vmlaq_n_s16(vd, esti, MOD);
}

template <int16_t MOD, int16_t COEF>
int16x8_t barret_mul_const(int16x8_t v) {
  return barret_mul<MOD>(v, COEF, gen_bar<int16_t, MOD>(COEF));
}

template <int16_t MOD, int16_t COEF>
void barret_mla_const(int16x8_t &vd, int16x8_t v1) {
  barret_mla<MOD>(vd, v1, COEF, gen_bar<int16_t, MOD>(COEF));
}

template <int16_t MOD, int16_t COEF>
void barret_mls_const(int16x8_t &vd, int16x8_t v1) {
  barret_mls<MOD>(vd, v1, COEF, gen_bar<int16_t, MOD>(COEF));
}

#endif // NEON_ARITH_H
