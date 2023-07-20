#ifndef NEON_ARITH_OPAQUE_H
#define NEON_ARITH_OPAQUE_H

#include <arm_neon.h>
#include <array>
#include <cstdint>

template <int LANE>
int16x8_t vqrdmulhq_laneq_s16_opaque(int16x8_t v1, int16x8_t v2) {
  int16x8_t res;
  asm ("sqrdmulh %[vd].8h, %[v1].8h, %[v2].h[%[lane]]"
      : [vd] "=w" (res)
      : [v1] "w" (v1), [v2] "x" (v2), [lane] "I" (LANE));
  return res;
}

template <int LANE>
int16x8_t vmulq_laneq_s16_opaque(int16x8_t v1, int16x8_t v2) {
  int16x8_t res;
  asm ("mul %[vd].8h, %[v1].8h, %[v2].h[%[lane]]"
      : [vd] "=w" (res)
      : [v1] "w" (v1), [v2] "x" (v2), [lane] "I" (LANE));
  return res;
}

template <int LANE>
void vmlaq_laneq_s16_opaque(int16x8_t &vd, int16x8_t v1, int16x8_t v2) {
  asm ("mla %[vd].8h, %[v1].8h, %[v2].h[%[lane]]"
      : [vd] "+w" (vd)
      : [v1] "w" (v1), [v2] "x" (v2), [lane] "I" (LANE));
}

template <int LANE>
void vmlsq_laneq_s16_opaque(int16x8_t &vd, int16x8_t v1, int16x8_t v2) {
  asm ("mls %[vd].8h, %[v1].8h, %[v2].h[%[lane]]"
      : [vd] "+w" (vd)
      : [v1] "w" (v1), [v2] "x" (v2), [lane] "I" (LANE));
}

template <int16_t MOD, int LANE>
int16x8_t barret_mul_laneq_opaque(int16x8_t v, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16_opaque<LANE>(v, bar);
  int16x8_t res = vmulq_laneq_s16_opaque<LANE>(v, coef);
  vmlsq_laneq_s16_opaque<7>(res, esti, vec_mod);
  return res;
}

template <int16_t MOD, int LANEC, int LANEB, int LANEM>
int16x8_t barret_mul_laneq_mix_opaque(int16x8_t v, int16x8_t coef_bar) {
  int16x8_t esti = vqrdmulhq_laneq_s16_opaque<LANEB>(v, coef_bar);
  int16x8_t res = vmulq_laneq_s16_opaque<LANEC>(v, coef_bar);
  vmlsq_laneq_s16_opaque<LANEM>(res, esti, coef_bar);
  return res;
}

template <int16_t MOD, int LANE>
void barret_mla_laneq_opaque(int16x8_t &vd, int16x8_t v1, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16_opaque<LANE>(v1, bar);
  vmlaq_laneq_s16_opaque<LANE>(vd, v1, coef);
  vmlsq_laneq_s16_opaque<7>(vd, esti, vec_mod);
}

template <int16_t MOD, int LANE>
void barret_mls_laneq_opaque(int16x8_t &vd, int16x8_t v1, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16_opaque<LANE>(v1, bar);
  vmlsq_laneq_s16_opaque<LANE>(vd, v1, coef);
  vmlaq_laneq_s16_opaque<7>(vd, esti, vec_mod);
}

template <int16_t MOD>
void barret_reduce_laneq_opaque(int16x8_t &v, int16x8_t vec_red, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16_opaque<7>(v, vec_red);
  vmlsq_laneq_s16_opaque<7>(v, esti, vec_mod);
}

template <int16_t MOD, int LANE>
int16x8_t barret_mul_2_laneq_opaque(int16x8_t v, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16_opaque<LANE>(v, bar);
  int16x8_t res = vshlq_n_s16(v, 1);
  vmlsq_laneq_s16_opaque<7>(res, esti, vec_mod);
  return res;
}

template <int16_t MOD, int LANE>
int16x8_t barret_mul_n2_laneq_opaque(int16x8_t v, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16_opaque<LANE>(v, bar);
  int16x8_t res = vshlq_n_s16(vnegq_s16(v), 1);
  vmlsq_laneq_s16_opaque<7>(res, esti, vec_mod);
  return res;
}

#endif // NEON_ARITH_OPAQUE_H
