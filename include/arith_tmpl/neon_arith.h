#ifndef NEON_ARITH_H
#define NEON_ARITH_H

#include <arm_neon.h>
#include <array>
#include <cstdint>

template <int16_t MOD>
constexpr std::array<int16_t, 8> back_mod(std::array<int16_t, 8> vec) {
  vec[7] = MOD;
  return vec;
}

template <int16_t MOD>
constexpr std::array<int16_t, 8> back_red(std::array<int16_t, 8> vec) {
  vec[7] = gen_bar<int16_t, MOD>(1);
  return vec;
}

constexpr int16x8_t arr_to_x8_t(std::array<int16_t, 8> arr) {
  return int16x8_t{arr[0], arr[1], arr[2], arr[3], arr[4], arr[5], arr[6], arr[7]};
}

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

template <int16_t MOD, int LANE>
int16x8_t barret_mul_laneq(int16x8_t v, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16(v, bar, LANE);
  int16x8_t res = vmulq_laneq_s16(v, coef, LANE);
  res = vmlsq_laneq_s16(res, esti, vec_mod, 7);
  return res;
}

template <int16_t MOD, int LANE>
void barret_mla_laneq(int16x8_t &vd, int16x8_t v1, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16(v1, bar, LANE);
  vd = vmlaq_laneq_s16(vd, v1, coef, LANE);
  vd = vmlsq_laneq_s16(vd, esti, vec_mod, 7);
}

template <int16_t MOD, int LANE>
void barret_mls_laneq(int16x8_t &vd, int16x8_t v1, int16x8_t coef, int16x8_t bar, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16(v1, bar, LANE);
  vd = vmlsq_laneq_s16(vd, v1, coef, LANE);
  vd = vmlaq_laneq_s16(vd, esti, vec_mod, 7);
}

template <int16_t MOD>
void barret_reduce_laneq(int16x8_t &v, int16x8_t vec_red, int16x8_t vec_mod) {
  int16x8_t esti = vqrdmulhq_laneq_s16(v, vec_red, 7);
  v = vmlsq_laneq_s16(v, esti, vec_mod, 7);
}

#endif // NEON_ARITH_H
