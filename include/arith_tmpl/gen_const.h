#ifndef GEN_CONST_H
#define GEN_CONST_H

#include <array>
#include <cmath>
#include <cstdint>

#include "arith_tmpl/arith.h"

template <typename INT_T, INT_T MOD>
constexpr INT_T gen_pow(INT_T base, int pow) {
  int64_t res = 1;
  for (int i = 0; i < pow; i++) {
    res = res * base % MOD;
  }
  return center_lift<int64_t, MOD>(res);
}

template <typename INT_T, size_t SZ, INT_T MOD>
constexpr std::array<INT_T, SZ> gen_pows(INT_T base) {
  std::array<INT_T, SZ> pows = {1};
  for (int i = 1; i < SZ; i++) {
    pows[i] = center_lift<int64_t, MOD>(int64_t(pows[i - 1]) * base);
  }
  return pows;
}

template <typename INT_T, INT_T MOD>
constexpr INT_T gen_bar(INT_T coef) {
  return std::round(double(coef) * (1 << 15) / MOD);
}

template <typename INT_T, size_t SZ, INT_T MOD>
constexpr std::array<INT_T, SZ> gen_bars(std::array<INT_T, SZ> arr) {
  std::array<INT_T, SZ> bars = {};
  for (int i = 0; i < SZ; i++) {
    bars[i] = gen_bar<INT_T, MOD>(arr[i]);
  }
  return bars;
}

#endif // GEN_CONST_H
