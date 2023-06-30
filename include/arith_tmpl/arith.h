#ifndef ARITH_H
#define ARITH_H

template <typename INT_T, INT_T MOD>
constexpr INT_T center_lift(INT_T val) {
  if (val >= 0) {
    val %= MOD;
  }
  else {
    val = (MOD - (-val) % MOD) % MOD;
  }
  if (val > (MOD - 1) / 2) {
    val -= MOD;
  }
  return val;
}

#endif // ARITH_H
