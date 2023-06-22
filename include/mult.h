#ifndef MULT_H
#define MULT_H

#include <cstdint>

constexpr int16_t P = 761;
constexpr int16_t Q = 4591;
constexpr int16_t CENTER_MAG = (Q - 1) / 2;

constexpr int64_t center_lift(int64_t val) {
  if (val >= 0) {
    val %= Q;
  }
  else {
    val = (Q - (-val) % Q) % Q;
  }
  if (val > (Q - 1) / 2) {
    val -= Q;
  }
  return val;
}

void mult(int16_t in1[], int16_t in2[], int16_t out[]);

#endif // MULT_H
