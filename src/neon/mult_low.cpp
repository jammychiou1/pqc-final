#include "neon/mult_low.h"

#include <arm_neon.h>
#include <array>
#include <cstdint>
#include <iostream>

#include "sntrup761.h"
#include "arith_tmpl/gen_const.h"
#include "arith_tmpl/neon_arith.h"
#include "utils/debug.h"

#include "arith_tmpl/arith.h"

void low_ntt_10(int16_t ntt[10][16], int16_t low[80]) {

}

void mult_low(int16_t in1_low[160], int16_t in2_low[160], int16_t out_low[81]) {
  int16_t in1_low_ntt[10][16];
  int16_t in2_low_ntt[10][16];
  int16_t out_low_ntt[10][16];

  for (int i = 0; i < 81; i++) {
    for (int j = 0; i + j < 81; j++) {
      out_low[i + j] = center_lift<int64_t, Q>(out_low[i + j] + int64_t(1) * in1_low[i] * in2_low[j]);
    }
  }
}

