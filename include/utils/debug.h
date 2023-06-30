#ifndef DEBUG_H
#define DEBUG_H

#include <arm_neon.h>
#include <iostream>

inline void debug_int16x8(int16x8_t val) {
  int16_t tmp[8];
  vst1q_s16(&tmp[0], val);
  for (int i = 0; i < 8; i++) {
    std::cerr << tmp[i] << " \n"[i == 7];
  }
}

inline void debug_int32x4(int32x4_t val) {
  int32_t tmp[4];
  vst1q_s32(&tmp[0], val);
  for (int i = 0; i < 4; i++) {
    std::cerr << tmp[i] << " \n"[i == 3];
  }
}

#endif // DEBUG_H
