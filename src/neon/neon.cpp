#include <arm_neon.h>
#include <iostream>
#include <cstring>
#include <cassert>
#include <array>

#include "mult.h"
#include "pretty.h"

template <size_t SZ>
constexpr std::array<int16_t, SZ> gen_pows(int16_t base) {
  std::array<int16_t, SZ> pows = {1};
  for (int i = 1; i < SZ; i++) {
    pows[i] = center_lift(int32_t(pows[i - 1]) * base);
  }
  return pows;
}

constexpr int32_t BARRET_Q = 467759;
constexpr int32_t BARRET_xW_9_Q = -751221681;
constexpr int32_t BARRET_xW_9_p2_Q = -423790064;
constexpr int32_t BARRET_xW_9_p4_Q = -445774759;
constexpr int32_t BARRET_TWIST_Q[3][3] = {
  {BARRET_Q, BARRET_Q, BARRET_Q},
  {BARRET_Q, BARRET_xW_9_Q, BARRET_xW_9_p2_Q},
  {BARRET_Q, BARRET_xW_9_p2_Q, BARRET_xW_9_p4_Q}
};

constexpr int ORD = 4590;
constexpr int16_t W_4590 = 11;

constexpr std::array<int16_t, ORD> ws = gen_pows<ORD>(W_4590);

constexpr int16_t W_10 = ws[ORD / 10];
constexpr int16_t W_5 = ws[ORD / 5];
constexpr int16_t W_9 = ws[ORD / 9];
constexpr int16_t W_3 = ws[ORD / 3];

constexpr std::array<int16_t, 17> w_5s = gen_pows<17>(W_5);
constexpr std::array<int16_t, 10> w_10s = gen_pows<10>(W_10);
constexpr std::array<int16_t, 17> w_9s = gen_pows<17>(W_9);
constexpr std::array<int16_t, 5> w_3s = gen_pows<5>(W_3);

void debug_int16x8(int16x8_t val) {
  int16_t tmp[8];
  vst1q_s16(&tmp[0], val);
  for (int i = 0; i < 8; i++) {
    std::cerr << tmp[i] << " \n"[i == 7];
  }
}

void debug_int32x4(int32x4_t val) {
  int32_t tmp[4];
  vst1q_s32(&tmp[0], val);
  for (int i = 0; i < 4; i++) {
    std::cerr << tmp[i] << " \n"[i == 3];
  }
}

void ntt_10(int16_t ntt[10][9][16], int16_t poly[1440]) {

  for (int j = 0; j < 9; j++) {
    {
      int32x4_t wide_front_low[5][2];
      int32x4_t wide_front_high[5][2];
      for (int i1 = 0; i1 < 5; i1++) {
        for (int i2 = 0; i2 < 2; i2++) {
          wide_front_low[i1][i2] = vdupq_n_s32(0);
          wide_front_high[i1][i2] = vdupq_n_s32(0);
        }
      }

      for (int i = j, t = 0, idx = j; t < (j <= 2 ? 6 : 5); t++, i = (i + 9) % 10, idx += 9) {
        int ii1 = i % 5;
        int i2 = i % 2;
        int16x8_t front = vld1q_s16(&poly[idx * 16]);
        if (i2 == 0) {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_front_low[i1][0] = vmlal_n_s16(wide_front_low[i1][0], vget_low_s16(front), twiddle);
            wide_front_high[i1][0] = vmlal_high_n_s16(wide_front_high[i1][0], front, twiddle);
          }
        }
        else {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_front_low[i1][1] = vmlal_n_s16(wide_front_low[i1][1], vget_low_s16(front), twiddle);
            wide_front_high[i1][1] = vmlal_high_n_s16(wide_front_high[i1][1], front, twiddle);
          }
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        int16x8_t tmp_front[2];
        int16x8_t result_front[2];
        for (int i2 = 0; i2 < 2; i2++) {
          int32x4_t esti;
          esti = vqrdmulhq_n_s32(wide_front_low[i1][i2], BARRET_Q);
          wide_front_low[i1][i2] = vmlsq_n_s32(wide_front_low[i1][i2], esti, Q);
          esti = vqrdmulhq_n_s32(wide_front_high[i1][i2], BARRET_Q);
          wide_front_high[i1][i2] = vmlsq_n_s32(wide_front_high[i1][i2], esti, Q);

          tmp_front[i2] = vuzp1q_s16(vreinterpretq_s16_s32(wide_front_low[i1][i2]), vreinterpretq_s16_s32(wide_front_high[i1][i2]));
        }
        result_front[0] = vaddq_s16(tmp_front[0], tmp_front[1]);
        result_front[1] = vsubq_s16(tmp_front[0], tmp_front[1]);
        vst1q_s16(&ntt[2 * i1][j][0], result_front[0]);
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][0], result_front[1]);
      }
    }

    {
      int32x4_t wide_back_low[5][2];
      int32x4_t wide_back_high[5][2];
      for (int i1 = 0; i1 < 5; i1++) {
        for (int i2 = 0; i2 < 2; i2++) {
          wide_back_low[i1][i2] = vdupq_n_s32(0);
          wide_back_high[i1][i2] = vdupq_n_s32(0);
        }
      }

      for (int i = j, t = 0, idx = j; t < (j <= 2 ? 6 : 5); t++, i = (i + 9) % 10, idx += 9) {
        int ii1 = i % 5;
        int i2 = i % 2;
        int16x8_t back = vld1q_s16(&poly[idx * 16 + 8]);
        if (i2 == 0) {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_back_low[i1][0] = vmlal_n_s16(wide_back_low[i1][0], vget_low_s16(back), twiddle);
            wide_back_high[i1][0] = vmlal_high_n_s16(wide_back_high[i1][0], back, twiddle);
          }
        }
        else {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_back_low[i1][1] = vmlal_n_s16(wide_back_low[i1][1], vget_low_s16(back), twiddle);
            wide_back_high[i1][1] = vmlal_high_n_s16(wide_back_high[i1][1], back, twiddle);
          }
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        int16x8_t tmp_back[2];
        int16x8_t result_back[2];
        for (int i2 = 0; i2 < 2; i2++) {
          int32x4_t esti;
          esti = vqrdmulhq_n_s32(wide_back_low[i1][i2], BARRET_Q);
          wide_back_low[i1][i2] = vmlsq_n_s32(wide_back_low[i1][i2], esti, Q);
          esti = vqrdmulhq_n_s32(wide_back_high[i1][i2], BARRET_Q);
          wide_back_high[i1][i2] = vmlsq_n_s32(wide_back_high[i1][i2], esti, Q);

          tmp_back[i2] = vuzp1q_s16(vreinterpretq_s16_s32(wide_back_low[i1][i2]), vreinterpretq_s16_s32(wide_back_high[i1][i2]));
        }
        result_back[0] = vaddq_s16(tmp_back[0], tmp_back[1]);
        result_back[1] = vsubq_s16(tmp_back[0], tmp_back[1]);
        vst1q_s16(&ntt[2 * i1][j][8], result_back[0]);
        vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][8], result_back[1]);
      }
    }

    // std::cerr << "ntt_10, j = " << j << '\n';
    // for (int i = 0; i < 10; i++) {
    //   std::cerr << "  i = " << i << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }
  }

}

void intt_10_x10(int16_t ntt[10][9][16], int16_t poly[1440]) {

  for (int j = 0; j < 9; j++) {
    {
      int32x4_t wide_front_low[5][2];
      int32x4_t wide_front_high[5][2];

      // i2 = 0
      for (int ii1 = 0; ii1 < 5; ii1++) {
        int16x8_t front = vld1q_s16(&ntt[ii1 * 4 % 10][j][0]);
        if (ii1 == 0) {
          for (int i1 = 0; i1 < 5; i1++) {
            wide_front_low[i1][0] = vmovl_s16(vget_low_s16(front));
            wide_front_high[i1][0] = vmovl_high_s16(front);
          }
        }
        else {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_front_low[i1][0] = vmlal_n_s16(wide_front_low[i1][0], vget_low_s16(front), twiddle);
            wide_front_high[i1][0] = vmlal_high_n_s16(wide_front_high[i1][0], front, twiddle);
          }
        }
      }

      // i2 = 1
      for (int ii1 = 0; ii1 < 5; ii1++) {
        int16x8_t front = vld1q_s16(&ntt[(ii1 * 4 + 5) % 10][j][0]);
        if (ii1 == 0) {
          for (int i1 = 0; i1 < 5; i1++) {
            wide_front_low[i1][1] = vmovl_s16(vget_low_s16(front));
            wide_front_high[i1][1] = vmovl_high_s16(front);
          }
        }
        else {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_front_low[i1][1] = vmlal_n_s16(wide_front_low[i1][1], vget_low_s16(front), twiddle);
            wide_front_high[i1][1] = vmlal_high_n_s16(wide_front_high[i1][1], front, twiddle);
          }
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        int16x8_t tmp_front[2];
        int16x8_t result_front[2];
        for (int i2 = 0; i2 < 2; i2++) {
          int32x4_t esti;
          esti = vqrdmulhq_n_s32(wide_front_low[i1][i2], BARRET_Q);
          wide_front_low[i1][i2] = vmlsq_n_s32(wide_front_low[i1][i2], esti, Q);
          esti = vqrdmulhq_n_s32(wide_front_high[i1][i2], BARRET_Q);
          wide_front_high[i1][i2] = vmlsq_n_s32(wide_front_high[i1][i2], esti, Q);

          tmp_front[i2] = vuzp1q_s16(vreinterpretq_s16_s32(wide_front_low[i1][i2]), vreinterpretq_s16_s32(wide_front_high[i1][i2]));
        }
        result_front[0] = vaddq_s16(tmp_front[0], tmp_front[1]);
        result_front[1] = vsubq_s16(tmp_front[0], tmp_front[1]);
        // vst1q_s16(&ntt[2 * i1][j][0], result_front[0]);
        // vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][0], result_front[1]);
        vst1q_s16(&poly[((2 * i1) * 81 + 10 * j) % 90 * 16], result_front[0]);
        vst1q_s16(&poly[((2 * i1 + 5) * 81 + 10 * j) % 90 * 16], result_front[1]);
        // std::cerr << i1 << ' ' << j << ": put to " << (((2 * i1) * 81 + 10 * j) % 90 * 16) << " and " << ((2 * i1 + 5) * 81 + 10 * j) % 90 * 16 << '\n';
      }
    }

    {
      int32x4_t wide_back_low[5][2];
      int32x4_t wide_back_high[5][2];

      // i2 = 0
      for (int ii1 = 0; ii1 < 5; ii1++) {
        int16x8_t back = vld1q_s16(&ntt[ii1 * 4 % 10][j][8]);
        if (ii1 == 0) {
          for (int i1 = 0; i1 < 5; i1++) {
            wide_back_low[i1][0] = vmovl_s16(vget_low_s16(back));
            wide_back_high[i1][0] = vmovl_high_s16(back);
          }
        }
        else {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_back_low[i1][0] = vmlal_n_s16(wide_back_low[i1][0], vget_low_s16(back), twiddle);
            wide_back_high[i1][0] = vmlal_high_n_s16(wide_back_high[i1][0], back, twiddle);
          }
        }
      }

      // i2 = 1
      for (int ii1 = 0; ii1 < 5; ii1++) {
        int16x8_t back = vld1q_s16(&ntt[(ii1 * 4 + 5) % 10][j][8]);
        if (ii1 == 0) {
          for (int i1 = 0; i1 < 5; i1++) {
            wide_back_low[i1][1] = vmovl_s16(vget_low_s16(back));
            wide_back_high[i1][1] = vmovl_high_s16(back);
          }
        }
        else {
          for (int i1 = 0; i1 < 5; i1++) {
            int16_t twiddle = w_5s[i1 * ii1];
            wide_back_low[i1][1] = vmlal_n_s16(wide_back_low[i1][1], vget_low_s16(back), twiddle);
            wide_back_high[i1][1] = vmlal_high_n_s16(wide_back_high[i1][1], back, twiddle);
          }
        }
      }

      for (int i1 = 0; i1 < 5; i1++) {
        int16x8_t tmp_back[2];
        int16x8_t result_back[2];
        for (int i2 = 0; i2 < 2; i2++) {
          int32x4_t esti;
          esti = vqrdmulhq_n_s32(wide_back_low[i1][i2], BARRET_Q);
          wide_back_low[i1][i2] = vmlsq_n_s32(wide_back_low[i1][i2], esti, Q);
          esti = vqrdmulhq_n_s32(wide_back_high[i1][i2], BARRET_Q);
          wide_back_high[i1][i2] = vmlsq_n_s32(wide_back_high[i1][i2], esti, Q);

          tmp_back[i2] = vuzp1q_s16(vreinterpretq_s16_s32(wide_back_low[i1][i2]), vreinterpretq_s16_s32(wide_back_high[i1][i2]));
        }
        result_back[0] = vaddq_s16(tmp_back[0], tmp_back[1]);
        result_back[1] = vsubq_s16(tmp_back[0], tmp_back[1]);
        // vst1q_s16(&ntt[2 * i1][j][8], result_back[0]);
        // vst1q_s16(&ntt[(2 * i1 + 5) % 10][j][8], result_back[1]);
        vst1q_s16(&poly[((2 * i1) * 81 + 10 * j) % 90 * 16 + 8], result_back[0]);
        vst1q_s16(&poly[((2 * i1 + 5) * 81 + 10 * j) % 90 * 16 + 8], result_back[1]);
      }
    }

    // std::cerr << "intt_10_x10, j = " << j << '\n';
    // for (int i = 0; i < 10; i++) {
    //   std::cerr << "  i = " << i << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

  }

}

void ntt_9(int16_t ntt[10][9][16]) {

  for (int i = 0; i < 10; i++) {
    // std::cerr << "ntt_9 input, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

    {
      int32x4_t wide_front_low[3][3];
      int32x4_t wide_front_high[3][3];

      for (int j2 = 0; j2 < 3; j2++) {
        for (int jj1 = 0; jj1 < 3; jj1++) {
          int16x8_t front = vld1q_s16(&ntt[i][3 * jj1 + j2][0]);
          if (jj1 == 0) {
            for (int j1 = 0; j1 < 3; j1++) {
              wide_front_low[j1][j2] = vmovl_s16(vget_low_s16(front));
              wide_front_high[j1][j2] = vmovl_high_s16(front);
            }
          }
          else {
            for (int j1 = 0; j1 < 3; j1++) {
              int16_t twiddle = w_3s[j1 * jj1];
              wide_front_low[j1][j2] = vmlal_n_s16(wide_front_low[j1][j2], vget_low_s16(front), twiddle);
              wide_front_high[j1][j2] = vmlal_high_n_s16(wide_front_high[j1][j2], front, twiddle);
            }
          }
        }
        for (int j1 = 1; j1 < 3; j1++) {
          int32x4_t esti;
          int32_t twiddle = w_9s[j1 * j2];

          esti = vqrdmulhq_n_s32(wide_front_low[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_front_low[j1][j2] = vmulq_n_s32(wide_front_low[j1][j2], twiddle);
          wide_front_low[j1][j2] = vmlsq_n_s32(wide_front_low[j1][j2], esti, Q);

          esti = vqrdmulhq_n_s32(wide_front_high[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_front_high[j1][j2] = vmulq_n_s32(wide_front_high[j1][j2], twiddle);
          wide_front_high[j1][j2] = vmlsq_n_s32(wide_front_high[j1][j2], esti, Q);
        }
      }

      for (int j1 = 0; j1 < 3; j1++) {
        int32x4_t result_front_low[3];
        int32x4_t result_front_high[3];

        for (int jj2 = 0; jj2 < 3; jj2++) {
          if (jj2 == 0) {
            for (int j2 = 0; j2 < 3; j2++) {
              result_front_low[j2] = wide_front_low[j1][jj2];
              result_front_high[j2] = wide_front_high[j1][jj2];
            }
          }
          else {
            for (int j2 = 0; j2 < 3; j2++) {
              int32_t twiddle = w_3s[j2 * jj2];
              result_front_low[j2] = vmlaq_n_s32(result_front_low[j2], wide_front_low[j1][jj2], twiddle);
              result_front_high[j2] = vmlaq_n_s32(result_front_high[j2], wide_front_high[j1][jj2], twiddle);
            }
          }
        }
        for (int j2 = 0; j2 < 3; j2++) {
          if (j2 > 0) {
            int32x4_t esti;

            esti = vqrdmulhq_n_s32(result_front_low[j2], BARRET_Q);
            result_front_low[j2] = vmlsq_n_s32(result_front_low[j2], esti, Q);

            esti = vqrdmulhq_n_s32(result_front_high[j2], BARRET_Q);
            result_front_high[j2] = vmlsq_n_s32(result_front_high[j2], esti, Q);
          }
          int16x8_t result_front = vuzp1q_s16(vreinterpretq_s16_s32(result_front_low[j2]), vreinterpretq_s16_s32(result_front_high[j2]));
          vst1q_s16(&ntt[i][j1 + 3 * j2][0], result_front);
        }
      }
    }

    {
      int32x4_t wide_back_low[3][3];
      int32x4_t wide_back_high[3][3];

      for (int j2 = 0; j2 < 3; j2++) {
        for (int jj1 = 0; jj1 < 3; jj1++) {
          int16x8_t back = vld1q_s16(&ntt[i][3 * jj1 + j2][8]);
          if (jj1 == 0) {
            for (int j1 = 0; j1 < 3; j1++) {
              wide_back_low[j1][j2] = vmovl_s16(vget_low_s16(back));
              wide_back_high[j1][j2] = vmovl_high_s16(back);
            }
          }
          else {
            for (int j1 = 0; j1 < 3; j1++) {
              int16_t twiddle = w_3s[j1 * jj1];
              wide_back_low[j1][j2] = vmlal_n_s16(wide_back_low[j1][j2], vget_low_s16(back), twiddle);
              wide_back_high[j1][j2] = vmlal_high_n_s16(wide_back_high[j1][j2], back, twiddle);
            }
          }
        }
        for (int j1 = 1; j1 < 3; j1++) {
          int32x4_t esti;
          int32_t twiddle = w_9s[j1 * j2];

          esti = vqrdmulhq_n_s32(wide_back_low[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_back_low[j1][j2] = vmulq_n_s32(wide_back_low[j1][j2], twiddle);
          wide_back_low[j1][j2] = vmlsq_n_s32(wide_back_low[j1][j2], esti, Q);

          esti = vqrdmulhq_n_s32(wide_back_high[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_back_high[j1][j2] = vmulq_n_s32(wide_back_high[j1][j2], twiddle);
          wide_back_high[j1][j2] = vmlsq_n_s32(wide_back_high[j1][j2], esti, Q);
        }
      }

      for (int j1 = 0; j1 < 3; j1++) {
        int32x4_t result_back_low[3];
        int32x4_t result_back_high[3];

        for (int jj2 = 0; jj2 < 3; jj2++) {
          if (jj2 == 0) {
            for (int j2 = 0; j2 < 3; j2++) {
              result_back_low[j2] = wide_back_low[j1][jj2];
              result_back_high[j2] = wide_back_high[j1][jj2];
            }
          }
          else {
            for (int j2 = 0; j2 < 3; j2++) {
              int32_t twiddle = w_3s[j2 * jj2];
              result_back_low[j2] = vmlaq_n_s32(result_back_low[j2], wide_back_low[j1][jj2], twiddle);
              result_back_high[j2] = vmlaq_n_s32(result_back_high[j2], wide_back_high[j1][jj2], twiddle);
            }
          }
        }
        for (int j2 = 0; j2 < 3; j2++) {
          if (j2 > 0) {
            int32x4_t esti;

            esti = vqrdmulhq_n_s32(result_back_low[j2], BARRET_Q);
            result_back_low[j2] = vmlsq_n_s32(result_back_low[j2], esti, Q);

            esti = vqrdmulhq_n_s32(result_back_high[j2], BARRET_Q);
            result_back_high[j2] = vmlsq_n_s32(result_back_high[j2], esti, Q);
          }
          int16x8_t result_back = vuzp1q_s16(vreinterpretq_s16_s32(result_back_low[j2]), vreinterpretq_s16_s32(result_back_high[j2]));
          vst1q_s16(&ntt[i][j1 + 3 * j2][8], result_back);
        }
      }
    }

    // std::cerr << "ntt_9 output, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }
  }

}

void intt_9_x9(int16_t ntt[10][9][16]) {

  for (int i = 0; i < 10; i++) {
    // std::cerr << "ntt_9 input, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }

    {
      int32x4_t wide_front_low[3][3];
      int32x4_t wide_front_high[3][3];

      for (int j2 = 0; j2 < 3; j2++) {
        for (int jj1 = 0; jj1 < 3; jj1++) {
          int16x8_t front = vld1q_s16(&ntt[i][(6 * jj1 + 8 * j2) % 9][0]);
          if (jj1 == 0) {
            for (int j1 = 0; j1 < 3; j1++) {
              wide_front_low[j1][j2] = vmovl_s16(vget_low_s16(front));
              wide_front_high[j1][j2] = vmovl_high_s16(front);
            }
          }
          else {
            for (int j1 = 0; j1 < 3; j1++) {
              int16_t twiddle = w_3s[j1 * jj1];
              wide_front_low[j1][j2] = vmlal_n_s16(wide_front_low[j1][j2], vget_low_s16(front), twiddle);
              wide_front_high[j1][j2] = vmlal_high_n_s16(wide_front_high[j1][j2], front, twiddle);
            }
          }
        }
        for (int j1 = 1; j1 < 3; j1++) {
          int32x4_t esti;
          int32_t twiddle = w_9s[j1 * j2];

          esti = vqrdmulhq_n_s32(wide_front_low[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_front_low[j1][j2] = vmulq_n_s32(wide_front_low[j1][j2], twiddle);
          wide_front_low[j1][j2] = vmlsq_n_s32(wide_front_low[j1][j2], esti, Q);

          esti = vqrdmulhq_n_s32(wide_front_high[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_front_high[j1][j2] = vmulq_n_s32(wide_front_high[j1][j2], twiddle);
          wide_front_high[j1][j2] = vmlsq_n_s32(wide_front_high[j1][j2], esti, Q);
        }
      }

      for (int j1 = 0; j1 < 3; j1++) {
        int32x4_t result_front_low[3];
        int32x4_t result_front_high[3];

        for (int jj2 = 0; jj2 < 3; jj2++) {
          if (jj2 == 0) {
            for (int j2 = 0; j2 < 3; j2++) {
              result_front_low[j2] = wide_front_low[j1][jj2];
              result_front_high[j2] = wide_front_high[j1][jj2];
            }
          }
          else {
            for (int j2 = 0; j2 < 3; j2++) {
              int32_t twiddle = w_3s[j2 * jj2];
              result_front_low[j2] = vmlaq_n_s32(result_front_low[j2], wide_front_low[j1][jj2], twiddle);
              result_front_high[j2] = vmlaq_n_s32(result_front_high[j2], wide_front_high[j1][jj2], twiddle);
            }
          }
        }
        for (int j2 = 0; j2 < 3; j2++) {
          if (j2 > 0) {
            int32x4_t esti;

            esti = vqrdmulhq_n_s32(result_front_low[j2], BARRET_Q);
            result_front_low[j2] = vmlsq_n_s32(result_front_low[j2], esti, Q);

            esti = vqrdmulhq_n_s32(result_front_high[j2], BARRET_Q);
            result_front_high[j2] = vmlsq_n_s32(result_front_high[j2], esti, Q);
          }
          int16x8_t result_front = vuzp1q_s16(vreinterpretq_s16_s32(result_front_low[j2]), vreinterpretq_s16_s32(result_front_high[j2]));
          vst1q_s16(&ntt[i][j1 + 3 * j2][0], result_front);
        }
      }
    }

    {
      int32x4_t wide_back_low[3][3];
      int32x4_t wide_back_high[3][3];

      for (int j2 = 0; j2 < 3; j2++) {
        for (int jj1 = 0; jj1 < 3; jj1++) {
          int16x8_t back = vld1q_s16(&ntt[i][(6 * jj1 + 8 * j2) % 9][8]);
          if (jj1 == 0) {
            for (int j1 = 0; j1 < 3; j1++) {
              wide_back_low[j1][j2] = vmovl_s16(vget_low_s16(back));
              wide_back_high[j1][j2] = vmovl_high_s16(back);
            }
          }
          else {
            for (int j1 = 0; j1 < 3; j1++) {
              int16_t twiddle = w_3s[j1 * jj1];
              wide_back_low[j1][j2] = vmlal_n_s16(wide_back_low[j1][j2], vget_low_s16(back), twiddle);
              wide_back_high[j1][j2] = vmlal_high_n_s16(wide_back_high[j1][j2], back, twiddle);
            }
          }
        }
        for (int j1 = 1; j1 < 3; j1++) {
          int32x4_t esti;
          int32_t twiddle = w_9s[j1 * j2];

          esti = vqrdmulhq_n_s32(wide_back_low[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_back_low[j1][j2] = vmulq_n_s32(wide_back_low[j1][j2], twiddle);
          wide_back_low[j1][j2] = vmlsq_n_s32(wide_back_low[j1][j2], esti, Q);

          esti = vqrdmulhq_n_s32(wide_back_high[j1][j2], BARRET_TWIST_Q[j1][j2]);
          wide_back_high[j1][j2] = vmulq_n_s32(wide_back_high[j1][j2], twiddle);
          wide_back_high[j1][j2] = vmlsq_n_s32(wide_back_high[j1][j2], esti, Q);
        }
      }

      for (int j1 = 0; j1 < 3; j1++) {
        int32x4_t result_back_low[3];
        int32x4_t result_back_high[3];

        for (int jj2 = 0; jj2 < 3; jj2++) {
          if (jj2 == 0) {
            for (int j2 = 0; j2 < 3; j2++) {
              result_back_low[j2] = wide_back_low[j1][jj2];
              result_back_high[j2] = wide_back_high[j1][jj2];
            }
          }
          else {
            for (int j2 = 0; j2 < 3; j2++) {
              int32_t twiddle = w_3s[j2 * jj2];
              result_back_low[j2] = vmlaq_n_s32(result_back_low[j2], wide_back_low[j1][jj2], twiddle);
              result_back_high[j2] = vmlaq_n_s32(result_back_high[j2], wide_back_high[j1][jj2], twiddle);
            }
          }
        }
        for (int j2 = 0; j2 < 3; j2++) {
          if (j2 > 0) {
            int32x4_t esti;

            esti = vqrdmulhq_n_s32(result_back_low[j2], BARRET_Q);
            result_back_low[j2] = vmlsq_n_s32(result_back_low[j2], esti, Q);

            esti = vqrdmulhq_n_s32(result_back_high[j2], BARRET_Q);
            result_back_high[j2] = vmlsq_n_s32(result_back_high[j2], esti, Q);
          }
          int16x8_t result_back = vuzp1q_s16(vreinterpretq_s16_s32(result_back_low[j2]), vreinterpretq_s16_s32(result_back_high[j2]));
          vst1q_s16(&ntt[i][j1 + 3 * j2][8], result_back);
        }
      }
    }

    // std::cerr << "ntt_9 output, i = " << i << '\n';
    // for (int j = 0; j < 9; j++) {
    //   std::cerr << "  j = " << j << '\n';
    //   std::cerr << "    ";
    //   for (int k = 0; k < 16; k++) {
    //     std::cerr << ntt[i][j][k]  << " \n"[k == 15];
    //   }
    // }
  }

}


void base_mul(int16_t in1_ntt[10][9][16], int16_t in2_ntt[10][9][16], int16_t out_ntt[10][9][16]) {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      // TODO outer product
      int64_t tmp[10][9][16] = {};
      int64_t twiddle = center_lift(int32_t(1) * w_10s[i] * w_9s[j]);
      for (int k1 = 0; k1 < 16; k1++) {
        for (int k2 = 0; k2 < 16; k2++) {
          if (k1 + k2 >= 16) {
            tmp[i][j][k1 + k2 - 16] += twiddle * in1_ntt[i][j][k1] * in2_ntt[i][j][k2];
          }
          else {
            tmp[i][j][k1 + k2] += int64_t(1) * in1_ntt[i][j][k1] * in2_ntt[i][j][k2];
          }
        }
      }
      for (int k = 0; k < 16; k++) {
        out_ntt[i][j][k] = center_lift(tmp[i][j][k]);
      }
    }
  }
}

void forward(int16_t in_poly[], int16_t out_ntt[10][9][16]) {
  ntt_10(out_ntt, in_poly);
  // print_ntt(out_ntt);
  ntt_9(out_ntt);
  // print_ntt(out_ntt);
}


constexpr int32_t INV90 = -51;
constexpr int32_t BARRET_INV90_Q = -23855732;
void div90_main(int16_t main_poly[1440]) {
  int16x8_t zeros = vdupq_n_s16(0);
  for (int i = 0; i < 1440; i += 8) {
    int16x8_t chunk = vld1q_s16(&main_poly[i]);
    int32x4_t wide_low = vmovl_s16(vget_low_s16(chunk));
    int32x4_t wide_high = vmovl_high_s16(chunk);

    int32x4_t esti;

    esti = vqrdmulhq_n_s32(wide_low, BARRET_INV90_Q);
    wide_low = vmulq_n_s32(wide_low, INV90);
    wide_low = vmlsq_n_s32(wide_low, esti, Q);

    esti = vqrdmulhq_n_s32(wide_high, BARRET_INV90_Q);
    wide_high = vmulq_n_s32(wide_high, INV90);
    wide_high = vmlsq_n_s32(wide_high, esti, Q);

    chunk = vuzp1q_s16(vreinterpretq_s16_s32(wide_low), vreinterpretq_s16_s32(wide_high));
    vst1q_s16(&main_poly[i], chunk);
  }
}
void backward(int16_t in_ntt[10][9][16], int16_t out_main[1440]) {
  intt_9_x9(in_ntt);
  intt_10_x10(in_ntt, out_main);
  div90_main(out_main);
  // for (int i = 0; i < 1440; i++) {
  //   int t = i >> 4;
  //   out_main[i] = center_lift(in_ntt[t % 10][t % 9][i % 16] * INV90);
  // }
}

void mult_low(int16_t in1_low[81], int16_t in2_low[81], int16_t out_low[81]) {
  std::memset(out_low, 0, sizeof(int16_t[81]));
  for (int i = 0; i < 81; i++) {
    for (int j = 0; i + j < 81; j++) {
      out_low[i + j] = center_lift(out_low[i + j] + int64_t(1) * in1_low[i] * in2_low[j]);
    }
  }
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
  // for (int i = 81; i < 761; i++) {
  //   poly[i] += main_poly[i];
  // }
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
  // for (int i = 0; i < 80; i++) {
  //   poly[i] += low[i];
  // }
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
  int16x8_t zeros = vdupq_n_s16(0);
  for (int i = 0; i < 768; i += 8) {
    int16x8_t chunk = vld1q_s16(&poly[i]);
    // debug_int16x8(chunk);
    int32x4_t wide_low = vmovl_s16(vget_low_s16(chunk));
    int32x4_t wide_high = vmovl_high_s16(chunk);
    // debug_int32x4(wide_low);
    // debug_int32x4(wide_high);

    int32x4_t esti;

    esti = vqrdmulhq_n_s32(wide_low, BARRET_Q);
    wide_low = vmlsq_n_s32(wide_low, esti, Q);

    esti = vqrdmulhq_n_s32(wide_high, BARRET_Q);
    wide_high = vmlsq_n_s32(wide_high, esti, Q);

    chunk = vuzp1q_s16(vreinterpretq_s16_s32(wide_low), vreinterpretq_s16_s32(wide_high));
    vst1q_s16(&poly[i], chunk);
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
  //   std::cout << out_main[i] << " \n"[i == 1439];
  // }
  crt(out_poly, out_main, out_low);
  center_poly(out_poly);
  // for (int i = 0; i < 761; i++) {
  //   out_poly[i] = center_lift(out_poly[i]);
  // }
}
