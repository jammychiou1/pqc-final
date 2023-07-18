#ifndef BASE_MUL_H
#define BASE_MUL_H

#include <cstdint>

void base_mul(int16_t in1_ntt[9][2][10][8], int16_t in2_ntt[9][2][10][8], int16_t out_ntt[9][2][10][8]);

#endif // BASE_MUL_H
