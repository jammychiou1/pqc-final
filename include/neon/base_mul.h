#ifndef BASE_MUL_H
#define BASE_MUL_H

#include <cstdint>

void base_mul(int16_t in1_ntt[10][9][16], int16_t in2_ntt[10][9][16], int16_t out_ntt[10][9][16]);

#endif // BASE_MUL_H
