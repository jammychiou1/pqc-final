#ifndef PRETTY_H
#define PRETTY_H

#include "mult.h"

inline void scan_poly(int16_t poly[]) {
  for (int i = 0; i < P; i++) {
    std::cin >> poly[i];
  }
}

inline void print_poly(int16_t poly[]) {
  for (int i = 0; i < P; i++) {
    std::cout << poly[i] << " \n"[i == P - 1];
  }
}

inline void print_ntt(int16_t ntt[10][9][16]) {
  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 9; j++) {
      for (int k = 0; k < 16; k++) {
        std::cout << ntt[i][j][k] << " \n"[k == 15];
      }
    }
  }
}

#endif // PRETTY_H
