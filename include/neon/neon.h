#ifndef NEON_H
#define NEON_H

#include <cstdint>

// inx_poly length must >= 768
// inx_poly[761 : 768] filled with zero
// out_poly length must >= 768
void mult(const int16_t in1_poly[], const int16_t in2_poly[], int16_t out_poly[]);

#endif // NEON_H
