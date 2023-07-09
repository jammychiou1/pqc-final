#include <stdlib.h>

#include "rader_outer/ntt_17.h"

int16_t ntt[1632];
int16_t poly[1632];

int main() {
  for (int i = 0; i < 1632; i++) {
    poly[i] = rand() % 4591 - 2295;
  }
  for (int t = 0; t < 1000000; t++) {
    poly[0] = rand() % 4591 - 2295;
    ntt_17(ntt, poly);
  }
}
