#include <stdlib.h>

#include "rader_outer/polymul.h"

int16_t poly_a[1632];
int16_t poly_b[1632];
int16_t poly_c[1632];

int main() {
    polymul(poly_a, poly_b, poly_c);
}
