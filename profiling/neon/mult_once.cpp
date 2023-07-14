#include "neon/neon.h"

int16_t a[800];
int16_t b[800];
int16_t c[768];

int main() {
    mult(a, b, c);
}
