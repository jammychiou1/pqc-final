// #include <arm_neon.h>
#include <inttypes.h>
#include <stdio.h>
#include <string.h>
const int32_t W_4590 = 7;
const int32_t W_10 = 2981;
const int32_t W_9 = 2985;
const int16_t Q = 4591;
#define P 761
const int NTT_SZ = sizeof(int16_t) * 10 * 9 * 16;
const int POLY_SZ = sizeof(int16_t) * P;

int32_t w_10s[10];
int32_t w_9s[9];

void init() {
    w_10s[0] = 1;
    for (int i = 1; i < 10; i++) {
        w_10s[i] = (w_10s[i - 1] * W_10) % Q;
    }
    w_9s[0] = 1;
    for (int i = 1; i < 9; i++) {
        w_9s[i] = (w_9s[i - 1] * W_9) % Q;
    }
}

void scan_poly(int16_t poly[]) {
    memset(poly, 0, POLY_SZ);
    for (int i = 0; i < P; i++) {
        scanf("%"PRId16" ", &poly[i]);
    }
}

void print_poly(int16_t poly[]) {
    for (int i = 0; i < P; i++) {
        printf("%"PRId16"%c", poly[i], " \n"[i == P - 1]);
    }
}

void print_ntt(int16_t ntt[10][9][16]) {
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 16; k++) {
                printf("%"PRId16"%c", ntt[i][j][k], " \n"[k == 15]);
            }
        }
        if (i != 9) {
            printf("\n");
        }
    }
}

int64_t tmp[10][9][16];
void ntt10(int16_t ntt[10][9][16]) {
    memset(tmp, 0, sizeof(tmp));
    for (int i = 0; i < 10; i++) {
        for (int ii = 0; ii < 10; ii++) {
            for (int j = 0; j < 9; j++) {
                for (int k = 0; k < 16; k++) {
                    tmp[i][j][k] += 1ll * ntt[ii][j][k] * w_10s[ii * i % 10];
                }
            }
        }
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 16; k++) {
                ntt[i][j][k] = tmp[i][j][k] % Q;
            }
        }
    }
}

void intt10(int16_t ntt[10][9][16]) {
    memset(tmp, 0, sizeof(tmp));
    for (int i = 0; i < 10; i++) {
        for (int ii = 0; ii < 10; ii++) {
            for (int j = 0; j < 9; j++) {
                for (int k = 0; k < 16; k++) {
                    tmp[i][j][k] += 1ll * ntt[ii][j][k] * w_10s[(10 - ii) * i % 10];
                }
            }
        }
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 16; k++) {
                ntt[i][j][k] = tmp[i][j][k] % Q;
            }
        }
    }
}
void ntt9(int16_t ntt[10][9][16]) {
    memset(tmp, 0, sizeof(tmp));
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int jj = 0; jj < 9; jj++) {
                for (int k = 0; k < 16; k++) {
                    tmp[i][j][k] += 1ll * ntt[i][jj][k] * w_9s[jj * j % 9];
                }
            }
        }
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 16; k++) {
                ntt[i][j][k] = tmp[i][j][k] % Q;
            }
        }
    }
}
void intt9(int16_t ntt[10][9][16]) {
    memset(tmp, 0, sizeof(tmp));
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int jj = 0; jj < 9; jj++) {
                for (int k = 0; k < 16; k++) {
                    tmp[i][j][k] += 1ll * ntt[i][jj][k] * w_9s[(9 - jj) * j % 9];
                }
            }
        }
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 16; k++) {
                ntt[i][j][k] = tmp[i][j][k] % Q;
            }
        }
    }
}
void base_mul(int16_t in1_ntt[10][9][16], int16_t in2_ntt[10][9][16], int16_t out_ntt[10][9][16]) {
    memset(tmp, 0, sizeof(tmp));
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            int64_t twiddle = 1ll * w_10s[i] * w_9s[j] % Q;
            for (int k1 = 0; k1 < 16; k1++) {
                for (int k2 = 0; k2 < 16; k2++) {
                    if (k1 + k2 >= 16) {
                        tmp[i][j][k1 + k2 - 16] += twiddle * in1_ntt[i][j][k1] * in2_ntt[i][j][k2];

                    }
                    else {
                        tmp[i][j][k1 + k2] += 1ll * in1_ntt[i][j][k1] * in2_ntt[i][j][k2];
                    }
                }
            }
        }
    }
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 9; j++) {
            for (int k = 0; k < 16; k++) {
                out_ntt[i][j][k] = tmp[i][j][k] % Q;
            }
        }
    }
}

void forward(int16_t in_poly[P], int16_t out_ntt[10][9][16]) {
    memset(out_ntt, 0, NTT_SZ);
    for (int i = 0; i < P; i++) {
        int t = i >> 4;
        out_ntt[t % 10][t % 9][i % 16] = in_poly[i];
    }
    ntt10(out_ntt);
    ntt9(out_ntt);
}

const int32_t INV90 = 4540;
void backward(int16_t in_ntt[10][9][16], int16_t out_poly[P]) {
    memset(out_poly, 0, POLY_SZ);
    intt9(in_ntt);
    intt10(in_ntt);
    for (int i = 0; i < 1440; i++) {
        int t = i >> 4;
        out_poly[i] = in_ntt[t % 10][t % 9][i % 16] * INV90 % Q;
    }
}

void mult_low(int16_t in1_poly[81], int16_t in2_poly[81], int16_t out_poly[81]) {
    memset(out_poly, 0, sizeof(int16_t) * 81);
    for (int i = 0; i < 81; i++) {
        for (int j = 0; i + j < 81; j++) {
            out_poly[i + j] += 1ll * in1_poly[i] * in2_poly[j] % Q;
            out_poly[i + j] %= Q;
        }
    }
}

int16_t in1_ntt[10][9][16];
int16_t in2_ntt[10][9][16];
int16_t out_ntt[10][9][16];
int16_t out_low[81];
int16_t out_main[1440];
void mult(int16_t in1_poly[P], int16_t in2_poly[P], int16_t out_poly[P]) {
    forward(in1_poly, in1_ntt);
    forward(in2_poly, in2_ntt);
    mult_low(in1_poly, in2_poly, out_low);
    base_mul(in1_ntt, in2_ntt, out_ntt);
    backward(out_ntt, out_main);
    for (int i = 0; i < 1440; i++) {
        printf("%"PRId16" ", out_main[i]);
    }
    printf("\n");
    for (int i = 0; i < 760; i++) {
        out_poly[i] = out_main[(i + 760) % 1440] + out_main[(i + 761) % 1440];
    }
    out_poly[0] -= out_main[760];
    out_poly[760] += out_main[80];
    for (int i = 81; i < 761; i++) {
        out_poly[i] += out_main[i];
    }
    for (int i = 680; i < 760; i++) {
        out_poly[i] += 2 * Q - out_low[i - 680] - out_low[i - 679];
    }
    out_poly[679] += Q - out_low[0];
    out_poly[760] += Q - out_low[80];
    for (int i = 0; i < 80; i++) {
        out_poly[i] += out_low[i];
    }
    out_poly[80] += out_low[80];
    for (int i = 0; i < 760; i++) {
        out_poly[i] %= Q;
    }
    out_poly[760] %= Q;
}

int16_t a_poly[P];
int16_t b_poly[P];
int16_t c_poly[P];
int main() {
    init();

    scan_poly(a_poly);
    scan_poly(b_poly);
    mult(a_poly, b_poly, c_poly);
    print_poly(c_poly);
}
