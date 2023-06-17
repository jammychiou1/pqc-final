#include <stdio.h>
#include <string.h>

const int Q = 4591;
const int P = 761;
const int ARR_SZ = sizeof(int) * P;

void print_poly(int p[P]) {
    for (int i = 0; i < P; i++) {
        printf("%d%c", p[i], " \n"[i == P - 1]);
    }
}

void mul1(int a[P], int b[P], int c1[P]) {
    for (int i = 0; i < P; i++) {
        for (int j = 0; i + j < P; j++) {
            c1[i + j] += a[i] * b[j];
            c1[i + j] %= Q;
        }
    }
}

void mul2(int a[P], int b[P], int c2[P]) {
    a[0] = (a[0] + a[760]) % Q;
    b[0] = (b[0] + b[760]) % Q;
    for (int i = 0; i < P - 1; i++) {
        for (int j = 0; j < P - 1; j++) {
            c2[(i + j) % (P - 1)] += a[i] * b[j];
            c2[(i + j) % (P - 1)] %= Q;
        }
    }
}

void mul(int a[P], int b[P], int c[P]) {
    int c1[P], c2[P];
    memset(c, 0, ARR_SZ);
    memset(c1, 0, ARR_SZ);
    memset(c2, 0, ARR_SZ);
    mul1(a, b, c1);
    // print_poly(c1);
    mul2(a, b, c2);
    // print_poly(c2);
    c[0] = (c1[0] - c1[1] + c2[1] + Q) % Q;
    for (int i = 1; i < P - 2; i++) {
        c[i] = (-c1[i + 1] + c2[i] + c2[i + 1] + Q) % Q;
    }
    c[P - 2] = (-c1[0] - c1[P - 1] + c2[0] + c2[P - 2] + 2 * Q) % Q;
    c[P - 1] = (-c1[0] + c2[0] + Q) % Q;
}

int main() {
    int a[P], b[P], c[P];
    for (int i = 0; i < P; i++) {
        scanf("%d ", &a[i]);
    }
    for (int i = 0; i < P; i++) {
        scanf("%d ", &b[i]);
    }
    mul(a, b, c);
    print_poly(c);
}
