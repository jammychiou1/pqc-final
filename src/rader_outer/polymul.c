#include "rader_outer/polymul.h"

#include <arm_neon.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "rader_outer/ntt_17.h"
#include "rader_outer/ntt_3_2.h"
#include "rader_outer/intt_17.h"
#include "rader_outer/intt_3_2.h"
#include "rader_outer/base_mul.h"

static int p = 761;
static int q = 4591;

int16_t Fq_freeze(int16_t x) {

    int64_t t;
    t = (int64_t)x * 935519;
    t = (t + (1LL << 31)) >> 32;

    return x - (int32_t)t * q;

}

// h = fg mod (Z_q, x^p - x - 1)
void Rq_reduce(int16_t *h, int16_t *fg) {
    int i;
    int16x8_t t0, t1;

    for (i = p; i < p + p - 1 - 8; i += 8) {
        t0 = vld1q_s16(fg + i);
        t1 = vld1q_s16(fg + i - p + 1);
        t0 = t0 + t1;
        t1 = vqdmulhq_n_s16(t0, 29235);
        t1 = vrshrq_n_s16(t1, 12);
        t0 = vmlsq_n_s16(t0, t1, q);
        vst1q_s16(fg + i - p + 1, t0);
    }
    for (; i < p + p - 1; i++) {
        fg[i - p + 1] = Fq_freeze(fg[i - p + 1] + fg[i]);
    }
    for (i = p; i < p + p - 1 - 8; i += 8) {
        t0 = vld1q_s16(fg + i);
        t1 = vld1q_s16(fg + i - p);
        t0 = t0 + t1;
        t1 = vqdmulhq_n_s16(t0, 29235);
        t1 = vrshrq_n_s16(t1, 12);
        t0 = vmlsq_n_s16(t0, t1, q);
        vst1q_s16(fg + i - p, t0);
    }
    for (; i < p + p - 1; i++) {
        fg[i - p] = Fq_freeze(fg[i - p] + fg[i]);
    }

    for (i = 0; i < p; i++) {
        h[i] = fg[i];
    }
}

static int16_t poly1_NTT[51 * 32], poly2_NTT[51 * 32];
static int16_t out_tmp[1632] = {};
void polymul(int16_t *des, const int16_t *src1, const int16_t *src2) {
    int16_t *res_NTT = poly1_NTT;

    ntt_17(poly1_NTT, src1);
    ntt_17(poly2_NTT, src2);

    ntt_3_2(poly1_NTT);
    ntt_3_2(poly2_NTT);

    base_mul(poly1_NTT, poly2_NTT, res_NTT);

    intt_3_2(res_NTT);

    intt_17(res_NTT, out_tmp);

    Rq_reduce(des, out_tmp);
}

// int16_t in1_poly[1632] = {};
// int16_t in2_poly[1632] = {};
// int16_t out_poly[1632] = {};
// int16_t out_red[1632] = {};

// int main() {
//     for (int i = 0; i < 100000; i++) {
//         polymul(out_poly, in1_poly, in2_poly);
//         Rq_reduce(out_red, out_poly);
//     }
// }

// int main() {
//     for (int i = 0; i < 761; i++) {
//         scanf("%hd\n", &in1_poly[i]);
//     }
//     for (int i = 0; i < 761; i++) {
//         scanf("%hd\n", &in2_poly[i]);
//     }
//     polymul(out_poly, in1_poly, in2_poly);
//     Rq_reduce(out_red, out_poly);
//     for (int i = 0; i < 761; i++) {
//         printf("%hd%c", out_red[i], " \n"[i == 760]);
//     }
// }

