#ifndef REDUCE_TBL_H
#define REDUCE_TBL_H

inline uint8x16_t pad_zip(uint8x8_t a, uint8x8_t b) {
    uint8x16_t res;
    asm ("zip1 %[vd].16b, %[v1].16b, %[v2].16b"
      : [vd] "=w" (res)
      : [v1] "w" (a), [v2] "w" (b));
    return res;
}

constexpr static uint8_t TBL_LOW[16] = {
    0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
    0x89, 0x9a, 0xab, 0xbc, 0xcd, 0xde, 0xef, 0x00,
};
constexpr static uint8_t TBL_HIGH[16] = {
    0x00, 0xee, 0xdc, 0xca, 0xb8, 0xa6, 0x94, 0x82,
    0x7d, 0x6b, 0x59, 0x47, 0x35, 0x23, 0x11, 0x00,
};

inline void reduce_tbl(int16x8_t &v) {

    uint8x16_t tbl_low = vld1q_u8(&TBL_LOW[0]);
    uint8x16_t tbl_high = vld1q_u8(&TBL_HIGH[0]);
    uint8x8_t idx = vmovn_u16(vshrq_n_u16(vreinterpretq_u16_s16(v), 12));
    uint8x8_t low = vqtbl1_u8(tbl_low, idx);
    uint8x8_t high = vqtbl1_u8(tbl_high, idx);
    int16x8_t comb = vreinterpretq_s16_u8(pad_zip(low, high));
    v = vaddq_s16(v, comb);
}

inline void reduce_tbl_x2(int16x8_t &v1, int16x8_t &v2) {
    uint8x16_t tbl_low = vld1q_u8(&TBL_LOW[0]);
    uint8x16_t tbl_high = vld1q_u8(&TBL_HIGH[0]);
    uint8x16_t idx = vshrq_n_u8(vuzp2q_u8(vreinterpretq_u8_s16(v1), vreinterpretq_u8_s16(v2)), 4);
    uint8x16_t low = vqtbl1q_u8(tbl_low, idx);
    uint8x16_t high = vqtbl1q_u8(tbl_high, idx);
    int16x8_t comb1 = vreinterpretq_s16_u8(vzip1q_u8(low, high));
    int16x8_t comb2 = vreinterpretq_s16_u8(vzip2q_u8(low, high));
    v1 = vaddq_s16(v1, comb1);
    v2 = vaddq_s16(v2, comb2);
}

#endif // REDUCE_TBL_H
