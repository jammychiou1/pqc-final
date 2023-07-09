.align 6
.global __asm_weighted16x16
.global ___asm_weighted16x16
__asm_weighted16x16:
___asm_weighted16x16:

    sub sp, sp, #(16*4)
    stp  d8,  d9, [sp, #16*0]
    stp d10, d11, [sp, #16*1]
    stp d12, d13, [sp, #16*2]
    stp d14, d15, [sp, #16*3]

    mov w4, #-15631
    mov v0.H[0], w4
    mov w5, #4591
    mov v0.H[1], w5


    ldr w6, [x3], #4
    mov v0.S[1], w6

    ldp q8, q9, [x1], #32
    ldp q4, q5, [x2], #32

    sqrdmulh  v11.8H, v9.8H, v0.H[3]
    sqrdmulh  v10.8H, v8.8H, v0.H[3]
    mul       v13.8H, v9.8H, v0.H[2]
    mul       v12.8H, v8.8H, v0.H[2]
    mls       v13.8H, v11.8H, v0.H[1]
    mls       v12.8H, v10.8H, v0.H[1]

    ext       v17.16B, v13.16B, v8.16B, #14
    ext       v18.16B, v13.16B, v8.16B, #12
    ext       v19.16B, v13.16B, v8.16B, #10
    ext       v20.16B, v13.16B, v8.16B, #8
    ext       v21.16B, v13.16B, v8.16B, #6
    ext       v22.16B, v13.16B, v8.16B, #4
    ext       v23.16B, v13.16B, v8.16B, #2

    smull     v28.4S, v8.4H, v4.H[0]
    ext       v14.16B, v12.16B, v13.16B, #14
    smull2    v29.4S, v8.8H, v4.H[0]
    ext       v15.16B, v12.16B, v13.16B, #12
    smull     v30.4S, v8.4H, v5.H[0]
    ext       v16.16B, v12.16B, v13.16B, #10
    smull2    v31.4S, v8.8H, v5.H[0]
    ext       v24.16B, v12.16B, v13.16B, #8
    smlal     v28.4S, v17.4H, v4.H[1]
    ext       v25.16B, v12.16B, v13.16B, #6
    smlal2    v29.4S, v17.8H, v4.H[1]
    ext       v26.16B, v12.16B, v13.16B, #4
    smlal     v30.4S, v17.4H, v5.H[1]
    ext       v27.16B, v12.16B, v13.16B, #2
    smlal2    v31.4S, v17.8H, v5.H[1]
    ext       v17.16B, v8.16B, v9.16B, #14
    smlal     v28.4S, v18.4H, v4.H[2]
    smlal2    v29.4S, v18.8H, v4.H[2]
    smlal     v30.4S, v18.4H, v5.H[2]
    smlal2    v31.4S, v18.8H, v5.H[2]
    ext       v18.16B, v8.16B, v9.16B, #12
    smlal     v28.4S, v19.4H, v4.H[3]
    smlal2    v29.4S, v19.8H, v4.H[3]
    smlal     v30.4S, v19.4H, v5.H[3]
    smlal2    v31.4S, v19.8H, v5.H[3]
    ext       v19.16B, v8.16B, v9.16B, #10
    smlal     v28.4S, v20.4H, v4.H[4]
    smlal2    v29.4S, v20.8H, v4.H[4]
    smlal     v30.4S, v20.4H, v5.H[4]
    smlal2    v31.4S, v20.8H, v5.H[4]
    ext       v20.16B, v8.16B, v9.16B, #8
    smlal     v28.4S, v21.4H, v4.H[5]
    smlal2    v29.4S, v21.8H, v4.H[5]
    smlal     v30.4S, v21.4H, v5.H[5]
    smlal2    v31.4S, v21.8H, v5.H[5]
    ext       v21.16B, v8.16B, v9.16B, #6
    smlal     v28.4S, v22.4H, v4.H[6]
    smlal2    v29.4S, v22.8H, v4.H[6]
    smlal     v30.4S, v22.4H, v5.H[6]
    smlal2    v31.4S, v22.8H, v5.H[6]
    ext       v22.16B, v8.16B, v9.16B, #4
    smlal     v28.4S, v23.4H, v4.H[7]
    smlal2    v29.4S, v23.8H, v4.H[7]
    smlal     v30.4S, v23.4H, v5.H[7]
    smlal2    v31.4S, v23.8H, v5.H[7]
    ext       v23.16B, v8.16B, v9.16B, #2

    ldr w6, [x3], #4
    mov v0.S[1], w6

    smlal     v28.4S, v13.4H, v5.H[0]
    smlal2    v29.4S, v13.8H, v5.H[0]
    smlal     v30.4S, v9.4H, v4.H[0]
    smlal2    v31.4S, v9.8H, v4.H[0]
    ldp q8, q9, [x1], #32
    smlal     v28.4S, v14.4H, v5.H[1]
    smlal2    v29.4S, v14.8H, v5.H[1]
    smlal     v30.4S, v17.4H, v4.H[1]
    smlal2    v31.4S, v17.8H, v4.H[1]
    smlal     v28.4S, v15.4H, v5.H[2]
    smlal2    v29.4S, v15.8H, v5.H[2]
    smlal     v30.4S, v18.4H, v4.H[2]
    smlal2    v31.4S, v18.8H, v4.H[2]
    smlal     v28.4S, v16.4H, v5.H[3]
    smlal2    v29.4S, v16.8H, v5.H[3]
    smlal     v30.4S, v19.4H, v4.H[3]
    smlal2    v31.4S, v19.8H, v4.H[3]
    smlal     v28.4S, v24.4H, v5.H[4]
    smlal2    v29.4S, v24.8H, v5.H[4]
    smlal     v30.4S, v20.4H, v4.H[4]
    smlal2    v31.4S, v20.8H, v4.H[4]
    smlal     v28.4S, v25.4H, v5.H[5]
    smlal2    v29.4S, v25.8H, v5.H[5]
    smlal     v30.4S, v21.4H, v4.H[5]
    smlal2    v31.4S, v21.8H, v4.H[5]
    smlal     v28.4S, v26.4H, v5.H[6]
    smlal2    v29.4S, v26.8H, v5.H[6]
    smlal     v30.4S, v22.4H, v4.H[6]
    smlal2    v31.4S, v22.8H, v4.H[6]
    smlal     v28.4S, v27.4H, v5.H[7]
    smlal2    v29.4S, v27.8H, v5.H[7]
    smlal     v30.4S, v23.4H, v4.H[7]
    smlal2    v31.4S, v23.8H, v4.H[7]

    ldp q4, q5, [x2], #32

    sqrdmulh  v11.8H, v9.8H, v0.H[3]
    uzp1      v14.8H, v28.8H, v29.8H
    sqrdmulh  v10.8H, v8.8H, v0.H[3]
    mul       v14.8H, v14.8H, v0.H[0]
    mul       v13.8H, v9.8H, v0.H[2]
    mul       v12.8H, v8.8H, v0.H[2]
    smlal     v28.4S, v14.4H, v0.H[1]
    smlal2    v29.4S, v14.8H, v0.H[1]
    mls       v13.8H, v11.8H, v0.H[1]
    mls       v12.8H, v10.8H, v0.H[1]
    uzp2      v14.8H, v28.8H, v29.8H

    str q14, [x0], #16

    uzp1      v15.8H, v30.8H, v31.8H
    ext       v17.16B, v13.16B, v8.16B, #14
    ext       v18.16B, v13.16B, v8.16B, #12
    mul       v15.8H, v15.8H, v0.H[0]
    ext       v19.16B, v13.16B, v8.16B, #10
    ext       v20.16B, v13.16B, v8.16B, #8
    smlal     v30.4S, v15.4H, v0.H[1]
    ext       v21.16B, v13.16B, v8.16B, #6
    smlal2    v31.4S, v15.8H, v0.H[1]
    ext       v22.16B, v13.16B, v8.16B, #4
    ext       v23.16B, v13.16B, v8.16B, #2
    uzp2      v15.8H, v30.8H, v31.8H

    str q15, [x0], #16

    mov w7, #50

    __asm_weighted16x16_loop:

    smull     v28.4S, v8.4H, v4.H[0]
    ext       v14.16B, v12.16B, v13.16B, #14
    smull2    v29.4S, v8.8H, v4.H[0]
    ext       v15.16B, v12.16B, v13.16B, #12
    smull     v30.4S, v8.4H, v5.H[0]
    ext       v16.16B, v12.16B, v13.16B, #10
    smull2    v31.4S, v8.8H, v5.H[0]
    ext       v24.16B, v12.16B, v13.16B, #8
    smlal     v28.4S, v17.4H, v4.H[1]
    ext       v25.16B, v12.16B, v13.16B, #6
    smlal2    v29.4S, v17.8H, v4.H[1]
    ext       v26.16B, v12.16B, v13.16B, #4
    smlal     v30.4S, v17.4H, v5.H[1]
    ext       v27.16B, v12.16B, v13.16B, #2
    smlal2    v31.4S, v17.8H, v5.H[1]
    ext       v17.16B, v8.16B, v9.16B, #14
    smlal     v28.4S, v18.4H, v4.H[2]
    smlal2    v29.4S, v18.8H, v4.H[2]
    smlal     v30.4S, v18.4H, v5.H[2]
    smlal2    v31.4S, v18.8H, v5.H[2]
    ext       v18.16B, v8.16B, v9.16B, #12
    smlal     v28.4S, v19.4H, v4.H[3]
    smlal2    v29.4S, v19.8H, v4.H[3]
    smlal     v30.4S, v19.4H, v5.H[3]
    smlal2    v31.4S, v19.8H, v5.H[3]
    ext       v19.16B, v8.16B, v9.16B, #10
    smlal     v28.4S, v20.4H, v4.H[4]
    smlal2    v29.4S, v20.8H, v4.H[4]
    smlal     v30.4S, v20.4H, v5.H[4]
    smlal2    v31.4S, v20.8H, v5.H[4]
    ext       v20.16B, v8.16B, v9.16B, #8
    smlal     v28.4S, v21.4H, v4.H[5]
    smlal2    v29.4S, v21.8H, v4.H[5]
    smlal     v30.4S, v21.4H, v5.H[5]
    smlal2    v31.4S, v21.8H, v5.H[5]
    ext       v21.16B, v8.16B, v9.16B, #6
    smlal     v28.4S, v22.4H, v4.H[6]
    smlal2    v29.4S, v22.8H, v4.H[6]
    smlal     v30.4S, v22.4H, v5.H[6]
    smlal2    v31.4S, v22.8H, v5.H[6]
    ext       v22.16B, v8.16B, v9.16B, #4
    smlal     v28.4S, v23.4H, v4.H[7]
    smlal2    v29.4S, v23.8H, v4.H[7]
    smlal     v30.4S, v23.4H, v5.H[7]
    smlal2    v31.4S, v23.8H, v5.H[7]
    ext       v23.16B, v8.16B, v9.16B, #2

    ldr w6, [x3], #4
    mov v0.S[1], w6

    smlal     v28.4S, v13.4H, v5.H[0]
    smlal2    v29.4S, v13.8H, v5.H[0]
    smlal     v30.4S, v9.4H, v4.H[0]
    smlal2    v31.4S, v9.8H, v4.H[0]
    ldp q8, q9, [x1], #32
    smlal     v28.4S, v14.4H, v5.H[1]
    smlal2    v29.4S, v14.8H, v5.H[1]
    smlal     v30.4S, v17.4H, v4.H[1]
    smlal2    v31.4S, v17.8H, v4.H[1]
    smlal     v28.4S, v15.4H, v5.H[2]
    smlal2    v29.4S, v15.8H, v5.H[2]
    smlal     v30.4S, v18.4H, v4.H[2]
    smlal2    v31.4S, v18.8H, v4.H[2]
    smlal     v28.4S, v16.4H, v5.H[3]
    smlal2    v29.4S, v16.8H, v5.H[3]
    smlal     v30.4S, v19.4H, v4.H[3]
    smlal2    v31.4S, v19.8H, v4.H[3]
    smlal     v28.4S, v24.4H, v5.H[4]
    smlal2    v29.4S, v24.8H, v5.H[4]
    smlal     v30.4S, v20.4H, v4.H[4]
    smlal2    v31.4S, v20.8H, v4.H[4]
    smlal     v28.4S, v25.4H, v5.H[5]
    smlal2    v29.4S, v25.8H, v5.H[5]
    smlal     v30.4S, v21.4H, v4.H[5]
    smlal2    v31.4S, v21.8H, v4.H[5]
    smlal     v28.4S, v26.4H, v5.H[6]
    smlal2    v29.4S, v26.8H, v5.H[6]
    smlal     v30.4S, v22.4H, v4.H[6]
    smlal2    v31.4S, v22.8H, v4.H[6]
    smlal     v28.4S, v27.4H, v5.H[7]
    smlal2    v29.4S, v27.8H, v5.H[7]
    smlal     v30.4S, v23.4H, v4.H[7]
    smlal2    v31.4S, v23.8H, v4.H[7]

    ldp q4, q5, [x2], #32

    sqrdmulh  v11.8H, v9.8H, v0.H[3]
    uzp1      v14.8H, v28.8H, v29.8H
    sqrdmulh  v10.8H, v8.8H, v0.H[3]
    mul       v14.8H, v14.8H, v0.H[0]
    mul       v13.8H, v9.8H, v0.H[2]
    mul       v12.8H, v8.8H, v0.H[2]
    smlal     v28.4S, v14.4H, v0.H[1]
    smlal2    v29.4S, v14.8H, v0.H[1]
    mls       v13.8H, v11.8H, v0.H[1]
    mls       v12.8H, v10.8H, v0.H[1]
    uzp2      v14.8H, v28.8H, v29.8H

    str q14, [x0], #16

    uzp1      v15.8H, v30.8H, v31.8H
    ext       v17.16B, v13.16B, v8.16B, #14
    ext       v18.16B, v13.16B, v8.16B, #12
    mul       v15.8H, v15.8H, v0.H[0]
    ext       v19.16B, v13.16B, v8.16B, #10
    ext       v20.16B, v13.16B, v8.16B, #8
    smlal     v30.4S, v15.4H, v0.H[1]
    ext       v21.16B, v13.16B, v8.16B, #6
    smlal2    v31.4S, v15.8H, v0.H[1]
    ext       v22.16B, v13.16B, v8.16B, #4
    ext       v23.16B, v13.16B, v8.16B, #2
    uzp2      v15.8H, v30.8H, v31.8H

    str q15, [x0], #16

    smull     v28.4S, v8.4H, v4.H[0]
    ext       v14.16B, v12.16B, v13.16B, #14
    smull2    v29.4S, v8.8H, v4.H[0]
    ext       v15.16B, v12.16B, v13.16B, #12
    smull     v30.4S, v8.4H, v5.H[0]
    ext       v16.16B, v12.16B, v13.16B, #10
    smull2    v31.4S, v8.8H, v5.H[0]
    ext       v24.16B, v12.16B, v13.16B, #8
    smlal     v28.4S, v17.4H, v4.H[1]
    ext       v25.16B, v12.16B, v13.16B, #6
    smlal2    v29.4S, v17.8H, v4.H[1]
    ext       v26.16B, v12.16B, v13.16B, #4
    smlal     v30.4S, v17.4H, v5.H[1]
    ext       v27.16B, v12.16B, v13.16B, #2
    smlal2    v31.4S, v17.8H, v5.H[1]
    ext       v17.16B, v8.16B, v9.16B, #14
    smlal     v28.4S, v18.4H, v4.H[2]
    smlal2    v29.4S, v18.8H, v4.H[2]
    smlal     v30.4S, v18.4H, v5.H[2]
    smlal2    v31.4S, v18.8H, v5.H[2]
    ext       v18.16B, v8.16B, v9.16B, #12
    smlal     v28.4S, v19.4H, v4.H[3]
    smlal2    v29.4S, v19.8H, v4.H[3]
    smlal     v30.4S, v19.4H, v5.H[3]
    smlal2    v31.4S, v19.8H, v5.H[3]
    ext       v19.16B, v8.16B, v9.16B, #10
    smlal     v28.4S, v20.4H, v4.H[4]
    smlal2    v29.4S, v20.8H, v4.H[4]
    smlal     v30.4S, v20.4H, v5.H[4]
    smlal2    v31.4S, v20.8H, v5.H[4]
    ext       v20.16B, v8.16B, v9.16B, #8
    smlal     v28.4S, v21.4H, v4.H[5]
    smlal2    v29.4S, v21.8H, v4.H[5]
    smlal     v30.4S, v21.4H, v5.H[5]
    smlal2    v31.4S, v21.8H, v5.H[5]
    ext       v21.16B, v8.16B, v9.16B, #6
    smlal     v28.4S, v22.4H, v4.H[6]
    smlal2    v29.4S, v22.8H, v4.H[6]
    smlal     v30.4S, v22.4H, v5.H[6]
    smlal2    v31.4S, v22.8H, v5.H[6]
    ext       v22.16B, v8.16B, v9.16B, #4
    smlal     v28.4S, v23.4H, v4.H[7]
    smlal2    v29.4S, v23.8H, v4.H[7]
    smlal     v30.4S, v23.4H, v5.H[7]
    smlal2    v31.4S, v23.8H, v5.H[7]
    ext       v23.16B, v8.16B, v9.16B, #2

    ldr w6, [x3], #4
    mov v0.S[1], w6

    smlal     v28.4S, v13.4H, v5.H[0]
    smlal2    v29.4S, v13.8H, v5.H[0]
    smlal     v30.4S, v9.4H, v4.H[0]
    smlal2    v31.4S, v9.8H, v4.H[0]
    ldp q8, q9, [x1], #32
    smlal     v28.4S, v14.4H, v5.H[1]
    smlal2    v29.4S, v14.8H, v5.H[1]
    smlal     v30.4S, v17.4H, v4.H[1]
    smlal2    v31.4S, v17.8H, v4.H[1]
    smlal     v28.4S, v15.4H, v5.H[2]
    smlal2    v29.4S, v15.8H, v5.H[2]
    smlal     v30.4S, v18.4H, v4.H[2]
    smlal2    v31.4S, v18.8H, v4.H[2]
    smlal     v28.4S, v16.4H, v5.H[3]
    smlal2    v29.4S, v16.8H, v5.H[3]
    smlal     v30.4S, v19.4H, v4.H[3]
    smlal2    v31.4S, v19.8H, v4.H[3]
    smlal     v28.4S, v24.4H, v5.H[4]
    smlal2    v29.4S, v24.8H, v5.H[4]
    smlal     v30.4S, v20.4H, v4.H[4]
    smlal2    v31.4S, v20.8H, v4.H[4]
    smlal     v28.4S, v25.4H, v5.H[5]
    smlal2    v29.4S, v25.8H, v5.H[5]
    smlal     v30.4S, v21.4H, v4.H[5]
    smlal2    v31.4S, v21.8H, v4.H[5]
    smlal     v28.4S, v26.4H, v5.H[6]
    smlal2    v29.4S, v26.8H, v5.H[6]
    smlal     v30.4S, v22.4H, v4.H[6]
    smlal2    v31.4S, v22.8H, v4.H[6]
    smlal     v28.4S, v27.4H, v5.H[7]
    smlal2    v29.4S, v27.8H, v5.H[7]
    smlal     v30.4S, v23.4H, v4.H[7]
    smlal2    v31.4S, v23.8H, v4.H[7]

    ldp q4, q5, [x2], #32

    sqrdmulh  v11.8H, v9.8H, v0.H[3]
    uzp1      v14.8H, v28.8H, v29.8H
    sqrdmulh  v10.8H, v8.8H, v0.H[3]
    mul       v14.8H, v14.8H, v0.H[0]
    mul       v13.8H, v9.8H, v0.H[2]
    mul       v12.8H, v8.8H, v0.H[2]
    smlal     v28.4S, v14.4H, v0.H[1]
    smlal2    v29.4S, v14.8H, v0.H[1]
    mls       v13.8H, v11.8H, v0.H[1]
    mls       v12.8H, v10.8H, v0.H[1]
    uzp2      v14.8H, v28.8H, v29.8H

    str q14, [x0], #16

    uzp1      v15.8H, v30.8H, v31.8H
    ext       v17.16B, v13.16B, v8.16B, #14
    ext       v18.16B, v13.16B, v8.16B, #12
    mul       v15.8H, v15.8H, v0.H[0]
    ext       v19.16B, v13.16B, v8.16B, #10
    ext       v20.16B, v13.16B, v8.16B, #8
    smlal     v30.4S, v15.4H, v0.H[1]
    ext       v21.16B, v13.16B, v8.16B, #6
    smlal2    v31.4S, v15.8H, v0.H[1]
    ext       v22.16B, v13.16B, v8.16B, #4
    ext       v23.16B, v13.16B, v8.16B, #2
    uzp2      v15.8H, v30.8H, v31.8H

    str q15, [x0], #16

    sub w7, w7, #1
    cbnz w7, __asm_weighted16x16_loop

    smull     v28.4S, v8.4H, v4.H[0]
    ext       v14.16B, v12.16B, v13.16B, #14
    smull2    v29.4S, v8.8H, v4.H[0]
    ext       v15.16B, v12.16B, v13.16B, #12
    smull     v30.4S, v8.4H, v5.H[0]
    ext       v16.16B, v12.16B, v13.16B, #10
    smull2    v31.4S, v8.8H, v5.H[0]
    ext       v24.16B, v12.16B, v13.16B, #8
    smlal     v28.4S, v17.4H, v4.H[1]
    ext       v25.16B, v12.16B, v13.16B, #6
    smlal2    v29.4S, v17.8H, v4.H[1]
    ext       v26.16B, v12.16B, v13.16B, #4
    smlal     v30.4S, v17.4H, v5.H[1]
    ext       v27.16B, v12.16B, v13.16B, #2
    smlal2    v31.4S, v17.8H, v5.H[1]
    ext       v17.16B, v8.16B, v9.16B, #14
    smlal     v28.4S, v18.4H, v4.H[2]
    smlal2    v29.4S, v18.8H, v4.H[2]
    smlal     v30.4S, v18.4H, v5.H[2]
    smlal2    v31.4S, v18.8H, v5.H[2]
    ext       v18.16B, v8.16B, v9.16B, #12
    smlal     v28.4S, v19.4H, v4.H[3]
    smlal2    v29.4S, v19.8H, v4.H[3]
    smlal     v30.4S, v19.4H, v5.H[3]
    smlal2    v31.4S, v19.8H, v5.H[3]
    ext       v19.16B, v8.16B, v9.16B, #10
    smlal     v28.4S, v20.4H, v4.H[4]
    smlal2    v29.4S, v20.8H, v4.H[4]
    smlal     v30.4S, v20.4H, v5.H[4]
    smlal2    v31.4S, v20.8H, v5.H[4]
    ext       v20.16B, v8.16B, v9.16B, #8
    smlal     v28.4S, v21.4H, v4.H[5]
    smlal2    v29.4S, v21.8H, v4.H[5]
    smlal     v30.4S, v21.4H, v5.H[5]
    smlal2    v31.4S, v21.8H, v5.H[5]
    ext       v21.16B, v8.16B, v9.16B, #6
    smlal     v28.4S, v22.4H, v4.H[6]
    smlal2    v29.4S, v22.8H, v4.H[6]
    smlal     v30.4S, v22.4H, v5.H[6]
    smlal2    v31.4S, v22.8H, v5.H[6]
    ext       v22.16B, v8.16B, v9.16B, #4
    smlal     v28.4S, v23.4H, v4.H[7]
    smlal2    v29.4S, v23.8H, v4.H[7]
    smlal     v30.4S, v23.4H, v5.H[7]
    smlal2    v31.4S, v23.8H, v5.H[7]
    ext       v23.16B, v8.16B, v9.16B, #2

    smlal     v28.4S, v13.4H, v5.H[0]
    smlal2    v29.4S, v13.8H, v5.H[0]
    smlal     v30.4S, v9.4H, v4.H[0]
    smlal2    v31.4S, v9.8H, v4.H[0]
    smlal     v28.4S, v14.4H, v5.H[1]
    smlal2    v29.4S, v14.8H, v5.H[1]
    smlal     v30.4S, v17.4H, v4.H[1]
    smlal2    v31.4S, v17.8H, v4.H[1]
    smlal     v28.4S, v15.4H, v5.H[2]
    smlal2    v29.4S, v15.8H, v5.H[2]
    smlal     v30.4S, v18.4H, v4.H[2]
    smlal2    v31.4S, v18.8H, v4.H[2]
    smlal     v28.4S, v16.4H, v5.H[3]
    smlal2    v29.4S, v16.8H, v5.H[3]
    smlal     v30.4S, v19.4H, v4.H[3]
    smlal2    v31.4S, v19.8H, v4.H[3]
    smlal     v28.4S, v24.4H, v5.H[4]
    smlal2    v29.4S, v24.8H, v5.H[4]
    smlal     v30.4S, v20.4H, v4.H[4]
    smlal2    v31.4S, v20.8H, v4.H[4]
    smlal     v28.4S, v25.4H, v5.H[5]
    smlal2    v29.4S, v25.8H, v5.H[5]
    smlal     v30.4S, v21.4H, v4.H[5]
    smlal2    v31.4S, v21.8H, v4.H[5]
    smlal     v28.4S, v26.4H, v5.H[6]
    smlal2    v29.4S, v26.8H, v5.H[6]
    smlal     v30.4S, v22.4H, v4.H[6]
    smlal2    v31.4S, v22.8H, v4.H[6]
    smlal     v28.4S, v27.4H, v5.H[7]
    smlal2    v29.4S, v27.8H, v5.H[7]
    smlal     v30.4S, v23.4H, v4.H[7]
    smlal2    v31.4S, v23.8H, v4.H[7]

    uzp1      v14.8H, v28.8H, v29.8H
    mul       v14.8H, v14.8H, v0.H[0]
    smlal     v28.4S, v14.4H, v0.H[1]
    smlal2    v29.4S, v14.8H, v0.H[1]
    uzp2      v14.8H, v28.8H, v29.8H

    str q14, [x0], #16

    uzp1      v15.8H, v30.8H, v31.8H
    mul       v15.8H, v15.8H, v0.H[0]
    smlal     v30.4S, v15.4H, v0.H[1]
    smlal2    v31.4S, v15.8H, v0.H[1]
    uzp2      v15.8H, v30.8H, v31.8H

    str q15, [x0], #16

    ldp  d8,  d9, [sp, #16*0]
    ldp d10, d11, [sp, #16*1]
    ldp d12, d13, [sp, #16*2]
    ldp d14, d15, [sp, #16*3]
    add sp, sp, #(16*4)

    br lr
