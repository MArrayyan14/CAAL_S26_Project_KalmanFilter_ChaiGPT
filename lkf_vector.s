# lkf_vector.s — Fully Vectorised Linear Kalman Filter (RISC-V RVV)
# Team ChaiGPT — Milestone 4

.option arch, +v

.equ WS_xtmp,      0
.equ WS_FP,        2208
.equ WS_FPFt,      611616
.equ WS_Ft,        1221024
.equ WS_Ht,        1830432
.equ WS_PHt,       1982784
.equ WS_PHtT,      2135136
.equ WS_Kt,        2287488
.equ WS_K,         2439840
.equ WS_S,         2592192
.equ WS_L,         2630280
.equ WS_y,         2668368
.equ WS_Hx,        2668920
.equ WS_Ky,        2669472
.equ WS_KH,        2671680
.equ WS_IKH,       3281088
.equ WS_IKHT,      3890496
.equ WS_IKH_P,     4499904
.equ WS_KR,        5109312
.equ WS_KTJ,       5261664
.equ WS_KRKt,      5414016
.equ LKF_WS_TOTAL, 6023424

.equ NN, 76176      # 276*276
.equ NM, 19044      # 276*69
.equ MM, 4761       # 69*69


.section .rodata
.align 3
ZERO_D: .double 0.0     # Used by cholesky_decomposition

.section .text

# run_lkf_asm 

.globl run_lkf_asm
.align 2
run_lkf_asm:
    addi sp,sp,-128
    sd ra,120(sp); sd s0,112(sp); sd s1,104(sp); sd s2,96(sp)
    sd s3,88(sp);  sd s4,80(sp);  sd s5,72(sp);  sd s6,64(sp)
    sd s7,56(sp);  sd s8,48(sp);  sd s9,40(sp);  sd s10,32(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    mv s4,a4; mv s5,a5; mv s6,a6

    li a0,LKF_WS_TOTAL; call malloc; mv s7,a0  # workspace
    li a0,2208;          call malloc; mv s8,a0  # x  (276 doubles)
    li a0,609408;        call malloc; mv s9,a0  # P  (276×276)

    # Zero x using RVV
    mv a0,s8; li a1,276; call vlkf_zero

    # Initialise x positions from first measurement frame
    li t0,0
.Lrlkf_ix:
    li t1,23; bge t0,t1,.Lrlkf_ix_d
    li t2,24; mul t3,t0,t2; add t3,s0,t3   # &meas[j*3]
    li t2,96; mul t4,t0,t2; add t4,s8,t4   # &x[j*12]
    fld ft0,0(t3);  fsd ft0,0(t4)          # px → x[j*12+0]
    fld ft1,8(t3);  fsd ft1,32(t4)         # py → x[j*12+4]
    fld ft2,16(t3); fsd ft2,64(t4)         # pz → x[j*12+8]
    addi t0,t0,1; j .Lrlkf_ix
.Lrlkf_ix_d:

    # P = I_276  (zero via RVV, diagonal set scalar)
    mv a0,s9; li a1,276; call init_identity

    # Frame loop
    li s10,0
.Lrlkf_loop:
    bge s10,s1,.Lrlkf_done
    li t0,552; mul t0,s10,t0; add t0,s0,t0   # z = &meas[k*69]
    mv a0,s8; mv a1,s9; mv a2,t0
    mv a3,s2; mv a4,s3; mv a5,s4; mv a6,s5; mv a7,s7
    call lkf_step
    # Store x into all_states[k*276..]
    li t0,2208; mul t0,s10,t0; add t0,s6,t0
    mv a0,t0; mv a1,s8; li a2,276; call vlkf_copy
    addi s10,s10,1; j .Lrlkf_loop

.Lrlkf_done:
    mv a0,s7; call free
    mv a0,s8; call free
    mv a0,s9; call free
    ld ra,120(sp); ld s0,112(sp); ld s1,104(sp); ld s2,96(sp)
    ld s3,88(sp);  ld s4,80(sp);  ld s5,72(sp);  ld s6,64(sp)
    ld s7,56(sp);  ld s8,48(sp);  ld s9,40(sp);  ld s10,32(sp)
    addi sp,sp,128; ret


# lkf_step 

.globl lkf_step
.align 2
lkf_step:
    addi sp,sp,-112
    sd ra,104(sp); sd s0,96(sp);  sd s1,88(sp);  sd s2,80(sp)
    sd s3,72(sp);  sd s4,64(sp);  sd s5,56(sp);  sd s6,48(sp)
    sd s7,40(sp);  sd s8,32(sp);  sd s9,24(sp);  sd s10,16(sp); sd s11,8(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    mv s4,a4; mv s5,a5; mv s6,a6; mv s7,a7

    # PREDICT 
    # x_tmp = F*x  → copy back to x
    li t0,WS_xtmp; add s8,s7,t0
    mv a0,s8; mv a1,s3; mv a2,s0; li a3,276; li a4,276; li a5,1
    call vlkf_matmul
    mv a0,s0; mv a1,s8; li a2,276; call vlkf_copy

    # FP = F*P
    li t0,WS_FP; add s8,s7,t0
    mv a0,s8; mv a1,s3; mv a2,s1; li a3,276; li a4,276; li a5,276
    call vlkf_matmul

    # Ft = F^T  (276×276)
    li t0,WS_Ft; add s9,s7,t0
    mv a0,s9; mv a1,s3; li a2,276; li a3,276
    call vlkf_transpose

    # FPFt = FP * Ft
    li t0,WS_FPFt; add s10,s7,t0
    mv a0,s10; mv a1,s8; mv a2,s9; li a3,276; li a4,276; li a5,276
    call vlkf_matmul

    # P = FPFt + Q
    mv a0,s1; mv a1,s10; mv a2,s4; li a3,NN
    call vlkf_add

    # UPDATE 
    # Hx = H*x 
    li t0,WS_Hx; add s8,s7,t0
    mv a0,s8; mv a1,s5; mv a2,s0; li a3,69; li a4,276; li a5,1
    call vlkf_matmul

    # y = z - Hx  
    li t0,WS_y; add s9,s7,t0
    mv a0,s9; mv a1,s2; mv a2,s8; li a3,69
    call vlkf_sub

    # Ht = H^T  
    li t0,WS_Ht; add s8,s7,t0
    mv a0,s8; mv a1,s5; li a2,69; li a3,276
    call vlkf_transpose

    # PHt = P * Ht  
    li t0,WS_PHt; add s10,s7,t0
    mv a0,s10; mv a1,s1; mv a2,s8; li a3,276; li a4,276; li a5,69
    call vlkf_matmul

    # S = H*PHt + R  
    li t0,WS_S; add s8,s7,t0
    mv a0,s8; mv a1,s5; mv a2,s10; li a3,69; li a4,276; li a5,69
    call vlkf_matmul
    mv a0,s8; mv a1,s8; mv a2,s6; li a3,MM
    call vlkf_add

    li t0,WS_L; add s11,s7,t0
    mv a0,s11; mv a1,s8; li a2,69
    call cholesky_decomposition


    li t0,WS_PHtT; add s8,s7,t0
    mv a0,s8; mv a1,s10; li a2,276; li a3,69
    call vlkf_transpose

    li t0,WS_Kt; add s9,s7,t0
    mv a0,s9; mv a1,s11; mv a2,s8; li a3,69; li a4,276
    call chol_solve_matrix

    # K = Kt^T  
    li t0,WS_K; add s8,s7,t0
    mv a0,s8; mv a1,s9; li a2,69; li a3,276
    call vlkf_transpose

    # Ky = K*y  
    li t0,WS_Ky; add s9,s7,t0
    li t1,WS_y;  add t1,s7,t1
    mv a0,s9; mv a1,s8; mv a2,t1; li a3,276; li a4,69; li a5,1
    call vlkf_matmul

    # x += Ky
    mv a0,s0; mv a1,s9; li a2,276
    call vlkf_add_inplace

    # KH = K*H  (276×69 * 69×276 = 276×276)
    li t0,WS_KH; add s9,s7,t0
    mv a0,s9; mv a1,s8; mv a2,s5; li a3,276; li a4,69; li a5,276
    call vlkf_matmul

    # IKH = I - KH
    li t0,WS_IKH; add s10,s7,t0
    mv a0,s10; mv a1,s9; li a2,276
    call vlkf_identity_minus

    # IKH_P = (I-KH)*P
    li t0,WS_IKH_P; add s9,s7,t0
    mv a0,s9; mv a1,s10; mv a2,s1; li a3,276; li a4,276; li a5,276
    call vlkf_matmul

    # IKHT = (I-KH)^T
    li t0,WS_IKHT; add s11,s7,t0
    mv a0,s11; mv a1,s10; li a2,276; li a3,276
    call vlkf_transpose

    
    li t0,WS_FPFt; add s10,s7,t0
    mv a0,s10; mv a1,s9; mv a2,s11; li a3,276; li a4,276; li a5,276
    call vlkf_matmul

    # KR = K*R  
    li t0,WS_KR;  add s9,s7,t0
    li t1,WS_K;   add t1,s7,t1
    mv a0,s9; mv a1,t1; mv a2,s6; li a3,276; li a4,69; li a5,69
    call vlkf_matmul

    # KTJ = K^T  
    li t0,WS_KTJ; add s11,s7,t0
    li t1,WS_K;   add t1,s7,t1
    mv a0,s11; mv a1,t1; li a2,276; li a3,69
    call vlkf_transpose


    li t0,WS_KRKt; add s8,s7,t0
    mv a0,s8; mv a1,s9; mv a2,s11; li a3,276; li a4,69; li a5,276
    call vlkf_matmul

    mv a0,s1; mv a1,s10; mv a2,s8; li a3,NN
    call vlkf_add

    ld ra,104(sp); ld s0,96(sp);  ld s1,88(sp);  ld s2,80(sp)
    ld s3,72(sp);  ld s4,64(sp);  ld s5,56(sp);  ld s6,48(sp)
    ld s7,40(sp);  ld s8,32(sp);  ld s9,24(sp);  ld s10,16(sp); ld s11,8(sp)
    addi sp,sp,112; ret


# vlkf_zero(ptr, n)

.globl vlkf_zero
.align 2
vlkf_zero:
    mv t0,a0; mv t1,a1
.Lvz_loop:
    beqz t1,.Lvz_done
    vsetvli t2,t1,e64,m8,ta,ma   # grant up to VLEN*8 / 64 lanes
    vmv.v.i v0,0                  # all lanes = integer 0 = +0.0 in IEEE 754
    vse64.v v0,(t0)
    slli t3,t2,3; add t0,t0,t3   # ptr  += vl * 8
    sub  t1,t1,t2                 # rem  -= vl
    j .Lvz_loop
.Lvz_done:
    ret

# vlkf_copy(dst, src, n)

.globl vlkf_copy
.align 2
vlkf_copy:
    mv t0,a0; mv t1,a1; mv t2,a2
.Lvc_loop:
    beqz t2,.Lvc_done
    vsetvli t3,t2,e64,m8,ta,ma
    vle64.v v0,(t1)
    vse64.v v0,(t0)
    slli t4,t3,3
    add t0,t0,t4; add t1,t1,t4
    sub t2,t2,t3
    j .Lvc_loop
.Lvc_done:
    ret


# vlkf_add(C, A, B, n)

.globl vlkf_add
.align 2
vlkf_add:
    mv t0,a0; mv t1,a1; mv t2,a2; mv t3,a3
.Lvadd_loop:
    beqz t3,.Lvadd_done
    vsetvli t4,t3,e64,m4,ta,ma
    vle64.v v0,(t1)
    vle64.v v4,(t2)
    vfadd.vv v8,v0,v4
    vse64.v v8,(t0)
    slli t5,t4,3
    add t0,t0,t5; add t1,t1,t5; add t2,t2,t5
    sub t3,t3,t4
    j .Lvadd_loop
.Lvadd_done:
    ret

# vlkf_sub(C, A, B, n)

.globl vlkf_sub
.align 2
vlkf_sub:
    mv t0,a0; mv t1,a1; mv t2,a2; mv t3,a3
.Lvsub_loop:
    beqz t3,.Lvsub_done
    vsetvli t4,t3,e64,m4,ta,ma
    vle64.v v0,(t1)
    vle64.v v4,(t2)
    vfsub.vv v8,v0,v4
    vse64.v v8,(t0)
    slli t5,t4,3
    add t0,t0,t5; add t1,t1,t5; add t2,t2,t5
    sub t3,t3,t4
    j .Lvsub_loop
.Lvsub_done:
    ret

# vlkf_add_inplace(a, b, n)

.globl vlkf_add_inplace
.align 2
vlkf_add_inplace:
    mv t0,a0; mv t1,a1; mv t2,a2
.Lvai_loop:
    beqz t2,.Lvai_done
    vsetvli t3,t2,e64,m4,ta,ma
    vle64.v v0,(t0)
    vle64.v v4,(t1)
    vfadd.vv v8,v0,v4
    vse64.v v8,(t0)
    slli t4,t3,3
    add t0,t0,t4; add t1,t1,t4
    sub t2,t2,t3
    j .Lvai_loop
.Lvai_done:
    ret

# vlkf_matmul(C, A, B, M, K, N)
# C = A * B
# a0=C a1=A a2=B a3=M a4=K a5=N
.globl vlkf_matmul
.align 2
vlkf_matmul:
    addi sp,sp,-96
    sd ra,88(sp); sd s0,80(sp); sd s1,72(sp); sd s2,64(sp)
    sd s3,56(sp); sd s4,48(sp); sd s5,40(sp); sd s6,32(sp)
    sd s7,24(sp); sd s8,16(sp); sd s9,8(sp)

    mv s0,a0; mv s1,a1; mv s2,a2
    mv s3,a3; mv s4,a4; mv s5,a5

    slli s6,s4,3          # s6 = K*8  (A row stride)
    slli s7,s5,3          # s7 = N*8  (B/C row stride)

    li   s8,0             # i = 0
.Lvmm_i:
    bge  s8,s3,.Lvmm_done

    # s9 = &A[i,0]
    mul  t0,s8,s6; add s9,s1,t0

    li   a6,0             # j = 0  (a6 is caller-saved, free to use)
.Lvmm_j:
    bge  a6,s5,.Lvmm_next_i

    # Set VL for this horizontal chunk (columns j..j+vl-1)
    sub  t0,s5,a6
    vsetvli t1,t0,e64,m8,ta,ma   # t1 = vl

    # Accumulator v0 = 0  (uses v0..v7 under m8)
    vmv.v.i v0,0

    mv   t2,s9            # &A[i,0]
    slli t0,a6,3
    add  t3,s2,t0         # &B[0,j]
    mv   t4,s4            # k = K

.Lvmm_k:
    beqz t4,.Lvmm_kend

    fld  ft0,0(t2)        # A[i,k] scalar
    vle64.v v8,(t3)       # B[k, j..j+vl-1]  (CONTIGUOUS)
    vfmacc.vf v0,ft0,v8   # v0 += A[i,k] * B[k, j..]

    addi t2,t2,8          # next A
    add  t3,t3,s7         # next B row
    addi t4,t4,-1
    j    .Lvmm_k

.Lvmm_kend:
    # Store C[i, j..j+vl-1]
    mul  t0,s8,s7; add t0,s0,t0
    slli t5,a6,3; add t0,t0,t5
    vse64.v v0,(t0)

    add  a6,a6,t1         # j += vl
    j    .Lvmm_j

.Lvmm_next_i:
    addi s8,s8,1
    j    .Lvmm_i

.Lvmm_done:
    ld ra,88(sp); ld s0,80(sp); ld s1,72(sp); ld s2,64(sp)
    ld s3,56(sp); ld s4,48(sp); ld s5,40(sp); ld s6,32(sp)
    ld s7,24(sp); ld s8,16(sp); ld s9,8(sp)
    addi sp,sp,96; ret

# vlkf_transpose(At, A, rows, cols)

.globl vlkf_transpose
.align 2
vlkf_transpose:
    addi sp,sp,-48
    sd ra,40(sp); sd s0,32(sp); sd s1,24(sp)
    sd s2,16(sp); sd s3,8(sp);  sd s4,0(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3

    slli s4,s2,3    # stride for scatter = rows * 8 bytes

    li t5,0         # i = row of A
.Lvtr_i:
    bge t5,s2,.Lvtr_done

    # src = &A[i, 0]  = A + i*cols*8
    mul t0,t5,s3; slli t0,t0,3; add t3,s1,t0
    # dst = &At[0, i] = At + i*8  (start of column i in At)
    slli t0,t5,3;               add t4,s0,t0

    mv t0,t3        # running src ptr
    mv t1,t4        # running dst ptr
    mv t2,s3        # remaining cols

.Lvtr_j:
    beqz t2,.Lvtr_ni
    vsetvli a0,t2,e64,m4,ta,ma  # a0 = vl
    vle64.v v0,(t0)              # load A[i, j..j+vl-1]  (contiguous)
    vsse64.v v0,(t1),s4          # scatter At[j,i], At[j+1,i], ...  stride=rows*8
    slli a1,a0,3; add t0,t0,a1  # src advance = vl * 8  (contiguous)
    mul  a2,a0,s4; add t1,t1,a2 # dst advance = vl * rows * 8  (strided)
    sub  t2,t2,a0
    j .Lvtr_j

.Lvtr_ni:
    addi t5,t5,1; j .Lvtr_i

.Lvtr_done:
    ld ra,40(sp); ld s0,32(sp); ld s1,24(sp)
    ld s2,16(sp); ld s3,8(sp);  ld s4,0(sp)
    addi sp,sp,48; ret

# vlkf_identity_minus(C, A, n)

.globl vlkf_identity_minus
.align 2
vlkf_identity_minus:
    addi sp,sp,-32
    sd ra,24(sp); sd s0,16(sp); sd s1,8(sp); sd s2,0(sp)
    mv s0,a0; mv s1,a1; mv s2,a2

    # Step 1: flat negate over n*n elements
    mul t0,s2,s2
    mv t1,s0; mv t2,s1; mv t3,t0
.Lvim_neg:
    beqz t3,.Lvim_diag
    vsetvli t4,t3,e64,m8,ta,ma
    vle64.v v0,(t2)
    vfneg.v v8,v0
    vse64.v v8,(t1)
    slli t5,t4,3
    add t1,t1,t5; add t2,t2,t5
    sub t3,t3,t4
    j .Lvim_neg

.Lvim_diag:
    # Step 2: diagonal += 1.0
    li t0,1; fcvt.d.w ft0,t0    # ft0 = 1.0
    li t1,0
.Lvim_d:
    bge t1,s2,.Lvim_done
    mul t2,t1,s2; add t2,t2,t1  # offset = i*n + i
    slli t2,t2,3; add t2,s0,t2
    fld ft1,0(t2)
    fadd.d ft1,ft1,ft0
    fsd ft1,0(t2)
    addi t1,t1,1; j .Lvim_d

.Lvim_done:
    ld ra,24(sp); ld s0,16(sp); ld s1,8(sp); ld s2,0(sp)
    addi sp,sp,32; ret

.globl init_identity
.align 2
init_identity:
    addi sp,sp,-32; sd ra,24(sp); sd s0,16(sp); sd s1,8(sp)
    mv s0,a0; mv s1,a1
    mul a1,s1,s1; mv a0,s0; call vlkf_zero   # zero n*n via RVV
    li t2,1; fcvt.d.w ft0,t2; li t0,0
.Lii_l:
    bge t0,s1,.Lii_d
    mul t1,t0,s1; add t1,t1,t0; slli t1,t1,3; add t1,s0,t1
    fsd ft0,0(t1)
    addi t0,t0,1; j .Lii_l
.Lii_d:
    ld ra,24(sp); ld s0,16(sp); ld s1,8(sp); addi sp,sp,32; ret

# chol_solve_matrix — solve LL^T * X = B  column by column

.globl chol_solve_matrix
.align 2
chol_solve_matrix:
    addi sp,sp,-80
    sd ra,72(sp); sd s0,64(sp); sd s1,56(sp); sd s2,48(sp)
    sd s3,40(sp); sd s4,32(sp); sd s5,24(sp); sd s6,16(sp); sd s7,8(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3; mv s4,a4

    # Allocate two n-element column buffers
    slli a0,s3,3; call malloc; mv s5,a0   # s5 = y_buf  (fwd solve output / col temp)
    slli a0,s3,3; call malloc; mv s6,a0   # s6 = x_buf  (bwd solve output)

    # stride for column access in B and X = m * 8 bytes
    slli s7,s4,3                           # s7 = stride_B_X = m*8

    li t6,0                                # col = 0
.Lcsm_col:
    bge t6,s4,.Lcsm_done


    # &B[0,col] = B + col*8
    slli t0,t6,3; add t0,s2,t0            # t0 = &B[0,col]
    mv a0,s5                              # dst = y_buf
    mv a1,t0                              # src = &B[0,col]
    mv a2,s7                              # stride = m*8
    mv a3,s3                              # element count = n
    call vlkf_strided_load                # y_buf[i] = B[i,col]

    mv a0,s6; mv a1,s1; mv a2,s5; mv a3,s3
    call forward_substitution


    mv a0,s5; mv a1,s1; mv a2,s6; mv a3,s3
    call backward_substitution


    slli t0,t6,3; add t0,s0,t0            # t0 = &X[0,col]
    mv a0,t0                              # dst = &X[0,col]
    mv a1,s5                              # src = y_buf
    mv a2,s7                              # stride = m*8
    mv a3,s3                              # element count = n
    call vlkf_strided_store               # X[i,col] = y_buf[i]

    addi t6,t6,1; j .Lcsm_col

.Lcsm_done:
    mv a0,s5; call free
    mv a0,s6; call free
    ld ra,72(sp); ld s0,64(sp); ld s1,56(sp); ld s2,48(sp)
    ld s3,40(sp); ld s4,32(sp); ld s5,24(sp); ld s6,16(sp); ld s7,8(sp)
    addi sp,sp,80; ret


# vlkf_strided_load(dst, src, stride, n)

.globl vlkf_strided_load
.align 2
vlkf_strided_load:
    mv t0,a0; mv t1,a1; mv t2,a3   # t2 = remaining elements
.Lvsl_loop:
    beqz t2,.Lvsl_done
    vsetvli t3,t2,e64,m4,ta,ma
    vlse64.v v0,(t1),a2            # load vl elements from src, stride a2
    vse64.v  v0,(t0)               # store contiguous to dst
    slli t4,t3,3; add t0,t0,t4    # dst  += vl * 8
    mul  t4,t3,a2; add t1,t1,t4   # src  += vl * stride
    sub  t2,t2,t3
    j .Lvsl_loop
.Lvsl_done:
    ret

# vlkf_strided_store(dst, src, stride, n)

.globl vlkf_strided_store
.align 2
vlkf_strided_store:
    mv t0,a0; mv t1,a1; mv t2,a3
.Lvss_loop:
    beqz t2,.Lvss_done
    vsetvli t3,t2,e64,m4,ta,ma
    vle64.v  v0,(t1)               # load vl elements contiguous from src
    vsse64.v v0,(t0),a2            # scatter to dst with stride a2
    slli t4,t3,3; add t1,t1,t4    # src  += vl * 8
    mul  t4,t3,a2; add t0,t0,t4   # dst  += vl * stride
    sub  t2,t2,t3
    j .Lvss_loop
.Lvss_done:
    ret

# SCALAR: cholesky_decomposition(L, A, n)

.globl cholesky_decomposition
.align 2
cholesky_decomposition:
    addi sp,sp,-96
    sd   ra,88(sp); sd s0,80(sp); sd s1,72(sp); sd s2,64(sp)
    sd   s3,56(sp); sd s4,48(sp); sd s5,40(sp)
    fsd  fs0,32(sp); fsd fs1,24(sp)

    mv s0,a0; mv s1,a1; mv s2,a2   # L, A, n
    li s3,0                         # i = 0

.Lchol_i:
    bge s3,s2,.Lchol_success

    li s4,0                         # j = 0
.Lchol_j:
    bgt s4,s3,.Lchol_next_i

    # fs0 = A[i,j]
    mul t0,s3,s2; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0
    fld fs0,0(t0)

    # fs0 -= sum_{k=0}^{j-1} L[i,k]*L[j,k]
    li s5,0
.Lchol_k:
    bge s5,s4,.Lchol_offdiag
    mul t1,s3,s2; add t1,t1,s5; slli t1,t1,3; add t1,s0,t1; fld ft0,0(t1)  # L[i,k]
    mul t2,s4,s2; add t2,t2,s5; slli t2,t2,3; add t2,s0,t2; fld ft1,0(t2)  # L[j,k]
    fnmsub.d fs0,ft0,ft1,fs0    # fs0 -= L[i,k]*L[j,k]
    addi s5,s5,1; j .Lchol_k

.Lchol_offdiag:
    beq s3,s4,.Lchol_diag

    # L[i,j] = fs0 / L[j,j]
    mul t1,s4,s2; add t1,t1,s4; slli t1,t1,3; add t1,s0,t1; fld fs1,0(t1)
    fdiv.d fs0,fs0,fs1
    mul t1,s3,s2; add t1,t1,s4; slli t1,t1,3; add t1,s0,t1; fsd fs0,0(t1)
    addi s4,s4,1; j .Lchol_j

.Lchol_diag:
    # L[i,i] = sqrt(fs0); check fs0 >= 0
    la t0,ZERO_D; fld fs1,0(t0)
    flt.d t1,fs0,fs1; bnez t1,.Lchol_fail
    fsqrt.d fs0,fs0
    mul t1,s3,s2; add t1,t1,s4; slli t1,t1,3; add t1,s0,t1; fsd fs0,0(t1)
    addi s4,s4,1; j .Lchol_j

.Lchol_next_i:
    # Zero upper triangle: L[i,j] = 0 for j > i
    addi s4,s3,1
.Lchol_zero_upper:
    bge s4,s2,.Lchol_inc_i
    mul t1,s3,s2; add t1,t1,s4; slli t1,t1,3; add t1,s0,t1
    la t0,ZERO_D; fld fs0,0(t0); fsd fs0,0(t1)
    addi s4,s4,1; j .Lchol_zero_upper
.Lchol_inc_i:
    addi s3,s3,1; j .Lchol_i

.Lchol_success:
    li a0,0; j .Lchol_epilog
.Lchol_fail:
    li a0,-1
.Lchol_epilog:
    ld ra,88(sp); ld s0,80(sp); ld s1,72(sp); ld s2,64(sp)
    ld s3,56(sp); ld s4,48(sp); ld s5,40(sp)
    fld fs0,32(sp); fld fs1,24(sp)
    addi sp,sp,96; ret

# SCALAR: forward_substitution(y, L, b, n)

.globl forward_substitution
.align 2
forward_substitution:
    addi sp,sp,-64
    sd ra,56(sp); sd s0,48(sp); sd s1,40(sp); sd s2,32(sp)
    sd s3,24(sp); sd s4,16(sp); sd s5,8(sp)
    fsd fs0,0(sp)

    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    li s4,0

.Lfwd_i:
    bge s4,s3,.Lfwd_done
    fcvt.d.w fs0,zero              # sum = 0
    li s5,0
.Lfwd_k:
    bge s5,s4,.Lfwd_end_k
    mul t0,s4,s3; add t0,t0,s5; slli t0,t0,3; add t0,s1,t0; fld ft0,0(t0)
    slli t1,s5,3; add t1,s0,t1; fld ft1,0(t1)
    fmadd.d fs0,ft0,ft1,fs0
    addi s5,s5,1; j .Lfwd_k
.Lfwd_end_k:
    slli t0,s4,3; add t0,s2,t0; fld ft0,0(t0)   # b[i]
    fsub.d ft0,ft0,fs0
    mul t0,s4,s3; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0; fld ft1,0(t0)
    fdiv.d ft0,ft0,ft1
    slli t0,s4,3; add t0,s0,t0; fsd ft0,0(t0)
    addi s4,s4,1; j .Lfwd_i

.Lfwd_done:
    ld ra,56(sp); ld s0,48(sp); ld s1,40(sp); ld s2,32(sp)
    ld s3,24(sp); ld s4,16(sp); ld s5,8(sp)
    fld fs0,0(sp); addi sp,sp,64; ret

# SCALAR: backward_substitution(x, L, y, n)

.globl backward_substitution
.align 2
backward_substitution:
    addi sp,sp,-64
    sd ra,56(sp); sd s0,48(sp); sd s1,40(sp); sd s2,32(sp)
    sd s3,24(sp); sd s4,16(sp); sd s5,8(sp)
    fsd fs0,0(sp)

    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    addi s4,s3,-1                  # i = n-1

.Lbwd_i:
    bltz s4,.Lbwd_done
    fcvt.d.w fs0,zero              # sum = 0
    addi s5,s4,1                   # k = i+1
.Lbwd_k:
    bge s5,s3,.Lbwd_end_k
    mul t0,s5,s3; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0; fld ft0,0(t0)  # L[k,i]
    slli t1,s5,3; add t1,s0,t1; fld ft1,0(t1)                                # x[k]
    fmadd.d fs0,ft0,ft1,fs0
    addi s5,s5,1; j .Lbwd_k
.Lbwd_end_k:
    slli t0,s4,3; add t0,s2,t0; fld ft0,0(t0)   # y[i]
    fsub.d ft0,ft0,fs0
    mul t0,s4,s3; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0; fld ft1,0(t0)
    fdiv.d ft0,ft0,ft1
    slli t0,s4,3; add t0,s0,t0; fsd ft0,0(t0)
    addi s4,s4,-1; j .Lbwd_i

.Lbwd_done:
    ld ra,56(sp); ld s0,48(sp); ld s1,40(sp); ld s2,32(sp)
    ld s3,24(sp); ld s4,16(sp); ld s5,8(sp)
    fld fs0,0(sp); addi sp,sp,64; ret

# buildF — Construct 276×276 state-transition matrix 

.globl buildF
.align 2
buildF:
    addi sp,sp,-64
    sd  ra,56(sp); sd s0,48(sp); sd s1,40(sp); sd s2,32(sp)
    fsd fs0,24(sp); fsd fs1,16(sp); fsd fs2,8(sp); fsd fs3,0(sp)

    mv   s0,a0; fmv.d fs0,fa0        # s0=F, fs0=dt

    fmul.d fs1,fs0,fs0               # fs1 = dt²
    fmul.d fs2,fs1,fs0               # fs2 = dt³
    li t0,2; fcvt.d.w ft0,t0; fdiv.d fs3,fs1,ft0    # fs3 = dt²/2
    li t0,6; fcvt.d.w ft0,t0; fdiv.d ft1,fs2,ft0
    fmv.d fs4,ft1                                     # fs4 = dt³/6

    # Zero F (276*276 = 76176 doubles) via RVV
    mv a0,s0; li a1,NN; call vlkf_zero

    li s1,0                           # joint = 0
.LbF_joint:
    li t0,23; bge s1,t0,.LbF_done
    li s2,0
.LbF_axis:
    li t0,3; bge s2,t0,.LbF_next_joint

    # base index: joint*12 + axis*4
    li t0,12; mul t1,s1,t0; slli t2,s2,2; add t1,t1,t2

    # Row 0 ptr: F + (base*276 + base)*8
    li t0,276; mul t2,t1,t0; add t2,t2,t1; slli t2,t2,3; add t2,s0,t2

    fcvt.d.w ft0,zero
    li t3,1; fcvt.d.w ft1,t3
    fsd ft1,0(t2); fsd fs0,8(t2); fsd fs3,16(t2); fsd fs4,24(t2)  # row 0
    li t3,2208; add t2,t2,t3
    fsd ft0,0(t2); fsd ft1,8(t2); fsd fs0,16(t2); fsd fs3,24(t2)  # row 1
    li t3,2208; add t2,t2,t3
    fsd ft0,0(t2); fsd ft0,8(t2); fsd ft1,16(t2); fsd fs0,24(t2)  # row 2
    li t3,2208; add t2,t2,t3
    fsd ft0,0(t2); fsd ft0,8(t2); fsd ft0,16(t2); fsd ft1,24(t2)  # row 3

    addi s2,s2,1; j .LbF_axis
.LbF_next_joint:
    addi s1,s1,1; j .LbF_joint
.LbF_done:
    ld  ra,56(sp); ld s0,48(sp); ld s1,40(sp); ld s2,32(sp)
    fld fs0,24(sp); fld fs1,16(sp); fld fs2,8(sp); fld fs3,0(sp)
    addi sp,sp,64; ret

# buildQ — Construct 276×276 process-noise covariance 

.globl buildQ
.align 2
buildQ:
    addi sp,sp,-96
    sd  ra,88(sp); sd s0,80(sp); sd s1,72(sp); sd s2,64(sp)
    fsd fs0,56(sp); fsd fs1,48(sp); fsd fs2,40(sp); fsd fs3,32(sp)
    fsd fs4,24(sp); fsd fs5,16(sp); fsd fs6,8(sp);  fsd fs7,0(sp)

    mv   s0,a0; fmv.d fs0,fa0; fmv.d fs1,fa1   # s0=Q, fs0=dt, fs1=sigma

    fmul.d fs2,fs0,fs0    # dt²
    fmul.d ft0,fs2,fs0    # dt³
    fmul.d ft1,ft0,fs0    # dt⁴
    fmul.d ft2,ft1,fs0    # dt⁵
    fmul.d ft3,ft2,fs0    # dt⁶
    fmul.d fs3,fs1,fs1    # sigma²

    # Zero Q (76176 doubles) via RVV
    mv a0,s0; li a1,NN; call vlkf_zero

    li s1,0
.LbQ_joint:
    li t0,23; bge s1,t0,.LbQ_done
    li s2,0
.LbQ_axis:
    li t0,3; bge s2,t0,.LbQ_next_joint

    li t0,12; mul t1,s1,t0; slli t2,s2,2; add t1,t1,t2

    # Compute Q sub-block entries
    li t0,36; fcvt.d.w ft4,t0; fdiv.d ft5,ft3,ft4; fmul.d fs4,ft5,fs3  # dt⁶/36*σ²
    li t0,12; fcvt.d.w ft4,t0; fdiv.d ft6,ft2,ft4; fmul.d fs5,ft6,fs3  # dt⁵/12*σ²
    li t0,6;  fcvt.d.w ft4,t0; fdiv.d ft7,ft1,ft4; fmul.d fs6,ft7,fs3  # dt⁴/6*σ²
              fdiv.d ft8,ft0,ft4; fmul.d fs7,ft8,fs3                     # dt³/6*σ²
    li t0,4;  fcvt.d.w ft4,t0; fdiv.d ft9,ft1,ft4; fmul.d fa2,ft9,fs3  # dt⁴/4*σ²
    li t0,2;  fcvt.d.w ft4,t0; fdiv.d ft10,ft0,ft4; fmul.d fa3,ft10,fs3 # dt³/2*σ²
              fdiv.d ft11,fs2,ft4; fmul.d fa4,ft11,fs3                   # dt²/2*σ²
    fmul.d fa5,fs2,fs3    # dt²*σ²
    fmul.d fa6,fs0,fs3    # dt*σ²
    fmv.d  fa7,fs3        # σ²

    li t0,276; mul t2,t1,t0; add t2,t2,t1; slli t2,t2,3; add t2,s0,t2

    fsd fs4,0(t2); fsd fs5,8(t2); fsd fs6,16(t2); fsd fs7,24(t2)  # row 0
    li t3,2208; add t2,t2,t3
    fsd fs5,0(t2); fsd fa2,8(t2); fsd fa3,16(t2); fsd fa4,24(t2)  # row 1
    li t3,2208; add t2,t2,t3
    fsd fs6,0(t2); fsd fa3,8(t2); fsd fa5,16(t2); fsd fa6,24(t2)  # row 2
    li t3,2208; add t2,t2,t3
    fsd fs7,0(t2); fsd fa4,8(t2); fsd fa6,16(t2); fsd fa7,24(t2)  # row 3

    addi s2,s2,1; j .LbQ_axis
.LbQ_next_joint:
    addi s1,s1,1; j .LbQ_joint
.LbQ_done:
    ld  ra,88(sp); ld s0,80(sp); ld s1,72(sp); ld s2,64(sp)
    fld fs0,56(sp); fld fs1,48(sp); fld fs2,40(sp); fld fs3,32(sp)
    fld fs4,24(sp); fld fs5,16(sp); fld fs6,8(sp);  fld fs7,0(sp)
    addi sp,sp,96; ret

# buildH — Construct 69×276 measurement matrix 
.globl buildH
.align 2
buildH:
    addi sp,sp,-32
    sd ra,24(sp); sd s0,16(sp); sd s1,8(sp); sd s2,0(sp)

    mv s0,a0

    # Zero H (69*276 = 19044 doubles) via RVV
    mv a0,s0; li a1,NM; call vlkf_zero

    li t0,1; fcvt.d.w ft1,t0    # ft1 = 1.0
    li s1,0
.LbH_joint:
    li t0,23; bge s1,t0,.LbH_done
    li s2,0
.LbH_axis:
    li t0,3; bge s2,t0,.LbH_next_joint

    # row = joint*3 + axis,  col = joint*12 + axis*4
    li t0,3;  mul t1,s1,t0; add t1,t1,s2
    li t0,12; mul t2,s1,t0; slli t3,s2,2; add t2,t2,t3

    li t0,276; mul t3,t1,t0; add t3,t3,t2; slli t3,t3,3; add t3,s0,t3
    fsd ft1,0(t3)   # H[row,col] = 1.0

    addi s2,s2,1; j .LbH_axis
.LbH_next_joint:
    addi s1,s1,1; j .LbH_joint
.LbH_done:
    ld ra,24(sp); ld s0,16(sp); ld s1,8(sp); ld s2,0(sp)
    addi sp,sp,32; ret

.globl matrix_zero
matrix_zero:         j vlkf_zero

.globl matrix_copy
matrix_copy:         j vlkf_copy

.globl matrix_add
matrix_add:          j vlkf_add

.globl matrix_subtract
matrix_subtract:     j vlkf_sub

.globl vector_add_inplace
vector_add_inplace:  j vlkf_add_inplace

.globl matrix_transpose
matrix_transpose:    j vlkf_transpose

.globl identity_minus
identity_minus:      j vlkf_identity_minus

.globl matrix_multiply
matrix_multiply:     j vlkf_matmul