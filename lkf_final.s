# lkf_final.s - Complete Standalone Linear Kalman Filter
# Team ChaiGPT - Milestone 3

# WORKSPACE LAYOUT AND CONSTANTS

.equ WS_xtmp,   0
.equ WS_FP,     2208
.equ WS_FPFt,   611616
.equ WS_Ft,     1221024
.equ WS_Ht,     1830432
.equ WS_PHt,    1982784
.equ WS_PHtT,   2135136
.equ WS_Kt,     2287488
.equ WS_K,      2439840
.equ WS_S,      2592192
.equ WS_L,      2630280
.equ WS_y,      2668368
.equ WS_Hx,     2668920
.equ WS_Ky,     2669472
.equ WS_KH,     2671680
.equ WS_IKH,    3281088
.equ WS_IKHT,   3890496
.equ WS_IKH_P,  4499904
.equ WS_KR,     5109312
.equ WS_KTJ,    5261664
.equ WS_KRKt,   5414016
.equ LKF_WS_TOTAL, 6023424

.equ NN, 76176
.equ NM, 19044
.equ MM, 4761

.section .text

# lkf_asm.s  –  Linear Kalman Filter in RISC-V scalar assembly  –  Milestone 3
# Team ChaiGPT: Laiba(31736) Ameer(30535) Usman(29126) Arrayyan(30557)
# Workspace layout (byte offsets inside the pre-malloc'd block):
#   WS_xtmp   =0        276 doubles   (x predict temp)
#   WS_FP     =2208     276×276       (F*P)
#   WS_FPFt   =611616   276×276       (F*P*F^T)
#   WS_Ft     =1221024  276×276       (F transposed)
#   WS_Ht     =1830432  276×69        (H transposed)
#   WS_PHt    =1982784  276×69
#   WS_PHtT   =2135136  69×276        (rhs for Cholesky solve)
#   WS_Kt     =2287488  69×276        (Cholesky solve result K^T)
#   WS_K      =2439840  276×69        (Kalman gain)
#   WS_S      =2592192  69×69         (innovation cov)
#   WS_L      =2630280  69×69         (Cholesky L)
#   WS_y      =2668368  69 doubles    (innovation z-Hx)
#   WS_Hx     =2668920  69 doubles
#   WS_Ky     =2669472  276 doubles   (K*y)
#   WS_KH     =2671680  276×276
#   WS_IKH    =3281088  276×276       (I-KH)
#   WS_IKHT   =3890496  276×276       ((I-KH)^T)
#   WS_IKH_P  =4499904  276×276       ((I-KH)*P)
#   WS_KR     =5109312  276×69        (K*R)
#   WS_KTJ    =5261664  69×276        (K^T for Joseph)
#   WS_KRKt   =5414016  276×276       (K*R*K^T)
#   TOTAL     =6023424  (~5.7 MB)

.equ WS_xtmp,   0
.equ WS_FP,     2208
.equ WS_FPFt,   611616
.equ WS_Ft,     1221024
.equ WS_Ht,     1830432
.equ WS_PHt,    1982784
.equ WS_PHtT,   2135136
.equ WS_Kt,     2287488
.equ WS_K,      2439840
.equ WS_S,      2592192
.equ WS_L,      2630280
.equ WS_y,      2668368
.equ WS_Hx,     2668920
.equ WS_Ky,     2669472
.equ WS_KH,     2671680
.equ WS_IKH,    3281088
.equ WS_IKHT,   3890496
.equ WS_IKH_P,  4499904
.equ WS_KR,     5109312
.equ WS_KTJ,    5261664
.equ WS_KRKt,   5414016
.equ LKF_WS_TOTAL, 6023424

.equ NN, 76176      # 276*276
.equ NM, 19044      # 276*69
.equ MM, 4761       # 69*69

.section .text

# Usage pattern: li t0,OFFSET; add REG,s_ws,t0
# run_lkf_asm(meas,num_frames,F,Q,H,R,all_states)
# a0=meas  a1=num_frames  a2=F  a3=Q  a4=H  a5=R  a6=all_states

.globl run_lkf_asm
.align 2
run_lkf_asm:
    addi sp,sp,-128
    sd ra,120(sp); sd s0,112(sp); sd s1,104(sp); sd s2,96(sp)
    sd s3,88(sp);  sd s4,80(sp);  sd s5,72(sp);  sd s6,64(sp)
    sd s7,56(sp);  sd s8,48(sp);  sd s9,40(sp);  sd s10,32(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    mv s4,a4; mv s5,a5; mv s6,a6

    # Workspace (reused every frame)
    li a0,LKF_WS_TOTAL; call malloc; mv s7,a0
    # State x (276 doubles)
    li a0,2208; call malloc; mv s8,a0
    # Covariance P (276×276 doubles)
    li a0,609408; call malloc; mv s9,a0

    # memset(x, 0, 2208)
    mv a0,s8; li a1,0; li a2,2208; call memset

    # Init x positions from first measurement
    li t0,0
.Lrlkf_ix:
    li t1,23; bge t0,t1,.Lrlkf_ix_d
    li t2,24; mul t3,t0,t2; add t3,s0,t3   # &meas[j*3]
    li t2,96; mul t4,t0,t2; add t4,s8,t4   # &x[j*12]
    fld ft0,0(t3);  fsd ft0,0(t4)    # px -> x[j*12+0]
    fld ft1,8(t3);  fsd ft1,32(t4)   # py -> x[j*12+4]
    fld ft2,16(t3); fsd ft2,64(t4)   # pz -> x[j*12+8]
    addi t0,t0,1; j .Lrlkf_ix
.Lrlkf_ix_d:

    # P = I_276
    mv a0,s9; li a1,276; call init_identity

    # Frame loop
    li s10,0
.Lrlkf_loop:
    bge s10,s1,.Lrlkf_done
    # z pointer: meas + k*69*8
    li t0,552; mul t0,s10,t0; add t0,s0,t0  # &meas[k*69]
    # lkf_step(x,P,z,F,Q,H,R,ws)
    mv a0,s8; mv a1,s9; mv a2,t0
    mv a3,s2; mv a4,s3; mv a5,s4; mv a6,s5; mv a7,s7
    call lkf_step
    # Store x into all_states[k*276..]
    li t0,2208; mul t0,s10,t0; add t0,s6,t0
    mv a0,t0; mv a1,s8; li a2,276; call matrix_copy
    addi s10,s10,1; j .Lrlkf_loop
.Lrlkf_done:
    mv a0,s7; call free
    mv a0,s8; call free
    mv a0,s9; call free
    ld ra,120(sp); ld s0,112(sp); ld s1,104(sp); ld s2,96(sp)
    ld s3,88(sp);  ld s4,80(sp);  ld s5,72(sp);  ld s6,64(sp)
    ld s7,56(sp);  ld s8,48(sp);  ld s9,40(sp);  ld s10,32(sp)
    addi sp,sp,128; ret

# lkf_step(x,P,z,F,Q,H,R,ws) — one predict+update step
# a0=x(276) a1=P(276×276) a2=z(69) a3=F(276×276) a4=Q a5=H(69×276) a6=R(69×69) a7=ws
#
# Register allocation (callee-saved):
#   s0=x  s1=P  s2=z  s3=F  s4=Q  s5=H  s6=R  s7=ws
#   s8..s11 = computed workspace sub-pointers

.globl lkf_step
.align 2
lkf_step:
    addi sp,sp,-112
    sd ra,104(sp); sd s0,96(sp);  sd s1,88(sp);  sd s2,80(sp)
    sd s3,72(sp);  sd s4,64(sp);  sd s5,56(sp);  sd s6,48(sp)
    sd s7,40(sp);  sd s8,32(sp);  sd s9,24(sp);  sd s10,16(sp); sd s11,8(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3; mv s4,a4; mv s5,a5; mv s6,a6; mv s7,a7

    # x = F*x  (276×276 * 276×1)
    li t0,WS_xtmp; add s8,s7,t0
    mv a0,s8; mv a1,s3; mv a2,s0; li a3,276; li a4,276; li a5,1
    call matrix_multiply
    mv a0,s0; mv a1,s8; li a2,276; call matrix_copy

    # FP = F*P  (276×276 * 276×276)
    li t0,WS_FP; add s8,s7,t0
    mv a0,s8; mv a1,s3; mv a2,s1; li a3,276; li a4,276; li a5,276
    call matrix_multiply

    # Ft = F^T
    li t0,WS_Ft; add s9,s7,t0
    mv a0,s9; mv a1,s3; li a2,276; li a3,276; call matrix_transpose

    # FPFt = FP * Ft
    li t0,WS_FPFt; add s10,s7,t0
    mv a0,s10; mv a1,s8; mv a2,s9; li a3,276; li a4,276; li a5,276
    call matrix_multiply

    # P = FPFt + Q
    mv a0,s1; mv a1,s10; mv a2,s4; li a3,NN; call matrix_add

    # Hx = H*x  (69×276 * 276×1)
    li t0,WS_Hx; add s8,s7,t0
    mv a0,s8; mv a1,s5; mv a2,s0; li a3,69; li a4,276; li a5,1
    call matrix_multiply

    # y = z - Hx
    li t0,WS_y; add s9,s7,t0
    mv a0,s9; mv a1,s2; mv a2,s8; li a3,69; call matrix_subtract

    # Ht = H^T  (276×69)
    li t0,WS_Ht; add s8,s7,t0
    mv a0,s8; mv a1,s5; li a2,69; li a3,276; call matrix_transpose

    # PHt = P * Ht  (276×276 * 276×69)
    li t0,WS_PHt; add s10,s7,t0
    mv a0,s10; mv a1,s1; mv a2,s8; li a3,276; li a4,276; li a5,69
    call matrix_multiply

    # S = H*PHt + R  (69×276 * 276×69 = 69×69)
    li t0,WS_S; add s8,s7,t0
    mv a0,s8; mv a1,s5; mv a2,s10; li a3,69; li a4,276; li a5,69
    call matrix_multiply
    mv a0,s8; mv a1,s8; mv a2,s6; li a3,MM; call matrix_add

    # L = cholesky(S)
    li t0,WS_L; add s11,s7,t0
    mv a0,s11; mv a1,s8; li a2,69; call cholesky_decomposition

    # PHtT = PHt^T  (69×276)
    li t0,WS_PHtT; add s8,s7,t0
    mv a0,s8; mv a1,s10; li a2,276; li a3,69; call matrix_transpose

    # Kt = chol_solve(L, PHtT)  solve LL^T * Kt = PHtT  (69×276)
    li t0,WS_Kt; add s9,s7,t0
    mv a0,s9; mv a1,s11; mv a2,s8; li a3,69; li a4,276
    call chol_solve_matrix

    # K = Kt^T  (276×69)
    li t0,WS_K; add s8,s7,t0
    mv a0,s8; mv a1,s9; li a2,69; li a3,276; call matrix_transpose

    # Ky = K*y  (276×69 * 69×1)
    li t0,WS_Ky; add s9,s7,t0
    li t1,WS_y; add t1,s7,t1
    mv a0,s9; mv a1,s8; mv a2,t1; li a3,276; li a4,69; li a5,1
    call matrix_multiply
    # x += Ky
    mv a0,s0; mv a1,s9; li a2,276; call vector_add_inplace

    # Joseph form: P = (I-KH)*P*(I-KH)^T + K*R*K^T 
    # KH = K*H  (276×69 * 69×276 = 276×276)
    li t0,WS_KH; add s9,s7,t0
    mv a0,s9; mv a1,s8; mv a2,s5; li a3,276; li a4,69; li a5,276
    call matrix_multiply

    # I_KH = I - KH
    li t0,WS_IKH; add s10,s7,t0
    mv a0,s10; mv a1,s9; li a2,276; call identity_minus

    # IKH_P = (I-KH)*P
    li t0,WS_IKH_P; add s9,s7,t0
    mv a0,s9; mv a1,s10; mv a2,s1; li a3,276; li a4,276; li a5,276
    call matrix_multiply

    # IKHT = (I-KH)^T
    li t0,WS_IKHT; add s11,s7,t0
    mv a0,s11; mv a1,s10; li a2,276; li a3,276; call matrix_transpose

    # P = IKH_P * IKHT  (into FPFt slot — safe, predict done)
    li t0,WS_FPFt; add s10,s7,t0
    mv a0,s10; mv a1,s9; mv a2,s11; li a3,276; li a4,276; li a5,276
    call matrix_multiply

    # KR = K*R  (276×69 * 69×69 = 276×69)
    li t0,WS_KR; add s9,s7,t0
    li t1,WS_K; add t1,s7,t1
    mv a0,s9; mv a1,t1; mv a2,s6; li a3,276; li a4,69; li a5,69
    call matrix_multiply

    # KTJ = K^T  (69×276)
    li t0,WS_KTJ; add s11,s7,t0
    li t1,WS_K; add t1,s7,t1
    mv a0,s11; mv a1,t1; li a2,276; li a3,69; call matrix_transpose

    # KRKt = KR * KTJ  (276×69 * 69×276 = 276×276)
    li t0,WS_KRKt; add s8,s7,t0
    mv a0,s8; mv a1,s9; mv a2,s11; li a3,276; li a4,69; li a5,276
    call matrix_multiply

    # P = P_part1 + KRKt
    mv a0,s1; mv a1,s10; mv a2,s8; li a3,NN; call matrix_add

    ld ra,104(sp); ld s0,96(sp);  ld s1,88(sp);  ld s2,80(sp)
    ld s3,72(sp);  ld s4,64(sp);  ld s5,56(sp);  ld s6,48(sp)
    ld s7,40(sp);  ld s8,32(sp);  ld s9,24(sp);  ld s10,16(sp); ld s11,8(sp)
    addi sp,sp,112; ret

# MATRIX MULTIPLY
# matrix_multiply.s 
# C = A*B
# a0=C a1=A a2=B a3=rows_A a4=cols_A a5=cols_B

.globl matrix_multiply
.type matrix_multiply,@function
matrix_multiply:
    addi sp,sp,-64
    sd ra,56(sp); sd s0,48(sp); sd s1,40(sp)
    sd s2,32(sp); sd s3,24(sp); sd s4,16(sp); sd s5,8(sp)
    mv s0,a0; mv s1,a1; mv s2,a2
    mv s3,a3; mv s4,a4; mv s5,a5
    li t0,0
.L_i: bge t0,s3,.L_end
    li t1,0
.L_j: bge t1,s5,.L_next_i
    fcvt.d.w fs0, zero
    li t2,0
.L_k: bge t2,s4,.L_store
    mul t3,t0,s4; add t3,t3,t2; slli t3,t3,3; add t3,s1,t3; fld ft0,0(t3)
    mul t4,t2,s5; add t4,t4,t1; slli t4,t4,3; add t4,s2,t4; fld ft1,0(t4)
    fmadd.d fs0,ft0,ft1,fs0
    addi t2,t2,1; j .L_k
.L_store:
    mul t3,t0,s5; add t3,t3,t1; slli t3,t3,3; add t3,s0,t3; fsd fs0,0(t3)
    addi t1,t1,1; j .L_j
.L_next_i:
    addi t0,t0,1; j .L_i
.L_end:
    ld ra,56(sp); ld s0,48(sp); ld s1,40(sp)
    ld s2,32(sp); ld s3,24(sp); ld s4,16(sp); ld s5,8(sp)
    addi sp,sp,64
    ret

# CHOLESKY DECOMPOSITION
# Computes L such that A = L * L^T for positive definite A

.data
.align 3
ZERO: .double 0.0

.text
.globl cholesky_decomposition
.type cholesky_decomposition,@function

cholesky_decomposition:
    # Prologue
    addi sp, sp, -96
    sd   ra, 88(sp)
    sd   s0, 80(sp)
    sd   s1, 72(sp)
    sd   s2, 64(sp)
    sd   s3, 56(sp)
    sd   s4, 48(sp)
    sd   s5, 40(sp)
    fsd  fs0, 32(sp)
    fsd  fs1, 24(sp)

    # Save arguments
    mv s0, a0      # s0 = L output matrix
    mv s1, a1      # s1 = A input matrix
    mv s2, a2      # s2 = n (dimension)

    li s3, 0       # i = 0

.Li:
    bge s3, s2, .L_success   # if i >= n, done

    li s4, 0       # j = 0
.Lj:
    bgt s4, s3, .L_next_i    # if j > i, skip to next i

    # sum = A[i,j]
    mul t0, s3, s2
    add t0, t0, s4
    slli t0, t0, 3
    add t0, s1, t0
    fld fs0, 0(t0)

    # Inner loop: sum -= L[i,k]*L[j,k] for k = 0..j-1
    li s5, 0       # k = 0
.Lk:
    bge s5, s4, .Ldiag_off

    # Load L[i,k]
    mul t1, s3, s2
    add t1, t1, s5
    slli t1, t1, 3
    add t1, s0, t1
    fld ft0, 0(t1)

    # Load L[j,k]
    mul t2, s4, s2
    add t2, t2, s5
    slli t2, t2, 3
    add t2, s0, t2
    fld ft1, 0(t2)

    # fs0 = fs0 - L[i,k]*L[j,k]
    fnmsub.d fs0, ft0, ft1, fs0

    addi s5, s5, 1
    j .Lk

.Ldiag_off:
    beq s3, s4, .Ldiag  # if i==j, diagonal

    # Off-diagonal element: L[i,j] = sum / L[j,j]
    mul t1, s4, s2
    add t1, t1, s4
    slli t1, t1, 3
    add t1, s0, t1
    fld fs1, 0(t1)      # fs1 = L[j,j]

    fdiv.d fs0, fs0, fs1

    # Store L[i,j]
    mul t1, s3, s2
    add t1, t1, s4
    slli t1, t1, 3
    add t1, s0, t1
    fsd fs0, 0(t1)

    addi s4, s4, 1
    j .Lj

.Ldiag:
    # Diagonal element: L[i,i] = sqrt(sum), check positive
    la t0, ZERO
    fld fs1, 0(t0)        # load 0.0 into fs1
    flt.d t1, fs0, fs1    # t1 = 1 if fs0 < 0.0
    bnez t1, .Lfail

    fsqrt.d fs0, fs0

    # Store L[i,i]
    mul t1, s3, s2
    add t1, t1, s4
    slli t1, t1, 3
    add t1, s0, t1
    fsd fs0, 0(t1)

    addi s4, s4, 1
    j .Lj

.L_next_i:
    # Zero upper triangle L[i,j] for j>i
    addi s4, s3, 1
.Lzero:
    bge s4, s2, .Linc_i
    mul t1, s3, s2
    add t1, t1, s4
    slli t1, t1, 3
    add t1, s0, t1
    la t0, ZERO
    fld fs0, 0(t0)
    fsd fs0, 0(t1)
    addi s4, s4, 1
    j .Lzero

.Linc_i:
    addi s3, s3, 1
    j .Li

.L_success:
    li a0, 0
    j .L_epilogue

.Lfail:
    li a0, -1

.L_epilogue:
    ld ra, 88(sp)
    ld s0, 80(sp)
    ld s1, 72(sp)
    ld s2, 64(sp)
    ld s3, 56(sp)
    ld s4, 48(sp)
    ld s5, 40(sp)
    fld fs0, 32(sp)
    fld fs1, 24(sp)
    addi sp, sp, 96
    ret

.size cholesky_decomposition, .-cholesky_decomposition

# TRIANGULAR SOLVE (Forward and Backward Substitution)
# triangular_solve.s - Forward and backward substitution for triangular systems

.text

# forward_substitution: Solve Ly = b for lower triangular L
# Arguments:
#   a0 = y (output vector, length n)
#   a1 = L (lower triangular matrix n×n)
#   a2 = b (input vector, length n)
#   a3 = n (dimension)
#
# Algorithm:
#   for i = 0 to n-1:
#     sum = 0
#     for k = 0 to i-1:
#       sum += L[i,k] * y[k]
#     y[i] = (b[i] - sum) / L[i,i]
#
# Register allocation:
#   s0 = y base
#   s1 = L base
#   s2 = b base
#   s3 = n
#   s4 = i (outer loop)
#   s5 = k (inner loop)
#   fs0 = sum accumulator
#   fs1 = L[i,i]

.globl forward_substitution
.type forward_substitution, @function

forward_substitution:
    # Prologue
    addi sp, sp, -64
    sd   ra, 56(sp)
    sd   s0, 48(sp)
    sd   s1, 40(sp)
    sd   s2, 32(sp)
    sd   s3, 24(sp)
    sd   s4, 16(sp)
    sd   s5, 8(sp)
    fsd  fs0, 0(sp)
    
    # Save arguments
    mv   s0, a0              # s0 = y
    mv   s1, a1              # s1 = L
    mv   s2, a2              # s2 = b
    mv   s3, a3              # s3 = n
    
    # Outer loop: i = 0 to n-1
    li   s4, 0
    
.L_fwd_loop_i:
    bge  s4, s3, .L_fwd_done
    
    # Initialize sum = 0.0
    fcvt.d.w fs0, zero
    
    # Inner loop: k = 0 to i-1
    li   s5, 0
    
.L_fwd_loop_k:
    bge  s5, s4, .L_fwd_end_k
    
    # Load L[i,k]
    mul  t0, s4, s3          # t0 = i * n
    add  t0, t0, s5          # t0 = i*n + k
    slli t0, t0, 3           # t0 = (i*n + k) * 8
    add  t0, s1, t0          # t0 = &L[i,k]
    fld  ft0, 0(t0)          # ft0 = L[i,k]
    
    # Load y[k]
    slli t1, s5, 3           # t1 = k * 8
    add  t1, s0, t1          # t1 = &y[k]
    fld  ft1, 0(t1)          # ft1 = y[k]
    
    # sum += L[i,k] * y[k]
    fmadd.d fs0, ft0, ft1, fs0
    
    addi s5, s5, 1
    j    .L_fwd_loop_k
    
.L_fwd_end_k:
    # Load b[i]
    slli t0, s4, 3           # t0 = i * 8
    add  t0, s2, t0          # t0 = &b[i]
    fld  ft0, 0(t0)          # ft0 = b[i]
    
    # ft0 = b[i] - sum
    fsub.d ft0, ft0, fs0
    
    # Load L[i,i]
    mul  t0, s4, s3          # t0 = i * n
    add  t0, t0, s4          # t0 = i*n + i
    slli t0, t0, 3
    add  t0, s1, t0
    fld  ft1, 0(t0)          # ft1 = L[i,i]
    
    # y[i] = (b[i] - sum) / L[i,i]
    fdiv.d ft0, ft0, ft1
    
    # Store y[i]
    slli t0, s4, 3
    add  t0, s0, t0
    fsd  ft0, 0(t0)
    
    addi s4, s4, 1
    j    .L_fwd_loop_i
    
.L_fwd_done:
    # Epilogue
    ld   ra, 56(sp)
    ld   s0, 48(sp)
    ld   s1, 40(sp)
    ld   s2, 32(sp)
    ld   s3, 24(sp)
    ld   s4, 16(sp)
    ld   s5, 8(sp)
    fld  fs0, 0(sp)
    addi sp, sp, 64
    ret

.size forward_substitution, .-forward_substitution

# backward_substitution: Solve L^T x = y for upper triangular L^T
# Arguments:
#   a0 = x (output vector, length n)
#   a1 = L (lower triangular matrix n×n, we solve with L^T)
#   a2 = y (input vector, length n)
#   a3 = n (dimension)
#
# Algorithm:
#   for i = n-1 down to 0:
#     sum = 0
#     for k = i+1 to n-1:
#       sum += L[k,i] * x[k]
#     x[i] = (y[i] - sum) / L[i,i]
#
# Register allocation:
#   s0 = x base
#   s1 = L base
#   s2 = y base
#   s3 = n
#   s4 = i (outer loop, counting down)
#   s5 = k (inner loop)
#   fs0 = sum accumulator

.globl backward_substitution
.type backward_substitution, @function

backward_substitution:
    # Prologue
    addi sp, sp, -64
    sd   ra, 56(sp)
    sd   s0, 48(sp)
    sd   s1, 40(sp)
    sd   s2, 32(sp)
    sd   s3, 24(sp)
    sd   s4, 16(sp)
    sd   s5, 8(sp)
    fsd  fs0, 0(sp)
    
    # Save arguments
    mv   s0, a0              # s0 = x
    mv   s1, a1              # s1 = L
    mv   s2, a2              # s2 = y
    mv   s3, a3              # s3 = n
    
    # Outer loop: i = n-1 down to 0
    addi s4, s3, -1          # s4 = n - 1
    
.L_bwd_loop_i:
    bltz s4, .L_bwd_done     # if i < 0, done
    
    # Initialize sum = 0.0
    fcvt.d.w fs0, zero
    
    # Inner loop: k = i+1 to n-1
    addi s5, s4, 1           # s5 = i + 1
    
.L_bwd_loop_k:
    bge  s5, s3, .L_bwd_end_k
    
    # Load L[k,i] (note: L^T[i,k] = L[k,i])
    mul  t0, s5, s3          # t0 = k * n
    add  t0, t0, s4          # t0 = k*n + i
    slli t0, t0, 3
    add  t0, s1, t0
    fld  ft0, 0(t0)          # ft0 = L[k,i]
    
    # Load x[k]
    slli t1, s5, 3
    add  t1, s0, t1
    fld  ft1, 0(t1)          # ft1 = x[k]
    
    # sum += L[k,i] * x[k]
    fmadd.d fs0, ft0, ft1, fs0
    
    addi s5, s5, 1
    j    .L_bwd_loop_k
    
.L_bwd_end_k:
    # Load y[i]
    slli t0, s4, 3
    add  t0, s2, t0
    fld  ft0, 0(t0)          # ft0 = y[i]
    
    # ft0 = y[i] - sum
    fsub.d ft0, ft0, fs0
    
    # Load L[i,i]
    mul  t0, s4, s3
    add  t0, t0, s4
    slli t0, t0, 3
    add  t0, s1, t0
    fld  ft1, 0(t0)          # ft1 = L[i,i]
    
    # x[i] = (y[i] - sum) / L[i,i]
    fdiv.d ft0, ft0, ft1
    
    # Store x[i]
    slli t0, s4, 3
    add  t0, s0, t0
    fsd  ft0, 0(t0)
    
    addi s4, s4, -1
    j    .L_bwd_loop_i
    
.L_bwd_done:
    # Epilogue
    ld   ra, 56(sp)
    ld   s0, 48(sp)
    ld   s1, 40(sp)
    ld   s2, 32(sp)
    ld   s3, 24(sp)
    ld   s4, 16(sp)
    ld   s5, 8(sp)
    fld  fs0, 0(sp)
    addi sp, sp, 64
    ret

.size backward_substitution, .-backward_substitution

# MATRIX BUILDERS (buildF, buildQ, buildH)
# build_matrices.s - Construct system matrices F, Q, H

.section .rodata
.align 3

.section .text

# buildF: Construct 276×276 state transition matrix
# Arguments:
#   a0 = F (output 276×276 matrix)
#   fa0 = dt (time step)

.globl buildF
.type buildF, @function

buildF:
    # Prologue
    addi sp, sp, -64
    sd   ra, 56(sp)
    sd   s0, 48(sp)
    sd   s1, 40(sp)
    sd   s2, 32(sp)
    fsd  fs0, 24(sp)
    fsd  fs1, 16(sp)
    fsd  fs2, 8(sp)
    fsd  fs3, 0(sp)
    
    mv   s0, a0              # s0 = F base
    fmv.d fs0, fa0           # fs0 = dt
    
    # Compute powers of dt
    fmul.d fs1, fs0, fs0     # fs1 = dt²
    fmul.d fs2, fs1, fs0     # fs2 = dt³
    
    # Compute F_axis entries
    li t0, 2
    fcvt.d.w ft0, t0
    fdiv.d fs3, fs1, ft0     # fs3 = dt²/2
    
    li t0, 6
    fcvt.d.w ft0, t0
    fdiv.d ft1, fs2, ft0     # ft1 = dt³/6 (store in fs4 instead)
    fmv.d fs4, ft1           # fs4 = dt³/6
    
    # Zero entire F
    li t0, 276
    mul t0, t0, t0
    slli t0, t0, 3
    add t1, s0, t0
    mv  t2, s0
    fcvt.d.w ft0, zero
    
.L_buildF_zero:
    bge  t2, t1, .L_buildF_zero_done
    fsd  ft0, 0(t2)
    addi t2, t2, 8
    j    .L_buildF_zero
    
.L_buildF_zero_done:
    # Loop over 23 joints
    li s1, 0
    
.L_buildF_joint:
    li t0, 23
    bge s1, t0, .L_buildF_done
    
    # Loop over 3 axes
    li s2, 0
    
.L_buildF_axis:
    li t0, 3
    bge s2, t0, .L_buildF_next_joint
    
    # Compute base index: joint*12 + axis*4
    li t0, 12
    mul t1, s1, t0
    slli t2, s2, 2
    add t1, t1, t2
    
    # Set F_axis block at [base:base+4, base:base+4]
    # Row 0: [1, dt, dt²/2, dt³/6]
    li t0, 276
    mul t2, t1, t0
    add t2, t2, t1
    slli t2, t2, 3
    add t2, s0, t2
    
    fcvt.d.w ft0, zero
    li t3, 1
    fcvt.d.w ft1, t3
    fsd ft1, 0(t2)           # F[base,base] = 1
    fsd fs0, 8(t2)           # F[base,base+1] = dt
    fsd fs3, 16(t2)          # F[base,base+2] = dt²/2
    fsd fs4, 24(t2)          # F[base,base+3] = dt³/6
    
    # Row 1: [0, 1, dt, dt²/2]
    # FIXED: Use li + add instead of addi for large immediate
    li t3, 2208              # 276*8 bytes
    add t2, t2, t3           # Next row
    fsd ft0, 0(t2)
    fsd ft1, 8(t2)
    fsd fs0, 16(t2)
    fsd fs3, 24(t2)
    
    # Row 2: [0, 0, 1, dt]
    li t3, 2208
    add t2, t2, t3
    fsd ft0, 0(t2)
    fsd ft0, 8(t2)
    fsd ft1, 16(t2)
    fsd fs0, 24(t2)
    
    # Row 3: [0, 0, 0, 1]
    li t3, 2208
    add t2, t2, t3
    fsd ft0, 0(t2)
    fsd ft0, 8(t2)
    fsd ft0, 16(t2)
    fsd ft1, 24(t2)
    
    addi s2, s2, 1
    j    .L_buildF_axis
    
.L_buildF_next_joint:
    addi s1, s1, 1
    j    .L_buildF_joint
    
.L_buildF_done:
    # Epilogue
    ld   ra, 56(sp)
    ld   s0, 48(sp)
    ld   s1, 40(sp)
    ld   s2, 32(sp)
    fld  fs0, 24(sp)
    fld  fs1, 16(sp)
    fld  fs2, 8(sp)
    fld  fs3, 0(sp)
    addi sp, sp, 64
    ret

.size buildF, .-buildF

# buildQ: Construct 276×276 process noise covariance matrix
# Arguments:
#   a0 = Q (output 276×276 matrix)
#   fa0 = dt
#   fa1 = sigma (jerk noise std dev)

.globl buildQ
.type buildQ, @function

buildQ:
    # Prologue - need more saved FP registers
    addi sp, sp, -96
    sd   ra, 88(sp)
    sd   s0, 80(sp)
    sd   s1, 72(sp)
    sd   s2, 64(sp)
    fsd  fs0, 56(sp)
    fsd  fs1, 48(sp)
    fsd  fs2, 40(sp)
    fsd  fs3, 32(sp)
    fsd  fs4, 24(sp)
    fsd  fs5, 16(sp)
    fsd  fs6, 8(sp)
    fsd  fs7, 0(sp)
    
    mv   s0, a0
    fmv.d fs0, fa0           # dt
    fmv.d fs1, fa1           # sigma
    
    # Compute dt powers
    fmul.d fs2, fs0, fs0     # dt²
    fmul.d ft0, fs2, fs0     # dt³
    fmul.d ft1, ft0, fs0     # dt⁴
    fmul.d ft2, ft1, fs0     # dt⁵
    fmul.d ft3, ft2, fs0     # dt⁶
    
    # Compute sigma²
    fmul.d fs3, fs1, fs1
    
    # Zero entire Q first
    li t0, 276
    mul t0, t0, t0
    slli t0, t0, 3
    add t1, s0, t0
    mv  t2, s0
    fcvt.d.w ft4, zero
    
.L_buildQ_zero:
    bge  t2, t1, .L_buildQ_zero_done
    fsd  ft4, 0(t2)
    addi t2, t2, 8
    j    .L_buildQ_zero
    
.L_buildQ_zero_done:
    # Loop over joints and axes
    li s1, 0
    
.L_buildQ_joint:
    li t0, 23
    bge s1, t0, .L_buildQ_done
    
    li s2, 0
    
.L_buildQ_axis:
    li t0, 3
    bge s2, t0, .L_buildQ_next_joint
    
    # Compute base
    li t0, 12
    mul t1, s1, t0
    slli t2, s2, 2
    add t1, t1, t2
    
    # Compute all Q_axis entries
    # Row 0
    li t0, 36
    fcvt.d.w ft4, t0
    fdiv.d ft5, ft3, ft4
    fmul.d fs4, ft5, fs3     # Q00 = dt⁶/36 * σ²
    
    li t0, 12
    fcvt.d.w ft4, t0
    fdiv.d ft6, ft2, ft4
    fmul.d fs5, ft6, fs3     # Q01 = dt⁵/12 * σ²
    
    li t0, 6
    fcvt.d.w ft4, t0
    fdiv.d ft7, ft1, ft4
    fmul.d fs6, ft7, fs3     # Q02 = dt⁴/6 * σ²
    
    fdiv.d ft8, ft0, ft4
    fmul.d fs7, ft8, fs3     # Q03 = dt³/6 * σ²
    
    # Row 1
    li t0, 4
    fcvt.d.w ft4, t0
    fdiv.d ft9, ft1, ft4
    fmul.d fa2, ft9, fs3     # Q11 = dt⁴/4 * σ²
    
    li t0, 2
    fcvt.d.w ft4, t0
    fdiv.d ft10, ft0, ft4
    fmul.d fa3, ft10, fs3    # Q12 = dt³/2 * σ²
    
    fdiv.d ft11, fs2, ft4
    fmul.d fa4, ft11, fs3    # Q13 = dt²/2 * σ²
    
    # Row 2
    fmul.d fa5, fs2, fs3     # Q22 = dt² * σ²
    fmul.d fa6, fs0, fs3     # Q23 = dt * σ²
    
    # Row 3
    fmv.d fa7, fs3           # Q33 = σ²
    
    # Calculate address and store
    li t0, 276
    mul t2, t1, t0
    add t2, t2, t1
    slli t2, t2, 3
    add t2, s0, t2
    
    # Store row 0
    fsd fs4, 0(t2)
    fsd fs5, 8(t2)
    fsd fs6, 16(t2)
    fsd fs7, 24(t2)
    
    # Store row 1 
    li t3, 2208
    add t2, t2, t3
    fsd fs5, 0(t2)
    fsd fa2, 8(t2)
    fsd fa3, 16(t2)
    fsd fa4, 24(t2)
    
    # Store row 2
    li t3, 2208
    add t2, t2, t3
    fsd fs6, 0(t2)
    fsd fa3, 8(t2)
    fsd fa5, 16(t2)
    fsd fa6, 24(t2)
    
    # Store row 3 
    li t3, 2208
    add t2, t2, t3
    fsd fs7, 0(t2)
    fsd fa4, 8(t2)
    fsd fa6, 16(t2)
    fsd fa7, 24(t2)
    
    addi s2, s2, 1
    j    .L_buildQ_axis
    
.L_buildQ_next_joint:
    addi s1, s1, 1
    j    .L_buildQ_joint
    
.L_buildQ_done:
    # Epilogue
    ld   ra, 88(sp)
    ld   s0, 80(sp)
    ld   s1, 72(sp)
    ld   s2, 64(sp)
    fld  fs0, 56(sp)
    fld  fs1, 48(sp)
    fld  fs2, 40(sp)
    fld  fs3, 32(sp)
    fld  fs4, 24(sp)
    fld  fs5, 16(sp)
    fld  fs6, 8(sp)
    fld  fs7, 0(sp)
    addi sp, sp, 96
    ret

.size buildQ, .-buildQ

# buildH: Construct 69×276 measurement matrix (sparse)
# Arguments:
#   a0 = H (output 69×276 matrix)

.globl buildH
.type buildH, @function

buildH:
    # Prologue
    addi sp, sp, -32
    sd   ra, 24(sp)
    sd   s0, 16(sp)
    sd   s1, 8(sp)
    sd   s2, 0(sp)
    
    mv   s0, a0
    
    # Zero entire H
    li t0, 69
    li t1, 276
    mul t0, t0, t1
    slli t0, t0, 3
    add t1, s0, t0
    mv  t2, s0
    fcvt.d.w ft0, zero
    
.L_buildH_zero:
    bge  t2, t1, .L_buildH_zero_done
    fsd  ft0, 0(t2)
    addi t2, t2, 8
    j    .L_buildH_zero
    
.L_buildH_zero_done:
    # Set ones: H[j*3+a, j*12+a*4] = 1.0
    li t0, 1
    fcvt.d.w ft1, t0
    
    li s1, 0
    
.L_buildH_joint:
    li t0, 23
    bge s1, t0, .L_buildH_done
    
    li s2, 0
    
.L_buildH_axis:
    li t0, 3
    bge s2, t0, .L_buildH_next_joint
    
    # Row index: j*3 + axis
    li t0, 3
    mul t1, s1, t0
    add t1, t1, s2
    
    # Col index: j*12 + axis*4
    li t0, 12
    mul t2, s1, t0
    slli t3, s2, 2
    add t2, t2, t3
    
    # Address: H + (row*276 + col)*8
    li t0, 276
    mul t3, t1, t0
    add t3, t3, t2
    slli t3, t3, 3
    add t3, s0, t3
    
    # H[row,col] = 1.0
    fsd ft1, 0(t3)
    
    addi s2, s2, 1
    j    .L_buildH_axis
    
.L_buildH_next_joint:
    addi s1, s1, 1
    j    .L_buildH_joint
    
.L_buildH_done:
    # Epilogue
    ld   ra, 24(sp)
    ld   s0, 16(sp)
    ld   s1, 8(sp)
    ld   s2, 0(sp)
    addi sp, sp, 32
    ret

.size buildH, .-buildH

# MATRIX OPERATIONS (add, subtract, transpose, etc.)
#   matrix_add, matrix_subtract, matrix_copy, matrix_zero,
#   matrix_transpose, identity_minus, init_identity,
#   vector_add_inplace, chol_solve_matrix
.section .text

# matrix_add(C,A,B,n): C[i]=A[i]+B[i]  a0=C a1=A a2=B a3=n
.globl matrix_add
.align 2
matrix_add:
    beqz a3,.Lma_d; li t0,0
.Lma_l: slli t1,t0,3
    add t2,a1,t1; fld ft0,0(t2)
    add t2,a2,t1; fld ft1,0(t2)
    fadd.d ft0,ft0,ft1
    add t2,a0,t1; fsd ft0,0(t2)
    addi t0,t0,1; blt t0,a3,.Lma_l
.Lma_d: ret

# matrix_subtract(C,A,B,n): C[i]=A[i]-B[i]
.globl matrix_subtract
.align 2
matrix_subtract:
    beqz a3,.Lms_d; li t0,0
.Lms_l: slli t1,t0,3
    add t2,a1,t1; fld ft0,0(t2)
    add t2,a2,t1; fld ft1,0(t2)
    fsub.d ft0,ft0,ft1
    add t2,a0,t1; fsd ft0,0(t2)
    addi t0,t0,1; blt t0,a3,.Lms_l
.Lms_d: ret

# matrix_copy(dst,src,n): copy n doubles
.globl matrix_copy
.align 2
matrix_copy:
    beqz a2,.Lmc_d; li t0,0
.Lmc_l: slli t1,t0,3
    add t2,a1,t1; ld t3,0(t2)
    add t2,a0,t1; sd t3,0(t2)
    addi t0,t0,1; blt t0,a2,.Lmc_l
.Lmc_d: ret

# matrix_zero(A,n): zero n doubles
.globl matrix_zero
.align 2
matrix_zero:
    beqz a1,.Lmz_d; li t0,0
.Lmz_l: slli t1,t0,3
    add t2,a0,t1; sd zero,0(t2)
    addi t0,t0,1; blt t0,a1,.Lmz_l
.Lmz_d: ret

# matrix_transpose(At,A,rows,cols): At[j*rows+i]=A[i*cols+j]
# a0=At a1=A a2=rows a3=cols
.globl matrix_transpose
.align 2
matrix_transpose:
    addi sp,sp,-56
    sd ra,48(sp); sd s0,40(sp); sd s1,32(sp)
    sd s2,24(sp); sd s3,16(sp); sd s4,8(sp); sd s5,0(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    li s4,0
.Lmt_i: bge s4,s2,.Lmt_d
    li s5,0
.Lmt_j: bge s5,s3,.Lmt_ni
    mul t0,s4,s3; add t0,t0,s5; slli t0,t0,3; add t0,s1,t0; ld t1,0(t0)
    mul t2,s5,s2; add t2,t2,s4; slli t2,t2,3; add t2,s0,t2; sd t1,0(t2)
    addi s5,s5,1; j .Lmt_j
.Lmt_ni: addi s4,s4,1; j .Lmt_i
.Lmt_d:
    ld ra,48(sp); ld s0,40(sp); ld s1,32(sp)
    ld s2,24(sp); ld s3,16(sp); ld s4,8(sp); ld s5,0(sp)
    addi sp,sp,56; ret

# identity_minus(C,A,n): C = I_n - A  (n×n matrices)
# a0=C a1=A a2=n
.globl identity_minus
.align 2
identity_minus:
    addi sp,sp,-40
    sd ra,32(sp); sd s0,24(sp); sd s1,16(sp); sd s2,8(sp); sd s3,0(sp)
    mv s0,a0; mv s1,a1; mv s2,a2
    li t2,1; fcvt.d.w ft2,t2          # ft2 = 1.0
    li s3,0
.Lim_i: bge s3,s2,.Lim_d
    li t6,0
.Lim_j: bge t6,s2,.Lim_ni
    mul t0,s3,s2; add t0,t0,t6; slli t0,t0,3
    add t1,s1,t0; fld ft0,0(t1)
    bne s3,t6,.Lim_off
    fsub.d ft1,ft2,ft0; j .Lim_st
.Lim_off: fneg.d ft1,ft0
.Lim_st: add t1,s0,t0; fsd ft1,0(t1)
    addi t6,t6,1; j .Lim_j
.Lim_ni: addi s3,s3,1; j .Lim_i
.Lim_d:
    ld ra,32(sp); ld s0,24(sp); ld s1,16(sp); ld s2,8(sp); ld s3,0(sp)
    addi sp,sp,40; ret

# init_identity(M,n): M = I_n
.globl init_identity
.align 2
init_identity:
    addi sp,sp,-32; sd ra,24(sp); sd s0,16(sp); sd s1,8(sp)
    mv s0,a0; mv s1,a1
    mul a1,s1,s1; mv a0,s0; call matrix_zero
    li t2,1; fcvt.d.w ft0,t2; li t0,0
.Lii_l: bge t0,s1,.Lii_d
    mul t1,t0,s1; add t1,t1,t0; slli t1,t1,3; add t1,s0,t1; fsd ft0,0(t1)
    addi t0,t0,1; j .Lii_l
.Lii_d:
    ld ra,24(sp); ld s0,16(sp); ld s1,8(sp); addi sp,sp,32; ret

# vector_add_inplace(a,b,n): a[i]+=b[i]
.globl vector_add_inplace
.align 2
vector_add_inplace:
    beqz a2,.Lvai_d; li t0,0
.Lvai_l: slli t1,t0,3
    add t2,a0,t1; fld ft0,0(t2)
    add t3,a1,t1; fld ft1,0(t3)
    fadd.d ft0,ft0,ft1; fsd ft0,0(t2)
    addi t0,t0,1; blt t0,a2,.Lvai_l
.Lvai_d: ret

# chol_solve_matrix(X,L,B,n,m): solve LL^T * X = B  column by column
# a0=X(n×m output)  a1=L(n×n)  a2=B(n×m)  a3=n  a4=m
# Uses forward_substitution + backward_substitution from triangular_solve.s
.globl chol_solve_matrix
.align 2
chol_solve_matrix:
    addi sp,sp,-80
    sd ra,72(sp); sd s0,64(sp); sd s1,56(sp); sd s2,48(sp)
    sd s3,40(sp); sd s4,32(sp); sd s5,24(sp); sd s6,16(sp); sd s7,8(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3; mv s4,a4
    # Alloc two n-element column buffers
    slli a0,s3,3; call malloc; mv s5,a0   # s5=y_buf (fwd output)
    slli a0,s3,3; call malloc; mv s6,a0   # s6=x_buf (bwd output)
    li s7,0                                # col=0
.Lcsm_col: bge s7,s4,.Lcsm_done
    # Extract B[:,col] into s5 (used as b_col for forward solve)
    li t0,0
.Lcsm_ex: bge t0,s3,.Lcsm_fwd
    mul t1,t0,s4; add t1,t1,s7; slli t1,t1,3; add t1,s2,t1; fld ft0,0(t1)
    slli t2,t0,3; add t2,s5,t2; fsd ft0,0(t2)
    addi t0,t0,1; j .Lcsm_ex
.Lcsm_fwd:
    # forward_substitution(y=s6, L=s1, b=s5, n=s3)
    mv a0,s6; mv a1,s1; mv a2,s5; mv a3,s3; call forward_substitution
    # backward_substitution(x=s5, L=s1, y=s6, n=s3)
    mv a0,s5; mv a1,s1; mv a2,s6; mv a3,s3; call backward_substitution
    # Write s5 back to X[:,col]
    li t0,0
.Lcsm_wr: bge t0,s3,.Lcsm_nc
    slli t1,t0,3; add t1,s5,t1; fld ft0,0(t1)
    mul t2,t0,s4; add t2,t2,s7; slli t2,t2,3; add t2,s0,t2; fsd ft0,0(t2)
    addi t0,t0,1; j .Lcsm_wr
.Lcsm_nc: addi s7,s7,1; j .Lcsm_col
.Lcsm_done:
    mv a0,s5; call free; mv a0,s6; call free
    ld ra,72(sp); ld s0,64(sp); ld s1,56(sp); ld s2,48(sp)
    ld s3,40(sp); ld s4,32(sp); ld s5,24(sp); ld s6,16(sp); ld s7,8(sp)
    addi sp,sp,80; ret
