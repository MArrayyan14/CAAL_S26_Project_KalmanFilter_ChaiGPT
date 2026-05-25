# ekf_vector.s — FULLY CORRECTED Extended Kalman Filter (RISC-V RVV)
# Team ChaiGPT - Milestone 4

.option arch, +v
.equ WS_xtmp,   0
.equ WS_FP,     2208
.equ WS_FPFt,   611616
.equ WS_Ft,     1221024
.equ WS_H,      1830432
.equ WS_Ht,     1982784
.equ WS_PHt,    2135136
.equ WS_PHtT,   2287488
.equ WS_Kt,     2439840
.equ WS_K,      2592192
.equ WS_S,      2744544
.equ WS_L,      2782632
.equ WS_y,      2820720
.equ WS_Hx,     2821272
.equ WS_Ky,     2821824
.equ WS_KH,     2824032
.equ WS_IKH,    3433440
.equ WS_IKHT,   4042848
.equ WS_IKH_P,  4652256
.equ WS_KR,     5261664
.equ WS_KTJ,    5414016
.equ WS_KRKt,   5566368
.equ EKF_WS_TOTAL, 6175776

.equ NN, 76176
.equ NM, 19044
.equ MM, 4761

.section .rodata
    .align 3
ZERO_D: .double 0.0 
PI_val:        .double 3.14159265358979323846
PI_half:       .double 1.57079632679489661923
PI_const:      .double 3.14159265358979323846
TWO_PI_const:  .double 6.28318530717958647692
EPS_R:         .double 1e-3
EPS_RXY:       .double 0.01
SKIP_THRESH:   .double 1e6
c0:            .double 0.999866
c1:            .double 0.3302995
c2:            .double 0.180141
c3:            .double 0.085133
c4:            .double 0.0208351

.section .text

# arctan2_manual (SCALAR ,unchanged from M3)
.globl arctan2_manual
.align 2
arctan2_manual:
    addi sp,sp,-80
    sd ra,72(sp)
    fsd fs0,64(sp); fsd fs1,56(sp); fsd fs2,48(sp)
    fsd fs3,40(sp); fsd fs4,32(sp); fsd fs5,24(sp)
    
    fmv.d fs0,fa0; fmv.d fs1,fa1
    fabs.d fs2,fs0; fabs.d fs3,fs1
    
    fcvt.d.w ft0,zero
    feq.d t0,fs2,ft0; feq.d t1,fs3,ft0; and t0,t0,t1
    bnez t0,.Lreturn_zero
    
.Lcontinue:
    flt.d t0,fs3,fs2; bnez t0,.Ly_bigger
    fdiv.d fs4,fs2,fs3; li t3,0; j .Lcalc
.Ly_bigger:
    fdiv.d fs4,fs3,fs2; li t3,1
    
.Lcalc:
    fmv.d fs5,fs4
    fmul.d ft0,fs4,fs4
    
    la t2,c4; fld ft1,0(t2)
    la t2,c3; fld ft2,0(t2); fmul.d ft3,ft0,ft1; fsub.d ft1,ft2,ft3
    la t2,c2; fld ft2,0(t2); fmul.d ft3,ft0,ft1; fsub.d ft1,ft2,ft3
    la t2,c1; fld ft2,0(t2); fmul.d ft3,ft0,ft1; fsub.d ft1,ft2,ft3
    la t2,c0; fld ft2,0(t2); fmul.d ft3,ft0,ft1; fsub.d ft1,ft2,ft3
    
    fmul.d fs4,fs5,ft1
    
    beqz t3,.Lquadrant
    la t2,PI_half; fld ft0,0(t2); fsub.d fs4,ft0,fs4
    
.Lquadrant:
    fcvt.d.w ft0,zero; flt.d t0,fs1,ft0; beqz t0,.Lsigny
    la t2,PI_val; fld ft0,0(t2); fsub.d fs4,ft0,fs4
    
.Lsigny:
    fcvt.d.w ft0,zero; flt.d t0,fs0,ft0; beqz t0,.Ldone
    fneg.d fs4,fs4
    
.Ldone:
    fmv.d fa0,fs4; j .Lepilog
    
.Lreturn_zero:
    fcvt.d.w fa0,zero
    
.Lepilog:
    ld ra,72(sp)
    fld fs0,64(sp); fld fs1,56(sp); fld fs2,48(sp)
    fld fs3,40(sp); fld fs4,32(sp); fld fs5,24(sp)
    addi sp,sp,80
    ret

# run_ekf_asm

.globl run_ekf_asm
.align 2
run_ekf_asm:
    addi sp,sp,-144
    sd ra,136(sp); sd s0,128(sp); sd s1,120(sp); sd s2,112(sp)
    sd s3,104(sp); sd s4,96(sp);  sd s5,88(sp);  sd s6,80(sp)
    sd s7,72(sp);  sd s8,64(sp);  sd s9,56(sp);  sd s10,48(sp); sd s11,40(sp)
    
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    mv s4,a4; mv s5,a5; mv s6,a6

    li a0,EKF_WS_TOTAL; call malloc; mv s7,a0
    li a0,2208;          call malloc; mv s8,a0
    li a0,609408;        call malloc; mv s9,a0
    li a0,152352;        call malloc; mv s11,a0

    mv a0,s8; li a1,276; call vekf_zero

    li t0,0
.Lrekf_ix:
    li t1,23; bge t0,t1,.Lrekf_ix_d
    li t2,24; mul t3,t0,t2; add t3,s0,t3
    li t2,96; mul t4,t0,t2; add t4,s8,t4
    fld ft0,0(t3);  fsd ft0,0(t4)
    fld ft1,8(t3);  fsd ft1,32(t4)
    fld ft2,16(t3); fsd ft2,64(t4)
    addi t0,t0,1; j .Lrekf_ix
.Lrekf_ix_d:

    mv a0,s9; li a1,276; call init_identity

    li s10,0
.Lrekf_loop:
    bge s10,s2,.Lrekf_done
    li t0,552; mul t0,s10,t0; add t0,s1,t0
    mv a0,s8; mv a1,s9; mv a2,t0
    mv a3,s3; mv a4,s4; mv a5,s5; mv a6,s11; mv a7,s7
    call ekf_step
    li t0,2208; mul t0,s10,t0; add t0,s6,t0
    mv a0,t0; mv a1,s8; li a2,276; call vekf_copy
    addi s10,s10,1; j .Lrekf_loop

.Lrekf_done:
    mv a0,s7; call free
    mv a0,s8; call free
    mv a0,s9; call free
    mv a0,s11; call free
    
    ld ra,136(sp); ld s0,128(sp); ld s1,120(sp); ld s2,112(sp)
    ld s3,104(sp); ld s4,96(sp);  ld s5,88(sp);  ld s6,80(sp)
    ld s7,72(sp);  ld s8,64(sp);  ld s9,56(sp);  ld s10,48(sp); ld s11,40(sp)
    addi sp,sp,144
    ret

# ekf_step

.globl ekf_step
.align 2
ekf_step:
    addi sp,sp,-112
    sd ra,104(sp); sd s0,96(sp);  sd s1,88(sp);  sd s2,80(sp)
    sd s3,72(sp);  sd s4,64(sp);  sd s5,56(sp);  sd s6,48(sp)
    sd s7,40(sp);  sd s8,32(sp);  sd s9,24(sp);  sd s10,16(sp); sd s11,8(sp)
    
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    mv s4,a4; mv s5,a5; mv s6,a6; mv s7,a7

    # PREDICT 
    li t0,WS_xtmp; add s8,s7,t0
    mv a0,s8; mv a1,s3; mv a2,s0; li a3,276; li a4,276; li a5,1
    call vekf_matmul
    mv a0,s0; mv a1,s8; li a2,276; call vekf_copy

    li t0,WS_FP; add s8,s7,t0
    mv a0,s8; mv a1,s3; mv a2,s1; li a3,276; li a4,276; li a5,276
    call vekf_matmul

    li t0,WS_Ft; add s9,s7,t0
    mv a0,s9; mv a1,s3; li a2,276; li a3,276
    call vekf_transpose

    li t0,WS_FPFt; add s10,s7,t0
    mv a0,s10; mv a1,s8; mv a2,s9; li a3,276; li a4,276; li a5,276
    call vekf_matmul

    mv a0,s1; mv a1,s10; mv a2,s4; li a3,NN
    call vekf_add

    # UPDATE 
    li t0,WS_H; add s8,s7,t0
    mv a0,s8; mv a1,s0; call computeJacobian

    li t0,WS_Hx; add s9,s7,t0
    mv a0,s9; mv a1,s0; call h_func

    li t0,WS_y; add s10,s7,t0
    mv a0,s10; mv a1,s2; mv a2,s9; li a3,69
    call vekf_sub

    # Angle wrapping
    la t0,PI_const;     fld fs0,0(t0)
    la t0,TWO_PI_const; fld fs1,0(t0)
    li s11,0
.Lekf_wrap:
    li t0,23; bge s11,t0,.Lekf_wrap_done

    li t0,24; mul t1,s11,t0; addi t1,t1,8
    add t1,s10,t1; fld ft0,0(t1)
    
.Lekf_wt:
    fle.d t2,ft0,fs0; bnez t2,.Lekf_wt_dn
    fsub.d ft0,ft0,fs1; j .Lekf_wt
.Lekf_wt_dn:
    fneg.d ft1,fs0
    fle.d t2,ft1,ft0; bnez t2,.Lekf_wt_ok
    fadd.d ft0,ft0,fs1; j .Lekf_wt_dn
.Lekf_wt_ok:
    fsd ft0,0(t1)

    addi t1,t1,8; fld ft0,0(t1)
    
.Lekf_wp:
    fle.d t2,ft0,fs0; bnez t2,.Lekf_wp_dn
    fsub.d ft0,ft0,fs1; j .Lekf_wp
.Lekf_wp_dn:
    fneg.d ft1,fs0
    fle.d t2,ft1,ft0; bnez t2,.Lekf_wp_ok
    fadd.d ft0,ft0,fs1; j .Lekf_wp_dn
.Lekf_wp_ok:
    fsd ft0,0(t1)
    
    addi s11,s11,1; j .Lekf_wrap
.Lekf_wrap_done:

    mv a0,s10; li a1,69; call vector_norm
    la t0,SKIP_THRESH; fld ft0,0(t0)
    flt.d t0,ft0,fa0
    bnez t0,.Lekf_skip

    li t0,WS_Ht; add s9,s7,t0
    mv a0,s9; mv a1,s8; li a2,69; li a3,276
    call vekf_transpose

    li t0,WS_PHt; add s10,s7,t0
    mv a0,s10; mv a1,s1; mv a2,s9; li a3,276; li a4,276; li a5,69
    call vekf_matmul

    li t0,WS_S; add s9,s7,t0
    mv a0,s9; mv a1,s8; mv a2,s10; li a3,69; li a4,276; li a5,69
    call vekf_matmul
    mv a0,s9; mv a1,s9; mv a2,s5; li a3,MM
    call vekf_add

    li t0,WS_L; add s11,s7,t0
    mv a0,s11; mv a1,s9; li a2,69
    call cholesky_decomposition
    bnez a0,.Lekf_skip

    li t0,WS_PHtT; add s9,s7,t0
    mv a0,s9; mv a1,s10; li a2,276; li a3,69
    call vekf_transpose

    li t0,WS_Kt; add s10,s7,t0
    mv a0,s10; mv a1,s11; mv a2,s9; li a3,69; li a4,276
    call chol_solve_matrix

    li t0,WS_K; add s9,s7,t0
    mv a0,s9; mv a1,s10; li a2,69; li a3,276
    call vekf_transpose

    li t0,WS_Ky; add s10,s7,t0
    li t1,WS_y;  add t1,s7,t1
    mv a0,s10; mv a1,s9; mv a2,t1; li a3,276; li a4,69; li a5,1
    call vekf_matmul

    mv a0,s0; mv a1,s10; li a2,276
    call vekf_add_inplace

    li t0,WS_KH; add s10,s7,t0
    mv a0,s10; mv a1,s9; mv a2,s8; li a3,276; li a4,69; li a5,276
    call vekf_matmul

    li t0,WS_IKH; add s11,s7,t0
    mv a0,s11; mv a1,s10; li a2,276
    call vekf_identity_minus

    li t0,WS_IKH_P; add s10,s7,t0
    mv a0,s10; mv a1,s11; mv a2,s1; li a3,276; li a4,276; li a5,276
    call vekf_matmul

    li t0,WS_IKHT; add s9,s7,t0
    mv a0,s9; mv a1,s11; li a2,276; li a3,276
    call vekf_transpose

    li t0,WS_FPFt; add s11,s7,t0
    mv a0,s11; mv a1,s10; mv a2,s9; li a3,276; li a4,276; li a5,276
    call vekf_matmul

    li t0,WS_KR;  add s10,s7,t0
    li t1,WS_K;   add t1,s7,t1
    mv a0,s10; mv a1,t1; mv a2,s5; li a3,276; li a4,69; li a5,69
    call vekf_matmul

    li t0,WS_KTJ; add s9,s7,t0
    li t1,WS_K;   add t1,s7,t1
    mv a0,s9; mv a1,t1; li a2,276; li a3,69
    call vekf_transpose

    li t0,WS_KRKt; add s8,s7,t0
    mv a0,s8; mv a1,s10; mv a2,s9; li a3,276; li a4,69; li a5,276
    call vekf_matmul

    mv a0,s1; mv a1,s11; mv a2,s8; li a3,NN
    call vekf_add

.Lekf_skip:
    ld ra,104(sp); ld s0,96(sp);  ld s1,88(sp);  ld s2,80(sp)
    ld s3,72(sp);  ld s4,64(sp);  ld s5,56(sp);  ld s6,48(sp)
    ld s7,40(sp);  ld s8,32(sp);  ld s9,24(sp);  ld s10,16(sp); ld s11,8(sp)
    addi sp,sp,112
    ret

# vector_norm (VECTOR)
.globl vector_norm
.align 2
vector_norm:
    addi sp,sp,-16
    sd ra,8(sp); sd s0,0(sp)

    mv s0,a0                   # save vector base
    mv t0,a1                   # remaining = n
    vmv.v.i v0,0               # v0[0] = 0.0  (reduction accumulator)

.Lvn_loop:
    beqz t0,.Lvn_done
    vsetvli t1,t0,e64,m8,ta,ma
    vle64.v v8,(s0)            # load chunk
    vfmul.vv v8,v8,v8          # square each lane
    vfredosum.vs v0,v8,v0      # v0[0] += sum(chunk)
    slli t2,t1,3
    add s0,s0,t2
    sub t0,t0,t1
    j .Lvn_loop

.Lvn_done:
    vfmv.f.s fa0,v0            # move scalar sum from v0[0] to fa0
    fsqrt.d fa0,fa0

    ld ra,8(sp); ld s0,0(sp)
    addi sp,sp,16
    ret

# ═════════════════════════════════════════════════════════════════════════════
# h_func (SCALAR - EXACT COPY FROM M3)
# ═════════════════════════════════════════════════════════════════════════════
.globl h_func
.align 2
h_func:
    addi sp,sp,-80
    sd ra,72(sp); sd s0,64(sp); sd s1,56(sp); sd s2,48(sp)
    fsd fs0,40(sp); fsd fs1,32(sp); fsd fs2,24(sp)
    fsd fs3,16(sp); fsd fs4,8(sp)
    
    mv s0,a0; mv s1,a1
    li s2,0
    
.Lhf_loop:
    li t0,23; bge s2,t0,.Lhf_done
    
    li t0,96; mul t1,s2,t0
    add t1,s1,t1
    fld fs0,0(t1)    # px
    fld fs1,32(t1)   # py
    fld fs2,64(t1)   # pz
    
    # r = sqrt(px^2 + py^2 + pz^2)
    fmul.d ft0,fs0,fs0
    fmul.d ft1,fs1,fs1
    fmadd.d ft0,fs2,fs2,ft0
    fadd.d ft0,ft0,ft1
    fsqrt.d fs3,ft0
    
    # theta = arctan2(py, px)
    fmv.d fa0,fs1
    fmv.d fa1,fs0
    call arctan2_manual
    fmv.d fs4,fa0
    
    # phi = arctan2(pz, sqrt(px^2 + py^2))
    fmul.d ft0,fs0,fs0
    fmul.d ft1,fs1,fs1
    fadd.d ft0,ft0,ft1
    fsqrt.d ft0,ft0
    fmv.d fa0,fs2
    fmv.d fa1,ft0
    call arctan2_manual
    
    # Store r, theta, phi
    li t0,24; mul t1,s2,t0
    add t1,s0,t1
    fsd fs3,0(t1)
    fsd fs4,8(t1)
    fsd fa0,16(t1)
    
    addi s2,s2,1
    j .Lhf_loop
    
.Lhf_done:
    ld ra,72(sp); ld s0,64(sp); ld s1,56(sp); ld s2,48(sp)
    fld fs0,40(sp); fld fs1,32(sp); fld fs2,24(sp)
    fld fs3,16(sp); fld fs4,8(sp)
    addi sp,sp,80
    ret

# computeJacobian 

.globl computeJacobian
.align 2
computeJacobian:
    addi sp,sp,-96
    sd ra,88(sp); sd s0,80(sp); sd s1,72(sp); sd s2,64(sp)
    fsd fs0,56(sp); fsd fs1,48(sp); fsd fs2,40(sp)
    fsd fs3,32(sp); fsd fs4,24(sp); fsd fs5,16(sp); fsd fs6,8(sp)
    
    mv s0,a0; mv s1,a1
    
    # Zero H first
    li a1,19044; mv a0,s0; call vekf_zero
    
    # Load epsilon constants
    la t0,EPS_R;   fld fs4,0(t0)
    la t0,EPS_RXY; fld fs5,0(t0)
    
    li s2,0
.Lcj_loop:
    li t0,23; bge s2,t0,.Lcj_done
    
    # Extract px, py, pz
    li t0,96; mul t1,s2,t0
    add t1,s1,t1
    fld fs0,0(t1)   # px
    fld fs1,32(t1)  # py
    fld fs2,64(t1)  # pz
    
    # rho_xy = sqrt(px^2 + py^2)
    fmul.d ft0,fs0,fs0
    fmul.d ft1,fs1,fs1
    fadd.d ft0,ft0,ft1
    fsqrt.d fs3,ft0
    
    # r = sqrt(rho_xy^2 + pz^2)
    fmul.d ft0,fs2,fs2
    fmadd.d ft0,fs3,fs3,ft0
    fsqrt.d fs6,ft0
    
    # Epsilon protection
    flt.d t0,fs6,fs4
    beqz t0,.Lcj_r_ok
    fmv.d fs6,fs4
.Lcj_r_ok:
    flt.d t0,fs3,fs5
    beqz t0,.Lcj_rxy_ok
    fmv.d fs3,fs5
.Lcj_rxy_ok:
    
    # Compute r2, rxy2
    fmul.d ft4,fs6,fs6
    fmul.d ft5,fs3,fs3
    
    # Row base: joint*3, Col base: joint*12
    li t0,3;  mul t3,s2,t0
    li t0,12; mul t4,s2,t0
    
    # ∂r/∂px = px/r → H[row,col+0]
    fdiv.d ft0,fs0,fs6
    li t0,276; mul t1,t3,t0; add t1,t1,t4; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    # ∂r/∂py = py/r → H[row,col+4]
    fdiv.d ft0,fs1,fs6
    li t0,276; mul t1,t3,t0; add t1,t1,t4; addi t1,t1,4; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    # ∂r/∂pz = pz/r → H[row,col+8]
    fdiv.d ft0,fs2,fs6
    li t0,276; mul t1,t3,t0; add t1,t1,t4; addi t1,t1,8; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    # ∂θ/∂px = -py/rxy2 → H[row+1,col+0]
    fneg.d ft0,fs1
    fdiv.d ft0,ft0,ft5
    addi t2,t3,1
    li t0,276; mul t1,t2,t0; add t1,t1,t4; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    # ∂θ/∂py = px/rxy2 → H[row+1,col+4]
    fdiv.d ft0,fs0,ft5
    li t0,276; mul t1,t2,t0; add t1,t1,t4; addi t1,t1,4; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    # ∂φ/∂px = -px*pz/(r2*rho_xy) → H[row+2,col+0]
    fmul.d ft0,fs0,fs2
    fmul.d ft1,ft4,fs3
    fdiv.d ft0,ft0,ft1
    fneg.d ft0,ft0
    addi t2,t3,2
    li t0,276; mul t1,t2,t0; add t1,t1,t4; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    # ∂φ/∂py = -py*pz/(r2*rho_xy) → H[row+2,col+4]
    fmul.d ft0,fs1,fs2
    fmul.d ft1,ft4,fs3
    fdiv.d ft0,ft0,ft1
    fneg.d ft0,ft0
    li t0,276; mul t1,t2,t0; add t1,t1,t4; addi t1,t1,4; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    # ∂φ/∂pz = rho_xy/r2 → H[row+2,col+8]
    fdiv.d ft0,fs3,ft4
    li t0,276; mul t1,t2,t0; add t1,t1,t4; addi t1,t1,8; slli t1,t1,3
    add t1,s0,t1; fsd ft0,0(t1)
    
    addi s2,s2,1
    j .Lcj_loop
    
.Lcj_done:
    ld ra,88(sp); ld s0,80(sp); ld s1,72(sp); ld s2,64(sp)
    fld fs0,56(sp); fld fs1,48(sp); fld fs2,40(sp)
    fld fs3,32(sp); fld fs4,24(sp); fld fs5,16(sp); fld fs6,8(sp)
    addi sp,sp,96
    ret

# VECTORIZED HELPER FUNCTIONS (from working lkf_vector.s)

.globl vekf_zero
.align 2
vekf_zero:
    mv t0,a0; mv t1,a1
.Lvez_loop:
    beqz t1,.Lvez_done
    vsetvli t2,t1,e64,m8,ta,ma
    vmv.v.i v0,0
    vse64.v v0,(t0)
    slli t3,t2,3; add t0,t0,t3
    sub t1,t1,t2; j .Lvez_loop
.Lvez_done:
    ret

.globl vekf_copy
.align 2
vekf_copy:
    mv t0,a0; mv t1,a1; mv t2,a2
.Lvec_loop:
    beqz t2,.Lvec_done
    vsetvli t3,t2,e64,m8,ta,ma
    vle64.v v0,(t1); vse64.v v0,(t0)
    slli t4,t3,3; add t0,t0,t4; add t1,t1,t4
    sub t2,t2,t3; j .Lvec_loop
.Lvec_done:
    ret

.globl vekf_add
.align 2
vekf_add:
    mv t0,a0; mv t1,a1; mv t2,a2; mv t3,a3
.Lvea_loop:
    beqz t3,.Lvea_done
    vsetvli t4,t3,e64,m4,ta,ma
    vle64.v v0,(t1); vle64.v v4,(t2)
    vfadd.vv v8,v0,v4; vse64.v v8,(t0)
    slli t5,t4,3
    add t0,t0,t5; add t1,t1,t5; add t2,t2,t5
    sub t3,t3,t4; j .Lvea_loop
.Lvea_done:
    ret

.globl vekf_sub
.align 2
vekf_sub:
    mv t0,a0; mv t1,a1; mv t2,a2; mv t3,a3
.Lves_loop:
    beqz t3,.Lves_done
    vsetvli t4,t3,e64,m4,ta,ma
    vle64.v v0,(t1); vle64.v v4,(t2)
    vfsub.vv v8,v0,v4; vse64.v v8,(t0)
    slli t5,t4,3
    add t0,t0,t5; add t1,t1,t5; add t2,t2,t5
    sub t3,t3,t4; j .Lves_loop
.Lves_done:
    ret

.globl vekf_add_inplace
.align 2
vekf_add_inplace:
    mv t0,a0; mv t1,a1; mv t2,a2
.Lveai_loop:
    beqz t2,.Lveai_done
    vsetvli t3,t2,e64,m8,ta,ma
    vle64.v v0,(t0); vle64.v v8,(t1)
    vfadd.vv v16,v0,v8; vse64.v v16,(t0)
    slli t4,t3,3; add t0,t0,t4; add t1,t1,t4
    sub t2,t2,t3; j .Lveai_loop
.Lveai_done:
    ret

# vlkf_matmul(C, A, B, M, K, N)
# C = A * B
# a0=C a1=A a2=B a3=M a4=K a5=N
.globl vlkf_matmul
.align 2
vekf_matmul:
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

.globl vekf_transpose
.align 2
vekf_transpose:
    addi sp,sp,-48
    sd ra,40(sp); sd s0,32(sp); sd s1,24(sp)
    sd s2,16(sp); sd s3,8(sp);  sd s4,0(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    slli s4,s2,3

    li t5,0
.Lvet_i:
    bge t5,s2,.Lvet_done
    mul t0,t5,s3; slli t0,t0,3; add t3,s1,t0
    slli t0,t5,3;               add t4,s0,t0
    mv t0,t3; mv t1,t4; mv t2,s3
.Lvet_j:
    beqz t2,.Lvet_ni
    vsetvli a0,t2,e64,m4,ta,ma
    vle64.v v0,(t0)
    vsse64.v v0,(t1),s4
    slli a1,a0,3; add t0,t0,a1
    mul  a2,a0,s4; add t1,t1,a2
    sub  t2,t2,a0; j .Lvet_j
.Lvet_ni:
    addi t5,t5,1; j .Lvet_i
.Lvet_done:
    ld ra,40(sp); ld s0,32(sp); ld s1,24(sp)
    ld s2,16(sp); ld s3,8(sp);  ld s4,0(sp)
    addi sp,sp,48
    ret

.globl vekf_identity_minus
.align 2
vekf_identity_minus:
    addi sp,sp,-32
    sd ra,24(sp); sd s0,16(sp); sd s1,8(sp); sd s2,0(sp)
    mv s0,a0; mv s1,a1; mv s2,a2

    mul t0,s2,s2; mv t1,s0; mv t2,s1; mv t3,t0
.Lveim_neg:
    beqz t3,.Lveim_diag
    vsetvli t4,t3,e64,m8,ta,ma
    vle64.v v0,(t2); vfneg.v v8,v0; vse64.v v8,(t1)
    slli t5,t4,3; add t1,t1,t5; add t2,t2,t5
    sub t3,t3,t4; j .Lveim_neg

.Lveim_diag:
    li t0,0
.Lveim_d_loop:
    bge t0,s2,.Lveim_d_done
    mul t1,t0,s2; add t1,t1,t0; slli t1,t1,3; add t1,s0,t1
    fld ft0,0(t1); li t2,1; fcvt.d.w ft1,t2
    fadd.d ft0,ft0,ft1; fsd ft0,0(t1)
    addi t0,t0,1; j .Lveim_d_loop
.Lveim_d_done:
    ld ra,24(sp); ld s0,16(sp); ld s1,8(sp); ld s2,0(sp)
    addi sp,sp,32
    ret

.globl init_identity
.align 2
init_identity:
    addi sp,sp,-32
    sd ra,24(sp); sd s0,16(sp); sd s1,8(sp)
    mv s0,a0; mv s1,a1

    mul a1,s1,s1; mv a0,s0; call vekf_zero

    li t0,0
.Linit_diag:
    bge t0,s1,.Linit_done
    mul t1,t0,s1; add t1,t1,t0; slli t1,t1,3; add t1,s0,t1
    li t2,1; fcvt.d.w ft0,t2; fsd ft0,0(t1)
    addi t0,t0,1; j .Linit_diag
.Linit_done:
    ld ra,24(sp); ld s0,16(sp); ld s1,8(sp)
    addi sp,sp,32
    ret

.globl chol_solve_matrix
.align 2
chol_solve_matrix:
    addi sp,sp,-80
    sd ra,72(sp); sd s0,64(sp); sd s1,56(sp); sd s2,48(sp)
    sd s3,40(sp); sd s4,32(sp); sd s5,24(sp); sd s6,16(sp); sd s7,8(sp)
    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3; mv s4,a4

    slli a0,s3,3; call malloc; mv s5,a0
    slli a0,s3,3; call malloc; mv s6,a0
    slli s7,s4,3

    li t6,0
.Lcsm_col:
    bge t6,s4,.Lcsm_done

    slli t0,t6,3; add t0,s2,t0
    mv a0,s5; mv a1,t0; mv a2,s7; mv a3,s3
    call vekf_strided_load

    mv a0,s6; mv a1,s1; mv a2,s5; mv a3,s3
    call forward_substitution

    mv a0,s5; mv a1,s1; mv a2,s6; mv a3,s3
    call backward_substitution

    slli t0,t6,3; add t0,s0,t0
    mv a0,t0; mv a1,s5; mv a2,s7; mv a3,s3
    call vekf_strided_store

    addi t6,t6,1; j .Lcsm_col
.Lcsm_done:
    mv a0,s5; call free
    mv a0,s6; call free
    ld ra,72(sp); ld s0,64(sp); ld s1,56(sp); ld s2,48(sp)
    ld s3,40(sp); ld s4,32(sp); ld s5,24(sp); ld s6,16(sp); ld s7,8(sp)
    addi sp,sp,80
    ret

.globl vekf_strided_load
.align 2
vekf_strided_load:
    mv t0,a0; mv t1,a1; mv t2,a3
.Lvsl_loop:
    beqz t2,.Lvsl_done
    vsetvli t3,t2,e64,m4,ta,ma
    vlse64.v v0,(t1),a2
    vse64.v  v0,(t0)
    slli t4,t3,3; add t0,t0,t4
    mul  t4,t3,a2; add t1,t1,t4
    sub  t2,t2,t3; j .Lvsl_loop
.Lvsl_done:
    ret

.globl vekf_strided_store
.align 2
vekf_strided_store:
    mv t0,a0; mv t1,a1; mv t2,a3
.Lvss_loop:
    beqz t2,.Lvss_done
    vsetvli t3,t2,e64,m4,ta,ma
    vle64.v  v0,(t1)
    vsse64.v v0,(t0),a2
    slli t4,t3,3; add t1,t1,t4
    mul  t4,t3,a2; add t0,t0,t4
    sub  t2,t2,t3; j .Lvss_loop
.Lvss_done:
    ret

# SCALAR FUNCTIONS (unchanged from M3)


.globl cholesky_decomposition
.align 2
cholesky_decomposition:
    addi sp,sp,-96
    sd  ra,88(sp); sd s0,80(sp); sd s1,72(sp); sd s2,64(sp)
    sd  s3,56(sp); sd s4,48(sp); sd s5,40(sp)
    fsd fs0,32(sp); fsd fs1,24(sp)

    mv s0,a0; mv s1,a1; mv s2,a2; li s3,0

.Lchol_i:
    bge s3,s2,.Lchol_succ
    li s4,0
.Lchol_j:
    bgt s4,s3,.Lchol_zero_upper

    mul t0,s3,s2; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0
    fld fs0,0(t0)

    li s5,0
.Lchol_k:
    bge s5,s4,.Lchol_k_done
    mul t0,s3,s2; add t0,t0,s5; slli t0,t0,3; add t0,s0,t0
    fld ft0,0(t0)
    mul t1,s4,s2; add t1,t1,s5; slli t1,t1,3; add t1,s0,t1
    fld ft1,0(t1)
    fnmsub.d fs0,ft0,ft1,fs0
    addi s5,s5,1; j .Lchol_k
.Lchol_k_done:

    bne s3,s4,.Lchol_offdiag

.Lchol_diag:
    fcvt.d.w ft0,zero
    fle.d t0,fs0,ft0
    bnez t0,.Lchol_fail
    fsqrt.d fs0,fs0
    mul t0,s3,s2; add t0,t0,s4; slli t0,t0,3; add t0,s0,t0
    fsd fs0,0(t0)
    j .Lchol_next_j

.Lchol_offdiag:
    mul t0,s4,s2; add t0,t0,s4; slli t0,t0,3; add t0,s0,t0
    fld fs1,0(t0)
    fdiv.d fs0,fs0,fs1
    mul t0,s3,s2; add t0,t0,s4; slli t0,t0,3; add t0,s0,t0
    fsd fs0,0(t0)

.Lchol_next_j:
    addi s4,s4,1; j .Lchol_j

.Lchol_zero_upper:
    addi s4,s3,1
.Lchol_zero_loop:
    bge s4,s2,.Lchol_next_i
    mul t0,s3,s2; add t0,t0,s4; slli t0,t0,3; add t0,s0,t0
    fcvt.d.w ft0,zero; fsd ft0,0(t0)
    addi s4,s4,1; j .Lchol_zero_loop

.Lchol_next_i:
    addi s3,s3,1; j .Lchol_i

.Lchol_succ:
    li a0,0; j .Lchol_end
.Lchol_fail:
    li a0,-1
.Lchol_end:
    ld  ra,88(sp); ld s0,80(sp); ld s1,72(sp); ld s2,64(sp)
    ld  s3,56(sp); ld s4,48(sp); ld s5,40(sp)
    fld fs0,32(sp); fld fs1,24(sp)
    addi sp,sp,96
    ret

.globl forward_substitution
.align 2
forward_substitution:
    addi sp,sp,-80
    sd ra,72(sp); sd s0,64(sp); sd s1,56(sp); sd s2,48(sp)
    sd s3,40(sp); sd s4,32(sp); sd s5,24(sp)
    fsd fs0,16(sp); fsd fs1,8(sp)

    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3; li s4,0

.Lfs_i:
    bge s4,s3,.Lfs_done

    slli t0,s4,3; add t0,s2,t0; fld fs0,0(t0)

    li s5,0
.Lfs_j:
    bge s5,s4,.Lfs_j_done
    mul t0,s4,s3; add t0,t0,s5; slli t0,t0,3; add t0,s1,t0
    fld ft0,0(t0)
    slli t1,s5,3; add t1,s0,t1; fld ft1,0(t1)
    fnmsub.d fs0,ft0,ft1,fs0
    addi s5,s5,1; j .Lfs_j
.Lfs_j_done:

    mul t0,s4,s3; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0
    fld fs1,0(t0)
    fdiv.d fs0,fs0,fs1

    slli t0,s4,3; add t0,s0,t0; fsd fs0,0(t0)

    addi s4,s4,1; j .Lfs_i
.Lfs_done:
    ld ra,72(sp); ld s0,64(sp); ld s1,56(sp); ld s2,48(sp)
    ld s3,40(sp); ld s4,32(sp); ld s5,24(sp)
    fld fs0,16(sp); fld fs1,8(sp)
    addi sp,sp,80
    ret

.globl backward_substitution
.align 2
backward_substitution:
    addi sp,sp,-80
    sd ra,72(sp); sd s0,64(sp); sd s1,56(sp); sd s2,48(sp)
    sd s3,40(sp); sd s4,32(sp); sd s5,24(sp)
    fsd fs0,16(sp); fsd fs1,8(sp)

    mv s0,a0; mv s1,a1; mv s2,a2; mv s3,a3
    mv s4,s3; addi s4,s4,-1

.Lbs_i:
    bltz s4,.Lbs_done

    slli t0,s4,3; add t0,s2,t0; fld fs0,0(t0)

    mv s5,s4; addi s5,s5,1
.Lbs_j:
    bge s5,s3,.Lbs_j_done
    mul t0,s5,s3; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0
    fld ft0,0(t0)
    slli t1,s5,3; add t1,s0,t1; fld ft1,0(t1)
    fnmsub.d fs0,ft0,ft1,fs0
    addi s5,s5,1; j .Lbs_j
.Lbs_j_done:

    mul t0,s4,s3; add t0,t0,s4; slli t0,t0,3; add t0,s1,t0
    fld fs1,0(t0)
    fdiv.d fs0,fs0,fs1

    slli t0,s4,3; add t0,s0,t0; fsd fs0,0(t0)

    addi s4,s4,-1; j .Lbs_i
.Lbs_done:
    ld ra,72(sp); ld s0,64(sp); ld s1,56(sp); ld s2,48(sp)
    ld s3,40(sp); ld s4,32(sp); ld s5,24(sp)
    fld fs0,16(sp); fld fs1,8(sp)
    addi sp,sp,80
    ret

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
    mv a0,s0; li a1,NN; call vekf_zero

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

# buildQ — Construct 276×276 process-noise covariance (run once)

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
    mv a0,s0; li a1,NN; call vekf_zero

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
    mv a0,s0; li a1,NM; call vekf_zero

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
matrix_zero:         j vekf_zero

.globl matrix_copy
matrix_copy:         j vekf_copy

.globl matrix_add
matrix_add:          j vekf_add

.globl matrix_subtract
matrix_subtract:     j vekf_sub

.globl vector_add_inplace
vector_add_inplace:  j vekf_add_inplace

.globl matrix_transpose
matrix_transpose:    j vekf_transpose

.globl identity_minus
identity_minus:      j vekf_identity_minus

.globl matrix_multiply
matrix_multiply:     j vekf_matmul