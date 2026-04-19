# Kalman Filter Implementation in RISC-V Scalar Assembly
This repository contains a high-performance implementation of **Linear (LKF)** and **Extended Kalman Filters (EKF)** written in **RISC-V Scalar Assembly**. Developed as part of the Computer Architecture and Assembly Language course at IBA Karachi, this project tracks full-body motion capture data consisting of 23 joints and 276 state variables.

## Overview
The project demonstrates the power of low-level optimization by translating complex state-estimation algorithms from C++ to hand-optimized assembly. By utilizing the **RV64D** (Double-Precision Floating-Point) extension, we achieve bit-exact accuracy compared to high-level implementations.

### Key Features
- **Scalable State Tracking:** Tracks 23 joints with 12 states each (Position, Velocity, Accel, Jerk across 3 axes).
- **Mathematical Precision:** Implements Double-Precision IEEE 754 arithmetic.
- **Robust Inversion:** Uses **Cholesky Decomposition (LLT)** and triangular solvers instead of explicit matrix inversion for numerical stability.
- **EKF Nonlinearity:** Custom `arctan2` implementation using a 4th-degree minimax polynomial for fast coordinate transformation.
- **Joseph Form:** Covariance updates are performed using the Joseph form to ensure the matrix remains positive semi-definite.

---

## System Specifications

| Parameter | Value |
| :--- | :--- |
| **ISA Extensions** | RV64I, RV64M, RV64F, RV64D |
| **State Dimension ($n$)** | 276 |
| **Measurement Dimension ($m$)** | 69 (3 per joint: $r, \theta, \phi$) |
| **Sampling Rate** | 100 Hz ($\Delta t = 0.01s$) |
| **Data Duration** | 3,040 Frames (30.4 seconds) |
| **Verification** | **Bit-Exact Matching** ($0.000$ error vs C++) |

---

## Technical Implementation

### Memory Management
The system utilizes a three-tier memory strategy to eliminate `malloc` overhead during the main loop:
1. **Stack:** 16-byte aligned frames for function calls.
2. **Heap:** Persistent storage for the State Vector ($2.2$ KB) and Covariance Matrix ($609$ KB).
3. **Workspace:** A pre-allocated **5.7 MB** buffer used as a scratchpad for 19 intermediate matrices, avoiding dynamic allocation fragmentation.

### Assembly Optimizations
- **Instruction Level Parallelism:** Extensive use of `fmadd.d` (Fused Multiply-Add) to reduce the floating-point operation count by 33% in matrix multiplications.
- **Register Allocation:** Optimized use of `s0-s11` and `fs0-fs11` registers to minimize memory access (load/store) within inner loops.
- **Loop Unrolling:** Applied in critical sections of the matrix-vector operations.

---

## 📈 Verification Results

The implementation was tested against the Milestone 2 C++ reference using over **839,040 state estimates**.

| Filter Type | Max Error | Status |
| :--- | :--- | :--- |
| **Linear Kalman Filter (LKF)** | $0.000 \times 10^0$ | **PASS** ✅ |
| **Extended Kalman Filter (EKF)** | $0.000 \times 10^0$ | **PASS** ✅ |

### Performance Visualization
The project includes visualization scripts that generate:
- **3D Trajectory Plots:** Overlaying True vs. Noisy vs. Estimated paths.
- **Error Heatmaps:** Visualizing accuracy across all 23 joints.
- **Derivative Analysis:** Smooth estimation of Velocity, Acceleration, and Jerk.

## Prerequisites
GCC Toolchain: riscv64-unknown-elf-gcc
Emulator: qemu-riscv64 (User mode)
Build and Run
code
Bash
### Clone the repository
git clone https://github.com/MArrayyan14/CAAL_S26_Project_KalmanFilter_ChaiGPT.git
cd CAAL_S26_Project_KalmanFilter_ChaiGPT/milestone-3

### Assemble and Link
riscv64-unknown-elf-gcc -march=rv64g -mabi=lp64d src/*.s -o kalman_filter

### Execute
qemu-riscv64 ./kalman_filter

---

## Team ChaiGPT
Laiba Irfan (31736)
Ameer Abdullah (30535)
Usman Iftikhar (29126)
Arrayyan Iqbal (30557)

Institution: Institute of Business Administration (IBA), Karachi
Course: Computer Architecture and Assembly Language (Spring 2026)
