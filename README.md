# Kalman Filter Implementation in RISC-V Vector Assembly

This repository contains the **Milestone 4** submission of our Kalman Filter project for **Computer Architecture and Assembly Language (Spring 2026)** at IBA Karachi.

## Team ChaiGPT
- Laiba Irfan (31736)
- Ameer Abdullah (30535)
- M Usman (29126)
- M Arrayyan Asad (30557)
---

# Overview

Milestone 4 extends our Milestone 3 scalar RISC-V assembly implementation by accelerating the core computational kernels using the **RISC-V Vector Extension (RVV)**.

The project demonstrates low-level performance optimisation through vectorisation of dense linear algebra routines used throughout the Kalman Filter pipeline. By replacing repeated scalar operations with RVV vector instructions, the implementation improves execution speed while preserving numerical correctness against the scalar reference.

Both the Linear Kalman Filter and Extended Kalman Filter were vectorised and verified against the Milestone 3 implementation.

---

# Key Features

### RVV Vectorised Computation
Major computational kernels were rewritten using the RISC-V Vector Extension for parallel execution across multiple double-precision elements.

### Scalable State Tracking
Tracks **23 joints** with **12 states per joint**:
- Position
- Velocity
- Acceleration
- Jerk

across the **x, y, and z** axes.

### Double-Precision Floating Point
Uses 64-bit floating-point arithmetic throughout to maintain numerical consistency with previous milestones.

### Vectorised Matrix Operations
Accelerates:
- Matrix multiplication
- Matrix-vector multiplication
- Matrix transpose
- Vector addition / subtraction
- Vector copy and zero initialisation
- Strided load/store operations

### EKF Nonlinear Measurement Model
The EKF preserves the nonlinear measurement update from previous milestones, including:
- custom `arctan2` approximation
- Jacobian computation
- angle wrapping

### Numerical Verification
Vectorised output is validated against the Milestone 3 scalar assembly reference to ensure numerical correctness.

---

# System Specifications

| Parameter | Value |
|---|---:|
| ISA Extensions | RV64I, RV64M, RV64F, RV64D, RVV |
| State Dimension (n) | 276 |
| Measurement Dimension (m) | 69 |
| Number of Joints | 23 |
| States per Joint | 12 |
| Sampling Rate | 100 Hz |
| Data Duration | 3,040 Frames |

---

# Technical Implementation

## Memory Management

The implementation follows the same memory strategy introduced in Milestone 3, extended for vectorised workloads.

### Stack
Used for function calls and register preservation.

### Heap
Persistent allocation for:
- State Vector (`x`)
- Covariance Matrix (`P`)

### Workspace
Pre-allocated workspace buffers are used for intermediate matrices required during prediction and update steps.

This avoids repeated dynamic allocation during the frame loop and improves runtime efficiency.

---

## RVV Optimisations

Key optimisations include:

### Vectorised Matrix Multiplication
Core matrix multiplication kernels are implemented using RVV instructions including:

- `vsetvli`
- `vle64.v`
- `vse64.v`
- `vfmacc.vf`

### Vector Transpose
Efficient transpose operations implemented using:

- `vlse64.v`
- `vsse64.v`

for strided memory access.

### Vector Arithmetic
Element-wise operations implemented using:

- `vfadd.vv`
- `vfsub.vv`

### Memory Operations
Vector copy / zero kernels implemented using:

- `vmv.v.i`
- `vle64.v`
- `vse64.v`

---

# Repository Files

| File | Description |
|---|---|
| `lkf_vector.s` | RVV implementation of Linear Kalman Filter |
| `ekf_vector.s` | RVV implementation of Extended Kalman Filter |
| `main_lkf_final.cpp` | Driver program for LKF |
| `main_ekf_final.cpp` | Driver program for EKF |
| `lkf_output_m4.csv` | Output of vectorised LKF |
| `ekf_output_m4.csv` | Output of vectorised EKF |
| `Makefile` | Build and execution targets |

---

# Build and Run

## Linear Kalman Filter

Build:

```bash
make lkf_vector
```

Run:

```bash
make run_lkf_vector
```

---

## Extended Kalman Filter

Build:

```bash
make ekf_vector
```

Run:

```bash
make run_ekf_vector
```


---

## Build Everything

```bash
make
```

---

# Verification Results

Verification condition:

```text
|x_vector - x_scalar| ≤ 1e-9
```

Outputs from both filters were checked across all frames to confirm that vectorisation preserved numerical correctness.

---

# Performance

Milestone 4 focuses on accelerating the computational bottlenecks of the Kalman Filter pipeline through RVV vectorisation.

Performance analysis includes:

- runtime comparison with Milestone 3 scalar assembly
- speedup analysis
- instruction count comparison
- kernel-level vector performance evaluation

The largest performance improvements come from vectorising repeated dense matrix operations used throughout both LKF and EKF.

---
