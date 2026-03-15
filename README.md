# CAAL S26 Project KalmanFilter Milestone 2

**Team ChaiGPT** | Laiba (31736) · Ameer (30535) · Usman (29126) · Arrayyan (30557)

Implementation of a Linear Kalman Filter (LKF) and Extended Kalman Filter (EKF) in C++ for 3D full-body human gait estimation across 23 joints.

---

## Files

| File | Description |
|------|-------------|
| `lkf.cpp` | Linear Kalman Filter — Cartesian measurement model |
| `ekf.cpp` | Extended Kalman Filter — spherical measurement model |
| `CMakeLists.txt` | CMake build configuration |
| `animate_skeleton.py` | Generates 4-panel 3D walking animation (True / Noisy / LKF / EKF) |
| `plot_results.py` | Generates all 8 analysis plots for a chosen joint |
| `gait_data_true.csv` | Ground truth positions (3040 × 69, with header) |
| `gait_data_noisy.csv` | Noisy measurements (3040 × 69, with header) |
| `lkf_output.csv` | LKF filtered state vectors (3040 × 276, no header) |
| `ekf_output.csv` | EKF filtered state vectors (3040 × 276, no header) |

---

## Build

Requires **Eigen3** and **CMake ≥ 3.10**.

```bash
mkdir build && cd build
cmake ..
make
```

This produces two executables: `lkf` and `ekf`.

---

## Run

```bash
# Linear Kalman Filter
./lkf ../gait_data_noisy.csv ../lkf_output.csv

# Extended Kalman Filter
./ekf ../gait_data_noisy.csv ../ekf_output.csv
```

---

## Visualisation

```bash
# Generate analysis plots (Joint 0)
python3 plot_results.py

# Generate 3D walking animation
python3 animate_skeleton.py \
    --true gait_data_true.csv \
    --noisy gait_data_noisy.csv \
    --lkf lkf_output.csv \
    --ekf ekf_output.csv \
    --output animation.mp4
```

Requires: `numpy`, `matplotlib`. For MP4 output, `ffmpeg` must be installed.

---

## State Vector Layout

Each joint has 12 states. All 23 joints are stacked into a single vector of size 276:

```
[px, vx, ax, jx,  py, vy, ay, jy,  pz, vz, az, jz]  ×  23 joints
```

In the output CSVs, positions for joint `j` are at columns `j*12`, `j*12+4`, `j*12+8`.

---

## Results (Joint 0)

| | Noise RMSE | LKF RMSE | EKF RMSE |
|--|--|--|--|
| X | 542.8 mm | 165.5 mm | 164.1 mm |
| Y | 310.4 mm | 99.7 mm | 101.1 mm |
| Z | 45.2 mm | 19.8 mm | 21.0 mm |
| **Average** | **299.5 mm** | **95.0 mm** | **95.4 mm** |

Both filters achieve ~68% noise reduction. See the full report PDF and animation in the [Drive link](https://drive.google.com/drive/folders/1m8FhrUTW-Cw1UNOjSTIPtesxCXueZNSV?usp=share_link).
