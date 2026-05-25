# ============================================================================
# Makefile for Kalman Filter Milestone 4 - 
# Team ChaiGPT
# ============================================================================

# ============================================================================
# Configuration
# ============================================================================
PREFIX = riscv64-linux-gnu-

AS  = $(PREFIX)as
CXX = $(PREFIX)g++

ARCH = rv64gv
ABI  = lp64d

ASFLAGS  = -march=$(ARCH) -mabi=$(ABI)
CXXFLAGS = -march=$(ARCH) -mabi=$(ABI) -O2 -std=c++11

LDFLAGS = -lm

# NEW
QEMU = qemu-riscv64 -cpu max,v=true,vlen=1024,elen=64 -tb-size 2048 -L /usr/riscv64-linux-gnu

# ============================================================================
# Targets
# ============================================================================
.PHONY: all clean clean_all help run_compare_lkf run_compare_ekf

all: lkf_scalar lkf_vector ekf_scalar ekf_vector

# ============================================================================
# LKF BUILD
# ============================================================================

lkf_scalar: lkf_final.o main_lkf_final.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✓ Built LKF Scalar"

lkf_vector: lkf_vector.o main_lkf_final.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✓ Built LKF Vector"

# ============================================================================
# EKF BUILD
# ============================================================================

ekf_scalar: ekf_final.o main_ekf_final.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✓ Built EKF Scalar"

ekf_vector: ekf_vector.o main_ekf_final.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✓ Built EKF Vector"

# ============================================================================
# Assembly compilation
# ============================================================================
%.o: %.s
	$(AS) $(ASFLAGS) -o $@ $<

# ============================================================================
# C++ compilation
# ============================================================================
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# ============================================================================
# RUN FULL DATASET
# ============================================================================

run_lkf_scalar:
	$(QEMU) ./lkf_scalar gait_data_noisy.csv lkf_scalar_out.csv

run_lkf_vector:
	$(QEMU) ./lkf_vector gait_data_noisy.csv lkf_vector_out.csv

run_ekf_scalar:
	$(QEMU) ./ekf_scalar gait_data_noisy.csv ekf_scalar_out.csv

run_ekf_vector:
	$(QEMU) ./ekf_vector gait_data_noisy.csv ekf_vector_out.csv

# ============================================================================
# BENCHMARK (FULL DATASET COMPARISON)
# ============================================================================

run_compare_lkf:
	@echo "===== LKF BENCHMARK ====="
	@echo "SCALAR:"
	/usr/bin/time -p $(QEMU) ./lkf_scalar gait_data_noisy.csv out1.csv
	@echo "VECTOR:"
	/usr/bin/time -p $(QEMU) ./lkf_vector gait_data_noisy.csv out2.csv

run_compare_ekf:
	@echo "===== EKF BENCHMARK ====="
	@echo "SCALAR:"
	/usr/bin/time -p $(QEMU) ./ekf_scalar gait_data_noisy.csv out1.csv
	@echo "VECTOR:"
	/usr/bin/time -p $(QEMU) ./ekf_vector gait_data_noisy.csv out2.csv

# ============================================================================
# CLEAN
# ============================================================================

clean:
	rm -f *.o lkf_scalar lkf_vector ekf_scalar ekf_vector temp_20.csv
	@echo "✓ Cleaned build files"

clean_all: clean
	rm -f *_out.csv
	@echo "✓ Cleaned outputs"

# ============================================================================
# HELP
# ============================================================================

help:
	@echo "Milestone 4 - Clean Benchmark System"
	@echo ""
	@echo "Build:"
	@echo "  make all"
	@echo ""
	@echo "Run:"
	@echo "  make run_lkf_scalar"
	@echo "  make run_lkf_vector"
	@echo "  make run_compare_lkf"
	@echo ""
	@echo "  make run_ekf_scalar"
	@echo "  make run_ekf_vector"
	@echo "  make run_compare_ekf"
	@echo ""