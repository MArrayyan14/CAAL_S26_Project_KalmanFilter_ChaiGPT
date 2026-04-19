# Makefile for Kalman Filter Milestone 3 
# Team ChaiGPT

# Configuration
PREFIX = riscv64-unknown-elf-
AS = $(PREFIX)as
CXX = $(PREFIX)g++

ARCH = rv64g
ABI = lp64d

ASFLAGS = -march=$(ARCH) -mabi=$(ABI)
CXXFLAGS = -march=$(ARCH) -mabi=$(ABI) -O2 -static -std=c++11
LDFLAGS = -static -lm

# QEMU (if using QEMU instead of Spike)
QEMU = qemu-riscv64

# Targets
.PHONY: all clean run_lkf run_ekf

all: lkf_final ekf_final

# Build LKF
lkf_final: lkf_final.o main_lkf_final.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✓ Built LKF: lkf_final"

# Build EKF
ekf_final: ekf_final.o main_ekf_final.o
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)
	@echo "✓ Built EKF: ekf_final"

# Assembly compilation
%.o: %.s
	$(AS) $(ASFLAGS) -o $@ $<

# C++ compilation
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Run targets
run_lkf: lkf_final
	@echo "Running LKF on gait_data_noisy.csv..."
	$(QEMU) ./lkf_final gait_data_noisy.csv lkf_output_m3.csv
	@echo "✓ Output saved to lkf_output_m3.csv"

run_ekf: ekf_final
	@echo "Running EKF on gait_data_noisy.csv..."
	$(QEMU) ./ekf_final gait_data_noisy.csv ekf_output_m3.csv
	@echo "✓ Output saved to ekf_output_m3.csv"

run_all: run_lkf run_ekf
	@echo "✓ Both filters complete!"

# Clean
clean:
	rm -f *.o lkf_final ekf_final
	@echo "✓ Cleaned"

clean_output:
	rm -f *_output_m3.csv
	@echo "✓ Cleaned output files"

clean_all: clean clean_output

# Help
help:
	@echo "Milestone 3 - Final Standalone Kalman Filters"
	@echo ""
	@echo "Targets:"
	@echo "  all         - Build both LKF and EKF"
	@echo "  lkf_final   - Build LKF only"
	@echo "  ekf_final   - Build EKF only"
	@echo "  run_lkf     - Run LKF on gait_data_noisy.csv"
	@echo "  run_ekf     - Run EKF on gait_data_noisy.csv"
	@echo "  run_all     - Run both filters"
	@echo "  clean       - Remove executables and objects"
	@echo "  clean_all   - Remove everything including outputs"
	@echo ""
	@echo "Usage examples:"
	@echo "  make all                    # Build both"
	@echo "  make run_lkf                # Run LKF"
	@echo "  make run_ekf                # Run EKF"
	@echo "  ./lkf_final input.csv out.csv   # Custom input"
