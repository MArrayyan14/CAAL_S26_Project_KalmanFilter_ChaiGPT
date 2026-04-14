#!/usr/bin/env python3
"""
Comparison Plots: Milestone 2 (C++) vs Milestone 3 (Assembly)
Team ChaiGPT - RISC-V Kalman Filter

Generates comprehensive comparison plots showing:
1. M2 vs M3 state trajectories
2. Error plots (should be zero everywhere)
3. Joint-wise comparisons
4. Statistical analysis plots
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import os

# Set style for publication-quality plots
plt.style.use('seaborn-v0_8-darkgrid')
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['axes.titlesize'] = 14
plt.rcParams['legend.fontsize'] = 10

def load_data(m2_file, m3_file):
    """Load both M2 and M3 outputs"""
    print(f"Loading {m2_file}...")
    m2_data = np.loadtxt(m2_file, delimiter=',')
    print(f"  Shape: {m2_data.shape}")
    
    print(f"Loading {m3_file}...")
    m3_data = np.loadtxt(m3_file, delimiter=',')
    print(f"  Shape: {m3_data.shape}")
    
    if m2_data.shape != m3_data.shape:
        print("ERROR: Shape mismatch!")
        sys.exit(1)
    
    return m2_data, m3_data

def plot_state_comparison(m2_data, m3_data, joint_idx, state_name, output_dir):
    """
    Plot comparison of specific state component across M2 and M3
    """
    num_frames = m2_data.shape[0]
    time = np.arange(num_frames) * 0.01  # dt = 0.01s
    
    # State indices mapping
    state_map = {
        'px': 0, 'vx': 1, 'ax': 2, 'jx': 3,
        'py': 4, 'vy': 5, 'ay': 6, 'jy': 7,
        'pz': 8, 'vz': 9, 'az': 10, 'jz': 11
    }
    
    state_idx = joint_idx * 12 + state_map[state_name]
    
    m2_state = m2_data[:, state_idx]
    m3_state = m3_data[:, state_idx]
    error = m3_state - m2_state
    
    fig, axes = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: Overlay
    axes[0].plot(time, m2_state, 'b-', label='M2 (C++)', linewidth=1.5, alpha=0.7)
    axes[0].plot(time, m3_state, 'r--', label='M3 (Assembly)', linewidth=1.5, alpha=0.7)
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel(f'{state_name}')
    axes[0].set_title(f'Joint {joint_idx} - {state_name}: M2 vs M3 Comparison')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Plot 2: Error
    axes[1].plot(time, error, 'g-', linewidth=1)
    axes[1].axhline(y=0, color='k', linestyle='--', alpha=0.3)
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Error (M3 - M2)')
    axes[1].set_title(f'Error: Joint {joint_idx} - {state_name}')
    axes[1].grid(True, alpha=0.3)
    
    # Add error statistics
    max_err = np.max(np.abs(error))
    mean_err = np.mean(np.abs(error))
    axes[1].text(0.02, 0.98, f'Max Error: {max_err:.2e}\nMean Error: {mean_err:.2e}',
                transform=axes[1].transAxes, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    filename = f'{output_dir}/comparison_joint{joint_idx}_{state_name}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

def plot_all_joints_single_state(m2_data, m3_data, state_name, output_dir):
    """
    Plot specific state component for all 23 joints
    """
    state_map = {
        'px': 0, 'vx': 1, 'ax': 2, 'jx': 3,
        'py': 4, 'vy': 5, 'ay': 6, 'jy': 7,
        'pz': 8, 'vz': 9, 'az': 10, 'jz': 11
    }
    
    num_frames = m2_data.shape[0]
    time = np.arange(num_frames) * 0.01
    
    fig, axes = plt.subplots(5, 5, figsize=(20, 16))
    axes = axes.flatten()
    
    for joint_idx in range(23):
        state_idx = joint_idx * 12 + state_map[state_name]
        
        m2_state = m2_data[:, state_idx]
        m3_state = m3_data[:, state_idx]
        
        axes[joint_idx].plot(time, m2_state, 'b-', label='M2', linewidth=1, alpha=0.7)
        axes[joint_idx].plot(time, m3_state, 'r--', label='M3', linewidth=1, alpha=0.7)
        axes[joint_idx].set_title(f'Joint {joint_idx}', fontsize=10)
        axes[joint_idx].grid(True, alpha=0.3)
        
        if joint_idx == 0:
            axes[joint_idx].legend(fontsize=8)
    
    # Hide unused subplots
    for idx in range(23, 25):
        axes[idx].axis('off')
    
    fig.suptitle(f'All Joints - {state_name}: M2 vs M3', fontsize=16)
    plt.tight_layout()
    
    filename = f'{output_dir}/all_joints_{state_name}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

def plot_error_heatmap(m2_data, m3_data, output_dir):
    """
    Create heatmap of errors across all joints and states
    """
    errors = np.abs(m3_data - m2_data)
    
    # Reshape to (frames, joints, states_per_joint)
    num_frames = errors.shape[0]
    errors_reshaped = errors.reshape(num_frames, 23, 12)
    
    # Average over frames
    avg_errors = np.mean(errors_reshaped, axis=0)  # (23, 12)
    
    fig, ax = plt.subplots(figsize=(14, 10))
    
    im = ax.imshow(avg_errors, cmap='hot', aspect='auto', interpolation='nearest')
    
    ax.set_xlabel('State Component', fontsize=12)
    ax.set_ylabel('Joint', fontsize=12)
    ax.set_title('Average Absolute Error: M3 vs M2\n(Averaged over 3040 frames)', fontsize=14)
    
    # Set ticks
    ax.set_xticks(np.arange(12))
    ax.set_xticklabels(['px', 'vx', 'ax', 'jx', 'py', 'vy', 'ay', 'jy', 'pz', 'vz', 'az', 'jz'])
    ax.set_yticks(np.arange(23))
    ax.set_yticklabels([f'J{i}' for i in range(23)])
    
    # Colorbar
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Average Absolute Error', rotation=270, labelpad=20)
    
    # Add text annotations
    max_error = np.max(avg_errors)
    ax.text(0.02, 0.98, f'Max Error: {max_error:.2e}',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    filename = f'{output_dir}/error_heatmap.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

def plot_error_distribution(m2_data, m3_data, output_dir):
    """
    Plot distribution of errors
    """
    errors = np.abs(m3_data - m2_data).flatten()
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # Histogram
    axes[0].hist(errors, bins=100, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Absolute Error')
    axes[0].set_ylabel('Frequency')
    axes[0].set_title('Error Distribution')
    axes[0].set_yscale('log')
    axes[0].grid(True, alpha=0.3)
    
    # Box plot by state component
    errors_by_state = []
    state_names = ['px', 'vx', 'ax', 'jx', 'py', 'vy', 'ay', 'jy', 'pz', 'vz', 'az', 'jz']
    
    for state_idx in range(12):
        indices = [joint * 12 + state_idx for joint in range(23)]
        state_errors = np.abs(m3_data[:, indices] - m2_data[:, indices]).flatten()
        errors_by_state.append(state_errors)
    
    bp = axes[1].boxplot(errors_by_state, labels=state_names, patch_artist=True)
    axes[1].set_xlabel('State Component')
    axes[1].set_ylabel('Absolute Error')
    axes[1].set_title('Error Distribution by State Component')
    axes[1].grid(True, alpha=0.3, axis='y')
    axes[1].set_yscale('log')
    
    # Color boxes
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    
    plt.tight_layout()
    filename = f'{output_dir}/error_distribution.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

def plot_trajectory_3d(m2_data, m3_data, joint_idx, output_dir):
    """
    3D trajectory plot for a specific joint
    """
    from mpl_toolkits.mplot3d import Axes3D
    
    # Extract position states
    px_m2 = m2_data[:, joint_idx * 12 + 0]
    py_m2 = m2_data[:, joint_idx * 12 + 4]
    pz_m2 = m2_data[:, joint_idx * 12 + 8]
    
    px_m3 = m3_data[:, joint_idx * 12 + 0]
    py_m3 = m3_data[:, joint_idx * 12 + 4]
    pz_m3 = m3_data[:, joint_idx * 12 + 8]
    
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Plot M2 trajectory
    ax.plot(px_m2, py_m2, pz_m2, 'b-', label='M2 (C++)', linewidth=1.5, alpha=0.6)
    
    # Plot M3 trajectory
    ax.plot(px_m3, py_m3, pz_m3, 'r--', label='M3 (Assembly)', linewidth=1.5, alpha=0.6)
    
    # Mark start and end points
    ax.scatter([px_m2[0]], [py_m2[0]], [pz_m2[0]], c='g', marker='o', s=100, label='Start')
    ax.scatter([px_m2[-1]], [py_m2[-1]], [pz_m2[-1]], c='k', marker='x', s=100, label='End')
    
    ax.set_xlabel('X Position (m)')
    ax.set_ylabel('Y Position (m)')
    ax.set_zlabel('Z Position (m)')
    ax.set_title(f'3D Trajectory - Joint {joint_idx}: M2 vs M3')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    filename = f'{output_dir}/trajectory_3d_joint{joint_idx}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

def plot_summary_dashboard(m2_data, m3_data, filter_name, output_dir):
    """
    Create comprehensive summary dashboard
    """
    errors = np.abs(m3_data - m2_data)
    
    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.3, wspace=0.3)
    
    # Plot 1: Overall error over time
    ax1 = fig.add_subplot(gs[0, :])
    max_errors_per_frame = np.max(errors, axis=1)
    mean_errors_per_frame = np.mean(errors, axis=1)
    time = np.arange(len(max_errors_per_frame)) * 0.01
    
    ax1.plot(time, max_errors_per_frame, 'r-', label='Max Error per Frame', linewidth=1)
    ax1.plot(time, mean_errors_per_frame, 'b-', label='Mean Error per Frame', linewidth=1)
    ax1.axhline(y=1e-9, color='g', linestyle='--', label='Threshold (1e-9)')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Error')
    ax1.set_title(f'{filter_name}: Error Evolution Over Time')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_yscale('log')
    
    # Plot 2: Error per joint
    ax2 = fig.add_subplot(gs[1, 0])
    errors_reshaped = errors.reshape(errors.shape[0], 23, 12)
    avg_errors_per_joint = np.mean(errors_reshaped, axis=(0, 2))
    
    ax2.bar(range(23), avg_errors_per_joint, color='steelblue', edgecolor='black')
    ax2.set_xlabel('Joint')
    ax2.set_ylabel('Average Error')
    ax2.set_title('Error by Joint')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Error per state component
    ax3 = fig.add_subplot(gs[1, 1])
    state_names = ['px', 'vx', 'ax', 'jx', 'py', 'vy', 'ay', 'jy', 'pz', 'vz', 'az', 'jz']
    avg_errors_per_state = []
    for state_idx in range(12):
        indices = [joint * 12 + state_idx for joint in range(23)]
        avg_errors_per_state.append(np.mean(errors[:, indices]))
    
    ax3.bar(range(12), avg_errors_per_state, color='coral', edgecolor='black')
    ax3.set_xticks(range(12))
    ax3.set_xticklabels(state_names, rotation=45)
    ax3.set_ylabel('Average Error')
    ax3.set_title('Error by State Component')
    ax3.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Error histogram
    ax4 = fig.add_subplot(gs[1, 2])
    ax4.hist(errors.flatten(), bins=100, edgecolor='black', alpha=0.7, color='green')
    ax4.set_xlabel('Absolute Error')
    ax4.set_ylabel('Frequency')
    ax4.set_title('Error Distribution')
    ax4.set_yscale('log')
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: Sample trajectories
    ax5 = fig.add_subplot(gs[2, :])
    joint_idx = 0  # Right ankle
    for state, offset, name in [(0, 0, 'px'), (4, 0.1, 'py'), (8, 0.2, 'pz')]:
        state_idx = joint_idx * 12 + state
        ax5.plot(time, m2_data[:, state_idx] + offset, 'b-', linewidth=1, alpha=0.6, label=f'M2 {name}')
        ax5.plot(time, m3_data[:, state_idx] + offset, 'r--', linewidth=1, alpha=0.6, label=f'M3 {name}')
    
    ax5.set_xlabel('Time (s)')
    ax5.set_ylabel('Position (m)')
    ax5.set_title(f'Sample Trajectories - Joint {joint_idx}')
    ax5.legend(ncol=3, fontsize=8)
    ax5.grid(True, alpha=0.3)
    
    # Add overall statistics text
    stats_text = f"""
    Total Frames: {m2_data.shape[0]}
    States per Frame: {m2_data.shape[1]}
    Max Error: {np.max(errors):.3e}
    Mean Error: {np.mean(errors):.3e}
    Threshold: 1.0e-09
    Status: {'PASS ✓' if np.max(errors) < 1e-9 or np.max(errors) == 0 else 'FAIL ✗'}
    """
    
    fig.text(0.02, 0.02, stats_text, fontsize=10, verticalalignment='bottom',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    fig.suptitle(f'{filter_name}: M2 vs M3 Comparison Dashboard', fontsize=16, fontweight='bold')
    
    filename = f'{output_dir}/dashboard_{filter_name.lower().replace(" ", "_")}.png'
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {filename}")

def main():
    if len(sys.argv) < 2:
        print("Usage:")
        print(f"  {sys.argv[0]} lkf    # Compare LKF outputs")
        print(f"  {sys.argv[0]} ekf    # Compare EKF outputs")
        print(f"  {sys.argv[0]} both   # Compare both")
        sys.exit(1)
    
    mode = sys.argv[1].lower()
    
    if mode in ['lkf', 'both']:
        print("\n" + "="*70)
        print("  GENERATING LKF COMPARISON PLOTS")
        print("="*70 + "\n")
        
        m2_lkf, m3_lkf = load_data('lkf_output_m2.csv', 'lkf_output_m3.csv')
        
        output_dir = 'plots_lkf_comparison'
        os.makedirs(output_dir, exist_ok=True)
        
        print("\n1. Generating summary dashboard...")
        plot_summary_dashboard(m2_lkf, m3_lkf, 'Linear Kalman Filter', output_dir)
        
        print("\n2. Generating error heatmap...")
        plot_error_heatmap(m2_lkf, m3_lkf, output_dir)
        
        print("\n3. Generating error distribution plots...")
        plot_error_distribution(m2_lkf, m3_lkf, output_dir)
        
        print("\n4. Generating individual state comparisons...")
        for state in ['px', 'vx', 'py', 'pz']:
            plot_state_comparison(m2_lkf, m3_lkf, 0, state, output_dir)
        
        print("\n5. Generating all joints comparison...")
        plot_all_joints_single_state(m2_lkf, m3_lkf, 'px', output_dir)
        
        print("\n6. Generating 3D trajectory...")
        plot_trajectory_3d(m2_lkf, m3_lkf, 0, output_dir)
        
        print(f"\n✓ LKF plots saved to {output_dir}/")
    
    if mode in ['ekf', 'both']:
        print("\n" + "="*70)
        print("  GENERATING EKF COMPARISON PLOTS")
        print("="*70 + "\n")
        
        m2_ekf, m3_ekf = load_data('ekf_output_m2.csv', 'ekf_output_m3.csv')
        
        output_dir = 'plots_ekf_comparison'
        os.makedirs(output_dir, exist_ok=True)
        
        print("\n1. Generating summary dashboard...")
        plot_summary_dashboard(m2_ekf, m3_ekf, 'Extended Kalman Filter', output_dir)
        
        print("\n2. Generating error heatmap...")
        plot_error_heatmap(m2_ekf, m3_ekf, output_dir)
        
        print("\n3. Generating error distribution plots...")
        plot_error_distribution(m2_ekf, m3_ekf, output_dir)
        
        print("\n4. Generating individual state comparisons...")
        for state in ['px', 'vx', 'py', 'pz']:
            plot_state_comparison(m2_ekf, m3_ekf, 0, state, output_dir)
        
        print("\n5. Generating all joints comparison...")
        plot_all_joints_single_state(m2_ekf, m3_ekf, 'px', output_dir)
        
        print("\n6. Generating 3D trajectory...")
        plot_trajectory_3d(m2_ekf, m3_ekf, 0, output_dir)
        
        print(f"\n✓ EKF plots saved to {output_dir}/")
    
    print("\n" + "="*70)
    print("  ALL COMPARISON PLOTS GENERATED SUCCESSFULLY")
    print("="*70 + "\n")

if __name__ == '__main__':
    main()