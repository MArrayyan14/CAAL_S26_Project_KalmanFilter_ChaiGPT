"""
Comprehensive Plotting Script - Milestone 2
Team ChaiGPT: Laiba (31736), Ameer (30535), Usman (29126), Arrayyan (30557)

Generates 8 plots for Joint 0:
1-4: Position, Velocity, Acceleration, Jerk comparisons
5-8: Error summaries for each
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from pathlib import Path
import sys

# Constants
STATES_PER_JOINT = 12
NUM_JOINTS = 23
DT = 0.01

def load_csv_skip_header(filename, expected_cols):
    """Load CSV and skip header row"""
    print(f"Loading {filename}...")
    try:
        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        print(f"  ✓ Loaded {data.shape[0]} data rows, {data.shape[1]} columns")
        
        if data.shape[0] != 3040:
            print(f"  ⚠ WARNING: Expected 3040 data rows, got {data.shape[0]}")
            
        return data
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        sys.exit(1)

def load_csv_no_skip(filename, expected_cols):
    """Load CSV without skipping header"""
    print(f"Loading {filename}...")
    try:
        data = np.loadtxt(filename, delimiter=',')
        print(f"  ✓ Loaded {data.shape[0]} rows, {data.shape[1]} columns")
        return data
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        sys.exit(1)

def extract_joint_states(states, joint_idx):
    """Extract all 12 states for a specific joint"""
    start = joint_idx * STATES_PER_JOINT
    return {
        'px': states[:, start + 0],  'vx': states[:, start + 1],
        'ax': states[:, start + 2],  'jx': states[:, start + 3],
        'py': states[:, start + 4],  'vy': states[:, start + 5],
        'ay': states[:, start + 6],  'jy': states[:, start + 7],
        'pz': states[:, start + 8],  'vz': states[:, start + 9],
        'az': states[:, start + 10], 'jz': states[:, start + 11]
    }

def extract_joint_positions(measurements, joint_idx):
    """Extract position measurements for a specific joint"""
    start = joint_idx * 3
    return {
        'px': measurements[:, start + 0],
        'py': measurements[:, start + 1],
        'pz': measurements[:, start + 2]
    }

# PLOT 1: POSITION - TRUE VS NOISY VS LKF VS EKF

def plot_position_comparison(joint_idx, true_pos, noisy_pos, lkf, ekf, time, output_dir):
    """Position: True vs Noisy vs LKF vs EKF"""
    print("  [1/8] Position - True vs Noisy vs LKF vs EKF...")
    
    fig, axes = plt.subplots(3, 1, figsize=(15, 11))
    fig.suptitle(f'Joint {joint_idx} - POSITION: True vs Noisy vs Estimated', 
                 fontsize=16, fontweight='bold')
    
    components = [
        ('px', 'X Position (m)', 'tab:blue'),
        ('py', 'Y Position (m)', 'tab:green'),
        ('pz', 'Z Position (m)', 'tab:red')
    ]
    
    for ax, (key, label, color) in zip(axes, components):
        ax.plot(time, true_pos[key], 'k-', label='True (Ground Truth)', 
                linewidth=2.5, alpha=0.9, zorder=4)
        ax.plot(time, noisy_pos[key], color='gray', label='Noisy (Measurements)', 
                linewidth=1, alpha=0.3, zorder=1)
        ax.plot(time, lkf[key], 'b-', label='LKF (Filtered)', 
                linewidth=2, alpha=0.85, zorder=3)
        ax.plot(time, ekf[key], 'r--', label='EKF (Filtered)', 
                linewidth=2, alpha=0.85, zorder=2)
        
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='upper right', fontsize=10, ncol=2, framealpha=0.9)
        
        noise_rmse = np.sqrt(np.mean((noisy_pos[key] - true_pos[key])**2)) * 1000
        lkf_rmse = np.sqrt(np.mean((lkf[key] - true_pos[key])**2)) * 1000
        ekf_rmse = np.sqrt(np.mean((ekf[key] - true_pos[key])**2)) * 1000
        lkf_reduction = ((noise_rmse - lkf_rmse) / noise_rmse) * 100
        ekf_reduction = ((noise_rmse - ekf_rmse) / noise_rmse) * 100
        
        textstr = (f'RMSE (mm) & Noise Reduction:\n'
                   f'Noise: {noise_rmse:.1f}\n'
                   f'LKF:   {lkf_rmse:.1f} ({lkf_reduction:+.1f}%)\n'
                   f'EKF:   {ekf_rmse:.1f} ({ekf_reduction:+.1f}%)')
        
        ax.text(0.02, 0.97, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    axes[-1].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / f'joint_{joint_idx}_1_position_comparison.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# PLOT 2: VELOCITY - LKF VS EKF

def plot_velocity_comparison(joint_idx, lkf, ekf, time, output_dir):
    """Velocity: LKF vs EKF"""
    print("  [2/8] Velocity - LKF vs EKF...")
    
    fig, axes = plt.subplots(3, 1, figsize=(15, 11))
    fig.suptitle(f'Joint {joint_idx} - VELOCITY: LKF vs EKF', 
                 fontsize=16, fontweight='bold')
    
    components = [
        ('vx', 'X Velocity (m/s)', 'tab:blue'),
        ('vy', 'Y Velocity (m/s)', 'tab:green'),
        ('vz', 'Z Velocity (m/s)', 'tab:red')
    ]
    
    for ax, (key, label, color) in zip(axes, components):
        ax.plot(time, lkf[key], 'b-', label='LKF', 
                linewidth=2, alpha=0.85, zorder=3)
        ax.plot(time, ekf[key], 'r--', label='EKF', 
                linewidth=2, alpha=0.85, zorder=2)
        
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='upper right', fontsize=10, ncol=2, framealpha=0.9)
        
        diff_rmse = np.sqrt(np.mean((lkf[key] - ekf[key])**2)) * 1000
        lkf_mean = np.mean(lkf[key]) * 1000
        ekf_mean = np.mean(ekf[key]) * 1000
        lkf_std = np.std(lkf[key]) * 1000
        ekf_std = np.std(ekf[key]) * 1000
        
        textstr = (f'Statistics (mm/s):\n'
                   f'LKF: μ={lkf_mean:.1f}, σ={lkf_std:.1f}\n'
                   f'EKF: μ={ekf_mean:.1f}, σ={ekf_std:.1f}\n'
                   f'Diff RMSE: {diff_rmse:.1f}')
        
        ax.text(0.02, 0.97, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    axes[-1].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / f'joint_{joint_idx}_2_velocity_comparison.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# PLOT 3: ACCELERATION - LKF VS EKF

def plot_acceleration_comparison(joint_idx, lkf, ekf, time, output_dir):
    """Acceleration: LKF vs EKF"""
    print("  [3/8] Acceleration - LKF vs EKF...")
    
    fig, axes = plt.subplots(3, 1, figsize=(15, 11))
    fig.suptitle(f'Joint {joint_idx} - ACCELERATION: LKF vs EKF', 
                 fontsize=16, fontweight='bold')
    
    components = [
        ('ax', 'X Acceleration (m/s²)', 'tab:blue'),
        ('ay', 'Y Acceleration (m/s²)', 'tab:green'),
        ('az', 'Z Acceleration (m/s²)', 'tab:red')
    ]
    
    for ax, (key, label, color) in zip(axes, components):
        ax.plot(time, lkf[key], 'b-', label='LKF', 
                linewidth=2, alpha=0.85, zorder=3)
        ax.plot(time, ekf[key], 'r--', label='EKF', 
                linewidth=2, alpha=0.85, zorder=2)
        
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='upper right', fontsize=10, ncol=2, framealpha=0.9)
        
        diff_rmse = np.sqrt(np.mean((lkf[key] - ekf[key])**2)) * 1000
        lkf_mean = np.mean(lkf[key]) * 1000
        ekf_mean = np.mean(ekf[key]) * 1000
        lkf_std = np.std(lkf[key]) * 1000
        ekf_std = np.std(ekf[key]) * 1000
        
        textstr = (f'Statistics (mm/s²):\n'
                   f'LKF: μ={lkf_mean:.1f}, σ={lkf_std:.1f}\n'
                   f'EKF: μ={ekf_mean:.1f}, σ={ekf_std:.1f}\n'
                   f'Diff RMSE: {diff_rmse:.1f}')
        
        ax.text(0.02, 0.97, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    axes[-1].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / f'joint_{joint_idx}_3_acceleration_comparison.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# PLOT 4: JERK - LKF VS EKF

def plot_jerk_comparison(joint_idx, lkf, ekf, time, output_dir):
    """Jerk: LKF vs EKF"""
    print("  [4/8] Jerk - LKF vs EKF...")
    
    fig, axes = plt.subplots(3, 1, figsize=(15, 11))
    fig.suptitle(f'Joint {joint_idx} - JERK: LKF vs EKF', 
                 fontsize=16, fontweight='bold')
    
    components = [
        ('jx', 'X Jerk (m/s³)', 'tab:blue'),
        ('jy', 'Y Jerk (m/s³)', 'tab:green'),
        ('jz', 'Z Jerk (m/s³)', 'tab:red')
    ]
    
    for ax, (key, label, color) in zip(axes, components):
        ax.plot(time, lkf[key], 'b-', label='LKF', 
                linewidth=2, alpha=0.85, zorder=3)
        ax.plot(time, ekf[key], 'r--', label='EKF', 
                linewidth=2, alpha=0.85, zorder=2)
        
        ax.set_ylabel(label, fontsize=12, fontweight='bold')
        ax.grid(True, alpha=0.3, linestyle='--')
        ax.legend(loc='upper right', fontsize=10, ncol=2, framealpha=0.9)
        
        diff_rmse = np.sqrt(np.mean((lkf[key] - ekf[key])**2)) * 1000
        lkf_mean = np.mean(lkf[key]) * 1000
        ekf_mean = np.mean(ekf[key]) * 1000
        lkf_std = np.std(lkf[key]) * 1000
        ekf_std = np.std(ekf[key]) * 1000
        
        textstr = (f'Statistics (mm/s³):\n'
                   f'LKF: μ={lkf_mean:.1f}, σ={lkf_std:.1f}\n'
                   f'EKF: μ={ekf_mean:.1f}, σ={ekf_std:.1f}\n'
                   f'Diff RMSE: {diff_rmse:.1f}')
        
        ax.text(0.02, 0.97, textstr, transform=ax.transAxes, fontsize=9,
                verticalalignment='top', family='monospace',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    axes[-1].set_xlabel('Time (s)', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_dir / f'joint_{joint_idx}_4_jerk_comparison.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# PLOT 5: POSITION ERROR SUMMARY

def plot_position_error_summary(joint_idx, lkf, ekf, true_pos, noisy_pos, time, output_dir):
    """Position Error Summary"""
    print("  [5/8] Position error summary...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    fig.suptitle(f'Joint {joint_idx} - POSITION Error Analysis', 
                 fontsize=16, fontweight='bold')
    
    ax1 = fig.add_subplot(gs[0, :])
    keys = ['px', 'py', 'pz']
    names = ['X', 'Y', 'Z']
    colors = ['tab:blue', 'tab:green', 'tab:red']
    
    for key, name, color in zip(keys, names, colors):
        lkf_err = np.abs(lkf[key] - true_pos[key]) * 1000
        ekf_err = np.abs(ekf[key] - true_pos[key]) * 1000
        
        ax1.plot(time, lkf_err, color=color, linestyle='-', 
                linewidth=1.5, alpha=0.7, label=f'{name}-LKF')
        ax1.plot(time, ekf_err, color=color, linestyle='--', 
                linewidth=1.5, alpha=0.7, label=f'{name}-EKF')
    
    ax1.set_ylabel('Absolute Error (mm)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Time (s)', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', ncol=3, fontsize=9)
    ax1.set_title('Position Errors Over Time', fontsize=13, fontweight='bold')
    
    ax2 = fig.add_subplot(gs[1, 0])
    
    noise_rmses = [np.sqrt(np.mean((noisy_pos[k] - true_pos[k])**2)) * 1000 for k in keys]
    lkf_rmses = [np.sqrt(np.mean((lkf[k] - true_pos[k])**2)) * 1000 for k in keys]
    ekf_rmses = [np.sqrt(np.mean((ekf[k] - true_pos[k])**2)) * 1000 for k in keys]
    
    x = np.arange(len(names))
    width = 0.25
    
    ax2.bar(x - width, noise_rmses, width, label='Noisy', color='gray', alpha=0.7)
    ax2.bar(x, lkf_rmses, width, label='LKF', color='tab:blue', alpha=0.8)
    ax2.bar(x + width, ekf_rmses, width, label='EKF', color='tab:red', alpha=0.8)
    
    ax2.set_ylabel('RMSE (mm)', fontsize=11, fontweight='bold')
    ax2.set_title('Position RMSE by Axis', fontsize=12, fontweight='bold')
    ax2.set_xticks(x)
    ax2.set_xticklabels(names)
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    summary = "POSITION ERROR STATISTICS\n" + "="*50 + "\n\n"
    for name, key in zip(names, keys):
        noise_rmse = np.sqrt(np.mean((noisy_pos[key] - true_pos[key])**2)) * 1000
        lkf_rmse = np.sqrt(np.mean((lkf[key] - true_pos[key])**2)) * 1000
        ekf_rmse = np.sqrt(np.mean((ekf[key] - true_pos[key])**2)) * 1000
        lkf_red = ((noise_rmse - lkf_rmse) / noise_rmse * 100)
        ekf_red = ((noise_rmse - ekf_rmse) / noise_rmse * 100)
        
        summary += f"{name}-Axis:\n"
        summary += f"  Noise:  {noise_rmse:6.2f} mm\n"
        summary += f"  LKF:    {lkf_rmse:6.2f} mm  ({lkf_red:5.1f}% reduction)\n"
        summary += f"  EKF:    {ekf_rmse:6.2f} mm  ({ekf_red:5.1f}% reduction)\n\n"
    
    overall_noise = np.mean(noise_rmses)
    overall_lkf = np.mean(lkf_rmses)
    overall_ekf = np.mean(ekf_rmses)
    
    summary += f"Overall Average:\n"
    summary += f"  Noise:  {overall_noise:.2f} mm\n"
    summary += f"  LKF:    {overall_lkf:.2f} mm\n"
    summary += f"  EKF:    {overall_ekf:.2f} mm\n"
    
    ax3.text(0.1, 0.95, summary, transform=ax3.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.savefig(output_dir / f'joint_{joint_idx}_5_position_error_summary.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# PLOT 6: VELOCITY ERROR SUMMARY

def plot_velocity_error_summary(joint_idx, lkf, ekf, time, output_dir):
    """Velocity Error Summary"""
    print("  [6/8] Velocity error summary...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    fig.suptitle(f'Joint {joint_idx} - VELOCITY Comparison Analysis', 
                 fontsize=16, fontweight='bold')
    
    ax1 = fig.add_subplot(gs[0, :])
    keys = ['vx', 'vy', 'vz']
    names = ['X', 'Y', 'Z']
    colors = ['tab:blue', 'tab:green', 'tab:red']
    
    for key, name, color in zip(keys, names, colors):
        diff = np.abs(lkf[key] - ekf[key]) * 1000
        ax1.plot(time, diff, color=color, linewidth=1.5, alpha=0.7, label=f'{name}-axis')
    
    ax1.set_ylabel('|LKF - EKF| Difference (mm/s)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Time (s)', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.set_title('Velocity Difference Between LKF and EKF', fontsize=13, fontweight='bold')
    
    ax2 = fig.add_subplot(gs[1, 0])
    
    diff_rmses = [np.sqrt(np.mean((lkf[k] - ekf[k])**2)) * 1000 for k in keys]
    
    ax2.bar(names, diff_rmses, color=colors, alpha=0.8)
    ax2.set_ylabel('Difference RMSE (mm/s)', fontsize=11, fontweight='bold')
    ax2.set_title('Velocity RMSE Difference by Axis', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    summary = "VELOCITY STATISTICS\n" + "="*50 + "\n\n"
    for name, key in zip(names, keys):
        lkf_mean = np.mean(lkf[key]) * 1000
        ekf_mean = np.mean(ekf[key]) * 1000
        lkf_std = np.std(lkf[key]) * 1000
        ekf_std = np.std(ekf[key]) * 1000
        diff_rmse = np.sqrt(np.mean((lkf[key] - ekf[key])**2)) * 1000
        correlation = np.corrcoef(lkf[key], ekf[key])[0, 1]
        
        summary += f"{name}-Axis (mm/s):\n"
        summary += f"  LKF:  μ={lkf_mean:6.2f}, σ={lkf_std:6.2f}\n"
        summary += f"  EKF:  μ={ekf_mean:6.2f}, σ={ekf_std:6.2f}\n"
        summary += f"  Diff RMSE: {diff_rmse:6.2f}\n"
        summary += f"  Correlation: {correlation:.4f}\n\n"
    
    ax3.text(0.1, 0.95, summary, transform=ax3.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.savefig(output_dir / f'joint_{joint_idx}_6_velocity_error_summary.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# PLOT 7: ACCELERATION ERROR SUMMARY

def plot_acceleration_error_summary(joint_idx, lkf, ekf, time, output_dir):
    """Acceleration Error Summary"""
    print("  [7/8] Acceleration error summary...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    fig.suptitle(f'Joint {joint_idx} - ACCELERATION Comparison Analysis', 
                 fontsize=16, fontweight='bold')
    
    ax1 = fig.add_subplot(gs[0, :])
    keys = ['ax', 'ay', 'az']
    names = ['X', 'Y', 'Z']
    colors = ['tab:blue', 'tab:green', 'tab:red']
    
    for key, name, color in zip(keys, names, colors):
        diff = np.abs(lkf[key] - ekf[key]) * 1000
        ax1.plot(time, diff, color=color, linewidth=1.5, alpha=0.7, label=f'{name}-axis')
    
    ax1.set_ylabel('|LKF - EKF| Difference (mm/s²)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Time (s)', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.set_title('Acceleration Difference Between LKF and EKF', fontsize=13, fontweight='bold')
    
    ax2 = fig.add_subplot(gs[1, 0])
    
    diff_rmses = [np.sqrt(np.mean((lkf[k] - ekf[k])**2)) * 1000 for k in keys]
    
    ax2.bar(names, diff_rmses, color=colors, alpha=0.8)
    ax2.set_ylabel('Difference RMSE (mm/s²)', fontsize=11, fontweight='bold')
    ax2.set_title('Acceleration RMSE Difference by Axis', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    summary = "ACCELERATION STATISTICS\n" + "="*50 + "\n\n"
    for name, key in zip(names, keys):
        lkf_mean = np.mean(lkf[key]) * 1000
        ekf_mean = np.mean(ekf[key]) * 1000
        lkf_std = np.std(lkf[key]) * 1000
        ekf_std = np.std(ekf[key]) * 1000
        diff_rmse = np.sqrt(np.mean((lkf[key] - ekf[key])**2)) * 1000
        correlation = np.corrcoef(lkf[key], ekf[key])[0, 1]
        
        summary += f"{name}-Axis (mm/s²):\n"
        summary += f"  LKF:  μ={lkf_mean:6.2f}, σ={lkf_std:6.2f}\n"
        summary += f"  EKF:  μ={ekf_mean:6.2f}, σ={ekf_std:6.2f}\n"
        summary += f"  Diff RMSE: {diff_rmse:6.2f}\n"
        summary += f"  Correlation: {correlation:.4f}\n\n"
    
    ax3.text(0.1, 0.95, summary, transform=ax3.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
    
    plt.savefig(output_dir / f'joint_{joint_idx}_7_acceleration_error_summary.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# PLOT 8: JERK ERROR SUMMARY

def plot_jerk_error_summary(joint_idx, lkf, ekf, time, output_dir):
    """Jerk Error Summary"""
    print("  [8/8] Jerk error summary...")
    
    fig = plt.figure(figsize=(16, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.3, wspace=0.3)
    fig.suptitle(f'Joint {joint_idx} - JERK Comparison Analysis', 
                 fontsize=16, fontweight='bold')
    
    ax1 = fig.add_subplot(gs[0, :])
    keys = ['jx', 'jy', 'jz']
    names = ['X', 'Y', 'Z']
    colors = ['tab:blue', 'tab:green', 'tab:red']
    
    for key, name, color in zip(keys, names, colors):
        diff = np.abs(lkf[key] - ekf[key]) * 1000
        ax1.plot(time, diff, color=color, linewidth=1.5, alpha=0.7, label=f'{name}-axis')
    
    ax1.set_ylabel('|LKF - EKF| Difference (mm/s³)', fontsize=12, fontweight='bold')
    ax1.set_xlabel('Time (s)', fontsize=11)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right', fontsize=10)
    ax1.set_title('Jerk Difference Between LKF and EKF', fontsize=13, fontweight='bold')
    
    ax2 = fig.add_subplot(gs[1, 0])
    
    diff_rmses = [np.sqrt(np.mean((lkf[k] - ekf[k])**2)) * 1000 for k in keys]
    
    ax2.bar(names, diff_rmses, color=colors, alpha=0.8)
    ax2.set_ylabel('Difference RMSE (mm/s³)', fontsize=11, fontweight='bold')
    ax2.set_title('Jerk RMSE Difference by Axis', fontsize=12, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    ax3 = fig.add_subplot(gs[1, 1])
    ax3.axis('off')
    
    summary = "JERK STATISTICS\n" + "="*50 + "\n\n"
    for name, key in zip(names, keys):
        lkf_mean = np.mean(lkf[key]) * 1000
        ekf_mean = np.mean(ekf[key]) * 1000
        lkf_std = np.std(lkf[key]) * 1000
        ekf_std = np.std(ekf[key]) * 1000
        diff_rmse = np.sqrt(np.mean((lkf[key] - ekf[key])**2)) * 1000
        correlation = np.corrcoef(lkf[key], ekf[key])[0, 1]
        
        summary += f"{name}-Axis (mm/s³):\n"
        summary += f"  LKF:  μ={lkf_mean:6.2f}, σ={lkf_std:6.2f}\n"
        summary += f"  EKF:  μ={ekf_mean:6.2f}, σ={ekf_std:6.2f}\n"
        summary += f"  Diff RMSE: {diff_rmse:6.2f}\n"
        summary += f"  Correlation: {correlation:.4f}\n\n"
    
    ax3.text(0.1, 0.95, summary, transform=ax3.transAxes,
            fontsize=9, verticalalignment='top', family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
    
    plt.savefig(output_dir / f'joint_{joint_idx}_8_jerk_error_summary.png', 
                dpi=300, bbox_inches='tight')
    plt.close()

# MAIN

def main():
    parser = argparse.ArgumentParser(description='Generate 8 plots for Milestone 2')
    
    parser.add_argument('--joint', type=int, default=0, help='Joint index (0-22)')
    parser.add_argument('--lkf', type=str, default='lkf_output.csv')
    parser.add_argument('--ekf', type=str, default='ekf_output.csv')
    parser.add_argument('--true', type=str, default='gait_data_true.csv')
    parser.add_argument('--noisy', type=str, default='gait_data_noisy.csv')
    parser.add_argument('--output', type=str, default='plots')
    
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("  MILESTONE 2 PLOTTING - Team ChaiGPT")
    print("="*70)
    print(f"\nJoint: {args.joint}")
    print(f"Output: {args.output}/\n")
    
    # Create output directory
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    print("Loading datasets...")
    true_all = load_csv_skip_header(args.true, 69)
    noisy_all = load_csv_skip_header(args.noisy, 69)
    lkf_all = load_csv_no_skip(args.lkf, 276)
    ekf_all = load_csv_no_skip(args.ekf, 276)
    
    # Verify row counts
    if not (true_all.shape[0] == noisy_all.shape[0] == lkf_all.shape[0] == ekf_all.shape[0]):
        print("\n ERROR: Row count mismatch!")
        sys.exit(1)
    
    print(f"  ✓ All files have {lkf_all.shape[0]} rows - MATCH!\n")
    
    # Extract joint data
    true_pos = extract_joint_positions(true_all, args.joint)
    noisy_pos = extract_joint_positions(noisy_all, args.joint)
    lkf = extract_joint_states(lkf_all, args.joint)
    ekf = extract_joint_states(ekf_all, args.joint)
    
    # Time vector
    time = np.arange(lkf_all.shape[0]) * DT
    
    # Generate plots
    print("="*70)
    print("Generating 8 plots...")
    print("="*70 + "\n")
    
    plot_position_comparison(args.joint, true_pos, noisy_pos, lkf, ekf, time, output_dir)
    plot_velocity_comparison(args.joint, lkf, ekf, time, output_dir)
    plot_acceleration_comparison(args.joint, lkf, ekf, time, output_dir)
    plot_jerk_comparison(args.joint, lkf, ekf, time, output_dir)
    plot_position_error_summary(args.joint, lkf, ekf, true_pos, noisy_pos, time, output_dir)
    plot_velocity_error_summary(args.joint, lkf, ekf, time, output_dir)
    plot_acceleration_error_summary(args.joint, lkf, ekf, time, output_dir)
    plot_jerk_error_summary(args.joint, lkf, ekf, time, output_dir)
    
    print("\n" + "="*70)
    print("ALL 8 PLOTS GENERATED SUCCESSFULLY!")
    print("="*70 + "\n")
    print("Output files:")
    print("  1. joint_0_1_position_comparison.png")
    print("  2. joint_0_2_velocity_comparison.png")
    print("  3. joint_0_3_acceleration_comparison.png")
    print("  4. joint_0_4_jerk_comparison.png")
    print("  5. joint_0_5_position_error_summary.png")
    print("  6. joint_0_6_velocity_error_summary.png")
    print("  7. joint_0_7_acceleration_error_summary.png")
    print("  8. joint_0_8_jerk_error_summary.png\n")

if __name__ == '__main__':
    main()
