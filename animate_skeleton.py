#!/usr/bin/env python3
"""
3D Full-Body Animation Script - Milestone 2
Team ChaiGPT: Laiba (31736), Ameer (30535), Usman (29126), Arrayyan (30557)
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import argparse
import sys

# Skeleton topology (joint connections)
SKELETON_CONNECTIONS = [
    (0, 1), (1, 2), (2, 3),           # Spine
    (3, 4), (4, 5), (5, 6), (6, 7),   # Right arm
    (3, 8), (8, 9), (9, 10), (10, 11), # Left arm
    (0, 12), (12, 13), (13, 14), (14, 15), # Right leg
    (0, 16), (16, 17), (17, 18), (18, 19), # Left leg
    (3, 20), (20, 21), (21, 22)       # Head
]

def load_csv_skip_header(filename, expected_cols):
    """Load CSV and skip header row"""
    print(f"Loading {filename}...")
    try:
        data = np.loadtxt(filename, delimiter=',', skiprows=1)
        print(f"  ✓ Loaded {data.shape[0]} data rows")
        
        if data.shape[0] != 3040:
            print(f"  ⚠ WARNING: Expected 3040 rows, got {data.shape[0]}")
        
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
        
        if data.shape[1] != expected_cols:
            print(f"  ⚠ WARNING: Expected {expected_cols} columns, got {data.shape[1]}")
            
        return data
    except Exception as e:
        print(f"  ✗ ERROR: {e}")
        sys.exit(1)

def extract_positions(data, is_state=False):
    """
    Extract 3D positions for all joints
    data: (frames, cols)
    Returns: (frames, 23, 3) array
    """
    num_frames = data.shape[0]
    positions = np.zeros((num_frames, 23, 3))
    
    if is_state:
        # State vector: extract px, py, pz from 12-state blocks
        for joint in range(23):
            state_base = joint * 12
            positions[:, joint, 0] = data[:, state_base + 0]  # px
            positions[:, joint, 1] = data[:, state_base + 4]  # py
            positions[:, joint, 2] = data[:, state_base + 8]  # pz
    else:
        # Measurement vector: already just positions
        for joint in range(23):
            meas_base = joint * 3
            positions[:, joint, 0] = data[:, meas_base + 0]
            positions[:, joint, 1] = data[:, meas_base + 1]
            positions[:, joint, 2] = data[:, meas_base + 2]
    
    return positions

def create_animation(true_pos, noisy_pos, lkf_pos, ekf_pos, output_file, fps=30, duration=10):
    """
    Create side-by-side animation of FOUR skeletons (True, Noisy, LKF, EKF)
    """
    num_frames = true_pos.shape[0]
    frame_step = max(1, num_frames // (fps * duration))
    frames_to_animate = num_frames // frame_step
    
    print(f"\nCreating animation:")
    print(f"  Total frames: {num_frames}")
    print(f"  Frames to animate: {frames_to_animate}")
    print(f"  FPS: {fps}")
    print(f"  Duration: ~{frames_to_animate/fps:.1f} seconds")
    
    fig = plt.figure(figsize=(24, 6))
    
    # Four subplots
    ax1 = fig.add_subplot(141, projection='3d')
    ax2 = fig.add_subplot(142, projection='3d')
    ax3 = fig.add_subplot(143, projection='3d')
    ax4 = fig.add_subplot(144, projection='3d')
    
    axes = [ax1, ax2, ax3, ax4]
    titles = ['True', 'Noisy', 'LKF', 'EKF']
    datasets = [true_pos, noisy_pos, lkf_pos, ekf_pos]
    colors = ['black', 'gray', 'blue', 'red']
    
    # Calculate axis limits from ALL datasets with padding
    all_data = np.concatenate([true_pos, noisy_pos, lkf_pos, ekf_pos], axis=0)
    
    x_min, x_max = all_data[:, :, 0].min(), all_data[:, :, 0].max()
    y_min, y_max = all_data[:, :, 1].min(), all_data[:, :, 1].max()
    z_min, z_max = all_data[:, :, 2].min(), all_data[:, :, 2].max()
    
    # Add 15% padding to make skeleton visible and centered
    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min
    
    x_pad = x_range * 0.15
    y_pad = y_range * 0.15
    z_pad = z_range * 0.15
    
    x_lim = [x_min - x_pad, x_max + x_pad]
    y_lim = [y_min - y_pad, y_max + y_pad]
    z_lim = [z_min - z_pad, z_max + z_pad]
    
    print(f"\nAxis limits (with 15% padding):")
    print(f"  X: [{x_lim[0]:.2f}, {x_lim[1]:.2f}] m")
    print(f"  Y: [{y_lim[0]:.2f}, {y_lim[1]:.2f}] m")
    print(f"  Z: [{z_lim[0]:.2f}, {z_lim[1]:.2f}] m")
    
    # Set up axes
    for ax, title in zip(axes, titles):
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title(title, fontsize=14, fontweight='bold')
        
        # Set dynamic limits for all axes
        ax.set_xlim(x_lim)
        ax.set_ylim(y_lim)
        ax.set_zlim(z_lim)
    
    # Initialize lines
    lines = []
    for ax, color in zip(axes, colors):
        line_objs = []
        for connection in SKELETON_CONNECTIONS:
            line, = ax.plot([], [], [], color=color, linewidth=2, marker='o', markersize=4)
            line_objs.append(line)
        lines.append(line_objs)
    
    def init():
        for ax_lines in lines:
            for line in ax_lines:
                line.set_data([], [])
                line.set_3d_properties([])
        return [line for ax_lines in lines for line in ax_lines]
    
    def animate(frame_idx):
        real_frame = frame_idx * frame_step
        
        for ax_idx, (ax_lines, positions) in enumerate(zip(lines, datasets)):
            for line_idx, (j1, j2) in enumerate(SKELETON_CONNECTIONS):
                x = [positions[real_frame, j1, 0], positions[real_frame, j2, 0]]
                y = [positions[real_frame, j1, 1], positions[real_frame, j2, 1]]
                z = [positions[real_frame, j1, 2], positions[real_frame, j2, 2]]
                
                ax_lines[line_idx].set_data(x, y)
                ax_lines[line_idx].set_3d_properties(z)
        
        fig.suptitle(f'Frame {real_frame}/{num_frames} - Time: {real_frame*0.01:.2f}s', 
                     fontsize=16, fontweight='bold')
        
        return [line for ax_lines in lines for line in ax_lines]
    
    print("\nRendering animation...")
    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                   frames=frames_to_animate, interval=1000/fps,
                                   blit=True)
    
    print(f"Saving to {output_file}...")
    
    if output_file.endswith('.mp4'):
        writer = animation.FFMpegWriter(fps=fps, bitrate=2000)
        anim.save(output_file, writer=writer)
    elif output_file.endswith('.gif'):
        anim.save(output_file, writer='pillow', fps=fps)
    else:
        print("⚠ Unknown format, saving as MP4")
        anim.save(output_file + '.mp4', writer=animation.FFMpegWriter(fps=fps, bitrate=2000))
    
    plt.close()
    print(f"  ✓ Animation saved successfully!")

def main():
    parser = argparse.ArgumentParser(description='Create 3D walking animation')
    parser.add_argument('--true', type=str, default='gait_data_true.csv')
    parser.add_argument('--noisy', type=str, default='gait_data_noisy.csv')
    parser.add_argument('--lkf', type=str, default='lkf_output.csv')
    parser.add_argument('--ekf', type=str, default='ekf_output.csv')
    parser.add_argument('--output', type=str, default='animation.mp4')
    parser.add_argument('--fps', type=int, default=15, help='Frames per second (default: 15 for slower animation)')
    parser.add_argument('--duration', type=int, default=20, help='Duration in seconds (default: 20 for slower animation)')
    parser.add_argument('--format', type=str, default='mp4', choices=['mp4', 'gif'])
    
    args = parser.parse_args()
    
    print("\n" + "="*70)
    print("  3D FULL-BODY ANIMATION - Milestone 2")
    print("  Team ChaiGPT")
    print("="*70)
    
    # Load data (skips headers)
    print("\nLoading datasets (skipping header rows)...")
    true_data = load_csv_skip_header(args.true, 69)
    noisy_data = load_csv_skip_header(args.noisy, 69)
    lkf_data = load_csv_no_skip(args.lkf, 276)
    ekf_data = load_csv_no_skip(args.ekf, 276)
    
    # Verify row counts
    if not (true_data.shape[0] == noisy_data.shape[0] == lkf_data.shape[0] == ekf_data.shape[0]):
        print("\n⚠ ERROR: Row count mismatch!")
        sys.exit(1)
    
    print(f"  ✓ All files have {true_data.shape[0]} rows - MATCH!")
    
    # Extract 3D positions
    print("\nExtracting 3D positions for 23 joints...")
    true_pos = extract_positions(true_data, is_state=False)
    noisy_pos = extract_positions(noisy_data, is_state=False)
    lkf_pos = extract_positions(lkf_data, is_state=True)
    ekf_pos = extract_positions(ekf_data, is_state=True)
    print(f"  ✓ Extracted positions: {true_pos.shape}")
    
    # Create animation
    output_file = f"{args.output.rsplit('.', 1)[0]}.{args.format}"
    create_animation(true_pos, noisy_pos, lkf_pos, ekf_pos, output_file, args.fps, args.duration)
    
    print("\n" + "="*70)
    print("ANIMATION COMPLETE!")
    print("="*70)
    print(f"\nGenerated file: {output_file}\n")

if __name__ == '__main__':
    main()
