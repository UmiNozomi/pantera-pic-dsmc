#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
H-alpha 光谱累积显示工具 - 将所有时间步的counts相加

用法：
    python plot_spectrum_cumulative.py
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import glob
import os

# 设置字体和样式
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans', 'Arial']
rcParams['axes.unicode_minus'] = False
rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 300

# H-alpha 静止波长
LAMBDA0 = 656.28  # nm

def plot_cumulative_spectrum():
    """读取并累积所有时间步的光谱"""
    
    # 查找所有光谱文件
    files = sorted(glob.glob('results/spectrum_halpha_*.dat'))
    
    if not files:
        print("✗ 未找到光谱文件！")
        print("  请检查 results/ 目录")
        return
    
    print(f"找到 {len(files)} 个光谱文件\n")
    
    # 初始化累积数组
    wavelength = None
    vlos = None
    cumulative_counts = None
    total_events = 0
    
    # 读取并累加每个文件
    for i, filename in enumerate(files):
        try:
            # 提取时间步
            import re
            match = re.search(r'spectrum_halpha_(\d+)\.dat', filename)
            timestep = int(match.group(1)) if match else i
            
            # 读取数据
            data = np.loadtxt(filename, comments='#')
            vlos_temp = data[:, 0]  # m/s
            wavelength_temp = data[:, 1]  # nm
            counts_temp = data[:, 2]
            
            # 读取总事件数
            events = 0
            with open(filename, 'r') as f:
                for line in f:
                    if 'Total events:' in line:
                        events = int(line.split(':')[1].strip())
                        break
            
            if events == 0:
                events = int(np.sum(counts_temp))
            
            # 初始化或累加
            if cumulative_counts is None:
                wavelength = wavelength_temp
                vlos = vlos_temp
                cumulative_counts = counts_temp.copy()
            else:
                cumulative_counts += counts_temp
            
            total_events += events
            
            print(f"✓ {os.path.basename(filename):40s} | t={timestep:6d} | events={events:8,d}")
            
        except Exception as e:
            print(f"✗ 读取失败: {filename}")
            print(f"  错误: {e}")
            continue
    
    if cumulative_counts is None:
        print("✗ 没有有效数据！")
        return
    
    # 归一化
    intensity = cumulative_counts / np.sum(cumulative_counts)
    
    print(f"\n累积总事件数: {total_events:,}")
    print(f"累积总counts: {np.sum(cumulative_counts):,.0f}\n")
    
    # 创建图形
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'H-alpha 累积光谱 (总事件数: {total_events:,})', 
                 fontsize=16, fontweight='bold')
    
    # 左图：波长空间
    ax1 = axes[0]
    ax1.plot(wavelength, intensity, 'b-', linewidth=2.5, label='PIC-MCC 模拟')
    ax1.axvline(LAMBDA0, color='red', linestyle='--', linewidth=2, 
                alpha=0.6, label=f'Hα = {LAMBDA0} nm')
    ax1.set_xlabel('Wavelength [nm]', fontsize=14, fontweight='bold')
    ax1.set_ylabel('Intensity [arb.u]', fontsize=14, fontweight='bold')
    ax1.set_title('波长光谱', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=12, loc='best', framealpha=0.9)
    ax1.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax1.tick_params(labelsize=12)
    ax1.set_xlim([wavelength.min(), wavelength.max()])
    
    # 右图：速度空间
    ax2 = axes[1]
    ax2.plot(vlos/1e3, intensity, 'b-', linewidth=2.5, label='PIC-MCC 模拟')
    ax2.axvline(0, color='red', linestyle='--', linewidth=2,
                alpha=0.6, label='v = 0')
    ax2.set_xlabel('Velocity [km/s]', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Intensity [arb.u]', fontsize=14, fontweight='bold')
    ax2.set_title('速度光谱', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=12, loc='best', framealpha=0.9)
    ax2.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
    ax2.tick_params(labelsize=12)
    ax2.set_xlim([vlos.min()/1e3, vlos.max()/1e3])
    
    plt.tight_layout()
    
    # 保存图片
    output_file = 'halpha_cumulative.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ 图片已保存: {output_file}")
    
    # 保存数据
    data_file = 'halpha_cumulative.txt'
    header = f'Cumulative H-alpha spectrum from {len(files)} timesteps\n'
    header += f'Total events: {total_events}\n'
    header += 'Velocity[m/s]  Wavelength[nm]  Counts  Intensity'
    np.savetxt(data_file, 
               np.column_stack([vlos, wavelength, cumulative_counts, intensity]),
               header=header, 
               fmt='%.6e  %.6f  %.1f  %.6e')
    print(f"✓ 数据已保存: {data_file}")
    
    plt.show()
    
    # 计算并输出一些统计信息
    print("\n" + "="*70)
    print("  光谱统计信息")
    print("="*70)
    
    peak_idx = np.argmax(intensity)
    peak_wavelength = wavelength[peak_idx]
    peak_velocity = vlos[peak_idx]
    
    # 计算FWHM
    half_max = intensity[peak_idx] / 2
    left_idx = peak_idx
    while left_idx > 0 and intensity[left_idx] > half_max:
        left_idx -= 1
    right_idx = peak_idx
    while right_idx < len(intensity) - 1 and intensity[right_idx] > half_max:
        right_idx += 1
    
    fwhm_nm = wavelength[right_idx] - wavelength[left_idx]
    fwhm_velocity = vlos[right_idx] - vlos[left_idx]
    
    print(f"\n峰值波长:     {peak_wavelength:.4f} nm")
    print(f"峰值速度:     {peak_velocity/1e3:.2f} km/s")
    print(f"波长频移:     {peak_wavelength - LAMBDA0:.4f} nm")
    print(f"FWHM (波长):  {fwhm_nm:.4f} nm")
    print(f"FWHM (速度):  {fwhm_velocity/1e3:.2f} km/s")
    print(f"\n处理文件数:   {len(files)}")
    print(f"累积事件数:   {total_events:,}")
    print()


if __name__ == '__main__':
    print("\n" + "="*70)
    print("  H-alpha 累积光谱分析")
    print("="*70 + "\n")
    
    plot_cumulative_spectrum()
    
    print("="*70)
    print("  完成！")
    print("="*70 + "\n")

