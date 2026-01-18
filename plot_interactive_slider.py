#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
交互式光谱滑块查看器
Interactive Spectrum Slider Viewer
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from scipy.interpolate import make_interp_spline
from scipy.signal import savgol_filter
from matplotlib.widgets import Slider

# 设置UTF-8编码以支持中文输出
if sys.platform == 'win32':
    os.environ['PYTHONIOENCODING'] = 'utf-8'
    sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def read_matrix_spectrum(filename='results/spectrum_halpha_cumulative.dat'):
    """读取矩阵格式的光谱数据"""
    data = pd.read_csv(filename, sep=r'\s+', comment='#')
    wavelengths = data.columns[1:].astype(float).values
    
    print(f"读取了 {len(data)} 个时间步")
    print(f"波长点数: {len(wavelengths)}")
    print(f"时间步范围: {data['timestep'].min()} - {data['timestep'].max()}")
    print(f"波长范围: {wavelengths[0]:.4f} - {wavelengths[-1]:.4f} nm")
    
    return data, wavelengths


def plot_interactive_spectrum_slider(filename='results/spectrum_halpha_cumulative.dat'):
    """
    创建交互式光谱图，可以通过滑块查看不同时间步的光谱
    """
    print("\n" + "=" * 70)
    print("绘制交互式时间滑块光谱图...")
    print("=" * 70)
    
    # 读取数据
    data, wavelengths = read_matrix_spectrum(filename)
    timesteps = data['timestep'].values
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 7))
    plt.subplots_adjust(bottom=0.25)  # 为滑块留出空间
    
    # 设置白色背景
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # 初始化：显示第一个时间步的累加光谱（就是第一步本身）
    initial_idx = 0
    # 累加：从第0步到initial_idx的总和
    counts = data.iloc[:initial_idx+1, 1:].sum(axis=0).values
    
    # 应用平滑滤波（window_length=21 提供更强的平滑效果）
    counts_filtered = savgol_filter(counts, window_length=11, polyorder=3)
    
    # 创建平滑曲线
    wavelengths_smooth = np.linspace(wavelengths.min(), wavelengths.max(), 
                                     len(wavelengths) * 10)
    spline = make_interp_spline(wavelengths, counts_filtered, k=3)
    counts_smooth = spline(wavelengths_smooth)
    
    # 绘制初始累加光谱线
    line, = ax.plot(wavelengths_smooth, counts_smooth, '-', 
                    color='teal', linewidth=1.5, 
                    label=f'Cumulative up to step: {timesteps[initial_idx]:.0f}')
    
    # 标记Hα中心线
    vline = ax.axvline(656.28, color='gray', linestyle='--', alpha=0.5, 
                       linewidth=1, label='Hα center (656.28 nm)')
    
    # 设置轴标签
    ax.set_xlabel('Wavelength [nm]', fontsize=12)
    ax.set_ylabel('Cumulative Intensity [arb.u.]', fontsize=12)
    ax.set_title('Hα Cumulative Spectrum Evolution - Interactive Slider', 
                 fontsize=14, fontweight='bold', pad=15)
    
    # 图例
    legend = ax.legend(fontsize=10, frameon=True, loc='upper right', 
                      framealpha=0.9, edgecolor='black')
    
    # 设置边框
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1)
    
    # 刻度样式
    ax.tick_params(axis='both', which='major', labelsize=10, 
                   direction='in', length=4, width=1)
    
    # 设置Y轴范围（固定，使用完全累加的最大值）
    final_cumulative = data.iloc[:, 1:].sum(axis=0).values
    final_filtered = savgol_filter(final_cumulative, window_length=21, polyorder=3)
    max_intensity = final_filtered.max()
    
    ax.set_ylim(0, max_intensity * 1.1)
    
    # 创建滑块轴
    ax_slider = plt.axes([0.15, 0.10, 0.7, 0.03], facecolor='lightgray')
    
    # 创建滑块
    slider = Slider(
        ax=ax_slider,
        label='Cumulative up to Time Step Index',
        valmin=0,
        valmax=len(timesteps) - 1,
        valinit=initial_idx,
        valstep=1,
        color='teal'
    )
    
    # 滑块更新函数
    def update(val):
        idx = int(slider.val)
        
        # 获取从第0步到第idx步的累加数据
        counts = data.iloc[:idx+1, 1:].sum(axis=0).values
        
        # 应用平滑（window_length=21 提供更强的平滑效果）
        counts_filtered = savgol_filter(counts, window_length=51, polyorder=3)
        
        # 样条插值
        spline = make_interp_spline(wavelengths, counts_filtered, k=3)
        counts_smooth = spline(wavelengths_smooth)
        
        # 更新曲线
        line.set_ydata(counts_smooth)
        
        # 更新图例标签
        line.set_label(f'Cumulative up to step: {timesteps[idx]:.0f}')
        legend = ax.legend(fontsize=10, frameon=True, loc='upper right', 
                          framealpha=0.9, edgecolor='black')
        
        # 重绘
        fig.canvas.draw_idle()
    
    # 连接滑块事件
    slider.on_changed(update)
    
    # 添加文本说明
    fig.text(0.5, 0.02, 'Drag slider to see cumulative spectrum growth over time', 
             ha='center', fontsize=10, style='italic', color='gray')
    
    print("\n✓ 交互式累加光谱图已创建")
    print("  提示: 拖动滑块查看从开始到某时间步的累加光谱")
    print("  关闭窗口以退出...")
    
    plt.show()
    
    return fig, ax, slider
    
    plt.show()
    
    return fig, ax, slider


if __name__ == '__main__':
    plot_interactive_spectrum_slider('results/spectrum_halpha_cumulative.dat')

