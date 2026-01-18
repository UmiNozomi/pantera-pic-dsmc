#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot Hα Spectral Diagnostics from PANTERA simulation
读取并可视化PANTERA模拟的Hα光谱诊断结果
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D

# 设置中文字体支持（可选）
plt.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

def read_spectrum_data(filename):
    """
    读取累积光谱文件
    
    Parameters:
    -----------
    filename : str
        光谱文件路径，例如 'results/spectrum_halpha_cumulative.dat'
    
    Returns:
    --------
    data : pandas.DataFrame
        包含 timestep, vLOS, wavelength, counts 的数据框
    """
    # 读取数据，跳过注释行
    data = pd.read_csv(filename, sep=r'\s+', comment='#', 
                       names=['timestep', 'vLOS', 'wavelength', 'counts'])
    
    print(f"读取了 {len(data)} 行数据")
    print(f"时间步范围: {data['timestep'].min()} - {data['timestep'].max()}")
    print(f"波长范围: {data['wavelength'].min():.2f} - {data['wavelength'].max():.2f} nm")
    
    return data


def plot_single_spectrum(data, timestep, save_fig=True):
    """
    绘制单个时间步的光谱
    
    Parameters:
    -----------
    data : pandas.DataFrame
        光谱数据
    timestep : int
        要绘制的时间步
    save_fig : bool
        是否保存图片
    """
    # 提取指定时间步的数据
    df = data[data['timestep'] == timestep]
    
    if len(df) == 0:
        print(f"警告: 时间步 {timestep} 没有数据!")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 绘制光谱
    ax.plot(df['wavelength'], df['counts'], 'b-', linewidth=1.5, label=f'Timestep {timestep}')
    ax.fill_between(df['wavelength'], df['counts'], alpha=0.3)
    
    # 标记Hα中心线
    ax.axvline(656.28, color='r', linestyle='--', alpha=0.5, label='Hα center (656.28 nm)')
    
    ax.set_xlabel('Wavelength [nm]', fontsize=12)
    ax.set_ylabel('Counts', fontsize=12)
    ax.set_title(f'Hα Spectrum - Timestep {timestep}', fontsize=14, fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_fig:
        output_file = f'halpha_spectrum_t{timestep:08d}.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"保存图片: {output_file}")
    
    plt.show()


def plot_time_evolution_heatmap(data, save_fig=True):
    """
    绘制光谱随时间演化的热图
    
    Parameters:
    -----------
    data : pandas.DataFrame
        光谱数据
    save_fig : bool
        是否保存图片
    """
    # 获取所有唯一的时间步
    timesteps = sorted(data['timestep'].unique())
    wavelengths = sorted(data['wavelength'].unique())
    
    # 创建2D数组
    n_time = len(timesteps)
    n_wave = len(wavelengths)
    intensity = np.zeros((n_time, n_wave))
    
    for i, t in enumerate(timesteps):
        df_t = data[data['timestep'] == t]
        intensity[i, :] = df_t['counts'].values
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 绘制热图
    im = ax.pcolormesh(wavelengths, timesteps, intensity, 
                       shading='auto', cmap='hot')
    
    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Counts', fontsize=12)
    
    # 标记Hα中心线
    ax.axvline(656.28, color='cyan', linestyle='--', alpha=0.7, 
               linewidth=2, label='Hα center')
    
    ax.set_xlabel('Wavelength [nm]', fontsize=12)
    ax.set_ylabel('Timestep', fontsize=12)
    ax.set_title('Hα Spectrum Evolution (Heatmap)', fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    if save_fig:
        output_file = 'halpha_evolution_heatmap.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"保存图片: {output_file}")
    
    plt.show()


def plot_time_evolution_3d(data, save_fig=True):
    """
    绘制光谱随时间演化的3D图
    
    Parameters:
    -----------
    data : pandas.DataFrame
        光谱数据
    save_fig : bool
        是否保存图片
    """
    # 获取所有唯一的时间步
    timesteps = sorted(data['timestep'].unique())
    wavelengths = sorted(data['wavelength'].unique())
    
    # 创建网格
    n_time = len(timesteps)
    n_wave = len(wavelengths)
    
    T, W = np.meshgrid(timesteps, wavelengths)
    I = np.zeros((n_wave, n_time))
    
    for i, t in enumerate(timesteps):
        df_t = data[data['timestep'] == t]
        I[:, i] = df_t['counts'].values
    
    fig = plt.figure(figsize=(14, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制3D表面
    surf = ax.plot_surface(W, T, I, cmap='hot', 
                           linewidth=0, antialiased=True, alpha=0.9)
    
    # 添加颜色条
    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)
    cbar.set_label('Counts', fontsize=12)
    
    ax.set_xlabel('Wavelength [nm]', fontsize=12)
    ax.set_ylabel('Timestep', fontsize=12)
    ax.set_zlabel('Counts', fontsize=12)
    ax.set_title('Hα Spectrum Evolution (3D)', fontsize=14, fontweight='bold')
    
    # 设置视角
    ax.view_init(elev=25, azim=45)
    
    plt.tight_layout()
    
    if save_fig:
        output_file = 'halpha_evolution_3d.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"保存图片: {output_file}")
    
    plt.show()


def plot_multiple_timesteps(data, timesteps, save_fig=True):
    """
    在一个图中比较多个时间步的光谱
    
    Parameters:
    -----------
    data : pandas.DataFrame
        光谱数据
    timesteps : list
        要比较的时间步列表
    save_fig : bool
        是否保存图片
    """
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # 使用不同颜色绘制各个时间步
    colors = cm.viridis(np.linspace(0, 1, len(timesteps)))
    
    for i, t in enumerate(timesteps):
        df_t = data[data['timestep'] == t]
        if len(df_t) > 0:
            ax.plot(df_t['wavelength'], df_t['counts'], 
                   color=colors[i], linewidth=1.5, label=f't={t}')
    
    # 标记Hα中心线
    ax.axvline(656.28, color='r', linestyle='--', alpha=0.5, label='Hα center')
    
    ax.set_xlabel('Wavelength [nm]', fontsize=12)
    ax.set_ylabel('Counts', fontsize=12)
    ax.set_title('Hα Spectrum Comparison', fontsize=14, fontweight='bold')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save_fig:
        output_file = 'halpha_comparison.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"保存图片: {output_file}")
    
    plt.show()


def main():
    """主函数 - 示例用法"""
    
    # ========== 配置 ==========
    # 修改为您的光谱文件路径
    spectrum_file = 'results/spectrum_halpha_cumulative.dat'
    
    # 读取数据
    print("=" * 60)
    print("读取光谱数据...")
    print("=" * 60)
    data = read_spectrum_data(spectrum_file)
    
    # 获取所有可用的时间步
    available_timesteps = sorted(data['timestep'].unique())
    print(f"\n可用的时间步: {available_timesteps}")
    
    # ========== 绘图示例 ==========
    
    # 1. 绘制单个时间步的光谱
    if len(available_timesteps) > 0:
        print("\n正在绘制单个时间步光谱...")
        plot_single_spectrum(data, available_timesteps[-1])  # 绘制最后一个时间步
    
    # 2. 绘制热图
    if len(available_timesteps) > 1:
        print("\n正在绘制时间演化热图...")
        plot_time_evolution_heatmap(data)
    
    # 3. 绘制3D图
    if len(available_timesteps) > 1:
        print("\n正在绘制3D演化图...")
        plot_time_evolution_3d(data)
    
    # 4. 比较多个时间步
    if len(available_timesteps) >= 3:
        print("\n正在绘制多时间步比较图...")
        # 选择几个有代表性的时间步进行比较
        timesteps_to_compare = [available_timesteps[0], 
                                available_timesteps[len(available_timesteps)//2],
                                available_timesteps[-1]]
        plot_multiple_timesteps(data, timesteps_to_compare)
    
    print("\n" + "=" * 60)
    print("绘图完成!")
    print("=" * 60)


if __name__ == '__main__':
    main()
