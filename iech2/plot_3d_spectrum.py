#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
简化版3D光谱绘制脚本 - 支持矩阵格式
Simplified 3D Spectrum Plotter - Matrix Format Support
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

def read_matrix_spectrum(filename='results/spectrum_halpha_cumulative.dat'):
    """
    读取矩阵格式的光谱数据
    
    参数:
    -----
    filename : str
        光谱数据文件路径
        
    返回:
    -----
    data : pandas.DataFrame
        数据框，第一列是timestep，其余列是各波长的计数
    wavelengths : numpy.array
        波长数组 [nm]
    """
    # 读取数据，跳过注释行
    data = pd.read_csv(filename, sep=r'\s+', comment='#')
    
    # 第一列是timestep，其余列名是波长
    wavelengths = data.columns[1:].astype(float).values
    
    print(f"读取了 {len(data)} 个时间步")
    print(f"波长点数: {len(wavelengths)}")
    print(f"时间步范围: {data['timestep'].min()} - {data['timestep'].max()}")
    print(f"波长范围: {wavelengths[0]:.4f} - {wavelengths[-1]:.4f} nm")
    
    return data, wavelengths


def plot_cumulative_spectrum(filename='results/spectrum_halpha_cumulative.dat', 
                             save=True):
    """
    绘制所有时间步累加的总光谱（时间积分光谱）
    
    参数:
    -----
    filename : str
        光谱数据文件路径
    save : bool
        是否保存图片
    """
    # 读取数据
    print("\n" + "="*70)
    print("绘制时间积分光谱...")
    print("="*70)
    
    data, wavelengths = read_matrix_spectrum(filename)
    
    # 计算所有时间步的累加和（对每列求和，跳过timestep列）
    total_counts = data.iloc[:, 1:].sum(axis=0).values
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(12, 7))
    
    # 绘制累加光谱
    ax.plot(wavelengths, total_counts, 'b-', linewidth=2, label='时间积分光谱')
    ax.fill_between(wavelengths, total_counts, alpha=0.3)
    
    # 标记Hα中心线
    ax.axvline(656.28, color='r', linestyle='--', alpha=0.7, 
               linewidth=2, label='Hα中心 (656.28 nm)')
    
    # 设置标签
    ax.set_xlabel('波长 (Wavelength) [nm]', fontsize=13)
    ax.set_ylabel('累积发光强度 (Total Counts)', fontsize=13)
    ax.set_title('Hα时间积分光谱\nTime-Integrated Hα Spectrum', 
                 fontsize=15, fontweight='bold')
    
    # 添加统计信息
    total_events = total_counts.sum()
    peak_wavelength = wavelengths[np.argmax(total_counts)]
    ax.text(0.02, 0.98, 
            f'总事件数: {int(total_events):,}\n峰值波长: {peak_wavelength:.4f} nm',
            transform=ax.transAxes, fontsize=11,
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    if save:
        output_file = 'halpha_cumulative_spectrum.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ 累加光谱图已保存: {output_file}")
    
    plt.show()


def plot_3d_spectrum(filename='results/spectrum_halpha_cumulative.dat', 
                     elevation=25, azimuth=45, save=True):
    """
    绘制Hα光谱的3D演化图
    
    参数:
    -----
    filename : str
        光谱数据文件路径
    elevation : float
        视角仰角 (0-90度)，推荐 20-30
    azimuth : float  
        视角方位角 (0-360度)，推荐 30-60
    save : bool
        是否保存图片
        
    交互操作:
    --------
    - 鼠标左键拖动: 旋转视角
    - 鼠标滚轮: 缩放
    - 右键拖动: 平移
    """
    
    # 读取数据
    print("\n" + "="*70)
    print("绘制3D演化图...")
    print("="*70)
    
    data, wavelengths = read_matrix_spectrum(filename)
    
    # 获取时间步
    timesteps = data['timestep'].values
    
    # 获取强度矩阵 (转置，使得行=波长，列=时间)
    intensity = data.iloc[:, 1:].T.values
    
    # 创建网格
    T, W = np.meshgrid(timesteps, wavelengths)
    
    # 创建3D图
    fig = plt.figure(figsize=(14, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制3D表面 - 使用'hot'配色（黑-红-黄-白）
    surf = ax.plot_surface(W, T, intensity, cmap='hot', 
                           linewidth=0, antialiased=True, 
                           alpha=0.9, edgecolor='none')
    
    # 添加颜色条
    cbar = fig.colorbar(surf, ax=ax, shrink=0.6, aspect=8, pad=0.1)
    cbar.set_label('发光强度 (Counts)', fontsize=13, labelpad=15)
    cbar.ax.tick_params(labelsize=11)
    
    # 设置坐标轴标签
    ax.set_xlabel('波长 (Wavelength) [nm]', fontsize=13, labelpad=10)
    ax.set_ylabel('时间步 (Timestep)', fontsize=13, labelpad=10)
    ax.set_zlabel('发光强度 (Counts)', fontsize=13, labelpad=10)
    
    # 设置标题
    ax.set_title('Hα光谱时间演化 (3D)\nHα Spectrum Time Evolution', 
                 fontsize=15, fontweight='bold', pad=20)
    
    # 设置视角
    ax.view_init(elev=elevation, azim=azimuth)
    
    # 优化网格线
    ax.grid(True, alpha=0.3)
    
    # 设置刻度字体大小
    ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    
    if save:
        output_file = 'halpha_3d_evolution.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ 3D图已保存: {output_file}")
    
    print("\n提示: 窗口打开后可以用鼠标拖动旋转视角")
    print("      鼠标滚轮可以缩放")
    
    plt.show()


def plot_time_evolution_heatmap(filename='results/spectrum_halpha_cumulative.dat',
                                save=True):
    """
    绘制光谱随时间演化的热图
    
    参数:
    -----
    filename : str
        光谱数据文件路径
    save : bool
        是否保存图片
    """
    # 读取数据
    print("\n" + "="*70)
    print("绘制热图...")
    print("="*70)
    
    data, wavelengths = read_matrix_spectrum(filename)
    
    # 获取时间步
    timesteps = data['timestep'].values
    
    # 获取强度矩阵
    intensity = data.iloc[:, 1:].values
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(14, 7))
    
    # 绘制热图
    im = ax.pcolormesh(wavelengths, timesteps, intensity, 
                       shading='auto', cmap='hot')
    
    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('发光强度 (Counts)', fontsize=12)
    
    # 标记Hα中心线
    ax.axvline(656.28, color='cyan', linestyle='--', alpha=0.7, 
               linewidth=2, label='Hα中心')
    
    ax.set_xlabel('波长 (Wavelength) [nm]', fontsize=12)
    ax.set_ylabel('时间步 (Timestep)', fontsize=12)
    ax.set_title('Hα光谱时间演化 (热图)\nHα Spectrum Time Evolution (Heatmap)', 
                 fontsize=14, fontweight='bold')
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    
    if save:
        output_file = 'halpha_evolution_heatmap.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ 热图已保存: {output_file}")
    
    plt.show()


if __name__ == '__main__':
    # ========== 使用示例 ==========
    
    print("=" * 70)
    print("Hα光谱可视化工具")
    print("=" * 70)
    
    filename = 'results/spectrum_halpha_cumulative.dat'
    
    # 1. 绘制时间积分累加光谱（所有时间步的总和）
    print("\n[1] 绘制时间积分累加光谱...")
    plot_cumulative_spectrum(filename, save=True)
    
    # 2. 绘制3D演化图
    print("\n[2] 绘制3D演化图...")
    plot_3d_spectrum(
        filename=filename,
        elevation=25,    # 仰角 (推荐 20-30)
        azimuth=45,      # 方位角 (推荐 30-60)  
        save=True
    )
    
    # 3. 绘制热图（取消注释以使用）
    # print("\n[3] 绘制热图...")
    # plot_time_evolution_heatmap(filename, save=True)
    
    print("\n" + "=" * 70)
    print("绘图完成!")
    print("=" * 70)

    """
    绘制Hα光谱的3D演化图
    
    参数:
    -----
    filename : str
        光谱数据文件路径
    elevation : float
        视角仰角 (0-90度)，推荐 20-30
    azimuth : float  
        视角方位角 (0-360度)，推荐 30-60
    save : bool
        是否保存图片
        
    交互操作:
    --------
    - 鼠标左键拖动: 旋转视角
    - 鼠标滚轮: 缩放
    - 右键拖动: 平移
    """
    
    # 读取数据
    print(f"正在读取: {filename}")
    data = pd.read_csv(filename, sep=r'\s+', comment='#', 
                       names=['timestep', 'vLOS', 'wavelength', 'counts'])
    
    print(f"数据行数: {len(data)}")
    
    # 获取唯一值
    timesteps = sorted(data['timestep'].unique())
    wavelengths = sorted(data['wavelength'].unique())
    
    print(f"时间步数: {len(timesteps)}")
    print(f"波长点数: {len(wavelengths)}")
    print(f"时间步范围: {timesteps[0]} - {timesteps[-1]}")
    print(f"波长范围: {wavelengths[0]:.2f} - {wavelengths[-1]:.2f} nm")
    
    # 创建网格
    n_time = len(timesteps)
    n_wave = len(wavelengths)
    
    T, W = np.meshgrid(timesteps, wavelengths)
    I = np.zeros((n_wave, n_time))
    
    # 填充强度数据
    for i, t in enumerate(timesteps):
        df_t = data[data['timestep'] == t]
        I[:, i] = df_t['counts'].values
    
    # 创建3D图
    fig = plt.figure(figsize=(14, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    # 绘制3D表面 - 使用'hot'配色（黑-红-黄-白）
    surf = ax.plot_surface(W, T, I, cmap='hot', 
                           linewidth=0, antialiased=True, 
                           alpha=0.9, edgecolor='none')
    
    # 添加颜色条
    cbar = fig.colorbar(surf, ax=ax, shrink=0.6, aspect=8, pad=0.1)
    cbar.set_label('发光强度 (Counts)', fontsize=13, labelpad=15)
    cbar.ax.tick_params(labelsize=11)
    
    # 设置坐标轴标签
    ax.set_xlabel('波长 (Wavelength) [nm]', fontsize=13, labelpad=10)
    ax.set_ylabel('时间步 (Timestep)', fontsize=13, labelpad=10)
    ax.set_zlabel('发光强度 (Counts)', fontsize=13, labelpad=10)
    
    # 设置标题
    ax.set_title('Hα光谱时间演化 (3D)\nHα Spectrum Time Evolution', 
                 fontsize=15, fontweight='bold', pad=20)
    
    # 设置视角
    ax.view_init(elev=elevation, azim=azimuth)
    
    # 优化网格线
    ax.grid(True, alpha=0.3)
    
    # 设置刻度字体大小
    ax.tick_params(labelsize=10)
    
    plt.tight_layout()
    
    if save:
        output_file = 'halpha_3d_evolution.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"\n✓ 图片已保存: {output_file}")
    
    print("\n提示: 窗口打开后可以用鼠标拖动旋转视角")
    print("      鼠标滚轮可以缩放")
    
    plt.show()


def plot_3d_multiple_views(filename='results/spectrum_halpha_cumulative.dat'):
    """
    在一个窗口中显示多个视角的3D图
    """
    # 读取数据
    data = pd.read_csv(filename, sep=r'\s+', comment='#', 
                       names=['timestep', 'vLOS', 'wavelength', 'counts'])
    
    timesteps = sorted(data['timestep'].unique())
    wavelengths = sorted(data['wavelength'].unique())
    
    n_time = len(timesteps)
    n_wave = len(wavelengths)
    
    T, W = np.meshgrid(timesteps, wavelengths)
    I = np.zeros((n_wave, n_time))
    
    for i, t in enumerate(timesteps):
        df_t = data[data['timestep'] == t]
        I[:, i] = df_t['counts'].values
    
    # 创建2x2子图，显示4个不同视角
    fig = plt.figure(figsize=(16, 12))
    
    # 定义4个不同的视角
    views = [
        (25, 45, '默认视角'),
        (10, 120, '侧视图'),
        (60, 45, '俯视图'),
        (25, 225, '反向视角')
    ]
    
    for idx, (elev, azim, title) in enumerate(views, 1):
        ax = fig.add_subplot(2, 2, idx, projection='3d')
        
        surf = ax.plot_surface(W, T, I, cmap='hot', 
                               linewidth=0, antialiased=True, 
                               alpha=0.9)
        
        ax.set_xlabel('Wavelength [nm]', fontsize=10)
        ax.set_ylabel('Timestep', fontsize=10)
        ax.set_zlabel('Counts', fontsize=10)
        ax.set_title(title, fontsize=12, fontweight='bold')
        ax.view_init(elev=elev, azim=azim)
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('halpha_3d_multiple_views.png', dpi=300, bbox_inches='tight')
    print("✓ 多视角图片已保存: halpha_3d_multiple_views.png")
    plt.show()


if __name__ == '__main__':
    # ========== 使用示例 ==========
    
    print("=" * 70)
    print("Hα光谱3D可视化")
    print("=" * 70)
    
    # 方法1: 绘制标准3D图（可交互）
    print("\n[1] 绘制交互式3D图...")
    plot_3d_spectrum(
        filename='results/spectrum_halpha_cumulative.dat',
        elevation=25,    # 仰角 (推荐 20-30)
        azimuth=45,      # 方位角 (推荐 30-60)  
        save=True
    )
    
    # 方法2: 绘制多视角对比图（取消注释以使用）
    # print("\n[2] 绘制多视角对比图...")
    # plot_3d_multiple_views('results/spectrum_halpha_cumulative.dat')
    
    print("\n" + "=" * 70)
    print("绘图完成!")
    print("=" * 70)
