#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
交互式光谱滑块查看器 - 支持按反应通道分析
Interactive Spectrum Slider Viewer with Per-Reaction Channel Support
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import os
from scipy.interpolate import make_interp_spline, interp1d
from scipy.signal import savgol_filter
from matplotlib.widgets import Slider, Button, CheckButtons
import matplotlib.colors as mcolors

# 设置UTF-8编码以支持中文输出
if sys.platform == 'win32':
    os.environ['PYTHONIOENCODING'] = 'utf-8'
    sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout, 'reconfigure') else None


def read_per_reaction_spectrum(filename='results/spectrum_halpha_cumulative.dat'):
    """
    读取按反应通道分离的光谱数据
    
    新格式: timestep, reaction_id, λ1, λ2, λ3, ...
    """
    # 读取文件头获取波长信息
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 跳过注释行找到波长头
    header_line = None
    for line in lines:
        if line.startswith('#'):
            continue
        header_line = line.strip()
        break
    
    if header_line is None:
        print("❌ 无法找到头部信息")
        return None, None, None
    
    # 解析头部
    parts = header_line.split()
    if parts[0] == 'timestep' and parts[1] == 'reaction_id':
        # 新格式: 有reaction_id列
        wavelengths = np.array([float(x) for x in parts[2:]])
        
        # 读取数据
        df = pd.read_csv(filename, sep=r'\s+', comment='#', skiprows=0)
        df.columns = ['timestep', 'reaction_id'] + [f'wl_{i}' for i in range(len(wavelengths))]
        
        print(f"\u68c0\u6d4b\u5230\u6309\u53cd\u5e94\u901a\u9053\u683c\u5f0f")
        print(f"\u65f6\u95f4\u6b65\u6570: {df['timestep'].nunique()}")
        print(f"\u53cd\u5e94\u901a\u9053\u6570: {df['reaction_id'].nunique()}")
        print(f"\u6ce2\u957f\u70b9\u6570: {len(wavelengths)}")
        print(f"\u53cd\u5e94ID\u5217\u8868: {sorted(df['reaction_id'].unique())}")
        
        # 保存波长列名列表
        wl_cols = [f'wl_{i}' for i in range(len(wavelengths))]
        return df, wavelengths, sorted(df['reaction_id'].unique()), wl_cols
    
    elif parts[0] == 'timestep':
        # 旧格式: 无reaction_id列
        wavelengths = np.array([float(x) for x in parts[1:]])
        df = pd.read_csv(filename, sep=r'\s+', comment='#')
        df.columns = ['timestep'] + [f'wl_{i}' for i in range(len(wavelengths))]
        
        print(f"\u68c0\u6d4b\u5230\u65e7\u7248\u683c\u5f0f (\u65e0\u53cd\u5e94\u5206\u79bb)")
        print(f"\u65f6\u95f4\u6b65\u6570: {df['timestep'].nunique()}")
        print(f"\u6ce2\u957f\u70b9\u6570: {len(wavelengths)}")
        
        # \u6ce2\u957f\u5217\u540d
        wl_cols = [f'wl_{i}' for i in range(len(wavelengths))]
        
        # \u6dfb\u52a0\u865a\u62df reaction_id \u5217\u4ee5\u517c\u5bb9
        df['reaction_id'] = 0
        return df, wavelengths, [0], wl_cols
    
    else:
        print(f"\u274c \u672a\u77e5\u7684\u6587\u4ef6\u683c\u5f0f")
        return None, None, None, None


def apply_symmetrization(wavelengths, intensities, center_wavelength=656.28):
    """对光谱进行对称化处理"""
    interp_func = interp1d(wavelengths, intensities, 
                           kind='linear', 
                           bounds_error=False, 
                           fill_value=0.0)
    
    symmetrized = np.zeros_like(intensities)
    for i, wl in enumerate(wavelengths):
        delta = wl - center_wavelength
        mirror_wl = center_wavelength - delta
        symmetrized[i] = intensities[i] + interp_func(mirror_wl)
    
    return symmetrized


def plot_per_reaction_spectrum(filename='results/spectrum_halpha_cumulative.dat'):
    """
    创建交互式光谱图，支持按反应通道显示
    - 时间滑块选择累加范围
    - 每个反应通道有独立复选框控制显示
    - 自动叠加所有选中的反应通道
    """
    print("\n" + "=" * 70)
    print("绘制按反应通道交互式光谱图...")
    print("=" * 70)
    
    # \u8bfb\u53d6\u6570\u636e
    df, wavelengths, reaction_ids, wl_cols = read_per_reaction_spectrum(filename)
    if df is None:
        print(f"\u274c \u9519\u8bef: \u65e0\u6cd5\u8bfb\u53d6\u6587\u4ef6 '{filename}'")
        return
    
    timesteps = sorted(df['timestep'].unique())
    n_reactions = len(reaction_ids)
    
    # 生成颜色
    colors = list(mcolors.TABLEAU_COLORS.values())
    while len(colors) < n_reactions:
        colors.extend(colors)
    
    # 创建图形
    fig, ax = plt.subplots(figsize=(14, 8))
    plt.subplots_adjust(bottom=0.20, right=0.75)
    
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # 初始化选中状态 (所有反应默认选中)
    selected = {rid: True for rid in reaction_ids}
    
    # 创建平滑波长网格
    wavelengths_smooth = np.linspace(wavelengths.min(), wavelengths.max(), len(wavelengths) * 5)
    
    # 每个反应创建一条线
    lines = {}
    for i, rid in enumerate(reaction_ids):
        line, = ax.plot(wavelengths_smooth, np.zeros_like(wavelengths_smooth), 
                       '-', color=colors[i], linewidth=1.5, alpha=0.7,
                       label=f'Reaction {rid}')
        lines[rid] = line
    
    # 创建总和线
    line_total, = ax.plot(wavelengths_smooth, np.zeros_like(wavelengths_smooth), 
                         '-', color='black', linewidth=2.5, 
                         label='Total (Sum)')
    
    # Hα中心线
    ax.axvline(656.28, color='gray', linestyle='--', alpha=0.5, linewidth=1)
    
    ax.set_xlabel('Wavelength [nm]', fontsize=12)
    ax.set_ylabel('Intensity [counts]', fontsize=12)
    ax.set_title('Hα Spectrum by Reaction Channel', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.2)
    ax.legend(loc='upper left', fontsize=8)
    
    # 计算全局Y轴范围
    all_counts = df.iloc[:, 2:-1].sum(axis=0).values if 'reaction_id' in df.columns else df.iloc[:, 1:].sum(axis=0).values
    global_max = all_counts.max() * 1.2 if all_counts.max() > 0 else 100
    ax.set_ylim(-0.02 * global_max, global_max)
    
    # 时间滑块
    ax_slider = plt.axes([0.15, 0.10, 0.55, 0.025], facecolor='lightgray')
    slider = Slider(ax_slider, 'Time Step', 0, max(1, len(timesteps)-1), 
                   valinit=max(0, len(timesteps)-1), valstep=1, color='teal')
    
    # 平滑窗口滑块
    ax_smooth = plt.axes([0.15, 0.06, 0.35, 0.025], facecolor='lightgray')
    smooth_slider = Slider(ax_smooth, 'Smooth Window', 3, 51, 
                          valinit=11, valstep=2, color='orange')
    
    # 多项式阶数滑块
    ax_poly = plt.axes([0.55, 0.06, 0.15, 0.025], facecolor='lightgray')
    poly_slider = Slider(ax_poly, 'Poly Order', 1, 5, 
                        valinit=3, valstep=1, color='green')
    
    # 反应选择复选框
    ax_check = plt.axes([0.77, 0.30, 0.20, 0.55], facecolor='white')
    ax_check.set_title('Reactions', fontsize=10)
    labels = [f'R{rid}' for rid in reaction_ids] + ['Total']
    actives = [True] * (n_reactions + 1)  # 默认全部选中
    check = CheckButtons(ax_check, labels, actives)
    
    # 设置复选框颜色
    for i, rid in enumerate(reaction_ids):
        check.labels[i].set_color(colors[i])
    check.labels[-1].set_color('black')
    check.labels[-1].set_fontweight('bold')
    
    # 平滑参数存储
    smooth_params = {'window': 11, 'polyorder': 3}
    
    # 对称化和Auto Scale按钮
    ax_btn_sym = plt.axes([0.77, 0.20, 0.10, 0.04])
    btn_sym = Button(ax_btn_sym, 'Symmetrize')
    symmetrize_enabled = [False]
    
    ax_btn_scale = plt.axes([0.77, 0.14, 0.10, 0.04])
    btn_scale = Button(ax_btn_scale, 'Auto Scale')
    
    def update(val=None):
        """更新所有曲线"""
        idx = int(slider.val)
        end_timestep = timesteps[idx]
        
        # 获取平滑参数
        window = int(smooth_slider.val)
        polyorder = int(poly_slider.val)
        # 确保window是奇数且大于polyorder
        if window % 2 == 0:
            window += 1
        if window <= polyorder:
            window = polyorder + 2
        
        # 筛选时间步 <= end_timestep 的数据
        df_subset = df[df['timestep'] <= end_timestep]
        
        total_counts = np.zeros(len(wavelengths))
        
        for i, rid in enumerate(reaction_ids):
            # 获取该反应在所有选定时间步的累加光谱
            df_r = df_subset[df_subset['reaction_id'] == rid]
            if len(df_r) == 0:
                counts = np.zeros(len(wavelengths))
            else:
                counts = df_r[wl_cols].sum(axis=0).values
            
            # 对称化（如果启用）
            if symmetrize_enabled[0]:
                counts = apply_symmetrization(wavelengths, counts)
            
            # 平滑 (使用滑块参数)
            try:
                actual_window = min(window, len(counts)//2*2-1)
                if actual_window > polyorder:
                    counts_filtered = savgol_filter(counts, window_length=actual_window, polyorder=polyorder)
                else:
                    counts_filtered = counts
                spline = make_interp_spline(wavelengths, counts_filtered, k=3)
                counts_smooth = np.maximum(spline(wavelengths_smooth), 0)
            except:
                counts_smooth = np.interp(wavelengths_smooth, wavelengths, counts)
            
            # 更新曲线
            lines[rid].set_ydata(counts_smooth)
            lines[rid].set_visible(selected.get(rid, True))
            
            # 累加到总和
            if selected.get(rid, True):
                total_counts += counts
        
        # 更新总和曲线
        if symmetrize_enabled[0]:
            total_counts = apply_symmetrization(wavelengths, total_counts) / 2
        
        try:
            actual_window = min(window, len(total_counts)//2*2-1)
            if actual_window > polyorder:
                total_filtered = savgol_filter(total_counts, window_length=actual_window, polyorder=polyorder)
            else:
                total_filtered = total_counts
            spline = make_interp_spline(wavelengths, total_filtered, k=3)
            total_smooth = np.maximum(spline(wavelengths_smooth), 0)
        except:
            total_smooth = np.interp(wavelengths_smooth, wavelengths, total_counts)
        
        line_total.set_ydata(total_smooth)
        line_total.set_visible(selected.get('total', True))
        
        ax.set_title(f'Hα Spectrum (up to step {end_timestep})', fontsize=14, fontweight='bold')
        fig.canvas.draw_idle()
    
    def toggle_reaction(label):
        """切换反应通道显示"""
        if label == 'Total':
            selected['total'] = not selected.get('total', True)
        else:
            rid = int(label[1:])  # 从 "R1" 提取 1
            selected[rid] = not selected.get(rid, True)
        update()
    
    def toggle_symmetrize(event):
        symmetrize_enabled[0] = not symmetrize_enabled[0]
        btn_sym.label.set_text('Sym: ON' if symmetrize_enabled[0] else 'Symmetrize')
        update()
    
    def autoscale(event):
        cur_max = max(line_total.get_ydata().max(), 
                     max(lines[rid].get_ydata().max() for rid in reaction_ids if selected.get(rid, True)))
        ax.set_ylim(-0.02 * cur_max, cur_max * 1.1)
        fig.canvas.draw_idle()
    
    slider.on_changed(update)
    smooth_slider.on_changed(update)
    poly_slider.on_changed(update)
    check.on_clicked(toggle_reaction)
    btn_sym.on_clicked(toggle_symmetrize)
    btn_scale.on_clicked(autoscale)
    
    # 初始更新
    update()
    
    print("\n✓ 交互式按反应通道光谱图已创建")
    print("  滑块说明：")
    print("    Time Step      - 选择累积时间范围")
    print("    Smooth Window  - Savgol滤波窗口大小 (3-51)")
    print("    Poly Order     - 多项式阶数 (1-5)")
    print("  右侧复选框：勾选/取消对应反应通道")
    plt.show()


if __name__ == '__main__':
    import argparse
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    parser = argparse.ArgumentParser(description='Interactive Per-Reaction Spectrum Viewer')
    parser.add_argument('file', nargs='?', default=None, help='Path to spectrum data file')
    args = parser.parse_args()
    
    candidates = [
        'spectrum_halpha_cumulative.dat',
    ]
    
    target_file = None
    
    if args.file:
        if os.path.exists(args.file):
            target_file = args.file
        else:
            print(f"⚠️ 指定的文件不存在: {args.file}")
    
    if target_file is None:
        search_dirs = [
            os.path.join(os.getcwd(), 'results'),
            os.path.join(script_dir, 'results')
        ]
        
        for d in search_dirs:
            for fname in candidates:
                fpath = os.path.join(d, fname)
                if os.path.exists(fpath):
                    target_file = fpath
                    print(f"自动找到文件: {target_file}")
                    break
            if target_file: 
                break
            
    if target_file:
        plot_per_reaction_spectrum(target_file)
    else:
        print("❌ 未找到光谱数据文件")
        print("请手动指定: python plot_interactive_slider.py path/to/file.dat")
