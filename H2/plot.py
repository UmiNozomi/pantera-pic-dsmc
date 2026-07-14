import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from tkinter import Tk, Label, Entry, Button, Frame, StringVar, ttk, LabelFrame, Toplevel, Checkbutton, BooleanVar, messagebox, filedialog
import tkinter as tk
from scipy import interpolate
from scipy.ndimage import gaussian_filter1d

# ── 配置区 ──
# 将 root_dir 设为脚本所在目录
root_dir = Path(__file__).parent   # 或者直接 Path('.')
file_stem = 'conservation_checks'   # 无后缀的文件"干名"

# 自动读取所有变量
def get_all_variables():
    """自动从CSV文件中读取所有变量名（除了time列）"""
    files = [f for f in root_dir.rglob('*') if f.is_file() and f.stem == file_stem]
    if not files:
        print("⚠️ 未找到任何匹配的文件")
        return []
    
    # 读取第一个文件来获取列名
    first_file = files[0]
    try:
        df = pd.read_csv(first_file, sep='\s+', header=0)
        # 排除time列，获取所有其他列名
        variables = [col for col in df.columns if col.lower() != 'time']
        print(f"📊 发现变量: {variables}")
        return variables
    except Exception as e:
        print(f"❌ 读取文件失败: {e}")
        return ['nreact_1', 'nreact_2', 'nreact_3']  # 默认值

variables = get_all_variables()
# ─────────────────

files = [f for f in root_dir.rglob('*') if f.is_file() and f.stem == file_stem]
data_list = [(f.parent.name, pd.read_csv(f, sep='\s+', header=0))
             for f in files]

# 变量选择界面
class VariableSelector:
    def __init__(self):
        self.root = Tk()
        self.root.title("📊 变量选择器")
        self.root.geometry("400x600")
        self.root.resizable(True, True)
        
        # 变量选择状态
        self.var_states = {}
        self.controllers = []  # 存储已打开的图表控制器
        
        self.setup_ui()
    
    def setup_ui(self):
        """设置用户界面"""
        # 主标题
        title_label = Label(self.root, text="🎯 选择要绘制的变量", 
                           font=('Arial', 16, 'bold'), pady=10)
        title_label.pack()
        
        # 创建中间容器框架，用于放置滚动区域
        middle_frame = Frame(self.root)
        middle_frame.pack(fill='both', expand=True, padx=10, pady=(0, 10))
        
        # 创建滚动框架
        self.canvas = tk.Canvas(middle_frame, highlightthickness=0)
        scrollbar = ttk.Scrollbar(middle_frame, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = ttk.Frame(self.canvas)
        
        # 配置滚动
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        )
        
        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=scrollbar.set)
        
        # 绑定鼠标滚轮事件
        self.bind_mousewheel()
        
        # 布局滚动区域
        self.canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # 变量选择区域
        var_frame = LabelFrame(self.scrollable_frame, text="📋 可用变量列表", font=('Arial', 12, 'bold'))
        var_frame.pack(pady=10, padx=10, fill='both', expand=True)
        
        # 添加全选/取消全选按钮
        button_frame = Frame(var_frame)
        button_frame.pack(fill='x', pady=5, padx=5)
        
        ttk.Button(button_frame, text="全选", command=self.select_all).pack(side='left', padx=5)
        ttk.Button(button_frame, text="取消全选", command=self.deselect_all).pack(side='left', padx=5)
        ttk.Button(button_frame, text="🔄 刷新变量", command=self.refresh_variables).pack(side='right', padx=5)
        
        # 变量复选框区域
        self.var_checkboxes_frame = Frame(var_frame)
        self.var_checkboxes_frame.pack(fill='both', expand=True, pady=5)
        
        self.create_variable_checkboxes()
        
        # 底部控制按钮区域
        control_frame = Frame(self.root)
        control_frame.pack(fill='x', pady=10, padx=10)
        
        # 选择统计信息
        self.info_label = Label(control_frame, text="已选择: 0 个变量", 
                               font=('Arial', 10))
        self.info_label.pack(pady=5)
        
        # 控制按钮
        button_control_frame = Frame(control_frame)
        button_control_frame.pack(fill='x', pady=5)
        
        ttk.Button(button_control_frame, text="📈 打开选中图表", 
                  command=self.open_selected_plots, 
                  style='Accent.TButton').pack(side='left', padx=5, fill='x', expand=True)
        
        ttk.Button(button_control_frame, text="📊 多变量合并图表", 
                  command=self.open_multi_variable_plot, 
                  style='Accent.TButton').pack(side='left', padx=5, fill='x', expand=True)
        
        ttk.Button(button_control_frame, text="❌ 关闭所有图表", 
                  command=self.close_all_plots).pack(side='left', padx=5, fill='x', expand=True)
        
        ttk.Button(button_control_frame, text="🚪 退出程序", 
                  command=self.exit_program).pack(side='left', padx=5, fill='x', expand=True)
        
        # 绑定选择变化事件
        self.update_info()
    
    def bind_mousewheel(self):
        """绑定鼠标滚轮事件"""
        def _on_mousewheel(event):
            self.canvas.yview_scroll(int(-1*(event.delta/120)), "units")
        
        def _bind_to_mousewheel(event):
            self.canvas.bind_all("<MouseWheel>", _on_mousewheel)
        
        def _unbind_from_mousewheel(event):
            self.canvas.unbind_all("<MouseWheel>")
        
        # 当鼠标进入canvas区域时绑定滚轮事件
        self.canvas.bind('<Enter>', _bind_to_mousewheel)
        # 当鼠标离开canvas区域时解绑滚轮事件
        self.canvas.bind('<Leave>', _unbind_from_mousewheel)
        
        # 直接绑定到滚动框架
        self.scrollable_frame.bind('<Enter>', _bind_to_mousewheel)
        self.scrollable_frame.bind('<Leave>', _unbind_from_mousewheel)
    
    def create_variable_checkboxes(self):
        """创建变量复选框"""
        # 清空现有复选框
        for widget in self.var_checkboxes_frame.winfo_children():
            widget.destroy()
        
        # 为每个变量创建复选框
        for var in variables:
            # 如果之前已经有状态，保持原状态，否则默认为False
            if var not in self.var_states:
                self.var_states[var] = BooleanVar(value=False)
            
            var_frame = Frame(self.var_checkboxes_frame)
            var_frame.pack(fill='x', pady=2, padx=10)
            
            cb = Checkbutton(var_frame, text=var, variable=self.var_states[var],
                            font=('Arial', 11), command=self.update_info)
            cb.pack(anchor='w')
        
        # 更新滚动区域
        self.root.update_idletasks()  # 确保所有组件都已渲染
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
    
    def select_all(self):
        """全选所有变量"""
        for var_state in self.var_states.values():
            var_state.set(True)
        self.update_info()
    
    def deselect_all(self):
        """取消选择所有变量"""
        for var_state in self.var_states.values():
            var_state.set(False)
        self.update_info()
    
    def refresh_variables(self):
        """刷新变量列表"""
        global variables
        variables = get_all_variables()
        
        # 重新加载数据
        global data_list
        files = [f for f in root_dir.rglob('*') if f.is_file() and f.stem == file_stem]
        data_list = [(f.parent.name, pd.read_csv(f, sep='\s+', header=0))
                     for f in files]
        
        # 重新创建复选框
        self.create_variable_checkboxes()
        self.update_info()
        
        # 更新滚动区域
        self.canvas.configure(scrollregion=self.canvas.bbox("all"))
        
        messagebox.showinfo("刷新完成", f"已刷新变量列表，发现 {len(variables)} 个变量")
    
    def update_info(self):
        """更新选择信息"""
        selected_count = sum(1 for var_state in self.var_states.values() if var_state.get())
        self.info_label.config(text=f"已选择: {selected_count} 个变量")
    
    def open_selected_plots(self):
        """打开选中的图表"""
        selected_vars = [var for var, var_state in self.var_states.items() if var_state.get()]
        
        if not selected_vars:
            messagebox.showwarning("未选择变量", "请至少选择一个变量来绘制图表")
            return
        
        # 关闭之前打开的图表
        self.close_all_plots()
        
        print(f"🚀 正在创建 {len(selected_vars)} 个图表...")
        
        # 为选中的变量创建图表控制器
        for var in selected_vars:
            controller = PlotController(var, master=self.root)
            self.controllers.append(controller)
        
        messagebox.showinfo("图表已打开", f"已成功打开 {len(selected_vars)} 个图表窗口")
    
    def open_multi_variable_plot(self):
        """打开多变量合并图表"""
        selected_vars = [var for var, var_state in self.var_states.items() if var_state.get()]
        
        if len(selected_vars) < 2:
            messagebox.showwarning("变量不足", "请至少选择2个变量来创建多变量合并图表")
            return
        
        print(f"🚀 正在创建多变量合并图表，包含 {len(selected_vars)} 个变量...")
        
        # 创建多变量图表控制器
        controller = MultiVariablePlotController(selected_vars, master=self.root)
        self.controllers.append(controller)
        
        messagebox.showinfo("多变量图表已打开", f"已成功创建包含 {len(selected_vars)} 个变量的合并图表")
    
    def close_all_plots(self):
        """关闭所有已打开的图表"""
        for controller in self.controllers:
            if not controller.is_closed:
                controller.on_closing()
        
        self.controllers.clear()
        
        # 清空全局控制器中的控制器列表
        global_controller.controllers.clear()
    
    def exit_program(self):
        """退出程序"""
        self.close_all_plots()
        self.root.quit()
        self.root.destroy()
    
    def run(self):
        """运行选择器"""
        self.root.mainloop()

# 全局控制器，用于协调所有图表
class GlobalController:
    def __init__(self):
        self.controllers = []
    
    def register_controller(self, controller):
        """注册一个图表控制器"""
        self.controllers.append(controller)
    
    def apply_to_all(self, action, *args, **kwargs):
        """对所有控制器执行指定操作"""
        for controller in self.controllers:
            if not controller.is_closed:
                action(controller, *args, **kwargs)

# 创建全局控制器实例
global_controller = GlobalController()

class PlotController:
    def __init__(self, var, master=None):
        # 创建主窗口
        self.root = Tk() if master is None else Toplevel(master)
        self.root.title(f"图表控制 - {var}")
        self.root.geometry("1200x950")  # 增加高度以适应新控件
        
        # 设置窗口关闭处理
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # 创建左右分栏
        self.left_frame = Frame(self.root)
        self.left_frame.pack(side='left', fill='both', expand=True, padx=10, pady=10)
        
        self.right_frame = Frame(self.root)
        self.right_frame.pack(side='right', fill='y', padx=10, pady=10)
        
        # 创建滚动条
        canvas = tk.Canvas(self.right_frame)
        scrollbar = ttk.Scrollbar(self.right_frame, orient="vertical", command=canvas.yview)
        self.scrollable_frame = ttk.Frame(canvas)
        
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # 创建图表
        self.fig = plt.figure(figsize=(10, 6))
        self.ax = self.fig.add_subplot(111)
        
        # 将matplotlib图表嵌入到tkinter窗口中
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.left_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        
        # 存储数据和图例信息
        self.data_series = []
        self.legend_vars = {}
        self.scale_vars = {}  # 存储倍率变量
        self.original_data = {}  # 存储原始数据
        self.visibility_vars = {}  # 存储显示/隐藏状态
        self.data_order = []  # 存储数据顺序
        self.series_colors = {}  # 存储每个数据系列的颜色
        self.point_size_var = StringVar(value="10")  # 数据点大小
        self.marker_style_var = StringVar(value="o")  # 数据点样式
        self.legend_size_var = StringVar(value="10")  # 图例标签文字大小
        self.legend_marker_scale_var = StringVar(value="1.0")  # 图例点缩放
        self.title_size_var = StringVar(value="12")  # 标题字体大小
        self.xlabel_size_var = StringVar(value="10")  # X轴标签字体大小
        self.ylabel_size_var = StringVar(value="10")  # Y轴标签字体大小
        self.label_position_var = StringVar(value="best")  # 图例位置
        self.x_unit_scale_var = StringVar(value="1")  # X轴单位缩放因子
        self.y_unit_scale_var = StringVar(value="1")  # Y轴单位缩放因子
        
        # 预定义颜色列表
        self.color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                          '#a6cee3', '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99']
        
        # 绘制数据
        color_index = 0
        for folder_name, df in data_list:
            if var not in df.columns:
                print(f"⚠️ {folder_name} 缺 {var}，跳过")
                continue
            # 存储原始数据
            self.original_data[folder_name] = (df['time'].values, df[var].values)
            
            # 为每个数据系列分配固定颜色
            color = self.color_list[color_index % len(self.color_list)]
            self.series_colors[folder_name] = color
            color_index += 1
            
            line = self.ax.scatter(df['time'], df[var], label=folder_name, s=10, alpha=0.7, color=color)
            self.data_series.append((line, folder_name, df))
            self.data_order.append(folder_name)  # 记录初始顺序
        
        self.ax.set_xlabel('time')
        self.ax.set_ylabel(var)
        self.ax.set_title(f'{var} vs time')
        # 只有在有带标签的图形元素时才创建图例
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:  # 如果有图例元素
            self.ax.legend()
        
        # 设置字体 [[memory:6787997]]
        plt.rcParams['font.family'] = ['Times New Roman', 'SimHei']  # 英文用Times New Roman，中文用SimHei
        
        # 设置样式
        self.apply_style()
        
        # 创建控制面板
        self.create_control_panel(var)
        
        # 存储变量名，用于保存文件
        self.var = var
        
        # 标记窗口状态
        self.is_closed = False
        
        # 注册到全局控制器
        global_controller.register_controller(self)
        
        # 绑定鼠标滚轮事件
        self.bind_mousewheel()

    def bind_mousewheel(self):
        """绑定鼠标滚轮事件"""
        def _on_mousewheel(event):
            self.scrollable_frame.master.yview_scroll(int(-1*(event.delta/120)), "units")
        
        def _bind_to_mousewheel(event):
            self.scrollable_frame.master.bind_all("<MouseWheel>", _on_mousewheel)
        
        def _unbind_from_mousewheel(event):
            self.scrollable_frame.master.unbind_all("<MouseWheel>")
        
        # 当鼠标进入滚动区域时绑定滚轮事件
        self.scrollable_frame.bind('<Enter>', _bind_to_mousewheel)
        self.scrollable_frame.bind('<Leave>', _unbind_from_mousewheel)

    def on_closing(self):
        """处理窗口关闭事件"""
        self.is_closed = True
        plt.close(self.fig)  # 关闭matplotlib图形
        self.root.destroy()  # 销毁窗口
        
    def apply_style(self):
        """应用图表样式"""
        self.ax.grid(True, alpha=0.3)
        self.ax.set_axisbelow(True)
        self.ax.spines['top'].set_visible(False)
        self.ax.spines['right'].set_visible(False)
        self.fig.tight_layout()
        
    def create_control_panel(self, var):
        """创建控制面板"""
        # 创建样式
        style = ttk.Style()
        style.configure('TButton', padding=5)
        style.configure('TEntry', padding=5)
        
        # 全局控制选项
        global_frame = LabelFrame(self.scrollable_frame, text="🌐 全局控制")
        global_frame.pack(pady=10, fill='x')
        
        # 应用到所有图表的复选框
        self.apply_to_all_var = BooleanVar(value=False)
        Checkbutton(global_frame, text="应用到所有图表", 
                   variable=self.apply_to_all_var,
                   font=('Arial', 10, 'bold')).pack(anchor='w', padx=5, pady=5)
        
        # 数据系列管理
        series_frame = LabelFrame(self.scrollable_frame, text="📊 数据系列管理")
        series_frame.pack(pady=10, fill='x')
        
        # 创建数据系列列表框架
        self.series_list_frame = Frame(series_frame)
        self.series_list_frame.pack(fill='both', expand=True, pady=5, padx=5)
        
        self.create_series_controls()
        
        # 标签设置
        label_frame = LabelFrame(self.scrollable_frame, text="🏷️ 标签设置")
        label_frame.pack(pady=10, fill='x')
        
        # 数据点大小设置
        size_frame = Frame(label_frame)
        size_frame.pack(fill='x', pady=5, padx=5)
        Label(size_frame, text="数据点大小:").pack(side='left')
        size_entry = Entry(size_frame, textvariable=self.point_size_var, width=10)
        size_entry.pack(side='left', padx=5)
        ttk.Button(size_frame, text="应用大小", command=self.update_point_size).pack(side='left', padx=5)
        
        Label(size_frame, text="样式:").pack(side='left', padx=(10, 0))
        style_combo = ttk.Combobox(size_frame, textvariable=self.marker_style_var, 
                                  values=['o', 's', '^', 'v', '<', '>', 'D', 'd', 'p', '*', 'h', 'H', '+', 'x'], 
                                  width=3, state='readonly')
        style_combo.pack(side='left', padx=5)
        style_combo.bind('<<ComboboxSelected>>', lambda e: self.redraw_plot())
        
        # 图例标签文字大小设置
        legend_size_frame = Frame(label_frame)
        legend_size_frame.pack(fill='x', pady=5, padx=5)
        Label(legend_size_frame, text="图例文字大小:").pack(side='left')
        legend_size_entry = Entry(legend_size_frame, textvariable=self.legend_size_var, width=10)
        legend_size_entry.pack(side='left', padx=5)
        ttk.Button(legend_size_frame, text="应用图例大小", command=self.update_legend_size).pack(side='left', padx=5)
        
        Label(legend_size_frame, text="点缩放:").pack(side='left', padx=(10, 0))
        scale_entry = Entry(legend_size_frame, textvariable=self.legend_marker_scale_var, width=5)
        scale_entry.pack(side='left', padx=5)
        ttk.Button(legend_size_frame, text="应用", command=self.update_legend).pack(side='left', padx=2)
        
        # 标题字体大小设置
        title_size_frame = Frame(label_frame)
        title_size_frame.pack(fill='x', pady=5, padx=5)
        Label(title_size_frame, text="标题字体大小:").pack(side='left')
        title_size_entry = Entry(title_size_frame, textvariable=self.title_size_var, width=10)
        title_size_entry.pack(side='left', padx=5)
        ttk.Button(title_size_frame, text="应用标题大小", command=self.update_title_size).pack(side='left', padx=5)
        
        # X轴标签字体大小设置
        xlabel_size_frame = Frame(label_frame)
        xlabel_size_frame.pack(fill='x', pady=5, padx=5)
        Label(xlabel_size_frame, text="X轴字体大小:").pack(side='left')
        xlabel_size_entry = Entry(xlabel_size_frame, textvariable=self.xlabel_size_var, width=10)
        xlabel_size_entry.pack(side='left', padx=5)
        ttk.Button(xlabel_size_frame, text="应用X轴大小", command=self.update_xlabel_size).pack(side='left', padx=5)
        
        # Y轴标签字体大小设置
        ylabel_size_frame = Frame(label_frame)
        ylabel_size_frame.pack(fill='x', pady=5, padx=5)
        Label(ylabel_size_frame, text="Y轴字体大小:").pack(side='left')
        ylabel_size_entry = Entry(ylabel_size_frame, textvariable=self.ylabel_size_var, width=10)
        ylabel_size_entry.pack(side='left', padx=5)
        ttk.Button(ylabel_size_frame, text="应用Y轴大小", command=self.update_ylabel_size).pack(side='left', padx=5)
        
        # 标签位置设置
        position_frame = Frame(label_frame)
        position_frame.pack(fill='x', pady=5, padx=5)
        Label(position_frame, text="标签位置:").pack(side='left')
        position_combo = ttk.Combobox(position_frame, textvariable=self.label_position_var, 
                                     values=['best', 'upper right', 'upper left', 'lower left', 
                                            'lower right', 'right', 'center left', 'center right', 
                                            'lower center', 'upper center', 'center'], width=15)
        position_combo.pack(side='left', padx=5)
        ttk.Button(position_frame, text="应用位置", command=self.update_label_position).pack(side='left', padx=5)
        
        # 标题控制
        title_frame = LabelFrame(self.scrollable_frame, text="📝 标题设置")
        title_frame.pack(pady=10, fill='x')
        self.title_var = StringVar(value=f'{var} vs time')
        title_entry = Entry(title_frame, textvariable=self.title_var)
        title_entry.pack(side='top', fill='x', pady=5, padx=5)
        ttk.Button(title_frame, text="更新标题", command=self.update_title).pack(side='top', fill='x', padx=5, pady=(0,5))
        
        # X轴标签控制
        xlabel_frame = LabelFrame(self.scrollable_frame, text="📊 X轴设置")
        xlabel_frame.pack(pady=10, fill='x')
        self.xlabel_var = StringVar(value='time')
        xlabel_entry = Entry(xlabel_frame, textvariable=self.xlabel_var)
        xlabel_entry.pack(side='top', fill='x', pady=5, padx=5)
        ttk.Button(xlabel_frame, text="更新X轴", command=self.update_xlabel).pack(side='top', fill='x', padx=5, pady=(0,5))
        
        # Y轴标签控制
        ylabel_frame = LabelFrame(self.scrollable_frame, text="📈 Y轴设置")
        ylabel_frame.pack(pady=10, fill='x')
        self.ylabel_var = StringVar(value=var)
        ylabel_entry = Entry(ylabel_frame, textvariable=self.ylabel_var)
        ylabel_entry.pack(side='top', fill='x', pady=5, padx=5)
        ttk.Button(ylabel_frame, text="更新Y轴", command=self.update_ylabel).pack(side='top', fill='x', padx=5, pady=(0,5))
        
        # 图例编辑
        legend_frame = LabelFrame(self.scrollable_frame, text="🏷️ 图例设置")
        legend_frame.pack(pady=10, fill='x')
        
        for line, folder_name, _ in self.data_series:
            series_frame = Frame(legend_frame)
            series_frame.pack(fill='x', pady=2, padx=5)
            self.legend_vars[folder_name] = StringVar(value=folder_name)
            Entry(series_frame, textvariable=self.legend_vars[folder_name]).pack(side='left', fill='x', expand=True)
        
        ttk.Button(legend_frame, text="更新图例", command=self.update_legend).pack(fill='x', padx=5, pady=5)
        
        # 数据倍率设置
        scale_frame = LabelFrame(self.scrollable_frame, text="⚖️ 数据倍率设置")
        scale_frame.pack(pady=10, fill='x')
        
        for line, folder_name, _ in self.data_series:
            series_frame = Frame(scale_frame)
            series_frame.pack(fill='x', pady=2, padx=5)
            
            # 显示数据系列名称
            Label(series_frame, text=f"{folder_name}:", width=15, anchor='w').pack(side='left')
            
            # 倍率输入框
            self.scale_vars[folder_name] = StringVar(value="1.0")
            scale_entry = Entry(series_frame, textvariable=self.scale_vars[folder_name], width=10)
            scale_entry.pack(side='left', padx=5)
            
            # 重置按钮
            ttk.Button(series_frame, text="重置", 
                      command=lambda fn=folder_name: self.reset_scale(fn)).pack(side='left', padx=2)
        
        ttk.Button(scale_frame, text="应用倍率", command=self.apply_scales).pack(fill='x', padx=5, pady=5)
        
        # Smooth化设置
        smooth_frame = LabelFrame(self.scrollable_frame, text="🌊 Smooth化设置")
        smooth_frame.pack(pady=10, fill='x')
        
        # 全局Smooth化控制
        global_smooth_control = Frame(smooth_frame)
        global_smooth_control.pack(fill='x', pady=5, padx=5)
        
        self.global_smoothing_var = BooleanVar(value=False)
        Checkbutton(global_smooth_control, text="启用全局Smooth化", 
                   variable=self.global_smoothing_var,
                   command=self.toggle_global_smoothing).pack(side='left', padx=5)
        
        # Smooth化方法选择
        method_frame = Frame(smooth_frame)
        method_frame.pack(fill='x', pady=2, padx=5)
        
        Label(method_frame, text="Smooth方法:").pack(side='left')
        self.smoothing_method_var = StringVar(value='gaussian')
        method_combo = ttk.Combobox(method_frame, textvariable=self.smoothing_method_var, 
                                   values=['gaussian', 'spline', 'moving_average'], 
                                   width=12, state='readonly')
        method_combo.pack(side='left', padx=5)
        method_combo.bind('<<ComboboxSelected>>', lambda e: self.update_smoothing_if_enabled())
        
        # Smooth化参数
        param_frame = Frame(smooth_frame)
        param_frame.pack(fill='x', pady=2, padx=5)
        
        Label(param_frame, text="参数:").pack(side='left')
        self.smoothing_param_var = StringVar(value="1.0")
        param_entry = Entry(param_frame, textvariable=self.smoothing_param_var, width=8)
        param_entry.pack(side='left', padx=5)
        
        Label(param_frame, text="插值点数:").pack(side='left', padx=(10,0))
        self.smoothing_points_var = StringVar(value="100")
        points_entry = Entry(param_frame, textvariable=self.smoothing_points_var, width=8)
        points_entry.pack(side='left', padx=5)
        
        ttk.Button(param_frame, text="应用Smooth化", command=self.apply_smoothing_settings).pack(side='left', padx=10)
        
        # 坐标轴类型控制
        axis_frame = LabelFrame(self.scrollable_frame, text="📐 坐标轴类型")
        axis_frame.pack(pady=10, fill='x')
        ttk.Button(axis_frame, text="切换X轴对数", command=lambda: self.toggle_log_scale('x')).pack(fill='x', pady=2, padx=5)
        ttk.Button(axis_frame, text="切换Y轴对数", command=lambda: self.toggle_log_scale('y')).pack(fill='x', pady=2, padx=5)
        
        # 坐标轴单位缩放
        unit_scale_frame = LabelFrame(self.scrollable_frame, text="📏 坐标轴单位缩放")
        unit_scale_frame.pack(pady=10, fill='x')
        
        # X轴单位缩放
        x_unit_frame = Frame(unit_scale_frame)
        x_unit_frame.pack(fill='x', pady=5, padx=5)
        Label(x_unit_frame, text="X轴缩放因子:").pack(side='left')
        x_unit_entry = Entry(x_unit_frame, textvariable=self.x_unit_scale_var, width=10)
        x_unit_entry.pack(side='left', padx=5)
        ttk.Button(x_unit_frame, text="应用", command=self.apply_unit_scale).pack(side='left', padx=5)
        Label(x_unit_frame, text="(例: 1e-5)").pack(side='left', padx=5)
        
        # Y轴单位缩放
        y_unit_frame = Frame(unit_scale_frame)
        y_unit_frame.pack(fill='x', pady=5, padx=5)
        Label(y_unit_frame, text="Y轴缩放因子:").pack(side='left')
        y_unit_entry = Entry(y_unit_frame, textvariable=self.y_unit_scale_var, width=10)
        y_unit_entry.pack(side='left', padx=5)
        ttk.Button(y_unit_frame, text="应用", command=self.apply_unit_scale).pack(side='left', padx=5)
        Label(y_unit_frame, text="(例: 1e-6)").pack(side='left', padx=5)
        
        # 说明文字
        help_frame = Frame(unit_scale_frame)
        help_frame.pack(fill='x', pady=2, padx=5)
        Label(help_frame, text="说明: 数据将除以缩放因子显示", 
              font=('Arial', 8), fg='gray').pack(anchor='w')
        
        # 保存按钮
        save_frame = LabelFrame(self.scrollable_frame, text="💾 保存")
        save_frame.pack(pady=10, fill='x')
        ttk.Button(save_frame, text="💾 保存图片（选择位置）", command=self.save_figure).pack(fill='x', pady=5, padx=5)

    def create_series_controls(self):
        """创建数据系列控制界面"""
        # 清空现有控件
        for widget in self.series_list_frame.winfo_children():
            widget.destroy()
        
        # 为每个数据系列创建控制行
        for i, folder_name in enumerate(self.data_order):
            # 创建系列控制框架
            series_control_frame = Frame(self.series_list_frame)
            series_control_frame.pack(fill='x', pady=2, padx=5)
            
            # 显示/隐藏复选框
            if folder_name not in self.visibility_vars:
                self.visibility_vars[folder_name] = BooleanVar(value=True)
            
            visibility_cb = Checkbutton(series_control_frame, text="显示", 
                                      variable=self.visibility_vars[folder_name],
                                      command=self.update_visibility)
            visibility_cb.pack(side='left', padx=(0, 5))
            
            # 系列名称标签
            name_label = Label(series_control_frame, text=folder_name, width=15, anchor='w')
            name_label.pack(side='left', padx=5)
            
            # 上移按钮
            up_btn = ttk.Button(series_control_frame, text="↑", width=3,
                               command=lambda idx=i: self.move_series_up(idx))
            up_btn.pack(side='right', padx=2)
            
            # 下移按钮
            down_btn = ttk.Button(series_control_frame, text="↓", width=3,
                                 command=lambda idx=i: self.move_series_down(idx))
            down_btn.pack(side='right', padx=2)
            
            # 禁用第一个的上移按钮和最后一个的下移按钮
            if i == 0:
                up_btn.config(state='disabled')
            if i == len(self.data_order) - 1:
                down_btn.config(state='disabled')

    def move_series_up(self, index):
        """将数据系列向上移动"""
        if index > 0:
            # 交换顺序
            self.data_order[index], self.data_order[index-1] = self.data_order[index-1], self.data_order[index]
            # 重新创建控制界面
            self.create_series_controls()
            # 重新绘制图表
            self.redraw_plot()

    def move_series_down(self, index):
        """将数据系列向下移动"""
        if index < len(self.data_order) - 1:
            # 交换顺序
            self.data_order[index], self.data_order[index+1] = self.data_order[index+1], self.data_order[index]
            # 重新创建控制界面
            self.create_series_controls()
            # 重新绘制图表
            self.redraw_plot()

    def update_visibility(self):
        """更新数据系列的显示/隐藏状态"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            global_controller.apply_to_all(self._update_visibility_single)
        else:
            # 仅应用到当前图表
            self._update_visibility_single()

    def _update_visibility_single(self, *args, **kwargs):
        """在单个图表上更新显示/隐藏状态"""
        self.redraw_plot()

    def update_point_size(self):
        """更新数据点大小"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            size_value = self.point_size_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_point_size_single(size_value))
        else:
            # 仅应用到当前图表
            self._update_point_size_single(self.point_size_var.get())

    def _update_point_size_single(self, size_value=None, *args, **kwargs):
        """在单个图表上更新数据点大小"""
        if size_value is None:
            size_value = self.point_size_var.get()
        else:
            self.point_size_var.set(size_value)
        
        try:
            size = int(size_value)
            # 只更新数据点大小，不重新绘制整个图表
            for folder_name in self.data_order:
                if self.visibility_vars[folder_name].get():
                    # 找到对应的散点图对象并更新大小
                    for line, series_name, _ in self.data_series:
                        if series_name == folder_name:
                            line.set_sizes([size] * len(line.get_offsets()))
                            break
            self.canvas.draw()
        except ValueError:
            messagebox.showerror("错误", "请输入有效的数字作为数据点大小")

    def update_legend_size(self):
        """更新图例标签文字大小"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            size_value = self.legend_size_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_legend_size_single(size_value))
        else:
            # 仅应用到当前图表
            self._update_legend_size_single(self.legend_size_var.get())

    def _update_legend_size_single(self, size_value=None, *args, **kwargs):
        """在单个图表上更新图例标签文字大小"""
        if size_value is None:
            size_value = self.legend_size_var.get()
        else:
            self.legend_size_var.set(size_value)
        
        try:
            size = int(size_value)
            # 重新创建图例以更新字体大小
            legend = self.ax.get_legend()
            if legend:
                # 移除旧图例
                legend.remove()
                # 重新创建图例（只有在有带标签的图形元素时）
                handles, labels = self.ax.get_legend_handles_labels()
                if handles:  # 如果有图例元素
                    self.ax.legend(loc=self.label_position_var.get(), prop={'size': size})
            self.canvas.draw()
        except ValueError:
            messagebox.showerror("错误", "请输入有效的数字作为图例文字大小")

    def update_title_size(self):
        """更新标题字体大小"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            size_value = self.title_size_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_title_size_single(size_value))
        else:
            # 仅应用到当前图表
            self._update_title_size_single(self.title_size_var.get())

    def _update_title_size_single(self, size_value=None, *args, **kwargs):
        """在单个图表上更新标题字体大小"""
        if size_value is None:
            size_value = self.title_size_var.get()
        else:
            self.title_size_var.set(size_value)
        
        try:
            size = int(size_value)
            # 重新设置标题字体大小
            self.ax.set_title(self.title_var.get(), fontsize=size)
            self.canvas.draw()
        except ValueError:
            messagebox.showerror("错误", "请输入有效的数字作为标题字体大小")

    def update_xlabel_size(self):
        """更新X轴标签字体大小"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            size_value = self.xlabel_size_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_xlabel_size_single(size_value))
        else:
            # 仅应用到当前图表
            self._update_xlabel_size_single(self.xlabel_size_var.get())

    def _update_xlabel_size_single(self, size_value=None, *args, **kwargs):
        """在单个图表上更新X轴标签字体大小"""
        if size_value is None:
            size_value = self.xlabel_size_var.get()
        else:
            self.xlabel_size_var.set(size_value)
        
        try:
            size = int(size_value)
            # 重新设置X轴标签字体大小
            self.ax.set_xlabel(self.xlabel_var.get(), fontsize=size)
            self.canvas.draw()
        except ValueError:
            messagebox.showerror("错误", "请输入有效的数字作为X轴字体大小")

    def update_ylabel_size(self):
        """更新Y轴标签字体大小"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            size_value = self.ylabel_size_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_ylabel_size_single(size_value))
        else:
            # 仅应用到当前图表
            self._update_ylabel_size_single(self.ylabel_size_var.get())

    def _update_ylabel_size_single(self, size_value=None, *args, **kwargs):
        """在单个图表上更新Y轴标签字体大小"""
        if size_value is None:
            size_value = self.ylabel_size_var.get()
        else:
            self.ylabel_size_var.set(size_value)
        
        try:
            size = int(size_value)
            # 重新设置Y轴标签字体大小
            self.ax.set_ylabel(self.ylabel_var.get(), fontsize=size)
            self.canvas.draw()
        except ValueError:
            messagebox.showerror("错误", "请输入有效的数字作为Y轴字体大小")

    def update_label_position(self):
        """更新标签位置"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            position_value = self.label_position_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_label_position_single(position_value))
        else:
            # 仅应用到当前图表
            self._update_label_position_single(self.label_position_var.get())

    def _update_label_position_single(self, position_value=None, *args, **kwargs):
        """在单个图表上更新标签位置"""
        if position_value is None:
            position_value = self.label_position_var.get()
        else:
            self.label_position_var.set(position_value)
        
        # 重新创建图例以更新位置
        legend = self.ax.get_legend()
        if legend:
            # 获取当前字体大小
            try:
                legend_size = int(self.legend_size_var.get())
            except ValueError:
                legend_size = 10
            
            # 移除旧图例
            legend.remove()
            # 重新创建图例（只有在有带标签的图形元素时）
            handles, labels = self.ax.get_legend_handles_labels()
            if handles:  # 如果有图例元素
                self.ax.legend(loc=position_value, prop={'size': legend_size})
        self.canvas.draw()

    def redraw_plot(self):
        """重新绘制图表"""
        # 清除当前图表
        self.ax.clear()
        
        # 获取单位缩放因子
        try:
            x_unit_scale = float(self.x_unit_scale_var.get())
            if x_unit_scale == 0:
                x_unit_scale = 1.0
        except (ValueError, ZeroDivisionError):
            x_unit_scale = 1.0
            self.x_unit_scale_var.set("1")
        
        try:
            y_unit_scale = float(self.y_unit_scale_var.get())
            if y_unit_scale == 0:
                y_unit_scale = 1.0
        except (ValueError, ZeroDivisionError):
            y_unit_scale = 1.0
            self.y_unit_scale_var.set("1")
        
        # 按顺序重新绘制可见的数据系列
        for folder_name in self.data_order:
            if self.visibility_vars[folder_name].get():  # 只绘制可见的系列
                # 获取原始数据
                x_data, y_data = self.original_data[folder_name]
                
                # 应用单位缩放
                x_data_scaled = x_data / x_unit_scale
                
                # 应用倍率
                try:
                    scale_factor = float(self.scale_vars[folder_name].get())
                except ValueError:
                    scale_factor = 1.0
                    self.scale_vars[folder_name].set("1.0")
                
                scaled_y_data = (y_data * scale_factor) / y_unit_scale
                
                # 获取图例名称
                legend_name = self.legend_vars[folder_name].get()
                
                # 使用固定的颜色绘制数据
                color = self.series_colors[folder_name]
                
                # 获取数据点大小
                try:
                    point_size = int(self.point_size_var.get())
                except ValueError:
                    point_size = 10
                
                marker = self.marker_style_var.get()
                
                # 应用Smooth化
                x_smooth, y_smooth = self.apply_single_smoothing(x_data_scaled, scaled_y_data)
                
                # 绘制数据
                if self.global_smoothing_var.get():
                    # 绘制原始数据点（较淡）
                    self.ax.scatter(x_data_scaled, scaled_y_data, s=point_size, alpha=0.4, color=color, edgecolors='none', marker=marker, label=f"{legend_name} (原始)")
                    # 生成Smooth曲线的颜色（更深或更亮的版本）
                    import matplotlib.colors as mcolors
                    # 将颜色转换为RGB，然后调整亮度
                    rgb = mcolors.to_rgb(color)
                    # 使颜色更深（乘以0.7）
                    smooth_color = tuple(max(0, c * 0.7) for c in rgb)
                    # 绘制Smooth曲线
                    self.ax.plot(x_smooth, y_smooth, label=f"{legend_name} (Smooth)", linewidth=2.5, alpha=0.9, color=smooth_color)
                else:
                    # 原始散点图
                    self.ax.scatter(x_data_scaled, scaled_y_data, label=legend_name, s=point_size, alpha=0.7, color=color, marker=marker)
        
        # 重新设置轴标签和标题
        try:
            title_size = int(self.title_size_var.get())
        except ValueError:
            title_size = 12
        
        try:
            xlabel_size = int(self.xlabel_size_var.get())
        except ValueError:
            xlabel_size = 10
        
        try:
            ylabel_size = int(self.ylabel_size_var.get())
        except ValueError:
            ylabel_size = 10
        
        self.ax.set_xlabel(self.xlabel_var.get(), fontsize=xlabel_size)
        self.ax.set_ylabel(self.ylabel_var.get(), fontsize=ylabel_size)
        self.ax.set_title(self.title_var.get(), fontsize=title_size)
        
        # 设置图例位置和大小
        try:
            legend_size = int(self.legend_size_var.get())
        except ValueError:
            legend_size = 10
        
        # 只有在有带标签的图形元素时才创建图例
        handles, labels = self.ax.get_legend_handles_labels()
        if handles:  # 如果有图例元素
            try:
                legend_scale = float(self.legend_marker_scale_var.get())
            except ValueError:
                legend_scale = 1.0
            self.ax.legend(loc=self.label_position_var.get(), prop={'size': legend_size}, markerscale=legend_scale)
        
        # 重新应用样式
        self.apply_style()
        
        # 刷新画布
        self.canvas.draw()

    def reset_scale(self, folder_name):
        """重置指定数据系列的倍率为1.0"""
        self.scale_vars[folder_name].set("1.0")

    def apply_scales(self):
        """应用所有数据倍率"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            global_controller.apply_to_all(self._apply_scales_single)
        else:
            # 仅应用到当前图表
            self._apply_scales_single()

    def _apply_scales_single(self, *args, **kwargs):
        """在单个图表上应用倍率"""
        self.redraw_plot()

    def update_legend(self):
        """更新图例标签"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            global_controller.apply_to_all(self._update_legend_single)
        else:
            # 仅应用到当前图表
            self._update_legend_single()

    def _update_legend_single(self, *args, **kwargs):
        """在单个图表上更新图例"""
        # 重新创建图例以更新标签
        legend = self.ax.get_legend()
        if legend:
            # 获取当前字体大小和位置
            try:
                legend_size = int(self.legend_size_var.get())
            except ValueError:
                legend_size = 10
            
            # 移除旧图例
            legend.remove()
            # 重新创建图例（只有在有带标签的图形元素时）
            handles, labels = self.ax.get_legend_handles_labels()
            if handles:  # 如果有图例元素
                try:
                    legend_scale = float(self.legend_marker_scale_var.get())
                except ValueError:
                    legend_scale = 1.0
                self.ax.legend(loc=self.label_position_var.get(), prop={'size': legend_size}, markerscale=legend_scale)
        
        self.canvas.draw()

    def update_title(self):
        if self.apply_to_all_var.get():
            # 应用到所有图表
            title_text = self.title_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_title_single(title_text))
        else:
            # 仅应用到当前图表
            self._update_title_single(self.title_var.get())

    def _update_title_single(self, title_text=None, *args, **kwargs):
        """在单个图表上更新标题"""
        if title_text is None:
            title_text = self.title_var.get()
        else:
            self.title_var.set(title_text)
        
        try:
            title_size = int(self.title_size_var.get())
        except ValueError:
            title_size = 12
        
        self.ax.set_title(title_text, fontsize=title_size)
        self.canvas.draw()
    
    def update_xlabel(self):
        if self.apply_to_all_var.get():
            # 应用到所有图表
            xlabel_text = self.xlabel_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_xlabel_single(xlabel_text))
        else:
            # 仅应用到当前图表
            self._update_xlabel_single(self.xlabel_var.get())

    def _update_xlabel_single(self, xlabel_text=None, *args, **kwargs):
        """在单个图表上更新X轴标签"""
        if xlabel_text is None:
            xlabel_text = self.xlabel_var.get()
        else:
            self.xlabel_var.set(xlabel_text)
        
        try:
            xlabel_size = int(self.xlabel_size_var.get())
        except ValueError:
            xlabel_size = 10
        
        self.ax.set_xlabel(xlabel_text, fontsize=xlabel_size)
        self.canvas.draw()

    def update_ylabel(self):
        if self.apply_to_all_var.get():
            # 应用到所有图表
            ylabel_text = self.ylabel_var.get()
            global_controller.apply_to_all(lambda controller: controller._update_ylabel_single(ylabel_text))
        else:
            # 仅应用到当前图表
            self._update_ylabel_single(self.ylabel_var.get())

    def _update_ylabel_single(self, ylabel_text=None, *args, **kwargs):
        """在单个图表上更新Y轴标签"""
        if ylabel_text is None:
            ylabel_text = self.ylabel_var.get()
        else:
            self.ylabel_var.set(ylabel_text)
        
        try:
            ylabel_size = int(self.ylabel_size_var.get())
        except ValueError:
            ylabel_size = 10
        
        self.ax.set_ylabel(ylabel_text, fontsize=ylabel_size)
        self.canvas.draw()
    
    def toggle_log_scale(self, axis):
        if self.apply_to_all_var.get():
            # 应用到所有图表
            global_controller.apply_to_all(lambda controller: controller._toggle_log_scale_single(axis))
        else:
            # 仅应用到当前图表
            self._toggle_log_scale_single(axis)

    def _toggle_log_scale_single(self, axis, *args, **kwargs):
        """在单个图表上切换对数坐标"""
        if axis == 'x':
            if self.ax.get_xscale() == 'log':
                self.ax.set_xscale('linear')
            else:
                self.ax.set_xscale('log')
        else:
            if self.ax.get_yscale() == 'log':
                self.ax.set_yscale('linear')
            else:
                self.ax.set_yscale('log')
        self.canvas.draw()
    
    def apply_unit_scale(self):
        """应用坐标轴单位缩放"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            x_scale_value = self.x_unit_scale_var.get()
            y_scale_value = self.y_unit_scale_var.get()
            global_controller.apply_to_all(lambda controller: controller._apply_unit_scale_single(x_scale_value, y_scale_value))
        else:
            # 仅应用到当前图表
            self._apply_unit_scale_single()

    def _apply_unit_scale_single(self, x_scale_value=None, y_scale_value=None, *args, **kwargs):
        """在单个图表上应用单位缩放"""
        if x_scale_value is not None:
            self.x_unit_scale_var.set(x_scale_value)
        if y_scale_value is not None:
            self.y_unit_scale_var.set(y_scale_value)
        
        self.redraw_plot()
    
    def toggle_global_smoothing(self):
        """切换全局Smooth化"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            global_controller.apply_to_all(self._toggle_global_smoothing_single)
        else:
            # 仅应用到当前图表
            self._toggle_global_smoothing_single()

    def _toggle_global_smoothing_single(self, *args, **kwargs):
        """在单个图表上切换全局Smooth化"""
        self.redraw_plot()

    def update_smoothing_if_enabled(self):
        """如果启用了Smooth化，则更新图表"""
        if self.global_smoothing_var.get():
            self.apply_smoothing_settings()

    def apply_smoothing_settings(self):
        """应用Smooth化设置"""
        if self.apply_to_all_var.get():
            # 应用到所有图表
            global_controller.apply_to_all(self._apply_smoothing_settings_single)
        else:
            # 仅应用到当前图表
            self._apply_smoothing_settings_single()

    def _apply_smoothing_settings_single(self, *args, **kwargs):
        """在单个图表上应用Smooth化设置"""
        self.redraw_plot()

    def apply_single_smoothing(self, x_data, y_data):
        """为单变量图表应用Smooth化处理"""
        if not self.global_smoothing_var.get():
            return x_data, y_data
        
        method = self.smoothing_method_var.get()
        
        try:
            if method == 'gaussian':
                # 高斯滤波
                sigma = float(self.smoothing_param_var.get())
                smooth_y = gaussian_filter1d(y_data, sigma=sigma)
                return x_data, smooth_y
            
            elif method == 'spline':
                # 样条插值
                num_points = int(self.smoothing_points_var.get())
                # 创建样条插值函数
                tck = interpolate.splrep(x_data, y_data, s=0)
                # 生成更密集的x点
                x_smooth = np.linspace(x_data.min(), x_data.max(), num_points)
                y_smooth = interpolate.splev(x_smooth, tck)
                return x_smooth, y_smooth
                
            elif method == 'moving_average':
                # 移动平均
                window = int(float(self.smoothing_param_var.get()))
                if window < 1:
                    window = 1
                # 计算移动平均
                smooth_y = np.convolve(y_data, np.ones(window)/window, mode='same')
                return x_data, smooth_y
        
        except (ValueError, TypeError):
            # 参数错误时返回原始数据
            return x_data, y_data
        
        return x_data, y_data

    def save_figure(self):
        """保存图片"""
        try:
            # 设置默认文件名
            default_filename = f"{self.var}_vs_time.png"
            
            # 打开文件保存对话框
            save_path = filedialog.asksaveasfilename(
                title="保存图片",
                defaultextension=".png",
                initialfile=default_filename,
                filetypes=[
                    ("PNG图片", "*.png"),
                    ("JPEG图片", "*.jpg"),
                    ("PDF文档", "*.pdf"),
                    ("SVG矢量图", "*.svg"),
                    ("所有文件", "*.*")
                ]
            )
            
            if save_path:
                # 根据文件扩展名确定保存格式
                file_ext = Path(save_path).suffix.lower()
                if file_ext == '.jpg' or file_ext == '.jpeg':
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='jpeg')
                elif file_ext == '.pdf':
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='pdf')
                elif file_ext == '.svg':
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='svg')
                else:
                    # 默认保存为PNG
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='png')
                
                print(f"图片已保存至: {save_path}")
                messagebox.showinfo("保存成功", f"图片已成功保存至:\n{save_path}")
            else:
                print("保存图片已取消")
        except Exception as e:
            print(f"保存图片时出错: {e}")
            messagebox.showerror("保存失败", f"保存图片时出错:\n{str(e)}")

# 多变量图表控制器
class MultiVariablePlotController:
    def __init__(self, variables, master=None):
        # 创建主窗口
        self.root = Tk() if master is None else Toplevel(master)
        self.root.title(f"多变量图表 - {', '.join(variables[:3])}")
        self.root.geometry("1400x950")
        
        # 设置窗口关闭处理
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # 存储变量
        self.variables = variables
        self.data_list = data_list
        
        # 创建左右分栏
        self.left_frame = Frame(self.root)
        self.left_frame.pack(side='left', fill='both', expand=True, padx=10, pady=10)
        
        self.right_frame = Frame(self.root)
        self.right_frame.pack(side='right', fill='y', padx=10, pady=10)
        
        # 创建滚动条
        canvas = tk.Canvas(self.right_frame)
        scrollbar = ttk.Scrollbar(self.right_frame, orient="vertical", command=canvas.yview)
        self.scrollable_frame = ttk.Frame(canvas)
        
        self.scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # 创建图表
        self.fig = plt.figure(figsize=(12, 8))
        self.ax1 = self.fig.add_subplot(111)
        self.ax2 = None  # 双Y轴的第二个轴
        
        # 将matplotlib图表嵌入到tkinter窗口中
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.left_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        
        # 存储变量状态
        self.variable_states = {}  # 存储每个变量的显示状态、Y轴分配、颜色等
        self.smoothing_states = {}  # 存储Smooth化状态
        self.smoothing_params = {}  # 存储Smooth化参数
        self.dual_y_enabled = BooleanVar(value=False)
        self.data_source_visibility = {}  # 存储每个数据源的显示状态
        
        # 坐标轴单位缩放因子
        self.x_unit_scale_var = StringVar(value="1")
        self.y_unit_scale_var = StringVar(value="1")
        self.left_y_unit_scale_var = StringVar(value="1")
        self.right_y_unit_scale_var = StringVar(value="1")
        self.legend_marker_scale_var = StringVar(value="1.0")
        
        # 预定义颜色列表
        self.color_list = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
                          '#a6cee3', '#fb9a99', '#fdbf6f', '#cab2d6', '#ffff99']
        
        # 初始化变量状态
        self.init_variable_states()
        
        # 创建控制面板
        self.create_control_panel()
        
        # 初始绘制
        self.redraw_plot()
    
        # 标记窗口状态
        self.is_closed = False
        
        # 设置字体 [[memory:6787997]]
        plt.rcParams['font.family'] = ['Times New Roman', 'SimHei']  # 英文用Times New Roman，中文用SimHei
        
        # 绑定鼠标滚轮事件
        self.bind_mousewheel()

    def init_variable_states(self):
        """初始化变量状态"""
        for i, var in enumerate(self.variables):
            self.variable_states[var] = {
                'visible': BooleanVar(value=True),
                'y_axis': StringVar(value='left'),  # 'left' or 'right'
                'color': self.color_list[i % len(self.color_list)],
                'marker': StringVar(value='o'),
                'smoothing': BooleanVar(value=False),
                'point_size': StringVar(value="10"),
                'line_width': StringVar(value="1.5"),
                'alpha': StringVar(value="0.7")
            }
            self.smoothing_states[var] = BooleanVar(value=False)
            self.smoothing_params[var] = {
                'method': StringVar(value='gaussian'),  # 'gaussian', 'spline', 'moving_average'
                'sigma': StringVar(value="1.0"),
                'points': StringVar(value="100")
            }
        
        # 初始化数据源可见性状态
        for folder_name, _ in self.data_list:
            self.data_source_visibility[folder_name] = BooleanVar(value=True)

    def bind_mousewheel(self):
        """绑定鼠标滚轮事件"""
        def _on_mousewheel(event):
            self.scrollable_frame.master.yview_scroll(int(-1*(event.delta/120)), "units")
        
        def _bind_to_mousewheel(event):
            self.scrollable_frame.master.bind_all("<MouseWheel>", _on_mousewheel)
        
        def _unbind_from_mousewheel(event):
            self.scrollable_frame.master.unbind_all("<MouseWheel>")
        
        self.scrollable_frame.bind('<Enter>', _bind_to_mousewheel)
        self.scrollable_frame.bind('<Leave>', _unbind_from_mousewheel)

    def create_control_panel(self):
        """创建控制面板"""
        # 双Y轴控制
        dual_y_frame = LabelFrame(self.scrollable_frame, text="📊 双Y轴设置")
        dual_y_frame.pack(pady=10, fill='x')
        
        Checkbutton(dual_y_frame, text="启用双Y轴", 
                   variable=self.dual_y_enabled,
                   command=self.toggle_dual_y_axis).pack(anchor='w', padx=5, pady=5)
        
        # 数据源管理区域
        data_source_frame = LabelFrame(self.scrollable_frame, text="📁 数据源管理")
        data_source_frame.pack(pady=10, fill='x')
        
        # 数据源全选/取消全选按钮
        data_source_btn_frame = Frame(data_source_frame)
        data_source_btn_frame.pack(fill='x', pady=5, padx=5)
        
        ttk.Button(data_source_btn_frame, text="全部显示", 
                  command=self.show_all_data_sources).pack(side='left', padx=5)
        ttk.Button(data_source_btn_frame, text="全部隐藏", 
                  command=self.hide_all_data_sources).pack(side='left', padx=5)
        
        # 为每个数据源创建复选框
        for folder_name, _ in self.data_list:
            ds_frame = Frame(data_source_frame)
            ds_frame.pack(fill='x', pady=2, padx=5)
            
            Checkbutton(ds_frame, text=folder_name, 
                       variable=self.data_source_visibility[folder_name],
                       command=self.redraw_plot).pack(anchor='w')
        
        # 变量控制区域
        vars_frame = LabelFrame(self.scrollable_frame, text="📈 变量设置")
        vars_frame.pack(pady=10, fill='x')
        
        # 为每个变量创建控制组
        for var in self.variables:
            self.create_variable_control_group(vars_frame, var)
        
        # Smooth化全局设置
        smooth_frame = LabelFrame(self.scrollable_frame, text="🌊 Smooth化设置")
        smooth_frame.pack(pady=10, fill='x')
        
        # 全局Smooth化控制
        global_smooth_frame = Frame(smooth_frame)
        global_smooth_frame.pack(fill='x', pady=5, padx=5)
        
        ttk.Button(global_smooth_frame, text="全部启用Smooth化", 
                  command=self.enable_all_smoothing).pack(side='left', padx=5)
        ttk.Button(global_smooth_frame, text="全部禁用Smooth化", 
                  command=self.disable_all_smoothing).pack(side='left', padx=5)
        
        # 标题和标签设置
        label_frame = LabelFrame(self.scrollable_frame, text="🏷️ 标题和标签设置")
        label_frame.pack(pady=10, fill='x')
        
        # 标题设置
        title_frame = Frame(label_frame)
        title_frame.pack(fill='x', pady=2, padx=5)
        Label(title_frame, text="图表标题:").pack(side='left')
        self.title_var = StringVar(value=f"多变量图表: {', '.join(self.variables)}")
        title_entry = Entry(title_frame, textvariable=self.title_var, width=30)
        title_entry.pack(side='left', padx=5, fill='x', expand=True)
        
        # X轴标签设置
        xlabel_frame = Frame(label_frame)
        xlabel_frame.pack(fill='x', pady=2, padx=5)
        Label(xlabel_frame, text="X轴标签:").pack(side='left')
        self.xlabel_var = StringVar(value="时间")
        xlabel_entry = Entry(xlabel_frame, textvariable=self.xlabel_var, width=20)
        xlabel_entry.pack(side='left', padx=5)
        
        # Y轴标签设置
        ylabel_frame = Frame(label_frame)
        ylabel_frame.pack(fill='x', pady=2, padx=5)
        
        # 单Y轴标签（当未启用双Y轴时显示）
        self.single_ylabel_frame = Frame(ylabel_frame)
        self.single_ylabel_frame.pack(fill='x')
        Label(self.single_ylabel_frame, text="Y轴标签:").pack(side='left')
        self.ylabel_var = StringVar(value="数值")
        ylabel_entry = Entry(self.single_ylabel_frame, textvariable=self.ylabel_var, width=20)
        ylabel_entry.pack(side='left', padx=5)
        
        # 双Y轴标签（当启用双Y轴时显示）
        self.dual_ylabel_frame = Frame(ylabel_frame)
        
        # 左Y轴标签
        left_ylabel_frame = Frame(self.dual_ylabel_frame)
        left_ylabel_frame.pack(fill='x', pady=1)
        Label(left_ylabel_frame, text="左Y轴标签:", fg='blue').pack(side='left')
        self.left_ylabel_var = StringVar(value="左Y轴数值")
        left_ylabel_entry = Entry(left_ylabel_frame, textvariable=self.left_ylabel_var, width=20)
        left_ylabel_entry.pack(side='left', padx=5)
        
        # 右Y轴标签
        right_ylabel_frame = Frame(self.dual_ylabel_frame)
        right_ylabel_frame.pack(fill='x', pady=1)
        Label(right_ylabel_frame, text="右Y轴标签:", fg='red').pack(side='left')
        self.right_ylabel_var = StringVar(value="右Y轴数值")
        right_ylabel_entry = Entry(right_ylabel_frame, textvariable=self.right_ylabel_var, width=20)
        right_ylabel_entry.pack(side='left', padx=5)
        
        # 字体大小设置
        font_frame = LabelFrame(self.scrollable_frame, text="🔤 字体大小设置")
        font_frame.pack(pady=10, fill='x')
        
        # 标题字体大小
        title_font_frame = Frame(font_frame)
        title_font_frame.pack(fill='x', pady=2, padx=5)
        Label(title_font_frame, text="标题字体:").pack(side='left')
        self.title_fontsize_var = StringVar(value="14")
        Entry(title_font_frame, textvariable=self.title_fontsize_var, width=6).pack(side='left', padx=5)
        
        # 轴标签字体大小
        axis_font_frame = Frame(font_frame)
        axis_font_frame.pack(fill='x', pady=2, padx=5)
        Label(axis_font_frame, text="轴标签字体:").pack(side='left')
        self.axis_fontsize_var = StringVar(value="12")
        Entry(axis_font_frame, textvariable=self.axis_fontsize_var, width=6).pack(side='left', padx=5)
        
        Label(axis_font_frame, text="刻度字体:").pack(side='left', padx=(10,0))
        self.tick_fontsize_var = StringVar(value="10")
        Entry(axis_font_frame, textvariable=self.tick_fontsize_var, width=6).pack(side='left', padx=5)
        
        # 图例字体大小
        legend_font_frame = Frame(font_frame)
        legend_font_frame.pack(fill='x', pady=2, padx=5)
        Label(legend_font_frame, text="图例字体:").pack(side='left')
        self.legend_fontsize_var = StringVar(value="10")
        Entry(legend_font_frame, textvariable=self.legend_fontsize_var, width=6).pack(side='left', padx=5)

        Label(legend_font_frame, text="点缩放:").pack(side='left', padx=(10,0))
        Entry(legend_font_frame, textvariable=self.legend_marker_scale_var, width=6).pack(side='left', padx=5)
        
        # 图例位置设置
        legend_pos_frame = LabelFrame(self.scrollable_frame, text="📌 图例位置设置")
        legend_pos_frame.pack(pady=10, fill='x')
        
        # 单Y轴图例位置
        self.single_legend_frame = Frame(legend_pos_frame)
        self.single_legend_frame.pack(fill='x', pady=2, padx=5)
        Label(self.single_legend_frame, text="图例位置:").pack(side='left')
        self.legend_position_var = StringVar(value="best")
        legend_combo = ttk.Combobox(self.single_legend_frame, textvariable=self.legend_position_var, 
                                   values=['best', 'upper right', 'upper left', 'lower left', 
                                          'lower right', 'right', 'center left', 'center right', 
                                          'lower center', 'upper center', 'center'], width=12, state='readonly')
        legend_combo.pack(side='left', padx=5)
        legend_combo.bind('<<ComboboxSelected>>', lambda e: self.redraw_plot())
        
        # 双Y轴图例位置
        self.dual_legend_frame = Frame(legend_pos_frame)
        
        # 左Y轴图例位置
        left_legend_frame = Frame(self.dual_legend_frame)
        left_legend_frame.pack(fill='x', pady=1, padx=5)
        Label(left_legend_frame, text="左Y轴图例:", fg='blue').pack(side='left')
        self.left_legend_position_var = StringVar(value="upper left")
        left_legend_combo = ttk.Combobox(left_legend_frame, textvariable=self.left_legend_position_var, 
                                        values=['best', 'upper right', 'upper left', 'lower left', 
                                               'lower right', 'right', 'center left', 'center right', 
                                               'lower center', 'upper center', 'center'], width=12, state='readonly')
        left_legend_combo.pack(side='left', padx=5)
        left_legend_combo.bind('<<ComboboxSelected>>', lambda e: self.redraw_plot())
        
        # 右Y轴图例位置
        right_legend_frame = Frame(self.dual_legend_frame)
        right_legend_frame.pack(fill='x', pady=1, padx=5)
        Label(right_legend_frame, text="右Y轴图例:", fg='red').pack(side='left')
        self.right_legend_position_var = StringVar(value="upper right")
        right_legend_combo = ttk.Combobox(right_legend_frame, textvariable=self.right_legend_position_var, 
                                         values=['best', 'upper right', 'upper left', 'lower left', 
                                                'lower right', 'right', 'center left', 'center right', 
                                                'lower center', 'upper center', 'center'], width=12, state='readonly')
        right_legend_combo.pack(side='left', padx=5)
        right_legend_combo.bind('<<ComboboxSelected>>', lambda e: self.redraw_plot())
        
        # 坐标轴类型控制
        axis_frame = LabelFrame(self.scrollable_frame, text="📐 坐标轴类型")
        axis_frame.pack(pady=10, fill='x')
        
        # X轴控制
        x_axis_frame = Frame(axis_frame)
        x_axis_frame.pack(fill='x', pady=2, padx=5)
        ttk.Button(x_axis_frame, text="切换X轴对数", 
                  command=lambda: self.toggle_log_scale('x')).pack(side='left', padx=5)
        
        # 单Y轴控制（当未启用双Y轴时显示）
        self.single_y_axis_frame = Frame(axis_frame)
        self.single_y_axis_frame.pack(fill='x', pady=2, padx=5)
        ttk.Button(self.single_y_axis_frame, text="切换Y轴对数", 
                  command=lambda: self.toggle_log_scale('y')).pack(side='left', padx=5)
        
        # 双Y轴控制（当启用双Y轴时显示）
        self.dual_y_axis_frame = Frame(axis_frame)
        
        # 左Y轴控制
        left_y_axis_frame = Frame(self.dual_y_axis_frame)
        left_y_axis_frame.pack(fill='x', pady=1, padx=5)
        ttk.Button(left_y_axis_frame, text="切换左Y轴对数", 
                  command=lambda: self.toggle_log_scale('left_y')).pack(side='left', padx=5)
        Label(left_y_axis_frame, text="(蓝色)", fg='blue').pack(side='left', padx=5)
        
        # 右Y轴控制
        right_y_axis_frame = Frame(self.dual_y_axis_frame)
        right_y_axis_frame.pack(fill='x', pady=1, padx=5)
        ttk.Button(right_y_axis_frame, text="切换右Y轴对数", 
                  command=lambda: self.toggle_log_scale('right_y')).pack(side='left', padx=5)
        Label(right_y_axis_frame, text="(红色)", fg='red').pack(side='left', padx=5)
        
        # 坐标轴范围设置
        range_frame = LabelFrame(self.scrollable_frame, text="📏 坐标轴范围设置")
        range_frame.pack(pady=10, fill='x')
        
        # X轴范围设置
        x_range_frame = Frame(range_frame)
        x_range_frame.pack(fill='x', pady=2, padx=5)
        Label(x_range_frame, text="X轴范围:").pack(side='left')
        Label(x_range_frame, text="最小值:").pack(side='left', padx=(10,0))
        self.x_min_var = StringVar(value="auto")
        Entry(x_range_frame, textvariable=self.x_min_var, width=8).pack(side='left', padx=2)
        Label(x_range_frame, text="最大值:").pack(side='left', padx=(5,0))
        self.x_max_var = StringVar(value="auto")
        Entry(x_range_frame, textvariable=self.x_max_var, width=8).pack(side='left', padx=2)
        ttk.Button(x_range_frame, text="应用X轴", command=self.apply_axis_limits).pack(side='left', padx=5)
        
        # 单Y轴范围设置
        self.single_y_range_frame = Frame(range_frame)
        self.single_y_range_frame.pack(fill='x', pady=2, padx=5)
        Label(self.single_y_range_frame, text="Y轴范围:").pack(side='left')
        Label(self.single_y_range_frame, text="最小值:").pack(side='left', padx=(10,0))
        self.y_min_var = StringVar(value="auto")
        Entry(self.single_y_range_frame, textvariable=self.y_min_var, width=8).pack(side='left', padx=2)
        Label(self.single_y_range_frame, text="最大值:").pack(side='left', padx=(5,0))
        self.y_max_var = StringVar(value="auto")
        Entry(self.single_y_range_frame, textvariable=self.y_max_var, width=8).pack(side='left', padx=2)
        ttk.Button(self.single_y_range_frame, text="应用Y轴", command=self.apply_axis_limits).pack(side='left', padx=5)
        
        # 双Y轴范围设置
        self.dual_y_range_frame = Frame(range_frame)
        
        # 左Y轴范围
        left_y_range_frame = Frame(self.dual_y_range_frame)
        left_y_range_frame.pack(fill='x', pady=1, padx=5)
        Label(left_y_range_frame, text="左Y轴:", fg='blue').pack(side='left')
        Label(left_y_range_frame, text="最小值:").pack(side='left', padx=(10,0))
        self.left_y_min_var = StringVar(value="auto")
        Entry(left_y_range_frame, textvariable=self.left_y_min_var, width=8).pack(side='left', padx=2)
        Label(left_y_range_frame, text="最大值:").pack(side='left', padx=(5,0))
        self.left_y_max_var = StringVar(value="auto")
        Entry(left_y_range_frame, textvariable=self.left_y_max_var, width=8).pack(side='left', padx=2)
        ttk.Button(left_y_range_frame, text="应用左Y轴", command=self.apply_axis_limits).pack(side='left', padx=5)
        
        # 右Y轴范围
        right_y_range_frame = Frame(self.dual_y_range_frame)
        right_y_range_frame.pack(fill='x', pady=1, padx=5)
        Label(right_y_range_frame, text="右Y轴:", fg='red').pack(side='left')
        Label(right_y_range_frame, text="最小值:").pack(side='left', padx=(10,0))
        self.right_y_min_var = StringVar(value="auto")
        Entry(right_y_range_frame, textvariable=self.right_y_min_var, width=8).pack(side='left', padx=2)
        Label(right_y_range_frame, text="最大值:").pack(side='left', padx=(5,0))
        self.right_y_max_var = StringVar(value="auto")
        Entry(right_y_range_frame, textvariable=self.right_y_max_var, width=8).pack(side='left', padx=2)
        ttk.Button(right_y_range_frame, text="应用右Y轴", command=self.apply_axis_limits).pack(side='left', padx=5)
        
        # 范围重置按钮
        reset_range_frame = Frame(range_frame)
        reset_range_frame.pack(fill='x', pady=5, padx=5)
        ttk.Button(reset_range_frame, text="🔄 重置所有范围为自动", 
                  command=self.reset_axis_limits).pack(side='left', padx=5)
        ttk.Button(reset_range_frame, text="📊 应用所有设置", 
                  command=self.apply_axis_limits).pack(side='left', padx=5)
        
        # 坐标轴单位缩放设置
        unit_scale_frame = LabelFrame(self.scrollable_frame, text="📏 坐标轴单位缩放")
        unit_scale_frame.pack(pady=10, fill='x')
        
        # 说明文字
        help_frame = Frame(unit_scale_frame)
        help_frame.pack(fill='x', pady=2, padx=5)
        Label(help_frame, text="说明: 数据将除以缩放因子显示 (例: 1e-5, 1e-6)", 
              font=('Arial', 8), fg='gray').pack(anchor='w')
        
        # X轴单位缩放
        x_unit_frame = Frame(unit_scale_frame)
        x_unit_frame.pack(fill='x', pady=2, padx=5)
        Label(x_unit_frame, text="X轴缩放因子:").pack(side='left')
        Entry(x_unit_frame, textvariable=self.x_unit_scale_var, width=10).pack(side='left', padx=5)
        
        # 单Y轴单位缩放
        self.single_y_unit_frame = Frame(unit_scale_frame)
        self.single_y_unit_frame.pack(fill='x', pady=2, padx=5)
        Label(self.single_y_unit_frame, text="Y轴缩放因子:").pack(side='left')
        Entry(self.single_y_unit_frame, textvariable=self.y_unit_scale_var, width=10).pack(side='left', padx=5)
        
        # 双Y轴单位缩放
        self.dual_y_unit_frame = Frame(unit_scale_frame)
        
        # 左Y轴单位缩放
        left_y_unit_frame = Frame(self.dual_y_unit_frame)
        left_y_unit_frame.pack(fill='x', pady=1, padx=5)
        Label(left_y_unit_frame, text="左Y轴缩放因子:", fg='blue').pack(side='left')
        Entry(left_y_unit_frame, textvariable=self.left_y_unit_scale_var, width=10).pack(side='left', padx=5)
        
        # 右Y轴单位缩放
        right_y_unit_frame = Frame(self.dual_y_unit_frame)
        right_y_unit_frame.pack(fill='x', pady=1, padx=5)
        Label(right_y_unit_frame, text="右Y轴缩放因子:", fg='red').pack(side='left')
        Entry(right_y_unit_frame, textvariable=self.right_y_unit_scale_var, width=10).pack(side='left', padx=5)
        
        # 应用按钮
        unit_apply_frame = Frame(unit_scale_frame)
        unit_apply_frame.pack(fill='x', pady=5, padx=5)
        ttk.Button(unit_apply_frame, text="📊 应用单位缩放", 
                  command=self.apply_unit_scale).pack(side='left', padx=5)
        ttk.Button(unit_apply_frame, text="🔄 重置为1", 
                  command=self.reset_unit_scale).pack(side='left', padx=5)
        
        # 图例文字编辑
        legend_text_frame = LabelFrame(self.scrollable_frame, text="✏️ 图例文字编辑")
        legend_text_frame.pack(pady=10, fill='x')
        
        # 存储图例文字自定义
        self.legend_text_vars = {}
        
        # 为每个变量创建图例文字编辑框
        for var in self.variables:
            var_legend_frame = Frame(legend_text_frame)
            var_legend_frame.pack(fill='x', pady=2, padx=5)
            
            Label(var_legend_frame, text=f"{var}:", width=15, anchor='w').pack(side='left')
            self.legend_text_vars[var] = StringVar(value=var)  # 默认使用变量名
            legend_text_entry = Entry(var_legend_frame, textvariable=self.legend_text_vars[var], width=30)
            legend_text_entry.pack(side='left', padx=5, fill='x', expand=True)
        
        ttk.Button(legend_text_frame, text="应用图例文字", command=self.redraw_plot).pack(fill='x', pady=5, padx=5)
        
        # 设置管理
        settings_frame = LabelFrame(self.scrollable_frame, text="⚙️ 设置管理")
        settings_frame.pack(pady=10, fill='x')
        
        settings_buttons_frame = Frame(settings_frame)
        settings_buttons_frame.pack(fill='x', pady=5, padx=5)
        
        ttk.Button(settings_buttons_frame, text="💾 保存设置", 
                  command=self.save_settings).pack(side='left', padx=5, fill='x', expand=True)
        ttk.Button(settings_buttons_frame, text="📂 加载设置", 
                  command=self.load_settings).pack(side='left', padx=5, fill='x', expand=True)
        ttk.Button(settings_buttons_frame, text="🔄 重置设置", 
                  command=self.reset_settings).pack(side='left', padx=5, fill='x', expand=True)
        
        # 图表控制
        control_frame = LabelFrame(self.scrollable_frame, text="🎯 图表控制")
        control_frame.pack(pady=10, fill='x')
        
        ttk.Button(control_frame, text="🔄 刷新图表", 
                  command=self.redraw_plot).pack(fill='x', pady=2, padx=5)
        ttk.Button(control_frame, text="💾 保存图片", 
                  command=self.save_figure).pack(fill='x', pady=2, padx=5)
        
        # 初始化控件显示状态
        self.update_dual_y_display()

    def create_variable_control_group(self, parent, var):
        """为单个变量创建控制组"""
        var_frame = LabelFrame(parent, text=f"📊 {var}")
        var_frame.pack(fill='x', pady=5, padx=5)
        
        # 第一行：显示、Y轴分配
        row1 = Frame(var_frame)
        row1.pack(fill='x', pady=2, padx=5)
        
        # 显示复选框
        Checkbutton(row1, text="显示", 
                   variable=self.variable_states[var]['visible'],
                   command=self.redraw_plot).pack(side='left', padx=5)
        
        # Y轴选择
        if self.dual_y_enabled.get():
            Label(row1, text="Y轴:").pack(side='left', padx=5)
            y_axis_combo = ttk.Combobox(row1, textvariable=self.variable_states[var]['y_axis'], 
                                       values=['left', 'right'], width=8, state='readonly')
            y_axis_combo.pack(side='left', padx=5)
            y_axis_combo.bind('<<ComboboxSelected>>', lambda e: self.redraw_plot())
        
        # 第二行：点大小、透明度
        row2 = Frame(var_frame)
        row2.pack(fill='x', pady=2, padx=5)
        
        Label(row2, text="点大小:").pack(side='left')
        Entry(row2, textvariable=self.variable_states[var]['point_size'], width=5).pack(side='left', padx=2)
        
        Label(row2, text="透明度:").pack(side='left', padx=(10,0))
        Entry(row2, textvariable=self.variable_states[var]['alpha'], width=5).pack(side='left', padx=2)
        
        Label(row2, text="样式:").pack(side='left', padx=(10,0))
        marker_combo = ttk.Combobox(row2, textvariable=self.variable_states[var]['marker'], 
                                   values=['o', 's', '^', 'v', '<', '>', 'D', 'd', 'p', '*', 'h', 'H', '+', 'x'], 
                                   width=3, state='readonly')
        marker_combo.pack(side='left', padx=2)
        marker_combo.bind('<<ComboboxSelected>>', lambda e: self.redraw_plot())
        
        # Smooth曲线线宽设置（仅在启用Smooth化时有用）
        Label(row2, text="曲线线宽:").pack(side='left', padx=(10,0))
        Entry(row2, textvariable=self.variable_states[var]['line_width'], width=5).pack(side='left', padx=2)
        
        # 第三行：Smooth化设置
        row3 = Frame(var_frame)
        row3.pack(fill='x', pady=2, padx=5)
        
        Checkbutton(row3, text="Smooth化", 
                   variable=self.smoothing_states[var],
                   command=self.redraw_plot).pack(side='left', padx=5)
        
        Label(row3, text="方法:").pack(side='left', padx=5)
        method_combo = ttk.Combobox(row3, textvariable=self.smoothing_params[var]['method'], 
                                   values=['gaussian', 'spline', 'moving_average'], 
                                   width=10, state='readonly')
        method_combo.pack(side='left', padx=2)
        method_combo.bind('<<ComboboxSelected>>', lambda e: self.redraw_plot())
        
        Label(row3, text="参数:").pack(side='left', padx=5)
        Entry(row3, textvariable=self.smoothing_params[var]['sigma'], width=6).pack(side='left', padx=2)

    def update_dual_y_display(self):
        """更新双Y轴控件的显示状态"""
        if self.dual_y_enabled.get():
            # 启用双Y轴：显示双Y轴控件，隐藏单Y轴控件
            self.single_ylabel_frame.pack_forget()
            self.dual_ylabel_frame.pack(fill='x')
            
            self.single_y_axis_frame.pack_forget()
            self.dual_y_axis_frame.pack(fill='x', pady=2, padx=5)
            
            # 图例位置控件
            self.single_legend_frame.pack_forget()
            self.dual_legend_frame.pack(fill='x', pady=2, padx=5)
            
            # 坐标轴范围控件
            self.single_y_range_frame.pack_forget()
            self.dual_y_range_frame.pack(fill='x', pady=2, padx=5)
            
            # 坐标轴单位缩放控件
            self.single_y_unit_frame.pack_forget()
            self.dual_y_unit_frame.pack(fill='x', pady=2, padx=5)
        else:
            # 禁用双Y轴：显示单Y轴控件，隐藏双Y轴控件
            self.dual_ylabel_frame.pack_forget()
            self.single_ylabel_frame.pack(fill='x')
            
            self.dual_y_axis_frame.pack_forget()
            self.single_y_axis_frame.pack(fill='x', pady=2, padx=5)
            
            # 图例位置控件
            self.dual_legend_frame.pack_forget()
            self.single_legend_frame.pack(fill='x', pady=2, padx=5)
            
            # 坐标轴范围控件
            self.dual_y_range_frame.pack_forget()
            self.single_y_range_frame.pack(fill='x', pady=2, padx=5)
            
            # 坐标轴单位缩放控件
            self.dual_y_unit_frame.pack_forget()
            self.single_y_unit_frame.pack(fill='x', pady=2, padx=5)

    def toggle_dual_y_axis(self):
        """切换双Y轴模式"""
        # 更新控件显示状态
        self.update_dual_y_display()
        
        # 重建变量控制组以显示/隐藏Y轴选择
        # 找到变量设置框架并重建
        for child in self.scrollable_frame.winfo_children():
            if isinstance(child, LabelFrame) and child.cget('text') == '📈 变量设置':
                # 清空现有控件
                for widget in child.winfo_children():
                    widget.destroy()
                # 重新创建变量控制组
                for var in self.variables:
                    self.create_variable_control_group(child, var)
                break
        
        self.redraw_plot()
    
    def show_all_data_sources(self):
        """显示所有数据源"""
        for folder_name in self.data_source_visibility:
            self.data_source_visibility[folder_name].set(True)
        self.redraw_plot()
    
    def hide_all_data_sources(self):
        """隐藏所有数据源"""
        for folder_name in self.data_source_visibility:
            self.data_source_visibility[folder_name].set(False)
        self.redraw_plot()
    
    def enable_all_smoothing(self):
        """启用所有变量的Smooth化"""
        for var in self.variables:
            self.smoothing_states[var].set(True)
        self.redraw_plot()
    
    def disable_all_smoothing(self):
        """禁用所有变量的Smooth化"""
        for var in self.variables:
            self.smoothing_states[var].set(False)
        self.redraw_plot()

    def apply_smoothing(self, x_data, y_data, var):
        """应用Smooth化处理"""
        if not self.smoothing_states[var].get():
            return x_data, y_data
        
        method = self.smoothing_params[var]['method'].get()
        
        try:
            if method == 'gaussian':
                # 高斯滤波
                sigma = float(self.smoothing_params[var]['sigma'].get())
                smooth_y = gaussian_filter1d(y_data, sigma=sigma)
                return x_data, smooth_y
            
            elif method == 'spline':
                # 样条插值
                num_points = int(self.smoothing_params[var]['points'].get())
                # 创建样条插值函数
                tck = interpolate.splrep(x_data, y_data, s=0)
                # 生成更密集的x点
                x_smooth = np.linspace(x_data.min(), x_data.max(), num_points)
                y_smooth = interpolate.splev(x_smooth, tck)
                return x_smooth, y_smooth
            
            elif method == 'moving_average':
                # 移动平均
                window = int(float(self.smoothing_params[var]['sigma'].get()))
                if window < 1:
                    window = 1
                # 计算移动平均
                smooth_y = np.convolve(y_data, np.ones(window)/window, mode='same')
                return x_data, smooth_y
        
        except (ValueError, TypeError):
            # 参数错误时返回原始数据
            return x_data, y_data
        
        return x_data, y_data

    def redraw_plot(self):
        """重新绘制图表"""
        # 清除当前图表
        self.ax1.clear()
        if self.ax2:
            self.ax2.clear()
        
        # 获取单位缩放因子
        try:
            x_unit_scale = float(self.x_unit_scale_var.get())
            if x_unit_scale == 0:
                x_unit_scale = 1.0
        except (ValueError, ZeroDivisionError):
            x_unit_scale = 1.0
        
        try:
            y_unit_scale = float(self.y_unit_scale_var.get())
            if y_unit_scale == 0:
                y_unit_scale = 1.0
        except (ValueError, ZeroDivisionError):
            y_unit_scale = 1.0
        
        try:
            left_y_unit_scale = float(self.left_y_unit_scale_var.get())
            if left_y_unit_scale == 0:
                left_y_unit_scale = 1.0
        except (ValueError, ZeroDivisionError):
            left_y_unit_scale = 1.0
        
        try:
            right_y_unit_scale = float(self.right_y_unit_scale_var.get())
            if right_y_unit_scale == 0:
                right_y_unit_scale = 1.0
        except (ValueError, ZeroDivisionError):
            right_y_unit_scale = 1.0
        
        # 创建双Y轴（如果需要）
        if self.dual_y_enabled.get() and self.ax2 is None:
            self.ax2 = self.ax1.twinx()
        elif not self.dual_y_enabled.get():
            if self.ax2:
                self.ax2.remove()
                self.ax2 = None
        
        # 设置双Y轴的样式
        if self.dual_y_enabled.get() and self.ax2:
            # 设置右Y轴的位置，稍微向右偏移
            self.ax2.spines['right'].set_position(('outward', 10))
            # 隐藏不需要的spine
            self.ax2.spines['top'].set_visible(False)
            self.ax1.spines['top'].set_visible(False)
            self.ax1.spines['right'].set_visible(False)
        
        # 准备图例信息
        left_lines = []
        left_labels = []
        right_lines = []
        right_labels = []
        
        # 绘制每个变量
        color_index = 0
        for var in self.variables:
            if not self.variable_states[var]['visible'].get():
                continue
            
            # 确定使用哪个Y轴
            if self.dual_y_enabled.get() and self.variable_states[var]['y_axis'].get() == 'right':
                current_ax = self.ax2
            else:
                current_ax = self.ax1
            
            if current_ax is None:
                continue
            
            # 获取样式参数
            try:
                point_size = int(self.variable_states[var]['point_size'].get())
                line_width = float(self.variable_states[var]['line_width'].get())
                alpha = float(self.variable_states[var]['alpha'].get())
                marker = self.variable_states[var]['marker'].get()
            except ValueError:
                point_size = 10
                line_width = 1.5
                alpha = 0.7
                marker = 'o'
            
            # 绘制每个数据源的数据
            for dataset_index, (folder_name, df) in enumerate(self.data_list):
                # 检查数据源是否可见
                if not self.data_source_visibility[folder_name].get():
                    continue
                
                if var not in df.columns:
                    continue
                
                # 为每个数据集+变量组合分配唯一颜色
                unique_color_index = color_index + dataset_index
                color = self.color_list[unique_color_index % len(self.color_list)]
                
                x_data = df['time'].values
                y_data = df[var].values
                
                # 应用单位缩放
                x_data_scaled = x_data / x_unit_scale
                
                # 根据Y轴分配应用相应的单位缩放
                if self.dual_y_enabled.get() and current_ax == self.ax2:
                    # 右Y轴
                    y_data_scaled = y_data / right_y_unit_scale
                elif self.dual_y_enabled.get():
                    # 左Y轴
                    y_data_scaled = y_data / left_y_unit_scale
                else:
                    # 单Y轴
                    y_data_scaled = y_data / y_unit_scale
                
                # 应用Smooth化
                x_smooth, y_smooth = self.apply_smoothing(x_data_scaled, y_data_scaled, var)
                
                # 使用自定义图例文字
                custom_var_name = self.legend_text_vars[var].get() if var in self.legend_text_vars else var
                label = f"{custom_var_name} ({folder_name})"
                
                # 只使用散点图，不绘制连接线
                if self.smoothing_states[var].get():
                    # Smooth化模式：绘制原始数据点和Smooth曲线
                    # 原始数据点（较淡）
                    scatter_orig = current_ax.scatter(x_data_scaled, y_data_scaled, s=point_size*0.6, alpha=alpha*0.4, 
                                                    color=color, edgecolors='none', marker=marker, label=f"{label} (原始)")
                    
                    # 生成Smooth曲线的颜色（更深的版本）
                    import matplotlib.colors as mcolors
                    rgb = mcolors.to_rgb(color)
                    smooth_color = tuple(max(0, c * 0.6) for c in rgb)
                    
                    # 绘制Smooth曲线
                    line = current_ax.plot(x_smooth, y_smooth, label=f"{label} (Smooth)", 
                                         linewidth=line_width*1.2, alpha=alpha, color=smooth_color)[0]
                else:
                    # 只绘制散点图
                    scatter = current_ax.scatter(x_data_scaled, y_data_scaled, s=point_size, alpha=alpha, 
                                              color=color, edgecolors='none', marker=marker, label=label)
                    line = scatter  # 将scatter赋值给line以便后续处理
                
                # 添加到对应的图例列表
                if current_ax == self.ax1:
                    left_lines.append(line)
                    if self.smoothing_states[var].get():
                        left_labels.append(f"{label} (Smooth)")
                    else:
                        left_labels.append(label)
                else:
                    right_lines.append(line)
                    if self.smoothing_states[var].get():
                        right_labels.append(f"{label} (Smooth)")
                    else:
                        right_labels.append(label)
            
            # 更新颜色索引，为下一个变量准备新的颜色基础
            color_index += len(self.data_list)
        
        # 获取字体大小设置
        try:
            title_fontsize = int(self.title_fontsize_var.get())
            axis_fontsize = int(self.axis_fontsize_var.get())
            tick_fontsize = int(self.tick_fontsize_var.get())
            legend_fontsize = int(self.legend_fontsize_var.get())
        except ValueError:
            title_fontsize = 14
            axis_fontsize = 12
            tick_fontsize = 10
            legend_fontsize = 10
        
        try:
            legend_scale = float(self.legend_marker_scale_var.get())
        except ValueError:
            legend_scale = 1.0
        
        # 设置轴标签
        xlabel_text = self.xlabel_var.get()
        self.ax1.set_xlabel(xlabel_text, fontsize=axis_fontsize, fontfamily='Times New Roman')
        
        if self.dual_y_enabled.get() and self.ax2:
            # 双Y轴标签
            left_vars = [var for var in self.variables 
                        if self.variable_states[var]['visible'].get() 
                        and self.variable_states[var]['y_axis'].get() == 'left']
            right_vars = [var for var in self.variables 
                         if self.variable_states[var]['visible'].get() 
                         and self.variable_states[var]['y_axis'].get() == 'right']
            
            if left_vars:
                ylabel_left = self.left_ylabel_var.get()
                self.ax1.set_ylabel(ylabel_left, fontsize=axis_fontsize, fontfamily='Times New Roman')
                # 设置左Y轴颜色和位置
                self.ax1.tick_params(axis='y', labelcolor='blue', labelsize=tick_fontsize)
                self.ax1.yaxis.label.set_color('blue')
                # 调整左Y轴标题位置，向左偏移
                self.ax1.yaxis.set_label_coords(-0.1, 0.5)
                
            if right_vars:
                ylabel_right = self.right_ylabel_var.get()
                self.ax2.set_ylabel(ylabel_right, fontsize=axis_fontsize, fontfamily='Times New Roman')
                # 设置右Y轴颜色和位置
                self.ax2.tick_params(axis='y', labelcolor='red', labelsize=tick_fontsize)
                self.ax2.yaxis.label.set_color('red')
                # 调整右Y轴标题位置，向右偏移
                self.ax2.yaxis.set_label_coords(1.1, 0.5)
        else:
            # 单Y轴标签
            ylabel_text = self.ylabel_var.get()
            self.ax1.set_ylabel(ylabel_text, fontsize=axis_fontsize, fontfamily='Times New Roman')
        
        # 设置X轴刻度字体大小
        self.ax1.tick_params(axis='x', labelsize=tick_fontsize)
        if not self.dual_y_enabled.get():
            self.ax1.tick_params(axis='y', labelsize=tick_fontsize)
        
        # 设置标题
        title_text = self.title_var.get()
        self.ax1.set_title(title_text, fontsize=title_fontsize, fontfamily='Times New Roman')
        
        # 创建图例
        if self.dual_y_enabled.get() and self.ax2:
            # 双Y轴模式：使用用户自定义的位置
            left_legend_pos = self.left_legend_position_var.get()
            right_legend_pos = self.right_legend_position_var.get()
            
            if left_lines:
                legend1 = self.ax1.legend(left_lines, left_labels, loc=left_legend_pos, 
                                        prop={'size': legend_fontsize, 'family': 'Times New Roman'}, markerscale=legend_scale)
                
            if right_lines:
                legend2 = self.ax2.legend(right_lines, right_labels, loc=right_legend_pos, 
                                        prop={'size': legend_fontsize, 'family': 'Times New Roman'}, markerscale=legend_scale)
        else:
            # 单Y轴模式：使用统一的图例位置设置
            legend_position = self.legend_position_var.get()
            if left_lines:
                legend1 = self.ax1.legend(left_lines, left_labels, loc=legend_position, 
                                        prop={'size': legend_fontsize, 'family': 'Times New Roman'}, markerscale=legend_scale)
        
        # 应用网格和样式
        self.ax1.grid(True, alpha=0.3)
        self.ax1.set_axisbelow(True)
        
        # 调整布局以避免Y轴标题重叠
        if self.dual_y_enabled.get() and self.ax2:
            # 双Y轴模式需要更多的左右边距
            self.fig.tight_layout(pad=2.0)
            # 手动调整子图位置，为Y轴标题留出空间
            self.fig.subplots_adjust(left=0.12, right=0.88, top=0.92, bottom=0.1)
        else:
            # 单Y轴模式使用标准布局
            self.fig.tight_layout(pad=1.5)
        
        # 重新应用用户设置的坐标轴范围（在布局调整后）
        self._reapply_axis_limits()
        
        self.canvas.draw()

    def save_figure(self):
        """保存图片"""
        try:
            # 设置默认文件名
            var_names = "_".join(self.variables[:3])  # 最多取前3个变量名
            if len(self.variables) > 3:
                var_names += f"_等{len(self.variables)}个变量"
            default_filename = f"多变量图表_{var_names}.png"
            
            # 打开文件保存对话框
            save_path = filedialog.asksaveasfilename(
                title="保存图片",
                defaultextension=".png",
                initialfile=default_filename,
                filetypes=[
                    ("PNG图片", "*.png"),
                    ("JPEG图片", "*.jpg"),
                    ("PDF文档", "*.pdf"),
                    ("SVG矢量图", "*.svg"),
                    ("所有文件", "*.*")
                ]
            )
            
            if save_path:
                # 根据文件扩展名确定保存格式
                file_ext = Path(save_path).suffix.lower()
                if file_ext == '.jpg' or file_ext == '.jpeg':
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='jpeg')
                elif file_ext == '.pdf':
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='pdf')
                elif file_ext == '.svg':
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='svg')
                else:
                    self.fig.savefig(save_path, dpi=300, bbox_inches='tight', format='png')
                
                print(f"图片已保存至: {save_path}")
                messagebox.showinfo("保存成功", f"图片已成功保存至:\n{save_path}")
            else:
                print("保存图片已取消")
        except Exception as e:
            print(f"保存图片时出错: {e}")
            messagebox.showerror("保存失败", f"保存图片时出错:\n{str(e)}")

    def apply_unit_scale(self):
        """应用坐标轴单位缩放"""
        self.redraw_plot()
    
    def reset_unit_scale(self):
        """重置坐标轴单位缩放为1"""
        self.x_unit_scale_var.set("1")
        self.y_unit_scale_var.set("1")
        self.left_y_unit_scale_var.set("1")
        self.right_y_unit_scale_var.set("1")
        self.redraw_plot()
    
    def toggle_log_scale(self, axis):
        """切换对数坐标"""
        if axis == 'x':
            # X轴对数切换（同时影响左右Y轴）
            if self.ax1.get_xscale() == 'log':
                self.ax1.set_xscale('linear')
                if self.ax2:
                    self.ax2.set_xscale('linear')
            else:
                self.ax1.set_xscale('log')
                if self.ax2:
                    self.ax2.set_xscale('log')
                    
        elif axis == 'y':
            # 单Y轴对数切换
            if self.ax1.get_yscale() == 'log':
                self.ax1.set_yscale('linear')
            else:
                self.ax1.set_yscale('log')
                
        elif axis == 'left_y':
            # 左Y轴对数切换
            if self.ax1.get_yscale() == 'log':
                self.ax1.set_yscale('linear')
            else:
                self.ax1.set_yscale('log')
                
        elif axis == 'right_y':
            # 右Y轴对数切换
            if self.ax2 and self.ax2.get_yscale() == 'log':
                self.ax2.set_yscale('linear')
            elif self.ax2:
                self.ax2.set_yscale('log')
        
        # 切换坐标轴类型后，重新应用用户设置的范围
        self._reapply_axis_limits()
        self.canvas.draw()

    def _reapply_axis_limits(self):
        """内部方法：重新应用用户设置的坐标轴范围（不触发重绘）"""
        try:
            # 应用X轴范围
            x_min_str = self.x_min_var.get().strip()
            x_max_str = self.x_max_var.get().strip()
            
            if x_min_str.lower() != 'auto' or x_max_str.lower() != 'auto':
                try:
                    x_min = float(x_min_str) if x_min_str.lower() != 'auto' else None
                    x_max = float(x_max_str) if x_max_str.lower() != 'auto' else None
                    
                    if x_min is not None and x_max is not None and x_min < x_max:
                        self.ax1.set_xlim(x_min, x_max)
                        if self.ax2:
                            self.ax2.set_xlim(x_min, x_max)
                    elif x_min is not None:
                        current_xlim = self.ax1.get_xlim()
                        self.ax1.set_xlim(x_min, current_xlim[1])
                        if self.ax2:
                            self.ax2.set_xlim(x_min, current_xlim[1])
                    elif x_max is not None:
                        current_xlim = self.ax1.get_xlim()
                        self.ax1.set_xlim(current_xlim[0], x_max)
                        if self.ax2:
                            self.ax2.set_xlim(current_xlim[0], x_max)
                except (ValueError, TypeError):
                    pass
            
            # 应用Y轴范围
            if self.dual_y_enabled.get() and self.ax2:
                # 双Y轴模式
                # 左Y轴范围
                left_y_min_str = self.left_y_min_var.get().strip()
                left_y_max_str = self.left_y_max_var.get().strip()
                
                if left_y_min_str.lower() != 'auto' or left_y_max_str.lower() != 'auto':
                    try:
                        y_min = float(left_y_min_str) if left_y_min_str.lower() != 'auto' else None
                        y_max = float(left_y_max_str) if left_y_max_str.lower() != 'auto' else None
                        
                        if y_min is not None and y_max is not None and y_min < y_max:
                            self.ax1.set_ylim(y_min, y_max)
                        elif y_min is not None:
                            current_ylim = self.ax1.get_ylim()
                            self.ax1.set_ylim(y_min, current_ylim[1])
                        elif y_max is not None:
                            current_ylim = self.ax1.get_ylim()
                            self.ax1.set_ylim(current_ylim[0], y_max)
                    except (ValueError, TypeError):
                        pass
                
                # 右Y轴范围
                right_y_min_str = self.right_y_min_var.get().strip()
                right_y_max_str = self.right_y_max_var.get().strip()
                
                if right_y_min_str.lower() != 'auto' or right_y_max_str.lower() != 'auto':
                    try:
                        y_min = float(right_y_min_str) if right_y_min_str.lower() != 'auto' else None
                        y_max = float(right_y_max_str) if right_y_max_str.lower() != 'auto' else None
                        
                        if y_min is not None and y_max is not None and y_min < y_max:
                            self.ax2.set_ylim(y_min, y_max)
                        elif y_min is not None:
                            current_ylim = self.ax2.get_ylim()
                            self.ax2.set_ylim(y_min, current_ylim[1])
                        elif y_max is not None:
                            current_ylim = self.ax2.get_ylim()
                            self.ax2.set_ylim(current_ylim[0], y_max)
                    except (ValueError, TypeError):
                        pass
            else:
                # 单Y轴模式
                y_min_str = self.y_min_var.get().strip()
                y_max_str = self.y_max_var.get().strip()
                
                if y_min_str.lower() != 'auto' or y_max_str.lower() != 'auto':
                    try:
                        y_min = float(y_min_str) if y_min_str.lower() != 'auto' else None
                        y_max = float(y_max_str) if y_max_str.lower() != 'auto' else None
                        
                        if y_min is not None and y_max is not None and y_min < y_max:
                            self.ax1.set_ylim(y_min, y_max)
                        elif y_min is not None:
                            current_ylim = self.ax1.get_ylim()
                            self.ax1.set_ylim(y_min, current_ylim[1])
                        elif y_max is not None:
                            current_ylim = self.ax1.get_ylim()
                            self.ax1.set_ylim(current_ylim[0], y_max)
                    except (ValueError, TypeError):
                        pass
        except Exception:
            pass  # 静默失败，不影响绘图
    
    def apply_axis_limits(self):
        """应用坐标轴范围设置（公开方法，会触发重绘）"""
        try:
            # 调用内部方法应用范围
            self._reapply_axis_limits()
            # 刷新显示
            self.canvas.draw()
        except Exception as e:
            messagebox.showerror("范围设置错误", f"应用坐标轴范围时出错:\n{str(e)}")

    def reset_axis_limits(self):
        """重置所有坐标轴范围为自动"""
        try:
            # 重置所有范围变量为auto
            self.x_min_var.set("auto")
            self.x_max_var.set("auto")
            self.y_min_var.set("auto")
            self.y_max_var.set("auto")
            self.left_y_min_var.set("auto")
            self.left_y_max_var.set("auto")
            self.right_y_min_var.set("auto")
            self.right_y_max_var.set("auto")
            
            # 重新绘制图表，让matplotlib自动计算范围
            self.redraw_plot()
            
            messagebox.showinfo("重置完成", "所有坐标轴范围已重置为自动")
            
        except Exception as e:
            messagebox.showerror("重置失败", f"重置坐标轴范围时出错:\n{str(e)}")

    def save_settings(self):
        """保存当前设置到文件"""
        import json
        from tkinter import filedialog
        
        try:
            # 收集所有设置
            settings = {
                # 基本设置
                'dual_y_enabled': self.dual_y_enabled.get(),
                'title': self.title_var.get(),
                'xlabel': self.xlabel_var.get(),
                'ylabel': self.ylabel_var.get() if hasattr(self, 'ylabel_var') else '',
                'left_ylabel': self.left_ylabel_var.get(),
                'right_ylabel': self.right_ylabel_var.get(),
                
                # 字体设置
                'title_fontsize': self.title_fontsize_var.get(),
                'axis_fontsize': self.axis_fontsize_var.get(),
                'tick_fontsize': self.tick_fontsize_var.get(),
                'legend_fontsize': self.legend_fontsize_var.get(),
                'legend_marker_scale': self.legend_marker_scale_var.get(),
                
                # 图例位置
                'legend_position': self.legend_position_var.get(),
                'left_legend_position': self.left_legend_position_var.get(),
                'right_legend_position': self.right_legend_position_var.get(),
                
                # 坐标轴范围
                'x_min': self.x_min_var.get(),
                'x_max': self.x_max_var.get(),
                'y_min': self.y_min_var.get(),
                'y_max': self.y_max_var.get(),
                'left_y_min': self.left_y_min_var.get(),
                'left_y_max': self.left_y_max_var.get(),
                'right_y_min': self.right_y_min_var.get(),
                'right_y_max': self.right_y_max_var.get(),
                
                # 坐标轴单位缩放
                'x_unit_scale': self.x_unit_scale_var.get(),
                'y_unit_scale': self.y_unit_scale_var.get(),
                'left_y_unit_scale': self.left_y_unit_scale_var.get(),
                'right_y_unit_scale': self.right_y_unit_scale_var.get(),
                
                # 变量设置
                'variable_states': {},
                'smoothing_states': {},
                'smoothing_params': {},
                'legend_texts': {},
                'data_source_visibility': {}
            }
            
            # 保存数据源可见性状态
            for folder_name in self.data_source_visibility:
                settings['data_source_visibility'][folder_name] = self.data_source_visibility[folder_name].get()
            
            # 保存每个变量的设置
            for var in self.variables:
                settings['variable_states'][var] = {
                    'visible': self.variable_states[var]['visible'].get(),
                    'y_axis': self.variable_states[var]['y_axis'].get(),
                    'point_size': self.variable_states[var]['point_size'].get(),
                    'line_width': self.variable_states[var]['line_width'].get(),
                    'alpha': self.variable_states[var]['alpha'].get(),
                    'marker': self.variable_states[var]['marker'].get()
                }
                settings['smoothing_states'][var] = self.smoothing_states[var].get()
                settings['smoothing_params'][var] = {
                    'method': self.smoothing_params[var]['method'].get(),
                    'sigma': self.smoothing_params[var]['sigma'].get(),
                    'points': self.smoothing_params[var]['points'].get()
                }
                if var in self.legend_text_vars:
                    settings['legend_texts'][var] = self.legend_text_vars[var].get()
            
            # 选择保存位置
            file_path = filedialog.asksaveasfilename(
                title="保存设置文件",
                defaultextension=".json",
                filetypes=[("JSON文件", "*.json"), ("所有文件", "*.*")]
            )
            
            if file_path:
                with open(file_path, 'w', encoding='utf-8') as f:
                    json.dump(settings, f, ensure_ascii=False, indent=2)
                messagebox.showinfo("保存成功", f"设置已保存到:\n{file_path}")
                
        except Exception as e:
            messagebox.showerror("保存失败", f"保存设置时出错:\n{str(e)}")

    def load_settings(self):
        """从文件加载设置"""
        import json
        from tkinter import filedialog
        
        try:
            # 选择设置文件
            file_path = filedialog.askopenfilename(
                title="加载设置文件",
                filetypes=[("JSON文件", "*.json"), ("所有文件", "*.*")]
            )
            
            if not file_path:
                return
                
            with open(file_path, 'r', encoding='utf-8') as f:
                settings = json.load(f)
            
            # 恢复基本设置
            if 'dual_y_enabled' in settings:
                self.dual_y_enabled.set(settings['dual_y_enabled'])
            if 'title' in settings:
                self.title_var.set(settings['title'])
            if 'xlabel' in settings:
                self.xlabel_var.set(settings['xlabel'])
            if 'ylabel' in settings and hasattr(self, 'ylabel_var'):
                self.ylabel_var.set(settings['ylabel'])
            if 'left_ylabel' in settings:
                self.left_ylabel_var.set(settings['left_ylabel'])
            if 'right_ylabel' in settings:
                self.right_ylabel_var.set(settings['right_ylabel'])
                
            # 恢复字体设置
            font_settings = ['title_fontsize', 'axis_fontsize', 'tick_fontsize', 'legend_fontsize']
            for font_setting in font_settings:
                if font_setting in settings:
                    getattr(self, f'{font_setting}_var').set(settings[font_setting])
            
            if 'legend_marker_scale' in settings:
                self.legend_marker_scale_var.set(settings['legend_marker_scale'])
            
            # 恢复图例位置
            legend_settings = ['legend_position', 'left_legend_position', 'right_legend_position']
            for legend_setting in legend_settings:
                if legend_setting in settings:
                    getattr(self, f'{legend_setting}_var').set(settings[legend_setting])
            
            # 恢复坐标轴范围
            axis_range_settings = ['x_min', 'x_max', 'y_min', 'y_max', 
                                 'left_y_min', 'left_y_max', 'right_y_min', 'right_y_max']
            for range_setting in axis_range_settings:
                if range_setting in settings:
                    getattr(self, f'{range_setting}_var').set(settings[range_setting])
            
            # 恢复坐标轴单位缩放
            unit_scale_settings = ['x_unit_scale', 'y_unit_scale', 'left_y_unit_scale', 'right_y_unit_scale']
            for scale_setting in unit_scale_settings:
                if scale_setting in settings:
                    getattr(self, f'{scale_setting}_var').set(settings[scale_setting])
            
            # 恢复变量设置
            if 'variable_states' in settings:
                for var, var_settings in settings['variable_states'].items():
                    if var in self.variable_states:
                        for key, value in var_settings.items():
                            if key in self.variable_states[var]:
                                self.variable_states[var][key].set(value)
            
            if 'smoothing_states' in settings:
                for var, state in settings['smoothing_states'].items():
                    if var in self.smoothing_states:
                        self.smoothing_states[var].set(state)
            
            if 'smoothing_params' in settings:
                for var, params in settings['smoothing_params'].items():
                    if var in self.smoothing_params:
                        for key, value in params.items():
                            if key in self.smoothing_params[var]:
                                self.smoothing_params[var][key].set(value)
            
            if 'legend_texts' in settings:
                for var, text in settings['legend_texts'].items():
                    if var in self.legend_text_vars:
                        self.legend_text_vars[var].set(text)
            
            # 恢复数据源可见性状态
            if 'data_source_visibility' in settings:
                for folder_name, visible in settings['data_source_visibility'].items():
                    if folder_name in self.data_source_visibility:
                        self.data_source_visibility[folder_name].set(visible)
            
            # 更新控件显示状态
            self.update_dual_y_display()
            
            # 重新绘制图表
            self.redraw_plot()
            
            messagebox.showinfo("加载成功", f"设置已从以下文件加载:\n{file_path}")
            
        except Exception as e:
            messagebox.showerror("加载失败", f"加载设置时出错:\n{str(e)}")

    def reset_settings(self):
        """重置所有设置为默认值"""
        try:
            # 重置基本设置
            self.dual_y_enabled.set(False)
            self.title_var.set(f"多变量图表: {', '.join(self.variables)}")
            self.xlabel_var.set("时间")
            if hasattr(self, 'ylabel_var'):
                self.ylabel_var.set("数值")
            self.left_ylabel_var.set("左Y轴数值")
            self.right_ylabel_var.set("右Y轴数值")
            
            # 重置字体设置
            self.title_fontsize_var.set("14")
            self.axis_fontsize_var.set("12")
            self.tick_fontsize_var.set("10")
            self.legend_fontsize_var.set("10")
            self.legend_marker_scale_var.set("1.0")
            
            # 重置图例位置
            self.legend_position_var.set("best")
            self.left_legend_position_var.set("upper left")
            self.right_legend_position_var.set("upper right")
            
            # 重置坐标轴范围
            self.x_min_var.set("auto")
            self.x_max_var.set("auto")
            self.y_min_var.set("auto")
            self.y_max_var.set("auto")
            self.left_y_min_var.set("auto")
            self.left_y_max_var.set("auto")
            self.right_y_min_var.set("auto")
            self.right_y_max_var.set("auto")
            
            # 重置坐标轴单位缩放
            self.x_unit_scale_var.set("1")
            self.y_unit_scale_var.set("1")
            self.left_y_unit_scale_var.set("1")
            self.right_y_unit_scale_var.set("1")
            
            # 重置变量设置
            for var in self.variables:
                self.variable_states[var]['visible'].set(True)
                self.variable_states[var]['y_axis'].set('left')
                self.variable_states[var]['point_size'].set("10")
                self.variable_states[var]['line_width'].set("1.5")
                self.variable_states[var]['alpha'].set("0.7")
                self.variable_states[var]['marker'].set('o')
                
                self.smoothing_states[var].set(False)
                self.smoothing_params[var]['method'].set('gaussian')
                self.smoothing_params[var]['sigma'].set("1.0")
                self.smoothing_params[var]['points'].set("100")
                
                if var in self.legend_text_vars:
                    self.legend_text_vars[var].set(var)
            
            # 重置数据源可见性状态
            for folder_name in self.data_source_visibility:
                self.data_source_visibility[folder_name].set(True)
            
            # 更新控件显示状态
            self.update_dual_y_display()
            
            # 重新绘制图表
            self.redraw_plot()
            
            messagebox.showinfo("重置完成", "所有设置已重置为默认值")
            
        except Exception as e:
            messagebox.showerror("重置失败", f"重置设置时出错:\n{str(e)}")

    def on_closing(self):
        """处理窗口关闭事件"""
        self.is_closed = True
        plt.close(self.fig)
        self.root.destroy()

def main():
    """主函数 - 启动变量选择器"""
    selector = VariableSelector()
    selector.run()

if __name__ == "__main__":
    main()