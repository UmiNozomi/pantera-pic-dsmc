#!/usr/bin/env python3
"""
Interactive Hα Spectral Analysis Tool
=====================================

Reads collision event data (CSV) and provides an interactive GUI to:
1. Explore spectra by reaction channel
2. Rotate line-of-sight (LOS) direction in real-time
3. Adjust binning resolution
4. Toggle channel visibility

Usage:
    python process_ha_events.py [path/to/csv_files]
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, CheckButtons, Button
from pathlib import Path
import sys
import glob

# Constants
C_LIGHT = 299792458.0  # m/s
LAMBDA_HALPHA = 656.28e-9  # m
NM_TO_M = 1e-9

# Channel Definitions
CHANNEL_NAMES = {
    7:  'H⁺ + H₂',
    13: 'H₂⁺ + H₂',
    18: 'H₃⁺ + H₂',
    22: 'H + H₂',
    26: 'H₂ + H₂'
}

CHANNEL_COLORS = {
    7:  '#e74c3c', # Red
    13: '#3498db', # Blue
    18: '#2ecc71', # Green
    22: '#f39c12', # Orange
    26: '#9b59b6'  # Purple
}

from scipy.signal import savgol_filter
from scipy.interpolate import interp1d
from scipy.optimize import minimize
import tkinter as tk
from tkinter import filedialog

class SpectralAnalyzer:
    def __init__(self, df):
        self.df = df
        self.fig, self.ax = plt.subplots(figsize=(12, 8))
        plt.subplots_adjust(left=0.1, bottom=0.35, right=0.8, top=0.95) # Reduce right margin for controls
        
        # Initial settings
        self.theta_min = 0.0   # Min polar angle (deg)
        self.theta_max = 180.0 # Max polar angle (deg)
        self.phi = 0.0         # Azimuthal angle (deg)
        self.n_angles = 4     # Number of angles to average
        self.nbins = 1600
        self.visible_channels = list(CHANNEL_NAMES.keys())
        self.show_total = True
        self.do_symmetrize = False
        self.do_smooth = True
        self.lock_y = False       # Lock Y-axis
        self.y_max_value = 1.0    # Manual Y-axis max
        
        # Data filtering
        self.reaction_ids = sorted(df['reaction_id'].unique())
        
        # Experimental data storage
        self.exp_wavelengths = None
        self.exp_intensities = None
        self.exp_loaded = False
        
        # Channel weights for fitting (scale factors)
        self.channel_weights = {rid: 1.0 for rid in self.reaction_ids}
        
        # Setup GUI
        self.setup_widgets()
        self.update_plot(None)
        
    def setup_widgets(self):
        # Sliders Area configuration
        ax_color = 'lightgoldenrodyellow'
        
        # 1. Theta Min Slider (0-180)
        ax_theta_min = plt.axes([0.1, 0.25, 0.35, 0.03], facecolor=ax_color)
        self.s_theta_min = Slider(ax_theta_min, 'θ min (°)', 0.0, 360.0, valinit=self.theta_min)
        self.s_theta_min.on_changed(self.update_plot)
        
        # 2. Theta Max Slider (0-180)
        ax_theta_max = plt.axes([0.1, 0.20, 0.35, 0.03], facecolor=ax_color)
        self.s_theta_max = Slider(ax_theta_max, 'θ max (°)', 0.0, 360.0, valinit=self.theta_max)
        self.s_theta_max.on_changed(self.update_plot)
        
        # 3. Phi Slider (0-360)
        ax_phi = plt.axes([0.1, 0.15, 0.35, 0.03], facecolor=ax_color)
        self.s_phi = Slider(ax_phi, 'Azimuth φ (°)', 0.0, 360.0, valinit=self.phi)
        self.s_phi.on_changed(self.update_plot)
        
        # 4. N Angles Slider
        ax_nangles = plt.axes([0.55, 0.25, 0.15, 0.03], facecolor=ax_color)
        self.s_nangles = Slider(ax_nangles, 'N θ', 1, 36, valinit=self.n_angles, valstep=1)
        self.s_nangles.on_changed(self.update_plot)
        
        # 5. Binning Slider
        ax_bins = plt.axes([0.1, 0.10, 0.35, 0.03], facecolor=ax_color)
        self.s_bins = Slider(ax_bins, 'Bins', 50, 5000, valinit=self.nbins, valstep=10)
        self.s_bins.on_changed(self.update_plot)
        
        # 6. Smooth Slider
        ax_smooth = plt.axes([0.1, 0.05, 0.35, 0.03], facecolor=ax_color)
        self.s_smooth = Slider(ax_smooth, 'Smooth Win', 2, 51, valinit=15, valstep=2)
        self.s_smooth.on_changed(self.update_plot)

        # 5. Checkboxes for Channels (Right Side)
        ax_check = plt.axes([0.82, 0.4, 0.15, 0.3], frameon=False)
        self.labels = [CHANNEL_NAMES.get(rid, f'R{int(rid)}') for rid in self.reaction_ids]
        self.labels.append('Total')
        
        actives = [True] * len(self.labels)
        self.check = CheckButtons(ax_check, self.labels, actives)
        self.check.on_clicked(self.toggle_visibility)
        ax_check.set_title("Channels")

        # 6. Feature Toggles (Symmetrize / Smooth / Lock Y)
        ax_features = plt.axes([0.82, 0.2, 0.15, 0.18], frameon=False)
        self.feature_labels = ['Symmetrize', 'Smooth', 'Lock Y']
        self.feature_actives = [False, False, False]
        self.check_features = CheckButtons(ax_features, self.feature_labels, self.feature_actives)
        self.check_features.on_clicked(self.toggle_features)
        ax_features.set_title("Processing")
        
        # Y Max Slider (only used when Lock Y is on)
        ax_ymax = plt.axes([0.55, 0.20, 0.15, 0.03], facecolor=ax_color)
        self.s_ymax = Slider(ax_ymax, 'Y Max', 0.1, 10.0, valinit=1.0, valstep=0.1)
        self.s_ymax.on_changed(self.update_plot)
        
        # 7. Reset Button
        ax_reset = plt.axes([0.55, 0.05, 0.08, 0.04])
        self.b_reset = Button(ax_reset, 'Reset', hovercolor='0.975')
        self.b_reset.on_clicked(self.reset)
        
        # 8. Load Experiment Button
        ax_load = plt.axes([0.55, 0.10, 0.08, 0.04])
        self.b_load = Button(ax_load, 'Load Exp', hovercolor='lightgreen')
        self.b_load.on_clicked(self.load_experiment)
        
        # 9. Auto Fit Button
        ax_fit = plt.axes([0.55, 0.15, 0.08, 0.04])
        self.b_fit = Button(ax_fit, 'Auto Fit', hovercolor='lightyellow')
        self.b_fit.on_clicked(self.auto_fit)

        # Info text
        self.txt_info = self.fig.text(0.1, 0.02, '', fontsize=10)

    def calculate_vlos(self, theta_deg, phi_deg):
        # Convert to radians
        theta = np.radians(theta_deg)
        phi = np.radians(phi_deg)
        
        # Cartesian components of LOS vector
        nx = np.sin(theta) * np.cos(phi)
        ny = np.sin(theta) * np.sin(phi)
        nz = np.cos(theta)
        
        # Calculate vLOS for all particles
        vlos = (self.df['vx'] * nx + 
                self.df['vy'] * ny + 
                self.df['vz'] * nz)
        return vlos

    def apply_symmetrization(self, wavelengths, intensities, center_wavelength=656.28):
        """Symmetrize spectrum around center wavelength."""
        # Create interpolation function
        interp_func = interp1d(wavelengths, intensities, 
                               kind='linear', 
                               bounds_error=False, 
                               fill_value=0.0)
        
        symmetrized = np.zeros_like(intensities)
        for i, wl in enumerate(wavelengths):
            delta = wl - center_wavelength
            mirror_wl = center_wavelength - delta
            # Average intensity at wl and its mirror point
            mirror_intensity = interp_func(mirror_wl)
            symmetrized[i] = (intensities[i] + mirror_intensity) / 2.0
            
        return symmetrized

    def update_plot(self, val):
        self.ax.clear()
        
        theta_min = self.s_theta_min.val
        theta_max = self.s_theta_max.val
        phi = self.s_phi.val
        n_angles = int(self.s_nangles.val)
        bins = int(self.s_bins.val)
        smooth_win = int(self.s_smooth.val)
        if smooth_win % 2 == 0: smooth_win += 1 # Ensure odd
        
        # Ensure theta_min <= theta_max
        if theta_min > theta_max:
            theta_min, theta_max = theta_max, theta_min

        # Update Info Text
        self.txt_info.set_text(f"LOS: θ=[{theta_min:.1f}°,{theta_max:.1f}°], φ={phi:.1f}°, N={n_angles} | Events: {len(self.df)}")
        
        # Generate angle array for averaging
        if n_angles == 1:
            theta_array = [(theta_min + theta_max) / 2]
        else:
            theta_array = np.linspace(theta_min, theta_max, n_angles)
        
        # Histogram parameters
        wl_min, wl_max = 634.0, 711.0
        hist_range = (wl_min, wl_max)
        
        # Combine data for plotting
        total_hist = np.zeros(bins)
        bin_edges = None
        
        max_y = 0
        
        # Check active status
        status = self.check.get_status()
        
        # Loop through data channels
        for i, rid in enumerate(self.reaction_ids):
            if not status[i]: continue # Skip if unchecked
            
            mask = self.df['reaction_id'] == rid
            if not mask.any(): continue
            
            weights = self.df.loc[mask, 'weight'].values
            
            # Accumulate histogram over all angles
            channel_hist = np.zeros(bins)
            
            for theta in theta_array:
                vlos = self.calculate_vlos(theta, phi)
                wavelengths_nm = LAMBDA_HALPHA * 1e9 * (1.0 + vlos / C_LIGHT)
                wl_data = wavelengths_nm[mask]
                
                hist, edges = np.histogram(wl_data, bins=bins, range=hist_range, weights=weights)
                channel_hist += hist
            
            centers = 0.5 * (edges[1:] + edges[:-1])
            bin_edges = edges
            
            # Normalize by number of angles
            channel_hist /= n_angles
            
            # Apply Processing
            if self.do_symmetrize:
                channel_hist = self.apply_symmetrization(centers, channel_hist)
                
            if self.do_smooth and len(channel_hist) > smooth_win:
                try:
                    channel_hist = savgol_filter(channel_hist, window_length=smooth_win, polyorder=3)
                    channel_hist[channel_hist < 0] = 0
                except:
                    pass
            
            # Add to plot
            label = CHANNEL_NAMES.get(rid, f'R{rid}')
            color = CHANNEL_COLORS.get(rid, 'gray')
            
            self.ax.plot(centers, channel_hist, label=label, color=color, alpha=0.8, linewidth=1.5)
            self.ax.fill_between(centers, channel_hist, color=color, alpha=0.1)
            
            total_hist += channel_hist
            max_y = max(max_y, channel_hist.max())

        # Use total_hist max for normalization (so Total peaks at 1.0)
        norm_factor = total_hist.max() if total_hist.max() > 0 else 1.0
        
        # Normalize all simulation data for display
        if norm_factor > 0:
            # Re-plot with normalized data
            self.ax.clear()
            
            for i, rid in enumerate(self.reaction_ids):
                if not status[i]: continue
                
                mask = self.df['reaction_id'] == rid
                if not mask.any(): continue
                
                weights = self.df.loc[mask, 'weight'].values
                channel_hist = np.zeros(bins)
                
                for theta in theta_array:
                    vlos = self.calculate_vlos(theta, phi)
                    wavelengths_nm = LAMBDA_HALPHA * 1e9 * (1.0 + vlos / C_LIGHT)
                    wl_data = wavelengths_nm[mask]
                    hist, edges = np.histogram(wl_data, bins=bins, range=hist_range, weights=weights)
                    channel_hist += hist
                
                centers = 0.5 * (edges[1:] + edges[:-1])
                channel_hist /= n_angles
                
                if self.do_symmetrize:
                    channel_hist = self.apply_symmetrization(centers, channel_hist)
                if self.do_smooth and len(channel_hist) > smooth_win:
                    try:
                        channel_hist = savgol_filter(channel_hist, window_length=smooth_win, polyorder=3)
                        channel_hist[channel_hist < 0] = 0
                    except:
                        pass
                
                # Normalize by total max
                channel_hist_norm = channel_hist / norm_factor
                
                label = CHANNEL_NAMES.get(rid, f'R{rid}')
                color = CHANNEL_COLORS.get(rid, 'gray')
                self.ax.plot(centers, channel_hist_norm, label=label, color=color, alpha=0.8, linewidth=1.5)
                self.ax.fill_between(centers, channel_hist_norm, color=color, alpha=0.1)
            
            # Plot normalized total (should peak at 1.0)
            if status[-1] and bin_edges is not None:
                centers = 0.5 * (bin_edges[1:] + bin_edges[:-1])
                total_hist_norm = total_hist / norm_factor
                # Apply smoothing to total if enabled
                if self.do_smooth and len(total_hist_norm) > smooth_win:
                    try:
                        total_hist_norm = savgol_filter(total_hist_norm, window_length=smooth_win, polyorder=3)
                        total_hist_norm[total_hist_norm < 0] = 0
                        # Renormalize after smoothing so peak = 1.0
                        if total_hist_norm.max() > 0:
                            total_hist_norm = total_hist_norm / total_hist_norm.max()
                    except:
                        pass
                self.ax.plot(centers, total_hist_norm, 'k--', label='Total', linewidth=2)

        # Plot experimental data if loaded (same axis, both normalized)
        if self.exp_loaded and self.exp_wavelengths is not None:
            self.ax.plot(self.exp_wavelengths, self.exp_intensities, 'ro-', 
                        label='Experiment', markersize=3, linewidth=1.5, alpha=0.8)

        # Decorate Plot
        self.ax.set_xlabel('Wavelength (nm)', fontsize=12)
        self.ax.set_ylabel('Normalized Intensity', fontsize=12)
        self.ax.set_title(f'Interactive Hα Spectrum Analysis (θ: {theta_min:.0f}°-{theta_max:.0f}°, N={n_angles})', fontsize=14)
        self.ax.axvline(656.28, color='gray', linestyle=':', label='Hα Rest')
        self.ax.set_xlim(wl_min, wl_max)
        
        # Y-axis: use locked value or auto-scale to 1.1
        if self.lock_y:
            self.ax.set_ylim(0, self.s_ymax.val)
        else:
            self.ax.set_ylim(0, 1.1)
        self.ax.legend(loc='upper right')
        self.ax.grid(True, alpha=0.3)
        
        self.fig.canvas.draw_idle()

    def toggle_visibility(self, label):
        self.update_plot(None)
        
    def toggle_features(self, label):
        index = self.feature_labels.index(label)
        # CheckButtons toggles state internally before callback, so get new status
        status = self.check_features.get_status()
        if label == 'Symmetrize':
            self.do_symmetrize = status[index]
        elif label == 'Smooth':
            self.do_smooth = status[index]
        elif label == 'Lock Y':
            self.lock_y = status[index]
        self.update_plot(None)

    def reset(self, event):
        self.s_theta_min.reset()
        self.s_theta_max.reset()
        self.s_phi.reset()
        self.s_nangles.reset()
        self.s_bins.reset()
        self.s_smooth.reset()

    def load_experiment(self, event):
        """Load experimental data from Excel file."""
        # Hide matplotlib window temporarily for file dialog
        root = tk.Tk()
        root.withdraw()
        
        filepath = filedialog.askopenfilename(
            title="Select Experimental Data",
            filetypes=[
                ("Excel files", "*.xlsx *.xls"),
                ("CSV files", "*.csv"),
                ("All files", "*.*")
            ]
        )
        root.destroy()
        
        if not filepath:
            return
            
        try:
            # Load data
            if filepath.endswith('.csv'):
                exp_df = pd.read_csv(filepath)
            else:
                exp_df = pd.read_excel(filepath)
            
            # Assume first column is wavelength, second is intensity
            cols = exp_df.columns.tolist()
            self.exp_wavelengths = exp_df[cols[0]].values
            self.exp_intensities = exp_df[cols[1]].values
            
            # Normalize experimental data
            self.exp_intensities = self.exp_intensities / np.max(self.exp_intensities)
            
            self.exp_loaded = True
            print(f"Loaded experimental data: {len(self.exp_wavelengths)} points from {filepath}")
            self.update_plot(None)
            
        except Exception as e:
            print(f"Error loading file: {e}")
            self.exp_loaded = False

    def compute_total_spectrum(self, theta_min, theta_max, phi, n_angles, bins, wl_min, wl_max):
        """Compute total spectrum for given parameters (for fitting)."""
        if n_angles == 1:
            theta_array = [(theta_min + theta_max) / 2]
        else:
            theta_array = np.linspace(theta_min, theta_max, n_angles)
        
        hist_range = (wl_min, wl_max)
        total_hist = np.zeros(bins)
        
        status = self.check.get_status()
        
        for i, rid in enumerate(self.reaction_ids):
            if not status[i]: continue
            
            mask = self.df['reaction_id'] == rid
            if not mask.any(): continue
            
            weights = self.df.loc[mask, 'weight'].values
            channel_hist = np.zeros(bins)
            
            for theta in theta_array:
                vlos = self.calculate_vlos(theta, phi)
                wavelengths_nm = LAMBDA_HALPHA * 1e9 * (1.0 + vlos / C_LIGHT)
                wl_data = wavelengths_nm[mask]
                hist, edges = np.histogram(wl_data, bins=bins, range=hist_range, weights=weights)
                channel_hist += hist
            
            channel_hist /= n_angles
            total_hist += channel_hist * self.channel_weights.get(rid, 1.0)
        
        centers = 0.5 * (edges[1:] + edges[:-1]) if 'edges' in dir() else np.linspace(wl_min, wl_max, bins)
        return centers, total_hist

    def auto_fit(self, event):
        """Automatically fit simulation to experimental data."""
        if not self.exp_loaded:
            print("Please load experimental data first!")
            return
        
        print("Starting auto-fit optimization...")
        
        # Get current wavelength range
        wl_min, wl_max = 652.0, 660.0
        bins = int(self.s_bins.val)
        
        # Interpolate experimental data to our bin centers
        sim_centers = np.linspace(wl_min, wl_max, bins)
        exp_interp = interp1d(self.exp_wavelengths, self.exp_intensities, 
                              kind='linear', bounds_error=False, fill_value=0.0)
        exp_target = exp_interp(sim_centers)
        
        # Optimization function
        def objective(params):
            theta_min, theta_max, phi, n_angles = params
            n_angles = max(1, int(n_angles))
            
            # Ensure valid theta range
            if theta_min > theta_max:
                theta_min, theta_max = theta_max, theta_min
            theta_min = np.clip(theta_min, 0, 180)
            theta_max = np.clip(theta_max, 0, 180)
            phi = phi % 360
            
            centers, sim_spec = self.compute_total_spectrum(
                theta_min, theta_max, phi, n_angles, bins, wl_min, wl_max)
            
            # Normalize simulation
            if np.max(sim_spec) > 0:
                sim_spec = sim_spec / np.max(sim_spec)
            
            # Mean squared error
            mse = np.mean((sim_spec - exp_target) ** 2)
            return mse
        
        # Initial guess from current slider values
        x0 = [
            self.s_theta_min.val,
            self.s_theta_max.val,
            self.s_phi.val,
            self.s_nangles.val
        ]
        
        # Bounds
        bounds = [
            (0, 180),    # theta_min
            (0, 180),    # theta_max  
            (0, 360),    # phi
            (1, 36)      # n_angles
        ]
        
        # Run optimization
        result = minimize(objective, x0, method='L-BFGS-B', bounds=bounds,
                         options={'maxiter': 100, 'disp': True})
        
        if result.success:
            theta_min, theta_max, phi, n_angles = result.x
            print(f"Fit complete! θ=[{theta_min:.1f}°, {theta_max:.1f}°], φ={phi:.1f}°, N={int(n_angles)}")
            print(f"Final MSE: {result.fun:.6f}")
            
            # Update sliders with fitted values
            self.s_theta_min.set_val(theta_min)
            self.s_theta_max.set_val(theta_max)
            self.s_phi.set_val(phi)
            self.s_nangles.set_val(int(n_angles))
        else:
            print(f"Optimization failed: {result.message}")


def main():
    # 1. Find Files
    if len(sys.argv) > 1:
        arg_path = sys.argv[1]
        if '*' in arg_path or '?' in arg_path:
            file_list = glob.glob(arg_path)
        elif Path(arg_path).is_dir():
            file_list = glob.glob(str(Path(arg_path) / 'ha_events*.csv'))
        else:
            file_list = [arg_path]
    else:
        possible_patterns = [
            'iech2/results/ha_events*.csv',
            'results/ha_events*.csv',
            '../results/ha_events*.csv'
        ]
        file_list = []
        for pat in possible_patterns:
            found = glob.glob(pat)
            if found:
                file_list = found
                break
    
    if not file_list:
        print("❌ Err: No 'ha_events*.csv' files found.")
        return

    print(f"📂 Found {len(file_list)} files. Loading...")
    
    # 2. Load Data
    dfs = []
    total_events = 0
    for f in file_list:
        try:
            df_chunk = pd.read_csv(f)
            dfs.append(df_chunk)
            total_events += len(df_chunk)
            # print(f"  + Loaded {len(df_chunk)} events from {Path(f).name}")
        except Exception as e:
            print(f"  ⚠️ Failed to load {f}: {e}")
    
    if not dfs:
        return

    df = pd.concat(dfs, ignore_index=True)
    
    # Filter out invalid reaction IDs (NaN)
    df = df.dropna(subset=['reaction_id'])
    df['reaction_id'] = df['reaction_id'].astype(int)
    
    print(f"✅ Loaded {len(df)} total events.")
    print(f"   Reactions found: {sorted(df['reaction_id'].unique())}")

    # 3. Launch Interactive Plot
    print("\n🚀 Launching Interactive Analyzer...")
    analyzer = SpectralAnalyzer(df)
    plt.show()

if __name__ == '__main__':
    main()
