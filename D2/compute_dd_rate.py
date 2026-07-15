#!/usr/bin/env python3
"""
compute_dd_rate.py - D-D Fusion Reaction Rate Post-Processor
=============================================================

Computes the volumetric D-D fusion reaction rate density R [reactions/(m³·s)]
from PANTERA PIC-DSMC VTK flowfield output using:

    R = n_beam(x) * n_target * <σv>(T(x))

where <σv>(T) is the D-D fusion reactivity from Bosch & Hale (1992),
Nuclear Fusion 32, 611.

For MCC background mode:
  - n_beam  : tracked species number density from VTK (e.g. D+, D, D2+, D2)
  - n_target: uniform MCC background density (from input file)
  - T       : ion/atom translational temperature from VTK

Usage:
  python compute_dd_rate.py results/dsmc_flowfield_10000.vtk
  python compute_dd_rate.py results/dsmc_flowfield_*.vtk --bg-density 1.5e20
  python compute_dd_rate.py results/dsmc_flowfield_10000.vtk --list-fields

Author: Generated for PANTERA PIC-DSMC project
"""

import numpy as np
import struct
import os
import sys
import glob
import re
import argparse
from pathlib import Path

# ==================== Physical Constants ====================
KB  = 1.380649e-23       # Boltzmann constant [J/K]
QE  = 1.602176634e-19    # Elementary charge [C]
AMU = 1.66053906660e-27  # Atomic mass unit [kg]


# ==================== Bosch-Hale D-D Fusion Reactivity ====================

def bosch_hale_sigmav(T_keV, reaction='dd_total'):
    """
    Compute D-D fusion reactivity <σv> [m³/s] using Bosch-Hale (1992).

    Parameters
    ----------
    T_keV : array_like
        Ion temperature in keV (1 keV ≈ 1.16e7 K)
    reaction : str
        'dd_n_he3' : D(d,n)³He branch
        'dd_p_t'   : D(d,p)T branch
        'dd_total' : sum of both branches (default)

    Returns
    -------
    sigmav : ndarray
        Reactivity <σv> in m³/s
    """
    T = np.atleast_1d(np.float64(T_keV))
    result = np.zeros_like(T)

    # Coefficients for D(d,n)³He
    if reaction in ('dd_n_he3', 'dd_total'):
        BG   = 31.3970
        mrc2 = 937814.0  # keV
        C = [5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0]

        mask = (T > 0.5) & (T < 550)
        Tm = T[mask]
        if len(Tm) > 0:
            theta = Tm / (1.0 - Tm * (C[1] + Tm * (C[3] + Tm * C[5]))
                          / (1.0 + Tm * (C[2] + Tm * (C[4] + Tm * C[6]))))
            xi = (BG**2 / (4.0 * theta))**(1.0 / 3.0)
            sv = C[0] * theta * np.sqrt(xi / (mrc2 * Tm**3)) * np.exp(-3.0 * xi)
            result[mask] += sv  # cm³/s

    # Coefficients for D(d,p)T
    if reaction in ('dd_p_t', 'dd_total'):
        BG   = 31.3970
        mrc2 = 937814.0
        C = [5.65718e-12, 3.41267e-3, 1.99167e-3, 0.0, 1.05060e-5, 0.0, 0.0]

        mask = (T > 0.5) & (T < 550)
        Tm = T[mask]
        if len(Tm) > 0:
            theta = Tm / (1.0 - Tm * (C[1] + Tm * (C[3] + Tm * C[5]))
                          / (1.0 + Tm * (C[2] + Tm * (C[4] + Tm * C[6]))))
            xi = (BG**2 / (4.0 * theta))**(1.0 / 3.0)
            sv = C[0] * theta * np.sqrt(xi / (mrc2 * Tm**3)) * np.exp(-3.0 * xi)
            result[mask] += sv

    # Convert cm³/s → m³/s
    result *= 1e-6
    return result


# ==================== VTK File Reader ====================

class VTKReader:
    """Reader for PANTERA Legacy VTK output files (ASCII and Binary)."""

    def __init__(self, filename):
        self.filename = filename
        self.fields = {}          # name -> (ncomp, ntuples, dtype_str, data)
        self.point_fields = {}
        self.grid_type = None
        self.dimensions = None
        self.x_coords = None
        self.y_coords = None
        self.z_coords = None
        self.points = None
        self.cells = None
        self.ncells = 0
        self.is_binary = False
        self._read()

    # ---- public API ----

    def get_field(self, name):
        """Return field data array by exact name, or None."""
        if name in self.fields:
            return self.fields[name][3]
        return None

    def find_field(self, *patterns):
        """Find the first field whose name contains ALL patterns (case-insensitive)."""
        for name in self.fields:
            nl = name.lower().strip()
            if all(p.lower() in nl for p in patterns):
                return self.fields[name][3], name
        return None, None

    def list_fields(self):
        return list(self.fields.keys())

    def get_cell_centers(self):
        if self.grid_type == 'RECTILINEAR_GRID' and self.x_coords is not None:
            nx = self.dimensions[0] - 1
            ny = max(self.dimensions[1] - 1, 1)
            cx = 0.5 * (self.x_coords[:-1] + self.x_coords[1:])
            cy = 0.5 * (self.y_coords[:-1] + self.y_coords[1:]) if ny > 1 else np.array([0.0])
            centers_x = np.zeros(self.ncells)
            centers_y = np.zeros(self.ncells)
            for j in range(ny):
                for i in range(nx):
                    idx = j * nx + i
                    if idx < self.ncells:
                        centers_x[idx] = cx[i]
                        centers_y[idx] = cy[j] if ny > 1 else 0.0
            return centers_x, centers_y
        elif self.points is not None and self.cells is not None:
            cx = np.zeros(self.ncells)
            cy = np.zeros(self.ncells)
            for ic, cell in enumerate(self.cells):
                if ic < self.ncells:
                    pts = self.points[cell]
                    cx[ic] = np.mean(pts[:, 0])
                    cy[ic] = np.mean(pts[:, 1])
            return cx, cy
        return None, None

    # ---- internal parsing ----

    def _read(self):
        with open(self.filename, 'rb') as f:
            f.readline()  # # vtk DataFile Version 3.0
            f.readline()  # title
            fmt_line = f.readline().decode('ascii', errors='ignore').strip()
            self.is_binary = 'BINARY' in fmt_line.upper()
            if self.is_binary:
                self._read_binary(f)
            else:
                content = f.read().decode('ascii', errors='ignore')
                self._read_ascii(content)

    def _read_ascii(self, content):
        lines = content.split('\n')
        i = 0
        target = self.fields  # switch between cell/point fields
        while i < len(lines):
            line = lines[i].strip()

            if line.startswith('DATASET'):
                self.grid_type = line.split()[1]
            elif line.startswith('DIMENSIONS'):
                self.dimensions = [int(x) for x in line.split()[1:4]]
            elif line.startswith('X_COORDINATES'):
                n = int(line.split()[1]); i += 1
                vals = self._collect_values(lines, i, n, float)
                self.x_coords = np.array(vals); i += self._lines_consumed; continue
            elif line.startswith('Y_COORDINATES'):
                n = int(line.split()[1]); i += 1
                vals = self._collect_values(lines, i, n, float)
                self.y_coords = np.array(vals); i += self._lines_consumed; continue
            elif line.startswith('Z_COORDINATES'):
                n = int(line.split()[1]); i += 1
                vals = self._collect_values(lines, i, n, float)
                self.z_coords = np.array(vals); i += self._lines_consumed; continue
            elif line.startswith('POINTS'):
                npts = int(line.split()[1]); i += 1
                vals = self._collect_values(lines, i, npts * 3, float)
                self.points = np.array(vals).reshape(npts, 3)
                i += self._lines_consumed; continue
            elif line.startswith('CELLS') and not line.startswith('CELL_'):
                nc = int(line.split()[1]); self.ncells = nc; i += 1
                cells = []
                for _ in range(nc):
                    parts = [int(x) for x in lines[i].strip().split()]
                    cells.append(parts[1:]); i += 1
                self.cells = cells; continue
            elif line.startswith('CELL_TYPES'):
                nc = int(line.split()[1]); i += 1
                for _ in range(nc): i += 1
                continue
            elif line.startswith('CELL_DATA'):
                self.ncells = int(line.split()[1])
                target = self.fields
            elif line.startswith('POINT_DATA'):
                target = self.point_fields
            elif line.startswith('FIELD'):
                nf = int(line.split()[2])
                for _ in range(nf):
                    i += 1
                    while i < len(lines) and lines[i].strip() == '': i += 1
                    if i >= len(lines): break
                    hdr = lines[i].strip().split()
                    # The field name may contain spaces from Fortran output
                    # Last 3 tokens are: ncomp ntuples dtype
                    if len(hdr) < 4: continue
                    dtype_s = hdr[-1]; nt = int(hdr[-2]); nc = int(hdr[-3])
                    fname = ' '.join(hdr[:-3]).strip()
                    total = nc * nt; i += 1
                    conv = int if dtype_s == 'integer' else float
                    vals = self._collect_values(lines, i, total, conv)
                    data = np.array(vals[:total])
                    target[fname] = (nc, nt, dtype_s, data)
                    i += self._lines_consumed
                continue
            i += 1

    def _collect_values(self, lines, start, n, conv):
        vals = []; lc = 0
        i = start
        while len(vals) < n and i < len(lines):
            for tok in lines[i].strip().split():
                try: vals.append(conv(tok))
                except ValueError: pass
            i += 1; lc += 1
        self._lines_consumed = lc
        return vals

    @staticmethod
    def _consume_newline(f):
        """Consume a trailing newline if present, otherwise seek back."""
        pos = f.tell()
        byte = f.read(1)
        if byte != b'\n':
            f.seek(pos)

    def _read_binary(self, f):
        target = self.fields
        while True:
            raw = f.readline()
            if not raw: break
            line = raw.decode('ascii', errors='ignore').strip()
            if not line: continue

            if line.startswith('DATASET'):
                self.grid_type = line.split()[1]
            elif line.startswith('DIMENSIONS'):
                self.dimensions = [int(x) for x in line.split()[1:4]]
            elif line.startswith('X_COORDINATES'):
                n = int(line.split()[1])
                self.x_coords = np.frombuffer(f.read(n * 8), dtype='>f8').copy()
                self._consume_newline(f)
            elif line.startswith('Y_COORDINATES'):
                n = int(line.split()[1])
                self.y_coords = np.frombuffer(f.read(n * 8), dtype='>f8').copy()
                self._consume_newline(f)
            elif line.startswith('Z_COORDINATES'):
                n = int(line.split()[1])
                self.z_coords = np.frombuffer(f.read(n * 8), dtype='>f8').copy()
                self._consume_newline(f)
            elif line.startswith('POINTS'):
                npts = int(line.split()[1])
                self.points = np.frombuffer(f.read(npts * 3 * 8), dtype='>f8').copy().reshape(npts, 3)
                self._consume_newline(f)
            elif line.startswith('CELLS') and not line.startswith('CELL_'):
                parts = line.split()
                nc = int(parts[1]); ti = int(parts[2]); self.ncells = nc
                raw_data = np.frombuffer(f.read(ti * 4), dtype='>i4')
                cells = []; idx = 0
                for _ in range(nc):
                    nn = raw_data[idx]; cells.append(raw_data[idx+1:idx+1+nn].tolist()); idx += 1 + nn
                self.cells = cells
                self._consume_newline(f)
            elif line.startswith('CELL_TYPES'):
                nc = int(line.split()[1]); f.read(nc * 4)
                self._consume_newline(f)
            elif line.startswith('CELL_DATA'):
                self.ncells = int(line.split()[1]); target = self.fields
            elif line.startswith('POINT_DATA'):
                target = self.point_fields
            elif line.startswith('FIELD'):
                nf = int(line.split()[2])
                for _ in range(nf):
                    fl = f.readline().decode('ascii', errors='ignore').strip()
                    while not fl:
                        fl = f.readline().decode('ascii', errors='ignore').strip()
                    hdr = fl.split()
                    if len(hdr) < 4: continue
                    dtype_s = hdr[-1]; nt = int(hdr[-2]); nc = int(hdr[-3])
                    fname = ' '.join(hdr[:-3]).strip()
                    total = nc * nt
                    if dtype_s == 'double':
                        data = np.frombuffer(f.read(total * 8), dtype='>f8').copy()
                    elif dtype_s == 'integer':
                        data = np.frombuffer(f.read(total * 4), dtype='>i4').copy()
                    else:
                        data = np.frombuffer(f.read(total * 4), dtype='>f4').copy()
                    self._consume_newline(f)
                    target[fname] = (nc, nt, dtype_s, data)

        # Infer ncells from field data if parser missed CELL_DATA header
        if self.ncells == 0 and self.fields:
            for name, (nc, nt, dt, data) in self.fields.items():
                if nt > 0:
                    self.ncells = nt
                    break


# ==================== VTK Writer ====================

def _sanitize_vtk_name(name):
    """Clean field name: remove spaces, replace special chars for VTK compatibility."""
    s = name.strip()
    # Collapse internal spaces (Fortran list-directed output adds spaces)
    s = '_'.join(s.split())
    # Replace characters that VTK parsers can't handle
    s = s.replace('+', 'plus').replace('-', 'minus')
    return s


def _write_vtk_data(f, data, ncells, vals_per_line=6):
    """Write data array in VTK-compatible format (multiple values per line)."""
    n = min(len(data), ncells)
    for i in range(n):
        f.write(f'{data[i]:.10e}')
        if (i + 1) % vals_per_line == 0 or i == n - 1:
            f.write('\n')
        else:
            f.write(' ')
    # Pad with zeros if data is shorter than ncells
    for i in range(n, ncells):
        f.write('0.0000000000e+00')
        if (i + 1) % vals_per_line == 0 or i == ncells - 1:
            f.write('\n')
        else:
            f.write(' ')


def write_vtk_ascii(reader, new_fields, output_file):
    """Write new VTK file with original grid + selected original fields + new computed fields."""

    # Determine consistent ncells from geometry
    ncells = reader.ncells
    if reader.cells is not None:
        ncells = len(reader.cells)
    reader.ncells = ncells  # ensure consistency

    with open(output_file, 'w') as f:
        f.write('# vtk DataFile Version 3.0\n')
        f.write('D-D Fusion Rate Post-Processing\n')
        f.write('ASCII\n')

        if reader.grid_type == 'RECTILINEAR_GRID':
            f.write('DATASET RECTILINEAR_GRID\n')
            f.write(f'DIMENSIONS {reader.dimensions[0]} {reader.dimensions[1]} {reader.dimensions[2]}\n')
            f.write(f'X_COORDINATES {reader.dimensions[0]} double\n')
            f.write(' '.join(f'{v:.15e}' for v in reader.x_coords) + '\n')
            f.write(f'Y_COORDINATES {reader.dimensions[1]} double\n')
            f.write(' '.join(f'{v:.15e}' for v in reader.y_coords) + '\n')
            f.write(f'Z_COORDINATES {reader.dimensions[2]} double\n')
            f.write(' '.join(f'{v:.15e}' for v in reader.z_coords) + '\n')
        elif reader.grid_type == 'UNSTRUCTURED_GRID' and reader.points is not None and reader.cells is not None:
            npts = len(reader.points)
            f.write('DATASET UNSTRUCTURED_GRID\n')
            f.write(f'POINTS {npts} double\n')
            for pt in reader.points:
                f.write(f'{pt[0]:.15e} {pt[1]:.15e} {pt[2]:.15e}\n')
            ti = sum(len(c) + 1 for c in reader.cells)
            f.write(f'CELLS {ncells} {ti}\n')
            for cell in reader.cells:
                f.write(f'{len(cell)} ' + ' '.join(str(n) for n in cell) + '\n')
            f.write(f'CELL_TYPES {ncells}\n')
            ctype = {2: 3, 3: 5, 4: 10}
            for cell in reader.cells:
                f.write(f'{ctype.get(len(cell), 5)}\n')
        else:
            print(f'  [ERROR] Cannot write VTK: grid data missing')
            return

        # Collect fields to write: original nrho/T/reaction_rate + new
        keep = {}
        for name in reader.fields:
            nl = name.lower().strip()
            if any(p in nl for p in ['nrho_mean', 'ttr_mean', 'reaction_rate']):
                clean_name = _sanitize_vtk_name(name)
                keep[clean_name] = reader.fields[name][3]

        all_fields = {**keep, **new_fields}
        f.write(f'CELL_DATA {ncells}\n')
        f.write(f'FIELD FieldData {len(all_fields)}\n')
        for name, data in all_fields.items():
            clean_name = _sanitize_vtk_name(name)
            f.write(f'{clean_name} 1 {ncells} double\n')
            _write_vtk_data(f, data, ncells)

    print(f'  Output VTK: {output_file}')


# ==================== Main Processing ====================

# Default configuration for Deuterium case
DEFAULT_CONFIG = {
    # Beam species contributing to fusion (name → per-nucleon factor)
    # D+  has 1 deuteron  → factor 1
    # D   has 1 deuteron  → factor 1
    # D2+ has 2 deuterons → factor 2 (each can fuse with background)
    # D2  has 2 deuterons → factor 2
    'beam_species': {
        'D+':  {'nucleon_factor': 1, 'mass_kg': 3.344e-27},
        'D':   {'nucleon_factor': 1, 'mass_kg': 3.344e-27},
        'D2+': {'nucleon_factor': 2, 'mass_kg': 6.689e-27},
        'D2':  {'nucleon_factor': 2, 'mass_kg': 6.689e-27},
    },
    'bg_density': 1.5e20,       # n(D2) background [1/m³]
    'bg_temp_K': 300.0,         # background temperature [K]
    'bg_nucleons_per_mol': 2,   # D2 has 2 deuterons
    'fusion_bias_factor': 1.0e8,  # Monte-Carlo sampling bias in R28-R31
}


def process_vtk(vtk_file, config, args):
    """Process a single VTK file and compute fusion rate distribution."""
    print(f'\n{"="*70}')
    print(f'  File: {vtk_file}')
    print(f'{"="*70}')

    reader = VTKReader(vtk_file)
    print(f'  Grid: {reader.grid_type}, cells: {reader.ncells}, binary: {reader.is_binary}')

    if args.list_fields:
        print(f'\n  Available fields ({len(reader.fields)}):')
        for name in sorted(reader.fields.keys()):
            info = reader.fields[name]
            print(f'    {name.strip():40s} ncomp={info[0]}, n={info[1]}, type={info[2]}')
        return None

    ncells = reader.ncells
    n_bg = config['bg_density']
    n_bg_nucleons = n_bg * config['bg_nucleons_per_mol']

    print(f'  Background: n_D2 = {n_bg:.3e} /m³ → n_D(nucleon) = {n_bg_nucleons:.3e} /m³')

    # Accumulate total fusion rate from all beam species
    # Defer R_total size until we know actual field size
    R_total = None
    new_fields = {}

    for sp_name, sp_info in config['beam_species'].items():
        # Find number density field
        nrho_data, nrho_key = reader.find_field('nrho_mean', sp_name)
        T_data, T_key = reader.find_field('ttr_mean', sp_name)

        if nrho_data is None:
            # Try alternative: Ttr (case variations)
            nrho_data, nrho_key = reader.find_field('nrho_mean', sp_name.replace('+', ''))
        if T_data is None:
            T_data, T_key = reader.find_field('ttr_mean', sp_name.replace('+', ''))

        if nrho_data is None:
            print(f'  [SKIP] {sp_name}: field nrho_mean not found')
            continue
        if T_data is None:
            print(f'  [SKIP] {sp_name}: field Ttr_mean not found')
            continue

        print(f'  Species: {sp_name}')
        print(f'    nrho field: {nrho_key}')
        print(f'    T    field: {T_key}')

        nrho = np.array(nrho_data, dtype=np.float64)
        T_K  = np.array(T_data, dtype=np.float64)

        # Convert molecular translational temperatures to the temperature of
        # their co-moving deuterons, then form the D-D relative-temperature
        # parameter used by the Bosch-Hale Maxwellian reactivity.
        m_beam = sp_info['mass_kg']
        m_target = 6.689e-27  # D2 mass
        m_deuteron = 3.344e-27
        T_bg = config['bg_temp_K']
        T_beam_D = T_K * (m_deuteron / m_beam)
        T_bg_D = T_bg * (m_deuteron / m_target)
        T_eff = 0.5 * (T_beam_D + T_bg_D)

        # Convert to keV
        T_keV = T_eff * KB / (QE * 1e3)

        # Compute <σv> (Bosch-Hale, per deuteron-deuteron pair)
        sv = bosch_hale_sigmav(T_keV, reaction='dd_total')

        # Number of deuterons in beam particle
        f_nucleon = sp_info['nucleon_factor']

        # Reaction rate density:
        # R = n_beam * f_nucleon * n_bg_nucleons * <σv>(T_eff)
        # Factor 1/2 NOT needed: beam-target geometry (not thermonuclear)
        R_sp = nrho * f_nucleon * n_bg_nucleons * sv

        # Initialize R_total on first species with correct size
        if R_total is None:
            R_total = np.zeros_like(R_sp)
            ncells = len(R_sp)  # update ncells from actual data

        # Statistics
        mask_valid = T_keV > 0.5  # Bosch-Hale valid range
        print(f'    n_max = {np.max(nrho):.3e} /m³')
        print(f'    T_max = {np.max(T_K):.1f} K = {np.max(T_keV):.4f} keV')
        print(f'    <σv>_max = {np.max(sv):.3e} m³/s')
        print(f'    R_max = {np.max(R_sp):.3e} reactions/(m³·s)')
        print(f'    Cells with T > 0.5 keV: {np.sum(mask_valid)} / {ncells}')

        new_fields[f'R_analytical_{sp_name}'] = R_sp
        new_fields[f'sigmav_{sp_name}'] = sv
        new_fields[f'T_eff_keV_{sp_name}'] = T_keV
        R_total += R_sp

    if R_total is None:
        R_total = np.zeros(ncells)
    new_fields['R_analytical_total'] = R_total
    reader.ncells = ncells  # ensure reader has correct ncells for output
    print(f'\n  Total analytical R_max = {np.max(R_total):.3e} reactions/(m³·s)')

    # Compare with simulation-counted rates if available
    sim_rates = {}
    for name in reader.fields:
        if 'reaction_rate' in name.lower().strip():
            sim_rates[name.strip()] = reader.fields[name][3]
    if sim_rates:
        print(f'\n  Simulation-counted reaction rates found:')
        bias_factor = config['fusion_bias_factor']
        for name, data in sim_rates.items():
            d = np.array(data, dtype=np.float64)
            d_physical = d / bias_factor
            print(f'    {name}: biased max = {np.max(d):.3e}, '
                  f'physical max = {np.max(d_physical):.3e}')

    # Write output VTK
    basename = os.path.splitext(os.path.basename(vtk_file))[0]
    outdir = os.path.dirname(vtk_file) or '.'
    out_vtk = os.path.join(outdir, f'{basename}_fusion_rate.vtk')
    write_vtk_ascii(reader, new_fields, out_vtk)

    # Generate static plots if matplotlib available
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        _make_plots(reader, new_fields, sim_rates, basename, outdir, config)
    except ImportError:
        print('  [INFO] matplotlib not available, skipping plots')

    # Open interactive window if requested
    if getattr(args, 'interactive', False):
        try:
            import matplotlib
            matplotlib.use('TkAgg', force=True)
            import importlib
            import matplotlib.pyplot as plt
            importlib.reload(plt)
            _interactive_axial_plot(reader, new_fields, basename, config)
        except Exception as e:
            print(f'  [WARN] Interactive plot failed: {e}')
            print('  [TIP] Try: pip install matplotlib')

    return new_fields


def _make_plots(reader, new_fields, sim_rates, basename, outdir, config):
    """Generate diagnostic plots."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm

    cx, cy = reader.get_cell_centers()
    if cx is None:
        print('  [INFO] Cannot determine cell centers, skipping plots')
        return

    R_total = new_fields.get('R_analytical_total', np.zeros(reader.ncells))

    # ---- Plot 1: Analytical rate 2D map ----
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))
    mask_pos = R_total > 0
    if np.any(mask_pos):
        sc = ax.scatter(cx * 100, cy * 100, c=R_total, s=2, cmap='hot',
                        norm=LogNorm(vmin=max(R_total[mask_pos].min(), 1e-10),
                                     vmax=R_total.max()))
        plt.colorbar(sc, ax=ax, label=r'$R$ [reactions/(m³·s)]')
    else:
        ax.scatter(cx * 100, cy * 100, c='gray', s=2)
        ax.set_title('No fusion reactions (T too low)')
    ax.set_xlabel('X [cm]')
    ax.set_ylabel('R [cm]')
    ax.set_title(f'D-D Fusion Rate Density (Analytical)\n{basename}')
    ax.set_aspect('equal')
    plot_file = os.path.join(outdir, f'{basename}_fusion_rate_map.png')
    fig.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Plot: {plot_file}')

    # ---- Plot 2: Compare analytical vs simulation-counted ----
    if sim_rates:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))

        ax = axes[0]
        if np.any(mask_pos):
            sc = ax.scatter(cx * 100, cy * 100, c=R_total, s=2, cmap='hot',
                            norm=LogNorm(vmin=max(R_total[mask_pos].min(), 1e-10),
                                         vmax=R_total.max()))
            plt.colorbar(sc, ax=ax, label=r'$R$ [reactions/(m³·s)]')
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('R [cm]')
        ax.set_title('Analytical: n·n·<σv>(T)')
        ax.set_aspect('equal')

        ax = axes[1]
        # Sum all simulation fusion rates (R28-R31)
        R_sim_total = np.zeros(reader.ncells)
        for name, data in sim_rates.items():
            # Extract reaction number
            m = re.search(r'(\d+)', name)
            if m:
                rnum = int(m.group(1))
                if rnum >= 28:  # Fusion reactions start at R28
                    R_sim_total += (np.array(data, dtype=np.float64)
                                    / config['fusion_bias_factor'])
        mask_sim = R_sim_total > 0
        if np.any(mask_sim):
            vmin = max(R_sim_total[mask_sim].min(), 1e-10)
            vmax = R_sim_total.max()
            sc = ax.scatter(cx * 100, cy * 100, c=R_sim_total, s=2, cmap='hot',
                            norm=LogNorm(vmin=vmin, vmax=vmax))
            plt.colorbar(sc, ax=ax, label=r'$R$ [reactions/(m³·s)]')
        ax.set_xlabel('X [cm]')
        ax.set_ylabel('R [cm]')
        ax.set_title('Simulation: de-biased FNUM × count / (V × T)')
        ax.set_aspect('equal')

        fig.suptitle(f'D-D Fusion Rate Comparison — {basename}', fontsize=14)
        fig.tight_layout()
        plot_file = os.path.join(outdir, f'{basename}_fusion_rate_comparison.png')
        fig.savefig(plot_file, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f'  Plot: {plot_file}')

    # ---- Plot 3: Reactivity curve ----
    fig, ax = plt.subplots(figsize=(8, 6))
    T_range = np.logspace(-1, 2.5, 500)  # 0.1 to ~300 keV
    sv_n = bosch_hale_sigmav(T_range, 'dd_n_he3')
    sv_p = bosch_hale_sigmav(T_range, 'dd_p_t')
    sv_t = bosch_hale_sigmav(T_range, 'dd_total')
    ax.loglog(T_range, sv_t, 'k-', linewidth=2, label='D-D total')
    ax.loglog(T_range, sv_n, 'b--', label=r'D(d,n)$^3$He')
    ax.loglog(T_range, sv_p, 'r--', label=r'D(d,p)T')

    # Mark range of temperatures in simulation
    for sp_name in config['beam_species']:
        key = f'T_eff_keV_{sp_name}'
        if key in new_fields:
            T_sp = new_fields[key]
            T_valid = T_sp[T_sp > 0.5]
            if len(T_valid) > 0:
                ax.axvspan(T_valid.min(), T_valid.max(), alpha=0.15, color='green',
                           label=f'{sp_name} T range')

    ax.set_xlabel('T [keV]')
    ax.set_ylabel(r'$\langle\sigma v\rangle$ [m³/s]')
    ax.set_title('D-D Fusion Reactivity (Bosch-Hale 1992)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_xlim(0.1, 300)
    ax.set_ylim(1e-35, 1e-20)
    plot_file = os.path.join(outdir, f'{basename}_reactivity_curve.png')
    fig.savefig(plot_file, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f'  Plot: {plot_file}')

    # ---- Plot 4: Axial profile R(x), smoothed and mirrored about x=0 ----
    R_total = new_fields.get('R_analytical_total', np.zeros(reader.ncells))
    if cx is not None and len(R_total) > 0 and np.any(R_total > 0):
        # Bin data along X axis, averaging R over the radial (Y) direction
        x_min, x_max = np.min(cx), np.max(cx)
        n_bins = min(200, max(50, int(np.sqrt(reader.ncells))))
        x_edges = np.linspace(x_min, x_max, n_bins + 1)
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])

        R_binned = np.zeros(n_bins)
        weight_sum = np.zeros(n_bins)

        for i in range(len(cx)):
            idx = np.searchsorted(x_edges, cx[i]) - 1
            idx = max(0, min(n_bins - 1, idx))
            # Volume-weight: for axisymmetric, weight by r (= cy) to account
            # for annular volume; for non-axi, weight uniformly
            w = max(cy[i], 1e-10) if np.any(cy > 0) else 1.0
            R_binned[idx] += R_total[i] * w
            weight_sum[idx] += w

        mask_valid = weight_sum > 0
        R_binned[mask_valid] /= weight_sum[mask_valid]

        # Gaussian smoothing
        try:
            from scipy.ndimage import gaussian_filter1d
            R_smooth = gaussian_filter1d(R_binned, sigma=2.0)
        except ImportError:
            # Simple moving average fallback
            kernel_size = 5
            kernel = np.ones(kernel_size) / kernel_size
            R_smooth = np.convolve(R_binned, kernel, mode='same')

        # Mirror about x=0: create symmetric profile [-x_max, x_max]
        x_cm = x_centers * 100  # convert to cm
        x_mirror = np.concatenate([-x_cm[::-1], x_cm])
        R_mirror = np.concatenate([R_smooth[::-1], R_smooth])

        fig, ax = plt.subplots(figsize=(12, 5))
        ax.plot(x_mirror, R_mirror, 'b-', linewidth=2, label='Analytical D-D rate')
        ax.fill_between(x_mirror, 0, R_mirror, alpha=0.15, color='blue')
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5, label='Axis of symmetry')
        ax.set_xlabel('X [cm]', fontsize=13)
        ax.set_ylabel(r'$R$ [reactions/(m³·s)]', fontsize=13)
        ax.set_title(f'D-D Fusion Rate — Axial Profile (Mirrored)\n{basename}', fontsize=14)
        ax.legend(fontsize=11)
        ax.grid(True, alpha=0.3)
        ax.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))

        # Use log scale if dynamic range is large
        if np.max(R_mirror) > 0 and np.max(R_mirror) / max(np.min(R_mirror[R_mirror > 0]), 1e-30) > 100:
            ax.set_yscale('log')
            ax.set_ylim(bottom=max(np.min(R_mirror[R_mirror > 0]) * 0.1, 1e-10))

        fig.tight_layout()
        plot_file = os.path.join(outdir, f'{basename}_axial_profile.png')
        fig.savefig(plot_file, dpi=200, bbox_inches='tight')
        plt.close(fig)
        print(f'  Plot: {plot_file}')


def _interactive_axial_plot(reader, new_fields, basename):
    """Open interactive matplotlib window with sliders for axial profile tuning."""
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider, CheckButtons, RadioButtons

    cx, cy = reader.get_cell_centers()
    R_total = new_fields.get('R_analytical_total', np.zeros(reader.ncells))

    if cx is None or not np.any(R_total > 0):
        print('  [INFO] No data for interactive plot')
        return

    x_raw = cx.copy()
    y_raw = cy.copy()
    R_raw = R_total.copy()
    has_radial = np.any(y_raw > 0)

    # ---- Fit model definitions ----
    def _gauss(x, A, mu, sigma):
        return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

    def _double_gauss(x, A1, mu1, s1, A2, mu2, s2):
        return A1 * np.exp(-0.5 * ((x - mu1) / s1) ** 2) + \
               A2 * np.exp(-0.5 * ((x - mu2) / s2) ** 2)

    def _lorentzian(x, A, x0, gamma):
        return A * gamma ** 2 / ((x - x0) ** 2 + gamma ** 2)

    def _try_fit(model_name, x, y):
        """Attempt curve fitting, returning (y_fit, param_str) or None."""
        try:
            from scipy.optimize import curve_fit
        except ImportError:
            return None, 'scipy not installed'

        y_max = np.max(y)
        x_peak = x[np.argmax(y)]
        half = y_max / 2
        above = np.where(y > half)[0]
        if len(above) > 1:
            w_est = (x[above[-1]] - x[above[0]]) / 2.355
        else:
            w_est = (x[-1] - x[0]) / 6

        try:
            if model_name == 'Gaussian':
                p0 = [y_max, x_peak, w_est]
                popt, _ = curve_fit(_gauss, x, y, p0=p0, maxfev=5000)
                y_fit = _gauss(x, *popt)
                txt = f'A = {popt[0]:.3e}\n'
                txt += f'μ = {popt[1]:.2f} cm\n'
                txt += f'σ = {popt[2]:.2f} cm\n'
                txt += f'FWHM = {abs(popt[2]) * 2.355:.2f} cm'
                return y_fit, txt

            elif model_name == 'Double Gaussian':
                mid = len(x) // 2
                x1 = x[np.argmax(y[:mid])] if mid > 0 else x_peak - w_est
                x2 = x[mid + np.argmax(y[mid:])] if mid < len(x) else x_peak + w_est
                p0 = [y_max, x1, w_est, y_max * 0.5, x2, w_est]
                popt, _ = curve_fit(_double_gauss, x, y, p0=p0, maxfev=10000)
                y_fit = _double_gauss(x, *popt)
                txt = f'Peak 1: A={popt[0]:.2e}, μ={popt[1]:.2f}, σ={popt[2]:.2f}\n'
                txt += f'Peak 2: A={popt[3]:.2e}, μ={popt[4]:.2f}, σ={popt[5]:.2f}'
                return y_fit, txt

            elif model_name == 'Lorentzian':
                p0 = [y_max, x_peak, w_est]
                popt, _ = curve_fit(_lorentzian, x, y, p0=p0, maxfev=5000)
                y_fit = _lorentzian(x, *popt)
                txt = f'A = {popt[0]:.3e}\n'
                txt += f'x₀ = {popt[1]:.2f} cm\n'
                txt += f'γ = {popt[2]:.2f} cm\n'
                txt += f'FWHM = {abs(popt[2]) * 2:.2f} cm'
                return y_fit, txt

            elif model_name == 'Polynomial':
                deg = 8
                coeffs = np.polyfit(x, y, deg)
                y_fit = np.polyval(coeffs, x)
                y_fit = np.maximum(y_fit, 0)
                txt = f'Degree: {deg}\n'
                txt += f'R² = {1 - np.sum((y - y_fit)**2) / np.sum((y - np.mean(y))**2):.6f}'
                return y_fit, txt

            else:
                return None, ''
        except Exception as e:
            return None, f'Fit failed: {e}'

    def compute_profile(n_bins, sigma):
        """Bin, average, smooth, and mirror the axial profile."""
        x_min, x_max = np.min(x_raw), np.max(x_raw)
        x_edges = np.linspace(x_min, x_max, n_bins + 1)
        x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])

        R_binned = np.zeros(n_bins)
        weight_sum = np.zeros(n_bins)

        for i in range(len(x_raw)):
            idx = np.searchsorted(x_edges, x_raw[i]) - 1
            idx = max(0, min(n_bins - 1, idx))
            w = max(y_raw[i], 1e-10) if has_radial else 1.0
            R_binned[idx] += R_raw[i] * w
            weight_sum[idx] += w

        mask = weight_sum > 0
        R_binned[mask] /= weight_sum[mask]

        if sigma > 0.1:
            try:
                from scipy.ndimage import gaussian_filter1d
                R_smooth = gaussian_filter1d(R_binned, sigma=sigma)
            except ImportError:
                k = max(3, int(sigma * 2 + 1))
                kernel = np.ones(k) / k
                R_smooth = np.convolve(R_binned, kernel, mode='same')
        else:
            R_smooth = R_binned.copy()

        x_cm = x_centers * 100
        x_mirror = np.concatenate([-x_cm[::-1], x_cm])
        R_mirror = np.concatenate([R_smooth[::-1], R_smooth])
        return x_mirror, R_mirror

    init_bins = min(200, max(50, int(np.sqrt(reader.ncells))))
    init_sigma = 2.0

    fig = plt.figure(figsize=(16, 8))
    ax = fig.add_axes([0.07, 0.22, 0.52, 0.68])
    ax_sigma = fig.add_axes([0.07, 0.08, 0.52, 0.03])
    ax_bins = fig.add_axes([0.07, 0.03, 0.52, 0.03])
    ax_check = fig.add_axes([0.65, 0.65, 0.15, 0.18])
    ax_radio = fig.add_axes([0.65, 0.35, 0.15, 0.25])
    ax_info = fig.add_axes([0.82, 0.35, 0.16, 0.55])
    ax_info.axis('off')

    x_plot, R_plot = compute_profile(init_bins, init_sigma)
    line, = ax.plot(x_plot, R_plot, 'b-', linewidth=2, label='Data')
    fill = ax.fill_between(x_plot, 0, R_plot, alpha=0.12, color='blue')
    fit_line, = ax.plot([], [], 'r-', linewidth=2, alpha=0.85, label='Fit')
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.4)
    ax.set_xlabel('X [cm]', fontsize=13)
    ax.set_ylabel(r'$R$ [reactions/(m³·s)]', fontsize=13)
    ax.set_title(f'D-D Fusion Rate — Axial Profile\n{basename}', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper left', fontsize=10)

    s_sigma = Slider(ax_sigma, 'Smoothing σ', 0.0, 10.0, valinit=init_sigma, valstep=0.1)
    s_bins = Slider(ax_bins, 'Bins', 20, 500, valinit=init_bins, valstep=1)

    check = CheckButtons(ax_check, ['Log Y-axis', 'Show fill', 'Show fit'],
                         [False, True, True])

    ax_radio.set_title('Fit Model', fontsize=10, fontweight='bold')
    radio = RadioButtons(ax_radio, ['Gaussian', 'Double Gaussian', 'Lorentzian', 'Polynomial'],
                         active=0)

    info_text = ax_info.text(0.0, 0.95, '', fontsize=9, family='monospace',
                             verticalalignment='top', transform=ax_info.transAxes)

    state = {'log': False, 'fill_visible': True, 'show_fit': True,
             'fit_model': 'Gaussian'}

    def update(val=None):
        sigma = s_sigma.val
        n_bins = int(s_bins.val)
        x_new, R_new = compute_profile(n_bins, sigma)

        line.set_xdata(x_new)
        line.set_ydata(R_new)

        nonlocal fill
        fill.remove()
        base = 0 if not state['log'] else max(np.min(R_new[R_new > 0]) * 0.5, 1e-20)
        fill = ax.fill_between(x_new, base, R_new,
                               alpha=0.12 if state['fill_visible'] else 0, color='blue')

        fit_txt = ''
        if state['show_fit']:
            y_fit, fit_txt = _try_fit(state['fit_model'], x_new, R_new)
            if y_fit is not None:
                fit_line.set_xdata(x_new)
                fit_line.set_ydata(y_fit)
                fit_line.set_visible(True)
                ss_res = np.sum((R_new - y_fit) ** 2)
                ss_tot = np.sum((R_new - np.mean(R_new)) ** 2)
                r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
                if 'R²' not in fit_txt:
                    fit_txt += f'\nR² = {r2:.6f}'
            else:
                fit_line.set_visible(False)
        else:
            fit_line.set_visible(False)

        ax.set_xlim(x_new.min(), x_new.max())
        R_pos = R_new[R_new > 0]
        if len(R_pos) > 0:
            if state['log']:
                ax.set_yscale('log')
                ax.set_ylim(R_pos.min() * 0.3, R_pos.max() * 3)
            else:
                ax.set_yscale('linear')
                ax.set_ylim(0, R_pos.max() * 1.15)

        R_pos_all = R_new[R_new > 0]
        info = '── Data ──\n'
        info += f'Max:  {np.max(R_new):.3e}\n'
        if len(R_pos_all) > 0:
            info += f'Mean: {np.mean(R_pos_all):.3e}\n'
        info += f'Bins: {n_bins}  σ: {sigma:.1f}\n'
        if state['show_fit'] and fit_txt:
            info += f'\n── {state["fit_model"]} Fit ──\n'
            info += fit_txt
        info_text.set_text(info)

        ax.legend(loc='upper left', fontsize=10)
        fig.canvas.draw_idle()

    def on_check(label):
        if label == 'Log Y-axis':
            state['log'] = not state['log']
        elif label == 'Show fill':
            state['fill_visible'] = not state['fill_visible']
        elif label == 'Show fit':
            state['show_fit'] = not state['show_fit']
        update()

    def on_radio(label):
        state['fit_model'] = label
        update()

    s_sigma.on_changed(update)
    s_bins.on_changed(update)
    check.on_clicked(on_check)
    radio.on_clicked(on_radio)

    update()
    plt.show()


# ==================== Entry Point ====================

def main():
    parser = argparse.ArgumentParser(
        description='D-D Fusion Reaction Rate Post-Processor for PANTERA PIC-DSMC',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python compute_dd_rate.py results/dsmc_flowfield_10000.vtk
  python compute_dd_rate.py results/dsmc_flowfield_*.vtk --bg-density 1.5e20
  python compute_dd_rate.py results/dsmc_flowfield_10000.vtk --list-fields

Note on cross-section scaling:
  The R28...R31 simulation tables intentionally use a 1e8 Monte-Carlo bias to
  raise fusion-event statistics. Simulation counts/rates must be divided by
  1e8 for physical comparison; analytical Bosch-Hale rates remain unscaled.
        """)
    parser.add_argument('vtk_files', nargs='+', help='VTK file(s) to process (glob patterns supported)')
    parser.add_argument('--bg-density', type=float, default=1.5e20,
                        help='MCC background number density n(D2) [1/m³] (default: 1.5e20)')
    parser.add_argument('--bg-temp', type=float, default=300.0,
                        help='MCC background temperature [K] (default: 300)')
    parser.add_argument('--list-fields', action='store_true',
                        help='Only list available VTK fields, do not compute')
    parser.add_argument('--species', nargs='*', default=None,
                        help='Beam species to include (default: D+ D D2+ D2)')
    parser.add_argument('--reaction', default='dd_total',
                        choices=['dd_total', 'dd_n_he3', 'dd_p_t'],
                        help='D-D reaction branch (default: dd_total)')
    parser.add_argument('-i', '--interactive', action='store_true',
                        help='Open interactive plot window with adjustable sliders')
    args = parser.parse_args()

    # Expand glob patterns
    files = []
    for pattern in args.vtk_files:
        expanded = sorted(glob.glob(pattern))
        if expanded:
            files.extend(expanded)
        else:
            files.append(pattern)

    # Configure
    config = DEFAULT_CONFIG.copy()
    config['bg_density'] = args.bg_density
    config['bg_temp_K'] = args.bg_temp

    if args.species:
        filtered = {}
        for sp in args.species:
            if sp in config['beam_species']:
                filtered[sp] = config['beam_species'][sp]
            else:
                print(f'  [WARN] Unknown species: {sp}')
        config['beam_species'] = filtered

    print('='*70)
    print('  D-D Fusion Rate Post-Processor for PANTERA PIC-DSMC')
    print('='*70)
    print(f'  Files to process: {len(files)}')
    print(f'  Background D2 density: {config["bg_density"]:.3e} /m³')
    print(f'  Background temperature: {config["bg_temp_K"]:.1f} K')
    print(f'  Beam species: {list(config["beam_species"].keys())}')

    for vtk_file in files:
        if not os.path.isfile(vtk_file):
            print(f'  [ERROR] File not found: {vtk_file}')
            continue
        process_vtk(vtk_file, config, args)

    print(f'\n{"="*70}')
    print('  Done.')
    print('='*70)


if __name__ == '__main__':
    main()
