"""
D-D Fusion Reaction Rate Post-Processing (ParaView Script)

Usage:
  1. Open VTK file in ParaView, click Apply
  2. Tools -> Python Shell -> Run Script -> select this file
  3. New pipeline nodes with fusion rate fields will appear

Cross-section includes x1e10 scaling factor!
"""

import paraview.simple as pvs
from paraview.simple import (GetActiveSource, ProgrammableFilter,
    RenameSource, Show, GetDisplayProperties, ColorBy,
    GetColorTransferFunction, GetActiveView, Render)
import numpy as np

# ============================================================
# Configuration
# ============================================================
N_D2_BACKGROUND = 1.5e20   # Background D2 density [m^-3]
M_D  = 3.344e-27           # D atom mass [kg]
M_D2 = 6.689e-27           # D2 molecule mass [kg]
MU_D_D2 = (M_D * M_D2) / (M_D + M_D2)  # Reduced mass [kg]
EV_TO_J = 1.602176634e-19

# Cross-section from dd_fusion.txt (with x1e10 scaling)
CROSS_SECTION_DATA = np.array([
    [0.000000e+00, 0.00000e+00],
    [1.000000e+02, 0.00000e+00],
    [5.000000e+02, 0.00000e+00],
    [1.000000e+03, 3.24796e-37],
    [2.000000e+03, 2.00051e-31],
    [3.000000e+03, 6.65795e-29],
    [4.000000e+03, 2.02748e-27],
    [5.000000e+03, 2.03152e-26],
    [6.000000e+03, 1.09400e-25],
    [7.000000e+03, 3.99875e-25],
    [8.000000e+03, 1.12634e-24],
    [9.000000e+03, 2.63697e-24],
    [1.000000e+04, 5.38484e-24],
    [1.100000e+04, 9.90663e-24],
    [1.200000e+04, 1.67952e-23],
    [1.300000e+04, 2.66716e-23],
    [1.400000e+04, 4.01592e-23],
    [1.500000e+04, 5.78625e-23],
    [1.600000e+04, 8.03506e-23],
    [1.700000e+04, 1.08145e-22],
    [1.800000e+04, 1.41714e-22],
    [1.900000e+04, 1.81462e-22],
    [2.000000e+04, 2.27738e-22],
    [2.500000e+04, 5.65050e-22],
    [3.000000e+04, 1.08632e-21],
    [3.500000e+04, 1.78384e-21],
    [4.000000e+04, 2.63708e-21],
    [4.500000e+04, 3.62066e-21],
    [5.000000e+04, 4.70886e-21],
    [6.000000e+04, 7.10690e-21],
    [7.000000e+04, 9.67646e-21],
    [8.000000e+04, 1.23085e-20],
    [9.000000e+04, 1.49306e-20],
    [1.000000e+05, 1.74967e-20],
    [1.200000e+05, 2.23612e-20],
    [1.500000e+05, 2.88562e-20],
    [2.000000e+05, 3.76332e-20],
    [2.500000e+05, 4.43744e-20],
    [3.000000e+05, 4.96501e-20],
    [3.500000e+05, 5.38764e-20],
    [4.000000e+05, 5.73408e-20],
    [4.500000e+05, 6.02423e-20],
    [5.000000e+05, 6.27207e-20],
    [5.240000e+05, 6.37906e-20],
])

CS_ENERGY = CROSS_SECTION_DATA[:, 0]
CS_SIGMA  = CROSS_SECTION_DATA[:, 1]


def main():
    source = GetActiveSource()
    if source is None:
        print("ERROR: Please open a VTK file first and click Apply!")
        return

    print("=" * 60)
    print("D-D Fusion Reaction Rate Calculation")
    print("=" * 60)

    source.UpdatePipeline()
    info = source.GetDataInformation()
    print("Cells: %d, Points: %d" % (info.GetNumberOfCells(), info.GetNumberOfPoints()))

    cell_info = info.GetCellDataInformation()
    available = []
    for i in range(cell_info.GetNumberOfArrays()):
        available.append(cell_info.GetArrayInformation(i).GetName())

    # Check fields
    for f in ['nrho_mean_D+', 'nrho_mean_D', 'Ttr_mean_D+', 'Ttr_mean_D',
              'Ttrx_mean_D+', 'Ttry_mean_D+', 'Ttrz_mean_D+',
              'vx_mean_D+', 'vy_mean_D+', 'vz_mean_D+']:
        status = "OK" if f in available else "MISSING"
        print("  %s: %s" % (f, status))

    has_R28 = 'reaction_rate_R28' in available
    has_R29 = 'reaction_rate_R29' in available
    print("Built-in R28: %s, R29: %s" % ("YES" if has_R28 else "NO",
                                          "YES" if has_R29 else "NO"))

    # ================================================================
    # Create Programmable Filter
    # ================================================================
    print("\nCreating Programmable Filter...")

    prog_filter = ProgrammableFilter(Input=source)
    prog_filter.OutputDataSetType = 'vtkUnstructuredGrid'

    script = '''
import numpy as np

# ====== Parameters ======
N_D2 = %e
M_D = %e
M_D2 = %e
MU = %e
EV_TO_J = %e
kB = 1.380649e-23  # Boltzmann constant [J/K]
kB_eV = 8.617333e-5  # [eV/K]

CS_E = np.array(%s)
CS_S = np.array(%s)

def sigma_interp(E_eV):
    E_eV = np.asarray(E_eV, dtype=float)
    result = np.zeros_like(E_eV)
    mask = E_eV >= CS_E[3]
    if not np.any(mask):
        return result
    high_mask = E_eV > CS_E[-1]
    result[high_mask] = CS_S[-1]
    interp_mask = mask & ~high_mask
    if np.any(interp_mask):
        valid = CS_S > 0
        log_E = np.log10(CS_E[valid])
        log_S = np.log10(CS_S[valid])
        result[interp_mask] = 10**np.interp(
            np.log10(np.maximum(E_eV[interp_mask], 1.0)), log_E, log_S)
    return result

# ====== Read VTK data ======
inp = self.GetInputDataObject(0, 0)
out = self.GetOutputDataObject(0)
out.ShallowCopy(inp)
n_cells = inp.GetNumberOfCells()

def get_array(name):
    arr = inp.GetCellData().GetArray(name)
    if arr is None:
        print("WARNING: %%s not found" %% name)
        return np.zeros(n_cells)
    return np.array([arr.GetValue(i) for i in range(n_cells)])

# ====== D+ fields ======
nrho_Dp = get_array("nrho_mean_D+")
vx_Dp   = get_array("vx_mean_D+")
vy_Dp   = get_array("vy_mean_D+")
vz_Dp   = get_array("vz_mean_D+")
Ttrx_Dp = get_array("Ttrx_mean_D+")  # directional temp [K]
Ttry_Dp = get_array("Ttry_mean_D+")
Ttrz_Dp = get_array("Ttrz_mean_D+")
Ttr_Dp  = get_array("Ttr_mean_D+")

# ====== D fields ======
nrho_D  = get_array("nrho_mean_D")
vx_D    = get_array("vx_mean_D")
vy_D    = get_array("vy_mean_D")
vz_D    = get_array("vz_mean_D")
Ttrx_D  = get_array("Ttrx_mean_D")
Ttry_D  = get_array("Ttry_mean_D")
Ttrz_D  = get_array("Ttrz_mean_D")
Ttr_D   = get_array("Ttr_mean_D")

# ====== DEBUG: print field statistics ======
print("\\n--- DEBUG: Field Statistics ---")
print("nrho_D+: min=%%e max=%%e mean=%%e" %% (nrho_Dp.min(), nrho_Dp.max(), nrho_Dp.mean()))
print("nrho_D:  min=%%e max=%%e mean=%%e" %% (nrho_D.min(), nrho_D.max(), nrho_D.mean()))
print("vx_D+:   min=%%e max=%%e" %% (vx_Dp.min(), vx_Dp.max()))
print("vy_D+:   min=%%e max=%%e" %% (vy_Dp.min(), vy_Dp.max()))
print("vz_D+:   min=%%e max=%%e" %% (vz_Dp.min(), vz_Dp.max()))
print("Ttr_D+:  min=%%e max=%%e K (=%%e - %%e eV)" %% (Ttr_Dp.min(), Ttr_Dp.max(), Ttr_Dp.min()*kB_eV, Ttr_Dp.max()*kB_eV))
print("Ttr_D:   min=%%e max=%%e K (=%%e - %%e eV)" %% (Ttr_D.min(), Ttr_D.max(), Ttr_D.min()*kB_eV, Ttr_D.max()*kB_eV))

# ====== Compute TOTAL kinetic energy per particle ======
# Total KE = drift KE + thermal KE
# E_total = 0.5 * m * v_drift^2 + 3/2 * kT
# This captures BOTH directed beam energy AND thermal energy

# Drift speed
v_drift_Dp = np.sqrt(vx_Dp**2 + vy_Dp**2 + vz_Dp**2)
v_drift_D  = np.sqrt(vx_D**2 + vy_D**2 + vz_D**2)

# Drift KE [J]
E_drift_Dp = 0.5 * M_D * v_drift_Dp**2
E_drift_D  = 0.5 * M_D * v_drift_D**2

# Thermal KE [J] = 3/2 * kB * T
E_thermal_Dp = 1.5 * kB * Ttr_Dp
E_thermal_D  = 1.5 * kB * Ttr_D

# Total KE per particle [J]
E_total_Dp = E_drift_Dp + E_thermal_Dp
E_total_D  = E_drift_D + E_thermal_D

# Effective speed from total KE: v_eff = sqrt(2 * E_total / m)
v_eff_Dp = np.sqrt(2.0 * E_total_Dp / M_D)
v_eff_D  = np.sqrt(2.0 * E_total_D / M_D)

# Convert to CM frame energy [eV]
# E_cm = (M_D2 / (M_D + M_D2)) * E_lab = 2/3 * E_lab
mass_ratio = M_D2 / (M_D + M_D2)  # = 2/3
E_cm_Dp_eV = mass_ratio * E_total_Dp / EV_TO_J
E_cm_D_eV  = mass_ratio * E_total_D / EV_TO_J

print("\\n--- Energy Statistics ---")
print("D+ drift energy:   min=%%e max=%%e eV" %% (E_drift_Dp.min()/EV_TO_J, E_drift_Dp.max()/EV_TO_J))
print("D+ thermal energy: min=%%e max=%%e eV" %% (E_thermal_Dp.min()/EV_TO_J, E_thermal_Dp.max()/EV_TO_J))
print("D+ total energy:   min=%%e max=%%e eV" %% (E_total_Dp.min()/EV_TO_J, E_total_Dp.max()/EV_TO_J))
print("D+ CM energy:      min=%%e max=%%e eV" %% (E_cm_Dp_eV.min(), E_cm_Dp_eV.max()))
print("D  drift energy:   min=%%e max=%%e eV" %% (E_drift_D.min()/EV_TO_J, E_drift_D.max()/EV_TO_J))
print("D  thermal energy: min=%%e max=%%e eV" %% (E_thermal_D.min()/EV_TO_J, E_thermal_D.max()/EV_TO_J))
print("D  CM energy:      min=%%e max=%%e eV" %% (E_cm_D_eV.min(), E_cm_D_eV.max()))

# Cross-section at CM energy
sigma_Dp = sigma_interp(E_cm_Dp_eV)
sigma_D  = sigma_interp(E_cm_D_eV)

print("\\n--- Cross-section ---")
print("sigma_D+: min=%%e max=%%e m^2" %% (sigma_Dp.min(), sigma_Dp.max()))
print("sigma_D:  min=%%e max=%%e m^2" %% (sigma_D.min(), sigma_D.max()))
print("Cells with sigma_D+ > 0: %%d / %%d" %% (np.sum(sigma_Dp > 0), n_cells))
print("Cells with sigma_D  > 0: %%d / %%d" %% (np.sum(sigma_D > 0), n_cells))

# Reaction rate: R = n_projectile * n_D2 * sigma(E_cm) * v_eff
R28 = nrho_Dp * N_D2 * sigma_Dp * v_eff_Dp
R29 = nrho_D  * N_D2 * sigma_D  * v_eff_D
R_total = R28 + R29

print("\\n--- Reaction Rates ---")
print("R28 (D+ + D2): max=%%e reactions/m3/s" %% R28.max())
print("R29 (D  + D2): max=%%e reactions/m3/s" %% R29.max())
print("R_total:       max=%%e reactions/m3/s" %% R_total.max())
print("Cells with R_total > 0: %%d / %%d" %% (np.sum(R_total > 0), n_cells))

# ====== Output arrays ======
from vtkmodules.vtkCommonCore import vtkDoubleArray

def add_array(name, data):
    arr = vtkDoubleArray()
    arr.SetName(name)
    arr.SetNumberOfTuples(n_cells)
    for i in range(n_cells):
        arr.SetValue(i, float(data[i]))
    out.GetCellData().AddArray(arr)

add_array("fusion_rate_R28_Dp_D2", R28)
add_array("fusion_rate_R29_D_D2", R29)
add_array("fusion_rate_total", R_total)
add_array("E_cm_Dp_eV", E_cm_Dp_eV)
add_array("E_cm_D_eV", E_cm_D_eV)
add_array("E_total_Dp_eV", E_total_Dp / EV_TO_J)
add_array("E_total_D_eV", E_total_D / EV_TO_J)
add_array("E_drift_Dp_eV", E_drift_Dp / EV_TO_J)
add_array("E_thermal_Dp_eV", E_thermal_Dp / EV_TO_J)
add_array("sigma_Dp_m2", sigma_Dp)
add_array("sigma_D_m2", sigma_D)
add_array("v_eff_Dp_ms", v_eff_Dp)
add_array("v_eff_D_ms", v_eff_D)

print("\\nDone! Filter output ready.")
''' % (N_D2_BACKGROUND, M_D, M_D2, MU_D_D2, EV_TO_J,
       repr(CS_ENERGY.tolist()), repr(CS_SIGMA.tolist()))

    prog_filter.Script = script
    RenameSource("DD_Fusion_Rate", prog_filter)
    prog_filter.UpdatePipeline()

    print("\nFilter created! New fields:")
    print("  - fusion_rate_R28_Dp_D2  (D+ + D2 rate [reactions/m3/s])")
    print("  - fusion_rate_R29_D_D2   (D  + D2 rate [reactions/m3/s])")
    print("  - fusion_rate_total      (Total rate [reactions/m3/s])")
    print("  - E_cm_Dp_eV, E_cm_D_eV (CM energy [eV])")
    print("  - E_total_*_eV           (Total KE = drift + thermal [eV])")
    print("  - E_drift_*_eV           (Drift KE only [eV])")
    print("  - E_thermal_*_eV         (Thermal KE only [eV])")
    print("  - sigma_*_m2             (Cross-section [m^2])")

    # Auto-display
    Show(prog_filter)
    display = GetDisplayProperties(prog_filter)
    try:
        ColorBy(display, ('CELLS', 'fusion_rate_total'))
        lut = GetColorTransferFunction('fusion_rate_total')
        lut.ApplyPreset('Cool to Warm', True)
        display.SetScalarBarVisibility(GetActiveView(), True)
        Render()
    except Exception as e:
        print("Auto-coloring note: %s" % str(e))

    print("\nCheck Python Shell output for DEBUG statistics.")
    print("If all energies < 1000 eV, no fusion occurs (threshold).")


main()
