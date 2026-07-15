"""
D-D Fusion Reaction Rate Post-Processing (ParaView Script)

Usage:
  1. Open VTK file in ParaView, click Apply
  2. Tools -> Python Shell -> Run Script -> select this file
  3. New pipeline nodes with fusion rate fields will appear

Simulation cross-sections include a 1e8 Monte-Carlo sampling bias.
Physical-rate fields are reported after division by the same bias.
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
FUSION_BIAS_FACTOR = 1.0e8

# Biased effective cross-section from D_fusion.txt; divide by FUSION_BIAS_FACTOR for physics
CROSS_SECTION_DATA = np.array([
    [0.000000e+00, 0.00000e+00],
    [2.000000e+02, 0.00000e+00],
    [1.000000e+03, 0.00000e+00],
    [2.000000e+03, 6.49592e-39],
    [4.000000e+03, 4.00102e-33],
    [6.000000e+03, 1.33159e-30],
    [8.000000e+03, 4.05496e-29],
    [1.000000e+04, 4.06304e-28],
    [1.200000e+04, 2.18800e-27],
    [1.400000e+04, 7.99750e-27],
    [1.600000e+04, 2.25268e-26],
    [1.800000e+04, 5.27394e-26],
    [2.000000e+04, 1.07697e-25],
    [2.200000e+04, 1.98133e-25],
    [2.400000e+04, 3.35904e-25],
    [2.600000e+04, 5.33432e-25],
    [2.800000e+04, 8.03184e-25],
    [3.000000e+04, 1.15725e-24],
    [3.200000e+04, 1.60701e-24],
    [3.400000e+04, 2.16290e-24],
    [3.600000e+04, 2.83428e-24],
    [3.800000e+04, 3.62924e-24],
    [4.000000e+04, 4.55476e-24],
    [5.000000e+04, 1.13010e-23],
    [6.000000e+04, 2.17264e-23],
    [7.000000e+04, 3.56768e-23],
    [8.000000e+04, 5.27416e-23],
    [9.000000e+04, 7.24132e-23],
    [1.000000e+05, 9.41772e-23],
    [1.200000e+05, 1.42138e-22],
    [1.400000e+05, 1.93529e-22],
    [1.600000e+05, 2.46170e-22],
    [1.800000e+05, 2.98612e-22],
    [2.000000e+05, 3.49934e-22],
    [2.400000e+05, 4.47224e-22],
    [3.000000e+05, 5.77124e-22],
    [4.000000e+05, 7.52664e-22],
    [5.000000e+05, 8.87488e-22],
    [6.000000e+05, 9.93002e-22],
    [7.000000e+05, 1.07753e-21],
    [8.000000e+05, 1.14682e-21],
    [9.000000e+05, 1.20485e-21],
    [1.000000e+06, 1.25441e-21],
    [1.048000e+06, 1.27581e-21],
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
FUSION_BIAS_FACTOR = %e
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

# D_fusion.txt is tabulated against projectile laboratory energy.
# For a D projectile and either D nucleus in stationary D2, E_DD,cm = E_lab/2.
# Its factor of two for the two target deuterons is already included in sigma_eff.
E_lab_Dp_eV = E_total_Dp / EV_TO_J
E_lab_D_eV  = E_total_D / EV_TO_J
E_cm_Dp_eV  = 0.5 * E_lab_Dp_eV
E_cm_D_eV   = 0.5 * E_lab_D_eV

print("\\n--- Energy Statistics ---")
print("D+ drift energy:   min=%%e max=%%e eV" %% (E_drift_Dp.min()/EV_TO_J, E_drift_Dp.max()/EV_TO_J))
print("D+ thermal energy: min=%%e max=%%e eV" %% (E_thermal_Dp.min()/EV_TO_J, E_thermal_Dp.max()/EV_TO_J))
print("D+ total energy:   min=%%e max=%%e eV" %% (E_total_Dp.min()/EV_TO_J, E_total_Dp.max()/EV_TO_J))
print("D+ CM energy:      min=%%e max=%%e eV" %% (E_cm_Dp_eV.min(), E_cm_Dp_eV.max()))
print("D  drift energy:   min=%%e max=%%e eV" %% (E_drift_D.min()/EV_TO_J, E_drift_D.max()/EV_TO_J))
print("D  thermal energy: min=%%e max=%%e eV" %% (E_thermal_D.min()/EV_TO_J, E_thermal_D.max()/EV_TO_J))
print("D  CM energy:      min=%%e max=%%e eV" %% (E_cm_D_eV.min(), E_cm_D_eV.max()))

# Biased effective D + D2 cross-section at projectile laboratory energy
sigma_Dp = sigma_interp(E_lab_Dp_eV)
sigma_D  = sigma_interp(E_lab_D_eV)
sigma_Dp_physical = sigma_Dp / FUSION_BIAS_FACTOR
sigma_D_physical  = sigma_D / FUSION_BIAS_FACTOR

print("\\n--- Cross-section ---")
print("sigma_D+ biased:   min=%%e max=%%e m^2" %% (sigma_Dp.min(), sigma_Dp.max()))
print("sigma_D biased:    min=%%e max=%%e m^2" %% (sigma_D.min(), sigma_D.max()))
print("sigma_D+ physical: min=%%e max=%%e m^2" %% (sigma_Dp_physical.min(), sigma_Dp_physical.max()))
print("sigma_D physical:  min=%%e max=%%e m^2" %% (sigma_D_physical.min(), sigma_D_physical.max()))
print("Cells with sigma_D+ > 0: %%d / %%d" %% (np.sum(sigma_Dp > 0), n_cells))
print("Cells with sigma_D  > 0: %%d / %%d" %% (np.sum(sigma_D > 0), n_cells))

# Existing field names retain the biased Monte-Carlo/visualization rate.
R28 = nrho_Dp * N_D2 * sigma_Dp * v_eff_Dp
R29 = nrho_D  * N_D2 * sigma_D  * v_eff_D
R_total = R28 + R29
R28_physical = R28 / FUSION_BIAS_FACTOR
R29_physical = R29 / FUSION_BIAS_FACTOR
R_total_physical = R_total / FUSION_BIAS_FACTOR

print("\\n--- Reaction Rates ---")
print("R28 biased (D+ + D2): max=%%e reactions/m3/s" %% R28.max())
print("R29 biased (D  + D2): max=%%e reactions/m3/s" %% R29.max())
print("R_total biased:       max=%%e reactions/m3/s" %% R_total.max())
print("R_total physical:     max=%%e reactions/m3/s" %% R_total_physical.max())
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
add_array("fusion_rate_physical_R28_Dp_D2", R28_physical)
add_array("fusion_rate_physical_R29_D_D2", R29_physical)
add_array("fusion_rate_physical_total", R_total_physical)
add_array("E_cm_Dp_eV", E_cm_Dp_eV)
add_array("E_cm_D_eV", E_cm_D_eV)
add_array("E_total_Dp_eV", E_total_Dp / EV_TO_J)
add_array("E_total_D_eV", E_total_D / EV_TO_J)
add_array("E_drift_Dp_eV", E_drift_Dp / EV_TO_J)
add_array("E_thermal_Dp_eV", E_thermal_Dp / EV_TO_J)
add_array("sigma_Dp_m2", sigma_Dp)
add_array("sigma_D_m2", sigma_D)
add_array("sigma_Dp_physical_m2", sigma_Dp_physical)
add_array("sigma_D_physical_m2", sigma_D_physical)
add_array("v_eff_Dp_ms", v_eff_Dp)
add_array("v_eff_D_ms", v_eff_D)

print("\\nDone! Filter output ready.")
''' % (N_D2_BACKGROUND, M_D, M_D2, MU_D_D2, EV_TO_J, FUSION_BIAS_FACTOR,
       repr(CS_ENERGY.tolist()), repr(CS_SIGMA.tolist()))

    prog_filter.Script = script
    RenameSource("DD_Fusion_Rate", prog_filter)
    prog_filter.UpdatePipeline()

    print("\nFilter created! New fields:")
    print("  - fusion_rate_R28_Dp_D2  (D+ + D2 rate [reactions/m3/s])")
    print("  - fusion_rate_R29_D_D2   (D  + D2 rate [reactions/m3/s])")
    print("  - fusion_rate_total      (Biased total rate [reactions/m3/s])")
    print("  - fusion_rate_physical_* (De-biased physical rates [reactions/m3/s])")
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
    print("If all projectile energies < 2000 eV, no fusion occurs (effective-table threshold).")


main()
