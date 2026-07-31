import hashlib
from pathlib import Path
import unittest


ROOT = Path(__file__).resolve().parents[1]


class P0StaticRegressionTests(unittest.TestCase):
    def test_d2_vibrational_table_and_reaction(self):
        table_path = ROOT / "D2" / "EXCITATION2.txt"
        normalized_bytes = table_path.read_bytes().replace(b"\r\n", b"\n")
        self.assertEqual(
            hashlib.sha256(normalized_bytes).hexdigest(),
            "d27fdfc9c57f6e778b49dc74a62e026beb16544fdd8271183f9da77708c63eb4",
        )
        lines = table_path.read_text(encoding="utf-8").splitlines()
        columns = next(i for i, line in enumerate(lines) if line.startswith("COLUMNS:"))

        rows = []
        for line in lines[columns + 2 :]:
            if line.startswith("-----"):
                break
            fields = line.split()
            self.assertEqual(len(fields), 2)
            rows.append((float(fields[0]), float(fields[1])))

        self.assertEqual(len(rows), 4963)
        self.assertAlmostEqual(rows[0][0], 0.38, places=12)
        self.assertAlmostEqual(rows[-1][0], 50.0, places=12)
        self.assertTrue(all(b[0] > a[0] for a, b in zip(rows, rows[1:])))
        self.assertTrue(all(cross_section >= 0.0 for _, cross_section in rows))
        peak_energy, peak_cross_section = max(rows, key=lambda row: row[1])
        self.assertGreaterEqual(peak_energy, 4.45)
        self.assertLessEqual(peak_energy, 4.47)
        self.assertAlmostEqual(peak_cross_section, 2.2898e-21, delta=1e-27)

        reaction = (ROOT / "D2" / "deuterium_modified.react").read_text(
            encoding="utf-8"
        )
        self.assertIn("lxcat 0.0 ELASTIC.txt elastic", reaction)
        self.assertIn("lxcat 0.3712 EXCITATION2.txt vibrational", reaction)

    def test_explicit_two_body_kinematics_bypass_bl_redistribution(self):
        source = (ROOT / "src" / "collisions.f90").read_text(encoding="utf-8")
        elastic_marker = (
            "IF (REACTIONS(JR)%KINEMATICS == REACTION_KIN_ELASTIC) THEN"
        )
        vibration_marker = (
            "ELSE IF (REACTIONS(JR)%KINEMATICS == "
            "REACTION_KIN_VIBRATIONAL) THEN"
        )
        statistical_marker = "ELSE IF (.NOT. REACTIONS(JR)%IS_CEX) THEN"

        elastic = source.split(elastic_marker, 1)[1].split(vibration_marker, 1)[0]
        vibration = source.split(vibration_marker, 1)[1].split(
            statistical_marker, 1
        )[0]

        self.assertIn("CALL RELATIVISTIC_TWO_BODY_SCATTER", elastic)
        self.assertIn("C1, C2, 0.d0", elastic)
        self.assertNotIn("CALL HS_SCATTER", elastic)
        self.assertNotIn("COLL_INTERNAL_ENERGY", elastic)
        self.assertIn("CALL RELATIVISTIC_TWO_BODY_SCATTER", vibration)
        self.assertIn("C1, C2, EA", vibration)
        self.assertIn("particles(JP2)%EVIB = 0.d0", vibration)
        self.assertIn("particles(JP2)%EVIB = particles(JP2)%EVIB + EA", vibration)
        self.assertNotIn("EA_CM", vibration)
        self.assertNotIn("ECOLL", vibration)
        self.assertNotIn("CALL HS_SCATTER", vibration)
        self.assertNotIn("COLL_INTERNAL_ENERGY", vibration)
        legacy = source.split(statistical_marker, 1)[1]
        self.assertIn(
            "ECOLL = MAX(0.d0, ETR - EA_CM + REACTIONS(JR)%Q_VALUE)", legacy
        )

    def test_d2_cathode_effective_absorption_matches_wall_and_energy_cll(self):
        wall_lines = (ROOT / "D2" / "cathode.wall").read_text(
            encoding="utf-8"
        ).splitlines()
        electron_absorption = []
        for index, line in enumerate(wall_lines):
            if line.strip() == "e- --> none":
                probability, minimum, maximum = map(
                    float, wall_lines[index + 1].split()
                )
                electron_absorption.append((probability, minimum, maximum))

        self.assertEqual(len(electron_absorption), 3)

        config_lines = (ROOT / "D2" / "input").read_text(
            encoding="utf-8"
        ).splitlines()
        energy_cll_line = next(
            line for line in config_lines
            if line.strip().startswith("cathode    energy_cll")
        )
        energy_cll_threshold = float(energy_cll_line.split()[-1])
        self.assertEqual(energy_cll_threshold, 100.0)

        def wall_absorption(energy):
            return sum(
                probability
                for probability, minimum, maximum in electron_absorption
                if minimum <= energy <= maximum
            )

        def effective_absorption(energy):
            wall_probability = wall_absorption(energy)
            if energy < energy_cll_threshold:
                return 1.0
            return wall_probability

        self.assertAlmostEqual(effective_absorption(90.0), 1.0)
        self.assertAlmostEqual(effective_absorption(100.0), 0.60)
        self.assertAlmostEqual(effective_absorption(1000.0), 0.60)

        species_line = next(
            line for line in
            (ROOT / "D2" / "D2.species").read_text(encoding="utf-8").splitlines()
            if line.split() and line.split()[0] == "e-"
        )
        species_fields = species_line.split()
        self.assertEqual(float(species_fields[8]), 1.0)
        self.assertEqual(float(species_fields[9]), -1.0)

    def test_d2_boundary_feedback_is_bounded(self):
        def high_energy_absorption(filename):
            lines = (ROOT / "D2" / filename).read_text(
                encoding="utf-8"
            ).splitlines()
            segments = [
                tuple(map(float, lines[index + 1].split()))
                for index, line in enumerate(lines)
                if line.strip() == "e- --> none"
            ]
            return segments[-1]

        self.assertEqual(
            high_energy_absorption("anode.wall"), (0.70, 80.0, 1.0e99)
        )
        self.assertEqual(
            high_energy_absorption("cathode.wall"), (0.60, 80.0, 1.0e99)
        )

        config = (ROOT / "D2" / "input").read_text(encoding="utf-8")
        wall_reaction_line = next(
            line for line in config.splitlines()
            if line.strip().startswith("wall     react")
        )
        active_wall_file = wall_reaction_line.split()[-1]
        self.assertEqual(active_wall_file, "wall.wall")
        self.assertEqual(
            high_energy_absorption(active_wall_file), (0.90, 80.0, 1.0e99)
        )

        active_wall = (ROOT / "D2" / active_wall_file).read_text(
            encoding="utf-8"
        ).splitlines()
        probability_sums = {}
        for index, line in enumerate(active_wall):
            if "-->" not in line or line.strip().startswith("!"):
                continue
            probability, minimum, maximum = map(
                float, active_wall[index + 1].split()
            )
            reactant = line.split("-->", 1)[0].strip()
            if reactant == "e-":
                continue
            key = (reactant, minimum, maximum)
            probability_sums[key] = probability_sums.get(key, 0.0) + probability
        self.assertTrue(probability_sums)
        for key, probability_sum in probability_sums.items():
            with self.subTest(wall_channel=key):
                self.assertAlmostEqual(probability_sum, 1.0)

        materials = {}
        for line in (ROOT / "D2" / "see_materials.dat").read_text(
            encoding="utf-8"
        ).splitlines():
            fields = line.split()
            if fields and fields[0] in {"SUS304", "Duran"}:
                materials[fields[0]] = list(map(float, fields[1:]))
        self.assertEqual(
            (
                materials["SUS304"][0],
                materials["SUS304"][4],
                materials["SUS304"][15],
            ),
            (1.50, 0.42, 0.045),
        )
        self.assertEqual(
            (materials["Duran"][0], materials["Duran"][4], materials["Duran"][15]),
            (1.30, 0.18, 0.022),
        )

    def test_boundary_current_is_tallied_once_after_final_wall_outcome(self):
        timecycle = (ROOT / "src" / "timecycle.f90").read_text(encoding="utf-8")
        self.assertNotIn("TIMESTEP_CHARGE_", timecycle)
        self.assertNotIn("SPICE_NODE_CURRENT", timecycle)
        self.assertEqual(timecycle.count("CALL TALLY_NET_BOUNDARY_CURRENT"), 1)

        energy_cll = timecycle.index(
            "IMPACT_ENERGY_EV < GRID_BC(FACE_PG)%ENERGY_THRESHOLD_EV"
        )
        deposit = timecycle.index("CALL DEPOSIT_NET_SURFACE_CHARGE")
        tally = timecycle.index("CALL TALLY_NET_BOUNDARY_CURRENT")
        reflected_flux = timecycle.index(
            "! Tally reflected particle fluxes", tally
        )
        self.assertLess(energy_cll, deposit)
        self.assertLess(timecycle.rfind("CALL WALL_REACT", 0, tally), deposit)
        self.assertLess(deposit, tally)
        self.assertLess(tally, reflected_flux)

    def test_controller_preserves_cv_then_cc_transition(self):
        fields = (ROOT / "src" / "fields.f90").read_text(encoding="utf-8")
        self.assertIn("CALL ERROR_DRIVEN_PID_STEP", fields)
        control_core = (ROOT / "src" / "boundary_control_core.f90").read_text(
            encoding="utf-8"
        )
        self.assertIn(
            "CURRENT_ERROR = TARGET_CURRENT - MEASURED_CURRENT", control_core
        )
        self.assertIn(
            "VOLTAGE_UNSAT = PREVIOUS_VOLTAGE + PID_SIGN*(", control_core
        )
        self.assertIn(
            "ABS(SMOOTHED_CURRENT) >= ABS(GRID_BC(IPG)%TARGET_CURRENT)", fields
        )
        config = (ROOT / "D2" / "input").read_text(encoding="utf-8")
        self.assertIn("cathode constant_current", config)
        self.assertIn("Boundary_voltage_ramp:", config)
        self.assertIn("cathode 0.0 1.0e-7", config)
        self.assertIn("SMOOTH_STARTUP_VOLTAGE", control_core)
        ramp_branch = fields.index("STARTUP_RAMP_ENABLED")
        target_switch = fields.index(
            "ABS(SMOOTHED_CURRENT) >= ABS(GRID_BC(IPG)%TARGET_CURRENT)"
        )
        self.assertLess(ramp_branch, target_switch)

        initialization = (ROOT / "src" / "initialization.f90").read_text(
            encoding="utf-8"
        )
        self.assertIn("GRID_TYPE /= UNSTRUCTURED", initialization)
        self.assertIn(
            "PIC_TYPE /= EXPLICIT .AND. PIC_TYPE /= EXPLICITLIMITED",
            initialization,
        )
        self.assertIn(
            "constant_current currently supports a negative cathode voltage",
            initialization,
        )
        self.assertIn(
            "constant_current boundary must retain a Dirichlet field condition",
            initialization,
        )
        self.assertIn(
            "Explicit LXCat kinematics currently require Collision_type VAHEDI_MCC",
            initialization,
        )
        self.assertIn("IEEE_IS_FINITE", initialization)

    def test_lxcat_arity_is_checked_before_parameter_access(self):
        source = (ROOT / "src" / "initialization.f90").read_text(
            encoding="utf-8"
        )
        block = source.split(
            "ELSE IF (STRARRAY(1) == 'lxcat' .OR. "
            "STRARRAY(1) == 'lxcat_fusion') THEN",
            1,
        )[1].split("OPEN(UNIT=in4", 1)[0]
        fusion_check = block.index("IF (N_STR /= 4)")
        first_threshold_read = block.index("READ(STRARRAY(2)")
        lxcat_check = block.index("IF (N_STR < 3 .OR. N_STR > 4)")
        second_threshold_read = block.index(
            "READ(STRARRAY(2)", first_threshold_read + 1
        )
        self.assertLess(fusion_check, first_threshold_read)
        self.assertLess(lxcat_check, second_threshold_read)

    def test_constant_current_diagnostics_preserve_existing_columns(self):
        source = (ROOT / "src" / "postprocess.f90").read_text(encoding="utf-8")
        old_columns = ["V_", "I_tot_", "I_ion_", "I_elec_", "I_see_"]
        new_columns = ["I_target_", "I_error_", "CC_limited_"]
        positions = [source.index(f"' {name}'") for name in old_columns + new_columns]
        self.assertEqual(positions, sorted(positions))
        self.assertIn("WRITE(54331,'(5(ES14.6))'", source)
        self.assertIn("WRITE(54331,'(3(ES14.6))'", source)
        self.assertIn("GRID_BC(IPG)%CC_CONTROL_LIMITED", source)

    def test_heavy_particle_lxcat_tables_stop_at_source_maximum(self):
        table_paths = sorted((ROOT / "D2").glob("H2-??-*.txt"))
        self.assertEqual(len(table_paths), 23)

        for table_path in table_paths:
            with self.subTest(table=table_path.name):
                lines = table_path.read_text(encoding="utf-8").splitlines()
                range_line = next(
                    line for line in lines if "Energy range:" in line
                )
                source_max_ev = float(range_line.split()[-2])
                columns = next(
                    i for i, line in enumerate(lines)
                    if line.startswith("COLUMNS:")
                )

                rows = []
                for line in lines[columns + 2 :]:
                    if line.startswith("-----"):
                        break
                    energy, cross_section = map(float, line.split())
                    rows.append((energy, cross_section))

                self.assertGreater(len(rows), 1)
                self.assertTrue(
                    all(b[0] > a[0] for a, b in zip(rows, rows[1:]))
                )
                self.assertTrue(
                    all(cross_section >= 0.0 for _, cross_section in rows)
                )
                self.assertAlmostEqual(
                    rows[-1][0],
                    source_max_ev,
                    delta=max(1.0e-12, source_max_ev * 1.0e-12),
                )
                self.assertTrue(
                    all(energy <= source_max_ev for energy, _ in rows)
                )

    def test_dplus_d2_non_dissociative_ionization_channel(self):
        table_path = ROOT / "D2" / "H2-09-IONIZATION.txt"
        lines = table_path.read_text(encoding="utf-8").splitlines()
        columns = next(i for i, line in enumerate(lines) if line.startswith("COLUMNS:"))

        rows = []
        for line in lines[columns + 2 :]:
            if line.startswith("-----"):
                break
            energy, cross_section = map(float, line.split())
            rows.append((energy, cross_section))

        self.assertGreater(len(rows), 2)
        self.assertAlmostEqual(rows[0][0], 0.1, places=12)
        self.assertAlmostEqual(rows[-1][0], 1.0e5, places=6)
        self.assertTrue(all(b[0] > a[0] for a, b in zip(rows, rows[1:])))

        a1, a2, a3, a4 = 1.864e-4, 1.216, 53.1, 0.897
        for energy_ev, cross_section in rows:
            if energy_ev <= 20.0:
                expected = 0.0
            else:
                x_kev = energy_ev / 1000.0 - 0.020
                expected = (
                    1.0e-20
                    * a1
                    * (x_kev / 0.01361) ** a2
                    / (1.0 + (x_kev / a3) ** (a2 + a4))
                )
            self.assertAlmostEqual(
                cross_section,
                expected,
                delta=max(1.0e-32, 5.0e-5 * abs(expected)),
            )

        reaction = (ROOT / "D2" / "deuterium_modified.react").read_text(
            encoding="utf-8"
        )
        self.assertEqual(
            reaction.count("D+ + D2 --> D+ + D2+ + e-"),
            1,
        )
        self.assertIn("lxcat 2.320500e+01 H2-09-IONIZATION.txt", reaction)
        self.assertNotIn("D+ + D2 --> D+ + D+ + D + e-", reaction)

        species = {}
        for line in (ROOT / "D2" / "D2.species").read_text(
            encoding="utf-8"
        ).splitlines():
            fields = line.split()
            if fields:
                species[fields[0]] = (float(fields[2]), float(fields[9]))
        self.assertEqual(species["D+"][1] + species["D2"][1], 1.0)
        self.assertEqual(
            species["D+"][1] + species["D2+"][1] + species["e-"][1],
            1.0,
        )
        ea_cm_ev = 23.205 * species["D2"][0] / (
            species["D+"][0] + species["D2"][0]
        )
        self.assertAlmostEqual(ea_cm_ev, 15.47, delta=0.005)




if __name__ == "__main__":
    unittest.main()
