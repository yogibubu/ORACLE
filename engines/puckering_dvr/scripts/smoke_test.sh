#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_DIR="${TMPDIR:-/tmp}/mw_path_dvr_smoke"
PYTHON="${PYTHON:-python}"

cd "$ROOT_DIR"
rm -rf "$TMP_DIR"
mkdir -p "$TMP_DIR"

"$PYTHON" -m compileall -q scripts
"$PYTHON" scripts/mw_path_dvr.py --help > "$TMP_DIR/help.txt"

cat > "$TMP_DIR/relaxed_scan_double_hybrid.log" <<'EOF'
 Charge = 0 Multiplicity = 1
 Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          1           0        0.000000    0.000000    0.000000
      2          1           0        0.000000    0.000000    0.700000
 ---------------------------------------------------------------------
 SCF Done:  E(RB3LYP) =  -1.0000000000     A.U. after    1 cycles
 E2(TestDH) = -0.1000000000D+00 E(TestDH) = -0.110000000000D+01
 Dipole        = 1.00000000D-01 2.00000000D-01 0.00000000D+00
 Rotational constants (GHZ):     10.0000000     9.0000000     8.0000000
 Step number   1 out of a maximum of  20 on scan point   1 out of   3
 Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          1           0        0.000000    0.000000    0.000000
      2          1           0        0.000000    0.000000    0.800000
 ---------------------------------------------------------------------
 SCF Done:  E(RB3LYP) =  -1.0100000000     A.U. after    1 cycles
 E2(TestDH) = -0.2000000000D+00 E(TestDH) = -0.120000000000D+01
 Dipole        = 2.00000000D-01 1.00000000D-01 0.00000000D+00
 Rotational constants (GHZ):      9.0000000     8.0000000     7.0000000
 Step number   3 out of a maximum of  20 on scan point   2 out of   3
 Input orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          1           0        0.000000    0.000000    0.000000
      2          1           0        0.000000    0.000000    0.900000
 ---------------------------------------------------------------------
 SCF Done:  E(RB3LYP) =  -1.0200000000     A.U. after    1 cycles
 E2(TestDH) = -0.1500000000D+00 E(TestDH) = -0.115000000000D+01
 Dipole        = 3.00000000D-01 0.00000000D+00 0.00000000D+00
 Rotational constants (GHZ):      8.0000000     7.0000000     6.0000000
 Step number   2 out of a maximum of  20 on scan point   3 out of   3
 Normal termination of Gaussian
EOF

cat > "$TMP_DIR/scan_vpt2_properties.csv" <<'EOF'
property,reference_value,vpt2_delta,total_perturbative_delta
A_MHz,9000.0,10.0,15.0
EOF

"$PYTHON" scripts/mw_path_dvr.py \
  --gaussian-log "$TMP_DIR/relaxed_scan_double_hybrid.log" \
  --log-selection all \
  --gaussian-energy post-scf \
  --energy-key energy_hartree --energy-unit hartree \
  --property A_MHz \
  --property dipole_debye \
  --vpt2-property-csv "$TMP_DIR/scan_vpt2_properties.csv" \
  --boundary nonperiodic \
  --solver sinc-dvr \
  --grid 31 \
  --levels 2 \
  --plot-max-state 1 \
  --plot-property-smooth-degree 1 \
  --outdir "$TMP_DIR/out_scan_post_scf" \
  --figdir "$TMP_DIR/figs_scan_post_scf" \
  --prefix scan_post_scf

"$PYTHON" - "$TMP_DIR/out_scan_post_scf/scan_post_scf_profile.csv" <<'PY'
import csv
import math
import sys

with open(sys.argv[1], newline="") as handle:
    rows = list(csv.DictReader(handle))

assert len(rows) == 3
relative = [float(row["relative_energy_cm-1"]) for row in rows]
assert math.isclose(relative[1], 0.0, abs_tol=1e-8)
assert relative[0] > relative[2] > relative[1]
assert math.isclose(float(rows[0]["energy_hartree"]), -1.1, abs_tol=1e-12)
assert float(rows[-1]["s_au"]) > 0.0
PY
test -s "$TMP_DIR/out_scan_post_scf/scan_post_scf_vpt2_property_correction.csv"
test -s "$TMP_DIR/out_scan_post_scf/scan_post_scf_property_comparison.csv"
grep -q "variational_minus_vpt2_1d" "$TMP_DIR/out_scan_post_scf/scan_post_scf_vpt2_property_correction.csv"
grep -q "dipole_debye" "$TMP_DIR/out_scan_post_scf/scan_post_scf_property_comparison.csv"
grep -q "local_taylor" "$TMP_DIR/out_scan_post_scf/scan_post_scf_property_comparison.csv"

cat > "$TMP_DIR/anharmonic_modes.log" <<'EOF'
 Harmonic frequencies (cm**-1)
                      1                      2                      3
 Frequencies --   -100.0000               300.0000               500.0000

 ........................................................
 :      QUADRATIC FORCE CONSTANTS IN NORMAL MODES       :
 :                                                      :
 : FI =  Frequency [cm-1]                               :
 :......................................................:

      I      J                  FI(I,J)       k(I,J)       K(I,J)

      1      1                  500.00000      1.00000      0.10000
      2      2                  300.00000      0.30000      0.03000
      3      3                 -100.00000     -0.10000     -0.01000

 Num. of 2nd derivatives larger than  0.100D-05: 3 over 6

 ........................................................
 :        CUBIC FORCE CONSTANTS IN NORMAL MODES         :
 :                                                      :
 : FI =  Reduced values [cm-1]  (default input)         :
 :......................................................:

      I      J      K          FI(I,J,K)     k(I,J,K)     K(I,J,K)

      2      2      2           -30.00000     -0.30000     -0.03000
      3      3      3            30.00000      0.30000      0.03000

 Num. of 3rd derivatives larger than  0.100D-05: 2 over 10

 ........................................................
 :                                                      :
 :       QUARTIC FORCE CONSTANTS IN NORMAL MODES        :
 :                                                      :
 : FI =  Reduced values [cm-1]  (default input)         :
 :......................................................:

      I      J      K      L  FI(I,J,K,L)   k(I,J,K,L)   K(I,J,K,L)

      1      1      1      1   -120.00000    -1.20000     -0.12000
      2      2      2      2   -120.00000    -1.20000     -0.12000
      3      3      3      3     50.00000     0.50000      0.05000

 Num. of 4th derivatives larger than  0.100D-05: 3 over 15
EOF

"$PYTHON" scripts/mw_path_dvr.py \
  --gaussian-log "$TMP_DIR/anharmonic_modes.log" \
  --anharmonic-mode 1 \
  --well-type double \
  --solver sinc-dvr \
  --grid 81 \
  --levels 3 \
  --save-states 2 \
  --outdir "$TMP_DIR/out_anharm_ts" \
  --figdir "$TMP_DIR/figs_anharm_ts" \
  --prefix anharm_ts

test ! -f "$TMP_DIR/out_anharm_ts/anharm_ts_anharmonic_vpt2_comparison.csv"

"$PYTHON" scripts/mw_path_dvr.py \
  --gaussian-log "$TMP_DIR/anharmonic_modes.log" \
  --anharmonic-mode 2 \
  --well-type single \
  --solver sinc-dvr \
  --grid 81 \
  --levels 3 \
  --save-states 2 \
  --outdir "$TMP_DIR/out_anharm_morse" \
  --figdir "$TMP_DIR/figs_anharm_morse" \
  --prefix anharm_morse

grep -q "handy_morse_reference" "$TMP_DIR/out_anharm_morse/anharm_morse_anharmonic_summary.txt"
test -f "$TMP_DIR/out_anharm_morse/anharm_morse_anharmonic_vpt2_comparison.csv"

"$PYTHON" scripts/mw_path_dvr.py \
  --gaussian-log "$TMP_DIR/anharmonic_modes.log" \
  --anharmonic-mode 3 \
  --well-type single \
  --solver sinc-dvr \
  --grid 81 \
  --levels 3 \
  --save-states 2 \
  --outdir "$TMP_DIR/out_anharm_gaussian" \
  --figdir "$TMP_DIR/figs_anharm_gaussian" \
  --prefix anharm_gaussian

grep -q "handy_gaussian_reference" "$TMP_DIR/out_anharm_gaussian/anharm_gaussian_anharmonic_summary.txt"
test -f "$TMP_DIR/out_anharm_gaussian/anharm_gaussian_anharmonic_vpt2_comparison.csv"
grep -q "variational_transition_cm-1" "$TMP_DIR/out_anharm_gaussian/anharm_gaussian_anharmonic_vpt2_comparison.csv"
! grep -Eq "morse_.*variational" "$TMP_DIR/out_anharm_gaussian/anharm_gaussian_anharmonic_vpt2_comparison.csv"

cat > "$TMP_DIR/tail_path.xyz" <<'EOF'
2
point=1
H 0.000000 0.000000 -0.350000
H 0.000000 0.000000  0.350000
2
point=2
H 0.000000 0.000000 -0.370000
H 0.000000 0.000000  0.370000
2
point=3
H 0.000000 0.000000 -0.395000
H 0.000000 0.000000  0.395000
2
point=4
H 0.000000 0.000000 -0.425000
H 0.000000 0.000000  0.425000
2
point=5
H 0.000000 0.000000 -0.460000
H 0.000000 0.000000  0.460000
EOF

cat > "$TMP_DIR/tail_props.csv" <<'EOF'
point,energy_cm-1,B_MHz
1,100.0,10000.0
2,55.0,9980.0
3,15.0,9940.0
4,0.0,9890.0
5,30.0,9820.0
EOF

cat > "$TMP_DIR/property_derivatives.csv" <<'EOF'
property,value,d1,d2,d3,origin,parity
odd_shift,10.0,2.0,,0.1,zero,auto
even_quad,5.0,,0.5,,zero,auto
EOF

"$PYTHON" scripts/mw_path_dvr.py \
  --xyz "$TMP_DIR/tail_path.xyz" \
  --properties-csv "$TMP_DIR/tail_props.csv" \
  --property-derivatives-csv "$TMP_DIR/property_derivatives.csv" \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --property B_MHz \
  --property odd_shift \
  --property even_quad \
  --temperature 50 \
  --temperature 300 \
  --boundary nonperiodic \
  --solver sinc-dvr \
  --grid 81 \
  --levels 3 \
  --save-states 2 \
  --path-symmetry half-even \
  --well-type double \
  --potential-extension repulsive-polynomial \
  --extension-degree 6 \
  --extension-length-au 10 \
  --extension-target-cm 1000 \
  --extension-points 8 \
  --potential-smoothing spline \
  --potential-spline-smoothing 0.0 \
  --outdir "$TMP_DIR/out_tail" \
  --figdir "$TMP_DIR/figs_tail" \
  --prefix tail_1d

"$PYTHON" scripts/mw_path_dvr.py \
  --xyz "$TMP_DIR/tail_path.xyz" \
  --properties-csv "$TMP_DIR/tail_props.csv" \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --property B_MHz \
  --boundary periodic \
  --solver fourier \
  --grid 80 \
  --levels 3 \
  --save-states 2 \
  --path-symmetry half-even-origin \
  --well-type double \
  --outdir "$TMP_DIR/out_periodic_sym" \
  --figdir "$TMP_DIR/figs_periodic_sym" \
  --prefix periodic_sym_1d

cat > "$TMP_DIR/single_path.xyz" <<'EOF'
2
point=1
H 0.000000 0.000000 -0.350000
H 0.000000 0.000000  0.350000
2
point=2
H 0.000000 0.000000 -0.370000
H 0.000000 0.000000  0.370000
2
point=3
H 0.000000 0.000000 -0.395000
H 0.000000 0.000000  0.395000
2
point=4
H 0.000000 0.000000 -0.425000
H 0.000000 0.000000  0.425000
2
point=5
H 0.000000 0.000000 -0.460000
H 0.000000 0.000000  0.460000
2
point=6
H 0.000000 0.000000 -0.500000
H 0.000000 0.000000  0.500000
2
point=7
H 0.000000 0.000000 -0.545000
H 0.000000 0.000000  0.545000
EOF

cat > "$TMP_DIR/single_props.csv" <<'EOF'
point,energy_cm-1,B_MHz
1,120.0,10000.0
2,45.0,9980.0
3,8.0,9940.0
4,0.0,9890.0
5,8.0,9820.0
6,25.0,9740.0
7,45.0,9650.0
EOF

"$PYTHON" scripts/mw_path_dvr.py \
  --xyz "$TMP_DIR/single_path.xyz" \
  --properties-csv "$TMP_DIR/single_props.csv" \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --property B_MHz \
  --boundary nonperiodic \
  --solver sinc-dvr \
  --grid 81 \
  --levels 3 \
  --save-states 2 \
  --well-type single \
  --potential-extension single-morse \
  --extension-length-au 8 \
  --extension-points 6 \
  --outdir "$TMP_DIR/out_single" \
  --figdir "$TMP_DIR/figs_single" \
  --prefix single_1d

cat > "$TMP_DIR/asym_path.xyz" <<'EOF'
2
point=1
H 0.000000 0.000000 -0.240000
H 0.000000 0.000000  0.240000
2
point=2
H 0.000000 0.000000 -0.280000
H 0.000000 0.000000  0.280000
2
point=3
H 0.000000 0.000000 -0.320000
H 0.000000 0.000000  0.320000
2
point=4
H 0.000000 0.000000 -0.360000
H 0.000000 0.000000  0.360000
2
point=5
H 0.000000 0.000000 -0.400000
H 0.000000 0.000000  0.400000
2
point=6
H 0.000000 0.000000 -0.440000
H 0.000000 0.000000  0.440000
2
point=7
H 0.000000 0.000000 -0.480000
H 0.000000 0.000000  0.480000
2
point=8
H 0.000000 0.000000 -0.520000
H 0.000000 0.000000  0.520000
2
point=9
H 0.000000 0.000000 -0.560000
H 0.000000 0.000000  0.560000
2
point=10
H 0.000000 0.000000 -0.600000
H 0.000000 0.000000  0.600000
2
point=11
H 0.000000 0.000000 -0.640000
H 0.000000 0.000000  0.640000
EOF

cat > "$TMP_DIR/asym_props.csv" <<'EOF'
point,energy_cm-1,B_MHz
1,22.4,10060.0
2,10.1,10030.0
3,12.0,10000.0
4,46.0,9970.0
5,81.8,9940.0
6,57.3,9900.0
7,12.9,9860.0
8,0.0,9820.0
9,6.9,9780.0
10,19.8,9740.0
11,36.4,9700.0
EOF

"$PYTHON" scripts/mw_path_dvr.py \
  --xyz "$TMP_DIR/asym_path.xyz" \
  --properties-csv "$TMP_DIR/asym_props.csv" \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --property B_MHz \
  --boundary nonperiodic \
  --solver sinc-dvr \
  --grid 91 \
  --levels 3 \
  --save-states 2 \
  --well-type double \
  --core-model asymmetric-parabola-gaussian \
  --potential-extension repulsive-polynomial \
  --extension-degree 8 \
  --extension-length-au 10 \
  --extension-target-cm 1000 \
  --extension-points 8 \
  --outdir "$TMP_DIR/out_asym" \
  --figdir "$TMP_DIR/figs_asym" \
  --prefix asym_1d

"$PYTHON" scripts/mw_path_dvr.py \
  --grid2d-csv examples/gaussian_outputs/metric_test_grid.csv \
  --grid2d-geom-xyz examples/xyz/metric_test_grid.xyz \
  --q1-key q1 --q2-key q2 \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --boundary1 nonperiodic --boundary2 nonperiodic \
  --metric-mode geometry \
  --metric-stencil 5 \
  --metric-smoothing auto \
  --solver1 sinc-dvr --solver2 sinc-dvr \
  --basis1 3 --basis2 3 \
  --levels 3 \
  --outdir "$TMP_DIR/out_geometry" \
  --figdir "$TMP_DIR/figs_geometry" \
  --prefix metric_geometry

"$PYTHON" scripts/mw_path_dvr.py \
  --grid2d-csv "$TMP_DIR/out_geometry/metric_geometry_2d_grid.csv" \
  --q1-key q1 --q2-key q2 \
  --energy-key relative_energy_cm-1 --energy-unit cm-1 \
  --boundary1 nonperiodic --boundary2 nonperiodic \
  --metric-mode csv \
  --metric-smoothing none \
  --g11-key g11 --g12-key g12 --g22-key g22 \
  --solver1 sinc-dvr --solver2 sinc-dvr \
  --basis1 3 --basis2 3 \
  --levels 3 \
  --outdir "$TMP_DIR/out_csv" \
  --figdir "$TMP_DIR/figs_csv" \
  --prefix metric_csv

"$PYTHON" scripts/mw_path_dvr.py \
  --grid2d-csv examples/gaussian_outputs/metric_test_grid.csv \
  --q1-key q1 --q2-key q2 \
  --energy-key energy_cm-1 --energy-unit cm-1 \
  --boundary1 nonperiodic --boundary2 nonperiodic \
  --metric-mode constant \
  --solver1 sinc-dvr --solver2 sinc-dvr \
  --basis1 3 --basis2 3 \
  --levels 3 \
  --outdir "$TMP_DIR/out_constant" \
  --figdir "$TMP_DIR/figs_constant" \
  --prefix metric_constant

test -s "$TMP_DIR/out_geometry/metric_geometry_2d_levels.csv"
test -s "$TMP_DIR/out_csv/metric_csv_2d_levels.csv"
test -s "$TMP_DIR/out_constant/metric_constant_2d_levels.csv"
test -s "$TMP_DIR/out_tail/tail_1d_levels.csv"
test -s "$TMP_DIR/out_tail/tail_1d_model_profile.csv"
test -s "$TMP_DIR/out_tail/tail_1d_thermal_expectations.csv"
test -s "$TMP_DIR/out_periodic_sym/periodic_sym_1d_levels.csv"
test -s "$TMP_DIR/out_single/single_1d_levels.csv"
test -s "$TMP_DIR/out_asym/asym_1d_levels.csv"
grep -q "left_tail_min_outward_step_cm-1" "$TMP_DIR/out_tail/tail_1d_summary.txt"
grep -q "double_well_criterion_status: ok" "$TMP_DIR/out_tail/tail_1d_summary.txt"
grep -q "property_odd_shift_expectation_rule: central_value_for_symmetric_potential" "$TMP_DIR/out_tail/tail_1d_summary.txt"
"$PYTHON" - "$TMP_DIR/out_tail/tail_1d_thermal_expectations.csv" <<'PY'
import csv
import sys

rows = list(csv.DictReader(open(sys.argv[1], newline="")))
temperatures = {float(row["temperature_K"]) for row in rows}
assert {0.0, 50.0, 300.0}.issubset(temperatures)
for row in rows:
    assert abs(float(row["odd_shift"]) - 10.0) < 1.0e-10
assert any(abs(float(row["even_quad"]) - 5.0) > 1.0e-6 for row in rows)
PY
grep -q "Boundary: periodic" "$TMP_DIR/out_periodic_sym/periodic_sym_1d_summary.txt"
grep -q "symmetry_reference_point: first" "$TMP_DIR/out_periodic_sym/periodic_sym_1d_summary.txt"
grep -q "single_morse_fit_rms_cm-1" "$TMP_DIR/out_single/single_1d_summary.txt"
grep -q "core_model_used: asymmetric-parabola-gaussian" "$TMP_DIR/out_asym/asym_1d_summary.txt"
grep -q "asym_pg_fit_rms_cm-1" "$TMP_DIR/out_asym/asym_1d_summary.txt"

echo "Smoke tests passed. Temporary outputs: $TMP_DIR"
