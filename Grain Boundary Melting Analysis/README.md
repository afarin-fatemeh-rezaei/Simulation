# Grain Boundary Premelting in 2D Lennard‑Jones Bicrystals  
**Part 1: Effect of Misorientation Angle**

This repository contains the codes and results for studying grain‑boundary premelting in a two‑dimensional triangular Lennard‑Jones crystal.  
In this first part, we investigate how the misorientation angle of a symmetric tilt grain boundary affects the onset of premelting and the disordering temperature.

---

## Folder Structure

```
.
├── convert.py                     # Converts .json structure files to .dat format
├── gb_gen.c                       # Bicrystal generation & relaxation (C, OpenMP)
├── gb_multi.c                     # Melting MD simulation for multiple input files (C, OpenMP)
├── generate_perfect_grain.py      # Generates perfect single crystal (triangular lattice)
├── part_one_results/              # Output plots (threshold crossing, excess disorder)
│   ├── analysis.py                    # Plots one melting curve vs perfect crystal
│   ├── analysis2.py                   # Threshold crossing & excess disorder analysis
│   ├── threshold_crossing_threshold0.5.png
│   ├── threshold_crossing_threshold0.6.png
│   ├── threshold_crossing_threshold0.7.png
│   ├── threshold_crossing_threshold0.8.png
│   ├── excess_disorder_threshold0.5.png
│   ├── excess_disorder_threshold0.6.png
│   ├── excess_disorder_threshold0.7.png
│   ├── excess_disorder_threshold0.8.png
│   ├── atoms_0deg.dat                   # Perfect single crystal (text format)
│   ├── si_scan_0deg_best.json           # Perfect single crystal (JSON)
│   ├── si_scan_0deg_best.png            # Visualisation of perfect crystal
│   ├── atoms_10deg.dat                  # Relaxed 10° bicrystal (DAT)
│   ├── bicrystal_scan_10deg_best.json   # Relaxed 10° bicrystal (JSON)
│   ├── bicrystal_scan_10deg_best.png    # Visualisation of 10° bicrystal
│   ├── atoms_15deg.dat                  # Relaxed 15° bicrystal (DAT)
│   ├── bicrystal_scan_15deg_best.json   # Relaxed 15° bicrystal (JSON)
│   ├── ...
│   ├── 0_3_psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 3 
│   ├── 0_4_psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 4 
│   ├── 0_5_psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 5 
│   ├── 0_6psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 6 
│   ├── 0_7_psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 7 
│   ├── 0_40_psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 40 
│   ├── 0_50_psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 50 
│   ├── 0_60_psi_GB_vs_T.csv              # Output of gb_multi for perfect single crystal with seed 60 
│   ├── 10_3_psi_GB_vs_T.csv              # Output of gb_multi for Relaxed 10° bicrystal with seed 3 
│   ├── 10_4_psi_GB_vs_T.csv              # Output of gb_multi for Relaxed 10° bicrystal with seed 4 
│   └── ...
```

---

## File Descriptions

### Structure Generation

| File | Description |
|------|-------------|
| `generate_perfect_grain.py` | Python script that generates a perfect triangular lattice in a rectangular box with periodic boundary conditions in \(y\). The box size is commensurate with the lattice. |
| `gb_gen.c` | C code for generating symmetric tilt bicrystals. It rotates two grains by ±θ/2, removes overlapping atoms, relaxes the structure via gradient descent, and searches for the optimal translation between grains. |
| `convert.py` | Converts structure files from `.json` format to `.dat` format (or vice versa) for use with the C codes. |

### Melting Simulation

| File | Description |
|------|-------------|
| `gb_multi.c` | Molecular dynamics code that heats a structure from \(T=0.05\) to \(T=1.0\) and outputs the average \(|\psi_6|\) in the grain‑boundary region (\(|x| < 3.0\)) at each temperature. It uses a shifted‑force Lennard‑Jones potential, velocity Verlet integration, Berendsen thermostat, and cell lists. The random seed is fixed as `SEED` in the code. |

### Analysis Scripts

| File | Description |
|------|-------------|
| `analysis.py` | Plots the \(\psi_6(T)\) curve of a single run against the perfect crystal reference. Useful for quick visual checks. |
| `analysis2.py` | Reads the CSV outputs from multiple seeds and calculates two quantities for each angle: (1) the threshold crossing temperature (where \(|\psi_6|\) first drops below a threshold), and (2) the maximum excess disorder \(\Delta\psi_6 = \psi_6^{\text{perfect}} - \psi_6^{\text{bicrystal}}\). It also computes error bars via bootstrap resampling and saves the plots. |

---

## Prerequisites

- **C compiler** with OpenMP support (e.g., GCC)
- **Python 3** with `numpy`, `pandas`, `matplotlib`
- On Linux/macOS: `gcc`, `python3`
- On Windows: MinGW‑w64 or similar

---

## Compilation

Compile the C codes with optimisation and OpenMP:

```bash
gcc -O3 -march=native -ffast-math -funroll-loops -fopenmp -o gb_gen gb_gen.c -lm
gcc -O3 -march=native -ffast-math -funroll-loops -fopenmp -o gb_multi gb_multi.c -lm
```

If you do not have OpenMP, remove `-fopenmp`. The code will still run, but slower.

---

## Usage

### 1. Generate the Perfect Single Crystal

Run:

```bash
python generate_perfect_grain.py --a 1.2 --theta 0 --xmin -120 --xmax 120 --ymin -36.373 --ymax 36.373 --output atoms_0deg.json --txt atoms_0deg.dat
```

This creates `atoms_0deg.json` and `atoms_0deg.dat`.

### 2. Generate a Bicrystal for a Desired Angle

Edit the angle in `gb_gen.c`:

```c
double theta_deg = 10.0;   // change to 10, 15, 20, 25, 30
```

Compile and run:

```bash
./gb_gen
```

The code scans translations, relaxes the structure, and outputs:

- `bicrystal_scan_<angle>deg_best.json`
- `best_structure.txt` (overwritten each time, so rename if needed)

### 3. Convert JSON to DAT (if needed)

```bash
python convert.py input.json output.dat
```

*(The exact usage of `convert.py` may vary; check the script.)*

### 4. Run Melting Simulations

The `gb_multi` code accepts multiple input files and processes them sequentially.  
For each angle, you must run the melting simulation with **8 different random seeds** to average out the effect of initial velocities.  
The seed is defined at the top of `gb_multi.c`:

```c
#define SEED        50
```

Change this value, recompile, and run. For example, to run the 10° bicrystal with seed 50:

```bash
gcc ... -o gb_multi gb_multi.c -lm
./gb_multi bicrystal_scan_10deg_best.dat
```

The output will be:

```
bicrystal_scan_10deg_best_seed50_psi_GB_vs_T.csv
```

Repeat for seeds 51, 52, …, 57 (or any 8 values) and for each angle.

### 5. Analyse the Results

Once you have all CSV files (for each angle and seed), run:

```bash
python analysis2.py
```

Make sure the script’s file‑name patterns match your files. The script will:

- Load all CSV files for each angle and each seed.
- Compute average \(\psi_6(T)\) curves.
- For thresholds \(0.8, 0.7, 0.6, 0.5\), find:
  - **Threshold crossing temperature** (first \(T\) where \(|\psi_6|\) drops below the threshold).
  - **Maximum excess disorder** \(\Delta\psi_6 = \psi_6^{\text{single}} - \psi_6^{\text{bicrystal}}\).
- Compute error bars via bootstrap resampling across seeds.
- Save plots in the current directory.

The plots will be named like:

```
threshold_crossing_with_bootstrap_errors.png
excess_disorder_with_bootstrap_errors.png
```

Move them to `part_one_results/` if desired.

---

## Parameters Used

| Parameter | Value |
|-----------|-------|
| Lattice constant \(a\) | \(1.2\) |
| Box \(x\)-range | \([-120, 120]\) |
| Box \(y\)-range | \([-36.373, 36.373]\) |
| Lennard‑Jones \(\sigma\) (generation/relaxation) | \(1.069 = a / 2^{1/6}\) |
| Lennard‑Jones \(\sigma\) (melting) | \(1.069\) (fixed to match) |
| Lennard‑Jones \(\varepsilon\) | \(1.0\) |
| Cutoff radius \(r_c\) | \(3.2\) |
| Timestep \(\Delta t\) | \(0.005\) |
| Equilibration steps \(N_{\text{eq}}\) | \(1000\) |
| Production steps \(N_{\text{prod}}\) | \(2000\) |
| Temperature range | \(0.05\) – \(1.0\) with step \(0.01\) |
| Thermostat time constant \(\tau\) | \(0.1\) |
| GB probe region | \(|x| < 3.0\) |
| Neighbour cutoff for \(\psi_6\) | \(1.8\) |
| Number of seeds | 8 per angle |

---

## Output Files

- **CSV files** from `gb_multi`: `T,psi_GB` columns, where `psi_GB` is the average \(|\psi_6|\) in the GB probe.
- **Plots** from `analysis2.py`: threshold crossing temperature and excess disorder vs misorientation angle, with error bars.
- **Structure files**: `.dat` and `.json` for each angle and for the perfect crystal.

---

## Notes

- The reflecting walls in \(x\) are **hard walls**, not free surfaces. This is a simplification; results are qualitative.
- The Berendsen thermostat does not produce a true canonical ensemble; it is used for simplicity.
- The simulation times per temperature are short (\(15\) time units); near the transition, longer runs would be needed for quantitative accuracy.
- The grain‑boundary region \(|x| < 3.0\) is measured in simulation length units, not lattice constants.
- The scripts assume all CSV files are in the same directory as the scripts. Adjust paths if needed.

---

## License

This project is licensed under the MIT License – see the LICENSE file for details.
