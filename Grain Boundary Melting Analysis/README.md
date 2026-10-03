# `README.md`

Save at `Grain Boundary Melting Analysis/README.md`.

```markdown
# Grain Boundary Melting Analysis

Molecular-dynamics study of grain-boundary premelting in a 2D Lennard-Jones
triangular crystal. We measure the grain-boundary (GB) energy as a function
of misorientation angle for a series of coincidence-site-lattice (CSL)
boundaries, and then quantify the degree to which each boundary disorders
before the bulk crystal melts.

**Main result.** Every boundary in the series disorders at least
$\Delta T \approx 0.02$ below the bulk melting temperature, and the
strength of the premelting signature decreases monotonically with decreasing
misorientation angle — from $\approx 0.038$ excess energy per atom at
$32.2^\circ$ down to $\approx 0.006$ at $4.8^\circ$.

---

## Table of contents

1. [Physical setup](#physical-setup)
2. [Theory](#theory)
   - [Lennard-Jones potential](#lennard-jones-potential)
   - [CSL misorientations in 2D](#csl-misorientations-in-2d)
   - [Grain-boundary energy](#grain-boundary-energy)
   - [Premelting](#premelting)
   - [Observables](#observables)
3. [Repository layout](#repository-layout)
4. [Pipeline](#pipeline)
5. [Results](#results)
6. [Requirements](#requirements)

---

## Physical setup

We study a two-dimensional Lennard-Jones system whose ground state is a
triangular lattice. A bicrystal is constructed by joining two semi-infinite
crystals rotated by an angle $\theta$ about a common axis perpendicular to
the plane. All simulations use a shifted-force cutoff at $r_c = 3.0\,\sigma$
so that the potential and force both go smoothly to zero at the cutoff.

Simulation parameters:

| Quantity | Value |
|---|---|
| Lattice constant $a_0$ | determined at $P = 0$ from a 40×40 perfect crystal |
| Reference energy $e_0$ | per-atom energy of the perfect crystal at $a_0$ |
| Cutoff $r_c$ | $3.0\,\sigma$ |
| Temperature range | $T = 0.28$ to $0.44$ (reduced units) |
| Thermostat | BAOAB Langevin, $\gamma = 1.0$ |
| Barostat | Berendsen, $P_{\rm target} = 0$, $\kappa = 10^{-2}$ |
| Timestep | $\Delta t = 0.005$ |
| Production length | $5 \times 10^5$ steps per temperature |

---

## Theory

### Lennard-Jones potential

The pairwise interaction is the standard LJ potential truncated and shifted
in the force (shifted-force scheme):

$$
U(r) = 4\epsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right],
$$

with the shifted-force correction

$$
U_{\rm SF}(r) = U(r) - U(r_c) + F(r_c)\,(r - r_c),
$$

$$
F_{\rm SF}(r) = -\frac{dU}{dr}\bigg|_{r} + F(r_c).
$$

The force is continuous and vanishes at $r = r_c$, so no long-range tail
correction is needed and energy is conserved at the level required for
NVT/NPT MD. All energies below are in units of $\epsilon$, lengths in
units of $\sigma$, and temperature in units of $\epsilon/k_B$.

### CSL misorientations in 2D

For a triangular (hexagonal) lattice, the coincidence-site-lattice (CSL)
misorientations form a discrete set. Any pair of coprime integers $(m, n)$
with $0 \leq n < m$ generates one CSL orientation:

$$
\Sigma = m^2 + mn + n^2,
\qquad
\theta = 2\,\arctan\!\left(\frac{\sqrt{3}\,n}{2m + n}\right).
$$

The true coincidence index (after dividing by the common factor 3 that
appears when $\Sigma$ is divisible by 3) is

$$
\Sigma_{\rm true} = \begin{cases} \Sigma/3 & \text{if } \Sigma \equiv 0 \pmod 3 \\ \Sigma & \text{otherwise.} \end{cases}
$$

Low-$\Sigma$ values are "special" orientations with short periodicities
and low GB energies. In our series we use
$\Sigma \in \{7, 39, 91, 139, 147, 421\}$, corresponding to
$\theta$ from $4.84^\circ$ to $32.20^\circ$.

The bicrystal is generated in a rectangular box that is periodic under
the CSL lattice:

$$
L_x = n_x\, a_0 \sqrt{\Sigma}, \qquad
L_y = n_y\, a_0 \sqrt{3\,\Sigma},
$$

with integer multipliers $n_x, n_y$ chosen so that (i) the box is at
least a few times larger than the lattice-dislocation spacing

$$
D = \frac{a_0}{\theta_{\rm rad}}
$$

at low angles, and (ii) the boundary does not interact with itself
across the periodic cell. The atom count is roughly
$N \approx 2\, n_x\, n_y\, \Sigma$.

### Grain-boundary energy

The grain-boundary energy per unit length is defined as the excess energy
of the bicrystal relative to the same number of atoms in the perfect
crystal at the same $a_0$:

$$
\sigma_{\rm GB} = \frac{E_{\rm bicrystal} - N\, e_0}{2\, L_y}.
$$

The factor of 2 accounts for the two GBs in the periodic cell. This
quantity is computed at $T = 0$ by minimising the potential energy with
FIRE, after a coarse-to-fine search over the rigid-body translations
$(s_x, s_y)$ that describe how the two grains slide relative to each
other. The minimum over translations is the equilibrium GB energy.

**Low-angle regime.** Below roughly $15^\circ$ the boundary is a periodic
array of lattice dislocations with spacing $D = a_0/\theta_{\rm rad}$.
The classical Read–Shockley result gives

$$
\sigma_{\rm GB}(\theta) = \sigma_0\,\theta\,(A - \ln\theta),
$$

which vanishes as $\theta \to 0$. This is the origin of the sharp drop
of $\sigma_{\rm GB}$ at low angle, and of the plateau that develops
above $15$–$20^\circ$ where the dislocation cores overlap and the
concept of isolated dislocations breaks down.

**Cusp structure.** Superimposed on this smooth curve are cusps at
every CSL angle: at $\Sigma = 7$, $\Sigma = 13$, and so on, the
boundary energy has a local minimum because the two lattices can adopt
a particularly commensurate, low-energy atomic arrangement. The depth
of each cusp scales roughly as $1/\Sigma$. These cusps are visible in
`sigma_cusps.png`.

### Premelting

Grain-boundary premelting is the phenomenon in which the boundary region
disorders and behaves liquid-like at temperatures **below** the bulk
melting point $T_m$ of the crystal. The physical picture is:

1. The GB core has broken bonds and excess free volume, so it has a
   lower local melting temperature than the bulk.
2. As $T$ approaches $T_m$ from below, the GB region becomes
   progressively more disordered and eventually forms a thin liquid-like
   film — a *premelting film*.
3. The film thickens as $T \to T_m$. Above $T_m$ the whole crystal melts.

Thermodynamically, GB premelting manifests as a melting-point depression
$\Delta T = T_m^{\rm bulk} - T_m^{\rm GB}$ and as a positive GB excess
energy $\Delta\mathrm{PE}(T)$ that grows with $T$.

In our simulations we do not measure $T_m$ directly. Instead we use the
layer-resolved bond-order parameter $\psi_6$ (see below) to detect the
temperature at which the bicrystal bulk disorders, and we compare it
to the same quantity in a perfect crystal of the same size and shape.
The difference $\Delta T$ between the two is the melting-point depression,
and it is ≈ 0.02 in reduced units for every angle in our series.

### Observables

**Bond-orientational order parameter.** For each atom $i$ we compute
the hexagonal Steinhardt parameter

$$
\psi_6^{(i)} = \left| \frac{1}{n_i} \sum_{j \in \mathcal{N}(i)} e^{6 i \theta_{ij}} \right|,
$$

where $\mathcal{N}(i)$ is the set of neighbours of $i$ within a shell
$[0.8\,a_{\rm nn}, 1.2\,a_{\rm nn}]$ and $\theta_{ij}$ is the angle of
the bond $\mathbf{r}_j - \mathbf{r}_i$. In a triangular crystal
$\psi_6 \to 1$; in a 2D liquid $\psi_6 \to 0$. The box-averaged value
$\langle \psi_6 \rangle$ is the standard structural indicator of
melting and is used to locate $T_m$ in each simulation box.

**Layer-resolved profiles.** The box is divided into `N_BINS = 60`
strips along $x$. Within each strip we accumulate the time-averaged
$\psi_6$, the coordination number $n_i$, and the local density
$\rho(x)$. This gives a profile $\psi_6(x)$ that shows the GB as a
local minimum of the order parameter and the bulk as a plateau.

**Grain-boundary and bulk bins.** For each temperature and seed, the
GB bin is identified as the minimum of the *smoothed* $\psi_6(x)$
profile (5-bin moving average to suppress per-bin noise). The bulk bin
is the maximum, at least 10 bins away (circularly) from the GB. The
same bin indices are used to evaluate the bulk reference in the
perfect crystal, so $\psi_6^{\rm GB}$ and $\psi_6^{\rm bulk}$ compare
the *same physical region* in the two boxes.

**GB excess energy.** The energy observable used throughout is

$$
\Delta\mathrm{PE}(T) = \frac{E_{\rm bicrystal}(T)}{N} - \frac{E_{\rm perfect}(T)}{N},
$$

the difference in potential energy per atom between the bicrystal and
the perfect crystal at the same temperature. It is positive because
the GB has broken bonds and strained coordination, and it grows with
$T$ because the boundary disorders progressively. $\Delta\mathrm{PE}$
does not depend on any bin identification — it is a whole-box
quantity — and so it is the most robust observable for the
angle-dependence study.

**GB structural disorder.** As a complementary observable we use

$$
d\psi_6^{\rm GB}(T) = \psi_6^{\rm GB,B}(T) - \psi_6^{\rm GB,P}(T),
$$

the difference in GB-bin order parameter between the bicrystal and the
perfect crystal at the same location. Since a perfect crystal has no
GB, this quantity measures how much the GB has lowered the local order
parameter relative to the bulk. It is negative and grows in magnitude
as $T \to T_m$.

---

## Repository layout

```
Grain Boundary Melting Analysis/
│
├── results/
│   ├── 4/                         θ = 4.84°,  Σ421
│   │   └── 5/                     (only one seed for this angle)
│   │       ├── energy_bicrystal_Sigma421_theta4.8381deg_T*_seed*.dat
│   │       ├── energy_perfect_Sigma421_theta4.8381deg_T*_seed*.dat
│   │       ├── profile_bicrystal_Sigma421_theta4.8381deg_T*_seed*.dat
│   │       ├── profile_perfect_Sigma421_theta4.8381deg_T*_seed*.dat
│   │       ├── rho_bicrystal_Sigma421_theta4.8381deg_T*_seed*.dat
│   │       └── rho_perfect_Sigma421_theta4.8381deg_T*_seed*.dat
│   │
│   ├── 10/                        θ = 10.42°, Σ91
│   │   ├── 5/
│   │   ├── 17/
│   │   └── 29/
│   │
│   ├── 16/                        θ = 16.43°, Σ147
│   │   ├── 5/  ├── 17/  └── 29/
│   ├── 21/                        θ = 21.79°, Σ7
│   │   ├── 5/  ├── 17/  └── 29/
│   ├── 25/                        θ = 25.46°, Σ139
│   │   ├── 5/  ├── 17/  └── 29/
│   ├── 32/                        θ = 32.20°, Σ39
│   │   ├── 5/  ├── 17/  └── 29/
│   │
│   ├── sweep_perfect_Sigma21_theta21.7868deg_vs_bicrystal_..._seed5.dat
│   ├── sweep_perfect_Sigma21_theta21.7868deg_vs_bicrystal_..._seed17.dat
│   ├── sweep_perfect_Sigma21_theta21.7868deg_vs_bicrystal_..._seed29.dat
│   ├── sweep_perfect_Sigma39_theta32.2042deg_vs_bicrystal_..._seed5.dat
│   ├── ...                                       (16 sweep files total)
│   │
│   ├── premelting_plots.py
│   ├── premelting_summary.txt
│   ├── fig1_melting_curves.png
│   ├── fig2_dPE_vs_T.png
│   ├── fig3_dpsi6_vs_T.png
│   ├── fig4_angle_dependence.png
│   └── fig5_dPE_per_angle.png
│
├── bicrystal.c                    Generate a bicrystal at a given CSL
├── perfect_grain.c                Generate a perfect triangular crystal
├── melt_scan.c                    Stage 2: temperature sweep for one angle
├── sigma_gb.c                     Stage 1: GB energy vs angle
│
├── bicrystal_Sigma21_theta21.7868deg.dat
├── bicrystal_Sigma39_theta32.2042deg.dat
├── bicrystal_Sigma91_theta10.4174deg.dat
├── bicrystal_Sigma139_theta25.4611deg.dat
├── bicrystal_Sigma147_theta16.4264deg.dat
├── bicrystal_Sigma421_theta4.8381deg.dat
│
├── perfect_Sigma21_theta21.7868deg.dat
├── perfect_Sigma39_theta32.2042deg.dat
├── perfect_Sigma91_theta10.4174deg.dat
├── perfect_Sigma139_theta25.4611deg.dat
├── perfect_Sigma147_theta16.4264deg.dat
├── perfect_Sigma421_theta4.8381deg.dat
│
├── check_structure.py
├── sigma_plot.py
├── sigma_gb_csl.csv
├── sigma_vs_theta.png
├── sigma_vs_sigma.png
├── sigma_cusps.png
└── readshockley_fit.png
```

### Per-temperature files inside each seed folder

For each seed folder (e.g. `results/21/17/`) the temperature sweep
$T = 0.28 \ldots 0.44$ with $\Delta T = 0.02$ (9 points) produces three
file families, each in a `_bicrystal_` and a `_perfect_` variant:

| File pattern | Columns |
|---|---|
| `energy_*_T<T>_seed<S>.dat` | step, PE, KE, $T_{\rm inst}$, coordination, $\psi_6$, $\psi_4$, P |
| `profile_*_T<T>_seed<S>.dat` | x, $\rho$, coordination, $\psi_6$, $\psi_4$ per bin |
| `rho_*_T<T>_seed<S>.dat` | x, $\rho$ |

The `_bicrystal_` and `_perfect_` variants are the two boxes compared
at each $T$ to compute the GB excess energy $\Delta\mathrm{PE}$.

---

## Pipeline

### Stage 1 — GB energy vs angle

`sigma_gb.c` generates each CSL bicrystal at $a = a_0$, removes
overlaps, and relaxes the atomic positions with FIRE. The GB energy is
computed via the definition above after a coarse-to-fine search over
rigid-body shifts $(s_x, s_y)$. An adaptive box sizing based on the
dislocation spacing $D = a_0 / \theta_{\rm rad}$ ensures the boundary
does not interact with itself across the periodic cell at low angles.

```bash
gcc -O2 -fopenmp -o sigma_gb sigma_gb.c -lm
OMP_NUM_THREADS=4 ./sigma_gb
python sigma_plot.py
```

**Output:** `sigma_gb_csl.csv` with columns
`theta_deg, sigma_gb, Sigma, Sigma_true, m, n, nx, ny, N_atoms,
best_sx_a, best_sy_a, pe_final, converged, dmin`.

### Stage 2 — Premelting scan

`melt_scan.c` reads a `perfect_Sigma*.dat` and a matching
`bicrystal_Sigma*.dat`, and runs a temperature sweep. At each $T$:

1. Runs the perfect and bicrystal boxes in parallel via OpenMP.
2. Uses a Berendsen barostat at $P_{\rm target} = 0$ to hold pressure.
3. Accumulates a layer-resolved $\psi_6$ profile over the production
   window.
4. Identifies the GB bin as the minimum of the smoothed $\psi_6$
   profile, and the bulk bin as the maximum, at least 10 bins away.
5. Reports $\psi_6$, $\psi_4$, coordination and density in the GB and
   bulk bins, and the GB excess energy $\Delta\mathrm{PE}(T)$.

```bash
gcc -O2 -fopenmp -o melt_scan melt_scan.c -lm
OMP_NUM_THREADS=4 ./melt_scan \
    perfect_Sigma21_theta21.7868deg.dat \
    bicrystal_Sigma21_theta21.7868deg.dat \
    0.28 0.44 0.02 500000 1.0 0.005 5 0.0 1
```

Arguments after the two file names:
`T_min T_max dT n_steps gamma dt seed P_target use_barostat`.

Each run writes the per-$T$ energy / profile / rho files into the
current directory and one `sweep_*.dat` summary at the end.

### Stage 3 — Analysis and figures

```bash
cd results
python premelting_plots.py
```

`premelting_plots.py` auto-discovers every `sweep_*.dat` in `results/`,
groups them by misorientation angle, averages across seeds, and writes
the five figures and the summary table.

| Figure | Content |
|---|---|
| `fig1_melting_curves.png` | $\langle\psi_6\rangle$ vs $T$, one subplot per angle |
| `fig2_dPE_vs_T.png` | $\Delta\mathrm{PE}$ vs $T$ for all six angles |
| `fig3_dpsi6_vs_T.png` | $d\psi_6^{\rm GB}$ vs $T$ |
| `fig4_angle_dependence.png` | **Headline.** $\Delta\mathrm{PE}$ and $d\psi_6^{\rm GB}$ vs $\theta$ at $T = 0.28$ |
| `fig5_dPE_per_angle.png` | Per-angle PE curves with $\Delta\mathrm{PE}$ overlaid |
| `premelting_summary.txt` | Mean ± std of the two observables at $T = 0.28, 0.32, 0.36$ |

---

## Results

### Melting-point depression

The box-averaged $\langle \psi_6 \rangle$ curve drops sharply when the
crystal melts. In every seed and every angle:

- the **bicrystal bulk** disorders between $T = 0.38$ and $T = 0.40$;
- the **perfect crystal** remains solid at $T = 0.40$ and disorders
  between $T = 0.40$ and $T = 0.42$.

The resulting melting-point depression is

$$
\Delta T \;=\; T_m^{\rm perfect} - T_m^{\rm bicrystal} \;\approx\; 0.02
$$

in reduced units, at every misorientation angle in the series. This is
the primary thermodynamic signature of GB premelting.

### Angle dependence

Both observables at $T = 0.28$ increase monotonically with
misorientation angle:

| $\theta$ (deg) | $\Sigma$ | seeds | $\Delta\mathrm{PE}$ | $d\psi_6^{\rm GB}$ |
|---:|---:|---:|---:|---:|
| 4.84  | 421 | 1 | 0.0058 | −0.013 |
| 10.42 | 91  | 3 | 0.0185 | −0.040 |
| 16.43 | 147 | 3 | 0.0304 | −0.067 |
| 21.79 | 7   | 3 | 0.0354 | −0.083 |
| 25.46 | 139 | 3 | 0.0354 | −0.093 |
| 32.20 | 39  | 3 | 0.0376 | −0.096 |

Both quantities vary by a factor of $\sim 6$–$7$ across the series.
The bulk-reference bin in the bicrystal matches the same bin in the
perfect crystal to within $\pm 0.01$ in $\psi_6$ at every angle and
temperature, confirming that the bulk bin is uncontaminated by the GB.

### Interpretation

Low-angle boundaries are arrays of lattice dislocations with narrow,
ordered cores and little free volume. As the angle increases, the
dislocation spacing $D = a_0/\theta_{\rm rad}$ shrinks, the cores
overlap, and the boundary becomes a broader disordered region with more
free volume and more broken bonds. This is reflected in:

- **larger GB excess energy** $\Delta\mathrm{PE}$ at high angle;
- **stronger structural disorder** $d\psi_6^{\rm GB}$ at high angle;
- **stronger premelting response** overall at high angle.

The 4.84° ($\Sigma 421$) boundary sits at the edge of detectability
within the signal-to-noise of our runs, while the 32.20° ($\Sigma 39$)
boundary gives the strongest signal. The monotonic trend across six
angles over a factor of six in strength, seen independently in both
a whole-box energy observable and a local structural observable, is
the central quantitative result of the study.

---

## Requirements

- C compiler with OpenMP support (`gcc`, `clang`)
- Python 3.8+
- Python packages: `numpy`, `matplotlib`

## Compilation

```bash
gcc -O2 -fopenmp -o sigma_gb      sigma_gb.c      -lm
gcc -O2 -fopenmp -o melt_scan     melt_scan.c     -lm
gcc -O2 -fopenmp -o bicrystal    bicrystal.c     -lm
gcc -O2          -o perfect_grain perfect_grain.c -lm
```

