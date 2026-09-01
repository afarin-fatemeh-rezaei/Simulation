# Equilibrium and Melting: Molecular Dynamics Simulations

**Lennard-Jones molecular dynamics simulations for thermal equilibration and solid-liquid phase transitions**

---

## 📋 Table of Contents

- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [The Lennard-Jones Potential](#the-lennard-jones-potential)
- [Simulation 1: System Reaching Equilibrium](#simulation-1-system-reaching-equilibrium)
  - [Objectives](#objectives)
  - [Methodology](#methodology)
  - [Key Results](#key-results)
  - [Outputs](#outputs)
- [Simulation 2: System Melting](#simulation-2-system-melting)
  - [Objectives](#objectives-1)
  - [Methodology](#methodology-1)
  - [Key Results](#key-results-1)
  - [Outputs](#outputs-1)
- [Dependencies](#dependencies)
- [How to Run](#how-to-run)
- [References](#references)
- [Acknowledgments](#acknowledgments)

---

## 📖 Overview

This repository contains two molecular dynamics (MD) simulation projects from a **Simulation course**, exploring fundamental concepts in statistical mechanics and condensed matter physics:

1. **System Reaching Equilibrium** – Study of thermal equilibration and velocity distributions in a 2D Lennard-Jones fluid
2. **System Melting** – Investigation of the solid-liquid phase transition in a 2D Lennard-Jones system

Both simulations use **reduced Lennard-Jones units** (σ = ε = k_B = 1) and the **velocity Verlet integration algorithm**.

---

## 📁 Repository Structure

```
Simulation/
│
├── System Reaching Equilibrium.py          # Thermal equilibration simulation
├── System Melting.py                       # Melting transition simulation
│
├── equilibrium/                            # Outputs from equilibration simulation
│   ├── System Reaching Equilibrium outputs/
│   │   └── ... (output files from simulation)
│
├── melting/                                # Outputs from melting simulation
│   ├── System Melting outputs/
│   │   └── ... (output files from simulation)
│
└── README.md                               # This file
```

---

## 🔬 The Lennard-Jones Potential

Both simulations use the **Lennard-Jones (LJ) potential**, a simple model for pairwise interactions:

\[
V(r) = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6 \right]
\]

In **reduced units** (σ = ε = k_B = 1):

\[
V(r) = 4 \left( r^{-12} - r^{-6} \right)
\]

The corresponding force is:

\[
\mathbf{F} = 24 \left( 2r^{-14} - r^{-8} \right) \mathbf{r}
\]

### Simulation Parameters (Common to Both)

| Parameter | Value | Meaning |
|-----------|-------|---------|
| **N** | 16 | Number of particles |
| **dt** | 0.005 | Integration time step |
| **cutoff** | 3.0 | Potential cutoff radius |
| **L** | 8.0 or 10.0 | Box side length |

---

## 🧪 Simulation 1: System Reaching Equilibrium

**File**: `System Reaching Equilibrium.py`

### Objectives

1. Simulate a 2D Lennard-Jones fluid reaching thermal equilibrium
2. Verify the **Maxwell-Boltzmann velocity distribution**
3. Verify the **Maxwell speed distribution** in 2D
4. Calculate the **coarse-grained entropy** evolution
5. Study the approach to equilibrium using multiple independent runs

### Methodology

**Initial Conditions**:
- **Positions**: 4×4 square lattice with small random perturbations
- **Velocities**: Fixed magnitude (v₀ = 1), random directions, zero total momentum

**Protocol**:
- Run 10 independent simulations with different random seeds
- Each run: 40 time units (N_steps = 8000)
- Binning: 0.5 time units per bin
- Track: Temperature, velocity distributions, coarse-grained entropy

**Coarse-Grained Entropy**:

The 2D velocity space is binned into a 20×20 grid. The entropy is calculated as:

\[
S_{cg} = -\sum_i p_i \ln p_i
\]

where \(p_i\) is the probability of finding a particle in velocity bin \(i\).

### Key Results

#### Temperature Equilibration

The system starts at a low temperature (from cold initial conditions) and equilibrates to a final temperature determined by the LJ interactions.

- **Running average** over 10 runs shows convergence
- **Error bars** (±1 standard deviation) indicate statistical fluctuations

#### Velocity Distributions

The \(v_x\) distribution follows a **Gaussian**:

\[
P(v_x) = \frac{1}{\sqrt{2\pi T}} \exp\left(-\frac{v_x^2}{2T}\right)
\]

The speed distribution follows the **2D Maxwell distribution**:

\[
P(v) = \frac{v}{T} \exp\left(-\frac{v^2}{2T}\right)
\]

#### Entropy Evolution

The coarse-grained entropy increases monotonically and saturates at equilibrium, consistent with the **Second Law of Thermodynamics**.

### Outputs

The outputs are stored in `equilibrium/System Reaching Equilibrium outputs/` 

---

## 🧊 Simulation 2: System Melting

**File**: `System Melting.py`

### Objectives

1. Simulate a 2D Lennard-Jones solid (triangular lattice)
2. Study the **solid-liquid phase transition** via heating
3. Track the **mean squared displacement (MSD)** as a melting criterion
4. Monitor **energy components** (kinetic, potential, total)
5. Visualize the **structural transformation** from solid to liquid

### Methodology

**Initial Conditions**:
- **Positions**: 4×4 square lattice with small random perturbations
- **Velocities**: Zero (cold start)

**Equilibration Phase (Steps 0-4000)**:
- System relaxes from square lattice to **triangular solid**
- Vibrations around equilibrium positions

**Heating Phase (Steps 0-6000)**:
- Gradual heating by scaling velocities every 50 steps
- Scale factor = 1.02 (2% increase)
- System eventually melts

**Heating Protocol**:
```
Every 50 steps:
    v → 1.02 × v
```

### Key Results

#### Structural Transformation

| Stage | Description | Snapshot Times |
|-------|-------------|----------------|
| **Initial** | Square lattice (ordered) | t = 0 |
| **Solid** | Triangular lattice with vibrational trails | t = 0.5, 1.0, 2.0 |
| **Melting** | Structural breakdown begins | t = 4.0, 7.5 |
| **Liquid** | Disordered, diffusive motion | t = 15.0, 25.0 |

#### Mean Squared Displacement (MSD)

The MSD is a key indicator of melting:

\[
\text{MSD}(t) = \langle |\mathbf{r}_i(t) - \mathbf{r}_i(0)|^2 \rangle
\]

- **Solid**: MSD plateaus (particles vibrate around fixed positions)
- **Liquid**: MSD grows linearly with time (diffusive motion)
- **Transition**: MSD shows a sharp increase at the melting point

#### Energy Evolution

| Energy Component | Behavior |
|------------------|----------|
| **Kinetic** | Increases with temperature (heating) |
| **Potential** | Decreases (becomes less negative) as structure breaks down |
| **Total** | Increases due to heating |

#### Snapshot Visualization

The melting process is visualized through:
1. **Vibrational trails** – Particle trajectories over time
2. **Snapshots** – Instantaneous configurations at key times
3. **Nearest-neighbor bonds** – Structural connectivity

### Outputs

The outputs are stored in `melting/System Melting outputs/` and include:


---

## 📊 Comparison of the Two Simulations

| Feature | System Reaching Equilibrium | System Melting |
|---------|----------------------------|----------------|
| **Initial State** | Square lattice, finite velocity | Square lattice, zero velocity |
| **Goal** | Thermal equilibration | Solid-liquid phase transition |
| **Heating** | None (system thermalizes internally) | Gradual velocity scaling |
| **Key Quantity** | Velocity distributions, entropy | MSD, energy components |
| **Final State** | Equilibrium fluid | Liquid |
| **Simulation Time** | 40 time units | ~50 time units |
| **Ensemble** | Multiple runs (10) | Single run |

---

## 📦 Dependencies

### Python
- Python 3.x
- NumPy
- Matplotlib

### Installation

```bash
pip install numpy matplotlib
```

---

## 🚀 How to Run

### System Reaching Equilibrium

```bash
python "System Reaching Equilibrium.py"
```

The script will:
1. Run 10 independent simulations
2. Update the plots interactively after each run
3. Display:
   - Temperature equilibration (with error bars)
   - \(v_x\) distribution (Gaussian fit)
   - Speed distribution (2D Maxwell fit)
   - Coarse-grained entropy evolution

### System Melting

```bash
python "System Melting.py"
```

The script will:
1. Equilibrate the system to a solid state
2. Record solid-state trajectories (500 steps)
3. Gradually heat the system
4. Generate and save:
   - `lattice_transformation.png`
   - `melting_snapshots.png`
   - `msd_melting.png`
   - `energy_temperature.png`

---

## 🔑 Key Concepts

### Statistical Mechanics
- **Equipartition theorem**: \(T = \langle E_{kin} \rangle / N\)
- **Maxwell-Boltzmann distribution**: Velocity distribution at equilibrium
- **Entropy**: \(S = -k_B \sum p_i \ln p_i\)
- **Second Law**: Entropy increases monotonically to equilibrium

### Molecular Dynamics
- **Velocity Verlet integration**: Symplectic, time-reversible
- **Periodic boundary conditions**: Minimum image convention
- **Lennard-Jones potential**: Model for noble gases
- **Reduced units**: σ = ε = k_B = 1

### Phase Transitions
- **Mean squared displacement**: Criterion for solid vs. liquid
- **Lindemann criterion**: Melting when MSD exceeds ~10% of nearest-neighbor distance
- **Structural disorder**: Loss of long-range order

---

## 📚 References

1. **Allen, M. P., & Tildesley, D. J.** (2017). *Computer Simulation of Liquids*, 2nd Edition. Oxford University Press.

2. **Frenkel, D., & Smit, B.** (2002). *Understanding Molecular Simulation: From Algorithms to Applications*, 2nd Edition. Academic Press.

3. **Rapaport, D. C.** (2004). *The Art of Molecular Dynamics Simulation*, 2nd Edition. Cambridge University Press.

4. **Hansen, J.-P., & McDonald, I. R.** (2013). *Theory of Simple Liquids*, 4th Edition. Academic Press.

5. **Landau, R. H., Páez, M. J., & Bordeianu, C. C.** (2015). *Computational Physics: Problem Solving with Python*, 3rd Edition. Wiley-VCH.

---

## 🙏 Acknowledgments

This repository contains exercises for the **Simulation** course, exploring fundamental computational methods in statistical mechanics and molecular dynamics.

---

**Happy Simulating!** ⚛️💻
