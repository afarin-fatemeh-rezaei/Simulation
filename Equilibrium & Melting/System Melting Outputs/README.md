## Outputs from System Melting Simulations

The System Melting simulation was run with **three different initial conditions** to study the sensitivity of the melting transition:

1. **Optimal** (baseline configuration)
2. **High Energy** (increased initial kinetic energy)
3. **High Initial Velocity** (higher initial velocities)
4. **High Randomness** (more disordered initial positions)
5. **Low Randomness** (more ordered initial positions)

Each simulation generated four diagnostic plots that reveal different aspects of the melting transition.

---

### Figure Sets

#### 1. Structure Plots: Initial Square Lattice → Triangular Solid

**Files**: `*_Structure.png` (e.g., `Optimal_Structure.png`, `HighEnergy_Structure.png`)

**Description**:
- **Left panel**: Initial 4×4 square lattice with nearest-neighbor bonds (gray lines)
- **Right panel**: Equilibrated solid configuration showing:
  - Vibrational trails (gray trajectories)
  - Time-averaged positions (blue points)
  - Nearest-neighbor bonds between averaged positions

**Physical Interpretation**:
- The system undergoes a structural transformation from a metastable square lattice to the ground-state **triangular lattice**
- Vibrational trails show particles oscillating around equilibrium positions
- The triangular lattice maximizes packing density for the Lennard-Jones potential

**Key Observations**:
- All initial conditions converge to the same triangular solid structure
- The degree of disorder in the initial state affects the equilibration time
- Vibrational amplitudes are larger at higher temperatures

---

#### 2. Melting Snapshots: Time-lapse with Vibrational Trails

**Files**: `*_Melting.png` (e.g., `Optimal_Melting.png`, `HighEnergy_Melting.png`)

**Description**:
Eight snapshots showing the system evolution during heating:

| Step | Time (t) | Stage |
|------|----------|-------|
| 0 | 0.00 | Initial square lattice (ordered) |
| 100 | 0.50 | Transition to triangular solid |
| 200 | 1.00 | Solid with vibrational motion |
| 400 | 2.00 | Solid, increasing vibrational amplitude |
| 800 | 4.00 | Pre-melting: enhanced vibrations, some disorder |
| 1500 | 7.50 | Melting onset: structural breakdown |
| 3000 | 15.00 | Partial melting: disordered regions |
| 5000 | 25.00 | Fully melted: liquid state |

**Physical Interpretation**:
- **Gray trails**: Particle trajectories over the last 40 steps (vibrational/translational motion)
- **Red points**: Instantaneous positions at the snapshot time
- **Solid phase**: Particles vibrate around fixed lattice positions (trails are localized)
- **Liquid phase**: Particles diffuse across the box (trails are extended)

**Key Observations by Condition**:

| Condition | Behavior |
|-----------|----------|
| **Optimal** | Clear solid→liquid transition at t ~ 7.5-15.0 |
| **High Energy** | Melts earlier due to higher initial kinetic energy |
| **High Initial Velocity** | Rapid melting, less defined solid phase |
| **High Randomness** | Melts slightly earlier than optimal |
| **Low Randomness** | Similar to optimal, slightly delayed melting |

---

#### 3. Mean Squared Displacement (MSD) Plots

**Files**: `*_r^2.png` (e.g., `Optimal_r^2.png`, `HighEnergy_r^2.png`)

**Description**:
Evolution of the mean squared displacement during heating:

\[
\text{MSD}(t) = \langle |\mathbf{r}_i(t) - \mathbf{r}_i(0)|^2 \rangle
\]

**Typical Y-axis Range**: 0.0 to 3.5 (reduced units)
**Typical X-axis Range**: 0 to 30 (reduced time units)

**Physical Interpretation**:
- **Solid region (t < 7.5)**: MSD plateaus at a small value (~0.2-0.5)
  - Particles vibrate around fixed positions
  - No long-range diffusion
  
- **Melting region (t ~ 7.5-15.0)**: MSD shows a steep increase
  - Particles begin to escape their lattice sites
  - Structural breakdown
  
- **Liquid region (t > 15.0)**: MSD grows linearly with time
  - Particles diffuse freely
  - Slope relates to diffusion coefficient: \(D = \frac{1}{4}\frac{d\langle r^2\rangle}{dt}\)

**Key Observations by Condition**:

| Condition | MSD Behavior | Melting Time |
|-----------|--------------|--------------|
| **Optimal** | Clear plateau → steep rise → linear growth | t ~ 10 |
| **High Energy** | Higher plateau, earlier rise | t ~ 7 |
| **High Initial Velocity** | Minimal plateau, rapid rise | t ~ 5 |
| **High Randomness** | Similar to optimal, slightly earlier | t ~ 8 |
| **Low Randomness** | Clear plateau, delayed rise | t ~ 12 |

---

#### 4. Energy Evolution Plots

**Files**: `*_Energy.png` (e.g., `Optimal_Energy.png`, `HighEnergy_Energy.png`)

**Description**:
Thermodynamic evolution during the heating process, showing:

1. **Temperature** (\(T = E_{kin}/N\)) evolution
2. **Kinetic energy** (blue line)
3. **Potential energy** (green line)
4. **Total energy** (black line)

**Typical Y-axis Range**: 0 to 800 (energy in ε units)
**Typical X-axis Range**: 0 to 30 (reduced time units)

**Physical Interpretation**:

| Energy Component | Behavior | Physical Meaning |
|------------------|----------|------------------|
| **Kinetic** | Increases during heating | Temperature rises |
| **Potential** | Becomes less negative (increases) | Structure breaks down, particles escape potential wells |
| **Total** | Steadily increases | Energy input from velocity scaling |
| **Temperature** | Follows kinetic energy | Direct measure of particle motion |

**Key Observations by Condition**:

| Condition | Energy Behavior |
|-----------|-----------------|
| **Optimal** | Smooth temperature increase; potential energy rises at melting |
| **High Energy** | Higher initial kinetic energy; melts earlier |
| **High Initial Velocity** | Very high kinetic energy; rapid melting |
| **High Randomness** | Similar to optimal, slight variations |
| **Low Randomness** | Delayed melting, sharper energy changes |

---

### Summary: Effect of Initial Conditions

| Initial Condition | Structure Evolution | Melting Time (MSD) | Energy Behavior |
|-------------------|---------------------|-------------------|-----------------|
| **Optimal** | Smooth solid → liquid | t ~ 10 | Smooth, gradual |
| **High Energy** | Faster transition, less ordered solid | t ~ 7 | Higher initial energy |
| **High Initial Velocity** | Almost immediate melting | t ~ 5 | Very high kinetic energy |
| **High Randomness** | Delayed equilibration | t ~ 8 | Slight variations |
| **Low Randomness** | Well-defined solid, delayed melting | t ~ 12 | Sharp transition |

---

## 📊 Key Insights for the Melting Simulation

### 1. Lindemann Criterion
The melting transition occurs when the **MSD exceeds ~0.3-0.5** (approximately 10-15% of the nearest-neighbor distance), consistent with the Lindemann melting criterion.

### 2. First-Order Nature
The melting transition appears **first-order** based on:
- Sharp increase in MSD
- Discontinuity in potential energy
- Coexistence of solid and liquid regions during transition

### 3. System Sensitivity
The melting behavior is **sensitive to initial conditions**:
- High initial velocity can cause immediate melting
- Low randomness produces a more stable solid
- High randomness accelerates the transition

### 4. Finite-Size Effects
With only 16 particles, finite-size effects are significant:
- No true phase transition (requires thermodynamic limit)
- Surface effects dominate
- Results are qualitative rather than quantitative

---

## 📚 Additional Output Files

The `System Melting outputs/` folder contains:

| File | Type | Description |
|------|------|-------------|
| `Optimal_Structure.png` | PNG | Lattice transformation (square → triangular) |
| `Optimal_Melting.png` | PNG | Melting snapshots with vibrational trails |
| `Optimal_r^2.png` | PNG | Mean squared displacement evolution |
| `Optimal_Energy.png` | PNG | Energy and temperature evolution |
| `HighEnergy_*.png` | PNG | Same plots for high-energy initial conditions |
| `HighInitialVelocity_*.png` | PNG | Same plots for high initial velocity |
| `HighRandomness_*.png` | PNG | Same plots for high randomness |
| `LowRandomness_*.png` | PNG | Same plots for low randomness |

---

## 🔑 Key Takeaways

1. **Structural transformation**: Square lattice → triangular solid → liquid
2. **Melting signature**: MSD plateau → rapid increase → linear diffusion
3. **Thermodynamic signature**: Potential energy rise during melting
4. **Sensitivity**: Melting behavior depends on initial conditions
5. **Lindemann criterion**: MSD melting threshold at ~0.3-0.5
6. **Finite-size effects**: Qualitative behavior, not quantitative phase diagram

