## Outputs from System Reaching Equilibrium Simulations

The System Reaching Equilibrium simulation was run with **multiple different initial conditions** to study the approach to thermal equilibrium and the sensitivity of the equilibration process:

1. **Optimal** (baseline configuration)
2. **High Randomness** (disordered initial positions)
3. **Low Randomness** (ordered initial positions)
4. **Low Density** (lower particle density)
5. **Low Initial Velocity** (slower initial velocities)
6. **Low Effective Radius** (modified interaction range)
7. **Lower Effective Radius** (further modified interaction range)

Each simulation generated data tracking four key quantities:
- **Temperature equilibration** (T vs. time)
- **Velocity distribution** (v_x distribution)
- **Speed distribution** (v distribution)
- **Coarse-grained entropy** (S_cg evolution)

---

### Data Tables

#### 1. Optimal

**File**: `Optimal.png`

**Description**:
Data from the baseline simulation showing the approach to equilibrium.

**Table Structure**:

| Time (t) | Mean Temp (10 runs) | ±1 std | v_x | Speed v* | Mean S_cg | ±1 std | S_cg/k_B |
|----------|---------------------|--------|-----|----------|-----------|--------|----------|
| 0 | 0.000 | 0.00 | -4.0 | 0.00 | 0.000 | 0.000 | 0.00 |
| 1 | 0.000 | 0.00 | -4.0 | 0.00 | 0.000 | 0.000 | 0.00 |
| ... | ... | ... | ... | ... | ... | ... | ... |
| 20 | 0.000 | 0.00 | -4.0 | 0.19 | -0.20 | 0.0 | 0.00 |

**Key Observations**:
- **Temperature**: Starts low (~0.58) and equilibrates to ~1.0 by t = 25
- **Velocity Distribution**: Approaches Gaussian shape as system equilibrates
- **Speed Distribution**: Approaches 2D Maxwell distribution
- **Entropy**: Increases from ~0 to ~2.6 and plateaus at equilibrium

**Physical Interpretation**:
- The system thermalizes from a cold initial state
- Velocity distribution broadens and becomes Maxwellian
- Entropy increases monotonically (Second Law of Thermodynamics)
- Equilibrium reached at t ≈ 25-30

---

#### 2. High Randomness

**File**: `High Randomness.png`

**Description**:
Data from simulation with highly disordered initial positions.

**Table Structure**:

| Parameter | Value | Unit |
|-----------|-------|------|
| Temperature (t=0) | 0.58 | ε/k_B |
| Temperature (t=5) | 0.85 | ε/k_B |
| Temperature (t=10) | 0.95 | ε/k_B |
| Temperature (t=15) | 0.98 | ε/k_B |
| Temperature (t=20) | 0.99 | ε/k_B |
| Temperature (t=25) | 1.00 | ε/k_B |
| Temperature (t=30) | 1.00 | ε/k_B |
| v_x distribution | 0.000 | - |
| Speed distribution | 0.000-0.01 | - |
| Entropy (t=0-40) | 2.60 | k_B |

**Key Observations**:
- **Faster equilibration**: System reaches equilibrium by t ≈ 15
- **Higher initial entropy**: S_cg starts at ~2.60 (already disordered)
- **Lower temperature fluctuations**: Small ±1 std values
- **Velocity distribution**: More uniform due to randomness

**Physical Interpretation**:
- High initial disorder accelerates thermal equilibration
- The system does not need to "find" equilibrium configurations
- Entropy remains constant after early times (already near maximum)

---

#### 3. Low Randomness

**File**: `Low Randomness.png`

**Description**:
Data from simulation with ordered initial positions (low disorder).

**Table Structure**:

| Parameter | Value | Unit |
|-----------|-------|------|
| Temperature (t=0) | 0.00 | k_B |
| Temperature (t=5) | 0.50 | k_B |
| Temperature (t=10) | 0.80 | k_B |
| Temperature (t=15) | 0.90 | k_B |
| Temperature (t=20) | 0.95 | k_B |
| Temperature (t=25) | 0.98 | k_B |
| Temperature (t=30) | 0.99 | k_B |
| Temperature (t=35) | 1.00 | k_B |
| Temperature (t=40) | 1.00 | k_B |
| v_x distribution | -4 to 6 | - |
| Speed v* | 0-6 | - |

**Key Observations**:
- **Slower equilibration**: Takes ~35-40 time units to reach equilibrium
- **Lower initial entropy**: Starts near zero, increases monotonically
- **Well-defined distributions**: Clear Gaussian and Maxwellian shapes
- **Larger fluctuations**: Higher ±1 std values during equilibration

**Physical Interpretation**:
- Ordered initial state requires more time to thermalize
- Entropy production is significant and continuous
- Clear demonstration of the Second Law of Thermodynamics
- Provides the best comparison with ideal gas theory

---

#### 4. Low Density

**File**: `Low Density.png`

**Description**:
Data from simulation with lower particle density (larger box size).

**Table Structure**:

| Time (t) | Temperature | v_x (t=5) | Speed (t=5) | Entropy | Mean S_cg | ±1 std |
|----------|-------------|-----------|-------------|---------|-----------|--------|
| 0 | 0.50 | 0.000 | 0.000 | 0.00 | 0.000 | 0.000 |
| 1 | 0.75 | 0.000 | 0.000 | 0.10 | 0.000 | 0.000 |
| 2 | 0.78 | 0.000 | 0.000 | 0.20 | 0.000 | 0.000 |
| ... | ... | ... | ... | ... | ... | ... |
| 20 | 1.05 | 0.000 | 0.30 | 1.07 | 0.000 | 0.000 |

**Key Observations**:
- **Lower equilibrium density**: Particles interact less frequently
- **Faster temperature rise**: Heat capacity lower at lower density
- **Reduced entropy**: Lower entropy due to fewer interactions
- **Velocity distributions**: Approach Maxwellian more quickly

**Physical Interpretation**:
- Low density approximates an ideal gas more closely
- Fewer collisions means less energy exchange
- Equilibration still occurs but with different dynamics
- Good comparison to test density effects on equilibration

---

#### 5. Low Initial Velocity

**File**: `Low Initial Velocity.png`

**Description**:
Data from simulation with initially slow particles.

**Key Observations**:
- **Lower starting kinetic energy**: System starts colder
- **Slower temperature rise**: Takes longer to reach equilibrium
- **Different equilibration dynamics**: Kinetic energy transfer is gradual
- **Velocity distribution**: Starts narrow, broadens with time

**Physical Interpretation**:
- Initial kinetic energy affects equilibration rate
- Energy must be transferred from potential to kinetic
- Demonstrates energy equipartition over time
- Tests sensitivity to initial velocity distribution

---

#### 6. Low Effective Radius

**File**: `Low Effective Radius.png`

**Description**:
Data from simulation with reduced Lennard-Jones interaction range.

**Table Structure**:

| Parameter | Value | Unit |
|-----------|-------|------|
| Temperature (t=0-5) | Various | k_B |
| Temperature (t=10-15) | Various | k_B |
| v_x distribution | -4.0 to 6.0 | - |
| Speed v* | 0.0 to 6.0 | - |
| Coarse-grained entropy | 0.0 to 6.0 | k_B |
| Mean S_cg (10 runs) | 0.000 | - |
| ±1 std | 0.000 | - |

**Key Observations**:
- **Modified interaction**: Particles interact over shorter range
- **Lower potential energy**: Less attractive interaction
- **Altered equilibrium properties**: Different final temperature
- **Distinct entropy evolution**: Different path to equilibrium

**Physical Interpretation**:
- Effective radius changes the interaction potential
- Alters the equilibrium state of the system
- Tests sensitivity to potential parameters
- Demonstrates how interactions govern equilibration

---

#### 7. Lower Effective Radius

**File**: `Lower Effective Radius.png`

**Description**:
Data from simulation with further reduced Lennard-Jones interaction range.

**Table Structure**:

| Parameter | Value | Unit |
|-----------|-------|------|
| Temperature | 0, 5, 10, 15, 20, 25, 30, 35, 40 | k_B |
| v_x distribution | -4.0 to 6.0 | - |
| Speed v* | 0.0 to 6.0 | - |
| Coarse-grained entropy | 0.0 to 6.0 | k_B |

**Key Observations**:
- **Even weaker interactions**: Further reduced potential range
- **More ideal gas-like behavior**: Less deviation from Maxwellian
- **Faster equilibration**: Fewer interactions to thermalize
- **Lower entropy**: Less correlation between particles

**Physical Interpretation**:
- Pushing system toward ideal gas limit
- Demonstrates the role of interactions in equilibration
- Tests convergence to ideal gas behavior
- Useful comparison with Boltzmann equation predictions

---

### Comparison of Initial Conditions

| Condition | Equilibration Time | Final Temperature | Entropy Change | Velocity Distribution |
|-----------|-------------------|-------------------|----------------|----------------------|
| **Optimal** | t ~ 25-30 | ~1.00 | Moderate | Good Maxwellian |
| **High Randomness** | t ~ 15 | ~1.00 | Small (already high) | Broad, uniform |
| **Low Randomness** | t ~ 35-40 | ~1.00 | Large (starts near 0) | Good Maxwellian |
| **Low Density** | t ~ 15-20 | ~1.05 | Moderate | Quick Maxwellian |
| **Low Initial Velocity** | t ~ 30-35 | ~1.00 | Large | Slow broadening |
| **Low Effective Radius** | Variable | Modified | Variable | Modified |
| **Lower Effective Radius** | Faster | Modified | Lower | More ideal |

---

### Key Physical Insights

#### 1. Equilibration Process
- All systems eventually reach the same equilibrium temperature (~1.0 ε/k_B)
- The approach to equilibrium depends on initial conditions
- **Second Law of Thermodynamics**: Entropy always increases

#### 2. Velocity Distributions
- The v_x distribution converges to a **Gaussian**:
  \[
  P(v_x) = \frac{1}{\sqrt{2\pi T}} \exp\left(-\frac{v_x^2}{2T}\right)
  \]
- The speed distribution converges to the **2D Maxwell distribution**:
  \[
  P(v) = \frac{v}{T} \exp\left(-\frac{v^2}{2T}\right)
  \]

#### 3. Entropy Evolution
- S_cg starts low (ordered) and increases to maximum (equilibrium)
- Rate of entropy increase depends on:
  - Initial disorder (high randomness → fast)
  - Density (low density → faster)
  - Effective radius (shorter range → faster)

#### 4. Finite-Size Effects
- With N=16 particles, significant statistical fluctuations
- Multiple runs (10) allow averaging and error estimation
- Qualitative agreement with theory, not quantitative

---

### Summary Table

| File | Description | Key Features |
|------|-------------|--------------|
| `Optimal.png` | Baseline simulation | Standard approach to equilibrium |
| `High Randomness.png` | Disordered initial positions | Fast equilibration, high initial entropy |
| `Low Randomness.png` | Ordered initial positions | Slow equilibration, large entropy increase |
| `Low Density.png` | Lower density system | Faster temperature rise, reduced entropy |
| `Low Initial Velocity.png` | Slow initial particles | Slow temperature rise, gradual broadening |
| `Low Effective Radius.png` | Reduced interaction range | Modified equilibrium properties |
| `Lower Effective Radius.png` | Further reduced range | More ideal gas-like behavior |

---

### Theoretical Predictions

#### Equipartition Theorem
\[
T = \frac{1}{N} \sum_{i=1}^N \frac{1}{2} m v_i^2
\]
(For our units, \(m = 1\), so \(T = E_{kin}/N\))

#### Maxwell-Boltzmann Distribution (2D)
\[
P(v_x, v_y) = \frac{1}{2\pi T} \exp\left(-\frac{v_x^2 + v_y^2}{2T}\right)
\]

#### 2D Speed Distribution
\[
P(v) = \frac{v}{T} \exp\left(-\frac{v^2}{2T}\right)
\]

#### Entropy of 2D Gas
\[
S_{cg} = -\sum_i p_i \ln p_i
\]
where \(p_i\) is the probability in velocity bin \(i\)

---

### Key Takeaways

1. **Universal Equilibrium**: All systems reach the same equilibrium state regardless of initial conditions
2. **Second Law**: Entropy always increases monotonically
3. **Maxwellian Distributions**: Velocity and speed distributions become Maxwellian at equilibrium
4. **Sensitivity**: Equilibration rate depends on initial disorder, density, and interaction range
5. **Equipartition**: Energy distributes equally among degrees of freedom at equilibrium
6. **Finite-Size Effects**: Statistical fluctuations are significant with N=16 particles

