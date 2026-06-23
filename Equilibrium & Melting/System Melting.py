import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

# ============================================================
# Simulation parameters (reduced units, carbon‑like)
# ============================================================
N = 16
L = 8.0          # box side → density = 1 particle/σ²
dt = 0.005
cutoff = 3.0
cutoff2 = cutoff**2

# ---------- helper functions ----------
def minimum_image(dr):
    """Minimum image convention using global L."""
    return dr - L * np.round(dr / L)

def apply_pbc(r):
    """Put a single particle back into the box."""
    return r - L * np.floor(r / L)

def compute_acc_and_pe(pos):
    """Lennard‑Jones acceleration (force, m=1) and potential energy."""
    acc = np.zeros_like(pos)
    pe = 0.0
    for i in range(N):
        for j in range(i+1, N):
            dr = minimum_image(pos[i] - pos[j])
            r2 = np.dot(dr, dr)
            if r2 < cutoff2 and r2 > 1e-12:
                r2i = 1.0 / r2
                r6i = r2i * r2i * r2i
                # Force: 24*(2/r^14 - 1/r^8) * dr
                ff = 24.0 * (2.0 * r6i * r6i - r6i) * r2i
                fij = ff * dr
                acc[i] += fij
                acc[j] -= fij
                # Potential: 4*(r^-12 - r^-6)
                pe += 4.0 * (r6i * r6i - r6i)
    return acc, pe

def compute_acc(pos):
    """Keep old interface for back‑compatibility (calls compute_acc_and_pe)."""
    acc, _ = compute_acc_and_pe(pos)
    return acc

def velocity_verlet(pos, vel, acc, dt):
    """One step of velocity‑Verlet integration."""
    vel_half = vel + 0.5 * dt * acc
    new_pos = pos + vel_half * dt
    new_pos = np.array([apply_pbc(r) for r in new_pos])
    new_acc, pe = compute_acc_and_pe(new_pos)
    new_vel = vel_half + 0.5 * dt * new_acc
    return new_pos, new_vel, new_acc, pe

# ---------- initial square lattice, zero velocity ----------
pos_init = np.array([[i+0.5, j+0.5] for i in range(4) for j in range(4)])
pos = pos_init + np.random.uniform(-0.05, 0.05, pos_init.shape)  # tiny noise
pos = np.array([apply_pbc(r) for r in pos])
vel = np.zeros((N, 2))
acc, pe = compute_acc_and_pe(pos)

# ---------- equilibration → triangular solid ----------
n_eq = 4000
traj_interval = 20
pos_traj_eq = []
for step in range(n_eq):
    if step % traj_interval == 0:
        pos_traj_eq.append(pos.copy())
    pos, vel, acc, pe = velocity_verlet(pos, vel, acc, dt)
pos_solid = pos.copy()
vel_solid = vel.copy()

# ---------- Record solid‑state trajectories (for vibration plot) ----------
n_record = 500
pos_record = np.zeros((n_record, N, 2))
pos = pos_solid.copy()
vel = vel_solid.copy()
acc, pe = compute_acc_and_pe(pos)
for k in range(n_record):
    pos_record[k] = pos.copy()
    pos, vel, acc, pe = velocity_verlet(pos, vel, acc, dt)

# ---------- Time‑averaged positions with PBC unwrapping ----------
ref = pos_record[0].copy()
unwrapped = np.zeros_like(pos_record)
unwrapped[0] = ref
for t in range(1, n_record):
    for i in range(N):
        dr = minimum_image(pos_record[t, i] - unwrapped[t-1, i])
        unwrapped[t, i] = unwrapped[t-1, i] + dr
pos_avg_unwrapped = unwrapped.mean(axis=0)
pos_avg = np.array([apply_pbc(r) for r in pos_avg_unwrapped])

# ---------- Quick neighbour check (instantaneous) ----------
print("Instantaneous neighbours (final solid snapshot):")
n_neigh = np.zeros(N, dtype=int)
for i in range(N):
    for j in range(N):
        if i == j: continue
        dr = minimum_image(pos_solid[i] - pos_solid[j])
        if np.dot(dr, dr) < 1.5**2:
            n_neigh[i] += 1
print("Neighbours per atom:", n_neigh)

# ---------- Melting: gradual heating ----------
n_heat = 6000
scale_every = 50
scale_factor = 1.02
ref_pos = pos_solid.copy()
msd = []
snapshot_times = [0, 100, 200, 400, 800, 1500, 3000, 5000]
traj_len = 40
snapshots_traj = {}

# Energy/temperature records
T_history = np.zeros(n_heat)
KE_history = np.zeros(n_heat)
PE_history = np.zeros(n_heat)
E_total_history = np.zeros(n_heat)

pos = pos_solid.copy()
vel = vel_solid.copy()
acc, pe = compute_acc_and_pe(pos)
rolling_buffer = []
for step in range(n_heat):
    rolling_buffer.append(pos.copy())
    if len(rolling_buffer) > traj_len:
        rolling_buffer.pop(0)
    if step in snapshot_times:
        snapshots_traj[step] = (np.array(rolling_buffer), pos.copy())
    
    # Thermodynamic measurements
    ke = 0.5 * np.sum(vel**2)
    T = ke / N               # reduced temperature, kB=1
    T_history[step] = T
    KE_history[step] = ke
    PE_history[step] = pe
    E_total_history[step] = ke + pe
    
    pos, vel, acc, pe = velocity_verlet(pos, vel, acc, dt)
    if step % scale_every == 0 and step > 0:
        vel *= scale_factor
    
    disp = np.array([minimum_image(pos[i] - ref_pos[i]) for i in range(N)])
    msd.append(np.mean(np.sum(disp**2, axis=1)))

time_heat = np.arange(n_heat) * dt

# ============================================================
# PLOTS
# ============================================================

# 1. Square lattice vs. solid with vibrations (two‑panel figure)
fig1, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5))

# --- Left: initial square lattice ---
ax1.scatter(pos_init[:,0], pos_init[:,1], s=50, c='blue', zorder=2)
segs = []
for i in range(N):
    for j in range(i+1, N):
        dr = minimum_image(pos_init[i] - pos_init[j])
        if np.dot(dr, dr) < 1.5**2:
            segs.append([pos_init[i], pos_init[j]])
lc = LineCollection(segs, colors='gray', linewidths=0.8, zorder=1)
ax1.add_collection(lc)
ax1.set_xlim(0, L); ax1.set_ylim(0, L)
ax1.set_aspect('equal')
ax1.set_title('Initial square lattice')
ax1.axis('off')

# --- Right: solid with vibrations ---
for i in range(N):
    ax2.plot(pos_record[:, i, 0], pos_record[:, i, 1], ',', color='gray', alpha=0.3)
ax2.scatter(pos_avg[:, 0], pos_avg[:, 1], s=40, c='blue', zorder=2)
for i in range(N):
    for j in range(i+1, N):
        dr = minimum_image(pos_avg[i] - pos_avg[j])
        if np.dot(dr, dr) < 1.5**2:
            ax2.plot([pos_avg[i,0], pos_avg[j,0]], [pos_avg[i,1], pos_avg[j,1]],
                     color='gray', lw=0.8)
ax2.set_xlim(0, L); ax2.set_ylim(0, L)
ax2.set_aspect('equal')
ax2.set_title('Solid: vibrations around triangular lattice')
ax2.axis('off')

fig1.tight_layout()
fig1.savefig('lattice_transformation.png', dpi=200)
print('Saved: lattice_transformation.png')

# 2. Melting snapshots with vibrational trails
fig2, axes = plt.subplots(2, 4, figsize=(12, 6))
axes = axes.flatten()
for idx, step in enumerate(snapshot_times):
    ax = axes[idx]
    traj, final_pos = snapshots_traj[step]
    for i in range(N):
        ax.plot(traj[:, i, 0], traj[:, i, 1], '-', color='gray', alpha=0.4, linewidth=0.5)
    ax.scatter(final_pos[:,0], final_pos[:,1], s=30, c='red', zorder=2)
    ax.set_xlim(0, L); ax.set_ylim(0, L)
    ax.set_aspect('equal')
    ax.set_title(f'step {step} (t = {step*dt:.2f})')
    ax.axis('off')
fig2.suptitle('Melting: time‑lapse snapshots with vibrational trails', fontsize=14)
fig2.tight_layout()
fig2.savefig('melting_snapshots.png', dpi=200)
print('Saved: melting_snapshots.png')

# 3. Mean‑squared displacement
fig3, ax = plt.subplots(figsize=(6,4))
ax.plot(time_heat, msd, 'b-', lw=1)
ax.set_xlabel('Time (reduced units)')
ax.set_ylabel('Mean squared displacement, ⟨Δr²⟩')
ax.set_title('Melting transition via MSD')
ax.grid(alpha=0.3)
fig3.tight_layout()
fig3.savefig('msd_melting.png', dpi=200)
print('Saved: msd_melting.png')

# 4. Energy and temperature during heating
fig4, (axT, axE) = plt.subplots(2, 1, figsize=(7, 6), sharex=True)

# Temperature
axT.plot(time_heat, T_history, 'r-', lw=1)
axT.set_ylabel('Temperature (ε/k_B)')
axT.set_title('Thermodynamic evolution during melting')
axT.grid(alpha=0.3)

# Energy components
axE.plot(time_heat, KE_history, 'b-', label='Kinetic', lw=1)
axE.plot(time_heat, PE_history, 'g-', label='Potential', lw=1)
axE.plot(time_heat, E_total_history, 'k-', label='Total', lw=1)
axE.set_xlabel('Time (reduced units)')
axE.set_ylabel('Energy (ε)')
axE.legend()
axE.grid(alpha=0.3)

fig4.tight_layout()
fig4.savefig('energy_temperature.png', dpi=200)
print('Saved: energy_temperature.png')

plt.show()