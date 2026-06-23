import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# Physical parameters (reduced units, dilute gas)
# ============================================================
N = 16
L = 10.0
dt = 0.005
cutoff = 3.0
t_total = 40.0
n_steps = int(t_total / dt)
bin_width = 0.5
steps_per_bin = int(bin_width / dt)
n_bins = int(t_total / bin_width)

v0 = 1
n_runs = 10

# Fixed 2D velocity bins for coarse‑grained entropy
v_bins_1d = np.linspace(-5, 5, 21)   # 20 bins in each direction
vx_edges = v_bins_1d
vy_edges = v_bins_1d

# ============================================================
# Lennard‑Jones force (correct)
# ============================================================
def compute_acc(pos):
    acc = np.zeros((N, 2))
    for i in range(N):
        for j in range(i+1, N):
            dr = pos[i] - pos[j]
            dr -= L * np.round(dr / L)          # minimum image
            r2 = np.dot(dr, dr)
            if r2 < cutoff*cutoff and r2 > 1e-12:
                r2i = 1.0 / r2
                r6i = r2i * r2i * r2i
                r12i = r6i * r6i
                inv_r = np.sqrt(r2i)
                # Force magnitude: 24 (2/r^13 - 1/r^7)
                f_mag = 24.0 * (2.0 * r12i * inv_r - r6i * inv_r)
                f = f_mag * dr
                acc[i] += f
                acc[j] -= f
    return acc

def verlet_step(pos, vel, acc, dt):
    vel_half = vel + 0.5 * dt * acc
    pos_new = pos + vel_half * dt
    pos_new %= L
    acc_new = compute_acc(pos_new)
    vel_new = vel_half + 0.5 * dt * acc_new
    return pos_new, vel_new, acc_new

# ============================================================
# Coarse‑grained entropy (Shannon) from 2D velocity histogram
# ============================================================
def coarse_grained_entropy(vel, vx_edges, vy_edges):
    # 2D histogram, normalized to probabilities
    H, xed, yed = np.histogram2d(vel[:,0], vel[:,1],
                                 bins=[vx_edges, vy_edges])
    # Flatten and compute probabilities
    p = H.ravel() / vel.shape[0]
    # Only non‑zero bins
    p = p[p > 0]
    return -np.sum(p * np.log(p))

# ============================================================
# Single simulation → binned T(t), S_cg(t), and velocity components
# ============================================================
def run_one(seed):
    rng = np.random.default_rng(seed)
    # initial positions: 4×4 grid, spacing 2.5σ + random 0.1σ
    grid = 1.5
    offset = (L - 3*grid) / 2.0
    pos = np.empty((N, 2))
    idx = 0
    for i in range(4):
        for j in range(4):
            pos[idx] = [offset + i*grid, offset + j*grid]
            idx += 1
    pos += rng.uniform(-0.5, 0.5, size=pos.shape)
    pos %= L

    # initial velocities: fixed magnitude, random direction
    angles = rng.uniform(0, 2*np.pi, N)
    vel = np.column_stack([v0*np.cos(angles), v0*np.sin(angles)])
    vel -= np.mean(vel, axis=0)   # zero total momentum

    acc = compute_acc(pos)

    bin_T, bin_S, bin_t = [], [], []
    temp_sum = 0.0
    S_sum = 0.0
    count = 0
    all_vx, all_vy = [], []

    for step in range(1, n_steps+1):
        ke = 0.5 * np.sum(vel**2)
        T = ke / N
        temp_sum += T

        # Coarse‑grained entropy at this instant
        S_inst = coarse_grained_entropy(vel, vx_edges, vy_edges)
        S_sum += S_inst
        count += 1

        for i in range(N):
            all_vx.append(vel[i,0])
            all_vy.append(vel[i,1])

        pos, vel, acc = verlet_step(pos, vel, acc, dt)

        if step % steps_per_bin == 0:
            bin_T.append(temp_sum / count)
            bin_S.append(S_sum / count)   # average entropy over the bin
            bin_t.append((step - steps_per_bin/2) * dt)
            temp_sum = 0.0
            S_sum = 0.0
            count = 0

    return (np.array(bin_t), np.array(bin_T), np.array(bin_S),
            np.array(all_vx), np.array(all_vy))

# ============================================================
# Interactive figure with 4 subplots (added entropy)
# ============================================================
plt.ion()
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 14))

all_T = []
all_S = []
all_t = None
vx_pool = []
vy_pool = []

# Axis labels
ax1.set_ylabel('Temperature (ε/k_B)')
ax1.set_title('Temperature equilibration (running average over runs)')
ax1.grid(alpha=0.3)

ax2.set_xlabel('v*_x')
ax2.set_ylabel('Probability density')
ax2.set_title('v_x distribution (pooled)')
ax2.grid(alpha=0.3)

ax3.set_xlabel('Speed v*')
ax3.set_ylabel('Probability density')
ax3.set_title('Speed distribution (pooled)')
ax3.grid(alpha=0.3)

ax4.set_xlabel('Time (reduced units)')
ax4.set_ylabel('S_cg / k_B')
ax4.set_title('Coarse‑grained entropy (2D velocity histogram)')
ax4.grid(alpha=0.3)

# Fixed binning for v_x and speed
vx_bins = np.linspace(-5, 5, 40)
speed_bins = np.linspace(0, 6, 40)

for run_idx in range(n_runs):
    seed = 42 + run_idx * 123
    t_arr, T_arr, S_arr, vx, vy = run_one(seed)

    if all_t is None:
        all_t = t_arr
    all_T.append(T_arr)
    all_S.append(S_arr)
    vx_pool.append(vx)
    vy_pool.append(vy)

    # Running averages and spread
    T_mean = np.mean(all_T, axis=0)
    T_std  = np.std(all_T, axis=0)
    S_mean = np.mean(all_S, axis=0)
    S_std  = np.std(all_S, axis=0)

    # Pool all velocities so far
    vx_all = np.concatenate(vx_pool)
    vy_all = np.concatenate(vy_pool)
    speeds_all = np.sqrt(vx_all**2 + vy_all**2)

    # Equilibrium temperature from last half of current mean
    half = len(T_mean)//2
    T_eq = np.mean(T_mean[half:])

    # ---- Update plots ----
    ax1.clear()
    ax1.plot(all_t, T_mean, 'b-', lw=1.5, label=f'Mean ({run_idx+1} runs)')
    ax1.fill_between(all_t, T_mean - T_std, T_mean + T_std,
                     alpha=0.3, label='±1 std')
    ax1.set_ylabel('Temperature (ε/k_B)')
    ax1.set_title('Temperature equilibration (running average)')
    ax1.legend()
    ax1.grid(alpha=0.3)

    ax2.clear()
    # v_x histogram
    hist_vx, _ = np.histogram(vx_all, bins=vx_bins, density=True)
    centers_vx = (vx_bins[:-1] + vx_bins[1:]) / 2
    ax2.plot(centers_vx, hist_vx, 'go', label='Simulation (pooled)')
    # Gaussian fit
    vx_th = np.linspace(-5, 5, 200)
    sigma = np.sqrt(T_eq)
    gauss = np.exp(-vx_th**2 / (2*sigma**2)) / np.sqrt(2*np.pi*sigma**2)
    ax2.plot(vx_th, gauss, 'r-', label=f'Gaussian (σ={sigma:.2f})')
    ax2.set_xlabel('v*_x')
    ax2.set_ylabel('Probability density')
    ax2.set_title('v_x distribution')
    ax2.legend()
    ax2.grid(alpha=0.3)

    ax3.clear()
    hist_sp, _ = np.histogram(speeds_all, bins=speed_bins, density=True)
    centers_sp = (speed_bins[:-1] + speed_bins[1:]) / 2
    ax3.plot(centers_sp, hist_sp, 'bo', label='Simulation (pooled)')
    v_th = np.linspace(0, 6, 200)
    p_th = (v_th / T_eq) * np.exp(-v_th**2 / (2*T_eq))
    p_th /= np.trapezoid(p_th, v_th)   # normalise
    ax3.plot(v_th, p_th, 'r-', label=f'2D Maxwell (T*={T_eq:.2f})')
    ax3.set_xlabel('Speed v*')
    ax3.set_ylabel('Probability density')
    ax3.set_title('Speed distribution')
    ax3.legend()
    ax3.grid(alpha=0.3)

    ax4.clear()
    ax4.plot(all_t, S_mean, 'm-', lw=1.5, label=f'Mean S_cg ({run_idx+1} runs)')
    ax4.fill_between(all_t, S_mean - S_std, S_mean + S_std,
                     alpha=0.3, label='±1 std')
    ax4.set_xlabel('Time (reduced units)')
    ax4.set_ylabel('S_cg / k_B')
    ax4.set_title('Coarse‑grained entropy')
    ax4.legend()
    ax4.grid(alpha=0.3)

    plt.tight_layout()
    plt.draw()
    plt.pause(0.5)

plt.ioff()
plt.show()

print(f"Final equilibrium temperature (averaged over {n_runs} runs): {T_eq:.3f}")