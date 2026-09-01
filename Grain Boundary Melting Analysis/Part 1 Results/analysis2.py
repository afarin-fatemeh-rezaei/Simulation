import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================
# Settings
# ============================================================
threshold = 0.8
T_min = 0.0
T_max = 1.00

angles = [0, 10, 15, 20, 25, 30]
seeds = [3, 4, 5, 6, 7, 40, 50, 60]
file_pattern = "{angle}_{seed}_psi_GB_vs_T.csv"

rng = np.random.default_rng(42)   # reproducible bootstrap
n_bootstrap = 200                 # number of bootstrap resamples

# ============================================================
# Load raw data
# ============================================================
raw = {}   # raw[angle][seed] = DataFrame with T, psi_GB

for angle in angles:
    raw[angle] = {}
    for seed in seeds:
        filename = file_pattern.format(angle=angle, seed=seed)
        if os.path.exists(filename):
            df = pd.read_csv(filename)
            df = df[(df["T"] >= T_min) & (df["T"] <= T_max)]
            df = df[["T", "psi_GB"]].sort_values("T").reset_index(drop=True)
            raw[angle][seed] = df
        else:
            print(f"Warning: {filename} not found")

# ============================================================
# Helper: average psi_GB across seeds for one angle
# ============================================================
def average_across_seeds(angle, seed_subset=None):
    if seed_subset is None:
        seed_subset = list(raw[angle].keys())

    dfs = []
    for seed in seed_subset:
        if seed in raw[angle]:
            dfs.append(raw[angle][seed])

    if not dfs:
        return None

    combined = pd.concat(dfs, ignore_index=True)
    avg = combined.groupby("T", as_index=False)["psi_GB"].mean()
    return avg

# ============================================================
# Helper: find first T where psi < threshold
# ============================================================
def find_crossing_T(df):
    below = df[df["psi_GB"] < threshold]
    if below.empty:
        return np.nan
    return below.iloc[0]["T"]

# ============================================================
# Helper: compute peak excess disorder for one angle
# ============================================================
def compute_excess_peak(angle, seed_subset=None):
    single_avg = average_across_seeds(0, seed_subset)
    if single_avg is None:
        return np.nan, np.nan

    bic_avg = average_across_seeds(angle, seed_subset)
    if bic_avg is None:
        return np.nan, np.nan

    merged = pd.merge(
        single_avg.rename(columns={"psi_GB": "psi_single"}),
        bic_avg.rename(columns={"psi_GB": "psi_bi"}),
        on="T",
        how="inner"
    )
    if merged.empty:
        return np.nan, np.nan

    merged["delta_psi6"] = merged["psi_single"] - merged["psi_bi"]
    peak_idx = merged["delta_psi6"].idxmax()
    return merged.loc[peak_idx, "T"], merged.loc[peak_idx, "delta_psi6"]

# ============================================================
# Central estimates using all seeds
# ============================================================
crossing_T_center = {}
excess_peak_T_center = {}
excess_peak_delta_center = {}

print("\n=== Central estimates (averaged curve first) ===")

for angle in angles:
    df_avg = average_across_seeds(angle)
    if df_avg is not None:
        crossing_T_center[angle] = find_crossing_T(df_avg)
    else:
        crossing_T_center[angle] = np.nan

    if angle != 0:
        T_peak, d_peak = compute_excess_peak(angle)
        excess_peak_T_center[angle] = T_peak
        excess_peak_delta_center[angle] = d_peak
        print(f"Angle {angle}°: crossing T = {crossing_T_center[angle]:.3f}, "
              f"peak Δψ6 = {d_peak:.3f} at T = {T_peak:.3f}")

# ============================================================
# Compute actual psi6 at each crossing temperature, with std
# ============================================================
crossing_psi_mean = {}
crossing_psi_std = {}

for angle in angles:
    if np.isnan(crossing_T_center[angle]):
        crossing_psi_mean[angle] = np.nan
        crossing_psi_std[angle] = np.nan
        continue

    Tcross = crossing_T_center[angle]
    vals = []
    for seed in raw[angle]:
        df = raw[angle][seed]
        # Find the row closest to Tcross
        idx = (df["T"] - Tcross).abs().idxmin()
        vals.append(df.loc[idx, "psi_GB"])

    if vals:
        crossing_psi_mean[angle] = np.mean(vals)
        crossing_psi_std[angle] = np.std(vals, ddof=1) if len(vals) > 1 else 0.0
    else:
        crossing_psi_mean[angle] = np.nan
        crossing_psi_std[angle] = np.nan

# ============================================================
# Bootstrap error bars
# ============================================================
crossing_T_boot = {a: [] for a in angles}
excess_peak_T_boot = {a: [] for a in angles if a != 0}
excess_peak_delta_boot = {a: [] for a in angles if a != 0}

common_seeds = list(set.intersection(*(set(raw[a].keys()) for a in angles if raw[a])))

for b in range(n_bootstrap):
    # Resample seeds with replacement
    seed_sample = list(rng.choice(common_seeds, size=len(common_seeds), replace=True))

    for angle in angles:
        df_avg = average_across_seeds(angle, seed_sample)
        if df_avg is not None:
            crossing_T_boot[angle].append(find_crossing_T(df_avg))
        else:
            crossing_T_boot[angle].append(np.nan)

        if angle != 0:
            T_peak, d_peak = compute_excess_peak(angle, seed_sample)
            excess_peak_T_boot[angle].append(T_peak)
            excess_peak_delta_boot[angle].append(d_peak)

# Convert to standard deviation
crossing_T_err = {a: np.nanstd(np.array(crossing_T_boot[a]), ddof=1) for a in angles}
excess_peak_T_err = {a: np.nanstd(np.array(excess_peak_T_boot[a]), ddof=1) for a in excess_peak_T_boot}
excess_peak_delta_err = {a: np.nanstd(np.array(excess_peak_delta_boot[a]), ddof=1) for a in excess_peak_delta_boot}

# ============================================================
# Plot 1: threshold crossing with error bars and connecting lines
# ============================================================
plt.figure(figsize=(8, 5))

plot_angles = [a for a in angles if not np.isnan(crossing_T_center[a])]
plot_T = [crossing_T_center[a] for a in plot_angles]
plot_err = [crossing_T_err[a] for a in plot_angles]

plt.errorbar(plot_angles, plot_T, yerr=plot_err,
             fmt='o-', capsize=5, linewidth=2, color='darkblue',
             markersize=8, elinewidth=1.5)

for a, T, err in zip(plot_angles, plot_T, plot_err):
    psi_mean = crossing_psi_mean[a]
    psi_std = crossing_psi_std[a]

    # Annotation now includes psi6 value ± std and T ± std
    plt.annotate(
        f"ψ6={psi_mean:.3f}±{psi_std:.3f}\nT={T:.3f}±{err:.3f}",
        xy=(a, T),
        xytext=(5, -10),
        textcoords="offset points",
        fontsize=8,
        color="black"
    )

plt.xlabel("Grain-boundary angle / single crystal angle (°)")
plt.ylabel(rf"$T_{{\mathrm{{cross}}}}$  (first $|\psi_6| < {threshold}$)")
plt.title("Disordering temperature vs grain-boundary angle")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{threshold}_threshold_crossing_with_bootstrap_errors.png", dpi=300)
plt.show()

# ============================================================
# Plot 2: excess disorder peak with error bars and connecting lines
# ============================================================
plt.figure(figsize=(8, 5))

plot_angles = [a for a in excess_peak_T_center if not np.isnan(excess_peak_T_center[a])]
plot_T = [excess_peak_T_center[a] for a in plot_angles]
plot_T_err = [excess_peak_T_err[a] for a in plot_angles]
plot_delta = [excess_peak_delta_center[a] for a in plot_angles]
plot_delta_err = [excess_peak_delta_err[a] for a in plot_angles]

plt.errorbar(plot_angles, plot_T, yerr=plot_T_err,
             fmt='s-', capsize=5, linewidth=2, color='firebrick',
             markersize=8, elinewidth=1.5)

for a, T, Terr, delta, derr in zip(plot_angles, plot_T, plot_T_err, plot_delta, plot_delta_err):
    plt.annotate(
        f"Δψ6={delta:.3f}±{derr:.3f}",
        xy=(a, T),
        xytext=(5, 10),
        textcoords="offset points",
        fontsize=8,
        color="black"
    )

plt.xlabel("Grain-boundary angle (°)")
plt.ylabel(r"$T_{\mathrm{peak}}$  (temperature of maximum excess disorder)")
plt.title("Peak grain-boundary premelting vs angle")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(f"{threshold}_excess_disorder_with_bootstrap_errors.png", dpi=300)
plt.show()