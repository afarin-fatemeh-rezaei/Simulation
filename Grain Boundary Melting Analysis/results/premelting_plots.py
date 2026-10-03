#!/usr/bin/env python3
"""
plot_premelting.py

Plots the premelting analysis for the 6-angle series.

Angles are always ordered LARGEST -> SMALLEST misorientation angle,
both in the fig1 subplot layout (left-to-right, top-to-bottom) and
in the fig2 / fig3 legend.

Reads sweep files of the form:
  sweep_perfect_Sigma<X>_theta<Y>deg_vs_bicrystal_Sigma<X>_theta<Y>deg_seed<Z>.dat

Produces:
  fig1_melting_curves.png    - per-angle psi6 vs T (perfect vs bicrystal)
  fig2_dPE_vs_T.png          - Delta PE vs T for all angles
  fig3_dpsi6_vs_T.png        - d_psi6_gb vs T for all angles
  fig4_angle_dependence.png  - headline: strength vs misorientation angle
  premelting_summary.txt     - tabulated numbers for the paper

Usage:
  python plot_premelting.py                        # auto-find sweep_*.dat
  python plot_premelting.py file1.dat file2.dat ...
"""

import sys
import re
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm


# ------------------------------------------------------------------
# Column layout (v6.5 sweep file)
# ------------------------------------------------------------------
C = {
    'T':           0,
    'psi6_avg_P':  1,  'psi6_gb_P':  2,  'psi6_bulk_P': 3,  'PE_P': 4,  'P_P': 5,
    'psi6_avg_B':  6,  'psi6_gb_B':  7,  'psi6_bulk_B': 8,  'PE_B': 9,  'P_B': 10,
    'gb_bin_B':   11,  'bulk_bin_B': 12,
    'd_psi6_avg': 13,  'd_psi6_gb': 14, 'd_psi6_bulk': 15, 'd_PE': 16,
}


# ------------------------------------------------------------------
# Filename parsing
# ------------------------------------------------------------------
FNAME_RE = re.compile(
    r"sweep_perfect_Sigma(?P<sigma>\d+)_theta(?P<theta>[\d.]+)deg"
    r"_vs_bicrystal_Sigma\d+_theta[\d.]+deg"
    r"_seed(?P<seed>\d+)\.dat"
)

def parse_fname(path):
    base = path.split('/')[-1].split('\\')[-1]
    m = FNAME_RE.search(base)
    if not m:
        return None
    return {
        'sigma': int(m.group('sigma')),
        'theta': float(m.group('theta')),
        'seed':  int(m.group('seed')),
    }


# ------------------------------------------------------------------
# Sweep file reader
# ------------------------------------------------------------------
def read_sweep(path):
    rows = []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.split()
            if len(parts) < 17:
                continue
            try:
                rows.append([float(x) for x in parts])
            except ValueError:
                continue
    if not rows:
        raise RuntimeError(f"No data in {path}")
    return np.array(rows)


# ------------------------------------------------------------------
# Collect all runs.  Angles sorted by DESCENDING theta.
# ------------------------------------------------------------------
def collect(paths):
    runs = []
    for p in paths:
        meta = parse_fname(p)
        if meta is None:
            print(f"  skipped (name pattern): {p}")
            continue
        data = read_sweep(p)
        meta['path'] = p
        meta['data'] = data
        runs.append(meta)
    if not runs:
        raise RuntimeError("No valid sweep files found")

    # Angles sorted by DESCENDING theta
    angles = sorted(set((r['sigma'], r['theta']) for r in runs),
                    key=lambda x: -x[1])

    print(f"\nFound {len(runs)} runs across {len(angles)} angles "
          f"(ordered largest -> smallest theta):")
    for sig, th in angles:
        seeds = sorted(r['seed'] for r in runs
                       if (r['sigma'], r['theta']) == (sig, th))
        print(f"  Sigma{sig:>4d}  theta = {th:6.4f} deg   seeds = {seeds}")
    return runs, angles


# ------------------------------------------------------------------
# Aggregate per angle (mean/std across seeds)
# ------------------------------------------------------------------
def aggregate(runs, sigma, theta):
    selected = [r for r in runs
                if r['sigma'] == sigma and r['theta'] == theta]
    T = selected[0]['data'][:, C['T']]
    for r in selected[1:]:
        assert np.allclose(r['data'][:, C['T']], T), \
            f"T grid mismatch for Sigma{sigma}"
    stack = np.stack([r['data'] for r in selected], axis=0)
    mean = stack.mean(axis=0)
    std  = stack.std(axis=0)
    return {
        'T': T,
        'mean': mean,
        'std': std,
        'n_seeds': len(selected),
        'seeds': sorted(r['seed'] for r in selected),
        'runs': selected,
    }


def find_T_index(T, target):
    return int(np.argmin(np.abs(T - target)))


# ------------------------------------------------------------------
# Figure 1: per-angle melting curves.
# Subplots are filled in the order of the `angles` list, which is
# already DESCENDING theta.  So position (1,1) is the largest angle,
# and the last position is the smallest angle.
# ------------------------------------------------------------------
def fig1_melting_curves(runs, angles, outname):
    n = len(angles)
    ncol = 3
    nrow = (n + ncol - 1) // ncol

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2*ncol, 3.6*nrow))
    axes = axes.ravel()

    for ax, (sig, th) in zip(axes, angles):
        agg = aggregate(runs, sig, th)
        T = agg['T']
        m = agg['mean']
        s = agg['std']

        ax.errorbar(T, m[:, C['psi6_avg_P']], yerr=s[:, C['psi6_avg_P']],
                    fmt='s-', color='#d62728', lw=2, capsize=3,
                    label='perfect, box avg')
        ax.errorbar(T, m[:, C['psi6_avg_B']], yerr=s[:, C['psi6_avg_B']],
                    fmt='o-', color='#1f77b4', lw=2, capsize=3,
                    label='bicrystal, box avg')
        ax.axvspan(0.38, 0.42, color='gray', alpha=0.10)
        ax.set_xlabel('T')
        ax.set_ylabel(r'$\psi_6$')
        ax.set_title(fr'$\Sigma${sig}  $\theta$ = {th:.4f}$\degree$  '
                     f'({agg["n_seeds"]} seed' +
                     ('s' if agg['n_seeds'] > 1 else '') + ')')
        ax.set_ylim(0.40, 1.00)
        ax.grid(alpha=0.3)
        if ax is axes[0]:
            ax.legend(fontsize=8, loc='lower left')

    for ax in axes[n:]:
        ax.axis('off')

    fig.suptitle('Melting-point depression: box-averaged '
                 r'$\psi_6$ vs T  (largest angle -> smallest, '
                 'reading left-to-right, top-to-bottom)',
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(outname, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved {outname}")


# ------------------------------------------------------------------
# Figure 2: dPE vs T for all angles.
# The loop over `angles` is in DESCENDING theta, so the legend is
# ordered largest -> smallest.
# ------------------------------------------------------------------
def fig2_dPE_vs_T(runs, angles, outname):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    cmap = cm.viridis

    thetas = [th for _, th in angles]
    tmin, tmax = min(thetas), max(thetas)

    # iterate in descending theta (already the order of `angles`)
    for sig, th in angles:
        agg = aggregate(runs, sig, th)
        T = agg['T']
        m = agg['mean']
        s = agg['std']
        color = cmap((th - tmin) / (tmax - tmin + 1e-9))
        label = fr'$\theta$ = {th:.2f}$\degree$  ($\Sigma${sig})'
        ax.errorbar(T, m[:, C['d_PE']], yerr=s[:, C['d_PE']],
                    fmt='o-', color=color, lw=2, ms=5, capsize=3,
                    label=label)

    ax.set_xlabel('T (reduced)')
    ax.set_ylabel(r'$\Delta$PE per atom  (bicrystal $-$ perfect)')
    ax.set_title('Grain-boundary excess energy vs temperature')
    ax.legend(fontsize=8, ncol=2, loc='upper left')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outname, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved {outname}")


# ------------------------------------------------------------------
# Figure 3: d_psi6_gb vs T for all angles.
# Same descending order as fig2.
# ------------------------------------------------------------------
def fig3_dpsi6_vs_T(runs, angles, outname):
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    cmap = cm.viridis

    thetas = [th for _, th in angles]
    tmin, tmax = min(thetas), max(thetas)

    for sig, th in angles:
        agg = aggregate(runs, sig, th)
        T = agg['T']
        m = agg['mean']
        s = agg['std']
        color = cmap((th - tmin) / (tmax - tmin + 1e-9))
        label = fr'$\theta$ = {th:.2f}$\degree$  ($\Sigma${sig})'
        ax.errorbar(T, m[:, C['d_psi6_gb']], yerr=s[:, C['d_psi6_gb']],
                    fmt='o-', color=color, lw=2, ms=5, capsize=3,
                    label=label)

    ax.axhline(0.0, color='k', ls='--', lw=1, alpha=0.6)
    ax.set_xlabel('T (reduced)')
    ax.set_ylabel(r'$\psi_6^{\rm GB} - \psi_6^{\rm perfect}$ at same bin')
    ax.set_title('GB structural disorder vs temperature')
    ax.legend(fontsize=8, ncol=2, loc='lower left')
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(outname, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved {outname}")


# ------------------------------------------------------------------
# Figure 4: angle dependence of premelting strength.
# The printed table is ordered largest -> smallest theta, and the
# markers on the plot are the same angles. No legend, so ordering
# is only for the printed summary.
# ------------------------------------------------------------------
def fig4_angle_dependence(runs, angles, outname, T_target=0.28):
    thetas = []
    dPE_mean = []; dPE_err = []
    dpsi_mean = []; dpsi_err = []
    nseeds = []

    print(f"\nExtracting strengths at T = {T_target}  "
          f"(largest angle first)")
    print(f"{'Sigma':>6}  {'theta':>8}  {'seeds':>5}  "
          f"{'dPE':>10}  {'+/-':>8}  {'d_psi6_gb':>12}  {'+/-':>8}")
    for sig, th in angles:
        agg = aggregate(runs, sig, th)
        i = find_T_index(agg['T'], T_target)
        thetas.append(th)
        dPE_mean.append(agg['mean'][i, C['d_PE']])
        dPE_err.append(agg['std'][i, C['d_PE']])
        dpsi_mean.append(agg['mean'][i, C['d_psi6_gb']])
        dpsi_err.append(agg['std'][i, C['d_psi6_gb']])
        nseeds.append(agg['n_seeds'])
        print(f"{sig:>6}  {th:>8.4f}  {agg['n_seeds']:>5}  "
              f"{agg['mean'][i, C['d_PE']]:>10.5f}  "
              f"{agg['std'][i, C['d_PE']]:>8.5f}  "
              f"{agg['mean'][i, C['d_psi6_gb']]:>12.5f}  "
              f"{agg['std'][i, C['d_psi6_gb']]:>8.5f}")

    thetas    = np.array(thetas)
    dPE_mean  = np.array(dPE_mean)
    dPE_err   = np.array(dPE_err)
    dpsi_mean = np.array(dpsi_mean)
    dpsi_err  = np.array(dpsi_err)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # --- panel A: Delta PE ---
    ax = axes[0]
    ax.errorbar(thetas, dPE_mean, yerr=dPE_err,
                fmt='o-', color='#1f77b4', lw=2.5, ms=10, capsize=6,
                markerfacecolor='white', markeredgewidth=2)
    for th, y, n in zip(thetas, dPE_mean, nseeds):
        ax.annotate(f'{n}', (th, y), xytext=(0, 8),
                    textcoords='offset points',
                    ha='center', fontsize=8, color='gray')
    ax.set_xlabel(r'misorientation angle $\theta$ ($\degree$)')
    ax.set_ylabel(r'$\Delta$PE per atom at T = ' + f'{T_target}')
    ax.set_title('(a) GB excess energy vs angle')
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 35)

    # --- panel B: d_psi6_gb ---
    ax = axes[1]
    ax.errorbar(thetas, dpsi_mean, yerr=dpsi_err,
                fmt='s-', color='#d62728', lw=2.5, ms=10, capsize=6,
                markerfacecolor='white', markeredgewidth=2)
    ax.axhline(0.0, color='k', ls='--', lw=1, alpha=0.6)
    for th, y, n in zip(thetas, dpsi_mean, nseeds):
        ax.annotate(f'{n}', (th, y), xytext=(0, 8),
                    textcoords='offset points',
                    ha='center', fontsize=8, color='gray')
    ax.set_xlabel(r'misorientation angle $\theta$ ($\degree$)')
    ax.set_ylabel(r'$\psi_6^{\rm GB} - \psi_6^{\rm perfect}$ at T = '
                  + f'{T_target}')
    ax.set_title('(b) GB structural disorder vs angle')
    ax.grid(alpha=0.3)
    ax.set_xlim(0, 35)

    fig.suptitle('Premelting strength vs misorientation angle',
                 fontsize=14, y=1.02)
    fig.tight_layout()
    fig.savefig(outname, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved {outname}")


# ------------------------------------------------------------------
# Figure 5: per-angle PE curves and their difference
# 6 subplots in the same layout as fig1.
# Each subplot shows:
#   - PE_P (perfect, box avg) on left axis
#   - PE_B (bicrystal, box avg) on left axis
#   - dPE = PE_B - PE_P on right axis (the excess energy)
# ------------------------------------------------------------------
def fig5_dPE_per_angle(runs, angles, outname):
    n = len(angles)
    ncol = 3
    nrow = (n + ncol - 1) // ncol

    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4*ncol, 3.8*nrow))
    axes = axes.ravel()

    for ax, (sig, th) in zip(axes, angles):
        agg = aggregate(runs, sig, th)
        T = agg['T']
        m = agg['mean']
        s = agg['std']

        # Absolute energies on the left axis
        ax.errorbar(T, m[:, C['PE_P']], yerr=s[:, C['PE_P']],
                    fmt='s-', color='#d62728', lw=1.6, ms=4, capsize=2,
                    alpha=0.7, label='perfect')
        ax.errorbar(T, m[:, C['PE_B']], yerr=s[:, C['PE_B']],
                    fmt='o-', color='#1f77b4', lw=1.6, ms=4, capsize=2,
                    alpha=0.7, label='bicrystal')
        ax.set_xlabel('T')
        ax.set_ylabel(r'PE per atom')
        ax.grid(alpha=0.3)

        # Excess energy on the right axis
        ax2 = ax.twinx()
        ax2.errorbar(T, m[:, C['d_PE']], yerr=s[:, C['d_PE']],
                     fmt='D-', color='#2ca02c', lw=2.2, ms=5, capsize=3,
                     label=r'$\Delta$PE')
        ax2.set_ylabel(r'$\Delta$PE per atom', color='#2ca02c')
        ax2.tick_params(axis='y', labelcolor='#2ca02c')
        ax2.set_ylim(bottom=0)   # ΔPE is always positive

        ax.set_title(fr'$\Sigma${sig}  $\theta$ = {th:.4f}$\degree$  '
                     f'({agg["n_seeds"]} seed' +
                     ('s' if agg['n_seeds'] > 1 else '') + ')')

        # Combined legend for both axes
        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        if ax is axes[0]:
            ax.legend(h1 + h2, l1 + l2, fontsize=7,
                      loc='lower left', framealpha=0.9)

    for ax in axes[n:]:
        ax.axis('off')

    fig.suptitle('Grain-boundary excess energy: absolute PE curves and '
                 r'$\Delta$PE per angle  (largest angle -> smallest)',
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(outname, bbox_inches='tight', dpi=150)
    plt.close(fig)
    print(f"Saved {outname}")


# ------------------------------------------------------------------
# Text summary (angles ordered largest -> smallest theta)
# ------------------------------------------------------------------
def write_summary(runs, angles, outname, T_refs=(0.28, 0.32, 0.36)):
    lines = []
    lines.append("Premelting summary for 2D LJ bicrystal series")
    lines.append("Angles ordered largest -> smallest misorientation angle")
    lines.append("=" * 70)
    lines.append("")
    for sig, th in angles:
        agg = aggregate(runs, sig, th)
        lines.append(f"Sigma{sig}  theta = {th:.4f} deg  "
                     f"seeds = {agg['seeds']}")
        lines.append("-" * 60)
        lines.append(f"{'T':>6}  {'dPE':>10}  {'+/-':>8}  "
                     f"{'d_psi6_gb':>12}  {'+/-':>8}  "
                     f"{'d_psi6_bulk':>13}  {'+/-':>8}")
        for T_ref in T_refs:
            i = find_T_index(agg['T'], T_ref)
            lines.append(
                f"{agg['T'][i]:>6.3f}  "
                f"{agg['mean'][i, C['d_PE']]:>10.5f}  "
                f"{agg['std'][i, C['d_PE']]:>8.5f}  "
                f"{agg['mean'][i, C['d_psi6_gb']]:>12.5f}  "
                f"{agg['std'][i, C['d_psi6_gb']]:>8.5f}  "
                f"{agg['mean'][i, C['d_psi6_bulk']]:>13.5f}  "
                f"{agg['std'][i, C['d_psi6_bulk']]:>8.5f}"
            )
        lines.append("")

    text = "\n".join(lines)
    with open(outname, "w") as f:
        f.write(text)
    print(f"Saved {outname}")
    print()
    print(text)


# ------------------------------------------------------------------
# Main
# ------------------------------------------------------------------
def main():
    if len(sys.argv) > 1:
        paths = sys.argv[1:]
    else:
        paths = sorted(glob.glob("sweep_*.dat"))
        if not paths:
            print("No sweep_*.dat files found.")
            sys.exit(1)
        print(f"Auto-found {len(paths)} sweep files")

    runs, angles = collect(paths)

    fig1_melting_curves(runs, angles, "fig1_melting_curves.png")
    fig2_dPE_vs_T(runs, angles, "fig2_dPE_vs_T.png")
    fig3_dpsi6_vs_T(runs, angles, "fig3_dpsi6_vs_T.png")
    fig4_angle_dependence(runs, angles, "fig4_angle_dependence.png",
                          T_target=0.28)
    fig5_dPE_per_angle(runs, angles, "fig5_dPE_per_angle.png")     # NEW
    write_summary(runs, angles, "premelting_summary.txt")

    print("\nAll figures generated.")
    print("Files written:")
    print("  fig1_melting_curves.png    per-angle psi6 melting curves")
    print("  fig2_dPE_vs_T.png          Delta PE vs T (all on one axis)")
    print("  fig3_dpsi6_vs_T.png        GB structural disorder vs T")
    print("  fig4_angle_dependence.png  headline angle-dependence figure")
    print("  fig5_dPE_per_angle.png     per-angle PE curves + Delta PE")  # NEW
    print("  premelting_summary.txt     numeric table for the paper")


if __name__ == "__main__":
    main()