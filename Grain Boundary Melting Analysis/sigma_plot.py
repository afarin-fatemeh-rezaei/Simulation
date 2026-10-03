#!/usr/bin/env python3
"""
Analyse sigma_gb_csl.csv produced by the 2D LJ CSL grain-boundary code.

Produces three plots:
  1. sigma_GB vs theta          -- the full curve
  2. sigma_GB/theta vs ln(theta) -- Read-Shockley linearisation + fit
  3. sigma_GB vs Sigma_true      -- CSL cusp structure

Also prints the fitted Read-Shockley parameters E0 and A.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

CSV = Path("sigma_gb_csl.csv")
if not CSV.exists():
    raise SystemExit(f"cannot find {CSV}")

# ------------------------------------------------------------------
# Load.  Column order from the C code:
# theta_deg,sigma_gb,Sigma,Sigma_true,m,n,nx,ny,N_atoms,best_sx_a,best_sy_a,pe_final,converged,dmin
# ------------------------------------------------------------------
data = np.genfromtxt(CSV, delimiter=",", names=True, dtype=None, encoding="utf-8")

theta_deg  = np.asarray(data["theta_deg"],  float)
sigma_gb   = np.asarray(data["sigma_gb"],   float)
Sigma      = np.asarray(data["Sigma"],      int)
Sigma_true = np.asarray(data["Sigma_true"], int)
dmin       = np.asarray(data["dmin"],       float)
converged  = np.asarray(data["converged"],  int)

# Filter out any NaNs (should not be any in your current run, but be safe)
good = np.isfinite(sigma_gb) & (converged == 1)
theta_deg  = theta_deg[good]
sigma_gb   = sigma_gb[good]
Sigma      = Sigma[good]
Sigma_true = Sigma_true[good]
dmin       = dmin[good]

# Sort by angle
idx        = np.argsort(theta_deg)
theta_deg  = theta_deg[idx]
sigma_gb   = sigma_gb[idx]
Sigma      = Sigma[idx]
Sigma_true = Sigma_true[idx]
dmin       = dmin[idx]

theta_rad  = np.deg2rad(theta_deg)

# ------------------------------------------------------------------
# 1) sigma_GB vs theta  -- the main curve
# ------------------------------------------------------------------
fig1, ax1 = plt.subplots(figsize=(7.5, 5.0))
ax1.plot(theta_deg, sigma_gb, "o-", color="C0", markersize=5, lw=1.2,
         label=r"$\sigma_{\rm GB}(\theta)$")
ax1.set_xlabel(r"misorientation angle  $\theta$  (deg)", fontsize=12)
ax1.set_ylabel(r"$\sigma_{\rm GB}$  (LJ units)", fontsize=12)
ax1.set_title(r"Grain-boundary energy vs misorientation angle")
ax1.grid(alpha=0.3)
ax1.legend()
fig1.tight_layout()
fig1.savefig("sigma_vs_theta.png", dpi=150)

# ------------------------------------------------------------------
# 2) Read-Shockley linearisation: sigma/theta vs ln(theta)
#    sigma/theta = E0 * (A - ln theta)
#    => y = c + m * x  with  x = ln theta, y = sigma/theta
#       m = -E0, c = E0*A
# ------------------------------------------------------------------
fig2, ax2 = plt.subplots(figsize=(7.5, 5.0))

# Only use small-angle region for the fit (Read-Shockley is a
# small-angle expansion).  Adjust THETA_FIT_MAX if you want.
THETA_FIT_MAX = 15.0
mask_fit = theta_deg <= THETA_FIT_MAX

x = np.log(theta_rad[mask_fit])
y = sigma_gb[mask_fit] / theta_rad[mask_fit]

ax2.plot(np.log(theta_rad), sigma_gb / theta_rad, "o", color="C1",
         label=r"all data")
ax2.plot(x, y, "s", color="C3", label=fr"fit range  $\theta \leq {THETA_FIT_MAX:.0f}^\circ$")

# Linear least-squares fit  y = m*x + c
m, c = np.polyfit(x, y, 1)
E0   = -m
A_rs = c / E0

xx = np.linspace(x.min(), x.max(), 100)
ax2.plot(xx, m*xx + c, "--", color="k",
         label=fr"fit: $E_0={E0:.4f}$, $A={A_rs:.3f}$")

ax2.set_xlabel(r"$\ln\theta$   ($\theta$ in radians)", fontsize=12)
ax2.set_ylabel(r"$\sigma_{\rm GB}/\theta$", fontsize=12)
ax2.set_title(r"Read–Shockley linearisation")
ax2.grid(alpha=0.3)
ax2.legend()
fig2.tight_layout()
fig2.savefig("readshockley_fit.png", dpi=150)

# Reconstruct the RS curve on the whole range for the first plot
theta_fine = np.linspace(theta_rad.min(), theta_rad.max(), 400)
sigma_rs   = E0 * theta_fine * (A_rs - np.log(theta_fine))
ax1.plot(np.rad2deg(theta_fine), sigma_rs, "--", color="C3",
         label=fr"Read–Shockley fit ($E_0={E0:.3f}$, $A={A_rs:.2f}$)")
ax1.legend()
fig1.savefig("sigma_vs_theta.png", dpi=150)

# ------------------------------------------------------------------
# 3) sigma_GB vs Sigma_true  -- CSL cusp structure
# ------------------------------------------------------------------
fig3, ax3 = plt.subplots(figsize=(7.5, 5.0))

order = np.argsort(Sigma_true)
ax3.plot(Sigma_true[order], sigma_gb[order], "o", color="C2", markersize=6)

# Annotate low-Sigma points
for s, t, sig in zip(Sigma_true, theta_deg, sigma_gb):
    if s <= 30:
        ax3.annotate(fr"$\Sigma={s}$" + "\n" + fr"${t:.1f}^\circ$",
                     xy=(s, sig), xytext=(5, 5), textcoords="offset points",
                     fontsize=8)

ax3.set_xscale("log")
ax3.set_xlabel(r"$\Sigma_{\rm true}$  (log scale)", fontsize=12)
ax3.set_ylabel(r"$\sigma_{\rm GB}$  (LJ units)", fontsize=12)
ax3.set_title(r"Grain-boundary energy vs CSL index")
ax3.grid(alpha=0.3, which="both")
fig3.tight_layout()
fig3.savefig("sigma_vs_sigma.png", dpi=150)

# ------------------------------------------------------------------
# 3b) Cusp structure:  residual sigma - smooth(theta)  vs  Sigma_true
# ------------------------------------------------------------------
# Build a smooth background by fitting a low-order polynomial in theta
# (or use a running mean / LOESS if you prefer).  A cubic in log-theta
# works well for this kind of curve.
from numpy.polynomial import polynomial as P

coef = np.polyfit(np.log(theta_rad), sigma_gb, 3)   # cubic in ln(theta)
sigma_smooth = np.polyval(coef, np.log(theta_rad))
residual = sigma_gb - sigma_smooth

fig3b, ax3b = plt.subplots(figsize=(7.5, 5.0))
ax3b.axhline(0.0, color="k", lw=0.8, ls="--")
ax3b.plot(Sigma_true, residual, "o", color="C2", markersize=6)

# Label the low-Sigma points so the cusps are visible
for s, t, r in zip(Sigma_true, theta_deg, residual):
    if s <= 40:
        ax3b.annotate(fr"$\Sigma={s}$" + "\n" + fr"${t:.1f}^\circ$",
                      xy=(s, r), xytext=(6, 6), textcoords="offset points",
                      fontsize=8)

ax3b.set_xscale("log")
ax3b.set_xlabel(r"$\Sigma_{\rm true}$  (log scale)", fontsize=12)
ax3b.set_ylabel(r"$\sigma_{\rm GB} - \bar\sigma(\theta)$", fontsize=12)
ax3b.set_title(r"CSL cusp structure (residual vs smooth background)")
ax3b.grid(alpha=0.3, which="both")
fig3b.tight_layout()
fig3b.savefig("sigma_cusps.png", dpi=150)

# ------------------------------------------------------------------
# Report
# ------------------------------------------------------------------
print("=" * 60)
print(f"Read–Shockley fit on  theta <= {THETA_FIT_MAX:.1f} deg")
print(f"  E0 = {E0:.5f}")
print(f"  A  = {A_rs:.5f}")
print(f"  theta* = exp(A-1) = {np.rad2deg(np.exp(A_rs-1)):.2f} deg"
      "  (maximum of the RS curve)")
print("=" * 60)
print(f"Number of points plotted: {len(theta_deg)}")
print(f"sigma_min = {sigma_gb.min():.4f} at theta = {theta_deg[np.argmin(sigma_gb)]:.3f} deg")
print(f"sigma_max = {sigma_gb.max():.4f} at theta = {theta_deg[np.argmax(sigma_gb)]:.3f} deg")
print(f"mean(sigma | theta>25 deg) = {sigma_gb[theta_deg>25].mean():.4f}"
      "   <-- plateau value")
print("=" * 60)

plt.show()