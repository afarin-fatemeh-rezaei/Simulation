#!/usr/bin/env python3
"""
Generate a perfect triangular lattice:
    a = 1.2
    x ∈ [-100a, 100a] = [-120, 120]
    y ∈ [-36.373, 36.373]   (commensurate with sqrt(3)*a)

Output: atoms_0deg.json

Also prints the 6-fold coordination check.
"""

import numpy as np
import json

# Parameters
a = 1.2
xmin = -100.0 * a
xmax =  100.0 * a
ymin = -36.373
ymax =  36.373

# Primitive vectors of the unrotated triangular lattice
e1 = np.array([a, 0.0])
e2 = np.array([a * 0.5, a * np.sqrt(3.0) / 2.0])

# Number of lattice points needed
N = int(np.ceil(max(xmax - xmin, ymax - ymin) / a)) + 5

points = []
for i in range(-N, N + 1):
    for j in range(-N, N + 1):
        p = i * e1 + j * e2
        x, y = p[0], p[1]
        if xmin <= x <= xmax and ymin <= y <= ymax:
            points.append((x, y))

# Save to JSON
output_file = "atoms_0deg.json"
with open(output_file, "w") as f:
    f.write("[\n")
    for idx, (x, y) in enumerate(points):
        comma = "," if idx < len(points) - 1 else ""
        f.write(f'  {{"x": {x:.12f}, "y": {y:.12f}}}{comma}\n')
    f.write("]\n")

print(f"Generated {len(points)} atoms")
print(f"Saved to {output_file}")
print(f"x range: [{min(p[0] for p in points):.6f}, {max(p[0] for p in points):.6f}]")
print(f"y range: [{min(p[1] for p in points):.6f}, {max(p[1] for p in points):.6f}]")

# 6-fold coordination check
x_arr = np.array([p[0] for p in points])
y_arr = np.array([p[1] for p in points])

Ly = ymax - ymin   # 72.746
cut = 1.8
neighbor_counts = []

for i in range(len(x_arr)):
    dx = x_arr - x_arr[i]
    dy = y_arr - y_arr[i]
    # periodic boundary in y only
    dy -= Ly * np.round(dy / Ly)
    r2 = dx*dx + dy*dy
    n = np.sum((r2 < cut**2) & (r2 > 1e-12))
    neighbor_counts.append(n)

counts = np.bincount(neighbor_counts)
print("\n6-fold coordination check (cutoff = 1.8):")
for k in range(len(counts)):
    if counts[k] > 0:
        print(f"  {k} neighbours: {counts[k]} atoms")