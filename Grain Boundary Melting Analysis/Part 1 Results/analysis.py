import pandas as pd
import matplotlib.pyplot as plt

b  = "30_"
a = "6_psi_GB_vs_T"
# Read the two CSV files
bi = pd.read_csv(f"{b}{a}.csv")
single = pd.read_csv(f"0_{a}.csv")

# Plot both lines
plt.figure(figsize=(8, 5))

plt.plot(bi["T"], bi["psi_GB"],
         label="Bicrystal GB", color="red", linewidth=2, marker="o")

plt.plot(single["T"], single["psi_GB"],
         label="Single crystal (same probe)", color="blue", linewidth=2, marker="s")

plt.xlabel("Temperature")
plt.ylabel(r"$|\psi_6|$")
plt.title(r"Grain-boundary region $|\psi_6|$ vs temperature")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()

# Save figure
plt.savefig(f"{b}{a}.png", dpi=300)

# Show plot
#plt.show()