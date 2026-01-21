#==========================================================================
# This script  will: 
#  Reads a .MDE file with MD step, temperature, and total energy.
# Skips blank/comment lines.
# Plots temperature and energy on dual Y-axes.
# Saves the plot as a PNG file.
#  The result will be a plot with:
#  X-axis: MD steps
#  Left Y-axis (red): Temperature in Kelvin
#  Right Y-axis (blue): Total energy in Rydbergs
#  Output file: mde_dual_axis_plot.png
# Feel free to drop your quiries :drmohanlv@gmail.com
#=========================================================================
import matplotlib.pyplot as plt
import sys

def read_mde_file(filename):
    steps = []
    temperatures = []
    energies = []

    with open(filename, 'r') as f:
        for line in f:
            if line.strip() == "" or line.startswith("#"):
                continue  # Skip empty or comment lines
            tokens = line.split()
            if len(tokens) >= 3:
                try:
                    steps.append(int(tokens[0]))
                    temperatures.append(float(tokens[1]))
                    energies.append(float(tokens[2]))
                except ValueError:
                    continue  # Skip malformed lines
    return steps, temperatures, energies

# === Handle command-line argument ===
filename = "peg6peo5nai.MDE"
if len(sys.argv) > 1:
    filename = sys.argv[1]

# === Load your data ===
steps, temperatures, energies = read_mde_file(filename)

# === Plot with dual Y-axes ===
fig, ax1 = plt.subplots(figsize=(8, 5))

# Temperature axis (left)
color1 = 'tab:red'
ax1.set_xlabel("MD Step")
ax1.set_ylabel("Temperature (K)", color=color1)
line1, = ax1.plot(steps, temperatures, color=color1, linestyle='-', marker='o', markersize=2, label="Temperature")
ax1.tick_params(axis='y', labelcolor=color1)

# Energy axis (right)
ax2 = ax1.twinx()
color2 = 'tab:blue'
ax2.set_ylabel("Total Energy (Ry)", color=color2)
line2, = ax2.plot(steps, energies, color=color2, linestyle='--', marker='x', markersize=2, label="Total Energy")
ax2.tick_params(axis='y', labelcolor=color2)

# Combined legend
lines = [line1, line2]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc='upper right')

# Title and layout
plt.title("AIMD: Temperature and Total Energy vs MD Step")
fig.tight_layout()
plt.grid(True)
plt.savefig("mde_dual_axis_plot.png")
plt.show()

print("✅ Combined plot saved as: mde_dual_axis_plot.png")

