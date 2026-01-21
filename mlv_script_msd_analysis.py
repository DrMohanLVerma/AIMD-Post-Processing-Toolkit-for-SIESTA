#================================================
# This script will : 
# Loads a trajectory from peg6peo5nai_from_ani.xyz
# Computes MSD:
# Averages over all atoms
# Plots and saves the MSD curve: MSD_vs_Time.png
#================================================
import numpy as np
from ase.io import read
import matplotlib.pyplot as plt

xyz_file = "peg6peo5nai_from_ani.xyz"
traj = read(xyz_file, index=":")
nframes = len(traj)
natoms = len(traj[0])

ref_pos = traj[0].get_positions()
msd = []

for i, atoms in enumerate(traj):
    disp = atoms.get_positions() - ref_pos
    squared = np.sum(disp**2, axis=1)
    msd.append(np.mean(squared))

plt.figure(figsize=(8,5))
plt.plot(range(nframes), msd, color="blue")
plt.xlabel("Time Step")
plt.ylabel("MSD (Å²)")
plt.title("Mean Squared Displacement (MSD)")
plt.grid(True)
plt.tight_layout()
plt.savefig("MSD_vs_Time.png")
plt.show()
print("✅ MSD plot saved as 'MSD_vs_Time.png'")
