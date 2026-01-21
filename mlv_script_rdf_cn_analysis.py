#============================================================================
#This script will:
# Convert .ANI to .XYZ
# Apply cell info from STRUCT_OUT (if missing)
# Compute RDFs for defined element pairs
# Save RDF data and a combined RDF plot
#===========================================================================
import os
import numpy as np
from ase import Atoms
from ase.io import read, write
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt
from tqdm import tqdm

# === USER CONFIGURATION ===
ani_file = "peg6peo5nai.ANI"
xyz_outfile = "peg6peo5nai_from_ani.xyz"
struct_out = "peg6peo5nai.STRUCT_OUT"
rdf_dir = "rdf_results"
coord_output = "coordination_numbers.txt"
element_pairs = [("O", "Na"), ("O", "H"), ("O", "O")]
nbins = 200
coordination_cutoff = 3.5  # Å

# === ATOMIC SYMBOLS ===
symbols_reference = ['C'] * 20 + ['O'] * 11 + ['H'] * 43 + ['Na'] * 2 + ['I'] * 2
natoms = len(symbols_reference)

# === 1. CONVERT .ANI TO .XYZ ===
print(f"\n📦 Converting {ani_file} to XYZ...")
frames = []
with open(ani_file, 'r') as f:
    lines = f.readlines()

i = 0
while i < len(lines):
    if not lines[i].strip():
        i += 1
        continue

    atom_count = int(lines[i])
    if atom_count != natoms:
        raise ValueError(f"❌ Atom count mismatch at line {i}: {atom_count} != {natoms}")
    i += 2  # skip comment

    positions, symbols = [], []
    for _ in range(natoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        positions.append(list(map(float, parts[1:4])))
        i += 1

    atoms = Atoms(symbols=symbols, positions=positions)
    frames.append(atoms)

write(xyz_outfile, frames)
print(f"✅ Wrote {len(frames)} frames to {xyz_outfile}")

# === 2. LOAD TRAJECTORY ===
traj = read(xyz_outfile, index=":", format="xyz")
natoms = len(traj[0])
print(f"✅ Loaded {len(traj)} frames, {natoms} atoms per frame")

# === 3. Load CELL IF MISSING ===
if traj[0].get_cell().volume == 0.0:
    print(f"📂 Adding cell info from {struct_out}")
    with open(struct_out, 'r') as f:
        lines = f.readlines()
    cell = np.array([[float(x) for x in line.split()] for line in lines[:3]])
    for atoms in traj:
        atoms.set_cell(cell)
        atoms.set_pbc(True)
    print(f"✅ Cell vectors set:\n{cell}")

# === 4. Set rmax from cell length ===
cell_lengths = traj[0].get_cell_lengths_and_angles()[:3]
rmax = 0.95 * min(cell_lengths) / 2
print(f"📏 Chosen rmax = {rmax:.2f} Å")

# === 5. Prepare outputs ===
os.makedirs(rdf_dir, exist_ok=True)
coord_results = []

# === 6. RDF + Coordination Numbers ===
all_rdfs = []
for el1, el2 in element_pairs:
    print(f"\n🔍 RDF + Coordination for {el1}–{el2}")
    try:
        analysis = Analysis(traj)
        rdf = analysis.get_rdf(rmax=rmax, nbins=nbins, elements=(el1, el2))[0]
        radii = np.linspace(0, rmax, nbins)
        all_rdfs.append((el1, el2, radii, rdf))

        # Save RDF
        outfile = os.path.join(rdf_dir, f"rdf_{el1}_{el2}.txt")
        np.savetxt(outfile, np.column_stack([radii, rdf]), header="r [Å]    g(r)")
        print(f"✅ Saved RDF: {outfile}")

        # Compute coordination number
        all_symbols = traj[0].get_chemical_symbols()
        count_el2 = all_symbols.count(el2)
        volume = traj[0].get_volume()
        rho_el2 = count_el2 / volume  # atoms per Å³

        dr = radii[1] - radii[0]
        idx = np.where(radii <= coordination_cutoff)[0][-1]
        r_cut = radii[:idx+1]
        g_cut = rdf[:idx+1]

        N_coord = 4 * np.pi * rho_el2 * np.sum(g_cut * r_cut**2 * dr)
        coord_results.append((el1, el2, N_coord))
        print(f"📊 Coordination number {el1}-{el2} (r<{coordination_cutoff} Å): {N_coord:.2f}")

    except Exception as e:
        print(f"❌ Error with pair {el1}-{el2}: {e}")

# === 7. Save coordination numbers ===
with open(coord_output, 'w') as f:
    f.write("# Coordination numbers (r < {:.2f} Å)\n".format(coordination_cutoff))
    f.write("# Pair      CN\n")
    for el1, el2, N in coord_results:
        f.write(f"{el1}-{el2:<6}  {N:.4f}\n")

print(f"\n📝 Coordination numbers saved to {coord_output}")

# === 8. Plot all RDFs ===
if all_rdfs:
    plt.figure(figsize=(8, 5))
    for el1, el2, radii, rdf in all_rdfs:
        plt.plot(radii, rdf, label=f"{el1}–{el2}")
    plt.xlabel("r (Å)")
    plt.ylabel("g(r)")
    plt.title("Radial Distribution Function (RDF)")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("rdf_all_pairs.png")
    print("📈 RDF plot saved: rdf_all_pairs.png")
else:
    print("⚠️ No RDFs computed.")

