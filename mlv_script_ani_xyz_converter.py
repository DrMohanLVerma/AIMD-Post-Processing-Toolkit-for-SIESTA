from ase import Atoms
from ase.io import write
import os

# === USER INPUT ===
ani_file = "peg6peo5nai.ANI"
xyz_outfile = "peg6peo5nai_from_ani.xyz"

# === ATOMIC SYMBOLS ===
# Your atom count: 20 C, 11 O, 43 H, 2 Na, 2 I => Total = 78
symbols_reference = ['C'] * 20 + ['O'] * 11 + ['H'] * 43 + ['Na'] * 2 + ['I'] * 2
natoms = len(symbols_reference)

# === READ ANI FILE ===
print(f"📦 Reading {ani_file}")
frames = []
with open(ani_file, 'r') as f:
    lines = f.readlines()

i = 0
frame_count = 0
while i < len(lines):
    line = lines[i].strip()
    if not line:
        i += 1
        continue

    try:
        atom_count = int(line)
    except ValueError:
        raise ValueError(f"❌ Expected number of atoms at line {i}, got: {line}")

    if atom_count != natoms:
        raise ValueError(f"❌ Atom count mismatch at line {i}: found {atom_count}, expected {natoms}")

    i += 1  # skip comment line
    comment = lines[i].strip()
    i += 1

    positions = []
    symbols = []

    for j in range(natoms):
        parts = lines[i].strip().split()
        if len(parts) != 4:
            raise ValueError(f"❌ Malformed atom line at line {i}: {lines[i].strip()}")

        symbol = parts[0]
        x, y, z = map(float, parts[1:4])
        positions.append([x, y, z])
        symbols.append(symbol)
        i += 1

    atoms = Atoms(symbols=symbols, positions=positions)
    frames.append(atoms)
    frame_count += 1

print(f"✅ Extracted {frame_count} frames with {natoms} atoms each.")
print(f"💾 Writing to XYZ file: {xyz_outfile}")
write(xyz_outfile, frames)

print("✅ Conversion complete.")

