#===================================================================
# This script will : 
# Reads SIESTA .ANI MD trajectories
# Computes mean square displacement (MSD) for specified ions (default: Na)
# Fits diffusion coefficient 𝐷
# Applies Yeh–Hummer correction
# Computes ionic conductivity using Nernst–Einstein relation
# Outputs a CSV with raw & corrected  𝐷 and 𝜎
#=============================================================================
import numpy as np
import matplotlib.pyplot as plt
import os, glob, csv

kB = 1.380649e-23  # J/K
e = 1.602176634e-19  # C
NA = 6.02214076e23  # 1/mol

def read_ani(filename):
    positions, species = [], []
    with open(filename, 'r') as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        try:
            natoms = int(lines[i].strip())
            i += 2  # skip comment
            frame, frame_species = [], []
            for _ in range(natoms):
                parts = lines[i].split()
                frame_species.append(parts[0])
                frame.append([float(x) for x in parts[1:4]])
                i += 1
            positions.append(frame)
            if not species:
                species = frame_species
        except:
            break
    return np.array(positions), species

def get_atom_indices_by_element(species_list, target_element):
    return [i for i, sym in enumerate(species_list) if sym.lower() == target_element.lower()]

def mean_square_displacement_per_atom(positions, atom_indices):
    nframes = positions.shape[0]
    msds = np.zeros((len(atom_indices), nframes))
    for i, idx in enumerate(atom_indices):
        r0 = positions[:, idx, :]
        for dt in range(1, nframes):
            dr = r0[dt:] - r0[:-dt]
            sq_disp = np.sum(dr**2, axis=1)
            msds[i, dt] = np.mean(sq_disp)
    return msds

def fit_diffusion_coefficients(msds, dt, dim=3, fit_fraction=0.5):
    nframes = msds.shape[1]
    times = np.arange(nframes) * dt
    fit_range = int(nframes * fit_fraction)
    Ds = []
    for msd in msds:
        if np.all(msd[1:fit_range] == 0):
            continue
        coeffs = np.polyfit(times[1:fit_range], msd[1:fit_range], 1)
        D = coeffs[0] * 1e-20 / (2 * dim)  # Convert Å²/s to m²/s
        Ds.append(D)
    Ds = np.array(Ds)
    return np.mean(Ds), np.std(Ds), Ds

def correct_diffusion_yeh_hummer(D, T, eta, box_length_m):
    xi = 2.837297  # constant for cubic box
    return D + (xi * kB * T) / (6 * np.pi * eta * box_length_m)

def compute_conductivity_fixed_concentration(D, c_mol_per_m3, z, T):
    q = z * e
    n = c_mol_per_m3 * NA
    return n * q**2 * D / (kB * T)

def process_all_ani_files_fixed_conc(directory='.', dt=1e-15, dim=3,
                                     target_element='Na', T=300, z=1,
                                     eta=0.00089, box_length_ang=25.0,
                                     c_mol_per_liter=1.0):
    results = []
    ani_files = glob.glob(os.path.join(directory, '*.ANI'))
    box_length_m = box_length_ang * 1e-10
    c_mol_per_m3 = c_mol_per_liter * 1000

    for file in ani_files:
        try:
            positions, species = read_ani(file)
            atom_indices = get_atom_indices_by_element(species, target_element)
            if not atom_indices:
                print(f"[!] No {target_element} atoms in {file}")
                continue
            msds = mean_square_displacement_per_atom(positions, atom_indices)
            D_raw, D_std, Ds = fit_diffusion_coefficients(msds, dt, dim)
            D_corr = correct_diffusion_yeh_hummer(D_raw, T, eta, box_length_m)
            sigma_raw = compute_conductivity_fixed_concentration(D_raw, c_mol_per_m3, z, T)
            sigma_corr = compute_conductivity_fixed_concentration(D_corr, c_mol_per_m3, z, T)

            results.append([os.path.basename(file), target_element,
                            D_raw, D_corr, D_std, sigma_raw, sigma_corr])
            print(f"[✓] {file}: D = {D_raw:.2e}, D_corr = {D_corr:.2e}, "
                  f"σ_raw = {sigma_raw:.2e}, σ_corr = {sigma_corr:.2e} S/m")

        except Exception as e:
            print(f"[✗] Error in {file}: {e}")
            continue

    with open('diffusion_results_fixed_conc.csv', 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['File', 'Element', 'D_raw (m^2/s)', 'D_corr (m^2/s)', 'D_std',
                         'Conductivity_raw (S/m)', 'Conductivity_corr (S/m)'])
        writer.writerows(results)

    print("\n✅ Results saved to 'diffusion_results_fixed_conc.csv'")

# === Main execution ===
if __name__ == "__main__":
    dt_fs = 1.0
    dt = dt_fs * 1e-15
    process_all_ani_files_fixed_conc(
        directory='.',
        dt=dt,
        dim=3,
        target_element='Na',
        T=300,
        z=1,
        eta=0.00089,
        box_length_ang=25.0,
        c_mol_per_liter=1.0
    )

