#===================================================================
# Ionic Dynamics Analysis from SIESTA .ANI Trajectory
# Includes:
#   • Mean Square Displacement (MSD)
#   • Velocity Auto-Correlation Function (VACF)
#   • Vibrational Density of States (VDOS)
#   • Diffusion coefficients (Einstein & Green–Kubo)
#   • Ionic conductivity (Nernst–Einstein)
#===================================================================

import numpy as np
import matplotlib.pyplot as plt

# ================= USER PARAMETERS =================
ani_file  = "peg6peo5nai.ANI"
cell_file = "result.out"
dt_fs = 1.0                 # timestep in femtoseconds
T = 300                     # temperature (K)
z = 1                       # ionic charge (Na+)
ion_symbol = "Na"           # mobile ion species
max_vacf_lag = 1000

# ================= CONSTANTS =================
dt = dt_fs * 1e-15          # timestep (s)
kB = 1.380649e-23           # J/K
e  = 1.602176634e-19        # C

# ================= FUNCTIONS =================

def read_ani(filename):
    positions, species = [], []
    with open(filename) as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        nat = int(lines[i])
        i += 2
        frame = []
        if not species:
            species = []
            for j in range(nat):
                parts = lines[i+j].split()
                species.append(parts[0])
        for j in range(nat):
            parts = lines[i+j].split()
            frame.append(list(map(float, parts[1:4])))
        positions.append(frame)
        i += nat
    return np.array(positions), species

def remove_com_drift(pos):
    com = np.mean(pos, axis=1, keepdims=True)
    return pos - com

def compute_velocities(pos, dt):
    return (pos[1:] - pos[:-1]) / dt   # Å/s

def compute_msd(pos):
    nframes, natoms, _ = pos.shape
    msd = np.zeros(nframes)
    for tau in range(1, nframes):
        disp = pos[tau:] - pos[:-tau]
        msd[tau] = np.mean(np.sum(disp**2, axis=2))
    return msd

def diffusion_from_msd(msd, dt):
    t = np.arange(len(msd)) * dt
    t1 = int(0.4 * len(t))
    t2 = int(0.9 * len(t))
    slope, _ = np.polyfit(t[t1:t2], msd[t1:t2], 1)
    D_ang2_s = slope / 6.0
    return D_ang2_s * 1e-20

def compute_vacf(vel, max_lag):
    nframes, natoms, _ = vel.shape
    vacf = np.zeros(max_lag)
    for lag in range(max_lag):
        dot = 0.0
        count = 0
        for t in range(nframes - lag):
            dot += np.sum(np.sum(vel[t] * vel[t+lag], axis=1))
            count += natoms
        vacf[lag] = dot / count
    return vacf

def diffusion_from_vacf(vacf, dt):
    D_ang2_s = np.trapz(vacf, dx=dt) / 3.0
    return D_ang2_s * 1e-20

def compute_vdos(vacf, dt):
    vacf -= np.mean(vacf[-50:])
    window = np.hanning(len(vacf))
    fft = np.fft.fft(vacf * window)
    freq = np.fft.fftfreq(len(vacf), d=dt)
    return freq[:len(freq)//2] * 1e-12, np.real(fft[:len(freq)//2])

def read_cell(filename):
    with open(filename) as f:
        lines = f.readlines()
    for i, l in enumerate(lines):
        if "unit cell vectors" in l.lower():
            a = list(map(float, lines[i+1].split()))
            b = list(map(float, lines[i+2].split()))
            c = list(map(float, lines[i+3].split()))
            return np.array([a, b, c])
    raise RuntimeError("Cell vectors not found.")

def conductivity(D, N, V, z, T):
    n = N / V
    q = z * e
    return n * q*q * D / (kB * T)

# ================= MAIN =================

print("📥 Reading trajectory...")
positions, species = read_ani(ani_file)

ion_indices = [i for i, s in enumerate(species) if s == ion_symbol]
positions = positions[:, ion_indices, :]
N_ions = len(ion_indices)

positions = remove_com_drift(positions)

print("📊 MSD analysis...")
msd = compute_msd(positions)
D_msd = diffusion_from_msd(msd, dt)

print("📈 VACF analysis...")
vel = compute_velocities(positions, dt)
vacf = compute_vacf(vel, min(max_vacf_lag, vel.shape[0]))
D_vacf = diffusion_from_vacf(vacf, dt)

print("🎵 VDOS...")
freqs, vdos = compute_vdos(vacf, dt)

cell = read_cell(cell_file)
V_m3 = abs(np.linalg.det(cell)) * 1e-30

sigma = conductivity(D_vacf, N_ions, V_m3, z, T)

# ================= OUTPUT =================

t_ps = np.arange(len(msd)) * dt * 1e12

plt.figure()
plt.plot(t_ps, msd)
plt.xlabel("Time (ps)")
plt.ylabel("MSD (Å²)")
plt.grid()
plt.savefig("msd.png")
plt.close()

plt.figure()
plt.plot(np.arange(len(vacf))*dt*1e12, vacf)
plt.xlabel("Time (ps)")
plt.ylabel("VACF (Å²/s²)")
plt.grid()
plt.savefig("vacf.png")
plt.close()

plt.figure()
plt.plot(freqs, vdos)
plt.xlabel("Frequency (THz)")
plt.ylabel("VDOS (a.u.)")
plt.grid()
plt.savefig("vdos.png")
plt.close()

np.savetxt("msd.txt", np.column_stack([t_ps, msd]), header="Time(ps) MSD(Å^2)")
np.savetxt("vacf.txt", np.column_stack([np.arange(len(vacf))*dt*1e12, vacf]), header="Time(ps) VACF(Å^2/s^2)")
np.savetxt("vdos.txt", np.column_stack([freqs, vdos]), header="Frequency(THz) VDOS")

print("\n✅ FINAL SUMMARY")
print(f" • Mobile ions          : {N_ions}")
print(f" • Diffusion (MSD)      : {D_msd:.3e} m²/s")
print(f" • Diffusion (VACF)     : {D_vacf:.3e} m²/s")
print(f" • Conductivity (NE)    : {sigma:.3e} S/m")
print("📁 Files saved: msd.png, vacf.png, vdos.png + .txt data")

