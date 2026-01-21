 
AIMD Post-Processing Toolkit for SIESTA
From Stability to Structure, Dynamics, and Transport
This repository contains a modular set of Python scripts for post-processing Ab Initio Molecular Dynamics (AIMD) simulations performed with SIESTA.
The toolkit enables a complete and physically consistent workflow:
MDE → RDF & Coordination → MSD/VACF → Diffusion → Ionic Conductivity
The scripts are designed for polymer electrolytes, solid electrolytes, and ion-conducting materials, and are suitable for research workflows, hands-on practice (HoP) sessions, PhD theses, and journal publications.

Repository Structure
.
├── mlv_script_ani_xyz_converter.py
├── mlv_script_plot_mde_dual_axis.py
├── mlv_script_rdf_cn_analysis.py
├── mlv_script_msd_analysis.py
├── mlv_script_vacf_vdos_diffusion_cond_final.py
├── mlv_script_diffusion_coeff_cond_atomic_area.py
│
├── example_inputs/
│   ├── peg6peo5nai.ANI
│   ├── peg6peo5nai.MDE
│   ├── peg6peo5nai.STRUCT_OUT
│   └── result.out
│
└── README.md

Requirements
    • Python 3.8 or newer
    • NumPy
    • Matplotlib
    • ASE (Atomic Simulation Environment)
    • tqdm
Install dependencies using:
pip install numpy matplotlib ase tqdm

Recommended Analysis Workflow
1. AIMD Stability Check (Mandatory)
Script:
mlv_script_plot_mde_dual_axis.py
Purpose
    • Reads the .MDE file
    • Plots temperature and total energy versus MD steps using dual Y-axes
Why
    • Validates correct NVE or NVT ensemble behavior
    • Ensures the AIMD trajectory is thermodynamically stable before further analysis
Usage
python mlv_script_plot_mde_dual_axis.py peg6peo5nai.MDE
Output
mde_dual_axis_plot.png

2. Trajectory Conversion (.ANI to .XYZ)
Script:
mlv_script_ani_xyz_converter.py
Purpose
    • Converts SIESTA .ANI trajectories into .xyz format
    • Preserves atom ordering and frame integrity
Usage
python mlv_script_ani_xyz_converter.py
Output
peg6peo5nai_from_ani.xyz

3. Structural Analysis: RDF and Coordination Numbers
Script:
mlv_script_rdf_cn_analysis.py
Purpose
    • Computes radial distribution functions (RDFs) for selected element pairs
    • Integrates RDFs to obtain coordination numbers
    • Automatically assigns simulation cell vectors from STRUCT_OUT
Typical Outputs
rdf_results/
├── rdf_O_Na.txt
├── rdf_O_H.txt
└── rdf_O_O.txt

rdf_all_pairs.png
coordination_numbers.txt
Physical Insight
    • RDFs describe local structural ordering
    • Coordination numbers quantify ion–environment binding

4. Mean Square Displacement (Quick MSD Check)
Script:
mlv_script_msd_analysis.py
Purpose
    • Computes MSD averaged over all atoms
    • Useful for quick validation and instructional purposes
Usage
python mlv_script_msd_analysis.py
Output
MSD_vs_Time.png

5. Full Ionic Dynamics and Transport Analysis
Script:
mlv_script_vacf_vdos_diffusion_cond_final.py
Includes
    • MSD using the Einstein relation
    • VACF using the Green–Kubo formalism
    • Vibrational density of states (VDOS)
    • Diffusion coefficients
    • Ionic conductivity using the Nernst–Einstein relation
Key Features
    • Ion-specific analysis (e.g., Li⁺ or Na⁺)
    • Center-of-mass drift removal
    • Consistent unit handling
    • Publication-ready plots
Outputs
msd.png
vacf.png
vdos.png
msd.txt
vacf.txt
vdos.txt

6. Diffusion and Conductivity with Yeh–Hummer Correction
Script:
mlv_script_diffusion_coeff_cond_atomic_area.py
Purpose
    • Computes per-ion MSD-based diffusion coefficients
    • Applies Yeh–Hummer finite-size correction
    • Computes conductivity at fixed concentration
    • Outputs a CSV summary file
Output
diffusion_results_fixed_conc.csv

Physical Interpretation Guide
Step
Key Question
MDE
Is the AIMD simulation stable?
RDF
How are ions locally coordinated?
CN
How strong is ion binding?
MSD/VACF
How fast do ions move?
Conductivity
How efficiently is charge transported?
Structure controls dynamics, and dynamics controls transport.

Important Notes
    • .ANI, .XYZ, and MD_CAR files do not contain lattice vectors
Always supply STRUCT_OUT or result.out
    • Reported conductivities represent Nernst–Einstein upper bounds
    • Ion–ion correlations (Haven ratio < 1) are expected, especially in polymers
    • Always discard equilibration frames before statistical analysis

Intended Use
    • Polymer electrolytes
    • Solid-state electrolytes
    • Hands-on practice and training courses
    • PhD and MS theses
    • Journal publications and supporting information

Citation and Attribution
If you use this toolkit in academic work, please acknowledge:
Dr. Mohan L. Verma
Email: drmohanlv@gmail.com

Contributions and Extensions
Contributions are welcome, including:
    • Li⁺ versus Na⁺ comparative workflows
    • Multivalent ion support
    • Error estimation and block averaging
    • Arrhenius analysis
    • Integration with NEB or advanced transport models

Summary
A complete and physically consistent AIMD post-processing toolkit for SIESTA, bridging stability, structure, dynamics, and transport.
 
