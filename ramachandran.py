import os
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PPBuilder

# 1. Point to your downloaded PDB file
PDB_FILE = "model_1.pdb"  # or model1.pdb from HDOCK/PEP-FOLD
OUTPUT_IMAGE = "AMP01_Ramachandran_Plot.png"

parser = PDBParser(QUIET=True)
structure = parser.get_structure("AMP01", PDB_FILE)
ppb = PPBuilder()

phi_psi_angles = []

for pp in ppb.build_peptides(structure):
    angles = pp.get_phi_psi_list()  # returns list of (phi, psi) in radians
    for phi, psi in angles:
        if phi is not None and psi is not None:
            phi_psi_angles.append((np.degrees(phi), np.degrees(psi)))

phi = [p[0] for p in phi_psi_angles]
psi = [p[1] for p in phi_psi_angles]

# 2. Plotting
fig, ax = plt.subplots(figsize=(7, 7), dpi=300)

# Draw reference grid & axes
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axvline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)

# Background typical favored zones (Alpha-helix & Beta-sheet)
ax.add_patch(plt.Rectangle((-180, 90), 120, 90, color="#d4edda", alpha=0.6, label="Favored (β-sheet)"))
ax.add_patch(plt.Rectangle((-120, -90), 80, 70, color="#cfe2ff", alpha=0.6, label="Favored (Right α-helix)"))
ax.add_patch(plt.Rectangle((40, 20), 40, 60, color="#fff3cd", alpha=0.6, label="Favored (Left α-helix)"))

# Scatter residue dihedral angles
ax.scatter(phi, psi, color="#d90429", edgecolors="black", s=65, zorder=5, label=f"Residues (N={len(phi)})")

ax.set_title("Ramachandran Plot: AMP01 Backbone Validation", fontsize=14, fontweight="bold", pad=15)
ax.set_xlabel(r"$\phi$ (Phi) Dihedral Angle [°]", fontsize=12)
ax.set_ylabel(r"$\psi$ (Psi) Dihedral Angle [°]", fontsize=12)
ax.set_xticks(np.arange(-180, 181, 60))
ax.set_yticks(np.arange(-180, 181, 60))
ax.grid(True, linestyle=":", alpha=0.5)
ax.legend(loc="upper right", frameon=True, fontsize=9)

plt.tight_layout()
plt.savefig(OUTPUT_IMAGE)
print(f"Ramachandran plot generated successfully -> {OUTPUT_IMAGE}")