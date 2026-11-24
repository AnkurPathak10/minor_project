import os
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

FASTA_FILE = "AMP50_JPred_Batch.fasta"

# 1. Parse sequences
sequences = []
current_seq = []
with open(FASTA_FILE, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_seq:
                sequences.append("".join(current_seq))
                current_seq = []
        else:
            current_seq.append(line)
    if current_seq:
        sequences.append("".join(current_seq))

# 20 standard amino acids
amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

# Calculate percentage for each peptide, then take mean across 50 AMPs
all_freqs = []
for seq in sequences:
    length = len(seq)
    counts = Counter(seq)
    freq = {aa: (counts[aa] / length) * 100 for aa in amino_acids}
    all_freqs.append(freq)

mean_freqs = {aa: np.mean([f[aa] for f in all_freqs]) for aa in amino_acids}
std_freqs = {aa: np.std([f[aa] for f in all_freqs]) for aa in amino_acids}

# Sort by mean frequency descending
sorted_aas = sorted(amino_acids, key=lambda x: mean_freqs[x], reverse=True)
sorted_means = [mean_freqs[aa] for aa in sorted_aas]
sorted_stds = [std_freqs[aa] for aa in sorted_aas]

# ----------------- PLOT 1: Overall Amino Acid Frequency -----------------
plt.figure(figsize=(10, 5), dpi=300)
colors = ['#1d3557' if aa in 'KRH' else '#457b9d' if aa in 'AVLIM' else '#e63946' if aa in 'DE' else '#2a9d8f' for aa in sorted_aas]

bars = plt.bar(sorted_aas, sorted_means, yerr=sorted_stds, capsize=4, color=colors, edgecolor='black', alpha=0.9)
plt.title("Mean Amino Acid Frequency Across 50 Screened AMPs", fontsize=13, fontweight='bold', pad=12)
plt.xlabel("Amino Acid Residue", fontsize=11, fontweight='bold')
plt.ylabel("Mean Composition (%)", fontsize=11, fontweight='bold')
plt.grid(axis='y', linestyle='--', alpha=0.5)

# Add custom legend
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#1d3557', edgecolor='black', label='Cationic (K, R, H)'),
    Patch(facecolor='#457b9d', edgecolor='black', label='Hydrophobic (A, V, L, I, M)'),
    Patch(facecolor='#e63946', edgecolor='black', label='Anionic (D, E)'),
    Patch(facecolor='#2a9d8f', edgecolor='black', label='Polar / Other')
]
plt.legend(handles=legend_elements, loc='upper right', frameon=True)
plt.tight_layout()
plt.savefig("AMP50_Amino_Acid_Composition_Bar.png")
print("Saved -> AMP50_Amino_Acid_Composition_Bar.png")

# ----------------- PLOT 2: Chemical Class Breakdown -----------------
classes = {
    'Hydrophobic\n(A,V,I,L,M)': sum(mean_freqs[aa] for aa in 'AVLIM'),
    'Cationic\n(K,R,H)': sum(mean_freqs[aa] for aa in 'KRH'),
    'Polar Neutral\n(S,T,N,Q)': sum(mean_freqs[aa] for aa in 'STNQ'),
    'Special/Small\n(G,P)': sum(mean_freqs[aa] for aa in 'GP'),
    'Aromatic\n(F,W,Y)': sum(mean_freqs[aa] for aa in 'FWY'),
    'Anionic\n(D,E)': sum(mean_freqs[aa] for aa in 'DE'),
    'Cysteine\n(C)': mean_freqs['C']
}

plt.figure(figsize=(9, 4.8), dpi=300)
class_names = list(classes.keys())
class_vals = list(classes.values())
class_colors = ['#2a9d8f', '#1d3557', '#a8dadc', '#e9c46a', '#f4a261', '#e76f51', '#6c757d']

bars2 = plt.bar(class_names, class_vals, color=class_colors, edgecolor='black', width=0.6)
for bar in bars2:
    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2.0, yval + 0.6, f"{yval:.1f}%", ha='center', va='bottom', fontsize=9, fontweight='bold')

plt.title("Functional Chemical Group Distribution (N = 50 AMPs)", fontsize=13, fontweight='bold', pad=12)
plt.ylabel("Cohort Composition (%)", fontsize=11, fontweight='bold')
plt.ylim(0, max(class_vals) + 6)
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("AMP50_Chemical_Group_Distribution.png")
print("Saved -> AMP50_Chemical_Group_Distribution.png")