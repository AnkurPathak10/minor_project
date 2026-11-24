import os
import shutil
from gradio_client import Client

FASTA_FILE = "AMP50_JPred_Batch.fasta"
OUTPUT_DIR = "tertiary_structures_pdb"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 1. Parse your FASTA file
peptides = []
current_id, current_seq = None, []

with open(FASTA_FILE, "r") as f:
    for line in f:
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                peptides.append((current_id, "".join(current_seq)))
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    if current_id:
        peptides.append((current_id, "".join(current_seq)))

print(f"Loaded {len(peptides)} AMP sequences.")

# 2. Connect to the ESMFold Space
print("Connecting to ESMFold server...")
client = Client("wwydmanski/esmfold")

# Predict for the first 5 AMPs
N_MODELS = 5

for i in range(N_MODELS):
    amp_id, seq = peptides[i]
    out_file = os.path.join(OUTPUT_DIR, f"{amp_id}.pdb")
    print(f"[{i+1}/{N_MODELS}] Folding {amp_id} ({len(seq)} aa)...", end=" ", flush=True)

    try:
        # Run structure prediction
        result_path = client.predict(seq, api_name="/pdb")
        
        # Copy generated PDB to output directory
        shutil.copy(result_path, out_file)
        print(f"Success! Saved -> {out_file}")
    except Exception as e:
        print(f"Error: {e}")

print("\nDone! Check the 'tertiary_structures_pdb' folder for your .pdb files.")