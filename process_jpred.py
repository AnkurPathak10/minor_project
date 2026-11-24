import os
import glob
import re
import pandas as pd
import matplotlib.pyplot as plt

RESULTS_DIR = "jp_batch_1795450214_JPred_RESTful_API_job_ALL_JOBS_ARCHIVE"

# Find all prediction files
jnet_files = glob.glob(os.path.join(RESULTS_DIR, "**", "*.jnet"), recursive=True)
if not jnet_files:
    jnet_files = glob.glob(os.path.join(RESULTS_DIR, "**", "*.concise"), recursive=True)

print(f"Found {len(jnet_files)} prediction files.")

records = []

for file_path in sorted(jnet_files):
    amp_id = None
    pred_seq = None
    orig_seq = None
    
    with open(file_path, "r", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                amp_id = line[1:].strip()
            elif "jnetpred:" in line:
                pred_seq = line.split("jnetpred:")[1].strip()
            elif line.startswith("jnetpred,"):
                pred_seq = line.split(",")[1].strip()
            elif "origseq:" in line:
                orig_seq = line.split("origseq:")[1].strip()
            elif line.startswith("origseq,"):
                orig_seq = line.split(",")[1].strip()

    if pred_seq:
        # Strip all commas, spaces, or extra punctuation
        clean_pred = re.sub(r'[^HEC\-]', '', pred_seq)
        
        # Strip commas from original sequence if present
        if orig_seq:
            orig_seq = re.sub(r'[^A-Z]', '', orig_seq.upper())

        total_len = len(clean_pred)
        if total_len > 0:
            h_count = clean_pred.count('H')
            e_count = clean_pred.count('E')
            c_count = clean_pred.count('-') + clean_pred.count('C')

            records.append({
                "Job_Folder": os.path.basename(os.path.dirname(file_path)),
                "Header_ID": amp_id,
                "Sequence": orig_seq,
                "Length": total_len,
                "Helix_Count": h_count,
                "Sheet_Count": e_count,
                "Coil_Count": c_count,
                "Alpha_Helix_%": round((h_count / total_len) * 100, 2),
                "Beta_Sheet_%": round((e_count / total_len) * 100, 2),
                "Random_Coil_%": round((c_count / total_len) * 100, 2),
                "Secondary_Structure_Prediction": clean_pred
            })

df_jpred = pd.DataFrame(records)

# 2. Merge with your ProtParam Excel sheet
excel_path = "AMP_50_ProtParam_Characterization.xlsx"
output_excel = "AMP_50_Characterization_with_SecondaryStructure.xlsx"

if os.path.exists(excel_path) and not df_jpred.empty:
    print(f"Reading {excel_path}...")
    df_existing = pd.read_excel(excel_path)
    
    # Clean sequences for matching
    df_existing["Seq_clean"] = df_existing["Sequence"].astype(str).str.replace(r'[^A-Z]', '', regex=True).str.upper()
    df_jpred["Seq_clean"] = df_jpred["Sequence"].astype(str).str.replace(r'[^A-Z]', '', regex=True).str.upper()
    
    # Merge based on sequence match
    cols_to_add = ["Seq_clean", "Alpha_Helix_%", "Beta_Sheet_%", "Random_Coil_%", "Secondary_Structure_Prediction"]
    
    # If sequences matched properly
    if df_jpred["Seq_clean"].isin(df_existing["Seq_clean"]).sum() > 30:
        df_merged = pd.merge(df_existing, df_jpred[cols_to_add], on="Seq_clean", how="left").drop(columns=["Seq_clean"])
    else:
        # Fallback to sequential row merge
        cols_direct = ["Alpha_Helix_%", "Beta_Sheet_%", "Random_Coil_%", "Secondary_Structure_Prediction"]
        df_merged = pd.concat([df_existing.drop(columns=["Seq_clean"], errors="ignore").reset_index(drop=True),
                               df_jpred[cols_direct].reset_index(drop=True)], axis=1)

    df_merged.to_excel(output_excel, index=False)
    print(f"Saved cleanly merged file to: {output_excel}")

# 3. Generate 100% Stacked Bar Chart with clean AMP01-AMP50 labels
if not df_jpred.empty:
    # Use Accession or sequential AMP labels for clean X-axis
    plot_labels = [f"AMP{i+1:02d}" for i in range(len(df_jpred))]
    df_jpred["Plot_Label"] = plot_labels

    plt.figure(figsize=(16, 7))
    plot_df = df_jpred.set_index("Plot_Label")[["Alpha_Helix_%", "Beta_Sheet_%", "Random_Coil_%"]]
    
    ax = plot_df.plot(kind="bar", stacked=True, figsize=(16, 7), 
                      color=["#2b8a3e", "#1971c2", "#e8590c"], width=0.8)
    
    plt.title("Secondary Structure Composition across 50 AMPs (JPred4)", fontsize=14, fontweight="bold", pad=15)
    plt.xlabel("Antimicrobial Peptides (AMP01 – AMP50)", fontsize=11, labelpad=10)
    plt.ylabel("Structural Composition (%)", fontsize=11)
    plt.legend(["% Alpha Helix (H)", "% Beta Sheet (E)", "% Random Coil (-)"], 
               loc="upper right", frameon=True)
    plt.xticks(rotation=90, fontsize=8)
    plt.ylim(0, 100)
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig("AMP50_Secondary_Structure_Distribution.png", dpi=300)
    plt.close()
    print("Regenerated clean plot: AMP50_Secondary_Structure_Distribution.png")