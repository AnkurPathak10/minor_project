import os
import time
import requests
import pandas as pd

# ----------------- CONFIGURATION -----------------
API_TOKEN = "4e2dc8a9da8fd8cff360cc4f27bae06467dc8fcd"  # Keep your token here
FASTA_FILE = "AMP50_JPred_Batch.fasta"
OUTPUT_DIR = "swissmodel_pdbs"
SUMMARY_CSV = "swissmodel_modelling_summary.csv"
BASE_URL = "https://swissmodel.expasy.org"
# -------------------------------------------------

os.makedirs(OUTPUT_DIR, exist_ok=True)
headers = {
    "Authorization": f"Token {API_TOKEN}",
    "Content-Type": "application/json"
}

def parse_fasta(filepath):
    peptides = []
    current_id = None
    current_seq = []
    with open(filepath, "r") as f:
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
    return peptides

peptides = parse_fasta(FASTA_FILE)
print(f"Loaded {len(peptides)} antimicrobial peptides from {FASTA_FILE}")

results = []

for amp_id, seq in peptides:
    print(f"\n--- Submitting {amp_id} ({len(seq)} aa) ---")
    
    # SWISS-MODEL requires 'target_sequences' as a list: [seq]
    payload = {
        "target_sequences": [seq],
        "project_title": f"AMP_50_{amp_id}"
    }

    try:
        response = requests.post(f"{BASE_URL}/automodel", json=payload, headers=headers)
        
        if response.status_code not in (200, 202):
            print(f"  [!] Failed to submit: {response.status_code} {response.text}")
            results.append({"AMP_ID": amp_id, "Status": f"Failed ({response.status_code})", "GMQE": None, "QMEAN": None, "PDB_File": None})
            continue

        res_json = response.json()
        project_id = res_json.get("project_id")
        print(f"  [+] Job accepted! Project ID: {project_id}")

        # Polling loop (max 3 minutes per peptide)
        status = "RUNNING"
        model_data = None
        max_attempts = 18  # 18 * 10s = 180 seconds
        attempts = 0

        while status in ("RUNNING", "SUBMITTED", "WAITING") and attempts < max_attempts:
            time.sleep(10)
            attempts += 1
            
            status_resp = requests.get(f"{BASE_URL}/project/{project_id}/models/summary/", headers=headers)
            if status_resp.status_code == 200:
                summary = status_resp.json()
                status = summary.get("status")
                print(f"      Polling status ({attempts * 10}s): {status}")
                
                if status == "COMPLETED":
                    # Fetch full model details to get download links & scores
                    detail_resp = requests.get(f"{BASE_URL}/project/{project_id}/models/", headers=headers)
                    if detail_resp.status_code == 200:
                        model_data = detail_resp.json()
                    break
                elif status in ("FAILED", "NO_TEMPLATES"):
                    print(f"  [-] Finished with status: {status}")
                    break
            else:
                print(f"      Status query returned {status_resp.status_code}")

        # Process model if completed
        if status == "COMPLETED" and model_data and model_data.get("models"):
            top_model = model_data["models"][0]
            gmqe = top_model.get("gmqe")
            qmean = top_model.get("qmean", {}).get("qmean4_norm")
            pdb_url = top_model.get("coordinates_url")

            pdb_filename = f"{amp_id}_swissmodel.pdb"
            if pdb_url:
                pdb_resp = requests.get(pdb_url)
                pdb_path = os.path.join(OUTPUT_DIR, pdb_filename)
                with open(pdb_path, "w") as f_pdb:
                    f_pdb.write(pdb_resp.text)
                print(f"  [✓] Model saved: {pdb_filename} | GMQE: {gmqe} | QMEAN: {qmean}")
                
            results.append({
                "AMP_ID": amp_id,
                "Status": "Success",
                "GMQE": gmqe,
                "QMEAN": qmean,
                "PDB_File": pdb_filename
            })
        else:
            final_status = status if status in ("NO_TEMPLATES", "FAILED") else "Timeout/No Model"
            print(f"  [-] {amp_id}: {final_status}")
            results.append({
                "AMP_ID": amp_id,
                "Status": final_status,
                "GMQE": None,
                "QMEAN": None,
                "PDB_File": None
            })

    except Exception as e:
        print(f"  [!] Exception for {amp_id}: {e}")
        results.append({"AMP_ID": amp_id, "Status": f"Error: {e}", "GMQE": None, "QMEAN": None, "PDB_File": None})

# Save cumulative summary
df_out = pd.DataFrame(results)
df_out.to_csv(SUMMARY_CSV, index=False)
print(f"\nCompleted! Results logged in {SUMMARY_CSV}")