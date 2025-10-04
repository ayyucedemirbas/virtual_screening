import os
import time
import csv
from rdkit import Chem
from lib.boltz_utils import write_boltz_input_yaml, run_boltz_predict, parse_boltz_output
from lib.docking_utils import dock_ligand_to_protein

# Configs
PROTEIN_FASTA = "data/protein.fasta"
LIGAND_SMI = "data/ligands.smi"
OUT_DIR = "outputs/boltz_with_docking"
TOP_K = 3  # how many top hits to dock/refine
EXTRA_BOLTZ_ARGS = ["--use_msa_server"]
VINA_EXE = "vina"  # or full path to vina binary

def load_protein_sequence(fasta_path: str) -> str:
    with open(fasta_path) as f:
        lines = f.readlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">"))
    return seq

def load_ligands(smiles_path: str):
    ligs = []
    with open(smiles_path) as f:
        for ln in f:
            parts = ln.strip().split()
            if not parts:
                continue
            smi = parts[0]
            ligs.append(smi)
    return ligs

def main():
    protein_seq = load_protein_sequence(PROTEIN_FASTA)
    ligands = load_ligands(LIGAND_SMI)
    os.makedirs(OUT_DIR, exist_ok=True)

    results = []
    for idx, smi in enumerate(ligands):
        print(f"[{idx}] Ligand: {smi}")
        yaml_in = os.path.join(OUT_DIR, f"input_{idx}.yaml")
        write_boltz_input_yaml(protein_seq, smi, yaml_in)

        subdir = os.path.join(OUT_DIR, f"run_{idx}")
        os.makedirs(subdir, exist_ok=True)

        try:
            json_out = run_boltz_predict(yaml_in, subdir, EXTRA_BOLTZ_ARGS)
            res = parse_boltz_output(json_out)
            res["ligand_smiles"] = smi
            res["boltz_json"] = json_out

            # If Boltz output includes a predicted complex PDB, store its path
            # Here we assume `res["structure_file"]` is a PDB or PDB-like path
            if res.get("structure_file"):
                res["complex_pdb"] = os.path.join(subdir, res["structure_file"])
            else:
                res["complex_pdb"] = None

            results.append(res)
        except Exception as e:
            print("Error running Boltz for ligand", smi, ":", e)
        time.sleep(0.1)

    # Filter valid predictions
    valid = [r for r in results if r.get("affinity_prob") is not None]
    # Sort by descending probability, then ascending affinity value
    valid.sort(key=lambda r: (-r["affinity_prob"], r.get("affinity_value", float("inf"))))

    # Save summary CSV (before docking)
    summary_csv = os.path.join(OUT_DIR, "summary_before_docking.csv")
    with open(summary_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand", "affinity_prob", "affinity_value", "complex_pdb", "boltz_json"])
        for r in valid:
            writer.writerow([
                r["ligand_smiles"],
                r["affinity_prob"],
                r.get("affinity_value", ""),
                r.get("complex_pdb", ""),
                r.get("boltz_json", "")
            ])

    # Take top K hits and run docking refinement
    top_hits = valid[:TOP_K]
    print(f"\nRefining top {len(top_hits)} hits by docking:")

    refined = []
    for rank_idx, hit in enumerate(top_hits):
        smi = hit["ligand_smiles"]
        pdb_in = hit["complex_pdb"]
        print(f"  Hit #{rank_idx}: {smi}, Boltz PDB: {pdb_in}")

        if pdb_in is None:
            print("    No complex structure from Boltz — skipping docking for this hit.")
            continue

        # Convert SMILES to RDKit Mol
        ligand_mol = Chem.MolFromSmiles(smi)
        if ligand_mol is None:
            print("    Failed to parse ligand SMILES into RDKit Mol — skipping.")
            continue

        # Define output PDB file path for docking result
        docked_pdb = os.path.join(OUT_DIR, f"docked_hit_{rank_idx}.pdb")

        try:
            print("    Running docking refinement...")
            docked = dock_ligand_to_protein(pdb_in, ligand_mol, docked_pdb, vina_exe=VINA_EXE)
            hit["docked_pdb"] = docked
            refined.append(hit)
            print("    Docking finished, output:", docked)
        except Exception as e:
            print("    Docking failed for hit:", smi, "error:", e)

    # Save refined results CSV
    refined_csv = os.path.join(OUT_DIR, "refined_hits.csv")
    with open(refined_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand", "affinity_prob", "affinity_value", "boltz_json", "complex_pdb", "docked_pdb"])
        for r in refined:
            writer.writerow([
                r["ligand_smiles"],
                r["affinity_prob"],
                r.get("affinity_value", ""),
                r.get("boltz_json", ""),
                r.get("complex_pdb", ""),
                r.get("docked_pdb", "")
            ])

    print("\nDone. Summary files:")
    print("  Before docking:", summary_csv)
    print("  After docking:", refined_csv)


if __name__ == "__main__":
    main()
