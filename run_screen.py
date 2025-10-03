import os
import time
from rdkit import Chem
from lib.boltz_utils import write_boltz_input_yaml, run_boltz_predict, parse_boltz_output

# Configs
PROTEIN_FASTA = "data/protein.fasta"
LIGAND_SMI = "data/ligands.smi"
OUT_DIR = "outputs/boltz_results"
TOP_K = 10  # keep top 10
EXTRA_BOLTZ_ARGS = ["--use_msa_server"]  # e.g. requiring MSA server

def load_protein_sequence(fasta_path: str) -> str:
    from Bio import SeqIO
    rec = next(SeqIO.parse(fasta_path, "fasta"))
    return str(rec.seq)

def load_ligands(smiles_path: str):
    ligs = []
    with open(smiles_path) as f:
        for line in f:
            smi = line.strip().split()[0]
            ligs.append(smi)
    return ligs

def main():
    protein_seq = load_protein_sequence(PROTEIN_FASTA)
    ligands = load_ligands(LIGAND_SMI)
    os.makedirs(OUT_DIR, exist_ok=True)

    results = []
    for idx, smi in enumerate(ligands):
        print(f"Processing ligand {idx}/{len(ligands)}: {smi}")
        yaml_in = os.path.join(OUT_DIR, f"input_{idx}.yaml")
        write_boltz_input_yaml(protein_seq, smi, yaml_in)
        sub_out = os.path.join(OUT_DIR, f"run_{idx}")
        os.makedirs(sub_out, exist_ok=True)
        try:
            out_json = run_boltz_predict(yaml_in, sub_out, EXTRA_BOLTZ_ARGS)
            res = parse_boltz_output(out_json)
            res["ligand"] = smi
            res["json"] = out_json
            results.append(res)
        except Exception as e:
            print("Boltz failed for ligand", smi, ":", e)
        time.sleep(0.1)  # small throttle

    # Sort by predicted affinity (lower = stronger binder, or by probability)
    # Note: depending on Boltz’s convention, affinity_pred_value could be something like ΔG (so smaller is better)
    # We'll sort by affinity_prob descending, then by affinity_value ascending.
    results = [r for r in results if r.get("affinity_prob") is not None]
    results.sort(key=lambda r: (-r["affinity_prob"], r["affinity_value"] if r["affinity_value"] is not None else float("inf")))

    top_hits = results[:TOP_K]
    print("Top hits:")
    for r in top_hits:
        print(f" Ligand {r['ligand']}: prob={r['affinity_prob']:.3f}, val={r['affinity_value']:.3f}")

    # Optionally: for top hits, take the output structure and refine via docking
    # (You’d need to extract the PDB file path from Boltz or reconstruct from JSON)
    # ...

    # Save results
    import csv
    with open(os.path.join(OUT_DIR, "summary.csv"), "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["ligand", "affinity_prob", "affinity_value", "json_path"])
        for r in results:
            writer.writerow([r["ligand"], r["affinity_prob"], r["affinity_value"], r["json"]])

if __name__ == "__main__":
    main()
