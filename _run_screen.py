import os
from lib.boltz_utils import write_boltz_input_yaml, run_boltz_predict, parse_boltz_output

PROTEIN_FASTA = "data/protein.fasta"
LIGAND_SMI = "data/ligands.smi"
OUT_DIR = "outputs/boltz_test"
EXTRA_BOLTZ_ARGS = ["--use_msa_server"]
TOP_K = 3

def load_protein_sequence(path):
    with open(path) as f:
        lines = f.readlines()
    seq = "".join(l.strip() for l in lines if not l.startswith(">"))
    return seq

def load_ligands(path):
    ligs = []
    with open(path) as f:
        for ln in f:
            parts = ln.strip().split()
            if len(parts) == 0:
                continue
            smi = parts[0]
            ligs.append(smi)
    return ligs

def main():
    protein_seq = load_protein_sequence(PROTEIN_FASTA)
    ligands = load_ligands(LIGAND_SMI)
    os.makedirs(OUT_DIR, exist_ok=True)

    all_res = []
    for i, smi in enumerate(ligands):
        print(f"Testing ligand {i}: {smi}")
        yaml_in = os.path.join(OUT_DIR, f"in_{i}.yaml")
        write_boltz_input_yaml(protein_seq, smi, yaml_in)
        sub = os.path.join(OUT_DIR, f"run_{i}")
        os.makedirs(sub, exist_ok=True)
        try:
            json_out = run_boltz_predict(yaml_in, sub, EXTRA_BOLTZ_ARGS)
            r = parse_boltz_output(json_out)
            r["ligand"] = smi
            all_res.append(r)
        except Exception as e:
            print("Error for", smi, ":", e)

    # Filter & sort
    valid = [r for r in all_res if r.get("affinity_prob") is not None]
    valid.sort(key=lambda r: (-r["affinity_prob"], r["affinity_value"] or 0))
    print("Top results:")
    for r in valid[:TOP_K]:
        print(r)

if __name__ == "__main__":
    main()
