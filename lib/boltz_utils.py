import os
import subprocess
import json
from pathlib import Path
from typing import List, Dict

def write_boltz_input_yaml(protein_seq: str, ligand_smiles: str, yaml_path: str):
    """
    Creates a minimal YAML input for Boltz-2 to predict binding affinity.
    """
    # This mirrors example from Boltz docs/blogs. :contentReference[oaicite:1]{index=1}
    data = {
        "version": 1,
        "sequences": [
            {
                "protein": {
                    "id": "P0",
                    "sequence": protein_seq
                }
            },
            {
                "ligand": {
                    "id": "L0",
                    "smiles": ligand_smiles
                }
            }
        ],
        "properties": [
            {
                "affinity": {
                    "binder": "L0"
                }
            }
        ]
    }
    import yaml
    with open(yaml_path, "w") as f:
        yaml.dump(data, f)

def run_boltz_predict(yaml_input: str, output_dir: str, extra_args: List[str] = None):
    """
    Calls the `boltz predict` CLI command and returns the path to output JSON.
    """
    cmd = ["boltz", "predict", yaml_input, "--out_dir", output_dir]
    if extra_args:
        cmd.extend(extra_args)
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    # We expect JSON under output_dir/predictions/…
    # Let's search for the .json file
    for root, dirs, files in os.walk(output_dir):
        for fn in files:
            if fn.endswith(".json"):
                return os.path.join(root, fn)
    raise RuntimeError("No JSON output found in Boltz output")

def parse_boltz_output(json_path: str) -> Dict:
    """
    Parses Boltz output JSON and returns relevant fields (affinity, probability, etc).
    """
    with open(json_path) as f:
        data = json.load(f)
    # Inspect the JSON structure; typical fields include `affinity_pred_value`, `affinity_probability_binary`
    # from Boltz docs. :contentReference[oaicite:2]{index=2}
    # Example:
    result = {
        "affinity_value": data.get("affinity_pred_value", None),
        "affinity_prob": data.get("affinity_probability_binary", None),
        "structure_file": data.get("complex_structure", None)  # hypothetical
    }
    return result
