# Boltz-2 Virtual Screening + Docking Pipeline

This repository implements a prototype **virtual screening pipeline** that uses **Boltz-2** for protein–ligand co-folding + binding affinity prediction, and optionally refines the top hits via molecular docking (e.g. AutoDock Vina).  

The pipeline covers:  
- YAML input generation for Boltz  
- Running Boltz predictions via CLI  
- Parsing outputs (affinity, probability, predicted complex structures)  
- Ranking/filtering ligands  
- Docking refinement on top hits  
- Output summaries (CSV, PDBs)  

---

## 📂 Project Structure

```

virtual_screening/
├── lib/
│   ├── boltz_utils.py        # helpers for Boltz input, execution, and parsing
│   ├── docking_utils.py      # docking / ligand–protein refinement routines
├── data/
│   ├── protein.fasta          # input protein sequence
│   └── ligands.smi            # SMILES list of ligand candidates
├── outputs/
│   ├── boltz_with_docking/    # output directory (JSONs, PDBs, CSV summaries)
│   └── …                       # other output subdirs per run
├── run_screen.py              # main driver script (Boltz + docking integration)
└── README.md                  # this file

```

---

## 🚀 Installation & Setup

### Prerequisites

- Python ≥ 3.8  
- GPU + CUDA (recommended) for efficient Boltz inference  
- RDKit (for handling ligands)  
- AutoDock Vina (or equivalent docking software) for refinement  
- Boltz package installed (with appropriate extensions)  

You can install Boltz (with GPU support) as:

```bash
pip install boltz[cuda] -U
````

Or from source (for the latest updates):

```bash
git clone https://github.com/jwohlwend/boltz.git
cd boltz
pip install -e .[cuda]
```

Make sure **vina** (or your docking tool) is in your `PATH`, or configure its full path in `run_screen.py`.

---

## 🧾 Usage

1. **Prepare input data**

   * `data/protein.fasta`: your target protein sequence in FASTA format
   * `data/ligands.smi`: one ligand SMILES per line (optionally with ligand names)

2. **Run the pipeline**

   ```bash
   python _run_screen.py
   ```

   This will:

   * Generate YAML inputs per ligand
   * Run `boltz predict` per ligand
   * Parse JSON outputs (affinity, probability, structure)
   * Rank ligands by predicted binding
   * For the top **K** hits, run docking refinement
   * Produce summary CSVs and output PDB files

3. **Inspect outputs**

   * `outputs/boltz_with_docking/summary_before_docking.csv`
   * `outputs/boltz_with_docking/refined_hits.csv`
   * Complex PDBs from Boltz, and docked PDBs from refinement

---

## 🔧 Configuration Options

In `run_screen.py` you can adjust:

| Parameter          | Purpose                                             | Default                |
| ------------------ | --------------------------------------------------- | ---------------------- |
| `TOP_K`            | Number of top hits to dock                          | 3                      |
| `EXTRA_BOLTZ_ARGS` | Extra CLI flags for Boltz (e.g. `--use_msa_server`) | `["--use_msa_server"]` |
| `VINA_EXE`         | Executable name or path for docking tool            | `"vina"`               |

You should also adjust docking box size, receptor/ligand preparation, constraints etc. in `docking_utils.py` depending on your system.

---

## 🧪 Example Toy Data

You can test the pipeline with the sample toy data included:

* `data/protein.fasta`:

  ```text
  >ExampleProtein1
  MKTAYIAKQRQISFVKSHFSRQDILDAVRDTDTVPQLVTAVLLAPEDYLMRLVASA
  ```

* `data/ligands.smi`:

  ```text
  CCO ethanol
  CCOCC diethyl_ether
  CC(=O)O acetic_acid
  C1=CC=CC=C1 benzene
  C1=CC=CN=C1 pyridine
  ```

This minimal set helps verify that the I/O, YAML generation, Boltz invocation, JSON parsing, ranking, and docking integration all run end-to-end (though binding predictions will likely not be meaningful with such toy data).

---

## ⚠️ Caveats & Limitations

* **Boltz output format changes**: The JSON keys used in parsing (e.g. `affinity_pred_value`, `affinity_probability_binary`) may differ depending on version; always inspect a sample output and adjust `parse_boltz_output()` accordingly.
* **Docking refinement is schematic**: `docking_utils.py` contains a simplified docking wrapper. You may need to adjust ligand/receptor prep, box definition, etc., for your specific target.
* **Performance and scaling**: Running Boltz per ligand can be slow and GPU-intensive. For large ligand libraries, you may want batching, parallelism, or hybrid heuristics (e.g. docking-guided filtering).
* **Biological realism**: The toy protein and ligands are artificial. For true screening, use realistic protein targets and curated ligand sets.
* **Prediction confidence**: Always inspect model confidence, visual quality of predicted complexes, and avoid overreliance on single metric.
* **Model & license**: Boltz-2 is open source under the MIT license.

---

## 🔍 Further Reading & Extensions

* **Boltz-2 technical paper & architecture** — see the original manuscript detailing how Boltz-2 co-folds and predicts affinity
* **Boltzina: Docking-guided screening** — a hybrid method that speeds up screening by combining docking poses and Boltz affinity prediction
* You can extend the pipeline by:

  * Adding **batch inference** via Boltz’s API
  * Incorporating **template / contact constraints** to guide structure prediction
  * Re-ranking with MM/GBSA or FEP methods on top hits
  * Visualizing results (e.g. PyMOL, seaborn affinity plots)
  * Integrating with generative molecular design workflows

---

## 🧮 Citation

If you use this pipeline or Boltz-2 results in your work, please cite:

> Passaro, S., Corso, G., Wohlwend, J. et al. *Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction.* (2025) (https://www.biorxiv.org/content/10.1101/2025.06.14.659707v1)

---

## 🧰 Acknowledgments

This project is a demonstration of integrating Boltz-2 with docking refinement, inspired by current trends in AI-driven molecular modeling. Thanks to the Boltz authors for open access to the model and documentation.

---
