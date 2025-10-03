from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
from pathlib import Path

def dock_ligand_to_protein(protein_pdb: str, ligand_mol: Chem.Mol, out_pdb: str, vina_exe: str = "vina"):
    """
    Runs docking (e.g. AutoDock Vina) to refine ligand binding pose.
    """
    # Convert ligand to PDBQT or relevant format
    ligand_pdb = out_pdb.replace(".pdb", "_ligand.pdb")
    Chem.MolToPDBFile(ligand_mol, ligand_pdb)
    # (You would need to convert to PDBQT, prepare receptor, set search box, etc.)
    # For simplicity, assume receptor and ligand PDBQT are ready:
    ligand_q = ligand_pdb.replace(".pdb", ".pdbqt")
    protein_q = protein_pdb.replace(".pdb", ".pdbqt")
    cmd = [vina_exe,
           "--receptor", protein_q,
           "--ligand", ligand_q,
           "--out", out_pdb,
           "--center_x", "0","--center_y","0","--center_z","0",  # dummy
           "--size_x", "20","--size_y","20","--size_z","20"]
    subprocess.run(cmd, check=True)
    return out_pdb
