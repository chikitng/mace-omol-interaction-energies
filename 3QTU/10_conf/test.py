#!/usr/bin/env python3

import csv
import re
import shutil
import subprocess
from pathlib import Path
from urllib.request import urlretrieve

from rdkit import Chem
from rdkit.Chem import AllChem


# =========================
# Settings
# =========================
PDB_ID = "3QTU"
TARGET_FOLDER = "009-CDK2"
PLREX_BASE = "https://raw.githubusercontent.com/Honza-R/PL-REX/main"

NUM_CONFS = 10
RANDOM_SEED = 2026
# PRUNE_RMS = 0.5
PRUNE_RMS = -1

BASE_DIR = Path(PDB_ID)
INPUT_DIR = BASE_DIR / "input"
CONF_DIR = BASE_DIR / "conformers"
MOPAC_DIR = BASE_DIR / "mopac"
RESULTS_DIR = BASE_DIR / "results"

INPUT_SDF = INPUT_DIR / "ligand.sdf"
CONF_SDF = CONF_DIR / "ligand_10confs.sdf"
RESULTS_CSV = RESULTS_DIR / "vacuum_results.csv"


# =========================
# Setup directories
# =========================
def setup_dirs():
    INPUT_DIR.mkdir(parents=True, exist_ok=True)
    CONF_DIR.mkdir(parents=True, exist_ok=True)
    MOPAC_DIR.mkdir(parents=True, exist_ok=True)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# =========================
# Download ligand from PL-REX
# =========================
def download_ligand():
    url = f"{PLREX_BASE}/{TARGET_FOLDER}/structures_pl-rex/{PDB_ID}/ligand.sdf"
    print(f"Downloading ligand from:\n  {url}")
    urlretrieve(url, INPUT_SDF)
    print(f"Saved to: {INPUT_SDF}")


# =========================
# Read charge
# =========================
def get_charge(mol: Chem.Mol) -> int:
    if mol.HasProp("charge"):
        return int(float(mol.GetProp("charge").strip()))
    return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())


# =========================
# Generate conformers
# =========================
def generate_conformers() -> int:
    suppl = Chem.SDMolSupplier(str(INPUT_SDF), removeHs=False)
    mol = suppl[0]
    if mol is None:
        raise ValueError(f"Could not read molecule from {INPUT_SDF}")

    charge = get_charge(mol)
    print(f"Using charge: {charge}")

    mol = Chem.AddHs(mol, addCoords=True)

    params = AllChem.ETKDGv3()
    params.randomSeed = RANDOM_SEED
    params.pruneRmsThresh = PRUNE_RMS
    params.useSmallRingTorsions = True
    params.enforceChirality = True

    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=NUM_CONFS, params=params))
    if not conf_ids:
        raise RuntimeError("No conformers were generated.")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        print("Optimizing conformers with MMFF...")
        for cid in conf_ids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid)
    else:
        print("MMFF unavailable, optimizing conformers with UFF...")
        for cid in conf_ids:
            AllChem.UFFOptimizeMolecule(mol, confId=cid)

    writer = Chem.SDWriter(str(CONF_SDF))
    for cid in conf_ids:
        writer.write(mol, confId=cid)
    writer.close()

    print(f"Generated {len(conf_ids)} conformers -> {CONF_SDF}")
    return charge


# =========================
# Write MOPAC input files
# Header:
# PM7 CHARGE=0
#
# Vacuum
# =========================
def build_mopac_keywords(charge: int) -> str:
    return f"PM7 CHARGE={charge}"


def mol_to_mopac_xyz_block(mol: Chem.Mol, conf_id: int = 0) -> str:
    conf = mol.GetConformer(conf_id)
    lines = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(
            f"{atom.GetSymbol():2s} {pos.x:12.6f} 1 {pos.y:12.6f} 1 {pos.z:12.6f} 1"
        )
    return "\n".join(lines)


def write_mopac_inputs(charge: int) -> list[Path]:
    mols = [m for m in Chem.SDMolSupplier(str(CONF_SDF), removeHs=False) if m is not None]
    keywords = build_mopac_keywords(charge)

    mop_files = []
    for i, mol in enumerate(mols, start=1):
        coord_block = mol_to_mopac_xyz_block(mol, conf_id=0)

        text = (
            f"{keywords}\n"
            f"\n"
            f"\n"
            f"{coord_block}\n"
            f"\n"
        )

        mop_path = MOPAC_DIR / f"conf_{i:02d}.mop"
        mop_path.write_text(text)
        mop_files.append(mop_path)

    print(f"Wrote {len(mop_files)} MOPAC input files -> {MOPAC_DIR}")
    return mop_files


# =========================
# Run MOPAC
# =========================
def find_mopac_executable() -> str:
    for exe in ["mopac", "MOPAC2016.exe", "MOPAC.exe"]:
        path = shutil.which(exe)
        if path:
            return path
    raise FileNotFoundError("Could not find MOPAC executable in PATH.")


def run_mopac_jobs(mop_files: list[Path]):
    mopac_exe = find_mopac_executable()
    print(f"Using MOPAC executable: {mopac_exe}")

    for mop_file in mop_files:
        print(f"Running {mop_file.name}")
        subprocess.run(
            [mopac_exe, mop_file.name],
            cwd=str(MOPAC_DIR),
            check=True
        )


# =========================
# Parse .out files
# =========================
def extract_heat_of_formation(text: str):
    patterns = [
        re.compile(r"FINAL\s+HEAT\s+OF\s+FORMATION\s*=\s*([-+0-9.Ee]+)", re.I),
        re.compile(r"HEAT\s+OF\s+FORMATION\s*=\s*([-+0-9.Ee]+)", re.I),
    ]
    for pat in patterns:
        m = pat.search(text)
        if m:
            return float(m.group(1))
    return None


def extract_total_energy(text: str):
    patterns = [
        re.compile(r"TOTAL\s+ENERGY\s*=\s*([-+0-9.Ee]+)", re.I),
        re.compile(r"ELECTRONIC\s+ENERGY\s*=\s*([-+0-9.Ee]+)", re.I),
    ]
    for pat in patterns:
        m = pat.search(text)
        if m:
            return float(m.group(1))
    return None


def extract_status(text: str) -> str:
    low = text.lower()
    if "job ended normally" in low or "== mopac done ==" in low:
        return "ok"
    if "unable to achieve self-consistence" in low:
        return "scf_failed"
    if "too many cycles" in low:
        return "too_many_cycles"
    if "error" in low:
        return "error"
    if "failed" in low:
        return "failed"
    return "unknown"


def parse_outputs():
    rows = []

    for out_file in sorted(MOPAC_DIR.glob("conf_*.out")):
        text = out_file.read_text(errors="ignore")
        status = extract_status(text)
        hof = extract_heat_of_formation(text)
        total_energy = extract_total_energy(text)
        rows.append([out_file.name, status, hof, total_energy])

    with open(RESULTS_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "status", "heat_of_formation", "total_energy"])
        writer.writerows(rows)

    print(f"Results written to: {RESULTS_CSV}")
    print("\nSummary:")
    for row in rows:
        print(
            f"{row[0]:15s} status={row[1]:10s} "
            f"heat_of_formation={row[2]} total_energy={row[3]}"
        )


# =========================
# Main
# =========================
def main():
    print(f"Running in current directory: {Path.cwd()}")
    setup_dirs()
    download_ligand()
    charge = generate_conformers()
    mop_files = write_mopac_inputs(charge)
    run_mopac_jobs(mop_files)
    parse_outputs()


if __name__ == "__main__":
    main()