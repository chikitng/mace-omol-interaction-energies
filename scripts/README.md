# Scripts

## mlp_minimise_restrained_torsions.py

Minimise one SDF conformer with strong restraints on all rotatable torsions and report the MLP energy with restraint forces excluded.

### Example

```bash
/home/campus.ncl.ac.uk/nfc78/miniforge3/envs/mace-omol/bin/python scripts/mlp_minimise_restrained_torsions.py \
  --sdf 3QTU/10_conf/3QTU/conformers/ligand_10confs.sdf \
  --conformer-index 0 \
  --platform CUDA \
  --minimised-sdf-out scripts/minimised_structure.sdf \
  --text-out scripts/minimisation_result.txt
```

### Useful options

- `--model`: MLP model name (default: `mace-off23-medium`).
- `--molecule-index`: Molecule index in the SDF (default: `0`).
- `--conformer-index`: Conformer index for the selected molecule (default: `0`).
- `--restraint-k`: Torsion restraint force constant in kcal mol-1 rad-2 (default: `100000.0`).
- `--max-iterations`: OpenMM minimisation max iterations (`0` means until convergence).
- `--platform`: OpenMM platform (`CUDA`, `OpenCL`, `HIP`, `CPU`, `Reference`).
- `--minimised-sdf-out`: Optional path for minimised-geometry SDF output.
- `--json-out`: Optional path for JSON output.
- `--text-out`: Optional path for plain-text output.
