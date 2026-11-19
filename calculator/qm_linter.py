"""
Pre-Flight Quantum Chemistry Input Linter & SLURM Job Script Generator.
Catches steric clashes (<0.65 A), charge/spin parity errors, and missing HPC resources
before users submit jobs to ORCA 6 or Gaussian 16.
"""
import re
from typing import Any, Dict, List, Tuple
import numpy as np

ATOMIC_NUMBERS = {
    "H": 1, "HE": 2, "LI": 3, "BE": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "NE": 10,
    "NA": 11, "MG": 12, "AL": 13, "SI": 14, "P": 15, "S": 16, "CL": 17, "AR": 18,
    "K": 19, "CA": 20, "FE": 26, "CO": 27, "NI": 28, "CU": 29, "ZN": 30, "BR": 35, "I": 53,
}

def get_z_number(sym: str) -> int:
    """Resolve atomic number for any element 1-118 via RDKit PeriodicTable, with fallback."""
    clean = sym.strip().capitalize()
    try:
        from rdkit import Chem
        return Chem.GetPeriodicTable().GetAtomicNumber(clean)
    except Exception:
        return ATOMIC_NUMBERS.get(clean.upper(), 6)


def lint_qm_geometry(
    atoms: List[Tuple[str, float, float, float]],
    charge: int = 0,
    multiplicity: int = 1,
    engine: str = "orca",
    nprocs: int = 8,
) -> Dict[str, Any]:
    """Validate geometry clashes, electron/spin parity, and estimate HPC memory."""
    warnings: List[str] = []
    errors: List[str] = []

    if not atoms:
        return {"valid": False, "errors": ["No atomic coordinates found."], "warnings": [], "slurm_script": ""}

    symbols = [a[0].strip().upper() for a in atoms]
    coords = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
    n_atoms = len(atoms)

    # 1. Spin / Charge Parity Check
    total_z = sum(get_z_number(sym) for sym in symbols)
    n_electrons = total_z - int(charge)
    if multiplicity < 1:
        errors.append(f"Invalid spin multiplicity ({multiplicity}): must be >= 1.")
    elif multiplicity > n_electrons + 1:
        errors.append(f"Spin multiplicity ({multiplicity}) exceeds maximum possible for {n_electrons} electrons.")
    elif (n_electrons % 2 == 0) == (multiplicity % 2 == 0):
        suggested = 1 if (n_electrons % 2 == 0) else 2
        errors.append(
            f"Spin/Charge Parity Mismatch: {n_electrons} electrons (charge {charge}) cannot have "
            f"multiplicity {multiplicity}. Suggested multiplicity: {suggested}."
        )

    # 2. Pairwise Interatomic Clash Detection (< 0.65 A causes SCF divergence)
    clashes = []
    if n_atoms <= 1500:
        diff = coords[:, None, :] - coords[None, :, :]
        dist_mat = np.sqrt(np.sum(diff ** 2, axis=-1))
        np.fill_diagonal(dist_mat, 999.0)
        i_idx, j_idx = np.where(np.triu(dist_mat) < 0.65)
        for i, j in zip(i_idx, j_idx):
            if i < j:
                clashes.append(f"{symbols[i]}{i+1}\u2013{symbols[j]}{j+1} ({dist_mat[i, j]:.3f} \u00c5)")
    if clashes:
        errors.append(f"Severe atomic clashes (< 0.65 \u00c5) detected: {', '.join(clashes[:5])}.")

    # 3. HPC Resource Estimation
    recommended_mem_gb = max(4, min(64, int(round(n_atoms * 0.25)) + 4))
    maxcore_mb = int((recommended_mem_gb * 1024 * 0.75) // max(1, nprocs))
    if n_atoms > 150:
        warnings.append(f"Large system ({n_atoms} atoms): consider RIJCOSX / def2-SVP pre-optimization before triple-zeta DFT.")

    slurm_script = build_slurm_script(engine=engine, nprocs=nprocs, mem_gb=recommended_mem_gb)

    return {
        "valid": len(errors) == 0,
        "atom_count": n_atoms,
        "electron_count": n_electrons,
        "recommended_mem_gb": recommended_mem_gb,
        "recommended_orca_maxcore_mb": maxcore_mb,
        "errors": errors,
        "warnings": warnings,
        "slurm_script": slurm_script,
    }


def build_slurm_script(engine: str = "orca", job_name: str = "ilovechem_qm", nprocs: int = 8, mem_gb: int = 16) -> str:
    safe_name = re.sub(r"[^A-Za-z0-9_\-]", "_", job_name)[:40] or "ilovechem_qm"
    if engine.lower() == "gaussian":
        run_cmd = f"g16 < {safe_name}.gjf > {safe_name}.log"
        mod_cmd = "module load gaussian/16"
    else:
        run_cmd = f"$(which orca) {safe_name}.inp > {safe_name}.out"
        mod_cmd = "module load orca/6.0 openmpi"

    return (
        "#!/bin/bash\n"
        f"#SBATCH --job-name={safe_name}\n"
        f"#SBATCH --ntasks={max(1, min(int(nprocs), 64))}\n"
        f"#SBATCH --mem={max(2, min(int(mem_gb), 256))}G\n"
        "#SBATCH --time=24:00:00\n"
        f"#SBATCH --output={safe_name}_%j.slurm.log\n\n"
        "set -euo pipefail\n"
        f"{mod_cmd}\n"
        f"{run_cmd}\n"
    )
