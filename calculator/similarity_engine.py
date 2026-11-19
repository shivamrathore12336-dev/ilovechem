"""
ECFP4 Morgan Tanimoto Similarity & Bemis-Murcko Scaffold-Hopping Engine.
"""
from typing import Any, Dict, List, Optional
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold


def _murcko_smiles(mol) -> str:
    try:
        core = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(core, canonical=True) if core and core.GetNumAtoms() > 0 else ""
    except Exception:
        return ""


def compare_library_to_reference(reference_smiles: str, candidate_lines: List[str], max_items: int = 100) -> Optional[Dict[str, Any]]:
    ref_mol = Chem.MolFromSmiles((reference_smiles or "").strip()[:2000])
    if ref_mol is None:
        return None

    ref_fp = AllChem.GetMorganFingerprintAsBitVect(ref_mol, radius=2, nBits=2048)
    ref_scaffold = _murcko_smiles(ref_mol)

    ranked = []
    for idx, line in enumerate(candidate_lines[:max_items], start=1):
        raw = (line or "").strip()
        if not raw or raw.startswith("#"):
            continue
        parts = raw.split(None, 1)
        smi = parts[0]
        name = parts[1].strip() if len(parts) > 1 else f"Candidate_{idx}"
        mol = Chem.MolFromSmiles(smi[:2000])
        if mol is None:
            continue
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        sim = round(float(DataStructs.TanimotoSimilarity(ref_fp, fp)), 4)
        scaff = _murcko_smiles(mol)
        same_scaffold = bool(ref_scaffold and scaff == ref_scaffold)

        if sim >= 0.85:
            tier = "Close Analog"
        elif sim >= 0.55 and not same_scaffold:
            tier = "Scaffold Hop Candidate"
        elif sim >= 0.55:
            tier = "Same-Core Series"
        else:
            tier = "Structurally Diverse"

        ranked.append({
            "name": name,
            "smiles": Chem.MolToSmiles(mol, canonical=True),
            "tanimoto_ecfp4": sim,
            "same_murcko_scaffold": same_scaffold,
            "murcko_scaffold": scaff,
            "classification": tier,
        })

    ranked.sort(key=lambda r: r["tanimoto_ecfp4"], reverse=True)
    return {
        "reference_smiles": Chem.MolToSmiles(ref_mol, canonical=True),
        "reference_scaffold": ref_scaffold,
        "compared_count": len(ranked),
        "results": ranked,
    }
