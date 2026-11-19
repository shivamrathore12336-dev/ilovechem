"""
Bioisostere & Medicinal Chemistry Lead Optimizer for iLoveChem.
Applies curated RDKit AllChem SMARTS reaction transforms to propose
higher-QED, metabolically stabilized, and bioisosteric drug analogs in RAM.
"""
from typing import Any, Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from calculator.medchem_engine import analyze_medchem_profile

BIOISOSTERE_TRANSFORMS = [
    (
        "Carboxylic Acid -> Tetrazole",
        "Improves membrane permeability & metabolic stability while retaining acidity (pKa ~4.9).",
        "[*:1][C](=O)[OH]>>[*:1]c1nnn[nH]1",
    ),
    (
        "Carboxylic Acid -> Primary Amide",
        "Removes negative charge to boost Blood-Brain Barrier (BBB) penetration.",
        "[C:1](=O)[OH]>>[C:1](=O)N",
    ),
    (
        "Ester -> Amide (Hydrolysis Shield)",
        "Replaces metabolically labile ester with plasma-stable amide bond.",
        "[C:1](=O)O[C:2]>>[C:1](=O)N[C:2]",
    ),
    (
        "Aromatic C-H -> C-F (Metabolic Blocker)",
        "Blocks CYP450 para-hydroxylation soft spots with minimal steric penalty.",
        "[cH:1]1[c:2][c:3][c:4][c:5][c:6]1>>[c:1](F)1[c:2][c:3][c:4][c:5][c:6]1",
    ),
    (
        "Phenyl -> Pyridine (Polarity Boost)",
        "Lowers LogP by ~1.0–1.5 units to improve aqueous solubility & reduce hERG risk.",
        "[cH:1]1[c:2][c:3][c:4][c:5][c:6]1>>[n:1]1[c:2][c:3][c:4][c:5][c:6]1",
    ),
    (
        "Aryl Methyl -> Trifluoromethyl (CF3)",
        "Prevents benzylic oxidation while increasing lipophilic pocket binding.",
        "[c:1][CH3]>>[c:1]C(F)(F)F",
    ),
    (
        "Phenol -> Fluorobenzene",
        "Eliminates Phase II glucuronidation/sulfation clearance site.",
        "[c:1][OH]>>[c:1]F",
    ),
]


def generate_bioisostere_analogs(smiles: str, max_analogs: int = 12) -> Optional[Dict[str, Any]]:
    parent = analyze_medchem_profile(smiles)
    if parent is None:
        return None

    parent_mol = Chem.MolFromSmiles(parent["canonical_smiles"])
    if parent_mol is None:
        return None

    seen = {parent["canonical_smiles"]}
    analogs: List[Dict[str, Any]] = []

    for rule_name, rationale, rxn_smarts in BIOISOSTERE_TRANSFORMS:
        try:
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((parent_mol,))
        except Exception:
            continue

        added_for_rule = 0
        for prod_tuple in products:
            if added_for_rule >= 2:
                break
            prod = prod_tuple[0]
            try:
                Chem.SanitizeMol(prod)
                can_smi = Chem.MolToSmiles(prod, canonical=True)
            except Exception:
                continue
            if can_smi in seen:
                continue
            seen.add(can_smi)
            prof = analyze_medchem_profile(can_smi)
            if prof is None:
                continue
            prof.pop("boiled_egg_svg", None)
            analogs.append({
                "rule": rule_name,
                "rationale": rationale,
                "smiles": can_smi,
                "formula": prof["formula"],
                "mw": prof["molecular_weight"],
                "logp": prof["logp"],
                "tpsa": prof["tpsa"],
                "qed": prof["qed_score"],
                "delta_qed": round(prof["qed_score"] - parent["qed_score"], 3),
                "delta_logp": round(prof["logp"] - parent["logp"], 2),
                "delta_tpsa": round(prof["tpsa"] - parent["tpsa"], 2),
                "pains_count": prof["alerts"]["pains_count"],
                "bbb_permeant": prof["boiled_egg"]["bbb_permeant"],
            })
            added_for_rule += 1

    # Rank analogs by clean PAINS status first, then highest QED gain
    analogs.sort(key=lambda x: (x["pains_count"] > 0, -x["qed"]))
    parent.pop("boiled_egg_svg", None)
    return {
        "parent": parent,
        "analog_count": len(analogs[:max_analogs]),
        "analogs": analogs[:max_analogs],
    }
