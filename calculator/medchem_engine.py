"""
Advanced Medicinal Chemistry, Structural Alert (PAINS/Brenk/NIH),
BOILED-Egg Permeability, and Batch Library Screening Engine for iLoveChem.
100% in-memory RDKit execution — zero external API calls or database storage.
"""
import html
import math
from typing import Any, Dict, List, Optional
from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors, FilterCatalog, QED, rdMolDescriptors

# Initialize RDKit FilterCatalogs once at module import (thread-safe for read-only GetMatches)
def _build_catalog(catalog_enum) -> FilterCatalog.FilterCatalog:
    params = FilterCatalog.FilterCatalogParams()
    params.AddCatalog(catalog_enum)
    return FilterCatalog.FilterCatalog(params)

_PAINS_CATALOG = _build_catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
_BRENK_CATALOG = _build_catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
_NIH_CATALOG = _build_catalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)


def classify_boiled_egg(tpsa: float, wlogp: float) -> Dict[str, Any]:
    """
    BOILED-Egg model (Daina & Zoete, ChemMedChem 2016) elliptical boundaries:
    - Yolk (BBB permeation): center (TPSA=38.117, WLogP=3.177), a=39.96, b=2.76 (approx ellipse)
    - White (HIA / GI absorption): center (TPSA=71.051, WLogP=2.292), a=71.0, b=4.15 (approx ellipse)
    """
    in_bbb = (((tpsa - 38.117) / 41.0) ** 2 + ((wlogp - 3.177) / 2.85) ** 2) <= 1.0 and tpsa <= 79.0
    in_hia = (((tpsa - 71.051) / 71.5) ** 2 + ((wlogp - 2.292) / 4.20) ** 2) <= 1.0 and tpsa <= 142.0

    if in_bbb:
        region = "BBB Permeant (Yolk)"
        gi_absorption = "High"
        bbb_permeant = True
    elif in_hia:
        region = "GI Absorbed (White)"
        gi_absorption = "High"
        bbb_permeant = False
    else:
        region = "Low Absorption / Out of Bounds"
        gi_absorption = "Low"
        bbb_permeant = False

    return {
        "region": region,
        "gi_absorption": gi_absorption,
        "bbb_permeant": bbb_permeant,
        "tpsa": round(tpsa, 2),
        "wlogp": round(wlogp, 2),
    }


def compute_ligand_efficiency(affinity_kcal: float, heavy_atoms: int, logp: Optional[float] = None) -> Dict[str, Any]:
    """
    Compute medicinal chemistry Ligand Efficiency (LE) and Lipophilic Ligand Efficiency (LLE).
    - LE = -DeltaG / N_heavy (kcal/mol per heavy atom; drug-like target >= 0.30)
    - pKi = -log10(Ki), where Ki = exp(DeltaG / (R * T)) at 298.15 K (RT = 0.59248 kcal/mol)
    - LLE = pKi - LogP (target >= 5.0)
    """
    if heavy_atoms <= 0:
        return {"le": None, "pki": None, "lle": None}
    le = round(-affinity_kcal / heavy_atoms, 3)
    rt = 0.5924847  # kcal/mol at 298.15 K
    pki = round(-affinity_kcal / (rt * math.log(10.0)), 2)
    lle = round(pki - logp, 2) if logp is not None else None
    return {
        "le": le,
        "le_pass": le >= 0.30,
        "pki": pki,
        "lle": lle,
        "lle_pass": (lle >= 5.0) if lle is not None else None,
    }


def analyze_medchem_profile(smiles: str) -> Optional[Dict[str, Any]]:
    """Run full drug-likeness, PAINS/Brenk/NIH alerts, and BOILED-Egg profile on a SMILES string."""
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return None
    clean_smi = smiles.strip()[:2000]
    mol = Chem.MolFromSmiles(clean_smi)
    if mol is None:
        return None

    mw = Descriptors.MolWt(mol)
    exact_mass = Descriptors.ExactMolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    heavy_atoms = mol.GetNumHeavyAtoms()
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    qed_val = QED.qed(mol)

    # Structural alerts
    pains_hits = [m.GetDescription() for m in _PAINS_CATALOG.GetMatches(mol)]
    brenk_hits = [m.GetDescription() for m in _BRENK_CATALOG.GetMatches(mol)]
    nih_hits = [m.GetDescription() for m in _NIH_CATALOG.GetMatches(mol)]

    # Medicinal chemistry rules
    lipinski_v = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    veber_pass = (rot_bonds <= 10) and (tpsa <= 140.0)
    pfizer_3_75_risk = (logp > 3.0) and (tpsa < 75.0)  # High promiscuity / toxicity risk
    gsk_4_400_pass = (mw <= 400.0) and (logp <= 4.0)
    boiled_egg = classify_boiled_egg(tpsa, logp)

    return {
        "smiles": clean_smi,
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "formula": rdMolDescriptors.CalcMolFormula(mol),
        "molecular_weight": round(mw, 2),
        "exact_mass": round(exact_mass, 4),
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 2),
        "h_bond_donors": hbd,
        "h_bond_acceptors": hba,
        "rotatable_bonds": rot_bonds,
        "heavy_atoms": heavy_atoms,
        "fraction_csp3": round(fsp3, 3),
        "rings": rings,
        "aromatic_rings": arom_rings,
        "qed_score": round(qed_val, 3),
        "lipinski_violations": lipinski_v,
        "lipinski_pass": lipinski_v <= 1,
        "veber_pass": veber_pass,
        "pfizer_3_75_pass": not pfizer_3_75_risk,
        "gsk_4_400_pass": gsk_4_400_pass,
        "boiled_egg": boiled_egg,
        "alerts": {
            "pains_count": len(pains_hits),
            "pains": pains_hits[:5],
            "brenk_count": len(brenk_hits),
            "brenk": brenk_hits[:5],
            "nih_count": len(nih_hits),
            "nih": nih_hits[:5],
            "clean": (len(pains_hits) == 0 and len(brenk_hits) == 0 and len(nih_hits) == 0),
        },
        "boiled_egg_svg": generate_boiled_egg_svg(tpsa, logp, label=rdMolDescriptors.CalcMolFormula(mol)),
    }


def batch_screen_smiles(lines: List[str], max_items: int = 100) -> Dict[str, Any]:
    """Screen up to 100 SMILES lines (optionally `SMILES <name>`) in memory."""
    results = []
    errors = []
    for idx, raw_line in enumerate(lines[:max_items], start=1):
        stripped = (raw_line or "").strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split(None, 1)
        smi = parts[0]
        name = parts[1].strip() if len(parts) > 1 else f"Compound_{idx}"
        profile = analyze_medchem_profile(smi)
        if profile is None:
            errors.append({"line": idx, "input": stripped[:80], "error": "Invalid SMILES"})
            continue
        profile["name"] = name
        profile.pop("boiled_egg_svg", None)  # Keep batch JSON compact
        results.append(profile)

    passed_all = sum(1 for r in results if r["lipinski_pass"] and r["alerts"]["pains_count"] == 0)
    return {
        "total_processed": len(results),
        "passed_lipinski_and_pains": passed_all,
        "failed_count": len(errors),
        "compounds": results,
        "errors": errors,
    }


def generate_boiled_egg_svg(tpsa: float, wlogp: float, label: str = "Query") -> str:
    """Render a crisp, publication-grade BOILED-Egg SVG with the query molecule plotted."""
    # Coordinate mapping: TPSA [0..180] -> X [50..370], WLogP [-2..8] -> Y [210..30]
    def tx(val: float) -> float:
        clamped = max(0.0, min(180.0, val))
        return 50.0 + (clamped / 180.0) * 320.0

    def ty(val: float) -> float:
        clamped = max(-2.0, min(8.0, val))
        return 210.0 - ((clamped + 2.0) / 10.0) * 180.0

    px, py = tx(tpsa), ty(wlogp)
    safe_label = html.escape(str(label)[:24])

    return (
        '<svg viewBox="0 0 400 250" width="100%" role="img" '
        'aria-label="BOILED-Egg BBB and GI absorption plot" xmlns="http://www.w3.org/2000/svg">'
        '<rect x="50" y="20" width="320" height="190" rx="6" fill="#f1f5f9" stroke="#cbd5e1" stroke-width="1"/>'
        # White Egg (HIA / GI Absorption)
        f'<ellipse cx="{tx(71.051):.1f}" cy="{ty(2.292):.1f}" rx="{(71.5/180.0)*320.0:.1f}" '
        f'ry="{(4.20/10.0)*180.0:.1f}" fill="#ffffff" stroke="#94a3b8" stroke-width="1.5"/>'
        # Yellow Yolk (BBB Permeation)
        f'<ellipse cx="{tx(38.117):.1f}" cy="{ty(3.177):.1f}" rx="{(41.0/180.0)*320.0:.1f}" '
        f'ry="{(2.85/10.0)*180.0:.1f}" fill="#fef08a" stroke="#eab308" stroke-width="1.5"/>'
        # Axis labels
        '<text x="210" y="238" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="11" fill="#475569">TPSA (Å²)</text>'
        '<text x="16" y="115" text-anchor="middle" font-family="IBM Plex Mono, monospace" font-size="11" fill="#475569" transform="rotate(-90 16 115)">WLogP</text>'
        '<text x="115" y="102" font-family="IBM Plex Sans, sans-serif" font-size="10" font-weight="600" fill="#854d0e">BBB (Yolk)</text>'
        '<text x="220" y="160" font-family="IBM Plex Sans, sans-serif" font-size="10" font-weight="600" fill="#475569">HIA (White)</text>'
        # Query Molecule Marker
        f'<circle cx="{px:.1f}" cy="{py:.1f}" r="6" fill="#0f766e" stroke="#ffffff" stroke-width="2"/>'
        f'<text x="{min(px + 9, 310):.1f}" y="{max(py - 8, 34):.1f}" font-family="IBM Plex Mono, monospace" '
        f'font-size="10.5" font-weight="600" fill="#0f172a">{safe_label}</text>'
        '</svg>'
    )
