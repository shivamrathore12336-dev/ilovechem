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
    esol_logs = compute_esol_logs(mol, mw, logp, rot_bonds)
    formula_label = rdMolDescriptors.CalcMolFormula(mol)

    return {
        "smiles": clean_smi,
        "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
        "formula": formula_label,
        "molecular_weight": round(mw, 2),
        "exact_mass": round(exact_mass, 4),
        "logp": round(logp, 2),
        "tpsa": round(tpsa, 2),
        "esol_logs": esol_logs,
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
        "boiled_egg_svg": generate_boiled_egg_svg(tpsa, logp, label=formula_label),
        "bioavailability_radar_svg": generate_bioavailability_radar_svg(mw, logp, tpsa, rot_bonds, fsp3, esol_logs, label=formula_label),
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


def compute_esol_logs(mol: Any, mw: float, logp: float, rot_bonds: int) -> float:
    """
    Delaney ESOL (Estimated Aqueous Solubility, LogS in mol/L).
    LogS = 0.16 - (0.63 * LogP) - (0.0062 * MW) + (0.066 * RB) - (0.74 * AP)
    where AP = aromatic proportion = (number of aromatic atoms / number of heavy atoms).
    """
    if mol is None:
        return 0.0
    try:
        heavy_atoms = mol.GetNumHeavyAtoms()
        arom_atoms = sum(1 for a in mol.GetAtoms() if a.GetIsAromatic())
        ap = arom_atoms / max(heavy_atoms, 1)
        logs = 0.16 - (0.63 * logp) - (0.0062 * mw) + (0.066 * rot_bonds) - (0.74 * ap)
        return round(float(logs), 2)
    except Exception:
        return 0.0


def generate_bioavailability_radar_svg(
    mw: float,
    logp: float,
    tpsa: float,
    rot_bonds: int,
    fsp3: float,
    logs: float,
    label: str = "Query",
) -> str:
    """
    SwissADME-style 6-axis polygonal Bioavailability Radar chart.
    Axes & Optimal Drug-like Boundaries:
    1. LIPO:   -0.7 <= XLogP <= 5.0    (Range: -2.0 .. 8.0)
    2. SIZE:   150 <= MW <= 500 Da     (Range: 50 .. 750)
    3. POLAR:  20 <= TPSA <= 130 Å²    (Range: 0 .. 180)
    4. INSOLU: -6.0 <= LogS <= 0.0     (Range: -10.0 .. 2.0)
    5. INSATU: 0.25 <= Fsp3 <= 1.0     (Range: 0.0 .. 1.0)
    6. FLEX:   0 <= RotB <= 9          (Range: 0 .. 16)
    """
    cx, cy = 160.0, 145.0
    r_max = 95.0

    axes = [
        ("LIPO", logp, -0.7, 5.0, -2.0, 8.0, ""),
        ("SIZE", mw, 150.0, 500.0, 50.0, 750.0, "Da"),
        ("POLAR", tpsa, 20.0, 130.0, 0.0, 180.0, "Å²"),
        ("INSOLU", logs, -6.0, 0.0, -10.0, 2.0, ""),
        ("INSATU", fsp3, 0.25, 1.0, 0.0, 1.0, ""),
        ("FLEX", float(rot_bonds), 0.0, 9.0, 0.0, 16.0, ""),
    ]

    def norm(val: float, v_min: float, v_max: float) -> float:
        clamped = max(v_min, min(v_max, val))
        return (clamped - v_min) / (v_max - v_min)

    angles = [-math.pi / 2.0 + i * (2.0 * math.pi / 6.0) for i in range(6)]

    # Concentric background polygon rings
    grid_svg = []
    for step in [0.25, 0.5, 0.75, 1.0]:
        pts = [f"{cx + step * r_max * math.cos(a):.1f},{cy + step * r_max * math.sin(a):.1f}" for a in angles]
        grid_svg.append(f'<polygon points="{" ".join(pts)}" fill="none" stroke="#e2e8f0" stroke-width="1"/>')

    # Spoke lines and text labels
    spokes_svg = []
    for i, ((name, val, _, _, _, _, unit), a) in enumerate(zip(axes, angles)):
        x_end = cx + r_max * math.cos(a)
        y_end = cy + r_max * math.sin(a)
        spokes_svg.append(
            f'<line x1="{cx:.1f}" y1="{cy:.1f}" x2="{x_end:.1f}" y2="{y_end:.1f}" '
            f'stroke="#cbd5e1" stroke-width="1" stroke-dasharray="2,2"/>'
        )

        lx = cx + (r_max + 18.0) * math.cos(a)
        ly = cy + (r_max + 18.0) * math.sin(a)
        anchor = "middle"
        if math.cos(a) > 0.3:
            anchor = "start"
        elif math.cos(a) < -0.3:
            anchor = "end"
        val_str = f"{val:.1f}{unit}" if isinstance(val, float) else f"{int(val)}"
        spokes_svg.append(
            f'<text x="{lx:.1f}" y="{ly:.1f}" text-anchor="{anchor}" '
            f'font-family="IBM Plex Sans, sans-serif" font-size="10" font-weight="700" fill="#334155">{name}'
            f'<tspan font-size="9" font-weight="400" fill="#64748b"> ({val_str})</tspan></text>'
        )

    # Optimal Zone polygon (Pink translucent sweet spot)
    opt_outer = []
    for (name, _, omin, omax, amin, amax, _), a in zip(axes, angles):
        r_out = (0.15 + 0.85 * norm(omax, amin, amax)) * r_max
        opt_outer.append((cx + r_out * math.cos(a), cy + r_out * math.sin(a)))

    opt_path_pts = " ".join([f"{x:.1f},{y:.1f}" for x, y in opt_outer])
    opt_poly = (
        f'<polygon points="{opt_path_pts}" fill="#fce7f3" fill-opacity="0.65" '
        f'stroke="#f472b6" stroke-width="1.5" stroke-dasharray="3,3"/>'
    )

    # Current molecule polygon & vertex markers
    mol_pts = []
    mol_dots = []
    violations = []
    for (name, val, omin, omax, amin, amax, unit), a in zip(axes, angles):
        if val < omin or val > omax:
            violations.append(name)
        r_val = (0.15 + 0.85 * norm(val, amin, amax)) * r_max
        mx = cx + r_val * math.cos(a)
        my = cy + r_val * math.sin(a)
        mol_pts.append(f"{mx:.1f},{my:.1f}")
        dot_color = "#dc2626" if (val < omin or val > omax) else "#2563eb"
        mol_dots.append(
            f'<circle cx="{mx:.1f}" cy="{my:.1f}" r="4" fill="{dot_color}" stroke="#ffffff" stroke-width="1.5"/>'
        )

    mol_poly = f'<polygon points="{" ".join(mol_pts)}" fill="#3b82f6" fill-opacity="0.30" stroke="#1d4ed8" stroke-width="2"/>'

    status_text = "Optimal Oral Space" if len(violations) == 0 else f"Outside: {', '.join(violations[:2])}"
    status_color = "#059669" if len(violations) == 0 else "#d97706"

    return (
        f'<svg viewBox="0 0 320 295" width="100%" role="img" '
        f'aria-label="Bioavailability Radar Chart" xmlns="http://www.w3.org/2000/svg">'
        f'<rect width="320" height="295" rx="12" fill="#ffffff" stroke="#e2e8f0" stroke-width="1"/>'
        f'<text x="16" y="24" font-family="IBM Plex Sans, sans-serif" font-size="12" font-weight="700" fill="#0f172a">Bioavailability Radar</text>'
        f'<text x="304" y="24" text-anchor="end" font-family="IBM Plex Sans, sans-serif" font-size="10.5" font-weight="700" fill="{status_color}">{status_text}</text>'
        f'{"".join(grid_svg)}'
        f'{opt_poly}'
        f'{"".join(spokes_svg)}'
        f'{mol_poly}'
        f'{"".join(mol_dots)}'
        f'<text x="160" y="285" text-anchor="middle" font-family="IBM Plex Sans, sans-serif" font-size="9.5" fill="#94a3b8">Pink zone: Oral Bioavailability Space (SwissADME criteria)</text>'
        f'</svg>'
    )
