"""
iLoveChem ChemDock v3.2 - Publication-Grade 2D Interaction Diagram & PyMOL Session Generator
PoseView / Maestro / LigPlot / PyMOL / ChimeraX Standard
100% Volatile in-memory execution (Pharma Zero-Storage Invariant)
"""

import math
import re
from typing import List, Dict, Any, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point2D


def generate_2d_interaction_svg(
    smiles_or_mol: Any,
    hbonds: List[Dict[str, Any]],
    hydrophobics: List[Dict[str, Any]],
    salt_bridges: List[Dict[str, Any]],
    width: int = 700,
    height: int = 540,
    theme: str = "dark"
) -> str:
    """
    Generates a publication-grade 2D Protein-Ligand Interaction Diagram (PoseView / Maestro style).
    Renders central 2D ligand scaffold surrounded by color-coded residue badges,
    dashed hydrogen bond vectors with distance readouts, and hydrophobic contact indicators.
    Includes bi-directional interactive attributes for clicking residues to highlight in 3D viewport.
    """
    if isinstance(smiles_or_mol, str):
        mol = Chem.MolFromSmiles(smiles_or_mol)
    else:
        mol = smiles_or_mol

    if not mol:
        return ""

    try:
        mol = Chem.Mol(mol)
        AllChem.Compute2DCoords(mol)
    except Exception:
        return ""

    pad = 95
    lig_w = max(100, width - 2 * pad)
    lig_h = max(100, height - 2 * pad)

    is_dark = (theme.lower() == "dark")
    bg_color = "#0B0F19" if is_dark else "#FFFFFF"
    text_color = "#F8FAFC" if is_dark else "#0F172A"
    subtext_color = "#94A3B8" if is_dark else "#64748B"
    node_bg = "#1E293B" if is_dark else "#F8FAFC"
    badge_bg = "#0F172A" if is_dark else "#FFFFFF"
    border_color = "#334155" if is_dark else "#E2E8F0"

    drawer = rdMolDraw2D.MolDraw2DSVG(lig_w, lig_h)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.bondLineWidth = 2.4
    if is_dark:
        opts.useBWAtomPalette()

    try:
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        lig_svg = drawer.GetDrawingText()
    except Exception:
        return ""

    conf = mol.GetConformer()
    atom_pts = {}
    center_x, center_y = width / 2.0, height / 2.0
    xs, ys = [], []
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        pt = drawer.GetDrawCoords(Point2D(pos.x, pos.y))
        px, py = pt.x + pad, pt.y + pad
        atom_pts[i] = (px, py)
        xs.append(px)
        ys.append(py)

    lig_cx = sum(xs) / len(xs) if xs else center_x
    lig_cy = sum(ys) / len(ys) if ys else center_y

    # Collect distinct interacting residues with priority: H-Bond > Salt Bridge > Hydrophobic
    residues = {}
    for hb in hbonds:
        r_name = hb.get("receptor_residue") or f"{hb.get('resn', 'RES')} {hb.get('resi', 0)}"
        l_idx = hb.get("ligand_atom_idx", hb.get("ligand_atom_index", 0))
        resn = hb.get("resn", "RES")
        resi = hb.get("resi", 0)
        chain = hb.get("chain", "A")
        if resn in {"ASP", "GLU"}:
            color = "#EF4444"  # Acidic Red
        elif resn in {"ARG", "LYS", "HIS"}:
            color = "#3B82F6"  # Basic Blue
        elif resn in {"SER", "THR", "ASN", "GLN", "CYS", "TYR"}:
            color = "#10B981"  # Polar Green
        else:
            color = "#F59E0B"  # Hydrophobic Amber

        if r_name not in residues or residues[r_name]["type"] != "hbond":
            residues[r_name] = {
                "type": "hbond",
                "res_label": r_name,
                "resn": resn,
                "resi": resi,
                "chain": chain,
                "dist": hb.get("distance_angstrom", 2.8),
                "ligand_atom_idx": l_idx,
                "color": color
            }

    for sb in salt_bridges:
        r_name = sb.get("receptor_residue") or f"{sb.get('resn', 'RES')} {sb.get('resi', 0)}"
        if r_name not in residues:
            residues[r_name] = {
                "type": "salt_bridge",
                "res_label": r_name,
                "resn": sb.get("resn", "RES"),
                "resi": sb.get("resi", 0),
                "chain": sb.get("chain", "A"),
                "dist": sb.get("distance_angstrom", 3.2),
                "ligand_atom_idx": sb.get("ligand_atom_idx", sb.get("ligand_atom_index", 0)),
                "color": "#EC4899"
            }

    for hp in hydrophobics:
        r_name = hp.get("receptor_residue") or f"{hp.get('resn', 'RES')} {hp.get('resi', 0)}"
        if r_name not in residues and len(residues) < 16:
            residues[r_name] = {
                "type": "hydrophobic",
                "res_label": r_name,
                "resn": hp.get("resn", "RES"),
                "resi": hp.get("resi", 0),
                "chain": hp.get("chain", "A"),
                "dist": hp.get("distance_angstrom", 3.8),
                "ligand_atom_idx": hp.get("ligand_atom_idx", hp.get("ligand_atom_index", 0)),
                "color": "#F59E0B"
            }

    svg_elements = []
    angles_used = []

    for r_name, r_data in residues.items():
        l_idx = r_data["ligand_atom_idx"]
        if l_idx in atom_pts:
            lx, ly = atom_pts[l_idx]
        else:
            lx, ly = lig_cx, lig_cy

        dx = lx - lig_cx
        dy = ly - lig_cy
        angle = math.atan2(dy, dx)
        if abs(dx) < 1e-3 and abs(dy) < 1e-3:
            angle = len(angles_used) * (2 * math.pi / max(1, len(residues)))
        angles_used.append(angle)

        dist_out = pad * 0.82
        rx = lx + dist_out * math.cos(angle)
        ry = ly + dist_out * math.sin(angle)

        rx = max(50, min(width - 60, rx))
        ry = max(45, min(height - 40, ry))

        # 1. Connectors
        if r_data["type"] == "hbond":
            svg_elements.append(
                f'<line x1="{lx:.1f}" y1="{ly:.1f}" x2="{rx:.1f}" y2="{ry:.1f}" '
                f'stroke="{r_data["color"]}" stroke-width="2.5" stroke-dasharray="5,3" opacity="0.95" />'
            )
            mx = (lx + rx) / 2.0
            my = (ly + ry) / 2.0
            svg_elements.append(
                f'<rect x="{mx-19:.1f}" y="{my-9:.1f}" width="38" height="18" rx="5" '
                f'fill="{badge_bg}" stroke="{r_data["color"]}" stroke-width="1.2" />'
            )
            svg_elements.append(
                f'<text x="{mx:.1f}" y="{my+4:.1f}" fill="{text_color}" font-size="9.5" '
                f'font-family="monospace" font-weight="bold" text-anchor="middle">{r_data["dist"]}Å</text>'
            )
        elif r_data["type"] == "salt_bridge":
            svg_elements.append(
                f'<line x1="{lx:.1f}" y1="{ly:.1f}" x2="{rx:.1f}" y2="{ry:.1f}" '
                f'stroke="#EC4899" stroke-width="2.5" stroke-dasharray="6,3" opacity="0.95" />'
            )
            mx = (lx + rx) / 2.0
            my = (ly + ry) / 2.0
            svg_elements.append(
                f'<rect x="{mx-24:.1f}" y="{my-9:.1f}" width="48" height="18" rx="5" '
                f'fill="{badge_bg}" stroke="#EC4899" stroke-width="1.2" />'
            )
            svg_elements.append(
                f'<text x="{mx:.1f}" y="{my+4:.1f}" fill="#EC4899" font-size="9" '
                f'font-family="monospace" font-weight="bold" text-anchor="middle">± ionic</text>'
            )
        else:
            svg_elements.append(
                f'<line x1="{lx:.1f}" y1="{ly:.1f}" x2="{rx:.1f}" y2="{ry:.1f}" '
                f'stroke="#F59E0B" stroke-width="1.8" stroke-dasharray="2,3" opacity="0.75" />'
            )

        # 2. Residue Node Badge with hover & click handler
        badge_w, badge_h = 72, 28
        bx = rx - badge_w / 2.0
        by = ry - badge_h / 2.0
        svg_elements.append(
            f'<g class="residue-interactive cursor-pointer group" '
            f'onclick="highlightResidueIn3D(\'{r_data["resn"]}\', {r_data["resi"]}, \'{r_data["chain"]}\')" '
            f'data-resn="{r_data["resn"]}" data-resi="{r_data["resi"]}" data-chain="{r_data["chain"]}">'
            f'<rect x="{bx:.1f}" y="{by:.1f}" width="{badge_w}" height="{badge_h}" rx="8" '
            f'fill="{node_bg}" stroke="{r_data["color"]}" stroke-width="2" '
            f'class="transition duration-200 group-hover:stroke-white group-hover:brightness-125" />'
            f'<text x="{rx:.1f}" y="{ry+4.5:.1f}" fill="{r_data["color"]}" font-size="10.5" '
            f'font-family="monospace" font-weight="bold" text-anchor="middle" '
            f'class="select-none">{r_data["res_label"]}</text>'
            f'</g>'
        )

    # Extract inner SVG body from lig_svg
    svg_body_match = re.search(r'<svg[^>]*>(.*?)</svg>', lig_svg, re.DOTALL)
    svg_inner = svg_body_match.group(1) if svg_body_match else lig_svg

    # Invert black strokes to high-contrast white in dark mode
    if is_dark:
        svg_inner = svg_inner.replace('stroke:#000000', 'stroke:#E2E8F0').replace('stroke:black', 'stroke:#E2E8F0')

    final_svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="100%" height="100%" class="select-none font-sans">
    <rect width="{width}" height="{height}" fill="{bg_color}" rx="14" stroke="{border_color}" stroke-width="1" />
    
    <!-- Title and Subtitle -->
    <g transform="translate(18, 22)">
        <text x="0" y="0" fill="{text_color}" font-size="12" font-weight="bold" font-family="system-ui, sans-serif">2D Interaction Map (PoseView &#8226; Maestro Standard)</text>
        <text x="0" y="15" fill="{subtext_color}" font-size="9.5" font-family="system-ui, sans-serif">Click any amino acid badge to zoom &amp; inspect residue in 3D viewport</text>
    </g>

    <!-- Legend Badges (Top Right) -->
    <g transform="translate({width - 250}, 16)">
        <circle cx="5" cy="5" r="4" fill="#EF4444" /><text x="13" y="8" fill="{subtext_color}" font-size="9" font-family="system-ui, sans-serif">Acidic</text>
        <circle cx="55" cy="5" r="4" fill="#3B82F6" /><text x="63" y="8" fill="{subtext_color}" font-size="9" font-family="system-ui, sans-serif">Basic</text>
        <circle cx="102" cy="5" r="4" fill="#10B981" /><text x="110" y="8" fill="{subtext_color}" font-size="9" font-family="system-ui, sans-serif">Polar</text>
        <circle cx="145" cy="5" r="4" fill="#F59E0B" /><text x="153" y="8" fill="{subtext_color}" font-size="9" font-family="system-ui, sans-serif">Hydrophobic</text>
    </g>
    
    <!-- Connectors and Amino Acid Badges -->
    {''.join(svg_elements)}
    
    <!-- Central 2D Ligand Scaffold -->
    <g transform="translate({pad}, {pad})">
        {svg_inner}
    </g>
</svg>'''
    return final_svg


def generate_pymol_script(complex_filename: str = "docked_complex.pdb", title: str = "iLoveChem ChemDock Complex") -> str:
    """
    Generates a publication-grade PyMOL (.pml) macro script.
    Allows researchers to load the docked complex directly into PyMOL with
    active-site sticks, ribbon transparency, H-bond dashed vectors, and residue labels.
    """
    script = f"""# ==============================================================================
# iLoveChem ChemDock v3.2 - PyMOL Publication Macro Script
# Scientific Macromolecular Visualization Session
# ==============================================================================
# Instructions:
# 1. Save docked complex PDB as '{complex_filename}' in your working directory.
# 2. Open PyMOL and run: @docking_session.pml
# ==============================================================================

# 1. Global Viewport & Rendering Settings
reinitialize
set antialias, 2
set ray_shadows, 0
set ray_trace_mode, 1
set ray_trace_fog, 0
bg_color white

# 2. Load Complex
load {complex_filename}, docked_complex
hide everything, docked_complex

# 3. Receptor Macromolecule Representation
select receptor, docked_complex and polymer
show cartoon, receptor
color slate, receptor
set cartoon_transparency, 0.4, receptor

# 4. Docked Ligand Small Molecule
select ligand, docked_complex and resn LIG
show sticks, ligand
color green, ligand and elem C
set stick_radius, 0.28, ligand

# 5. Active Site Pocket Residues (within 5.0 A of ligand)
select pocket, byres (receptor within 5.0 of ligand)
show sticks, pocket
color gray70, pocket and elem C
set stick_radius, 0.20, pocket

# 6. Hydrogen Bond Polar Contacts
distance hbonds, ligand, pocket, 3.5, mode=2
set dash_color, yellow
set dash_width, 2.5
set dash_gap, 0.3
set dash_length, 0.2

# 7. Active Site Residue Labels
label pocket and name CA, "%s%s" % (resn, resi)
set label_size, 14
set label_font_id, 7
set label_color, black
set label_position, [0.0, 1.2, 1.0]

# 8. Focus Camera on Pocket
zoom pocket, buffer=3.0
orient ligand

print(">>> iLoveChem PyMOL Publication Session Ready! <<<")
"""
    return script
