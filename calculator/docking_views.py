import functools
import logging
import math
import re
import ssl
import urllib.error
import urllib.request

from django.shortcuts import render
from django.views.decorators.csrf import ensure_csrf_cookie
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status

from .docking_core import ChemDockEngine
from .interaction_diagram import generate_2d_interaction_svg, generate_pymol_script

logger = logging.getLogger(__name__)

# Modern PDB IDs are 4 characters: a digit 1-9 followed by three letters/digits (ASCII only).
PDB_ID_RE = re.compile(r"^[1-9][A-Za-z0-9]{3}$")

MAX_RECEPTOR_CHARS = 10 * 1024 * 1024
MAX_LIGAND_CHARS = 1 * 1024 * 1024
GRID_SIZE_RANGE = (5.0, 60.0)  # Angstroms


@ensure_csrf_cookie
def docking_studio_view(request):
    """
    Renders the dedicated Inbuilt Molecular Docking Studio UI.
    """
    return render(request, "calculator/docking.html")


MAX_RCSB_BYTES = 25 * 1024 * 1024

@functools.lru_cache(maxsize=2)
def _fetch_rcsb_cached(clean_id: str) -> str:
    """In-memory LRU cache for RCSB PDB structures. Prevents redundant remote downloads and caps buffer to 25 MB."""
    rcsb_url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    ctx = ssl.create_default_context()
    req = urllib.request.Request(rcsb_url, headers={'User-Agent': 'iLoveChem-Docking-Engine/1.0'})
    try:
        with urllib.request.urlopen(req, context=ctx, timeout=15) as resp:  # nosec: B310
            raw_bytes = resp.read(MAX_RCSB_BYTES + 1)
            if len(raw_bytes) > MAX_RCSB_BYTES:
                raise ValueError(f"PDB entry '{clean_id}' exceeds the 25 MB in-memory safety limit.")
    except urllib.error.HTTPError as exc:
        # urlopen raises on 4xx/5xx; it never returns a non-200 response object.
        if exc.code == 404:
            raise FileNotFoundError(clean_id) from exc
        raise
    return raw_bytes.decode('utf-8', errors='ignore')


@api_view(['GET'])
def api_fetch_rcsb_pdb(request, pdb_id: str):
    """
    Fetches official macromolecular crystal structure from RCSB Protein Data Bank.
    Auto-detects co-crystallized ligand and calculates active site pocket centroid.
    """
    clean_id = pdb_id.strip().upper()
    if not PDB_ID_RE.match(clean_id):
        return Response(
            {"success": False, "error": f"Invalid PDB ID format: '{clean_id}'. PDB IDs must be 4 alphanumeric characters (e.g. 1HSG, 6LU7, 3ERT)."},
            status=status.HTTP_400_BAD_REQUEST
        )

    try:
        pdb_content = _fetch_rcsb_cached(clean_id)
    except FileNotFoundError:
        return Response(
            {"success": False, "error": f"PDB ID '{clean_id}' not found on RCSB Protein Data Bank."},
            status=status.HTTP_404_NOT_FOUND
        )
    except Exception:
        logger.exception("RCSB fetch failed for %s", clean_id)
        return Response(
            {"success": False, "error": "Could not reach the RCSB Protein Data Bank. Try again in a moment."},
            status=status.HTTP_502_BAD_GATEWAY
        )

    # Sniff co-crystallized ligand and calculate pocket centroid
    pocket_info = ChemDockEngine.detect_co_crystallized_ligand_and_centroid(pdb_content)
    pocket_info["found_co_crystallized_ligand"] = pocket_info.get("has_co_crystal_ligand", False)
    pocket_info["co_crystallized_ligand_resname"] = pocket_info.get("ligand_resname")
    pocket_info["co_crystallized_atoms_count"] = pocket_info.get("atom_count", 0)
    pocket_info["centroid_x"] = pocket_info.get("center_x", 0.0)
    pocket_info["centroid_y"] = pocket_info.get("center_y", 0.0)
    pocket_info["centroid_z"] = pocket_info.get("center_z", 0.0)

    return Response({
        "success": True,
        "pdb_id": clean_id,
        "pdb_content": pdb_content,
        "pocket": pocket_info
    }, status=status.HTTP_200_OK)


@api_view(['POST'])
def api_run_docking(request):
    """
    Executes in-memory molecular docking using the native ChemDock engine.
    POST Payload:
      - receptor_pdb (string): PDB macromolecule content
      - ligand_input (string): SMILES, Molfile, or PDB string
      - ligand_format (string, optional): 'smiles' or 'pdb'
      - center_x, center_y, center_z (float): Pocket grid centroid
      - size_x, size_y, size_z (float, optional): Pocket box dimensions (default: 20 Å)
      - exhaustiveness (int, optional): Search depth (default: 6)
      - num_modes (int, optional): Number of ranked poses (default: 9)
    """
    data = request.data
    if not isinstance(data, dict):
        return Response(
            {"success": False, "error": "Invalid JSON object payload."},
            status=status.HTTP_400_BAD_REQUEST
        )
    receptor_pdb = str(data.get("receptor_pdb") or "").strip()
    ligand_input = str(data.get("ligand_input") or "").strip()
    ligand_format = str(data.get("ligand_format") or "smiles").strip().lower()

    if not receptor_pdb:
        return Response(
            {"success": False, "error": "Receptor macromolecule (PDB structure) is required."},
            status=status.HTTP_400_BAD_REQUEST
        )

    if not ligand_input:
        return Response(
            {"success": False, "error": "Ligand small molecule (SMILES or structure) is required."},
            status=status.HTTP_400_BAD_REQUEST
        )

    if len(receptor_pdb) > MAX_RECEPTOR_CHARS or len(ligand_input) > MAX_LIGAND_CHARS:
        return Response(
            {"success": False, "error": "Input is too large for the in-browser docking engine."},
            status=status.HTTP_413_REQUEST_ENTITY_TOO_LARGE
        )

    try:
        center_x = float(data.get("center_x", 0.0))
        center_y = float(data.get("center_y", 0.0))
        center_z = float(data.get("center_z", 0.0))
        size_x = float(data.get("size_x", 20.0))
        size_y = float(data.get("size_y", 20.0))
        size_z = float(data.get("size_z", 20.0))
        exhaustiveness = max(1, min(8, int(data.get("exhaustiveness", 8))))
        num_modes = max(1, min(15, int(data.get("num_modes", 9))))
    except (ValueError, TypeError):
        return Response(
            {"success": False, "error": "Grid centre and size must be numbers."},
            status=status.HTTP_400_BAD_REQUEST
        )

    if not all(math.isfinite(v) for v in (center_x, center_y, center_z, size_x, size_y, size_z)):
        return Response(
            {"success": False, "error": "Grid values must be finite numbers."},
            status=status.HTTP_400_BAD_REQUEST
        )
    lo, hi = GRID_SIZE_RANGE
    size_x, size_y, size_z = (min(max(v, lo), hi) for v in (size_x, size_y, size_z))

    # Step 1: Prepare rigid receptor PDBQT
    try:
        rec_pdbqt, rec_meta = ChemDockEngine.prepare_receptor_pdbqt(receptor_pdb, strip_co_crystal_ligands=True, pocket_center=(center_x, center_y, center_z), pocket_radius=max(size_x, size_y, size_z) + 12.0)
    except Exception:
        logger.exception("Receptor preparation failed")
        return Response(
            {"success": False, "error": "Receptor preparation failed. Check that the PDB file is valid."},
            status=status.HTTP_422_UNPROCESSABLE_ENTITY
        )

    # Step 2: Prepare flexible ligand PDBQT
    try:
        lig_pdbqt, lig_meta = ChemDockEngine.prepare_ligand_pdbqt(ligand_input, input_format=ligand_format)
    except Exception:
        logger.exception("Ligand preparation failed")
        return Response(
            {"success": False, "error": "Ligand preparation failed. Check the SMILES or structure."},
            status=status.HTTP_422_UNPROCESSABLE_ENTITY
        )

    # Step 3: Run ChemDock Native Monte Carlo & Empirical Scoring
    try:
        dock_results = ChemDockEngine.run_docking(
            receptor_pdbqt=rec_pdbqt,
            ligand_pdbqt=lig_pdbqt,
            center_x=center_x,
            center_y=center_y,
            center_z=center_z,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z,
            exhaustiveness=exhaustiveness,
            num_modes=num_modes
        )
    except Exception:
        logger.exception("Docking run failed")
        return Response(
            {"success": False, "error": "The docking calculation failed."},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR
        )

    # Step 4: Create consolidated Receptor-Ligand PDB complex for top pose
    consolidated_complex_pdb = ""
    if dock_results.get("poses"):
        top_pose_pdbqt = dock_results["poses"][0]["pdbqt_block"]
        consolidated_complex_pdb = _merge_receptor_and_pose_to_pdb(receptor_pdb, top_pose_pdbqt)

    lig_smiles = lig_meta.get("smiles", "") or ligand_input
    poses = dock_results.get("poses", [])

    # Wire Ligand Efficiency (LE & LLE) into docking poses
    from calculator.medchem_engine import analyze_medchem_profile, compute_ligand_efficiency
    lig_profile = analyze_medchem_profile(lig_smiles) if lig_smiles else None
    heavy_atoms = lig_profile.get("heavy_atoms", 0) if lig_profile else 0
    lig_logp = lig_profile.get("logp", None) if lig_profile else None

    for idx, p in enumerate(poses):
        aff = p.get("affinity")
        if isinstance(aff, (int, float)) and heavy_atoms > 0:
            p["efficiency"] = compute_ligand_efficiency(aff, heavy_atoms, lig_logp)
        p["rmsd_lower_bound"] = p.get("rmsd_lb", 0.0)
        p["rmsd_upper_bound"] = p.get("rmsd_ub", 0.0)
        p_hb = p.get("hbond_contacts", [])
        p_hp = p.get("hydrophobic_contacts", [])
        p_sb = p.get("salt_bridges", [])
        p["contacts"] = p_hb
        p["num_h_bonds"] = len(p_hb)
        # Generate dynamic 2D interaction diagram for each pose
        p["interaction_diagram_svg"] = generate_2d_interaction_svg(
            lig_smiles, p_hb, p_hp, p_sb, width=700, height=540, theme="dark"
        )
        p["interaction_diagram_svg_white"] = generate_2d_interaction_svg(
            lig_smiles, p_hb, p_hp, p_sb, width=700, height=540, theme="white"
        )

    # Deliverables for top pose
    top_hb = dock_results.get("hbond_contacts", [])
    top_hp = dock_results.get("hydrophobic_contacts", [])
    top_sb = dock_results.get("salt_bridges", [])
    top_svg_dark = generate_2d_interaction_svg(lig_smiles, top_hb, top_hp, top_sb, width=700, height=540, theme="dark")
    top_svg_white = generate_2d_interaction_svg(lig_smiles, top_hb, top_hp, top_sb, width=700, height=540, theme="white")
    pymol_script = generate_pymol_script("docked_complex.pdb")

    if poses and "efficiency" in poses[0]:
        dock_results["top_efficiency"] = poses[0]["efficiency"]

    dock_results["num_poses"] = dock_results.get("poses_count", len(poses))
    dock_results["exhaustiveness"] = exhaustiveness
    dock_results["rotatable_bonds"] = lig_meta.get("rotatable_bonds", 0)
    dock_results["receptor_meta"] = rec_meta
    dock_results["ligand_meta"] = lig_meta
    dock_results["consolidated_complex_pdb"] = consolidated_complex_pdb
    dock_results["interaction_diagram_svg"] = top_svg_dark
    dock_results["interaction_diagram_svg_white"] = top_svg_white
    dock_results["pymol_macro_script"] = pymol_script

    return Response(dock_results, status=status.HTTP_200_OK)


def _merge_receptor_and_pose_to_pdb(receptor_pdb: str, pose_pdbqt: str) -> str:
    """Merges receptor protein lines with docked ligand coordinates into a standard PDB complex."""
    out = []
    max_serial = 0
    for line in (receptor_pdb or "").splitlines():
        if line.startswith(("ATOM  ", "HETATM", "CRYST1", "TER")):
            out.append(line.rstrip())
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 11:
                try:
                    max_serial = max(max_serial, int(line[6:11].strip()))
                except ValueError:
                    pass
    out.append("TER")
    atom_idx = max(max_serial + 1, 1)
    for line in (pose_pdbqt or "").splitlines():
        if (line.startswith("ATOM  ") or line.startswith("HETATM")) and len(line) >= 54:
            serial_str = f"{min(atom_idx, 99999):5d}"
            name = line[12:16].strip()
            x = line[30:38]
            y = line[38:46]
            z = line[46:54]
            element = line[77:79].strip() if len(line) >= 79 else ""
            if not element:
                element = name[:1]
            out.append(f"HETATM{serial_str} {name:<4s} LIG Z 900    {x}{y}{z}  1.00 20.00          {element:>2s}")
            atom_idx += 1
    out.append("END")
    return "\n".join(out) + "\n"
