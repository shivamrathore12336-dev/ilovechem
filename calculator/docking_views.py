import os
import json
import urllib.request
import ssl
from typing import Dict, Any

from django.shortcuts import render
from django.http import JsonResponse, HttpResponse
from django.views.decorators.csrf import ensure_csrf_cookie
from rest_framework.decorators import api_view
from rest_framework.response import Response
from rest_framework import status

from .docking_core import ChemDockEngine


@ensure_csrf_cookie
def docking_studio_view(request):
    """
    Renders the dedicated Inbuilt Molecular Docking Studio UI.
    """
    return render(request, "calculator/docking.html")


@api_view(['GET'])
def api_fetch_rcsb_pdb(request, pdb_id: str):
    """
    Fetches official macromolecular crystal structure from RCSB Protein Data Bank.
    Auto-detects co-crystallized ligand and calculates active site pocket centroid.
    """
    clean_id = pdb_id.strip().upper()
    if not (clean_id.isalnum() and len(clean_id) == 4):
        return Response(
            {"success": False, "error": f"Invalid PDB ID format: '{clean_id}'. PDB IDs must be 4 alphanumeric characters (e.g. 1HSG, 6LU7, 3ERT)."},
            status=status.HTTP_400_BAD_REQUEST
        )

    rcsb_url = f"https://files.rcsb.org/download/{clean_id}.pdb"
    ctx = ssl.create_default_context()

    try:
        req = urllib.request.Request(rcsb_url, headers={'User-Agent': 'iLoveChem-Docking-Engine/1.0'})
        with urllib.request.urlopen(req, context=ctx, timeout=15) as resp:
            if resp.status != 200:
                return Response(
                    {"success": False, "error": f"PDB ID '{clean_id}' not found on RCSB Protein Data Bank."},
                    status=status.HTTP_404_NOT_FOUND
                )
            pdb_content = resp.read().decode('utf-8', errors='ignore')
    except Exception as e:
        return Response(
            {"success": False, "error": f"Failed to retrieve PDB '{clean_id}' from RCSB: {str(e)}"},
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
    receptor_pdb = data.get("receptor_pdb", "").strip()
    ligand_input = data.get("ligand_input", "").strip()
    ligand_format = data.get("ligand_format", "smiles").strip().lower()

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

    try:
        center_x = float(data.get("center_x", 0.0))
        center_y = float(data.get("center_y", 0.0))
        center_z = float(data.get("center_z", 0.0))
        size_x = float(data.get("size_x", 20.0))
        size_y = float(data.get("size_y", 20.0))
        size_z = float(data.get("size_z", 20.0))
        exhaustiveness = max(1, min(64, int(data.get("exhaustiveness", 8))))
        num_modes = max(1, min(15, int(data.get("num_modes", 9))))
    except (ValueError, TypeError) as e:
        return Response(
            {"success": False, "error": f"Invalid numerical grid coordinate or dimension: {str(e)}"},
            status=status.HTTP_400_BAD_REQUEST
        )

    # Step 1: Prepare rigid receptor PDBQT
    try:
        rec_pdbqt, rec_meta = ChemDockEngine.prepare_receptor_pdbqt(receptor_pdb, strip_co_crystal_ligands=True, pocket_center=(center_x, center_y, center_z), pocket_radius=max(size_x, size_y, size_z) + 12.0)
    except Exception as e:
        return Response(
            {"success": False, "error": f"Receptor preparation failed: {str(e)}"},
            status=status.HTTP_422_UNPROCESSABLE_ENTITY
        )

    # Step 2: Prepare flexible ligand PDBQT
    try:
        lig_pdbqt, lig_meta = ChemDockEngine.prepare_ligand_pdbqt(ligand_input, input_format=ligand_format)
    except Exception as e:
        return Response(
            {"success": False, "error": f"Ligand 3D conformation & PDBQT preparation failed: {str(e)}"},
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
    except Exception as e:
        return Response(
            {"success": False, "error": f"ChemDock docking calculation failed: {str(e)}"},
            status=status.HTTP_500_INTERNAL_SERVER_ERROR
        )

    # Step 4: Create consolidated Receptor-Ligand PDB complex for top pose
    consolidated_complex_pdb = ""
    if dock_results.get("poses"):
        top_pose_pdbqt = dock_results["poses"][0]["pdbqt_block"]
        consolidated_complex_pdb = _merge_receptor_and_pose_to_pdb(receptor_pdb, top_pose_pdbqt)

    hbond_contacts = dock_results.get("hbond_contacts", [])
    poses = dock_results.get("poses", [])
    for idx, p in enumerate(poses):
        p["rmsd_lower_bound"] = p.get("rmsd_lb", 0.0)
        p["rmsd_upper_bound"] = p.get("rmsd_ub", 0.0)
        if idx == 0:
            p["contacts"] = hbond_contacts
            p["num_h_bonds"] = len(hbond_contacts)
        else:
            p["contacts"] = []
            p["num_h_bonds"] = 0

    dock_results["num_poses"] = dock_results.get("poses_count", len(poses))
    dock_results["receptor_meta"] = rec_meta
    dock_results["ligand_meta"] = lig_meta
    dock_results["consolidated_complex_pdb"] = consolidated_complex_pdb

    return Response(dock_results, status=status.HTTP_200_OK)


def _merge_receptor_and_pose_to_pdb(receptor_pdb: str, pose_pdbqt: str) -> str:
    """Merges receptor protein lines with docked ligand coordinates into a standard PDB complex."""
    out_lines = []
    for line in receptor_pdb.splitlines():
        rec = line[0:6].strip()
        if rec in ["CRYST1", "ATOM"]:
            out_lines.append(line)

    out_lines.append("TER")

    atom_idx = 1000
    for line in pose_pdbqt.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atom_idx += 1
            # Format as HETATM LIG X 900
            name = line[12:16].strip()
            x = line[30:38]
            y = line[38:46]
            z = line[46:54]
            element = line[77:79].strip() or name[:1]
            pdb_line = f"HETATM{atom_idx:5d} {name:<4s} LIG Z 900    {x}{y}{z}  1.00 20.00          {element:>2s}"
            out_lines.append(pdb_line)

    out_lines.append("END")
    return "\n".join(out_lines) + "\n"
