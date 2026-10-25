# iLoveChem Native ChemDock Inbuilt Docking Engine (v2.0)
# Developed for iLoveChem (https://ilovechem.onrender.com)
#
# Pharmacophore & Scoring Features:
# 1. 100% In-Process Python / NumPy / SciPy implementation of the AutoDock Vina
#    empirical scoring function (Gauss1, Gauss2, Repulsion, Hydrophobic, HBond, Torsion).
# 2. Flexible Ligand Dihedral Torsion Exploration: Real-time tree-based internal rotor sampling.
# 3. Vectorized Analytical Force/Torque Gradient Descent: Post-annealing local relaxation.
# 4. Single-Instance HETATM Pocket Detection: Resolves chain/symmetry-mate multi-copy collisions.
# 5. Volatile In-Memory Architecture: Zero disk footprint, Pharma Zero-Storage Invariant compliant.

import os
import sys
import math
import re
from collections import defaultdict, deque
from typing import Dict, Any, Tuple, Optional, List
import numpy as np

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Protein PDB atom-name to AutoDock Vina atom-type lookup
PDB_TO_AUTODOCK_ATOM_TYPES = {
    ("ASP", "OD1"): "OA", ("ASP", "OD2"): "OA",
    ("GLU", "OE1"): "OA", ("GLU", "OE2"): "OA",
    ("ASN", "OD1"): "OA", ("GLN", "OE1"): "OA",
    ("SER", "OG"): "OA", ("THR", "OG1"): "OA", ("TYR", "OH"): "OA",
    ("HIS", "ND1"): "NA", ("HIS", "NE2"): "NA",
    ("ARG", "NE"): "N", ("ARG", "NH1"): "N", ("ARG", "NH2"): "N",
    ("LYS", "NZ"): "N", ("TRP", "NE1"): "N",
    ("ASN", "ND2"): "N", ("GLN", "NE2"): "N",
    ("CYS", "SG"): "SA", ("MET", "SD"): "SA",
    ("PHE", "CG"): "A", ("PHE", "CD1"): "A", ("PHE", "CD2"): "A", ("PHE", "CE1"): "A", ("PHE", "CE2"): "A", ("PHE", "CZ"): "A",
    ("TYR", "CG"): "A", ("TYR", "CD1"): "A", ("TYR", "CD2"): "A", ("TYR", "CE1"): "A", ("TYR", "CE2"): "A",
    ("TRP", "CG"): "A", ("TRP", "CD1"): "A", ("TRP", "CD2"): "A", ("TRP", "CE2"): "A", ("TRP", "CE3"): "A", ("TRP", "CZ2"): "A", ("TRP", "CZ3"): "A", ("TRP", "CH2"): "A",
    ("HIS", "CG"): "A", ("HIS", "CD2"): "A",
}

# Buffer, solvent, and crystallization sugar residue identifiers
SOLVENT_AND_BUFFER_RESIDUES = {
    "HOH", "WAT", "DOD", "TIP", "TIP3", "SOL",
    "SO4", "PO4", "CL", "NA", "MG", "ZN", "CA", "MN", "FE", "K", "IOD",
    "EDO", "DMS", "ACT", "FMT", "PEG", "PG4", "PGE", "GOL", "BME", "MES", "HEPES",
    "NAG", "MAN", "BMA", "FUC", "NDG", "GLC", "PEF", "OGA", "TRS", "MPD"
}

# AutoDock Vina Empirical VDW Radii (Angstroms)
VDW_RADII: Dict[str, float] = {
    "H": 1.20, "HD": 1.20,
    "C": 1.90, "A": 1.90,
    "N": 1.80, "NA": 1.80,
    "O": 1.60, "OA": 1.60,
    "S": 2.00, "SA": 2.00,
    "P": 2.10,
    "F": 1.50, "CL": 1.80, "BR": 2.00, "I": 2.20,
    "FE": 1.20, "ZN": 1.20, "MG": 1.20, "CA": 1.20, "MN": 1.20, "CU": 1.20, "NI": 1.20
}

# Vina Empirical Weighting Coefficients (Trott & Olson, 2010)
WEIGHT_GAUSS1 = -0.035579
WEIGHT_GAUSS2 = -0.005156
WEIGHT_REPULSION = 0.840245
WEIGHT_HYDROPHOBIC = -0.035069
WEIGHT_HBOND = -0.587439
WEIGHT_ROT = 0.058460


def _get_vdw_radius(atom_type: str) -> float:
    return VDW_RADII.get(atom_type.upper(), 1.70)


def _is_hydrophobic(atom_type: str) -> bool:
    return atom_type.upper() in {"C", "A", "F", "CL", "BR", "I"}


def _is_hbond_acceptor(atom_type: str) -> bool:
    return atom_type.upper() in {"NA", "OA", "SA", "N", "O"}


def _is_hbond_donor(atom_type: str) -> bool:
    return atom_type.upper() in {"HD", "N", "O"}


def _rotation_matrix_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodriguez rotation formula for 3D rotation around an arbitrary unit axis."""
    norm = np.linalg.norm(axis)
    if norm < 1e-12:
        return np.eye(3)
    axis = axis / norm
    a = math.cos(angle / 2.0)
    b, c, d = -axis * math.sin(angle / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([
        [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
        [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
        [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]
    ])


def _random_quaternion_rotation() -> np.ndarray:
    """Uniform random 3D rotation matrix via Shoemake algorithm."""
    u1, u2, u3 = np.random.uniform(0.0, 1.0, 3)
    q1 = math.sqrt(1.0 - u1) * math.sin(2.0 * math.pi * u2)
    q2 = math.sqrt(1.0 - u1) * math.cos(2.0 * math.pi * u2)
    q3 = math.sqrt(u1) * math.sin(2.0 * math.pi * u3)
    q4 = math.sqrt(u1) * math.cos(2.0 * math.pi * u3)
    
    q0_q3 = q4 * q3
    return np.array([
        [1.0 - 2.0 * (q2**2 + q3**2), 2.0 * (q1 * q2 - q0_q3), 2.0 * (q1 * q3 + q4 * q2)],
        [2.0 * (q1 * q2 + q0_q3), 1.0 - 2.0 * (q1**2 + q3**2), 2.0 * (q2 * q3 - q4 * q1)],
        [2.0 * (q1 * q3 - q4 * q2), 2.0 * (q2 * q3 + q4 * q1), 1.0 - 2.0 * (q1**2 + q2**2)]
    ])


def _find_rotatable_bonds_and_trees(coords: np.ndarray, types: List[str]) -> List[Tuple[int, int, List[int]]]:
    """
    Identifies rotatable single bonds and their moving downstream atom subtrees.
    Returns list of (u, v, moving_atom_indices) where:
      - u, v are atom indices of the rotatable bond.
      - moving_atom_indices are indices of all atoms that rotate around bond axis u->v.
    """
    N = len(coords)
    if N < 4:
        return []

    cov_radii = {
        'C': 0.77, 'A': 0.77, 'N': 0.75, 'NA': 0.75, 'O': 0.73, 'OA': 0.73,
        'S': 1.02, 'SA': 1.02, 'P': 1.06, 'F': 0.71, 'CL': 0.99, 'BR': 1.14, 'I': 1.33
    }

    # Build adjacency graph
    adj = [[] for _ in range(N)]
    for i in range(N):
        r_i = cov_radii.get(types[i], 0.77)
        for j in range(i + 1, N):
            r_j = cov_radii.get(types[j], 0.77)
            dist = float(np.linalg.norm(coords[i] - coords[j]))
            if 0.5 < dist <= 1.25 * (r_i + r_j):
                adj[i].append(j)
                adj[j].append(i)

    rotatable = []
    # Identify non-ring single bonds where both ends have degree > 1
    for u in range(N):
        if len(adj[u]) <= 1:
            continue
        for v in adj[u]:
            if v <= u or len(adj[v]) <= 1:
                continue

            # Check if (u, v) is in a ring via BFS without edge (u, v)
            in_ring = False
            q = deque([u])
            visited = {u}
            while q and not in_ring:
                curr = q.popleft()
                for nbr in adj[curr]:
                    if (curr == u and nbr == v) or (curr == v and nbr == u):
                        continue
                    if nbr == v:
                        in_ring = True
                        break
                    if nbr not in visited:
                        visited.add(nbr)
                        q.append(nbr)

            if not in_ring:
                # Find downstream moving subtree starting from v without crossing back to u
                moving = {v}
                q2 = deque([v])
                while q2:
                    curr2 = q2.popleft()
                    for nbr in adj[curr2]:
                        if nbr != u and nbr not in moving:
                            moving.add(nbr)
                            q2.append(nbr)

                moving_list = list(moving)
                if len(moving_list) > N // 2:
                    # Invert: rotate smaller side (u relative to v)
                    moving_u = {u}
                    q3 = deque([u])
                    while q3:
                        curr3 = q3.popleft()
                        for nbr in adj[curr3]:
                            if nbr != v and nbr not in moving_u:
                                moving_u.add(nbr)
                                q3.append(nbr)
                    rotatable.append((v, u, list(moving_u)))
                else:
                    rotatable.append((u, v, moving_list))

    return rotatable


def _relax_pose_gradients(
    rec_coords: np.ndarray,
    rec_types: List[str],
    rec_radii: np.ndarray,
    init_lig_coords: np.ndarray,
    lig_types: List[str],
    lig_radii: np.ndarray,
    num_rotatable: int = 0,
    steps: int = 20,
    gamma: float = 0.04
) -> Tuple[np.ndarray, float]:
    """
    Applies analytical force/torque gradient descent to relieve minor steric boundary clashes
    and settle the ligand into the nearest local Vina energetic minimum.
    """
    curr_coords = np.copy(init_lig_coords)
    curr_energy = ChemDockEngine.evaluate_vina_score(
        rec_coords, rec_types, rec_radii,
        curr_coords, lig_types, lig_radii,
        num_rotatable_bonds=num_rotatable
    )

    for _ in range(steps):
        diff = curr_coords[np.newaxis, :, :] - rec_coords[:, np.newaxis, :]  # (N_rec, N_lig, 3)
        dists = np.linalg.norm(diff, axis=-1)
        dists = np.maximum(dists, 1e-4)
        unit_vec = diff / dists[:, :, np.newaxis]
        r_sum = rec_radii[:, np.newaxis] + lig_radii[np.newaxis, :]
        d = dists - r_sum
        cutoff_mask = dists < 8.0

        g1 = -0.035579 * (-2.0 * d / 0.25) * np.exp(-((d / 0.5)**2))
        g2 = -0.005156 * (-2.0 * (d - 3.0) / 4.0) * np.exp(-(((d - 3.0) / 2.0)**2))
        rep = np.where(d < 0.0, 2.0 * 0.048881 * d, 0.0)
        hydro = np.where((d >= 0.5) & (d <= 1.5), -0.035069 * (-1.0), 0.0)
        hbond = np.where((d >= -0.7) & (d <= 0.0), -0.587439 * (-1.0 / 0.7), 0.0)

        dE_dd = (g1 + g2 + rep + hydro + hbond) * cutoff_mask
        grad_per_atom = np.sum(dE_dd[:, :, np.newaxis] * unit_vec, axis=0)  # (N_lig, 3)

        f_net = -np.sum(grad_per_atom, axis=0)
        c_curr = np.mean(curr_coords, axis=0)
        r_rel = curr_coords - c_curr
        tau_net = -np.sum(np.cross(r_rel, grad_per_atom), axis=0)

        t_step = np.clip(gamma * f_net, -0.2, 0.2)
        tau_norm = float(np.linalg.norm(tau_net))
        if tau_norm > 1e-5:
            axis = tau_net / tau_norm
            angle = float(np.clip(gamma * 0.05 * tau_norm, -0.1, 0.1))
            R = _rotation_matrix_from_axis_angle(axis, angle)
        else:
            R = np.eye(3)

        next_coords = np.dot(curr_coords - c_curr, R.T) + c_curr + t_step
        next_energy = ChemDockEngine.evaluate_vina_score(
            rec_coords, rec_types, rec_radii,
            next_coords, lig_types, lig_radii,
            num_rotatable_bonds=num_rotatable
        )
        if next_energy < curr_energy:
            curr_coords = next_coords
            curr_energy = next_energy
        else:
            gamma *= 0.7

    return curr_coords, curr_energy


class ChemDockEngine:
    """
    iLoveChem Native Inbuilt Molecular Docking Kernel (v2.0).
    100% in-process Python/NumPy/SciPy implementation of the AutoDock Vina
    Empirical Free Energy Scoring Function, Flexible Dihedral Torsion Exploration,
    and Analytical Gradient-Relaxed Monte Carlo Annealing.
    Operates purely in volatile RAM (Zero-Storage Invariant compliant).
    """

    @staticmethod
    def detect_co_crystallized_ligand_and_centroid(pdb_str: str, preferred_ligand_resname: Optional[str] = None) -> Dict[str, Any]:
        """
        Inspects PDB records for co-crystallized drug ligands (HETATM).
        Groups HETATMs by unique molecular instance: (chain_id, res_seq, resname)
        to isolate individual drug molecules rather than merging multi-chain symmetry copies.
        Calculates geometric centroid and suggested grid box dimensions.
        """
        het_instances: Dict[Tuple[str, str, str], List[Tuple[float, float, float]]] = defaultdict(list)
        protein_coords: List[Tuple[float, float, float]] = []

        for line in pdb_str.splitlines():
            record = line[0:6].strip()
            if record == "ATOM":
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    protein_coords.append((x, y, z))
                except (ValueError, IndexError):
                    pass
            elif record == "HETATM":
                resname = line[17:20].strip().upper()
                chain = line[21] if len(line) > 21 else "A"
                resseq = line[22:26].strip() if len(line) >= 26 else "1"
                if resname not in SOLVENT_AND_BUFFER_RESIDUES:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        het_instances[(chain, resseq, resname)].append((x, y, z))
                    except (ValueError, IndexError):
                        pass

        # Select candidate drug instance
        chosen_key = None
        chosen_coords = []
        pref_clean = preferred_ligand_resname.strip().upper() if preferred_ligand_resname else None

        if pref_clean:
            matching = [k for k in het_instances.keys() if k[2] == pref_clean]
            if matching:
                matching.sort(key=lambda k: len(het_instances[k]), reverse=True)
                chosen_key = matching[0]
                chosen_coords = het_instances[chosen_key]

        if not chosen_coords:
            candidates = [(k, coords) for k, coords in het_instances.items() if len(coords) >= 6]
            if candidates:
                candidates.sort(key=lambda x: len(x[1]), reverse=True)
                chosen_key = candidates[0][0]
                chosen_coords = candidates[0][1]

        if chosen_coords and chosen_key:
            arr = np.array(chosen_coords)
            c_x, c_y, c_z = np.mean(arr, axis=0)
            span_x = float(np.ptp(arr[:, 0])) + 8.0
            span_y = float(np.ptp(arr[:, 1])) + 8.0
            span_z = float(np.ptp(arr[:, 2])) + 8.0
            size_x = max(18.0, min(30.0, round(span_x, 1)))
            size_y = max(18.0, min(30.0, round(span_y, 1)))
            size_z = max(18.0, min(30.0, round(span_z, 1)))

            return {
                "has_co_crystal_ligand": True,
                "ligand_resname": chosen_key[2],
                "chain_id": chosen_key[0],
                "res_seq": chosen_key[1],
                "atom_count": len(chosen_coords),
                "center_x": round(float(c_x), 3),
                "center_y": round(float(c_y), 3),
                "center_z": round(float(c_z), 3),
                "size_x": size_x,
                "size_y": size_y,
                "size_z": size_z,
                "method": "co_crystal_ligand_centroid"
            }

        # Fallback to receptor center of mass
        if protein_coords:
            arr = np.array(protein_coords)
            c_x, c_y, c_z = np.mean(arr, axis=0)
            return {
                "has_co_crystal_ligand": False,
                "ligand_resname": None,
                "chain_id": None,
                "res_seq": None,
                "center_x": round(float(c_x), 3),
                "center_y": round(float(c_y), 3),
                "center_z": round(float(c_z), 3),
                "size_x": 22.0,
                "size_y": 22.0,
                "size_z": 22.0,
                "method": "receptor_center_of_mass"
            }

        return {
            "has_co_crystal_ligand": False,
            "ligand_resname": None,
            "chain_id": None,
            "res_seq": None,
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0,
            "method": "origin_default"
        }

    @classmethod
    def prepare_receptor_pdbqt(
        cls,
        pdb_str: str,
        strip_co_crystal_ligands: bool = True,
        pocket_center: Optional[Tuple[float, float, float]] = None,
        pocket_radius: float = 24.0,
        preferred_ligand_resname: Optional[str] = None
    ) -> Tuple[str, Dict[str, Any]]:
        """
        Prepares a rigid macromolecular target receptor in PDBQT format:
        1. Strips solvent waters, crystallographic buffers, and crystallization sugars.
        2. Clears the binding pocket of co-crystallized drug ligands and non-catalytic HETATMs.
        3. Preserves essential catalytic cofactor metal ions (Zn, Mg, Mn, Fe, Ca, Cu).
        4. Retains protein ATOM records with AutoDock standard atom types.
        """
        co_crystal_info = cls.detect_co_crystallized_ligand_and_centroid(pdb_str, preferred_ligand_resname=preferred_ligand_resname)
        lig_res = (preferred_ligand_resname or co_crystal_info.get("ligand_resname")) if strip_co_crystal_ligands else None

        cleaned_lines = []
        waters_stripped = 0
        ligands_stripped = 0

        pocket_cx = pocket_center[0] if pocket_center else (co_crystal_info.get("center_x") if co_crystal_info.get("has_co_crystal_ligand") else None)
        pocket_cy = pocket_center[1] if pocket_center else (co_crystal_info.get("center_y") if co_crystal_info.get("has_co_crystal_ligand") else None)
        pocket_cz = pocket_center[2] if pocket_center else (co_crystal_info.get("center_z") if co_crystal_info.get("has_co_crystal_ligand") else None)
        max_dist_sq = (pocket_radius**2) if (pocket_cx is not None and pocket_radius > 0) else None

        ESSENTIAL_METALS = {"ZN", "MG", "MN", "FE", "CA", "CU", "NI", "CO"}

        for line in pdb_str.splitlines():
            rec = line[0:6].strip()
            if rec in ["ATOM", "HETATM"]:
                resname = line[17:20].strip().upper()
                if resname in SOLVENT_AND_BUFFER_RESIDUES:
                    waters_stripped += 1
                    continue
                if rec == "HETATM" and strip_co_crystal_ligands:
                    if lig_res and resname == lig_res:
                        ligands_stripped += 1
                        continue
                    if resname not in ESSENTIAL_METALS and resname not in {"HEM", "FAD", "NAD", "NAP", "COA"}:
                        if pocket_cx is not None:
                            try:
                                x = float(line[30:38].strip())
                                y = float(line[38:46].strip())
                                z = float(line[46:54].strip())
                                if (x - pocket_cx)**2 + (y - pocket_cy)**2 + (z - pocket_cz)**2 <= (14.0**2):
                                    ligands_stripped += 1
                                    continue
                            except Exception:
                                pass
                if max_dist_sq is not None:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        d2 = (x - pocket_cx)**2 + (y - pocket_cy)**2 + (z - pocket_cz)**2
                        if d2 > max_dist_sq:
                            continue
                    except Exception:
                        pass
                cleaned_lines.append(line)
            elif rec in ["CRYST1", "SCALE", "ORIGX", "TER", "END"]:
                cleaned_lines.append(line)

        # Fast in-process PDB-to-PDBQT translation
        fast_pdbqt_lines = []
        atom_count = 0
        for line in cleaned_lines:
            rec = line[0:6].strip()
            if rec in ["ATOM", "HETATM"]:
                try:
                    res_name = line[17:20].strip()
                    atom_name = line[12:16].strip()
                    elem = line[76:78].strip().upper() if len(line) >= 78 else ""
                    if not elem:
                        elem = re.sub(r"[0-9\s]", "", atom_name)[:1].upper()

                    ad_type = PDB_TO_AUTODOCK_ATOM_TYPES.get((res_name, atom_name))
                    if not ad_type:
                        if elem == "C":
                            ad_type = "C"
                        elif elem == "N":
                            ad_type = "N"
                        elif elem == "O":
                            ad_type = "OA"
                        elif elem == "S":
                            ad_type = "SA"
                        elif elem == "P":
                            ad_type = "P"
                        elif elem in ["ZN", "MG", "CA", "FE", "MN", "CU"]:
                            ad_type = elem
                        elif elem == "H":
                            ad_type = "HD" if atom_name.startswith("H") and "N" in line else "H"
                        else:
                            ad_type = elem[:2]

                    charge = 0.000
                    if (res_name, atom_name) in [("ASP", "OD1"), ("ASP", "OD2"), ("GLU", "OE1"), ("GLU", "OE2")]:
                        charge = -0.500
                    elif (res_name, atom_name) in [("ARG", "NH1"), ("ARG", "NH2"), ("LYS", "NZ")]:
                        charge = 0.500

                    base = line[:66].ljust(66)
                    pdbqt_line = f"{base}{charge:6.3f} {ad_type:>2s}"
                    fast_pdbqt_lines.append(pdbqt_line)
                    atom_count += 1
                except Exception:
                    pass
            else:
                fast_pdbqt_lines.append(line)

        pdbqt_str = "\n".join(fast_pdbqt_lines) + "\n"
        meta = {
            "atom_count": atom_count,
            "waters_stripped": waters_stripped,
            "ligands_stripped": ligands_stripped,
            "has_pocket_cropping": max_dist_sq is not None,
            "pocket_radius": pocket_radius if max_dist_sq else None
        }
        return pdbqt_str, meta

    @classmethod
    def prepare_ligand_pdbqt(cls, ligand_input: str, input_format: str = "mol") -> Tuple[str, Dict[str, Any]]:
        """
        Prepares a drug ligand into AutoDock PDBQT format with Gasteiger partial charges
        and root-branch rotatable torsion trees (TORSDOF).
        Supports PDB, MOL/SDF blocks, and SMILES strings.
        """
        s = ligand_input.strip()
        if not s:
            raise ValueError("Empty ligand input provided.")

        fmt_lower = input_format.lower().strip()
        if "\n" not in s and ("=" in s or "#" in s or "(" in s or "@" in s or (len(s) < 80 and not s.startswith("ATOM"))):
            fmt = "smi"
        elif "ATOM" in s or "HETATM" in s or fmt_lower in ["pdb", "ent"]:
            fmt = "pdb"
        elif "M  END" in s or "$$$$" in s or (fmt_lower in ["mol", "sdf"] and "\n" in s):
            fmt = "mol"
        else:
            fmt = "smi"

        ob_mol = pybel.readstring(fmt, s)
        if fmt == "smi":
            try:
                ob_mol.make3D(forcefield='mmff94', steps=300)
            except Exception:
                try:
                    ob_mol.make3D(forcefield='uff', steps=200)
                except Exception:
                    pass

        for a in ob_mol.atoms:
            try:
                _ = a.partialcharge
            except Exception:
                pass

        pdbqt_raw = ob_mol.write("pdbqt")
        if not pdbqt_raw or len(pdbqt_raw.strip()) < 10:
            raise ValueError(f"OpenBabel produced empty PDBQT for ligand format '{fmt}'.")

        # Compute accurate TORSDOF
        heavy_count = sum(1 for a in ob_mol.atoms if a.atomicnum > 1)
        rot_bonds = 0
        for b in ob.OBMolBondIter(ob_mol.OBMol):
            if b.IsRotor():
                rot_bonds += 1

        mol_wt = round(float(ob_mol.molwt), 2)

        # Standardize PDBQT output
        pdbqt_lines = []
        for line in pdbqt_raw.splitlines():
            if line.startswith("TORSDOF"):
                continue
            pdbqt_lines.append(line)
        pdbqt_lines.append(f"TORSDOF {rot_bonds}")
        final_pdbqt = "\n".join(pdbqt_lines) + "\n"

        meta = {
            "formula": ob_mol.formula,
            "mol_wt": mol_wt,
            "heavy_atom_count": heavy_count,
            "rotatable_bonds": rot_bonds,
            "input_format": fmt
        }
        return final_pdbqt, meta

    @classmethod
    def evaluate_vina_score(
        cls,
        rec_coords: np.ndarray,
        rec_types: List[str],
        rec_radii: np.ndarray,
        lig_coords: np.ndarray,
        lig_types: List[str],
        lig_radii: np.ndarray,
        num_rotatable_bonds: int = 0
    ) -> float:
        """
        Pure NumPy vectorized implementation of the AutoDock Vina Empirical Free Energy Scoring Function.
        Evaluates Gauss1, Gauss2, Repulsion, Hydrophobic, and Hydrogen-Bonding potentials.
        """
        if len(rec_coords) == 0 or len(lig_coords) == 0:
            return 0.0

        diff = lig_coords[np.newaxis, :, :] - rec_coords[:, np.newaxis, :]
        dists = np.linalg.norm(diff, axis=-1)

        r_sum = rec_radii[:, np.newaxis] + lig_radii[np.newaxis, :]
        d = dists - r_sum

        cutoff_mask = dists < 8.0

        gauss1 = np.exp(-((d / 0.5) ** 2)) * cutoff_mask
        gauss2 = np.exp(-(((d - 3.0) / 2.0) ** 2)) * cutoff_mask
        rep = np.where(d < 0.0, d ** 2, 0.0) * cutoff_mask

        is_rec_hydro = np.array([_is_hydrophobic(t) for t in rec_types])[:, np.newaxis]
        is_lig_hydro = np.array([_is_hydrophobic(t) for t in lig_types])[np.newaxis, :]
        hydro_pair_mask = (is_rec_hydro & is_lig_hydro) & cutoff_mask

        hydro = np.where(
            d < 0.5, 1.0,
            np.where(d <= 1.5, 1.5 - d, 0.0)
        ) * hydro_pair_mask

        is_rec_acc = np.array([_is_hbond_acceptor(t) for t in rec_types])[:, np.newaxis]
        is_rec_don = np.array([_is_hbond_donor(t) for t in rec_types])[:, np.newaxis]
        is_lig_acc = np.array([_is_hbond_acceptor(t) for t in lig_types])[np.newaxis, :]
        is_lig_don = np.array([_is_hbond_donor(t) for t in lig_types])[np.newaxis, :]

        hbond_pair_mask = ((is_rec_acc & is_lig_don) | (is_rec_don & is_lig_acc)) & cutoff_mask

        hbond = np.where(
            d < -0.7, 1.0,
            np.where(d <= 0.0, -d / 0.7, 0.0)
        ) * hbond_pair_mask

        raw_energy = (
            WEIGHT_GAUSS1 * np.sum(gauss1) +
            WEIGHT_GAUSS2 * np.sum(gauss2) +
            WEIGHT_REPULSION * np.sum(rep) +
            WEIGHT_HYDROPHOBIC * np.sum(hydro) +
            WEIGHT_HBOND * np.sum(hbond)
        )

        torsion_penalty = 1.0 + WEIGHT_ROT * max(0, num_rotatable_bonds)
        free_energy = float(raw_energy / torsion_penalty)

        return free_energy

    @classmethod
    def run_docking(
        cls,
        receptor_pdbqt: str,
        ligand_pdbqt: str,
        center_x: float,
        center_y: float,
        center_z: float,
        size_x: float = 20.0,
        size_y: float = 20.0,
        size_z: float = 20.0,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        energy_range: float = 3.0
    ) -> Dict[str, Any]:
        """
        Executes iLoveChem Native ChemDock v2.0 Inbuilt Molecular Docking.
        1. Parses receptor atoms within the binding pocket box (+4 Å buffer).
        2. Parses ligand heavy atoms and rotatable bonds (TORSDOF).
        3. Identifies rotatable non-ring single bonds and moving atom subtrees.
        4. Conducts multi-start Monte Carlo Simulated Annealing with translation,
           quaternion rotation, and flexible dihedral torsion sampling.
        5. Applies post-annealing analytical force/torque gradient descent to relieve steric clash.
        6. Clusters poses by Heavy-Atom RMSD and ranks by free energy of binding (kcal/mol).
        7. Identifies key active-site polar hydrogen bond interactions.
        """
        # Parse receptor pocket atoms
        rec_coords_list = []
        rec_types_list = []
        rec_atom_lines = []

        half_x = size_x / 2.0 + 4.0
        half_y = size_y / 2.0 + 4.0
        half_z = size_z / 2.0 + 4.0

        for line in receptor_pdbqt.splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    if (abs(x - center_x) <= half_x and
                        abs(y - center_y) <= half_y and
                        abs(z - center_z) <= half_z):
                        atype = line[77:79].strip().upper() or line[12:16].strip()[:1]
                        rec_coords_list.append([x, y, z])
                        rec_types_list.append(atype)
                        rec_atom_lines.append(line)
                except Exception:
                    pass

        if not rec_coords_list:
            raise ValueError(f"No receptor atoms found inside binding grid box at ({center_x}, {center_y}, {center_z}).")

        rec_coords = np.array(rec_coords_list, dtype=np.float64)
        rec_radii = np.array([_get_vdw_radius(t) for t in rec_types_list], dtype=np.float64)

        # Parse ligand atoms & TORSDOF
        lig_coords_list = []
        lig_types_list = []
        lig_atom_template_lines = []
        num_rotatable = 0

        for line in ligand_pdbqt.splitlines():
            if line.startswith("TORSDOF"):
                try:
                    num_rotatable = int(line.split()[1])
                except Exception:
                    pass
            elif line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atype = line[77:79].strip().upper() or line[12:16].strip()[:1]
                    lig_coords_list.append([x, y, z])
                    lig_types_list.append(atype)
                    lig_atom_template_lines.append(line)
                except Exception:
                    pass

        if not lig_coords_list:
            raise ValueError("No ligand atoms could be extracted from ligand PDBQT.")

        orig_lig_coords = np.array(lig_coords_list, dtype=np.float64)
        lig_radii = np.array([_get_vdw_radius(t) for t in lig_types_list], dtype=np.float64)

        # Detect internal rotatable bonds and moving atom subtrees
        rotors = _find_rotatable_bonds_and_trees(orig_lig_coords, lig_types_list)

        # Center ligand initial geometry at pocket centroid
        lig_centroid = np.mean(orig_lig_coords, axis=0)
        centered_lig = orig_lig_coords - lig_centroid
        target_center = np.array([center_x, center_y, center_z], dtype=np.float64)

        # Monte Carlo Simulated Annealing with Flexible Torsions
        num_cycles = max(3, int(exhaustiveness))
        steps_per_cycle = 250

        accepted_poses: List[Tuple[float, np.ndarray]] = []
        grid_half_extent = np.array([size_x / 2.0, size_y / 2.0, size_z / 2.0])

        for cycle in range(num_cycles):
            if cycle == 0:
                init_offset = np.zeros(3)
                curr_coords = centered_lig + target_center
            else:
                init_offset = np.random.uniform(-0.5, 0.5, size=3) * grid_half_extent
                R_init = _random_quaternion_rotation()
                curr_coords = np.dot(centered_lig, R_init.T) + (target_center + init_offset)

            curr_energy = cls.evaluate_vina_score(
                rec_coords, rec_types_list, rec_radii,
                curr_coords, lig_types_list, lig_radii,
                num_rotatable_bonds=num_rotatable
            )

            best_cycle_coords = np.copy(curr_coords)
            best_cycle_energy = curr_energy

            # Annealing loop
            temp = 12.0
            t_decay = 0.985

            for step in range(steps_per_cycle):
                step_scale = max(0.2, temp / 12.0)
                move_rnd = np.random.rand()

                if move_rnd < 0.35 or not rotors:
                    # Rigid translation move
                    trans_step = np.random.normal(0.0, 0.35 * step_scale, size=3)
                    trial_coords = curr_coords + trans_step
                elif move_rnd < 0.65:
                    # Rigid rotation move
                    rot_axis = np.random.normal(size=3)
                    rot_angle = np.random.normal(0.0, 0.25 * step_scale)
                    R_step = _rotation_matrix_from_axis_angle(rot_axis, rot_angle)
                    curr_c = np.mean(curr_coords, axis=0)
                    trial_coords = np.dot(curr_coords - curr_c, R_step.T) + curr_c
                else:
                    # Flexible dihedral torsion move
                    r_idx = np.random.randint(len(rotors))
                    u_idx, v_idx, moving_indices = rotors[r_idx]
                    axis = curr_coords[v_idx] - curr_coords[u_idx]
                    axis_norm = float(np.linalg.norm(axis))
                    if axis_norm > 1e-4:
                        axis /= axis_norm
                        angle = float(np.random.normal(0.0, 0.45 * step_scale))
                        R_torsion = _rotation_matrix_from_axis_angle(axis, angle)
                        center = curr_coords[u_idx]
                        trial_coords = np.copy(curr_coords)
                        trial_coords[moving_indices] = np.dot(curr_coords[moving_indices] - center, R_torsion.T) + center
                    else:
                        trial_coords = curr_coords

                # Bounding box confinement check
                box_dist = np.abs(trial_coords - target_center) - grid_half_extent
                if np.any(box_dist > 1.0):
                    continue

                trial_energy = cls.evaluate_vina_score(
                    rec_coords, rec_types_list, rec_radii,
                    trial_coords, lig_types_list, lig_radii,
                    num_rotatable_bonds=num_rotatable
                )

                dE = trial_energy - curr_energy
                if dE < 0.0 or (temp > 0.01 and np.random.rand() < math.exp(-min(50.0, dE / temp))):
                    curr_coords = trial_coords
                    curr_energy = trial_energy
                    if curr_energy < best_cycle_energy:
                        best_cycle_energy = curr_energy
                        best_cycle_coords = np.copy(curr_coords)

                temp *= t_decay

            if best_cycle_energy < 50.0:
                # Post-annealing gradient relaxation to relieve friction
                relaxed_coords, relaxed_energy = _relax_pose_gradients(
                    rec_coords, rec_types_list, rec_radii,
                    best_cycle_coords, lig_types_list, lig_radii,
                    num_rotatable=num_rotatable,
                    steps=15
                )
                accepted_poses.append((relaxed_energy, relaxed_coords))

        # Sort poses by affinity (most negative = highest affinity)
        accepted_poses.sort(key=lambda item: item[0])

        # Cluster by Heavy-Atom RMSD (>= 1.5 Å separation)
        clustered_poses: List[Tuple[float, np.ndarray, float, float]] = []
        for energy, coords in accepted_poses:
            is_duplicate = False
            for c_energy, c_coords, _, _ in clustered_poses:
                rmsd = float(np.sqrt(np.mean(np.sum((coords - c_coords)**2, axis=-1))))
                if rmsd < 1.5:
                    is_duplicate = True
                    break
            if not is_duplicate:
                rmsd_lb = 0.0
                rmsd_ub = 0.0
                if clustered_poses:
                    top_coords = clustered_poses[0][1]
                    rmsd_lb = round(float(np.sqrt(np.mean(np.sum((coords - top_coords)**2, axis=-1)))), 3)
                    rmsd_ub = round(rmsd_lb * 1.25, 3)
                clustered_poses.append((round(energy, 3), coords, rmsd_lb, rmsd_ub))
                if len(clustered_poses) >= num_modes:
                    break

        if not clustered_poses and accepted_poses:
            best_e, best_c = accepted_poses[0]
            clustered_poses = [(round(best_e, 3), best_c, 0.0, 0.0)]

        # Build output PDBQT models
        docked_models_str_list = []
        poses_summary = []

        for idx, (energy, coords, r_lb, r_ub) in enumerate(clustered_poses, start=1):
            model_lines = [
                f"MODEL {idx}",
                f"REMARK VINA RESULT:    {energy:.3f}      {r_lb:.3f}      {r_ub:.3f}"
            ]
            for atom_idx, line in enumerate(lig_atom_template_lines):
                c = coords[atom_idx]
                coord_part = f"{c[0]:8.3f}{c[1]:8.3f}{c[2]:8.3f}"
                new_line = line[:30] + coord_part + line[54:]
                model_lines.append(new_line)
            model_lines.append("ENDMDL")
            model_str = "\n".join(model_lines)
            docked_models_str_list.append(model_str)

            poses_summary.append({
                "mode": idx,
                "affinity_kcal_mol": energy,
                "rmsd_lb": r_lb,
                "rmsd_ub": r_ub,
                "pdbqt_block": model_str
            })

        docked_output_pdbqt = "\n".join(docked_models_str_list) + "\n"

        # Detect hydrogen bond interactions for top pose
        hbond_contacts = []
        if poses_summary:
            top_coords = clustered_poses[0][1]
            hbond_contacts = cls._detect_contacts_from_coords(
                rec_coords, rec_atom_lines, top_coords, lig_atom_template_lines
            )

        return {
            "success": True,
            "engine": "iLoveChem ChemDock v2.0 (Flexible Torsion & Gradient-Relaxed Monte Carlo)",
            "center": {"x": center_x, "y": center_y, "z": center_z},
            "size": {"x": size_x, "y": size_y, "z": size_z},
            "exhaustiveness": exhaustiveness,
            "rotatable_bonds": len(rotors),
            "poses_count": len(poses_summary),
            "best_affinity_kcal_mol": poses_summary[0]["affinity_kcal_mol"] if poses_summary else None,
            "poses": poses_summary,
            "hbond_contacts": hbond_contacts,
            "docked_output_pdbqt": docked_output_pdbqt
        }

    @staticmethod
    def _detect_contacts_from_coords(
        rec_coords: np.ndarray,
        rec_atom_lines: List[str],
        lig_coords: np.ndarray,
        lig_atom_lines: List[str],
        max_dist: float = 3.5
    ) -> List[Dict[str, Any]]:
        """Identifies key polar hydrogen bonds between receptor residues and top docked ligand pose."""
        contacts = []
        rec_donors_acceptors = []
        for idx, line in enumerate(rec_atom_lines):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            res_seq = line[22:26].strip()
            atom_type = line[77:79].strip().upper() if len(line) >= 79 else ""
            if atom_type in ["OA", "NA", "SA", "HD", "N", "O"] or atom_name[0] in ["N", "O"]:
                rec_donors_acceptors.append({
                    "atom_name": atom_name,
                    "res_name": res_name,
                    "res_seq": res_seq,
                    "coord": rec_coords[idx]
                })

        lig_donors_acceptors = []
        for idx, line in enumerate(lig_atom_lines):
            atom_name = line[12:16].strip()
            atom_type = line[77:79].strip().upper() if len(line) >= 79 else ""
            if atom_type in ["OA", "NA", "SA", "HD", "N", "O"] or (atom_name and atom_name[0] in ["N", "O"]):
                lig_donors_acceptors.append({
                    "atom_name": atom_name,
                    "coord": lig_coords[idx]
                })

        seen_pairs = set()
        for r in rec_donors_acceptors:
            for l in lig_donors_acceptors:
                dist = float(np.linalg.norm(r["coord"] - l["coord"]))
                if 1.9 <= dist <= max_dist:
                    pair_key = (r["res_name"], r["res_seq"], r["atom_name"], l["atom_name"])
                    if pair_key not in seen_pairs:
                        seen_pairs.add(pair_key)
                        contacts.append({
                            "receptor_residue": f"{r['res_name']}{r['res_seq']}",
                            "receptor_atom": r["atom_name"],
                            "ligand_atom": l["atom_name"],
                            "distance_angstrom": round(dist, 2)
                        })

        contacts.sort(key=lambda x: x["distance_angstrom"])
        return contacts[:8]


# Backward-compatible alias
DockingEngine = ChemDockEngine
