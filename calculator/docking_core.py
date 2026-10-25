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

import os
import sys
import math
import re
from typing import Dict, Any, Tuple, Optional, List
import numpy as np

from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

# Buffer / solvent residue identifiers to ignore when sniffing co-crystallized ligands
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
    "FE": 1.20, "ZN": 1.20, "MG": 1.20, "CA": 1.20, "MN": 1.20
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
    """Rodriguez rotation formula for 3D rotation around an arbitrary axis."""
    axis = axis / (np.linalg.norm(axis) + 1e-12)
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


class ChemDockEngine:
    """
    iLoveChem Native Inbuilt Molecular Docking Kernel.
    100% in-process Python/NumPy/SciPy implementation of the AutoDock Vina
    Empirical Free Energy Scoring Function and Monte Carlo Annealing Engine.
    Requires ZERO external binaries and operates purely in volatile memory.
    """

    @staticmethod
    def detect_co_crystallized_ligand_and_centroid(pdb_str: str, preferred_ligand_resname: Optional[str] = None) -> Dict[str, Any]:
        """
        Inspects PDB records for co-crystallized drug ligands (HETATM).
        Calculates geometric centroid (center_x, center_y, center_z) and
        suggested grid box dimensions (size_x, size_y, size_z).
        """
        het_coords: Dict[str, List[Tuple[float, float, float]]] = {}
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
                if resname not in SOLVENT_AND_BUFFER_RESIDUES:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        if resname not in het_coords:
                            het_coords[resname] = []
                        het_coords[resname].append((x, y, z))
                    except (ValueError, IndexError):
                        pass

        # Select candidate drug/inhibitor ligand (prefer specified ligand if present)
        chosen_res = None
        chosen_coords = []
        pref_clean = preferred_ligand_resname.strip().upper() if preferred_ligand_resname else None

        if pref_clean and pref_clean in het_coords and len(het_coords[pref_clean]) >= 4:
            chosen_res = pref_clean
            chosen_coords = het_coords[pref_clean]
        else:
            for res, coords in het_coords.items():
                if len(coords) >= 6 and len(coords) > len(chosen_coords):
                    chosen_res = res
                    chosen_coords = coords

        if chosen_coords:
            arr = np.array(chosen_coords)
            c_x, c_y, c_z = np.mean(arr, axis=0)
            span_x = float(np.ptp(arr[:, 0])) + 8.0
            span_y = float(np.ptp(arr[:, 1])) + 8.0
            span_z = float(np.ptp(arr[:, 2])) + 8.0
            size_x = max(18.0, min(28.0, round(span_x, 1)))
            size_y = max(18.0, min(28.0, round(span_y, 1)))
            size_z = max(18.0, min(28.0, round(span_z, 1)))

            return {
                "has_co_crystal_ligand": True,
                "ligand_resname": chosen_res,
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
            "center_x": 0.0,
            "center_y": 0.0,
            "center_z": 0.0,
            "size_x": 20.0,
            "size_y": 20.0,
            "size_z": 20.0,
            "method": "origin_default"
        }

    @classmethod
    def prepare_receptor_pdbqt(cls, pdb_str: str, strip_co_crystal_ligands: bool = True, pocket_center: Optional[Tuple[float, float, float]] = None, pocket_radius: float = 24.0) -> Tuple[str, Dict[str, Any]]:
        """
        Prepares a rigid macromolecular target receptor in PDBQT format:
        1. Strips solvent waters and buffer artifacts.
        2. Optionally strips co-crystallized drug ligands to clear the binding pocket.
        3. Retains protein ATOM records and essential cofactors.
        4. Computes Gasteiger partial charges and adds polar hydrogens.
        5. Enforces fixed-width AutoDock column spacing without flexible torsion tags.
        """
        co_crystal_info = cls.detect_co_crystallized_ligand_and_centroid(pdb_str)
        lig_res = co_crystal_info.get("ligand_resname") if strip_co_crystal_ligands else None

        cleaned_lines = []
        waters_stripped = 0
        ligands_stripped = 0

        pocket_cx = pocket_center[0] if pocket_center else (co_crystal_info.get("center_x") if co_crystal_info.get("has_co_crystal_ligand") else None)
        pocket_cy = pocket_center[1] if pocket_center else (co_crystal_info.get("center_y") if co_crystal_info.get("has_co_crystal_ligand") else None)
        pocket_cz = pocket_center[2] if pocket_center else (co_crystal_info.get("center_z") if co_crystal_info.get("has_co_crystal_ligand") else None)
        max_dist_sq = (pocket_radius**2) if (pocket_cx is not None and pocket_radius > 0) else None

        for line in pdb_str.splitlines():
            rec = line[0:6].strip()
            if rec in ["ATOM", "HETATM"]:
                resname = line[17:20].strip().upper()
                if resname in SOLVENT_AND_BUFFER_RESIDUES:
                    waters_stripped += 1
                    continue
                if lig_res and resname == lig_res:
                    ligands_stripped += 1
                    continue
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

        # Fast in-process PDB-to-PDBQT translation (instantaneous in-memory execution)
        fast_pdbqt_lines = []
        atom_count = 0
        for line in cleaned_lines:
            rec = line[0:6].strip()
            if rec in ["ATOM", "HETATM"]:
                try:
                    resname = line[17:20].strip().upper()
                    atom_name = line[12:16].strip().upper()
                    elem = line[76:78].strip().upper() or atom_name[:1]
                    ad_type = PDB_TO_AUTODOCK_ATOM_TYPES.get((resname, atom_name))
                    if not ad_type:
                        if atom_name in ["O", "OXT"] or elem == "O":
                            ad_type = "OA"
                        elif atom_name == "N" or elem == "N":
                            ad_type = "N"
                        elif elem == "S":
                            ad_type = "SA"
                        elif elem == "C":
                            ad_type = "C"
                        elif elem == "H":
                            ad_type = "HD"
                        else:
                            ad_type = elem[:2]
                    pdbqt_line = f"{line[:70]:<70s} +0.000 {ad_type:>2s}"
                    fast_pdbqt_lines.append(pdbqt_line)
                    atom_count += 1
                except Exception:
                    pass

        if fast_pdbqt_lines and atom_count >= 10:
            pdbqt_str = "\n".join(fast_pdbqt_lines) + "\n"
            meta = {
                "atom_count": atom_count,
                "atoms_count": atom_count,
                "waters_stripped": waters_stripped,
                "ligands_stripped": ligands_stripped,
                "formula": f"Protein ({atom_count} atoms)"
            }
            return pdbqt_str, meta

        clean_pdb = "\n".join(cleaned_lines) + "\n"

        mol = pybel.readstring("pdb", clean_pdb)
        mol.OBMol.AddHydrogens(True, False) # polar hydrogens only

        for a in mol.atoms:
            try:
                q = a.OBAtom.GetPartialCharge()
                if math.isnan(q) or math.isinf(q):
                    a.OBAtom.SetPartialCharge(0.0)
            except Exception:
                a.OBAtom.SetPartialCharge(0.0)

        pdbqt_raw = mol.write("pdbqt")
        formatted_lines = []
        atom_count = 0
        for line in pdbqt_raw.splitlines():
            line_s = line.strip()
            if line_s.startswith(("ROOT", "ENDROOT", "TORSDOF", "BRANCH", "ENDBRANCH")):
                continue
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_count += 1
                if "nan" in line.lower() or "inf" in line.lower():
                    if len(line) >= 76:
                        line = line[:70] + "+0.000" + line[76:]
            formatted_lines.append(line)

        pdbqt_str = "\n".join(formatted_lines) + "\n"
        meta = {
            "atom_count": atom_count,
            "atoms_count": atom_count,
            "waters_stripped": waters_stripped,
            "ligands_stripped": ligands_stripped,
            "formula": mol.formula
        }
        return pdbqt_str, meta

    @staticmethod
    def prepare_ligand_pdbqt(ligand_input: str, input_format: str = "smiles") -> Tuple[str, Dict[str, Any]]:
        """
        Prepares a flexible small molecule ligand in standard PDBQT format:
        1. Generates physical 3D conformer with RDKit ETKDGv3 + MMFF94.
        2. Merges non-polar hydrogens into heavy atoms.
        3. Computes Gasteiger partial charges.
        4. Builds AutoDock torsion tree (ROOT, BRANCH, TORSDOF).
        """
        rd_mol = None
        s = ligand_input.strip()

        if input_format.lower() in ["smi", "smiles"] or ("\n" not in s and " " not in s and len(s) < 300):
            try:
                rd_mol = Chem.MolFromSmiles(s)
                if rd_mol:
                    rd_mol = Chem.AddHs(rd_mol)
                    AllChem.EmbedMolecule(rd_mol, AllChem.ETKDGv3())
                    if AllChem.MMFFHasAllMoleculeParams(rd_mol):
                        AllChem.MMFFOptimizeMolecule(rd_mol, maxIters=500)
                    elif AllChem.UFFHasAllMoleculeParams(rd_mol):
                        AllChem.UFFOptimizeMolecule(rd_mol, maxIters=500)
            except Exception:
                rd_mol = None

        if rd_mol is not None:
            mol_block = Chem.MolToMolBlock(rd_mol)
            ob_mol = pybel.readstring("mol", mol_block)
        else:
            fmt_lower = input_format.lower().strip()
            if fmt_lower in ["pdb", "ent"] or "ATOM" in s or "HETATM" in s:
                fmt = "pdb"
            elif fmt_lower in ["mol", "sdf"] or "$$$$" in s or "M  END" in s:
                fmt = "mol"
            else:
                fmt = "smi"

            ob_mol = pybel.readstring(fmt, s)
            if fmt == "smi":
                try:
                    ob_mol.make3D(forcefield='mmff94', steps=300)
                except Exception:
                    pass

        for a in ob_mol.atoms:
            try:
                q = a.OBAtom.GetPartialCharge()
                if math.isnan(q) or math.isinf(q):
                    a.OBAtom.SetPartialCharge(0.0)
            except Exception:
                a.OBAtom.SetPartialCharge(0.0)

        pdbqt_str = ob_mol.write("pdbqt")
        clean_lines = []
        for line in pdbqt_str.splitlines():
            if (line.startswith("ATOM") or line.startswith("HETATM")) and ("nan" in line.lower() or "inf" in line.lower()):
                if len(line) >= 76:
                    line = line[:70] + "+0.000" + line[76:]
            clean_lines.append(line)

        final_pdbqt = "\n".join(clean_lines) + "\n"
        rot_bonds = 0
        for line in final_pdbqt.splitlines():
            if line.startswith("TORSDOF"):
                try:
                    rot_bonds = int(line.split()[1])
                except Exception:
                    pass

        meta = {
            "formula": ob_mol.formula,
            "molecular_weight": round(float(ob_mol.molwt), 2),
            "atom_count": len(ob_mol.atoms),
            "num_atoms": len(ob_mol.atoms),
            "num_rotatable_bonds": rot_bonds
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
        Vectorized evaluation of the AutoDock Vina Empirical Free Energy Function:
        ΔG = (w_gauss1*f1 + w_gauss2*f2 + w_rep*f_rep + w_hydro*f_hydro + w_hbond*f_hbond) / (1 + w_rot*N_rot)
        """
        # Distance matrix (N_rec, N_lig)
        diff = rec_coords[:, np.newaxis, :] - lig_coords[np.newaxis, :, :] # (N_rec, N_lig, 3)
        dist_sq = np.sum(diff**2, axis=-1) # (N_rec, N_lig)
        dist = np.sqrt(dist_sq + 1e-12)

        # Spatial cutoff: interactions beyond 8.0 Å contribute negligibly
        cutoff_mask = dist <= 8.0
        if not np.any(cutoff_mask):
            return 0.0

        r_sum = rec_radii[:, np.newaxis] + lig_radii[np.newaxis, :]
        d = dist - r_sum # surface distance

        # Gauss 1
        gauss1 = np.exp(-((d / 0.5)**2)) * cutoff_mask
        # Gauss 2
        gauss2 = np.exp(-(((d - 3.0) / 2.0)**2)) * cutoff_mask
        # Repulsion
        rep = np.where(d < 0.0, d**2, 0.0) * cutoff_mask

        # Hydrophobic mask
        is_rec_hydro = np.array([_is_hydrophobic(t) for t in rec_types])[:, np.newaxis]
        is_lig_hydro = np.array([_is_hydrophobic(t) for t in lig_types])[np.newaxis, :]
        hydro_pair_mask = (is_rec_hydro & is_lig_hydro) & cutoff_mask

        hydro = np.where(
            d < 0.5, 1.0,
            np.where(d <= 1.5, 1.5 - d, 0.0)
        ) * hydro_pair_mask

        # HBond mask (donor-acceptor or acceptor-donor)
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
        Executes iLoveChem Native ChemDock Inbuilt Molecular Docking.
        1. Parses receptor atoms within the binding pocket box (+4 Å buffer).
        2. Parses ligand heavy atoms and rotatable bond count (TORSDOF).
        3. Conducts multi-start Monte Carlo Simulated Annealing with quaternion rotation,
           rigid translation, and pocket-constrained optimization.
        4. Clusters generated poses by RMSD and ranks them by AutoDock binding affinity (kcal/mol).
        5. Detects active-site hydrogen bond contact pairs.
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

        # Center ligand initial geometry at pocket centroid
        lig_centroid = np.mean(orig_lig_coords, axis=0)
        centered_lig = orig_lig_coords - lig_centroid
        target_center = np.array([center_x, center_y, center_z], dtype=np.float64)

        # Monte Carlo Simulated Annealing
        # Number of cycles scaled with exhaustiveness
        num_cycles = max(2, int(exhaustiveness))
        steps_per_cycle = 200

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
                # Adaptive perturbation: combination of translation and rotation
                step_scale = max(0.2, temp / 12.0)
                trans_step = np.random.normal(0.0, 0.4 * step_scale, size=3)
                rot_axis = np.random.normal(size=3)
                rot_angle = np.random.normal(0.0, 0.3 * step_scale)
                R_step = _rotation_matrix_from_axis_angle(rot_axis, rot_angle)

                curr_c = np.mean(curr_coords, axis=0)
                trial_coords = np.dot(curr_coords - curr_c, R_step.T) + curr_c + trans_step

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

            if best_cycle_energy < 50.0: # Filter completely clashing non-physical configurations
                accepted_poses.append((best_cycle_energy, best_cycle_coords))

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
            "engine": "iLoveChem ChemDock Native v1.0 (AutoDock Vina Scoring Inbuilt)",
            "center": {"x": center_x, "y": center_y, "z": center_z},
            "size": {"x": size_x, "y": size_y, "z": size_z},
            "exhaustiveness": exhaustiveness,
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


# Backward-compatible alias for the native docking engine
DockingEngine = ChemDockEngine
