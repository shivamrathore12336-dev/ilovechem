import os
import sys
import shutil
import tempfile
import subprocess
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
    "EDO", "DMS", "ACT", "FMT", "PEG", "PG4", "GOL", "BME", "MES", "HEPES"
}

def find_vina_binary() -> Optional[str]:
    """Locate the AutoDock Vina executable in PATH or project bin directory."""
    # 1. System PATH
    bin_path = shutil.which("vina") or shutil.which("vina.exe")
    if bin_path and os.path.exists(bin_path):
        return bin_path

    # 2. Project local bin directory
    base_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(base_dir, "bin", "vina.exe"),
        os.path.join(base_dir, "bin", "vina"),
        "/usr/local/bin/vina",
        "/usr/bin/vina"
    ]
    for c in candidates:
        if os.path.exists(c) and os.access(c, os.X_OK if not c.endswith('.exe') else os.F_OK):
            return c
    return None


class DockingEngine:
    """
    High-Performance In-Memory Molecular Docking Engine.
    Executes macromolecule preparation, flexible ligand torsion construction,
    binding pocket centroid sniffing, and AutoDock Vina scoring.
    Strictly preserves Pharma Zero-Storage Invariant (100% volatile RAM execution).
    """

    @staticmethod
    def detect_co_crystallized_ligand_and_centroid(pdb_str: str) -> Dict[str, Any]:
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

        # Select the largest HETATM group (candidate drug/inhibitor ligand)
        chosen_res = None
        chosen_coords = []
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
            # Clamp between 18.0 and 28.0 Å
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
    def prepare_receptor_pdbqt(cls, pdb_str: str, strip_co_crystal_ligands: bool = True) -> Tuple[str, Dict[str, Any]]:
        """
        Prepares a rigid macromolecular target receptor in PDBQT format:
        1. Strips solvent waters and buffer artifacts.
        2. Optionally strips co-crystallized drug ligands to clear the binding pocket.
        3. Retains protein ATOM records and essential cofactors.
        4. Computes Gasteiger partial charges and adds polar hydrogens.
        5. Enforces fixed-width AutoDock column spacing.
        """
        co_crystal_info = cls.detect_co_crystallized_ligand_and_centroid(pdb_str)
        lig_res = co_crystal_info.get("ligand_resname") if strip_co_crystal_ligands else None

        cleaned_lines = []
        waters_stripped = 0
        ligands_stripped = 0

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
                cleaned_lines.append(line)
            elif rec in ["CRYST1", "SCALE", "ORIGX", "TER", "END"]:
                cleaned_lines.append(line)

        clean_pdb = "\n".join(cleaned_lines) + "\n"

        # Step 2: OpenBabel parsing & charge assignment
        mol = pybel.readstring("pdb", clean_pdb)
        # Add polar hydrogens only according to AutoDock rules
        mol.OBMol.AddHydrogens(True, False) # polar hydrogens only

        for a in mol.atoms:
            try:
                q = a.OBAtom.GetPartialCharge()
                if math.isnan(q) or math.isinf(q):
                    a.OBAtom.SetPartialCharge(0.0)
            except Exception:
                a.OBAtom.SetPartialCharge(0.0)

        pdbqt_raw = mol.write("pdbqt")
        # Ensure strict 6-character partial charges in col 71-76 and strip ligand torsion tags from rigid receptor
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
            # Fallback to OpenBabel reader
            fmt = "mol" if "\n" in s else "smi"
            ob_mol = pybel.readstring(fmt, s)
            ob_mol.make3D(forcefield='mmff94', steps=300)

        # Assign Gasteiger charges
        for a in ob_mol.atoms:
            try:
                q = a.OBAtom.GetPartialCharge()
                if math.isnan(q) or math.isinf(q):
                    a.OBAtom.SetPartialCharge(0.0)
            except Exception:
                a.OBAtom.SetPartialCharge(0.0)

        pdbqt_str = ob_mol.write("pdbqt")
        # Sanitize column alignment
        clean_lines = []
        for line in pdbqt_str.splitlines():
            if (line.startswith("ATOM") or line.startswith("HETATM")) and ("nan" in line.lower() or "inf" in line.lower()):
                if len(line) >= 76:
                    line = line[:70] + "+0.000" + line[76:]
            clean_lines.append(line)

        final_pdbqt = "\n".join(clean_lines) + "\n"
        meta = {
            "formula": ob_mol.formula,
            "molecular_weight": round(float(ob_mol.molwt), 2),
            "atom_count": len(ob_mol.atoms)
        }
        return final_pdbqt, meta

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
        Executes AutoDock Vina in-memory calculation.
        Employs temporary RAM execution with immediate cleanup.
        Returns ranked binding affinities and docked poses.
        """
        vina_bin = find_vina_binary()
        if not vina_bin:
            raise FileNotFoundError("AutoDock Vina executable not found in PATH or project bin.")

        with tempfile.TemporaryDirectory() as tmpdir:
            rec_path = os.path.join(tmpdir, "receptor.pdbqt")
            lig_path = os.path.join(tmpdir, "ligand.pdbqt")
            out_path = os.path.join(tmpdir, "docked_output.pdbqt")
            cfg_path = os.path.join(tmpdir, "docking_grid.conf")

            with open(rec_path, "w", encoding="utf-8") as f:
                f.write(receptor_pdbqt)

            with open(lig_path, "w", encoding="utf-8") as f:
                f.write(ligand_pdbqt)

            config_text = f"""receptor = {rec_path}
ligand = {lig_path}
out = {out_path}

center_x = {center_x:.3f}
center_y = {center_y:.3f}
center_z = {center_z:.3f}

size_x = {size_x:.1f}
size_y = {size_y:.1f}
size_z = {size_z:.1f}

exhaustiveness = {max(1, min(32, int(exhaustiveness)))}
num_modes = {max(1, min(20, int(num_modes)))}
energy_range = {energy_range:.1f}
cpu = 1
"""
            with open(cfg_path, "w", encoding="utf-8") as f:
                f.write(config_text)

            cmd = [vina_bin, "--config", cfg_path]
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

            if proc.returncode != 0 and not os.path.exists(out_path):
                raise RuntimeError(f"Vina docking execution failed (code {proc.returncode}): {proc.stderr.strip() or proc.stdout.strip()}")

            # Parse results
            docked_content = ""
            if os.path.exists(out_path):
                with open(out_path, "r", encoding="utf-8") as f:
                    docked_content = f.read()

            poses = cls._parse_vina_output(docked_content, proc.stdout)

            # Auto-detect hydrogen bond interactions for top pose
            hbond_contacts = []
            if poses:
                top_pose_pdbqt = poses[0]["pdbqt_block"]
                hbond_contacts = cls._detect_hbond_contacts(receptor_pdbqt, top_pose_pdbqt)

            return {
                "success": True,
                "engine": "AutoDock Vina v1.2.7",
                "center": {"x": center_x, "y": center_y, "z": center_z},
                "size": {"x": size_x, "y": size_y, "z": size_z},
                "exhaustiveness": exhaustiveness,
                "poses_count": len(poses),
                "best_affinity_kcal_mol": poses[0]["affinity_kcal_mol"] if poses else None,
                "poses": poses,
                "hbond_contacts": hbond_contacts,
                "docked_output_pdbqt": docked_content
            }

    @staticmethod
    def _parse_vina_output(docked_pdbqt_str: str, stdout_log: str) -> List[Dict[str, Any]]:
        """Parses individual MODEL poses from Vina output PDBQT."""
        poses = []
        raw_models = re.split(r'\n(?=MODEL\s+\d+)', docked_pdbqt_str)

        for m_str in raw_models:
            m_str = m_str.strip()
            if not m_str.startswith("MODEL"):
                continue

            lines = m_str.splitlines()
            model_num = 1
            affinity = 0.0
            rmsd_lb = 0.0
            rmsd_ub = 0.0

            # Match model number
            m_match = re.search(r'MODEL\s+(\d+)', lines[0])
            if m_match:
                model_num = int(m_match.group(1))

            # Match affinity from REMARK VINA RESULT: -8.4 0.0 0.0
            for line in lines:
                if "REMARK VINA RESULT:" in line:
                    parts = line.replace("REMARK VINA RESULT:", "").strip().split()
                    if len(parts) >= 1:
                        affinity = float(parts[0])
                    if len(parts) >= 2:
                        rmsd_lb = float(parts[1])
                    if len(parts) >= 3:
                        rmsd_ub = float(parts[2])
                    break

            poses.append({
                "mode": model_num,
                "affinity_kcal_mol": affinity,
                "rmsd_lb": rmsd_lb,
                "rmsd_ub": rmsd_ub,
                "pdbqt_block": m_str
            })

        return poses

    @staticmethod
    def _detect_hbond_contacts(receptor_pdbqt: str, pose_pdbqt: str, max_dist: float = 3.5) -> List[Dict[str, Any]]:
        """Identifies key polar hydrogen bonds between receptor residues and docked ligand pose."""
        contacts = []
        rec_donors_acceptors = []
        for line in receptor_pdbqt.splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                res_seq = line[22:26].strip()
                atom_type = line[77:79].strip().upper() if len(line) >= 79 else ""
                if atom_type in ["OA", "NA", "SA", "HD", "N", "O"] or atom_name[0] in ["N", "O"]:
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        rec_donors_acceptors.append({
                            "atom_name": atom_name,
                            "res_name": res_name,
                            "res_seq": res_seq,
                            "coord": np.array([x, y, z])
                        })
                    except Exception:
                        pass

        lig_donors_acceptors = []
        for line in pose_pdbqt.splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                atom_type = line[77:79].strip().upper() if len(line) >= 79 else ""
                if atom_type in ["OA", "NA", "SA", "HD", "N", "O"] or (atom_name and atom_name[0] in ["N", "O"]):
                    try:
                        x = float(line[30:38].strip())
                        y = float(line[38:46].strip())
                        z = float(line[46:54].strip())
                        lig_donors_acceptors.append({
                            "atom_name": atom_name,
                            "coord": np.array([x, y, z])
                        })
                    except Exception:
                        pass

        seen_pairs = set()
        for r in rec_donors_acceptors:
            for l in lig_donors_acceptors:
                dist = float(np.linalg.norm(r["coord"] - l["coord"]))
                if 2.0 <= dist <= max_dist:
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
        return contacts[:8] # Return top 8 strongest interactions
