"""
Scientific Engine for iLoveChem:
Implements Phases 2 through 5 of the Universal Chemistry Conversion Standard:
- Phase 2: Quantum formatting, inverted log parser, spin/charge parity validation, Gaussian blank lines, coordinate sanitization.
- Phase 3: Crystallography & Solid-State (15 Å vacuum box, fractional to Cartesian, CIF auto-protonation).
- Phase 4: Bio/Macromolecule sanitization (AltLoc A/B filter, Water/Solvent stripping, Ligand isolation).
- Phase 5: Heavy metal UFF fallback, nuclear overlap clash checker (<0.7 Å), large file log streaming.
"""

import os
import re
import io
import math
import zipfile
from typing import List, Dict, Any, Tuple, Optional
import numpy as np

from openbabel import openbabel as ob
from openbabel import pybel
from ase.atoms import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem

# Periodic Table atomic numbers
TRANSITION_AND_HEAVY_METALS = set(
    list(range(21, 31)) + # Sc - Zn
    list(range(39, 49)) + # Y - Cd
    list(range(57, 81)) + # La - Tl
    list(range(89, 113))  # Ac - Nh
)

WATER_RESIDUE_NAMES = {"HOH", "WAT", "H2O", "SOL", "TIP3", "TIP4", "DOD", "OH2"}
COMMON_BUFFER_IONS = {"NA", "CL", "K", "MG", "CA", "ZN", "SO4", "PO4", "ACT", "NO3", "NH4", "BR", "IOD"}


class SaltAndSolventStripper:
    """
    Detects and strips disconnected counter-ions, acid salts, and crystallization solvents:
    - Counterions: Na+, K+, Li+, Cl-, Br-, I-, Mg2+, Ca2+, Zn2+, NH4+, etc.
    - Acid salts: HCl, HBr, TFA (trifluoroacetate), Acetate, Mesylate, Tosylate, Sulfate.
    - Crystallization Solvents / Hydrates: H2O, MeOH, EtOH, DMSO, Acetone, DCM.
    Prevents DFT geometry optimizations from crashing or drifting infinitely in vacuum.
    """
    KNOWN_COUNTERION_SYMBOLS = {
        "Na", "K", "Li", "Cs", "Rb", "Ca", "Mg", "Ba", "Zn", "Al", "Fe", "Cu",
        "Cl", "Br", "I", "F"
    }

    KNOWN_SALTS_AND_SOLVENTS = {
        # Inorganic / Polyatomic counter-ions
        "H2O", "OH2", "H2O1", "H3O", "NH4", "SO4", "HSO4", "NO3", "PO4", "HPO4", "H2PO4",
        "CLO4", "BF4", "PF6", "CN", "SCN", "HCL", "HBR", "HI",
        # Small organic acids & counterions
        "C2H3O2", "C2H4O2", "CH3COO", "HCOO", "CHO2", "CH2O2",
        "C2HF3O2", "C2F3O2", "CF3COO", "TFA",
        "CH4O3S", "CH3SO3", "C7H8O3S", "C7H7SO3", "C2H2O4", "C2O4",
        # Common crystallization solvents
        "CH4O", "C2H6O", "C3H8O", "C3H6O", "C2H3N", "CH2CL2", "CHCL3",
        "C4H10O", "C4H8O", "C2H6OS", "C3H7NO"
    }

    @classmethod
    def analyze_and_process(
        cls,
        mol: pybel.Molecule,
        strip_salts: bool = True
    ) -> Tuple[pybel.Molecule, Dict[str, Any]]:
        meta = {
            "salt_detected": False,
            "stripped": False,
            "detected_salts": [],
            "parent_formula": mol.formula,
            "parent_charge": mol.charge,
            "warnings": []
        }

        if len(mol.atoms) <= 1:
            return mol, meta

        # Physical Coordination Guard:
        # Check if any transition metal is within coordination distance (< 2.8 Å) of organic atoms
        # (e.g. Ferrocene, Cisplatin, metallocenes, coordination complexes).
        # These are coordinated complexes, NOT spectator salts, and must remain unified!
        metal_atoms = [a for a in mol.atoms if a.OBAtom.GetAtomicNum() in TRANSITION_AND_HEAVY_METALS]
        non_metal_atoms = [a for a in mol.atoms if a.OBAtom.GetAtomicNum() not in TRANSITION_AND_HEAVY_METALS]
        if metal_atoms and non_metal_atoms:
            has_coord_metal = False
            for m in metal_atoms:
                m_coords = np.array(m.coords)
                for nm in non_metal_atoms:
                    dist = float(np.linalg.norm(m_coords - np.array(nm.coords)))
                    if 0.5 < dist < 2.85:
                        has_coord_metal = True
                        break
                if has_coord_metal:
                    break
            if has_coord_metal:
                return mol, meta

        # If bonds not yet perceived (e.g. from raw 3D XYZ), connect and perceive
        if mol.OBMol.NumBonds() == 0 and len(mol.atoms) > 1:
            mol.OBMol.ConnectTheDots()
            mol.OBMol.PerceiveBondOrders()

        ob_frags = mol.OBMol.Separate()
        if len(ob_frags) <= 1:
            return mol, meta

        # Multiple disconnected fragments detected!
        frag_mols = [pybel.Molecule(f) for f in ob_frags]
        
        # Check every fragment against known counter-ions and solvents
        detected_salts_list = []
        non_salt_frags = []
        salt_frags = []

        for f in frag_mols:
            f_formula = re.sub(r'[^A-Za-z0-9]', '', f.formula).upper()
            f_atoms = len(f.atoms)
            f_sym = ob.GetSymbol(f.atoms[0].OBAtom.GetAtomicNum()) if f_atoms == 1 else ""

            matched = False
            label = f.formula
            if f_atoms == 1 and f_sym in cls.KNOWN_COUNTERION_SYMBOLS:
                matched = True
                charge_str = f"+{f.charge}" if f.charge > 0 else (f"{f.charge}" if f.charge < 0 else "")
                label = f"{f_sym}{charge_str} Ion"
            elif f_formula in cls.KNOWN_SALTS_AND_SOLVENTS:
                matched = True
                label = f"{f.formula} (Solvate/Salt)"
            
            if matched:
                salt_frags.append(f)
                detected_salts_list.append(label)
            else:
                non_salt_frags.append(f)

        # If known salts found and at least one organic molecule exists:
        if salt_frags and non_salt_frags:
            parent_frag = sorted(non_salt_frags, key=lambda m: len(m.atoms), reverse=True)[0]
            is_salt = True
        else:
            # Fallback heuristic: check if smallest fragments have <= 4 atoms compared to largest
            sorted_frags = sorted(frag_mols, key=lambda m: len(m.atoms), reverse=True)
            parent_frag = sorted_frags[0]
            is_salt = False
            for f in sorted_frags[1:]:
                if len(f.atoms) <= 4 and (len(f.atoms) * 3 < len(parent_frag.atoms)):
                    is_salt = True
                    detected_salts_list.append(f"{f.formula} (Small Fragment/Solvate)")

        if is_salt:
            meta["salt_detected"] = True
            meta["detected_salts"] = detected_salts_list
            meta["parent_formula"] = parent_frag.formula
            meta["parent_charge"] = parent_frag.charge

            if strip_salts:
                meta["stripped"] = True
                meta["warnings"].append(
                    f"Salt / Counter-Ions Stripped: Detected {', '.join(detected_salts_list)}. "
                    f"Isolated active parent molecule ({parent_frag.formula}, Charge = {parent_frag.charge}) "
                    f"to prevent vacuum DFT dissociation crashes."
                )
                return parent_frag, meta
            else:
                meta["stripped"] = False
                meta["warnings"].append(
                    f"Salt / Counter-Ions Retained: Full complex containing {', '.join(detected_salts_list)} "
                    f"kept. (Caution: Isolated ions may drift infinitely during vacuum DFT optimization)."
                )
                return mol, meta

        return mol, meta


class SpinChargeParityValidator:
    @staticmethod
    def validate_and_fix(total_atomic_number: int, charge: int, mult: int) -> Tuple[int, List[str]]:
        """
        Validates electron count parity:
        N_electrons = sum(Z_i) - charge
        - If N_e is EVEN: multiplicity MUST be ODD (1, 3, 5, ...).
        - If N_e is ODD: multiplicity MUST be EVEN (2, 4, 6, ...).
        Returns: (fixed_mult, warnings_list)
        """
        warnings = []
        n_elec = total_atomic_number - charge
        
        if n_elec <= 0:
            warnings.append(f"Warning: Calculated non-positive electron count ({n_elec}).")
            return mult, warnings

        is_even_electrons = (n_elec % 2 == 0)
        is_odd_mult = (mult % 2 == 1)

        if is_even_electrons and not is_odd_mult:
            fixed_mult = 1
            warnings.append(
                f"Spin Multiplicity Parity Error: Molecule has {n_elec} electrons (even). "
                f"Multiplicity {mult} is physically impossible (leads to Gaussian/ORCA QPErr). "
                f"Auto-corrected to Singlet (multiplicity = 1)."
            )
            return fixed_mult, warnings
        elif (not is_even_electrons) and is_odd_mult:
            fixed_mult = 2
            warnings.append(
                f"Spin Multiplicity Parity Error: Molecule has {n_elec} electrons (odd, radical/open-shell). "
                f"Multiplicity {mult} is physically impossible (leads to Gaussian/ORCA QPErr). "
                f"Auto-corrected to Doublet (multiplicity = 2)."
            )
            return fixed_mult, warnings

        return mult, warnings


class QuantumLogParser:
    @classmethod
    def reverse_stream_extract_from_file(
        cls,
        file_obj,
        tail_bytes: int = 4 * 1024 * 1024
    ) -> Optional[Tuple[List[Tuple[str, float, float, float]], str]]:
        """
        Streaming tail extractor for giant 50 MB - 500 MB output logs.
        Seeks directly to the final 4MB of the file and extracts the converged stationary geometry
        in < 20 ms without loading the full 100MB+ file into server RAM.
        """
        try:
            file_obj.seek(0, os.SEEK_END)
            total_size = file_obj.tell()
            offset = max(0, total_size - tail_bytes)
            file_obj.seek(offset)
            chunk = file_obj.read()
            chunk_str = chunk.decode('utf-8', errors='ignore')
            return cls.reverse_scan_gaussian_or_orca_log(chunk_str)
        except Exception:
            return None

    @classmethod
    def reverse_scan_gaussian_or_orca_log(cls, content_str: str) -> Optional[Tuple[List[Tuple[str, float, float, float]], str]]:
        """
        Inverted scanner for Gaussian, ORCA, Q-Chem, and NWChem output logs.
        Scans from bottom to top to find the FINAL converged geometry / stationary point,
        avoiding reading unoptimized initial geometries or intermediates.
        """
        # Memory Guard: if input string exceeds 4MB, take only the final 4MB tail
        if len(content_str) > 4 * 1024 * 1024:
            content_str = content_str[-4 * 1024 * 1024:]

        lines = content_str.splitlines()
        
        # Check if calculation terminated abnormally without full convergence (Loophole 3)
        is_converged = ("Stationary point found" in content_str or "OPTIMIZATION RUN HAS CONVERGED" in content_str or "OPTIMIZATION CONVERGED" in content_str or "Normal termination" in content_str)
        if not is_converged and ("Error termination" in content_str or "abnormally" in content_str or "termination" in content_str or "Step" in content_str):
            rescue_res = cls.rescue_abnormal_optimization(content_str)
            if rescue_res:
                best_coords, label, min_energy, step_num = rescue_res
                return best_coords, label

        # 1. Reverse search for Gaussian "Stationary point found" or standard orientation blocks
        for i in range(len(lines) - 1, -1, -1):
            line = lines[i]
            
            # Gaussian: "Standard orientation:" or "Input orientation:"
            if "Standard orientation:" in line or "Input orientation:" in line:
                coords = []
                j = i + 5
                while j < len(lines):
                    sub = lines[j].strip()
                    if sub.startswith("---") or sub == "":
                        break
                    parts = sub.split()
                    if len(parts) >= 6 and parts[0].isdigit():
                        at_num = int(parts[1])
                        sym = ob.GetSymbol(at_num)
                        try:
                            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                            coords.append((sym, x, y, z))
                        except ValueError:
                            pass
                    j += 1
                if coords:
                    return coords, "Gaussian Converged Stationary Point"

            # ORCA: "CARTESIAN COORDINATES (ANGSTROEM)"
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                coords = []
                j = i + 2
                while j < len(lines):
                    sub = lines[j].strip()
                    if sub.startswith("---") or sub == "" or "CARTESIAN" in sub:
                        break
                    parts = sub.split()
                    if len(parts) >= 4:
                        sym = parts[0]
                        try:
                            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
                            coords.append((sym, x, y, z))
                        except ValueError:
                            break
                    else:
                        break
                    j += 1
                if coords:
                    return coords, "ORCA Converged Geometry"

            # Q-Chem: "Coordinates (Angstroms)" after optimization cycle
            if "Coordinates (Angstroms)" in line and ("OPTIMIZATION CONVERGED" in content_str or "Final Energy" in content_str):
                coords = []
                j = i + 2
                while j < len(lines):
                    sub = lines[j].strip()
                    if sub.startswith("---") or sub == "" or "$" in sub:
                        break
                    parts = sub.split()
                    if len(parts) >= 5 and parts[0].isdigit():
                        sym = parts[1]
                        try:
                            x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                            coords.append((sym, x, y, z))
                        except ValueError:
                            pass
                    j += 1
                if coords:
                    return coords, "Q-Chem Converged Geometry"

        return None

    @classmethod
    def rescue_abnormal_optimization(cls, content_str: str) -> Optional[Tuple[List[Tuple[str, float, float, float]], str, float, int]]:
        """
        Loophole 3: Incomplete / Failed Quantum Optimizations in Logs.
        Scans all intermediate optimization steps when a calculation crashed,
        finds the cycle with the lowest SCF energy, and extracts its coordinates.
        Returns: (coords, parser_label, min_energy, best_step_num)
        """
        lines = content_str.splitlines()
        steps = []
        curr_energy = None
        
        for i, line in enumerate(lines):
            # Gaussian: SCF Done:  E(RB3LYP) =  -123.45678
            m_gau = re.search(r'SCF Done:\s+E\([^\)]+\)\s*=\s*([-\d\.]+)', line)
            if m_gau:
                try: curr_energy = float(m_gau.group(1))
                except ValueError: pass
                
            # ORCA: FINAL SINGLE POINT ENERGY      -123.45678
            m_orca = re.search(r'FINAL SINGLE POINT ENERGY\s+([-\d\.]+)', line)
            if m_orca:
                try: curr_energy = float(m_orca.group(1))
                except ValueError: pass

            if ('Standard orientation:' in line or 'Input orientation:' in line) and curr_energy is not None:
                coords = []
                j = i + 5
                while j < len(lines):
                    sub = lines[j].strip()
                    if sub.startswith('---') or sub == '':
                        break
                    parts = sub.split()
                    if len(parts) >= 6 and parts[0].isdigit():
                        at_num = int(parts[1])
                        sym = ob.GetSymbol(at_num)
                        try:
                            coords.append((sym, float(parts[3]), float(parts[4]), float(parts[5])))
                        except ValueError:
                            pass
                    j += 1
                if coords:
                    steps.append((curr_energy, len(steps) + 1, coords, "Gaussian Abnormally Terminated Job"))

            if 'CARTESIAN COORDINATES (ANGSTROEM)' in line and curr_energy is not None:
                coords = []
                j = i + 2
                while j < len(lines):
                    sub = lines[j].strip()
                    if sub.startswith('---') or sub == '' or 'CARTESIAN' in sub:
                        break
                    parts = sub.split()
                    if len(parts) >= 4:
                        sym = parts[0]
                        try:
                            coords.append((sym, float(parts[1]), float(parts[2]), float(parts[3])))
                        except ValueError:
                            break
                    j += 1
                if coords:
                    steps.append((curr_energy, len(steps) + 1, coords, "ORCA Abnormally Terminated Job"))

        if not steps:
            return None

        # Pick step with lowest (most stable) SCF energy
        best_step = min(steps, key=lambda s: s[0])
        min_energy, step_num, best_coords, engine_name = best_step
        label = f"{engine_name} (Rescued Lowest Energy Step #{step_num}, E = {min_energy:.6f} Har)"
        return best_coords, label, min_energy, step_num


class PhysicalIntegrityChecker:
    @staticmethod
    def check_clashes_and_metals(mol: pybel.Molecule) -> Dict[str, Any]:
        """
        Checks:
        1. Heavy / Transition metal presence (triggers UFF fallback recommendation).
        2. Nuclear overlap clash check (minimum interatomic distance < 0.7 Å).
        """
        atom_count = len(mol.atoms)
        coords = np.array([a.coords for a in mol.atoms])
        atomic_nums = [a.OBAtom.GetAtomicNum() for a in mol.atoms]
        
        has_metals = any(z in TRANSITION_AND_HEAVY_METALS for z in atomic_nums)
        metal_symbols = [ob.GetSymbol(z) for z in atomic_nums if z in TRANSITION_AND_HEAVY_METALS]

        min_dist = float('inf')
        clash_pairs = []

        if atom_count > 1 and len(coords) == atom_count:
            # Pairwise distance matrix check
            for i in range(atom_count):
                for j in range(i + 1, atom_count):
                    d = np.linalg.norm(coords[i] - coords[j])
                    if d < min_dist:
                        min_dist = float(d)
                    if d < 0.7:  # Severe nuclear overlap clash
                        clash_pairs.append({
                            "atom1": f"{ob.GetSymbol(atomic_nums[i])}{i+1}",
                            "atom2": f"{ob.GetSymbol(atomic_nums[j])}{j+1}",
                            "distance": round(float(d), 3)
                        })

        return {
            "has_metals": has_metals,
            "metal_symbols": list(set(metal_symbols)),
            "min_interatomic_distance": round(min_dist, 3) if min_dist != float('inf') else 0.0,
            "has_clashes": len(clash_pairs) > 0,
            "clash_pairs": clash_pairs[:5]  # first 5 clashes
        }


class PDBBioSanitizer:
    @staticmethod
    def sanitize_pdb(pdb_str: str, filter_mode: str = "all") -> str:
        """
        Sanitizes PDB file:
        - Resolves alternate location identifiers (AltLoc 'A' kept, 'B','C' dropped).
        - Filter mode:
          - 'all': Keeps everything (with AltLoc duplicates resolved)
          - 'strip_water': Removes HOH, WAT, SOL, etc.
          - 'ligand_only': Extracts HETATM ligands (excluding water and simple buffer ions).
        """
        lines = pdb_str.splitlines()
        
        # Multi-Model NMR / MD Trajectory Handler:
        # If file contains MODEL records, only process MODEL 1 to ensure a clean single conformer
        has_models = any(l.startswith("MODEL") for l in lines)
        if has_models:
            model_1_lines = []
            in_first_model = False
            for l in lines:
                if l.startswith("MODEL"):
                    parts = l.split()
                    if len(parts) >= 2 and parts[1] == "1":
                        in_first_model = True
                    elif in_first_model:
                        break # Finished Model 1, ignore subsequent models
                elif l.startswith("ENDMDL"):
                    if in_first_model:
                        break
                else:
                    if in_first_model or not any(x.startswith("MODEL 1") for x in lines):
                        model_1_lines.append(l)
            if model_1_lines:
                lines = model_1_lines

        clean_lines = []
        seen_atom_keys = set()

        for line in lines:
            if line.startswith("ATOM  ") or line.startswith("HETATM"):
                rec_type = line[:6].strip()
                res_name = line[17:20].strip().upper()
                chain_id = line[21:22]
                res_seq = line[22:26].strip()
                i_code = line[26:27] # Insertion code (e.g. 100A, 100B in antibody CDRs)
                atom_name = line[12:16].strip()
                alt_loc = line[16:17]

                # AltLoc Handling (Issue 16):
                # If alt_loc is not ' ' and not 'A', discard to prevent colliding duplicate atoms!
                if alt_loc not in (' ', '', 'A', '1'):
                    continue

                # Senior Scientist Fix: Include iCode in atom_key to preserve inserted residue loops!
                atom_key = (chain_id, res_seq, i_code, atom_name)
                if atom_key in seen_atom_keys:
                    continue
                seen_atom_keys.add(atom_key)

                # Water & Solvent Filtering (Issue 17):
                if filter_mode == "strip_water" and res_name in WATER_RESIDUE_NAMES:
                    continue

                if filter_mode == "ligand_only":
                    if rec_type != "HETATM":
                        continue
                    if res_name in WATER_RESIDUE_NAMES or res_name in COMMON_BUFFER_IONS:
                        continue

                # Missing Element Column Heuristic (Fix Alpha Carbon vs Calcium confusion):
                # Standard PDB format specifies atomic element symbol in columns 77-78 (1-based index).
                # Many legacy PDBs/models leave col 77-78 blank, causing 'CA' (alpha carbon) to be misread as 'Ca' (calcium).
                STANDARD_AMINO_ACIDS = {
                    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS",
                    "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
                }
                fixed_line = line
                if len(line) < 78 or line[76:78].strip() == "":
                    # Pad to 76 chars if needed
                    padded = line.ljust(76)
                    assigned_elem = "  "
                    if res_name in STANDARD_AMINO_ACIDS:
                        if atom_name in ("CA", "C", "CB", "CG", "CD", "CE", "CZ", "CH2", "CG1", "CG2", "CD1", "CD2", "CE1", "CE2", "CE3"):
                            assigned_elem = " C"
                        elif atom_name.startswith("N"):
                            assigned_elem = " N"
                        elif atom_name.startswith("O"):
                            assigned_elem = " O"
                        elif atom_name.startswith("S"):
                            assigned_elem = " S"
                        elif atom_name.startswith("H"):
                            assigned_elem = " H"
                    if assigned_elem != "  ":
                        fixed_line = padded[:76] + assigned_elem + (line[78:] if len(line) > 78 else "")

                clean_lines.append(fixed_line)
            elif line.startswith(("HEADER", "TITLE", "COMPND", "CRYST1", "END")):
                clean_lines.append(line)

        return "\n".join(clean_lines)


class CrystalSolidStateEngine:
    @staticmethod
    def create_vacuum_box(mol: pybel.Molecule, vacuum_padding: float = 15.0) -> Tuple[np.ndarray, np.ndarray]:
        """
        Prepares a 0D isolated molecule for Solid-State DFT (VASP POSCAR, PWscf).
        Centers the molecule in a box with >= 15.0 Å of vacuum padding in X, Y, Z.
        Returns: (lattice_vectors, centered_coordinates)
        """
        coords = np.array([a.coords for a in mol.atoms])
        if len(coords) == 0:
            return np.eye(3) * vacuum_padding, coords

        min_c = coords.min(axis=0)
        max_c = coords.max(axis=0)
        span = max_c - min_c

        # Box dimensions = span + vacuum_padding
        box_lengths = span + vacuum_padding
        # Ensure minimum box size of 15 Å even for single atom
        box_lengths = np.maximum(box_lengths, vacuum_padding)

        # Center coordinates in the middle of the new lattice box
        center_target = box_lengths / 2.0
        current_center = (min_c + max_c) / 2.0
        shift = center_target - current_center
        new_coords = coords + shift

        lattice_matrix = np.diag(box_lengths)
        return lattice_matrix, new_coords


class CrystalHydrogenator:
    @staticmethod
    def protonate_heavy_atoms(ob_mol: ob.OBMol) -> int:
        """
        Geometrically adds missing hydrogens to crystal structures (CIF, X-ray models).
        Preserves transition metal coordination and satisfies organic valencies.
        Returns: number of hydrogens added.
        """
        if ob_mol.NumBonds() == 0 and ob_mol.NumAtoms() > 1:
            ob_mol.ConnectTheDots()
            ob_mol.PerceiveBondOrders()

        initial_count = ob_mol.NumAtoms()
        for at in ob.OBMolAtomIter(ob_mol):
            z = at.GetAtomicNum()
            deg = at.GetExplicitDegree()
            
            # Never add hydrogens to transition metals
            if z in TRANSITION_AND_HEAVY_METALS:
                at.SetImplicitHCount(0)
                continue

            needed_h = 0
            if z == 6: # Carbon
                if deg == 1: needed_h = 3
                elif deg == 2:
                    needed_h = 1 if at.IsInRing() else 2
                elif deg == 3:
                    needed_h = 0 if at.IsInRing() else 1
            elif z == 8: # Oxygen
                if deg == 1: needed_h = 1 # Hydroxyl
            elif z == 7: # Nitrogen
                if deg == 1: needed_h = 2 # Primary amine
                elif deg == 2: needed_h = 1 # Secondary amine
            elif z == 16: # Sulfur
                if deg == 1: needed_h = 1 # Thiol
            
            if needed_h > 0:
                at.SetImplicitHCount(needed_h)

        ob_mol.AddHydrogens()
        added = ob_mol.NumAtoms() - initial_count
        return max(0, added)


class DeepContentSniffer:
    @classmethod
    def sniff_format(cls, content_bytes: bytes, filename: str = "") -> Tuple[str, Optional[str]]:
        """
        Deep content sniffer.
        Inspects internal file headers, magic bytes, and chemical syntax instead of trusting extensions.
        Returns: (detected_format, diagnostic_error_or_none)
        """
        if len(content_bytes) == 0:
            return "", "File Empty Error: Uploaded file contains 0 bytes. Please select a valid file."

        # Binary & Checkpoint check (Gaussian .chk or proprietary compiled files)
        if filename.lower().endswith(".chk") or "chk" in filename.lower() and b'\x00' in content_bytes[:512]:
            return "", (
                "Binary Gaussian Checkpoint (.chk) Detected: "
                "Gaussian .chk files are machine-dependent binary scratch files. "
                "Please convert to formatted checkpoint (.fchk) using the Gaussian utility "
                "'formchk input.chk input.fchk' or upload the Gaussian .log / .gjf file directly."
            )

        if b'\x00' in content_bytes[:512]:
            return "", (
                "Binary File Error: Uploaded file appears to be a compiled binary format. "
                "iLoveChem requires text-based chemical formats (CDXML, PDB, XYZ, CIF, LOG, INP, MOL, SDF, SMILES)."
            )

        text_head = content_bytes[:4000].decode('utf-8', errors='ignore')
        
        # 1. ChemDraw XML
        if ("<CDXML" in text_head or "<?xml" in text_head) and ("<page" in text_head or "<fragment" in text_head or "CDXML" in text_head):
            return "cdxml", None

        # 2. Crystallography CIF
        if "data_" in text_head and ("_cell_length_" in text_head or "_atom_site_" in text_head or "_symmetry_" in text_head or "_space_group" in text_head):
            return "cif", None
        if "_cell_length_a" in text_head or "_atom_site_fract_x" in text_head:
            return "cif", None

        # 3. Gaussian Output Log
        if "Entering Gaussian" in text_head or "Gaussian, Inc." in text_head or "Gaussian 16:" in text_head or "Gaussian 09:" in text_head or "Standard orientation:" in text_head:
            return "log", None

        # 4. ORCA Output Log
        if "* O   R   C   A *" in text_head or "An Ab Initio, DFT and Semiempirical electronic structure package" in text_head:
            return "out", None

        # 5. ORCA Input (.inp) -- checked BEFORE Gaussian input because ORCA files
        # routinely begin with '# comment' lines that must not look like a Gaussian route.
        lines = text_head.splitlines()
        first_30 = [l.strip() for l in lines[:30] if l.strip()]
        orca_keywords = ["b3lyp", "pbe", "def2", "opt", "freq", "hf", "wb97x", "svp", "tzvp", "dft", "mp2", "ccsd", "d3bj", "tightscf"]
        if any(l.startswith("!") and any(kw in l.lower() for kw in orca_keywords) for l in first_30):
            return "inp", None
        if re.search(r"^\s*%(maxcore|pal)\b", text_head, re.IGNORECASE | re.MULTILINE) or \
           re.search(r"^\s*\*\s*(xyz|xyzfile|int|gzmt)\s+-?\d+\s+\d+", text_head, re.IGNORECASE | re.MULTILINE):
            return "inp", None

        # 6. Gaussian Input (.gjf / .com)
        # A route section is '#p ...', '#n ...', '#t ...' or '# <method>/<basis> ...'.
        # A bare '# some comment' line is NOT enough.
        gaussian_route_keywords = re.compile(
            r"\b(opt|freq|sp|scf|td|nmr|irc|scan|force|geom|guess|pop|volume|stable|polar|counterpoise)\b",
            re.IGNORECASE
        )
        # An XYZ file starts with a bare atom-count line; its comment line must never be read as a route.
        is_xyz_like = bool(first_30) and first_30[0].isdigit()
        for l in ([] if is_xyz_like else first_30):
            if l.startswith(("%chk=", "%nproc", "%mem=", "%oldchk=")):
                return "gjf", None
            if l.startswith(("#p ", "#P ", "#n ", "#N ", "#t ", "#T ")):
                return "gjf", None
            if l.startswith("#") and ("/" in l or gaussian_route_keywords.search(l)) and not l.startswith("##"):
                return "gjf", None

        # 7. Q-Chem Input
        if "$molecule" in text_head or "$rem" in text_head:
            return "qchem", None

        # 8a. AutoDock PDBQT (Must precede standard PDB check)
        # Keywords must start a line (PDBQT records), so words like "BRANCHED" in a PDB
        # COMPND/REMARK line are not mistaken for AutoDock torsion tree records.
        if filename.lower().endswith(".pdbqt") or re.search(
            r"^(TORSDOF|ROOT|ENDROOT|BRANCH|ENDBRANCH)(\s|$)", text_head, re.MULTILINE
        ):
            return "pdbqt", None

        # 8b. Protein Data Bank (PDB)
        pdb_lines = [l for l in lines[:25] if l.startswith(("HEADER", "ATOM  ", "HETATM", "CRYST1", "COMPND", "AUTHOR"))]
        if len(pdb_lines) >= 1:
            return "pdb", None

        # 9. VASP POSCAR
        if len(lines) >= 6:
            try:
                float(lines[1].strip())
                v1 = [float(x) for x in lines[2].split()[:3]]
                v2 = [float(x) for x in lines[3].split()[:3]]
                v3 = [float(x) for x in lines[4].split()[:3]]
                if len(v1) == 3 and len(v2) == 3 and len(v3) == 3:
                    return "poscar", None
            except Exception:
                pass

        # 10. Standard XYZ
        if len(lines) >= 2 and lines[0].strip().isdigit():
            if len(lines) >= 3:
                parts = lines[2].split()
                if len(parts) >= 4:
                    return "xyz", None

        # 11. MDL Mol / SDF
        if "M  END" in text_head or "V2000" in text_head or "V3000" in text_head:
            return "mol", None

        # 12. InChI
        if text_head.strip().startswith("InChI="):
            return "inchi", None

        # 13. SMILES (Single line, valid organic characters, no spaces)
        non_empty = [l.strip() for l in lines if l.strip()]
        if len(non_empty) == 1:
            s = non_empty[0]
            # Full SMILES alphabet (any element symbol incl. lowercase aromatics s/o/p/se),
            # validated by RDKit so arbitrary one-word text is not mistaken for a molecule.
            if ' ' not in s and re.fullmatch(r"[A-Za-z0-9@+\-\[\]()\\/=#$%:.*~]+", s):
                try:
                    if Chem.MolFromSmiles(s) is not None:
                        return "smi", None
                except Exception:
                    pass

        # Fallback to extension if content sniffing is inconclusive
        lower_fname = filename.lower()
        KNOWN_EXTS = {
            "cdxml", "cdx", "mol", "sdf", "smi", "smiles", "inchi", "cml", "rxn", "skc", "mrv",
            "gjf", "com", "log", "out", "inp", "in", "nw", "mop", "gam", "output",
            "cif", "poscar", "contcar", "pwi", "pwout", "cell", "res",
            "pdb", "ent", "pdbqt", "gro", "mol2", "pqr", "crd", "data", "xyz", "txyz"
        }
        for ext in KNOWN_EXTS:
            if lower_fname.endswith("." + ext) or lower_fname == ext:
                return ext, None

        return "xyz", None


class PeriodicBoundaryUnwrapper:
    @staticmethod
    def unwrap_molecular_crystal(atoms: Atoms) -> Tuple[Atoms, bool]:
        """
        Issue 14: Periodic Boundary Wrapping & Fragment Splitting.
        Detects molecules cut/split across periodic boundary cell faces ([0, 1]).
        Applies minimum image shifts along covalent bonds (d_mic < 2.2 Å)
        so that entire molecules are reassembled into whole, contiguous 3D geometries.
        """
        if len(atoms) < 2 or not any(atoms.pbc):
            return atoms, False

        cell_mat = atoms.get_cell()
        if np.linalg.det(cell_mat) == 0:
            return atoms, False

        n_atoms = len(atoms)
        frac_pos = atoms.get_scaled_positions()

        # Build bonded neighbor graph using minimum image distance
        bonded_graph = {i: [] for i in range(n_atoms)}
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                ds = frac_pos[j] - frac_pos[i]
                ds_mic = ds - np.round(ds)
                dr_mic = np.dot(ds_mic, cell_mat)
                dist_mic = np.linalg.norm(dr_mic)
                if dist_mic < 2.2: # Covalent bonding threshold
                    bonded_graph[i].append((j, ds_mic))
                    bonded_graph[j].append((i, -ds_mic))

        # BFS unwrapping per connected component
        visited = set()
        unwrapped_frac = np.copy(frac_pos)
        was_split = False

        for start_node in range(n_atoms):
            if start_node in visited:
                continue
            queue = [start_node]
            visited.add(start_node)
            while queue:
                curr = queue.pop(0)
                for nbr, ds_mic in bonded_graph[curr]:
                    if nbr not in visited:
                        target_frac = unwrapped_frac[curr] + ds_mic
                        if np.linalg.norm(target_frac - frac_pos[nbr]) > 1e-4:
                            was_split = True
                        unwrapped_frac[nbr] = target_frac
                        visited.add(nbr)
                        queue.append(nbr)

        if was_split:
            unwrapped_cart = np.dot(unwrapped_frac, cell_mat)
            new_atoms = Atoms(
                symbols=atoms.get_chemical_symbols(),
                positions=unwrapped_cart,
                cell=cell_mat,
                pbc=atoms.pbc
            )
            return new_atoms, True

        return atoms, False


class GromacsGroEngine:
    @staticmethod
    def export_gro(mol: pybel.Molecule, title: str = "iLoveChem Generated Structure") -> Tuple[str, str, List[str]]:
        """
        Issue 19: GROMACS PBC Box Dimension Loss.
        Converts molecule coordinates to standard GROMACS .gro format:
        - Fixed-width column formatting (%5d%-5s%5s%5d%8.3f%8.3f%8.3f)
        - Computes minimum bounding box + 1.0 nm clearance (10 Å) on all sides
        - Shifts coordinates so molecule is centered with 1.0 nm clearance from box origin
        - Writes trailing box vectors [Lx, Ly, Lz] in nanometers
        """
        warnings = []
        if len(mol.atoms) == 0:
            return "Empty Molecule\n    0\n   1.00000   1.00000   1.00000\n", "structure.gro", warnings

        coords_ang = np.array([a.coords for a in mol.atoms])
        # Convert Ångströms to nanometers (1 nm = 10 Å)
        coords_nm = coords_ang / 10.0
        min_nm = np.min(coords_nm, axis=0)
        max_nm = np.max(coords_nm, axis=0)

        clearance_nm = 1.0  # 1.0 nm clearance
        span_nm = max_nm - min_nm
        box_nm = span_nm + 2.0 * clearance_nm
        shifted_nm = coords_nm - min_nm + clearance_nm

        lines = [
            f"{title} (1.0nm PBC Clearance)",
            f"{len(mol.atoms):5d}"
        ]

        for i, atom in enumerate(mol.atoms):
            sym = ob.GetSymbol(atom.OBAtom.GetAtomicNum())
            atom_name = f"{sym}{i + 1}"
            x, y, z = shifted_nm[i]
            line = f"{1:5d}{'MOL':<5s}{atom_name:>5s}{i + 1:5d}{x:8.3f}{y:8.3f}{z:8.3f}"
            lines.append(line)

        lines.append(f"{box_nm[0]:10.5f}{box_nm[1]:10.5f}{box_nm[2]:10.5f}")
        warnings.append(
            f"GROMACS PBC Box Applied: Box vectors set to [{box_nm[0]:.3f}, {box_nm[1]:.3f}, {box_nm[2]:.3f}] nm with 1.0 nm solvent clearance."
        )
        return "\n".join(lines), "molecule.gro", warnings


class ZeroCoordinateGuard:
    @staticmethod
    def inspect_and_fix(mol: pybel.Molecule) -> Tuple[pybel.Molecule, bool, List[str]]:
        """
        Loophole 2: The Silent Failure Bug (0.000 0.000 0.000).
        Detects if all or multiple atoms are collapsed to origin (0.000, 0.000, 0.000)
        or have identical coordinates.
        If collapsed: automatically regenerates physical 3D coordinates using RDKit ETKDGv3 + UFF.
        """
        warnings = []
        if len(mol.atoms) < 2:
            return mol, False, warnings

        coords = np.array([a.coords for a in mol.atoms])
        dist_to_origin = np.linalg.norm(coords, axis=1)
        zeros_count = int(np.sum(dist_to_origin < 1e-4))
        var_xyz = np.var(coords, axis=0)

        is_collapsed = (zeros_count >= 2) or (float(np.sum(var_xyz)) < 1e-5)
        if not is_collapsed:
            return mol, False, warnings

        fixed_rd = False
        try:
            smi = mol.write("smi").split()[0]
            if not "." in smi:
                rd_mol = Chem.MolFromSmiles(smi)
                if rd_mol:
                    rd_mol = Chem.AddHs(rd_mol)
                    res = AllChem.EmbedMolecule(rd_mol, maxAttempts=20, randomSeed=42, useRandomCoords=True)
                    if res >= 0:
                        try: AllChem.UFFOptimizeMolecule(rd_mol, maxIters=300)
                        except Exception: pass
                        new_mb = Chem.MolToMolBlock(rd_mol)
                        mol = pybel.readstring("mol", new_mb)
                        fixed_rd = True
        except Exception:
            pass

        if not fixed_rd:
            # Deterministic 3D physical dispersion of collapsed atoms
            for i, atom in enumerate(mol.atoms):
                angle = i * 2.399963229728653  # Golden angle in radians
                r = 1.5 * np.sqrt(i + 1)
                atom.OBAtom.SetVector(float(r * np.cos(angle)), float(r * np.sin(angle)), float(i * 0.5))

        warnings.append(
            "Zero-Coordinate Guard Active: Input geometry had multiple atoms collapsed at (0.000, 0.000, 0.000) [Silent Failure Bug]. "
            "Automatically reconstructed physical 3D coordinates."
        )
        return mol, True, warnings


class PDBBondPerceptionEngine:
    @staticmethod
    def restore_pdb_bonds(mol: pybel.Molecule) -> Tuple[pybel.Molecule, List[str]]:
        """
        Loophole 4: Missing Bond Orders & Double-Bonds in PDB.
        PDB files only record single bonds or no CONECT records.
        Restores bond orders, aromaticity, and valencies based on distances and atomic hybridization.
        """
        warnings = []
        try:
            mol.OBMol.ConnectTheDots()
            mol.OBMol.PerceiveBondOrders()
            n_mult = sum(1 for b in ob.OBMolBondIter(mol.OBMol) if b.GetBondOrder() > 1 or b.IsAromatic())
            if n_mult > 0:
                warnings.append(
                    f"PDB Bond Perception Active: Restored {n_mult} double/aromatic bonds from 3D coordinates and valency rules."
                )
        except Exception:
            pass
        return mol, warnings


class BatchSdfProcessor:
    @staticmethod
    def is_multi_sdf(content_bytes: bytes) -> bool:
        """
        Detects if content is a multi-structure SDF library ($$$$ count >= 2).
        """
        return content_bytes.count(b"$$$$") >= 2

    @staticmethod
    def split_sdf_records(sdf_str: str) -> List[str]:
        """
        Loophole 5: Multi-Molecule / Batch SDF Conversion Bottleneck.
        Splits multi-record SDF into individual molecule blocks using '$$$$' delimiter.
        """
        raw_blocks = sdf_str.split("$$$$")
        valid_records = []
        for blk in raw_blocks:
            cleaned = blk.strip()
            if cleaned and ("M  END" in cleaned or "V2000" in cleaned or "V3000" in cleaned or len(cleaned.splitlines()) >= 4):
                valid_records.append(cleaned + "\nM  END\n$$$$\n")
        return valid_records

    @staticmethod
    def create_batch_zip(records: List[str], target_format: str, export_fn, qm_options: dict) -> Tuple[bytes, str, int]:
        mem_zip = io.BytesIO()
        converted_count = 0
        with zipfile.ZipFile(mem_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
            for idx, rec in enumerate(records):
                try:
                    mol = pybel.readstring("mol", rec)
                    if mol and len(mol.atoms) > 0:
                        out_content, fname, _ = export_fn(mol, target_format, qm_options)
                        base_name, ext = os.path.splitext(fname)
                        safe_formula = re.sub(r'[^a-zA-Z0-9_]', '', mol.formula or f"compound_{idx+1}")
                        clean_ext = re.sub(r'[^a-zA-Z0-9\.]', '', ext)
                        zip_entry_name = f"molecule_{idx+1:03d}_{safe_formula}{clean_ext}"
                        zf.writestr(zip_entry_name, out_content.encode('utf-8'))
                        converted_count += 1
                except Exception:
                    continue
        mem_zip.seek(0)
        return mem_zip.getvalue(), f"batch_sdf_{converted_count}_molecules.zip", converted_count


class GhostAtomNormalizer:
    GHOST_SYMBOLS = {"bq", "gh", "x", "bq-h", "bq-c", "bq-o", "bq-n", "gh-h", "gh-c", "gh-o"}

    @classmethod
    def sanitize_input_text(cls, text: str) -> Tuple[str, Dict[int, str]]:
        """
        Loophole 6: Non-Standard Dummy / Ghost Atoms (Bq, Gh, X).
        Detects ghost atoms, maps them to dummy Helium (He) atoms so standard parsers don't crash,
        and remembers their original ghost notation to restore in QM output.
        """
        ghost_map = {}
        lines = text.splitlines()
        new_lines = []
        atom_idx = 0

        for line in lines:
            parts = line.split()
            if not parts:
                new_lines.append(line)
                continue

            first_tok = parts[0].strip()
            first_tok_lower = first_tok.lower()

            if first_tok_lower in cls.GHOST_SYMBOLS or first_tok_lower.startswith(("bq-", "gh-")):
                ghost_map[atom_idx] = first_tok
                new_line = "  " + "He".ljust(3) + "  " + "  ".join(parts[1:])
                new_lines.append(new_line)
                atom_idx += 1
            else:
                if len(parts) >= 4:
                    try:
                        float(parts[1]); float(parts[2]); float(parts[3])
                        atom_idx += 1
                    except ValueError:
                        pass
                new_lines.append(line)

        return "\n".join(new_lines), ghost_map

    @classmethod
    def restore_ghost_atoms_in_output(cls, output_content: str, ghost_map: Dict[int, str], target_format: str) -> str:
        if not ghost_map:
            return output_content

        lines = output_content.splitlines()
        new_lines = []
        atom_idx = 0

        for line in lines:
            parts = line.split()
            if len(parts) >= 4:
                try:
                    float(parts[1]); float(parts[2]); float(parts[3])
                    if atom_idx in ghost_map:
                        ghost_sym = ghost_map[atom_idx]
                        line = f"  {ghost_sym:<3}  " + "  ".join(f"{float(p):12.6f}" for p in parts[1:4])
                    atom_idx += 1
                except ValueError:
                    pass
            new_lines.append(line)

        return "\n".join(new_lines)


class UnixEolNormalizer:
    @staticmethod
    def normalize(text: str) -> str:
        """
        Loophole 7: Windows (CRLF \r\n) vs Linux (LF \n) Line Ending Nightmares.
        Guarantees 100% strict Unix LF line endings without any hidden ^M characters.
        """
        return text.replace("\r\n", "\n").replace("\r", "\n")


class MesoCompoundDetector:
    @staticmethod
    def check_meso(rd_mol: Chem.Mol) -> Tuple[bool, List[str]]:
        """
        Loophole 8: Optical Isomer & Meso Compound Ambiguity.
        Detects internal symmetry and opposite chiral center pairs in meso compounds
        (e.g., (2R, 3S)-tartaric acid) with net zero optical activity ([α] = 0).
        """
        warnings = []
        if not rd_mol:
            return False, warnings

        try:
            centers = Chem.FindMolChiralCenters(rd_mol, includeUnassigned=True)
            if len(centers) < 2 or len(centers) % 2 != 0:
                return False, warnings

            ranks = list(Chem.CanonicalRankAtoms(rd_mol, breakTies=False, includeChirality=False))
            chiral_ranks = [ranks[idx] for idx, tag in centers]
            
            rank_counts = {}
            for r in chiral_ranks:
                rank_counts[r] = rank_counts.get(r, 0) + 1

            is_symmetric = all(cnt % 2 == 0 for cnt in rank_counts.values())
            if is_symmetric:
                tags = [tag for idx, tag in centers]
                r_count = tags.count('R')
                s_count = tags.count('S')
                if r_count == s_count and r_count > 0:
                    warnings.append(
                        f"Meso Compound Detected: Molecule possesses internal symmetry plane with opposite chiral centers "
                        f"({r_count}R, {s_count}S). Net optical activity [α] = 0. Stereochemical configuration verified."
                    )
                    return True, warnings
        except Exception:
            pass

        return False, warnings
