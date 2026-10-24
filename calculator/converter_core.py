import os
import re
import io
from typing import Dict, Any, Tuple, Optional, List
import numpy as np

from openbabel import openbabel as ob
from openbabel import pybel
import cclib
from ase import io as ase_io
from ase.atoms import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem

from .format_registry import COMPLETE_TARGET_FORMATS
from .phase1_engine import (
    sanitize_and_extract_all_fragments,
    build_3d_mol_from_fragment,
    create_zip_archive_of_all_fragments
)
from .scientific_engine import (
    SpinChargeParityValidator,
    QuantumLogParser,
    PhysicalIntegrityChecker,
    PDBBioSanitizer,
    CrystalSolidStateEngine,
    TRANSITION_AND_HEAVY_METALS,
    SaltAndSolventStripper,
    CrystalHydrogenator,
    DeepContentSniffer,
    PeriodicBoundaryUnwrapper,
    GromacsGroEngine,
    ZeroCoordinateGuard,
    PDBBondPerceptionEngine,
    BatchSdfProcessor,
    GhostAtomNormalizer,
    UnixEolNormalizer,
    MesoCompoundDetector
)

KNOWN_INPUT_EXTENSIONS = {
    "cdxml", "cdx", "mol", "sdf", "smi", "smiles", "inchi", "cml", "rxn", "skc", "mrv",
    "gjf", "com", "log", "out", "inp", "in", "nw", "mop", "gam", "output",
    "cif", "poscar", "contcar", "pwi", "pwout", "cell", "res",
    "pdb", "ent", "pdbqt", "gro", "mol2", "pqr", "crd", "data", "xyz", "txyz"
}

def detect_format(filename: str, content_str: str) -> str:
    """
    Deep content sniffer: Inspects chemical headers and signatures first,
    falling back to extension only if content is ambiguous.
    """
    fmt, _ = DeepContentSniffer.sniff_format(content_str.encode('utf-8', errors='ignore'), filename)
    return fmt or "xyz"


def parse_and_extract_atoms(
    input_bytes: bytes,
    in_format: str,
    filename: str = "",
    fragment_index: int = 0,
    pdb_filter_mode: str = "all",
    add_hydrogens: bool = True,
    strip_salts: bool = True
) -> Tuple[Optional[pybel.Molecule], Optional[dict]]:
    """
    Universal chemical reader with automated scientific guardrails:
    - Phase 1: ChemDraw Reaction schemes & multi-fragments
    - Phase 2: Inverted output log parser (Gaussian / ORCA stationary points)
    - Phase 3: Crystallography (CIF / POSCAR) & fractional handling
    - Phase 4: PDB AltLoc deduplication & water/ligand isolation
    - Phase 5: Heavy metal UFF fallback & nuclear clash auditing
    """
    meta_info = {
        "is_2d": False,
        "is_periodic": False,
        "parsed_by": "",
        "formula": "",
        "num_atoms": 0,
        "input_format": in_format,
        "multiple_fragments": False,
        "fragment_count": 1,
        "fragments_list": [],
        "filtered_elements": {},
        "warnings": [],
        "has_metals": False,
        "metal_symbols": [],
        "has_clashes": False,
        "min_interatomic_distance": 0.0,
        "salt_detected": False,
        "salts_stripped": False,
        "detected_salts": [],
        "charge": 0
    }

    # 1. ChemDraw XML Direct Engine with Reaction Scheme Handling
    if in_format == "cdxml" or b"CDXML" in input_bytes[:1000]:
        try:
            cdxml_str = input_bytes.decode('utf-8', errors='ignore')
            fragments, stats = sanitize_and_extract_all_fragments(cdxml_str)
            meta_info["filtered_elements"] = stats
            
            if len(fragments) > 0:
                meta_info["is_2d"] = True
                meta_info["multiple_fragments"] = (len(fragments) > 1)
                meta_info["fragment_count"] = len(fragments)
                meta_info["fragments_list"] = [
                    {
                        "index": f["index"],
                        "formula": f["formula"],
                        "atom_count": f["atom_count"],
                        "is_largest": f["is_largest"]
                    }
                    for f in fragments
                ]
                
                # Pick requested fragment (or default to largest product)
                target_frag = fragments[0]
                if fragment_index < len(fragments):
                    target_frag = fragments[fragment_index]
                else:
                    target_frag = sorted(fragments, key=lambda x: x['atom_count'], reverse=True)[0]

                py_mol = build_3d_mol_from_fragment(target_frag)
                if py_mol:
                    meta_info["parsed_by"] = "ChemDraw Reaction Engine + 3D Conformer"
                    meta_info["formula"] = py_mol.formula
                    meta_info["num_atoms"] = len(py_mol.atoms)
                    py_mol = _apply_integrity_audit(py_mol, meta_info, strip_salts=strip_salts)
                    return py_mol, meta_info
        except Exception as e:
            meta_info["warnings"].append(f"ChemDraw parser note: {str(e)}")

    # 1b. InChI Native Support via RDKit
    if in_format in ["inchi"] or (isinstance(input_bytes, (bytes, bytearray)) and input_bytes.strip().startswith(b"InChI=")):
        try:
            inchi_str = input_bytes.decode('utf-8', errors='ignore').strip()
            rd_mol = Chem.MolFromInchi(inchi_str)
            if rd_mol:
                if add_hydrogens:
                    rd_mol = Chem.AddHs(rd_mol)
                AllChem.EmbedMolecule(rd_mol, AllChem.ETKDGv3())
                AllChem.MMFFOptimizeMolecule(rd_mol)
                molblock = Chem.MolToMolBlock(rd_mol)
                py_mol = pybel.readstring("mol", molblock)
                meta_info["parsed_by"] = "RDKit InChI Engine (3D Conformer)"
                meta_info["formula"] = py_mol.formula
                meta_info["num_atoms"] = len(py_mol.atoms)
                py_mol = _apply_integrity_audit(py_mol, meta_info, strip_salts=strip_salts)
                return py_mol, meta_info
        except Exception as e:
            meta_info["warnings"].append(f"InChI parser note: {str(e)}")

    # 2. Output logs via Inverted QuantumLogParser (Issue 7) or cclib
    if in_format in ["log", "out", "output", "pwout"]:
        try:
            content_str = input_bytes.decode('utf-8', errors='ignore')
            
            # Step 2a: Try reverse-scan parser for exact stationary/converged geometry
            scan_res = QuantumLogParser.reverse_scan_gaussian_or_orca_log(content_str)
            if scan_res:
                rev_coords, parser_label = scan_res
                ob_mol = ob.OBMol()
                for sym, x, y, z in rev_coords:
                    at = ob_mol.NewAtom()
                    at.SetAtomicNum(ob.GetAtomicNum(sym))
                    at.SetVector(x, y, z)
                py_mol = pybel.Molecule(ob_mol)
                meta_info["parsed_by"] = f"Streaming Tail Parser ({parser_label})"
                meta_info["formula"] = py_mol.formula
                meta_info["num_atoms"] = len(py_mol.atoms)
                py_mol = _apply_integrity_audit(py_mol, meta_info, strip_salts=strip_salts)
                return py_mol, meta_info

            # Step 2b: Fallback to cclib, fed from an in-memory stream so the upload is
            # never written to disk (zero-disk invariant). A cclib failure must not
            # skip the rescue parser below, so it is isolated in its own try block.
            data = None
            try:
                parser = cclib.io.ccopen(io.StringIO(content_str))
                data = parser.parse()
            except Exception as cclib_err:
                meta_info["warnings"].append(f"cclib could not parse this log ({cclib_err}); trying rescue parser.")

            if data is not None and hasattr(data, 'atomcoords') and len(data.atomcoords) > 0:
                final_coords = data.atomcoords[-1]
                atom_nos = data.atomnos
                
                ob_mol = ob.OBMol()
                for i, at_num in enumerate(atom_nos):
                    at = ob_mol.NewAtom()
                    at.SetAtomicNum(int(at_num))
                    coord = final_coords[i]
                    at.SetVector(float(coord[0]), float(coord[1]), float(coord[2]))
                
                py_mol = pybel.Molecule(ob_mol)
                pkg = getattr(data, 'metadata', {}).get('package', 'Quantum Log')
                meta_info["parsed_by"] = f"cclib ({pkg}) [Final Geometry]"
                meta_info["formula"] = py_mol.formula
                meta_info["num_atoms"] = len(py_mol.atoms)
                py_mol = _apply_integrity_audit(py_mol, meta_info, strip_salts=strip_salts)
                return py_mol, meta_info

            # Step 2c: Loophole 3 — Rescue Incomplete/Failed Quantum Optimizations
            rescue_res = QuantumLogParser.rescue_abnormal_optimization(content_str)
            if rescue_res:
                res_coords, parser_label, min_energy, step_num = rescue_res
                ob_mol = ob.OBMol()
                for sym, x, y, z in res_coords:
                    at = ob_mol.NewAtom()
                    at.SetAtomicNum(ob.GetAtomicNum(sym))
                    at.SetVector(x, y, z)
                py_mol = pybel.Molecule(ob_mol)
                meta_info["parsed_by"] = parser_label
                meta_info["formula"] = py_mol.formula
                meta_info["num_atoms"] = len(py_mol.atoms)
                meta_info["warnings"].append(
                    f"Abnormal Termination Rescued: Calculation terminated before complete convergence. "
                    f"Extracted lowest-energy intermediate geometry (Step #{step_num}, E = {min_energy:.6f} Har) so computing time is not lost."
                )
                py_mol = _apply_integrity_audit(py_mol, meta_info, strip_salts=strip_salts)
                return py_mol, meta_info
        except Exception as e:
            meta_info["warnings"].append(f"Quantum log parser note: {str(e)}")

    # 3a. AutoDock PDBQT Ligand Support
    if in_format in ["pdbqt"]:
        try:
            content_str = input_bytes.decode('utf-8', errors='ignore')
            mol = pybel.readstring("pdbqt", content_str)
            meta_info["parsed_by"] = "AutoDock PDBQT Ligand Engine"
            meta_info["formula"] = mol.formula
            meta_info["num_atoms"] = len(mol.atoms)
            mol = _apply_integrity_audit(mol, meta_info, strip_salts=strip_salts)
            return mol, meta_info
        except Exception as e:
            meta_info["warnings"].append(f"PDBQT parser note: {str(e)}")

    # 3b. PDB / Macromolecules Sanitization (Issues 16 & 17 + Loophole 4)
    if in_format in ["pdb", "ent"]:
        try:
            raw_pdb = input_bytes.decode('utf-8', errors='ignore')
            sanitized_pdb = PDBBioSanitizer.sanitize_pdb(raw_pdb, filter_mode=pdb_filter_mode)
            mol = pybel.readstring("pdb", sanitized_pdb)
            # Loophole 4: Restore double and aromatic bonds in PDB
            mol, pdb_bond_warns = PDBBondPerceptionEngine.restore_pdb_bonds(mol)
            if pdb_bond_warns:
                meta_info["warnings"].extend(pdb_bond_warns)
            if add_hydrogens:
                mol.addh()
            meta_info["parsed_by"] = f"PDB Bio-Sanitizer (Mode: {pdb_filter_mode})"
            meta_info["formula"] = mol.formula
            meta_info["num_atoms"] = len(mol.atoms)
            mol = _apply_integrity_audit(mol, meta_info, strip_salts=strip_salts)
            return mol, meta_info
        except Exception as e:
            meta_info["warnings"].append(f"PDB sanitizer note: {str(e)}")

    # 4. Materials & Crystals via ASE (Issues 12, 13, 14)
    if in_format in ["cif", "poscar", "contcar", "pwi", "cell", "res"]:
        try:
            content_str = input_bytes.decode('utf-8', errors='ignore')
            f_in = io.StringIO(content_str)
            fmt_map = {"cif": "cif", "poscar": "vasp", "contcar": "vasp", "pwi": "espresso-in"}
            atoms = ase_io.read(f_in, format=fmt_map.get(in_format, "cif"))

            # Issue 14: Periodic Boundary Wrapping & Split Molecule Unwrapper
            atoms, was_split = PeriodicBoundaryUnwrapper.unwrap_molecular_crystal(atoms)
            if was_split:
                meta_info["warnings"].append(
                    "Periodic Boundary Unwrapping Active: Reassembled molecules cut across unit cell boundary faces into contiguous 3D chemical structures."
                )
            
            ob_mol = ob.OBMol()
            for at in atoms:
                o_at = ob_mol.NewAtom()
                o_at.SetAtomicNum(int(at.number))
                o_at.SetVector(float(at.position[0]), float(at.position[1]), float(at.position[2]))
            
            has_h = any(at.number == 1 for at in atoms)
            if add_hydrogens and in_format == "cif" and not has_h:
                added = CrystalHydrogenator.protonate_heavy_atoms(ob_mol)
                if added > 0:
                    meta_info["warnings"].append(
                        f"CIF Auto-Hydrogenation Active: Input crystal structure lacked hydrogen coordinates (common in X-ray crystallography). "
                        f"Geometrically added {added} missing hydrogen atoms to complete valence shells."
                    )
            
            py_mol = pybel.Molecule(ob_mol)
            meta_info["parsed_by"] = "ASE Crystallography Engine (Fractional to Cartesian)"
            meta_info["is_periodic"] = True
            meta_info["formula"] = str(atoms.symbols)
            meta_info["num_atoms"] = len(atoms)
            py_mol = _apply_integrity_audit(py_mol, meta_info, strip_salts=strip_salts)
            return py_mol, meta_info
        except Exception as e:
            meta_info["warnings"].append(f"ASE crystal parser note: {str(e)}")

    # 5. Standard OpenBabel Formats with Heavy Metal UFF Detection (Issue 20 & Loophole 6)
    try:
        ob_fmt = in_format
        if in_format in ["gjf", "com"]: ob_fmt = "gjf"
        elif in_format in ["inp"]: ob_fmt = "orca"

        raw_text = input_bytes.decode('utf-8', errors='ignore')
        # Loophole 6: Ghost / Dummy Atom Normalization
        clean_text, ghost_map = GhostAtomNormalizer.sanitize_input_text(raw_text)
        meta_info["ghost_map"] = ghost_map
        if ghost_map:
            meta_info["warnings"].append(
                f"Ghost Atoms Preserved: Detected {len(ghost_map)} BSSE counterpoise ghost/dummy centers ({list(ghost_map.values())[:3]})."
            )

        mol = pybel.readstring(ob_fmt, clean_text)
        
        # Loophole 2: Inspect if input was a 3D format with collapsed (0,0,0) coordinates
        mol, was_fixed, zero_warns = ZeroCoordinateGuard.inspect_and_fix(mol)
        if was_fixed:
            meta_info["warnings"].extend(zero_warns)
            meta_info["silent_failure_prevented"] = True

        has_3d = False
        coords = [a.coords for a in mol.atoms]
        if coords and not was_fixed:
            z_vals = [c[2] for c in coords]
            if any(abs(z) > 1e-4 for z in z_vals):
                has_3d = True
        elif was_fixed:
            has_3d = True

        if not has_3d or in_format in ["cdx", "smi", "smiles", "mol"]:
            meta_info["is_2d"] = True
            
            # Pre-strip salts and counter-ions before 3D conformer embedding
            mol, salt_meta = SaltAndSolventStripper.analyze_and_process(mol, strip_salts=strip_salts)
            meta_info["salt_detected"] = salt_meta["salt_detected"]
            meta_info["detected_salts"] = salt_meta["detected_salts"]
            meta_info["salts_stripped"] = salt_meta["stripped"]
            meta_info["salt_processed"] = True
            if salt_meta["warnings"]:
                meta_info["warnings"].extend(salt_meta["warnings"])
            if salt_meta["stripped"]:
                meta_info["formula"] = salt_meta["parent_formula"]
                meta_info["charge"] = salt_meta["parent_charge"]

            if add_hydrogens:
                mol.addh()
            
            # Heavy metal & hypervalent check for force field selection (Issue 20)
            contains_metal = any(a.OBAtom.GetAtomicNum() in TRANSITION_AND_HEAVY_METALS for a in mol.atoms)
            is_hypervalent = any(
                (a.OBAtom.GetAtomicNum() == 16 and a.OBAtom.GetExplicitValence() > 4) or # Hypervalent S(VI)
                (a.OBAtom.GetAtomicNum() == 15 and a.OBAtom.GetExplicitValence() > 4) or # Hypervalent P(V)
                (a.OBAtom.GetAtomicNum() in [53, 35, 17] and a.OBAtom.GetExplicitValence() > 1) # Hypervalent Halogens (e.g. Periodinane)
                for a in mol.atoms
            )
            forcefield = 'uff' if (contains_metal or is_hypervalent) else 'mmff94'
            try:
                mol.make3D(forcefield=forcefield, steps=300)
            except Exception:
                try: mol.make3D(forcefield='uff', steps=300)
                except Exception: pass

        meta_info["parsed_by"] = "OpenBabel 3.1"
        meta_info["formula"] = mol.formula
        meta_info["num_atoms"] = len(mol.atoms)
        mol = _apply_integrity_audit(mol, meta_info, strip_salts=strip_salts)
        return mol, meta_info

    except Exception:
        try:
            rd_mol = Chem.MolFromMolBlock(input_bytes.decode('utf-8', errors='ignore'))
            if not rd_mol:
                rd_mol = Chem.MolFromSmiles(input_bytes.decode('utf-8', errors='ignore').strip())
            if rd_mol:
                if add_hydrogens:
                    rd_mol = Chem.AddHs(rd_mol)
                AllChem.EmbedMolecule(rd_mol, AllChem.ETKDGv3())
                try:
                    AllChem.MMFFOptimizeMolecule(rd_mol)
                except Exception:
                    AllChem.UFFOptimizeMolecule(rd_mol)
                mb = Chem.MolToMolBlock(rd_mol)
                mol = pybel.readstring("mol", mb)
                meta_info["parsed_by"] = "RDKit ETKDGv3"
                meta_info["formula"] = mol.formula
                meta_info["num_atoms"] = len(mol.atoms)
                mol = _apply_integrity_audit(mol, meta_info, strip_salts=strip_salts)
                return mol, meta_info
        except Exception:
            pass

    return None, None


def _apply_integrity_audit(mol: pybel.Molecule, meta_info: dict, strip_salts: bool = True) -> pybel.Molecule:
    """Audits molecular integrity: salts/counter-ions, transition metals, nuclear overlap clashes, and silent failures."""
    # Loophole 2: Zero-Coordinate Guard (Silent Failure Bug interceptor)
    mol, was_fixed, zero_warns = ZeroCoordinateGuard.inspect_and_fix(mol)
    if was_fixed:
        meta_info["warnings"].extend(zero_warns)
        meta_info["silent_failure_prevented"] = True

    # 1. Salt & Counter-Ion Stripper (if not already handled)
    if not meta_info.get("salt_processed", False):
        mol, salt_meta = SaltAndSolventStripper.analyze_and_process(mol, strip_salts=strip_salts)
        meta_info["salt_detected"] = salt_meta["salt_detected"]
        meta_info["detected_salts"] = salt_meta["detected_salts"]
        meta_info["salts_stripped"] = salt_meta["stripped"]
        meta_info["salt_processed"] = True
        if salt_meta["warnings"]:
            meta_info["warnings"].extend(salt_meta["warnings"])
        if salt_meta["stripped"]:
            meta_info["formula"] = salt_meta["parent_formula"]
            meta_info["num_atoms"] = len(mol.atoms)
            meta_info["charge"] = salt_meta["parent_charge"]

    # 2. Transition Metals & Clashes
    audit = PhysicalIntegrityChecker.check_clashes_and_metals(mol)
    meta_info["has_metals"] = audit["has_metals"]
    meta_info["metal_symbols"] = audit["metal_symbols"]
    meta_info["has_clashes"] = audit["has_clashes"]
    meta_info["min_interatomic_distance"] = audit["min_interatomic_distance"]
    if audit["has_clashes"]:
        meta_info["warnings"].append(
            f"Nuclear Overlap Warning: Detected atom clash with minimum distance {audit['min_interatomic_distance']} Å (< 0.7 Å). "
            f"Clashing pairs: {audit['clash_pairs']}"
        )

    # Loophole 8: Meso Compound Detector & Optical Inactivity Verification
    try:
        smi = mol.write("smi").split()[0]
        rd_mol = Chem.MolFromSmiles(smi)
        if rd_mol:
            is_meso, meso_warns = MesoCompoundDetector.check_meso(rd_mol)
            meta_info["is_meso"] = is_meso
            if meso_warns:
                meta_info["warnings"].extend(meso_warns)
    except Exception:
        pass

    # 3. Valence Electron, Formal Charge & Spin Multiplicity Calculation
    detected_charge = mol.charge
    total_z = sum(a.OBAtom.GetAtomicNum() for a in mol.atoms)
    n_electrons = total_z - detected_charge
    suggested_mult = 2 if (n_electrons % 2 == 1) else 1

    # Special Physical Invariant: Dioxygen (O2) ground state is a Triplet (3Sigma_g-)
    is_o2 = (len(mol.atoms) == 2 and all(a.OBAtom.GetAtomicNum() == 8 for a in mol.atoms) and detected_charge == 0)
    if is_o2:
        suggested_mult = 3
        electronic_state = f"Ground-State Triplet Diradical (O₂, {n_electrons} e⁻)"
    elif suggested_mult == 2:
        electronic_state = f"Radical (Doublet, {n_electrons} e⁻)"
    elif detected_charge < 0:
        electronic_state = f"Anion (Charge {detected_charge}, Singlet, {n_electrons} e⁻)"
    elif detected_charge > 0:
        electronic_state = f"Cation (Charge +{detected_charge}, Singlet, {n_electrons} e⁻)"
    else:
        electronic_state = f"Closed-Shell Neutral (Singlet, {n_electrons} e⁻)"

    meta_info["charge"] = detected_charge
    meta_info["detected_charge"] = detected_charge
    meta_info["total_electrons"] = n_electrons
    meta_info["suggested_mult"] = suggested_mult
    meta_info["electronic_state"] = electronic_state
    return mol


def export_to_target_format(
    mol: pybel.Molecule,
    target_format: str,
    qm_options: Optional[Dict[str, Any]] = None
) -> Tuple[str, str, List[str]]:
    """
    Universal target exporter with scientific guardrails:
    - Spin multiplicity & charge parity validator (Issue 8)
    - Guaranteed trailing blank lines for Gaussian (Issue 9)
    - Strict Cartesian coordinate sanitization (Issue 11)
    - 15 Å vacuum padding for VASP POSCAR (Issue 15)
    - AutoDock PDBQT Gasteiger charges (Issue 18)
    Returns: (output_content, filename, list_of_audit_warnings)
    """
    if qm_options is None:
        qm_options = {}

    def _clean_param(val: Any, default: str) -> str:
        s = str(val if val is not None else default).strip()
        # Security Guard: Strip all newlines, carriage returns, null bytes to prevent QM header injection
        s = re.sub(r'[\r\n\0]', '', s)
        # Allow only safe alphanumeric and standard chemical basis/method characters
        cleaned = re.sub(r'[^\w\s\.\-\/\+\*\(\)\,]', '', s).strip()
        return cleaned[:60] if cleaned else default

    method = _clean_param(qm_options.get("method"), "B3LYP")
    basis = _clean_param(qm_options.get("basis"), "def2-SVP")
    calc_type = _clean_param(qm_options.get("calc_type"), "opt freq")
    charge = int(qm_options.get("charge", 0))
    requested_mult = int(qm_options.get("mult", 1))
    nprocs = max(1, min(128, int(qm_options.get("nprocs", 8))))
    mem_gb = max(1, min(1024, int(qm_options.get("mem_gb", 16))))
    warnings = []

    # Senior Scientist Guard: Heavy Element vs Pople Basis Set Incompatibility (BIOVIA/Schrödinger Standard)
    pople_patterns = ["3-21g", "6-31g", "6-311g", "sto-3g", "sto-6g"]
    if any(p in basis.lower() for p in pople_patterns):
        heavy_atoms = [ob.GetSymbol(a.OBAtom.GetAtomicNum()) for a in mol.atoms if a.OBAtom.GetAtomicNum() > 36]
        if heavy_atoms:
            unique_heavy = sorted(set(heavy_atoms))
            warnings.append(
                f"Basis Set Compatibility Notice: Structure contains heavy element(s) ({', '.join(unique_heavy)}, Z > 36). "
                f"Pople basis set '{basis}' does not support elements beyond Krypton and will fail in Gaussian/ORCA. "
                f"Recommendation: Switch basis set to 'def2-SVP', 'def2-TZVP', or 'LANL2DZ' with effective core potentials (ECP)."
            )

    # Calculate total atomic number for parity validation (Issue 8)
    total_z = sum(atom.OBAtom.GetAtomicNum() for atom in mol.atoms)
    mult, parity_warnings = SpinChargeParityValidator.validate_and_fix(total_z, charge, requested_mult)
    warnings.extend(parity_warnings)

    # Dispersion: only append D3BJ for functionals that need it and support it in ORCA.
    # HF, MP2, coupled cluster, wB97X-D3/D4, B97-3c, r2SCAN-3c, etc. either do not use it
    # or already include dispersion, and ORCA rejects or double-counts it.
    D3BJ_FUNCTIONALS = {"b3lyp", "pbe0", "pbe", "bp86", "blyp", "tpss", "tpssh", "revpbe", "b3pw91", "pw91"}
    orca_dispersion = ""
    if method.lower().replace(" ", "") in D3BJ_FUNCTIONALS:
        orca_dispersion = " D3BJ"

    ghost_map = qm_options.get("ghost_map", {})
    def _finalize(content: str, fname: str) -> Tuple[str, str, List[str]]:
        # Loophole 6: Restore ghost/dummy atoms
        if ghost_map:
            content = GhostAtomNormalizer.restore_ghost_atoms_in_output(content, ghost_map, target_format)
        # Loophole 7: Universal Unix EOL normalization
        content = UnixEolNormalizer.normalize(content)
        return content, fname, warnings

    # 1. ORCA INP (Issues 8, 10, 11)
    if target_format in ["orca", "inp"]:
        lines = [
            f"# Generated by iLoveChem Studio",
            f"# Electron count: {total_z - charge} | Parity verified",
            f"! {method} {basis} {calc_type}{orca_dispersion}",
            f"%maxcore {int((mem_gb / nprocs) * 1024)}",
            f"%pal nprocs {nprocs} end",
            f"",
            f"* xyz {charge} {mult}",
        ]
        # Strict Cartesian coordinate formatting (Issue 11)
        for atom in mol.atoms:
            sym = ob.GetSymbol(atom.OBAtom.GetAtomicNum())
            x, y, z = atom.coords
            lines.append(f"  {sym:<3}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
        lines.append("*")
        return _finalize("\n".join(lines), "calculation.inp")

    # 2. Gaussian GJF / COM (Issues 8, 9, 11)
    if target_format in ["gjf", "com", "gaussian"]:
        lines = [
            f"%nprocshared={nprocs}",
            f"%mem={mem_gb}GB",
            f"%chk=calculation.chk",
            f"#p {method}/{basis} {calc_type}",
            f"",
            f"Generated by iLoveChem (Total electrons: {total_z - charge})",
            f"",
            f"{charge} {mult}",
        ]
        # Strict Cartesian coordinate formatting (Issue 11)
        for atom in mol.atoms:
            sym = ob.GetSymbol(atom.OBAtom.GetAtomicNum())
            x, y, z = atom.coords
            lines.append(f"{sym:<3}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
        # Guaranteed trailing 2 blank lines to prevent Gaussian QPErr / EOF crash (Issue 9)
        lines.append("")
        lines.append("")
        out_fname = "molecule.gjf" if target_format == "gaussian" else f"molecule.{target_format}"
        return _finalize("\n".join(lines), out_fname)

    # 3. Standard 3D XYZ
    if target_format == "xyz":
        lines = [
            str(len(mol.atoms)),
            f"Generated by iLoveChem | Formula: {mol.formula} | Charge: {charge}"
        ]
        for atom in mol.atoms:
            sym = ob.GetSymbol(atom.OBAtom.GetAtomicNum())
            x, y, z = atom.coords
            lines.append(f"{sym:<3}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
        return _finalize("\n".join(lines), "structure.xyz")

    # 4. Q-Chem
    if target_format in ["qchem", "in"]:
        lines = [
            "$molecule",
            f"{charge} {mult}"
        ]
        for atom in mol.atoms:
            sym = ob.GetSymbol(atom.OBAtom.GetAtomicNum())
            x, y, z = atom.coords
            lines.append(f"{sym:<3}  {x:12.6f}  {y:12.6f}  {z:12.6f}")
        lines.append("$end\n")
        lines.append("$rem")
        # Honor the requested calculation type (Q-Chem runs one job type per section;
        # a combined 'opt freq' starts with the optimization).
        calc_lc = calc_type.lower()
        if "opt" in calc_lc:
            qchem_job = "OPT"
        elif "freq" in calc_lc:
            qchem_job = "FREQ"
        else:
            qchem_job = "SP"
        lines.append(f"  JOBTYPE       {qchem_job}\n  METHOD        {method}\n  BASIS         {basis}\n$end\n")
        return _finalize("\n".join(lines), "qchem.in")

    # 5. Solid-State: VASP POSCAR with 15 Å Vacuum Box (Issue 15)
    if target_format == "poscar":
        lattice_mat, centered_coords = CrystalSolidStateEngine.create_vacuum_box(mol, vacuum_padding=15.0)
        
        # Group atoms by element for VASP 5.x convention
        elem_groups = {}
        for idx, atom in enumerate(mol.atoms):
            sym = ob.GetSymbol(atom.OBAtom.GetAtomicNum())
            if sym not in elem_groups:
                elem_groups[sym] = []
            elem_groups[sym].append(centered_coords[idx])

        poscar_lines = [
            f"iLoveChem Box (15A vacuum padding) - {mol.formula}",
            "1.0",
            f"  {lattice_mat[0,0]:14.8f}  {lattice_mat[0,1]:14.8f}  {lattice_mat[0,2]:14.8f}",
            f"  {lattice_mat[1,0]:14.8f}  {lattice_mat[1,1]:14.8f}  {lattice_mat[1,2]:14.8f}",
            f"  {lattice_mat[2,0]:14.8f}  {lattice_mat[2,1]:14.8f}  {lattice_mat[2,2]:14.8f}",
            "  " + "  ".join(elem_groups.keys()),
            "  " + "  ".join(str(len(v)) for v in elem_groups.values()),
            "Cartesian"
        ]
        for sym, c_list in elem_groups.items():
            for c in c_list:
                poscar_lines.append(f"  {c[0]:14.8f}  {c[1]:14.8f}  {c[2]:14.8f}")

        return _finalize("\n".join(poscar_lines), "POSCAR")

    # 6. Solid-State: CIF / PWscf via ASE with Vacuum Box
    if target_format in ["cif", "pwi"]:
        lattice_mat, centered_coords = CrystalSolidStateEngine.create_vacuum_box(mol, vacuum_padding=15.0)
        symbols = [ob.GetSymbol(a.OBAtom.GetAtomicNum()) for a in mol.atoms]
        atoms = Atoms(symbols=symbols, positions=centered_coords, cell=lattice_mat, pbc=True)
        buf = io.StringIO()
        
        if target_format == "cif":
            byte_buf = io.BytesIO()
            ase_io.write(byte_buf, atoms, format="cif")
            cif_content = byte_buf.getvalue().decode('utf-8')
            return _finalize(cif_content, "crystal.cif")
        elif target_format == "pwi":
            pseudos = {sym: f"{sym}.pbe-rrkjus.UPF" for sym in set(symbols)}
            ase_io.write(buf, atoms, format="espresso-in", pseudopotentials=pseudos)
            return _finalize(buf.getvalue(), "pwscf.pwi")

    # 7. AutoDock PDBQT (Issue 18)
    if target_format in ["pdbqt", "vina"]:
        try:
            pdbqt_str = mol.write("pdbqt")
            return _finalize(pdbqt_str, "docking_ligand.pdbqt")
        except Exception:
            pass

    # 8. GROMACS GRO with 1.0 nm clearance (Issue 19)
    if target_format in ["gro", "gromacs"]:
        gro_str, gro_fname, gro_warns = GromacsGroEngine.export_gro(mol)
        warnings.extend(gro_warns)
        return _finalize(gro_str, gro_fname)

    # 9. General OpenBabel Formats
    try:
        out_str = mol.write(target_format)
        ext_map = {
            "mop": ".mop", "nw": ".nw", "gam": ".inp", "fhiaims": ".in",
            "lmpdat": ".data", "xtl": ".xtl", "pdbqt": ".pdbqt", "pdb": ".pdb",
            "mol2": ".mol2", "gro": ".gro", "pqr": ".pqr", "cdxml": ".cdxml",
            "mol": ".mol", "sdf": ".sdf", "smi": ".smi", "inchi": ".txt",
            "inchikey": ".txt", "cml": ".cml", "txyz": ".xyz", "svg": ".svg"
        }
        ext = ext_map.get(target_format, f".{target_format}")
        return _finalize(out_str, f"converted{ext}")
    except Exception:
        # Never fail silently: tell the user the requested format could not be written.
        warnings.append(
            f"Could not write the requested format '{target_format}'. "
            f"The structure was exported as standard XYZ instead."
        )
        out_str = mol.write("xyz")
        return _finalize(out_str, "converted.xyz")
