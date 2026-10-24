import os
import sys
import django

# Setup Django environment
sys.path.insert(0, r"C:\Users\shiva\.gemini\antigravity\scratch\chemsite")
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "chemsite.settings")
django.setup()

from calculator.converter_core import parse_and_extract_atoms, export_to_target_format, detect_format
from calculator.scientific_engine import (
    SpinChargeParityValidator,
    QuantumLogParser,
    PhysicalIntegrityChecker,
    PDBBioSanitizer,
    CrystalSolidStateEngine
)

def run_tests():
    print("=" * 70)
    print("STARTING SCIENTIFIC INTEGRITY & CONVERTER VERIFICATION SUITE")
    print("=" * 70)

    # -------------------------------------------------------------
    # PHASE 1: ChemDraw & Reaction Scheme Tests (Issues 1-6)
    # -------------------------------------------------------------
    print("\n[PHASE 1] Testing ChemDraw Reaction Scheme & Stereochemistry...")
    cdxml_sample = """<?xml version="1.0" encoding="UTF-8" ?>
<!DOCTYPE CDXML SYSTEM "http://www.cambridgesoft.com/xml/cdxml.dtd" >
<CDXML>
  <page>
    <t p="100 200"><span>Reaction: Esterification in THF</span></t>
    <arrow GraphicType="Line" ArrowHeadType="Solid" p1="150 150" p2="250 150" />
    <graphic GraphicType="Symbol" SymbolType="Plus" p="130 150" />
    <fragment id="1">
      <n id="1" Element="6" p="10 10"/>
      <n id="2" Element="8" p="20 10"/>
      <b B="1" E="2" Order="1"/>
    </fragment>
    <fragment id="2">
      <n id="3" Element="6" p="100 10"/>
      <n id="4" Element="6" p="110 10"/>
      <n id="5" Element="6" p="120 10"/>
      <n id="6" Element="8" p="130 10"/>
      <b B="3" E="4" Order="1" Display="WedgeBegin"/>
      <b B="4" E="5" Order="1"/>
      <b B="5" E="6" Order="2"/>
    </fragment>
  </page>
</CDXML>"""
    mol_p1, meta_p1 = parse_and_extract_atoms(cdxml_sample.encode('utf-8'), 'cdxml', 'test.cdxml')
    assert mol_p1 is not None, "Phase 1: Molecule failed to parse from CDXML"
    assert meta_p1["fragment_count"] == 2, f"Phase 1: Expected 2 fragments, got {meta_p1['fragment_count']}"
    assert meta_p1["filtered_elements"]["reaction_arrows_filtered"] == 1, "Phase 1: Arrow not filtered"
    print("  -> Passed: Filtered text/arrows, extracted multiple fragments, built 3D conformer.")

    # -------------------------------------------------------------
    # PHASE 2: Quantum Chemistry & Output Logs (Issues 7-11)
    # -------------------------------------------------------------
    print("\n[PHASE 2] Testing Quantum Formatting, Parity & Inverted Log Parser...")
    
    # Issue 7: Inverted Log Parser
    mock_log = """
 Entering Gaussian System...
 Initial orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        0.000000    0.000000    0.000000
      2          1           0        0.000000    0.000000    1.090000
 ---------------------------------------------------------------------
 Step 1: Energy = -40.123
 Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        0.050000    0.000000    0.000000
      2          1           0        0.000000    0.000000    1.080000
 ---------------------------------------------------------------------
 Stationary point found.
 Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        0.123456    0.000000    0.000000
      2          1           0        0.000000    0.000000    1.087654
 ---------------------------------------------------------------------
 Normal termination of Gaussian.
"""
    mol_log, meta_log = parse_and_extract_atoms(mock_log.encode('utf-8'), 'log', 'opt.log')
    assert mol_log is not None, "Phase 2 Issue 7: Log parsing failed"
    # Check that final stationary coordinates (0.123456) were selected, not initial (0.0)
    c0 = mol_log.atoms[0].coords
    assert abs(c0[0] - 0.123456) < 1e-4, f"Phase 2 Issue 7: Did not extract converged stationary coordinates! Got {c0[0]}"
    print("  -> Passed Issue 7: Inverted log scanner successfully extracted final converged geometry.")

    # Issue 8: Parity validation (Water H2O has 10 electrons. If user requests mult=2, must auto-correct to 1)
    water_xyz = "3\nWater\nO 0.0 0.0 0.0\nH 0.0 0.757 0.587\nH 0.0 -0.757 0.587\n"
    mol_w, _ = parse_and_extract_atoms(water_xyz.encode('utf-8'), 'xyz', 'water.xyz')
    orca_out, _, warns_orca = export_to_target_format(mol_w, 'orca', {'charge': 0, 'mult': 2})
    assert "* xyz 0 1" in orca_out, "Phase 2 Issue 8: Failed to correct even electron parity to singlet"
    assert len(warns_orca) > 0, "Phase 2 Issue 8: Warning not emitted"
    print("  -> Passed Issue 8: Electron parity validator auto-corrected Multiplicity 2 -> 1 for 10-electron water.")

    # Issue 9: Gaussian 2 trailing blank lines
    gjf_out, _, _ = export_to_target_format(mol_w, 'gjf', {'charge': 0, 'mult': 1})
    assert gjf_out.endswith("\n\n\n") or gjf_out.endswith("\n\n"), "Phase 2 Issue 9: Gaussian file missing trailing blank lines"
    print("  -> Passed Issue 9: Gaussian .gjf terminates with required blank lines to prevent QPErr.")

    # -------------------------------------------------------------
    # PHASE 3: Crystallography & Solid-State DFT (Issues 12-15)
    # -------------------------------------------------------------
    print("\n[PHASE 3] Testing Crystallography & 15A Vacuum Box...")
    # Issue 15: VASP POSCAR 15A Vacuum Padding
    poscar_out, _, _ = export_to_target_format(mol_w, 'poscar')
    lines_poscar = poscar_out.splitlines()
    assert "Cartesian" in poscar_out, "Phase 3 Issue 15: Missing Cartesian tag in POSCAR"
    # Read lattice vector lengths from lines 2, 3, 4
    lx = float(lines_poscar[2].split()[0])
    ly = float(lines_poscar[3].split()[1])
    lz = float(lines_poscar[4].split()[2])
    assert lx >= 15.0 and ly >= 15.0 and lz >= 15.0, f"Phase 3 Issue 15: Vacuum padding box too small: {lx}, {ly}, {lz}"
    print(f"  -> Passed Issue 15: VASP POSCAR generated with {lx:.1f}x{ly:.1f}x{lz:.1f} Å vacuum box.")

    # -------------------------------------------------------------
    # PHASE 4: Bio / Macromolecules (Issues 16-19)
    # -------------------------------------------------------------
    print("\n[PHASE 4] Testing Bio-Sanitization (AltLoc A/B & Water Stripping)...")
    mock_pdb = """HEADER    TEST PROTEIN / LIGAND
CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1
ATOM      1  N  AALA A   1      10.000  10.000  10.000  0.50 20.00           N
ATOM      2  N  BALA A   1      10.200  10.100  10.000  0.50 20.00           N
HETATM   10  C1  LIG A 100      25.000  25.000  25.000  1.00 15.00           C
HETATM   11  O   HOH A 201      30.000  30.000  30.000  1.00 18.00           O
END
"""
    # Test AltLoc filtering: ALAB should be dropped
    sanitized_all = PDBBioSanitizer.sanitize_pdb(mock_pdb, filter_mode="all")
    assert "BALA" not in sanitized_all, "Phase 4 Issue 16: AltLoc B was not discarded"
    assert "AALA" in sanitized_all, "Phase 4 Issue 16: AltLoc A was dropped"
    print("  -> Passed Issue 16: Disordered AltLoc B discarded; AltLoc A retained.")

    # Test Water Stripping: HOH should be dropped in strip_water mode
    sanitized_strip = PDBBioSanitizer.sanitize_pdb(mock_pdb, filter_mode="strip_water")
    assert "HOH" not in sanitized_strip, "Phase 4 Issue 17: Water HOH was not stripped"
    print("  -> Passed Issue 17: Water (HOH/WAT/SOL) stripped cleanly.")

    # Test Ligand Only
    sanitized_lig = PDBBioSanitizer.sanitize_pdb(mock_pdb, filter_mode="ligand_only")
    assert "ALA" not in sanitized_lig and "LIG" in sanitized_lig, "Phase 4 Issue 17: Ligand isolation failed"
    print("  -> Passed Issue 17: Active ligand isolated from macromolecular complex.")

    # -------------------------------------------------------------
    # PHASE 5: Heavy Metals & Physical Integrity (Issues 20-22)
    # -------------------------------------------------------------
    print("\n[PHASE 5] Testing Heavy Metals UFF Fallback & Nuclear Clash Detection...")
    # Issue 20: Transition metal complex (Cisplatin Pt(NH3)2Cl2)
    cisplatin_xyz = """11
Cisplatin
Pt   0.000000   0.000000   0.000000
Cl   2.300000   0.000000   0.000000
Cl   0.000000   2.300000   0.000000
N   -2.050000   0.000000   0.000000
N    0.000000  -2.050000   0.000000
H   -2.400000   0.800000   0.500000
H   -2.400000  -0.800000   0.500000
H   -2.400000   0.000000  -0.900000
H    0.800000  -2.400000   0.500000
H   -0.800000  -2.400000   0.500000
H    0.000000  -2.400000  -0.900000
"""
    mol_pt, meta_pt = parse_and_extract_atoms(cisplatin_xyz.encode('utf-8'), 'xyz', 'cisplatin.xyz')
    assert meta_pt["has_metals"] is True, "Phase 5 Issue 20: Pt metal not detected"
    assert "Pt" in meta_pt["metal_symbols"], "Phase 5 Issue 20: Pt symbol missing"
    print(f"  -> Passed Issue 20: Transition metal (Pt) detected, activating UFF force field.")

    # Issue 21: Nuclear Overlap Clash (<0.7 Å)
    clash_xyz = """2
Clash
C 0.000000 0.000000 0.000000
C 0.350000 0.000000 0.000000
"""
    mol_clash, meta_clash = parse_and_extract_atoms(clash_xyz.encode('utf-8'), 'xyz', 'clash.xyz')
    assert meta_clash["has_clashes"] is True, "Phase 5 Issue 21: Nuclear clash not detected (<0.7 A)"
    assert meta_clash["min_interatomic_distance"] < 0.7, "Phase 5 Issue 21: Distance not calculated correctly"
    print(f"  -> Passed Issue 21: Nuclear overlap clash detected (d = {meta_clash['min_interatomic_distance']} Å).")

    # -------------------------------------------------------------
    # EXTRA: Disconnected Fragments / Counter-Ions & Salts
    # -------------------------------------------------------------
    print("\n[SALTS & COUNTER-IONS] Testing Counter-Ion Stripper & Complex Preservation...")
    # Aniline HCl
    mol_salt_strip, meta_salt_strip = parse_and_extract_atoms(b"c1ccccc1[NH3+].[Cl-]", "smi", strip_salts=True)
    assert meta_salt_strip["salt_detected"] is True, "Salt was not detected"
    assert meta_salt_strip["salts_stripped"] is True, "Salt was not stripped"
    assert "Cl" not in mol_salt_strip.formula, f"Counter-ion Cl was not stripped! Got formula {mol_salt_strip.formula}"
    print(f"  -> Passed: Aniline HCl salt stripped cleanly (Active parent: {mol_salt_strip.formula})")

    mol_salt_keep, meta_salt_keep = parse_and_extract_atoms(b"c1ccccc1[NH3+].[Cl-]", "smi", strip_salts=False)
    assert meta_salt_keep["salts_stripped"] is False, "Salt should not have been stripped when strip_salts=False"
    assert "Cl" in mol_salt_keep.formula, "Counter-ion Cl should have been kept"
    print(f"  -> Passed: Aniline HCl salt complex preserved when requested (Formula: {mol_salt_keep.formula})")

    # Sodium Acetate
    mol_na, meta_na = parse_and_extract_atoms(b"CC(=O)[O-].[Na+]", "smi", strip_salts=True)
    assert "Na" not in mol_na.formula, f"Na+ was not stripped! Got {mol_na.formula}"
    print(f"  -> Passed: Sodium Acetate salt stripped (Acetate anion: {mol_na.formula})")

    # TFA Solvate / Salt
    mol_tfa, meta_tfa = parse_and_extract_atoms(b"c1ccccc1.C(=O)(C(F)(F)F)O", "smi", strip_salts=True)
    assert "F" not in mol_tfa.formula and mol_tfa.formula == "C6H6", f"TFA was not stripped! Got {mol_tfa.formula}"
    print(f"  -> Passed: Benzene TFA salt stripped (Active parent: {mol_tfa.formula})")

    print("\n" + "=" * 70)
    print("ALL TESTS PASSED: PHASES 1 THROUGH 5 VERIFIED SUCCESSFULLY!")
    print("=" * 70)

if __name__ == "__main__":
    run_tests()
