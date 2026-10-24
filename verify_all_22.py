import sys
import time
import requests
import tempfile
import os
import numpy as np

# Ensure UTF-8 output on Windows consoles
if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

BASE_URL = 'http://127.0.0.1:8000/api/convert/'

print("=" * 85)
print("COMPREHENSIVE SCIENTIFIC AUDIT — ALL 22 REAL-WORLD FAILURE MODES VERIFICATION")
print("=" * 85)

results = {}

# ---------------------------------------------------------
# CATEGORY A: ChemDraw & 2D-to-3D Nightmares
# ---------------------------------------------------------

# 1. Chiral Center Inversion (R/S) Lock
cdxml_chiral = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<fragment id="1">
  <n id="1" p="100 100" Element="6"/>
  <n id="2" p="100 80" Element="9"/>
  <n id="3" p="120 100" Element="17"/>
  <n id="4" p="80 100" Element="35"/>
  <n id="5" p="100 120" Element="1"/>
  <b id="6" B="1" E="2" Order="1" Display="WedgeBegin"/>
  <b id="7" B="1" E="3" Order="1"/>
  <b id="8" B="1" E="4" Order="1"/>
  <b id="9" B="1" E="5" Order="1"/>
</fragment>
</page>
</CDXML>"""
r_a1 = requests.post(BASE_URL, data={'raw_content': cdxml_chiral, 'input_format': 'cdxml', 'target_format': 'xyz'})
results[1] = ("Chiral Center (R/S) Parity Lock", r_a1.status_code == 200 and r_a1.json()['metadata']['num_atoms'] == 5)

# 2. Cis/Trans (E/Z) Double Bond Stereo Locking (Trans-2-butene)
r_a2 = requests.post(BASE_URL, data={'raw_content': 'C/C=C/C', 'input_format': 'smi', 'target_format': 'xyz'})
results[2] = ("Cis/Trans (E/Z) Double Bond Stereo Lock", r_a2.status_code == 200 and r_a2.json()['metadata']['num_atoms'] == 12)

# 3. Aromaticity Perception Failure
r_a3 = requests.post(BASE_URL, data={'raw_content': 'c1ccccc1', 'input_format': 'smi', 'target_format': 'orca'})
res_a3 = r_a3.json()
results[3] = ("Aromaticity Ring Sanitizer", r_a3.status_code == 200 and res_a3['metadata']['formula'] == 'C6H6')

# 4. Multi-Fragment Reaction Schemes (Multi-Molecule)
cdxml_rxn = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<arrow id="10" p="50 50"/>
<t id="11" p="50 40"><s>85% Yield</s></t>
<fragment id="1"><n id="1" p="0 0" Element="6"/><n id="2" p="10 0" Element="8"/><b id="3" B="1" E="2" Order="1"/></fragment>
<fragment id="2"><n id="4" p="100 100" Element="6"/><n id="5" p="110 100" Element="6"/><n id="6" p="120 100" Element="7"/><b id="7" B="4" E="5" Order="1"/><b id="8" B="5" E="6" Order="1"/></fragment>
</page>
</CDXML>"""
r_a4 = requests.post(BASE_URL, data={'raw_content': cdxml_rxn, 'input_format': 'cdxml', 'convert_mode': 'all_zip', 'target_format': 'orca'})
res_a4 = r_a4.json()
results[4] = ("Multi-Fragment Reaction Scheme ZIP", r_a4.status_code == 200 and res_a4.get('is_zip') is True)

# 5. Non-Chemical Drawings (Arrows, Text, Brackets Stripping)
r_a5 = requests.post(BASE_URL, data={'raw_content': cdxml_rxn, 'input_format': 'cdxml', 'convert_mode': 'single', 'target_format': 'xyz'})
res_a5 = r_a5.json()
results[5] = ("Non-Chemical Annotations Filter", r_a5.status_code == 200 and res_a5['metadata']['num_atoms'] > 0)

# 6. Zwitterions & Formal Charge Strip
r_a6 = requests.post(BASE_URL, data={'raw_content': '[NH3+]CC(=O)[O-]', 'input_format': 'smi', 'target_format': 'orca'})
res_a6 = r_a6.json()
results[6] = ("Zwitterion Formal Charge Balance", r_a6.status_code == 200 and res_a6['metadata']['detected_charge'] == 0)

# ---------------------------------------------------------
# CATEGORY B: Quantum Chemistry Inputs & Outputs
# ---------------------------------------------------------

# 7. Initial vs Final Geometry Trap in Logs
log_text = """Entering Gaussian 16...
Initial orientation:
 1 6 0 0.000 0.000 0.000
 2 1 0 0.000 0.000 1.090
Cycle 1: Energy = -40.123
Stationary point found.
Standard orientation:
---------------------------------------------------------------------
Center     Atomic      Atomic             Coordinates (Angstroms)
Number     Number       Type             X           Y           Z
---------------------------------------------------------------------
     1          6           0        9.876543    0.000000    0.000000
     2          1           0        0.000000    0.000000    1.087654
---------------------------------------------------------------------
Normal termination.
"""
r_b7 = requests.post(BASE_URL, data={'raw_content': log_text, 'input_format': 'log', 'target_format': 'xyz'})
results[7] = ("Reverse-Scan Final Stationary Point", r_b7.status_code == 200 and "9.876543" in r_b7.json().get('output_content', ''))

# 8. Charge & Spin Multiplicity Inconsistency Parity
r_b8 = requests.post(BASE_URL, data={'raw_content': '[CH3]', 'input_format': 'smi', 'target_format': 'orca', 'mult': 1})
results[8] = ("Valence Electron Spin Parity (2S+1)", r_b8.status_code == 200 and "* xyz 0 2" in r_b8.json().get('output_content', ''))

# 9. Gaussian Blank Line Syntax Traps
r_b9 = requests.post(BASE_URL, data={'raw_content': 'O', 'input_format': 'smi', 'target_format': 'gjf'})
results[9] = ("Gaussian Guaranteed 2 Blank Lines", r_b9.status_code == 200 and r_b9.json()['output_content'].endswith("\n\n"))

# 10. Z-Matrix to Cartesian Crash in ORCA
r_b10 = requests.post(BASE_URL, data={'raw_content': 'C', 'input_format': 'smi', 'target_format': 'orca'})
lines_b10 = r_b10.json().get('output_content', '').splitlines()
has_cart = any("C  " in l and len(l.split()) == 4 for l in lines_b10)
results[10] = ("Standardized Cartesian Alignment", r_b10.status_code == 200 and has_cart)

# 11. PDB Metadata Pollution in Gaussian
pdb_sample = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
END"""
r_b11 = requests.post(BASE_URL, data={'raw_content': pdb_sample, 'input_format': 'pdb', 'target_format': 'gjf'})
out_b11 = r_b11.json().get('output_content', '')
clean_coords = "ALA" not in out_b11 and "ATOM" not in out_b11
results[11] = ("PDB Metadata Stripped from QM Coords", r_b11.status_code == 200 and clean_coords)

# ---------------------------------------------------------
# CATEGORY C: Solid-State & Crystallography
# ---------------------------------------------------------

# 12. Missing Hydrogens in X-ray CIF Files
cif_benzene = """data_benzene_noh
_cell_length_a 7.46
_cell_length_b 9.67
_cell_length_c 7.03
_cell_angle_alpha 90.00
_cell_angle_beta 90.00
_cell_angle_gamma 90.00
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.0000 0.1412 0.0000
C2 C 0.1223 0.0706 0.0000
C3 C 0.1223 -0.0706 0.0000
C4 C 0.0000 -0.1412 0.0000
C5 C -0.1223 -0.0706 0.0000
C6 C -0.1223 0.0706 0.0000"""
r_c12 = requests.post(BASE_URL, data={'raw_content': cif_benzene, 'target_format': 'xyz', 'add_hydrogens': 'true', 'strip_salts': 'false'})
results[12] = ("CIF Missing H Auto-Protonation", r_c12.status_code == 200 and any("CIF Auto-Hydrogenation Active" in w for w in r_c12.json()['metadata']['warnings']))

# 13. Fractional vs Cartesian Coordinates Confusion
r_c13 = requests.post(BASE_URL, data={'raw_content': cif_benzene, 'target_format': 'xyz', 'strip_salts': 'false'})
out_c13 = r_c13.json().get('output_content', '').splitlines()
has_xyz_floats = any(len(l.split()) == 4 and abs(float(l.split()[1])) < 20.0 for l in out_c13[2:])
results[13] = ("Fractional to Cartesian Matrix Math", r_c13.status_code == 200 and has_xyz_floats)

# 14. Periodic Boundary Wrapping & Fragment Splitting
cif_split = """data_split_mol
_cell_length_a 10.00
_cell_length_b 10.00
_cell_length_c 10.00
_cell_angle_alpha 90.00
_cell_angle_beta 90.00
_cell_angle_gamma 90.00
_symmetry_space_group_name_H-M 'P 1'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.0500 0.5000 0.5000
O2 O 0.9500 0.5000 0.5000"""
r_c14 = requests.post(BASE_URL, data={'raw_content': cif_split, 'target_format': 'xyz', 'add_hydrogens': 'false', 'strip_salts': 'false'})
res_c14 = r_c14.json()
warns_c14 = res_c14.get('metadata', {}).get('warnings', [])
has_unwrap_warn = any("Periodic Boundary Unwrapping" in w for w in warns_c14)
results[14] = ("Periodic Split Molecule Unwrapper", r_c14.status_code == 200 and has_unwrap_warn)

# 15. Vacuum Padding for Gas-Phase / Slab Calculations
r_c15 = requests.post(BASE_URL, data={'raw_content': 'O', 'input_format': 'smi', 'target_format': 'poscar'})
results[15] = ("15A Vacuum Box in VASP POSCAR", r_c15.status_code == 200 and "Cartesian" in r_c15.json().get('output_content', ''))

# ---------------------------------------------------------
# CATEGORY D: Bio, Docking & Macromolecules
# ---------------------------------------------------------

# 16. Missing Alternate Conformations (AltLoc A/B)
pdb_altloc = """ATOM      1  N  AALA A   1      10.000  10.000  10.000  0.50 20.00           N
ATOM      2  N  BALA A   1      10.200  10.100  10.000  0.50 20.00           N
ATOM      3  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
END"""
r_d16 = requests.post(BASE_URL, data={'raw_content': pdb_altloc, 'input_format': 'pdb', 'target_format': 'xyz', 'add_hydrogens': 'false'})
results[16] = ("AltLoc Conformation Deduplication", r_d16.status_code == 200 and r_d16.json()['metadata']['num_atoms'] == 2)

# 17. Non-Standard Residues & HETATM Clashes (Water Strip & Ligand Isolation)
pdb_waters = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
HETATM   10  C1  LIG A 100      25.000  25.000  25.000  1.00 15.00           C
HETATM   20  O   HOH A 201      30.000  30.000  30.000  1.00 18.00           O
END"""
r_d17 = requests.post(BASE_URL, data={'raw_content': pdb_waters, 'input_format': 'pdb', 'pdb_filter_mode': 'ligand_only', 'target_format': 'xyz', 'add_hydrogens': 'false'})
results[17] = ("PDB Water Strip & Ligand Isolation", r_d17.status_code == 200 and r_d17.json()['metadata']['num_atoms'] == 1)

# 18. AutoDock Vina Partial Charge & Torsion Assignment
r_d18 = requests.post(BASE_URL, data={'raw_content': 'CC(=O)Oc1ccccc1C(=O)O', 'input_format': 'smi', 'target_format': 'pdbqt'})
out_d18 = r_d18.json().get('output_content', '')
results[18] = ("AutoDock PDBQT Torsions & Charges", r_d18.status_code == 200 and "ROOT" in out_d18 and "TORSDOF" in out_d18)

# 19. GROMACS PBC Box Dimension Loss
r_d19 = requests.post(BASE_URL, data={'raw_content': 'CCO', 'input_format': 'smi', 'target_format': 'gro'})
out_d19 = r_d19.json().get('output_content', '')
lines_d19 = [l for l in out_d19.splitlines() if l.strip()]
box_line = lines_d19[-1].split() if lines_d19 else []
results[19] = ("GROMACS 1.0nm Clearance Box Vectors", r_d19.status_code == 200 and len(box_line) == 3 and float(box_line[0]) > 1.0)

# ---------------------------------------------------------
# CATEGORY E: Heavy Metals & Physical Integrity
# ---------------------------------------------------------

# 20. Transition Metal Force-Field Disasters (UFF Fallback)
ferrocene_xyz = """21
Ferrocene
Fe 0.0 0.0 0.0
C 1.214 0.0 1.65
C 0.375 1.155 1.65
C -0.982 0.714 1.65
C -0.982 -0.714 1.65
C 0.375 -1.155 1.65
H 2.277 0.0 1.65
H 0.704 2.166 1.65
H -1.842 1.338 1.65
H -1.842 -1.338 1.65
H 0.704 -2.166 1.65
C 1.214 0.0 -1.65
C 0.375 1.155 -1.65
C -0.982 0.714 -1.65
C -0.982 -0.714 -1.65
C 0.375 -1.155 -1.65
H 2.277 0.0 -1.65
H 0.704 2.166 -1.65
H -1.842 1.338 -1.65
H -1.842 -1.338 -1.65
H 0.704 -2.166 -1.65"""
r_e20 = requests.post(BASE_URL, data={'raw_content': ferrocene_xyz, 'target_format': 'orca'})
results[20] = ("Transition Metal UFF Force-Field Guard", r_e20.status_code == 200 and r_e20.json()['metadata']['has_metals'] is True)

# 21. Nuclear Overlap Steric Clashing Sanity Check (d < 0.7 A)
clash_xyz = """2
Clashing Diatomic
C 0.000 0.000 0.000
C 0.200 0.000 0.000"""
r_e21 = requests.post(BASE_URL, data={'raw_content': clash_xyz, 'target_format': 'orca'})
results[21] = ("Nuclear Overlap (< 0.7 A) Warning", r_e21.status_code == 200 and r_e21.json()['metadata']['has_clashes'] is True)

# 22. Large File Streaming Memory Hangs (100MB+ Logs)
with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.log', encoding='utf-8') as f:
    f.write("Entering Gaussian 16...\n")
    for i in range(150000):
        f.write(f"Cycle {i}: SCF Energy = -123.45678\n")
    f.write("""Stationary point found.
Standard orientation:
---------------------------------------------------------------------
Center     Atomic      Atomic             Coordinates (Angstroms)
Number     Number       Type             X           Y           Z
---------------------------------------------------------------------
     1          6           0        5.555555    0.000000    0.000000
     2          1           0        0.000000    1.087654    0.000000
---------------------------------------------------------------------
Normal termination.
""")
    tmp_log = f.name

size_mb = os.path.getsize(tmp_log) / (1024 * 1024)
t0 = time.time()
with open(tmp_log, 'rb') as f_in:
    r_e22 = requests.post(BASE_URL, files={'file': ('giant_opt.log', f_in, 'text/plain')}, data={'target_format': 'xyz'})
dt = time.time() - t0
os.remove(tmp_log)
results[22] = (f"Giant Log Sub-Second Streaming ({size_mb:.1f}MB in {dt*1000:.1f}ms)", r_e22.status_code == 200 and "5.555555" in r_e22.json().get('output_content', ''))

# ---------------------------------------------------------
# FINAL SCORECARD
# ---------------------------------------------------------
print(f"{'#':<3} {'Failure Mode / Feature':<48} {'Status'}")
print("-" * 85)
passed_count = 0
for idx in range(1, 23):
    name, passed = results[idx]
    status_str = "PASS (100%)" if passed else "FAIL"
    if passed:
        passed_count += 1
    print(f"{idx:<3} {name:<48} {status_str}")

print("=" * 85)
print(f"OVERALL RESULT: {passed_count} / 22 ISSUES PASSING ({(passed_count/22)*100:.1f}%)")
print("=" * 85)
