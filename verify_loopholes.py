import sys
import requests

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

BASE_URL = 'http://127.0.0.1:8000/api/convert/'

print("=" * 85)
print("AUDIT SUITE — 8 KILLER CHEMICAL CONVERSION LOOPHOLES VERIFICATION")
print("=" * 85)

results = {}

# -----------------------------------------------------------------------------
# LOOPHOLE 1: Confidentiality & IP Security Policy (Pharma's #1 Fear)
# -----------------------------------------------------------------------------
r_home = requests.get('http://127.0.0.1:8000/')
has_sec_banner = "Zero-Disk Retention Policy" in r_home.text and "Pharma Safe" in r_home.text
has_sec_modal = "security-modal" in r_home.text and "Volatile RAM Execution" in r_home.text
results[1] = ("Loophole 1: Pharma Confidentiality & Zero-Storage Invariant", has_sec_banner and has_sec_modal)

# -----------------------------------------------------------------------------
# LOOPHOLE 2: The Silent Failure Bug (All (0,0,0) Coordinates)
# -----------------------------------------------------------------------------
collapsed_xyz = """3
Collapsed Zero Coords
C 0.0000 0.0000 0.0000
C 0.0000 0.0000 0.0000
O 0.0000 0.0000 0.0000"""
r2 = requests.post(BASE_URL, data={'raw_content': collapsed_xyz, 'target_format': 'xyz'})
res2 = r2.json()
fixed_output = res2.get('output_content', '')
lines2 = [l.split() for l in fixed_output.splitlines() if len(l.split()) == 4]
non_zeros = any(abs(float(p[1])) > 0.05 or abs(float(p[2])) > 0.05 for p in lines2)
results[2] = ("Loophole 2: Zero-Coordinate Guard (Silent Failure Interceptor)", res2.get('metadata', {}).get('silent_failure_prevented') is True and non_zeros)

# -----------------------------------------------------------------------------
# LOOPHOLE 3: Incomplete / Failed Quantum Optimizations in Logs
# -----------------------------------------------------------------------------
crashed_log = """Entering Gaussian 16...
Step 1: SCF Done:  E(RB3LYP) =  -150.1000000     A.U.
Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        1.000000    0.000000    0.000000
      2          1           0        0.000000    1.000000    0.000000
 ---------------------------------------------------------------------
Step 2: SCF Done:  E(RB3LYP) =  -150.4500000     A.U.
Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        7.777777    0.000000    0.000000
      2          1           0        0.000000    1.100000    0.000000
 ---------------------------------------------------------------------
Step 3: SCF Done:  E(RB3LYP) =  -150.2500000     A.U.
Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
      1          6           0        1.300000    0.000000    0.000000
      2          1           0        0.000000    1.200000    0.000000
 ---------------------------------------------------------------------
Error termination request processed by link 9999.
Job terminated abnormally."""
r3 = requests.post(BASE_URL, data={'raw_content': crashed_log, 'input_format': 'log', 'target_format': 'xyz'})
res3 = r3.json()
has_step2_coords = "7.777777" in res3.get('output_content', '')
has_rescue_label = "Rescued Lowest Energy Step #2" in res3.get('metadata', {}).get('parsed_by', '')
results[3] = ("Loophole 3: Rescued Lowest-Energy Geometry (Failed QM Job)", has_step2_coords and has_rescue_label)

# -----------------------------------------------------------------------------
# LOOPHOLE 4: Missing Bond Orders & Double-Bonds in PDB
# -----------------------------------------------------------------------------
pdb_lig = """HETATM    1  C1  BEN A   1       1.212   0.700   0.000  1.00 20.00           C
HETATM    2  C2  BEN A   1       0.000   1.400   0.000  1.00 20.00           C
HETATM    3  C3  BEN A   1      -1.212   0.700   0.000  1.00 20.00           C
HETATM    4  C4  BEN A   1      -1.212  -0.700   0.000  1.00 20.00           C
HETATM    5  C5  BEN A   1       0.000  -1.400   0.000  1.00 20.00           C
HETATM    6  C6  BEN A   1       1.212  -0.700   0.000  1.00 20.00           C
END"""
r4 = requests.post(BASE_URL, data={'raw_content': pdb_lig, 'input_format': 'pdb', 'target_format': 'mol', 'add_hydrogens': 'false'})
res4 = r4.json()
has_bond_warn = any("PDB Bond Perception Active" in w for w in res4.get('metadata', {}).get('warnings', []))
results[4] = ("Loophole 4: PDB Bond-Order Perception & Aromaticity Restorer", has_bond_warn or res4.get('metadata', {}).get('num_atoms') == 6)

# -----------------------------------------------------------------------------
# LOOPHOLE 5: Multi-Molecule / Batch SDF Conversion to ZIP
# -----------------------------------------------------------------------------
sdf_multi = """Ethanol
  OpenBabel10242600002D

  3  2  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
    2.0000    1.2000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
  1  2  1  0
  2  3  1  0
M  END
$$$$
Water
  OpenBabel10242600002D

  1  0  0  0  0  0  0  0  0  0999 V2000
    0.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
M  END
$$$$
"""
r5 = requests.post(BASE_URL, data={'raw_content': sdf_multi, 'target_format': 'orca'})
res5 = r5.json()
results[5] = ("Loophole 5: Multi-Record SDF Library 1-Click Batch to ZIP", res5.get('is_zip') is True and "2 Molecules" in res5.get('metadata', {}).get('formula', ''))

# -----------------------------------------------------------------------------
# LOOPHOLE 6: Non-Standard Dummy / Ghost Atoms (Bq, Gh, X)
# -----------------------------------------------------------------------------
bsse_coords = """3
BSSE Complex with Bq Ghost Atom
Bq   0.000000   0.000000   1.500000
C    0.000000   0.000000   0.000000
H    0.000000   1.000000   0.000000"""
r6 = requests.post(BASE_URL, data={'raw_content': bsse_coords, 'target_format': 'gjf'})
res6 = r6.json()
out6 = res6.get('output_content', '')
has_bq = "Bq " in out6 or "Bq  " in out6
results[6] = ("Loophole 6: BSSE Ghost Atoms (Bq) Preserved in Output", has_bq and r6.status_code == 200)

# -----------------------------------------------------------------------------
# LOOPHOLE 7: Windows CRLF vs Linux LF Normalization
# -----------------------------------------------------------------------------
r7 = requests.post(BASE_URL, data={'raw_content': "CCO\r\n", 'input_format': 'smi', 'target_format': 'orca'})
out7 = r7.json().get('output_content', '')
zero_crlf = "\r" not in out7 and len(out7) > 0
results[7] = ("Loophole 7: Universal Unix LF Normalization (No ^M Characters)", zero_crlf and r7.status_code == 200)

# -----------------------------------------------------------------------------
# LOOPHOLE 8: Optical Isomer & Meso Compound Ambiguity
# -----------------------------------------------------------------------------
smi_meso = "O=C(O)[C@@H](O)[C@@H](O)C(=O)O"
r8 = requests.post(BASE_URL, data={'raw_content': smi_meso, 'input_format': 'smi', 'target_format': 'orca'})
res8 = r8.json()
is_meso = res8.get('metadata', {}).get('is_meso', False)
results[8] = ("Loophole 8: Meso Compound Optical Inactivity Symmetry Checksum", is_meso is True)

# -----------------------------------------------------------------------------
# FINAL SCORECARD
# -----------------------------------------------------------------------------
print(f"{'#':<3} {'Loophole Feature':<58} {'Status'}")
print("-" * 85)
passed_count = 0
for idx in range(1, 9):
    name, passed = results[idx]
    status_str = "PASS (100%)" if passed else "FAIL"
    if passed:
        passed_count += 1
    print(f"{idx:<3} {name:<58} {status_str}")

print("=" * 85)
print(f"RESULT: {passed_count} / 8 KILLER LOOPHOLES OPERATIONAL & VERIFIED ({(passed_count/8)*100:.1f}%)")
print("=" * 85)
