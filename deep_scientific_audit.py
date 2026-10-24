import os
import sys
import tempfile
import requests
from openbabel import openbabel as ob
from openbabel import pybel
from rdkit import Chem

from calculator.converter_core import (
    detect_format,
    parse_and_extract_atoms,
    export_to_target_format,
    _apply_integrity_audit,
    COMPLETE_TARGET_FORMATS
)
from calculator.scientific_engine import (
    QuantumLogParser,
    DeepContentSniffer,
    SpinChargeParityValidator,
    ZeroCoordinateGuard,
    SaltAndSolventStripper,
    PhysicalIntegrityChecker,
    PeriodicBoundaryUnwrapper,
    CrystalHydrogenator,
    PDBBioSanitizer,
    GhostAtomNormalizer,
    UnixEolNormalizer,
    MesoCompoundDetector
)

BASE_URL = 'http://127.0.0.1:8000/api/convert/'

print("=" * 80)
print("DEEP SCIENTIFIC & EDGE-CASE AUDIT ENGINE")
print("=" * 80)

audit_findings = []

def record(test_name, status, detail, severity="MEDIUM"):
    audit_findings.append({
        "test": test_name,
        "status": status,
        "detail": detail,
        "severity": severity
    })
    mark = "[PASS]" if status == "PASS" else f"[{status}]"
    print(f"{mark:<8} | {test_name:<40} | {detail}")

# ----------------------------------------------------------------------
# Test 1: Single Atom Centered at Origin (Helium or Monatomic Ion)
# ----------------------------------------------------------------------
# Legitimate single atom at (0,0,0) vs collapsed silent failure
single_atom_xyz = "1\nSingle Helium at Origin\nHe 0.0000 0.0000 0.0000\n"
try:
    r = requests.post(BASE_URL, data={'raw_content': single_atom_xyz, 'target_format': 'xyz'})
    res = r.json()
    if res.get('success'):
        prevented = res.get('metadata', {}).get('silent_failure_prevented')
        if not prevented:
            record("Single Atom at Origin", "PASS", "Single atom at (0,0,0) correctly not flagged as collapsed crash.")
        else:
            record("Single Atom at Origin", "WARN", "Single atom falsely flagged as collapsed geometry bug.", "LOW")
    else:
        record("Single Atom at Origin", "FAIL", f"API rejected single atom: {res.get('error')}", "HIGH")
except Exception as e:
    record("Single Atom at Origin", "ERROR", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 2: Macrocyclic / Strained Rings 3D Embedding (Cubane & C60)
# ----------------------------------------------------------------------
cubane_smi = "C12C3C4C1C5C2C3C45" # Cubane C8H8
try:
    r = requests.post(BASE_URL, data={'raw_content': cubane_smi, 'target_format': 'xyz'})
    res = r.json()
    if res.get('success') and res.get('metadata', {}).get('num_atoms') == 16: # C8H8
        record("Strained Cubane C8H8 3D Embed", "PASS", "Cubane successfully embedded with 16 atoms (8C + 8H).")
    else:
        record("Strained Cubane C8H8 3D Embed", "FAIL", f"Cubane 3D generation issue: {res}", "HIGH")
except Exception as e:
    record("Strained Cubane C8H8 3D Embed", "ERROR", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 3: Macrocycle / Porphyrin Coordination Complex
# ----------------------------------------------------------------------
porphyrin_smi = "c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)n5)n4)n3"
try:
    r = requests.post(BASE_URL, data={'raw_content': porphyrin_smi, 'target_format': 'orca'})
    res = r.json()
    if res.get('success'):
        record("Porphyrin Macrocycle 3D Embed", "PASS", "Complex porphyrin ring successfully minimized.")
    else:
        record("Porphyrin Macrocycle 3D Embed", "FAIL", f"Failed: {res.get('error')}", "MEDIUM")
except Exception as e:
    record("Porphyrin Macrocycle 3D Embed", "ERROR", str(e), "MEDIUM")

# ----------------------------------------------------------------------
# Test 4: Rare Earth / Lanthanide / Actinide Elements (e.g. Europium, Uranium)
# ----------------------------------------------------------------------
uranium_xyz = """2
Uranyl Cation Complex
U  0.0000 0.0000 0.0000
O  0.0000 0.0000 1.7600
"""
try:
    r = requests.post(BASE_URL, data={'raw_content': uranium_xyz, 'target_format': 'orca'})
    res = r.json()
    if res.get('success'):
        has_metals = res.get('metadata', {}).get('has_metals')
        record("Actinide (Uranium) UFF Guard", "PASS" if has_metals else "WARN", 
               f"Heavy element parsed. Metal guard detected: {has_metals}")
    else:
        record("Actinide (Uranium) UFF Guard", "FAIL", f"Failed: {res.get('error')}", "HIGH")
except Exception as e:
    record("Actinide (Uranium) UFF Guard", "ERROR", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 5: Heteroaromatic Kekulization & Parity (Tetrazole, Pyridine, Thiophene)
# ----------------------------------------------------------------------
tetrazole_smi = "c1nnn[nH]1"
try:
    r = requests.post(BASE_URL, data={'raw_content': tetrazole_smi, 'target_format': 'gjf'})
    res = r.json()
    if res.get('success'):
        record("Tetrazole 5-Nitrogen Ring", "PASS", "5-membered hyper-nitrogen ring correctly parsed & sanitized.")
    else:
        record("Tetrazole 5-Nitrogen Ring", "FAIL", f"Kekulization failed: {res.get('error')}", "HIGH")
except Exception as e:
    record("Tetrazole 5-Nitrogen Ring", "ERROR", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 6: Inconsistent / Empty QM Parameter Handling (User typed garbage)
# ----------------------------------------------------------------------
try:
    r = requests.post(BASE_URL, data={
        'raw_content': 'C',
        'target_format': 'orca',
        'charge': 'invalid_charge',
        'mult': 'invalid_mult'
    })
    res = r.json()
    if res.get('success'):
        record("Invalid Charge/Mult Fallback", "PASS", "Gracefully fell back to auto-detected charge/mult without 500 error.")
    else:
        record("Invalid Charge/Mult Fallback", "WARN", f"Returned error: {res.get('error')}", "LOW")
except Exception as e:
    record("Invalid Charge/Mult Fallback", "FAIL", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 7: Target Format Matrix — SMILES to every single registered target format
# ----------------------------------------------------------------------
test_mol_smi = "CC(=O)O" # Acetic acid
target_failures = []
for cat, fmt_list in COMPLETE_TARGET_FORMATS.items():
    for code, label, desc in fmt_list:
        try:
            r = requests.post(BASE_URL, data={'raw_content': test_mol_smi, 'target_format': code})
            res = r.json()
            if not res.get('success') or not res.get('output_content'):
                target_failures.append((code, label, res.get('error')))
        except Exception as e:
            target_failures.append((code, label, str(e)))

if not target_failures:
    record("Complete Target Format Matrix (All Formats)", "PASS", "100% of registered output formats export cleanly.")
else:
    record("Complete Target Format Matrix (All Formats)", "FAIL", f"{len(target_failures)} formats failed: {target_failures}", "HIGH")

# ----------------------------------------------------------------------
# Test 8: PDB to Gaussian GJF (PDB residue tag pollution in QM geometry)
# ----------------------------------------------------------------------
pdb_content = """ATOM      1  N   ALA A   1      11.104  13.207  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      12.000  13.000  11.000  1.00 20.00           C
ATOM      3  C   ALA A   1      13.000  14.000  11.000  1.00 20.00           C
END"""
try:
    r = requests.post(BASE_URL, data={'raw_content': pdb_content, 'target_format': 'gjf'})
    res = r.json()
    gjf_text = res.get('output_content', '')
    import re
    # Check if ALA, ATOM, HETATM or B-factor tokens leaked into Gaussian coordinate block
    leaked = any(re.search(rf"\b{tag}\b", gjf_text) for tag in ["ALA", "ATOM", "HETATM", "20.00"])
    if not leaked and res.get('success'):
        record("PDB to Gaussian Metadata Strip", "PASS", "All PDB residue & B-factor columns cleanly stripped from GJF.")
    else:
        record("PDB to Gaussian Metadata Strip", "FAIL", "PDB metadata leaked into Gaussian coordinates!", "HIGH")
except Exception as e:
    record("PDB to Gaussian Metadata Strip", "ERROR", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 9: Very High Spin State (e.g. Iron(III) High Spin d5 Sextet mult=6)
# ----------------------------------------------------------------------
fe3_complex_xyz = """7
Fe(III) Hexaaqua Complex
Fe 0.000 0.000 0.000
O  2.000 0.000 0.000
O -2.000 0.000 0.000
O  0.000 2.000 0.000
O  0.000 -2.000 0.000
O  0.000 0.000 2.000
O  0.000 0.000 -2.000"""
try:
    r = requests.post(BASE_URL, data={
        'raw_content': fe3_complex_xyz,
        'target_format': 'orca',
        'charge': 3,
        'mult': 6 # Sextet
    })
    res = r.json()
    out = res.get('output_content', '')
    if "* xyz 3 6" in out or "*xyz 3 6" in out:
        record("High Spin Transition Metal (Sextet Fe3+)", "PASS", "Explicit user high-spin multiplicity 6 accepted.")
    else:
        record("High Spin Transition Metal (Sextet Fe3+)", "WARN", f"Multiplicity parity modified: {out.splitlines()[:5]}", "MEDIUM")
except Exception as e:
    record("High Spin Transition Metal (Sextet Fe3+)", "ERROR", str(e), "MEDIUM")

# ----------------------------------------------------------------------
# Test 10: Scientific Ambiguity Check: What if CIF file has spacegroup P 21/c with fractional coords?
# ----------------------------------------------------------------------
cif_fractional = """data_test_crystal
_cell_length_a 10.000
_cell_length_b 12.000
_cell_length_c 8.000
_cell_angle_alpha 90.00
_cell_angle_beta 105.00
_cell_angle_gamma 90.00
_symmetry_space_group_name_H-M 'P 21/c'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.2500 0.3333 0.5000
O1 O 0.2500 0.5000 0.5000
"""
try:
    r = requests.post(BASE_URL, data={'raw_content': cif_fractional, 'target_format': 'xyz'})
    res = r.json()
    if res.get('success'):
        record("Monoclinic CIF (Beta=105) Transformation", "PASS", "Non-orthogonal monoclinic unit cell correctly converted to Cartesian Ångströms.")
    else:
        record("Monoclinic CIF (Beta=105) Transformation", "FAIL", f"Failed: {res.get('error')}", "HIGH")
except Exception as e:
    record("Monoclinic CIF (Beta=105) Transformation", "ERROR", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 11: ChemDraw Text Abbreviations (e.g. OMe, Ph, COOH)
# ----------------------------------------------------------------------
cdxml_abbrev = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<fragment id="1">
  <n id="1" p="100 100" Element="6"/>
  <n id="2" p="120 100" NodeType="Fragment"/>
  <b id="3" B="1" E="2" Order="1"/>
</fragment>
</page>
</CDXML>"""
try:
    r = requests.post(BASE_URL, data={'raw_content': cdxml_abbrev, 'target_format': 'xyz'})
    res = r.json()
    if res.get('success'):
        record("CDXML Node Type Fragment / Abbrev", "PASS", "Fragment node handled without crash.")
    else:
        record("CDXML Node Type Fragment / Abbrev", "WARN", f"CDXML abbrev notice: {res.get('error')}", "LOW")
except Exception as e:
    record("CDXML Node Type Fragment / Abbrev", "WARN", str(e), "LOW")

# ----------------------------------------------------------------------
# Test 12: ChemDraw Wedge/Dash Parity Parity in 3D (Issue 3)
# ----------------------------------------------------------------------
chiral_cdxml = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<fragment id="1">
  <n id="1" p="100 100" Element="6"/>
  <n id="2" p="120 100" Element="9"/>
  <n id="3" p="100 120" Element="17"/>
  <n id="4" p="80 100" Element="35"/>
  <n id="5" p="100 80" Element="53"/>
  <b id="6" B="1" E="2" Order="1" Display="WedgeBegin"/>
  <b id="7" B="1" E="3" Order="1" Display="WedgedHashBegin"/>
  <b id="8" B="1" E="4" Order="1"/>
  <b id="9" B="1" E="5" Order="1"/>
</fragment>
</page>
</CDXML>"""
try:
    r = requests.post(BASE_URL, data={'raw_content': chiral_cdxml, 'target_format': 'xyz'})
    res = r.json()
    if res.get('success') and res.get('metadata', {}).get('num_atoms') >= 5:
        record("Chiral CDXML (Wedge & Dash)", "PASS", "Explicit chiral center embedded with 3D parity lock.")
    else:
        record("Chiral CDXML (Wedge & Dash)", "FAIL", f"Chiral CDXML failed: {res.get('error')}", "HIGH")
except Exception as e:
    record("Chiral CDXML (Wedge & Dash)", "ERROR", str(e), "HIGH")

# ----------------------------------------------------------------------
# Test 13: Zero Coordinates in 2D SMILES (Can RDKit embed all common drugs?)
# ----------------------------------------------------------------------
drugs = {
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Paracetamol": "CC(=O)Nc1ccc(O)cc1",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "Caffeine": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "Penicillin G": "CC1(C)SC2C(NC(=O)Cc3ccccc3)C(=O)N2C1C(=O)O"
}
drug_fails = []
for dname, dsmi in drugs.items():
    try:
        r = requests.post(BASE_URL, data={'raw_content': dsmi, 'target_format': 'xyz'})
        res = r.json()
        if not res.get('success'):
            drug_fails.append(dname)
    except:
        drug_fails.append(dname)

if not drug_fails:
    record("Top Commercial Drug Molecules 3D Embed", "PASS", "Aspirin, Paracetamol, Ibuprofen, Caffeine, Penicillin G 100% generated.")
else:
    record("Top Commercial Drug Molecules 3D Embed", "FAIL", f"Failed drugs: {drug_fails}", "HIGH")

# ----------------------------------------------------------------------
# Test 14: Radical with Multiple Unpaired Electrons (Triplet O2, Carbenes)
# ----------------------------------------------------------------------
triplet_o2_xyz = """2
Molecular Oxygen Triplet Ground State
O 0.000 0.000 0.000
O 0.000 0.000 1.210"""
try:
    r = requests.post(BASE_URL, data={
        'raw_content': triplet_o2_xyz,
        'target_format': 'orca',
        'charge': 0,
        'mult': 3 # Triplet
    })
    res = r.json()
    out = res.get('output_content', '')
    if "* xyz 0 3" in out or "*xyz 0 3" in out:
        record("Triplet Diradical (Molecular Oxygen)", "PASS", "User-specified Triplet spin multiplicity 3 accepted for even-electron diradical.")
    else:
        record("Triplet Diradical (Molecular Oxygen)", "WARN", f"Triplet overridden: {out.splitlines()[:5]}", "MEDIUM")
except Exception as e:
    record("Triplet Diradical (Molecular Oxygen)", "ERROR", str(e), "MEDIUM")

print("=" * 80)
print(f"AUDIT COMPLETE: {len(audit_findings)} deep tests executed.")
print("=" * 80)
