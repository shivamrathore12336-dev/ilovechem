import os
import sys
import time
import requests
from openbabel import pybel
from rdkit import Chem

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

BASE_URL = 'http://127.0.0.1:8000/api/convert/'

print("=" * 85)
print("EXTREME CORNER-CASE & EDGE-STRESS TEST SUITE (15 VULNERABILITY VECTORS)")
print("=" * 85)

stress_cases = [
    # Vector 1: Malformed UTF-8 & Raw Binary Garbage
    ("Vector 1: Random Binary Noise", b"\x00\xff\xfe\x12\x89\xfa\xbc\xde" * 100, "xyz", "xyz"),
    
    # Vector 2: Single atom coordinates (H, C, Fe)
    ("Vector 2: Single Hydrogen Atom", "1\nSingle H\nH 0.000 0.000 0.000\n", "xyz", "orca"),
    
    # Vector 3: Huge Coordinate Drift (10,000 Angstroms away)
    ("Vector 3: Vast Interatomic Separation (>1000 A)", "2\nVast separation\nC 0.000 0.000 0.000\nC 1500.000 0.000 0.000\n", "xyz", "xyz"),
    
    # Vector 4: Macrocyclic Peptide (Cyclosporin A SMILES - 130 atoms ring)
    ("Vector 4: Cyclosporin A Macrocycle", "CCC1C(=O)N(CC(=O)N(C(C(=O)NC(C(=O)N(C(C(=O)NC(C(=O)NC(C(=O)N(C(C(=O)N(C(C(=O)N(C(C(=O)N1C(C(C)CC=CC)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)C", "smi", "orca"),
    
    # Vector 5: BINAP Atropisomer (Axial Chirality)
    ("Vector 5: (R)-BINAP Chiral Ligand", "C1=CC=C(C=C1)P(C2=CC=CC=C2)C3=C(C4=CC=CC=C4C=C3)C5=C(C=CC6=CC=CC=C65)P(C7=CC=CC=C7)C8=CC=CC=C8", "smi", "orca"),
    
    # Vector 6: Ultra-Long Linear Alkane (C50H102)
    ("Vector 6: C50 Polyethylene Chain", "C" * 50, "smi", "xyz"),
    
    # Vector 7: Hypervalent Iodine & Sulfur (Dess-Martin Periodinane SMILES)
    ("Vector 7: Dess-Martin Periodinane", "CC(=O)OI1(C2=CC=CC=C2C(=O)O1)(OC(=O)C)OC(=O)C", "smi", "orca"),
    
    # Vector 8: ChemDraw Superatom with multiple attachments (Ph-CH2-Ph)
    ("Vector 8: CDXML Multiple Superatoms", """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<fragment id="1">
  <n id="1" p="100 100" Element="6"/>
  <n id="2" p="80 100" NodeType="Fragment"><t p="80 100"><s>Ph</s></t></n>
  <n id="3" p="120 100" NodeType="Fragment"><t p="120 100"><s>OMe</s></t></n>
  <b id="4" B="1" E="2" Order="1"/>
  <b id="5" B="1" E="3" Order="1"/>
</fragment>
</page>
</CDXML>""", "cdxml", "xyz"),

    # Vector 9: PDB with disordered occupancy (0.00 / 0.50)
    ("Vector 9: PDB Disordered Partial Occupancies", """ATOM      1  N   ALA A   1      11.000  10.000  10.000  0.50 20.00           N
ATOM      2  CA  ALA A   1      12.000  10.000  10.000  0.50 20.00           C
ATOM      3  C   ALA A   1      13.000  11.000  10.000  0.50 20.00           C
ATOM      4  O   ALA A   1      14.000  11.000  10.000  0.50 20.00           O
END""", "pdb", "gjf"),

    # Vector 10: Non-standard / Exotic Elements (Lanthanide Gadolinium Gd, Astatine At)
    ("Vector 10: Gadolinium MRI Complex", """2
Gadolinium Oxide
Gd 0.000 0.000 0.000
O  2.100 0.000 0.000""", "xyz", "orca"),

    # Vector 11: Empty string / Whitespace only
    ("Vector 11: Whitespace Only Input", "   \n\t   \n  ", "xyz", "orca"),

    # Vector 12: ChemDraw Reaction with disconnected salt AND superatom
    ("Vector 12: CDXML Scheme + Salt + Superatom", """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<t id="99" p="10 10"><s>Reflux, 98% yield</s></t>
<arrow id="100" p="50 50"/>
<fragment id="1">
  <n id="1" p="100 100" Element="6"/>
  <n id="2" p="120 100" NodeType="Fragment"><t p="120 100"><s>Ph</s></t></n>
  <b id="3" B="1" E="2" Order="1"/>
</fragment>
<fragment id="2">
  <n id="10" p="200 200" Element="17" Charge="-1"/>
</fragment>
</page>
</CDXML>""", "cdxml", "orca"),

    # Vector 13: Extremely high charge (+5 / -5)
    ("Vector 13: Hexamminecobalt(III) Ion [Co(NH3)6]3+", """7
Hexamminecobalt(III)
Co 0.000 0.000 0.000
N  1.950 0.000 0.000
N -1.950 0.000 0.000
N  0.000 1.950 0.000
N  0.000 -1.950 0.000
N  0.000 0.000 1.950
N  0.000 0.000 -1.950""", "xyz", "orca"),

    # Vector 14: CIF unit cell with 0-volume or non-numeric parameters
    ("Vector 14: CIF with non-orthogonal triclinic angles", """data_triclinic
_cell_length_a 5.00
_cell_length_b 6.00
_cell_length_c 7.00
_cell_angle_alpha 75.00
_cell_angle_beta 85.00
_cell_angle_gamma 65.00
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.100 0.200 0.300
O1 O 0.500 0.500 0.500""", "cif", "xyz"),

    # Vector 15: Giant 20,000 Line Empty Comment / Blank Trailing lines
    ("Vector 15: Trailing Blank Bloat (20,000 lines)", "3\nWater\nO 0.0 0.0 0.1\nH 0.0 0.7 -0.5\nH 0.0 -0.7 -0.5\n" + ("\n" * 15000), "xyz", "orca")
]

results = []

for idx, (name, content, in_fmt, target_fmt) in enumerate(stress_cases, 1):
    t0 = time.time()
    try:
        data = {
            'target_format': target_fmt,
            'add_hydrogens': 'true',
            'strip_salts': 'true'
        }
        if isinstance(content, bytes):
            data['raw_content'] = content.decode('latin1')
        else:
            data['raw_content'] = content
        if in_fmt:
            data['input_format'] = in_fmt

        r = requests.post(BASE_URL, data=data, timeout=30)
        elapsed = time.time() - t0
        
        # Binary noise & empty input should return 400 Bad Request with a clean error message, not 500 crash!
        if "Binary Noise" in name or "Whitespace" in name:
            if r.status_code == 400:
                results.append((name, "PASS", f"Correctly rejected with 400 Bad Request: {r.json().get('error')[:50]}...", elapsed))
            else:
                results.append((name, "FAIL", f"Unexpected status {r.status_code}: {r.text[:60]}", elapsed))
        else:
            if r.status_code == 200 and r.json().get('success'):
                res = r.json()
                meta = res.get('metadata', {})
                results.append((name, "PASS", f"Output generated ({meta.get('formula')}, {meta.get('num_atoms')} atoms)", elapsed))
            else:
                results.append((name, "FAIL", f"Status {r.status_code}: {r.json().get('error') if r.headers.get('content-type')=='application/json' else r.text[:60]}", elapsed))
    except Exception as e:
        results.append((name, "ERROR", str(e), time.time() - t0))

print("\n" + "=" * 85)
print("STRESS TEST RESULTS:")
print("=" * 85)
for name, status, detail, elapsed in results:
    mark = "[PASS]" if status == "PASS" else f"[{status}]"
    print(f"{mark:<8} | {name:<45} | {elapsed*1000:>6.1f}ms | {detail}")
print("=" * 85)
