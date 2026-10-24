"""
Senior Scientist Audit Suite (BIOVIA / Schrödinger / PyMOL Standard)
Probing deep chemical, crystallographic, and quantum invariants.
"""
import os, sys, math, json, io, requests
sys.stdout.reconfigure(encoding='utf-8')
BASE = "http://127.0.0.1:8000/api/convert/"

def post(data=None, files=None):
    r = requests.post(BASE, data=data or {}, files=files, timeout=60)
    try: return r.status_code, r.json()
    except Exception: return r.status_code, {"_text": r.text[:300]}

print("=" * 80)
print("SENIOR SCIENTIFIC BENCHMARK & VULNERABILITY AUDIT")
print("=" * 80)

# Probe 1: [1.1.1]Propellane (Inverted bridgehead carbon tetrahedral geometry)
print("\n[PROBE 1] Highly Strained [1.1.1]Propellane (C5H6)...")
sc, j = post({"raw_content": "C1C23CC12C3", "target_format": "xyz"})
print(f"Status: {sc}, Formula: {j.get('metadata', {}).get('formula')}, Atoms: {j.get('metadata', {}).get('num_atoms')}")
assert sc == 200 and j.get('metadata', {}).get('num_atoms') == 11 # 5 C + 6 H

# Probe 2: Azulene (Non-benzenoid 5-7 fused aromatic ring system)
print("\n[PROBE 2] Azulene C10H8 (Aromaticity Perception across 5-7 fused rings)...")
sc, j = post({"raw_content": "c1cccc2cccc-2c1", "target_format": "orca"})
print(f"Status: {sc}, Formula: {j.get('metadata', {}).get('formula')}, Atoms: {j.get('metadata', {}).get('num_atoms')}")
assert sc == 200 and j.get('metadata', {}).get('formula') == 'C10H8'

# Probe 3: Pentavalent Nitro Group representation (N(=O)=O vs [N+](=O)[O-])
print("\n[PROBE 3] Pentavalent Nitro Group SMILES Sanitization (c1ccc(cc1)N(=O)=O)...")
sc, j = post({"raw_content": "c1ccc(cc1)N(=O)=O", "target_format": "xyz"})
print(f"Status: {sc}, Formula: {j.get('metadata', {}).get('formula')}, Parse Error?: {j.get('error')}")

# Probe 4: TEMPO Free Radical (Nitroxide Doublet Radical)
print("\n[PROBE 4] TEMPO Stable Nitroxide Free Radical (Doublet Spin Multiplicity)...")
sc, j = post({"raw_content": "CC1(C)CCCC(C)(C)N1[O]", "target_format": "orca"})
print(f"Status: {sc}, Mult: {j.get('metadata', {}).get('suggested_mult')}, State: {j.get('metadata', {}).get('electronic_state')}")

# Probe 5: Heavy Element (Gold / Platinum) with Pople Basis Set Check
print("\n[PROBE 5] Heavy Transition Metal (Cisplatin Pt) with 6-31G* basis set...")
cisplatin_xyz = "11\nCisplatin\nPt 0 0 0\nCl 2.3 0 0\nCl 0 2.3 0\nN -2.0 0 0\nN 0 -2.0 0\nH -2.3 0.8 0.5\nH -2.3 -0.8 0.5\nH -2.3 0 -0.9\nH 0.8 -2.3 0.5\nH -0.8 -2.3 0.5\nH 0 -2.3 -0.9"
sc, j = post({"raw_content": cisplatin_xyz, "target_format": "gjf", "basis": "6-31G*", "method": "B3LYP"})
warns = j.get('metadata', {}).get('warnings', [])
print(f"Status: {sc}, Basis Warnings: {warns}")

# Probe 6: PDB Insertion Codes (Antibody CDR loops e.g. Residue 100A, 100B)
print("\n[PROBE 6] PDB with Insertion Codes (100A, 100B)...")
pdb_insertion = """ATOM      1  N   ALA A 100      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A 100      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A 100      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A 100      11.000  12.000  10.000  1.00 20.00           O
ATOM      5  N   GLY A 100A     12.500  11.000  10.000  1.00 20.00           N
ATOM      6  CA  GLY A 100A     13.000  12.000  10.000  1.00 20.00           C
ATOM      7  C   GLY A 100A     14.000  12.000  10.000  1.00 20.00           C
ATOM      8  O   GLY A 100A     14.500  11.000  10.000  1.00 20.00           O
ATOM      9  N   SER A 100B     14.500  13.000  10.000  1.00 20.00           N
ATOM     10  CA  SER A 100B     15.500  13.000  10.000  1.00 20.00           C
END
"""
sc, j = post({"raw_content": pdb_insertion, "input_format": "pdb", "target_format": "xyz"})
print(f"Status: {sc}, Formula: {j.get('metadata', {}).get('formula')}, Atoms: {j.get('metadata', {}).get('num_atoms')}")

# Probe 7: Carbon Monoxide (CO) Formal Charge Balance in InChI / SMILES
print("\n[PROBE 7] Carbon Monoxide [C-]#[O+] (Zwitterionic Dipole)...")
sc, j = post({"raw_content": "[C-]#[O+]", "target_format": "orca"})
print(f"Status: {sc}, Charge: {j.get('metadata', {}).get('detected_charge')}, Mult: {j.get('metadata', {}).get('suggested_mult')}")

print("\n" + "=" * 80)
