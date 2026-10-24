from calculator.scientific_engine import PDBBioSanitizer
from openbabel import pybel

# PDB line without columns 77-78 where atom is CA in ALA
legacy_pdb = """ATOM      1  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00
ATOM      2  CB  ALA A   1      12.000  10.000  10.000  1.00 20.00
ATOM      3  N   ALA A   1      10.000  10.000  10.000  1.00 20.00
END"""

sanitized = PDBBioSanitizer.sanitize_pdb(legacy_pdb)
print("Sanitized PDB lines:")
for l in sanitized.splitlines():
    print(repr(l))

mol = pybel.readstring("pdb", sanitized)
atomic_nums = [a.OBAtom.GetAtomicNum() for a in mol.atoms]
print("Atomic numbers parsed:", atomic_nums)
has_calcium = 20 in atomic_nums
has_carbon = 6 in atomic_nums
print(f"Has Calcium (Ca=20): {has_calcium} | Has Carbon (C=6): {has_carbon}")
assert not has_calcium, "Error: Alpha carbon was misidentified as Calcium!"
print("SUCCESS: Element column correctly enforced as Carbon (atomic #6)!")
