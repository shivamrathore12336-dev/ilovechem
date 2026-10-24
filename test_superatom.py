from rdkit import Chem
from rdkit.Chem import AllChem

# Kekulized SMILES with dummy atom (*) at attachment point
SUPERATOMS = {
    "PH": "*c1ccccc1",              # attachment at dummy * (neighbor is aromatic C)
    "OME": "*OC",
    "OET": "*OCC",
    "ME": "*C",
    "ET": "*CC",
    "PR": "*CCC",
    "IPR": "*C(C)C",
    "BU": "*CCCC",
    "TBU": "*C(C)(C)C",
    "AC": "*C(=O)C",
    "BOC": "*C(=O)OC(C)(C)C",
    "OTS": "*S(=O)(=O)c1ccc(C)cc1",
    "OMS": "*S(=O)(=O)C",
    "TFA": "*C(=O)C(F)(F)F",
    "BN": "*Cc1ccccc1",
    "OBN": "*OCc1ccccc1",
    "COOH": "*C(=O)O",
    "NO2": "*[N+](=O)[O-]",
    "NH2": "*N",
    "OH": "*O",
    "CN": "*C#N"
}

def attach_group(base_mol, attach_atom_idx, group_name):
    group_smi = SUPERATOMS.get(group_name.upper())
    if not group_smi:
        return base_mol
    group_mol = Chem.MolFromSmiles(group_smi)
    # Combine base_mol and group_mol
    combo = Chem.CombineMols(base_mol, group_mol)
    rw = Chem.RWMol(combo)
    
    # In combo, the dummy atom (*) has atomic num 0
    dummy_idx = None
    dummy_neighbor_idx = None
    
    for atom in rw.GetAtoms():
        if atom.GetAtomicNum() == 0:
            dummy_idx = atom.GetIdx()
            neighbors = atom.GetNeighbors()
            if neighbors:
                dummy_neighbor_idx = neighbors[0].GetIdx()
            break
            
    if dummy_idx is not None and dummy_neighbor_idx is not None:
        rw.AddBond(attach_atom_idx, dummy_neighbor_idx, Chem.BondType.SINGLE)
        rw.RemoveAtom(dummy_idx)
        
    res = rw.GetMol()
    Chem.SanitizeMol(res)
    return res

# Test: Benzene ring with a bond to Ph -> Biphenyl
benzene = Chem.MolFromSmiles("c1ccccc1")
biphenyl = attach_group(benzene, 0, "PH")
print("Biphenyl SMILES:", Chem.MolToSmiles(biphenyl))

# Test: Amine with Boc -> Boc-protected amine
amine = Chem.MolFromSmiles("CC(C)N")
boc_amine = attach_group(amine, 3, "BOC")
print("Boc-Amine SMILES:", Chem.MolToSmiles(boc_amine))
