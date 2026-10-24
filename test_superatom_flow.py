from rdkit import Chem

SUPERATOMS = {
    "PH": "*c1ccccc1",              # phenyl ring
    "OME": "*OC",                   # methoxy
    "OET": "*OCC",                  # ethoxy
    "ME": "*C",                     # methyl
    "ET": "*CC",                    # ethyl
    "PR": "*CCC",                   # propyl
    "IPR": "*C(C)C",                # isopropyl
    "BU": "*CCCC",                  # butyl
    "TBU": "*C(C)(C)C",             # tert-butyl
    "AC": "*C(=O)C",                # acetyl
    "BOC": "*C(=O)OC(C)(C)C",       # tert-butoxycarbonyl
    "OTS": "*S(=O)(=O)c1ccc(C)cc1", # tosyl
    "OMS": "*S(=O)(=O)C",           # mesyl
    "TFA": "*C(=O)C(F)(F)F",        # trifluoroacetyl
    "BN": "*Cc1ccccc1",             # benzyl
    "OBN": "*OCc1ccccc1",           # benzyloxy
    "COOH": "*C(=O)O",              # carboxyl
    "NO2": "*[N+](=O)[O-]",         # nitro
    "NH2": "*N",                    # amino
    "OH": "*O",                     # hydroxyl
    "CN": "*C#N"                    # cyano
}

def build_with_superatoms():
    # Node 1 is Carbon (id='1')
    # Node 2 is Ph superatom (id='2', text='Ph')
    # Bond between 1 and 2
    nodes = {
        '1': {'elem': 6, 'charge': 0, 'superatom': None, 'x': 5.0, 'y': 5.0},
        '2': {'elem': 0, 'charge': 0, 'superatom': 'PH', 'x': 6.0, 'y': 5.0}
    }
    bonds = [{'start': '1', 'end': '2', 'order': 1}]

    # Step 1: separate real atoms from superatom nodes
    real_nodes = {nid: d for nid, d in nodes.items() if not d.get('superatom')}
    super_nodes = {nid: d for nid, d in nodes.items() if d.get('superatom')}

    rw = Chem.RWMol()
    id_to_idx = {}
    for nid, d in real_nodes.items():
        atom = Chem.Atom(d['elem'])
        idx = rw.AddAtom(atom)
        id_to_idx[nid] = idx

    # Step 2: For each superatom node, find its connection to real atoms
    for s_nid, s_data in super_nodes.items():
        group_smi = SUPERATOMS.get(s_data['superatom'].upper())
        if not group_smi:
            continue
        # Find which real node it connects to
        connected_real_nid = None
        for b in bonds:
            if b['start'] == s_nid and b['end'] in real_nodes:
                connected_real_nid = b['end']
                break
            elif b['end'] == s_nid and b['start'] in real_nodes:
                connected_real_nid = b['start']
                break
        
        if connected_real_nid is not None:
            group_mol = Chem.MolFromSmiles(group_smi)
            # Combine current rw and group_mol
            base_mol = rw.GetMol()
            combo = Chem.CombineMols(base_mol, group_mol)
            rw = Chem.RWMol(combo)
            
            # Find the dummy atom (*)
            dummy_idx = None
            dummy_neighbor = None
            for a in rw.GetAtoms():
                if a.GetAtomicNum() == 0:
                    dummy_idx = a.GetIdx()
                    neighs = a.GetNeighbors()
                    if neighs:
                        dummy_neighbor = neighs[0].GetIdx()
                    break
            
            if dummy_idx is not None and dummy_neighbor is not None:
                real_attach_idx = id_to_idx[connected_real_nid]
                rw.AddBond(real_attach_idx, dummy_neighbor, Chem.BondType.SINGLE)
                rw.RemoveAtom(dummy_idx)
                # Re-index id_to_idx if dummy was before real_attach_idx (dummy is always at the end in CombineMols, but safe check)

    res_mol = rw.GetMol()
    Chem.SanitizeMol(res_mol)
    print("Result SMILES:", Chem.MolToSmiles(res_mol))

build_with_superatoms()
