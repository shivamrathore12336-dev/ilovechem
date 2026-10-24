import os
import re
import io
import zipfile
import tempfile
import xml.etree.ElementTree as ET
from typing import List, Dict, Any, Tuple, Optional
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel
from openbabel import openbabel as ob

# Standard ChemDraw & Organic Chemistry Superatom Expansion Dictionary
# Maps common abbreviations to SMILES with attachment dummy (*) atom
SUPERATOM_EXPANSIONS = {
    "PH": "*c1ccccc1",              # Phenyl (C6H5)
    "OME": "*OC",                   # Methoxy (OCH3)
    "OET": "*OCC",                  # Ethoxy (OCH2CH3)
    "ME": "*C",                     # Methyl (CH3)
    "ET": "*CC",                    # Ethyl (C2H5)
    "PR": "*CCC",                   # Propyl (n-Pr)
    "IPR": "*C(C)C",                # Isopropyl (i-Pr)
    "BU": "*CCCC",                  # Butyl (n-Bu)
    "TBU": "*C(C)(C)C",             # tert-Butyl (t-Bu)
    "AC": "*C(=O)C",                # Acetyl (COCH3)
    "BOC": "*C(=O)OC(C)(C)C",       # tert-Butoxycarbonyl
    "OTS": "*S(=O)(=O)c1ccc(C)cc1", # Tosyl (p-toluenesulfonyl)
    "OMS": "*S(=O)(=O)C",           # Mesyl (methanesulfonyl)
    "TFA": "*C(=O)C(F)(F)F",        # Trifluoroacetyl
    "BN": "*Cc1ccccc1",             # Benzyl (CH2Ph)
    "OBN": "*OCc1ccccc1",           # Benzyloxy (OCH2Ph)
    "COOH": "*C(=O)O",              # Carboxyl
    "CO2ME": "*C(=O)OC",            # Methyl ester
    "CO2ET": "*C(=O)OCC",           # Ethyl ester
    "NO2": "*[N+](=O)[O-]",         # Nitro
    "NH2": "*N",                    # Amino
    "OH": "*O",                     # Hydroxyl
    "CN": "*C#N",                   # Cyano
    "TMS": "*[Si](C)(C)C",          # Trimethylsilyl
    "TBS": "*[Si](C)(C)C(C)(C)C",   # tert-Butyldimethylsilyl (TBDMS)
    "TBDMS": "*[Si](C)(C)C(C)(C)C", # tert-Butyldimethylsilyl
    "CBZ": "*C(=O)OCc1ccccc1",      # Carboxybenzyl
    "FMOC": "*C(=O)OCC1c2ccccc2-c3ccccc13" # Fluorenylmethyloxycarbonyl
}


def sanitize_and_extract_all_fragments(raw_text: str) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Parses full ChemDraw Reaction Schemes:
    - Filters out non-chemical elements (<arrow>, <t>, <graphic>, <bracketedgroup>)
    - Extracts EVERY molecule (<fragment>) on the page independently
    - Detects and tags superatom text abbreviations (Ph, OMe, Boc, OTs, Ac, etc.)
    - Auto-computes formula, molecular weight, and 2D topology for each
    Returns: (list_of_molecules_data, inspection_summary)
    """
    if len(raw_text) > 50 * 1024 * 1024:
        raise ValueError("CDXML file exceeds maximum allowed size (50 MB).")
    # Security Guard: Defuse XXE / Entity injection attacks (handles multiline bracketed DOCTYPE)
    clean_xml = re.sub(r'<!DOCTYPE\s+[^\[>]*(\[[^\]]*\])?\s*>', '', raw_text, flags=re.DOTALL | re.IGNORECASE)
    clean_xml = re.sub(r'<!ENTITY[^>]*>', '', clean_xml, flags=re.IGNORECASE)
    clean_xml = re.sub(r'<!DOCTYPE[^>]*>', '', clean_xml, flags=re.IGNORECASE)
    root = ET.fromstring(clean_xml)

    stats = {
        "text_labels_filtered": 0,
        "reaction_arrows_filtered": 0,
        "graphics_filtered": 0,
        "total_molecules_found": 0,
        "superatoms_detected": []
    }

    for t in root.iter('t'): stats["text_labels_filtered"] += 1
    for a in root.iter('arrow'): stats["reaction_arrows_filtered"] += 1
    for g in root.iter('graphic'): stats["graphics_filtered"] += 1
    for b in root.iter('bracketedgroup'): stats["graphics_filtered"] += 1

    fragments_data = []

    for frag_idx, frag in enumerate(root.iter('fragment')):
        nodes = {}
        bonds = []

        for n in frag.iter('n'):
            nid = n.get('id')
            if not nid: continue
            
            # Inspect if node has superatom text label (e.g. <s font="1">Ph</s>)
            node_type = n.get('NodeType', '')
            superatom_label = None
            text_el = n.find('.//s')
            if text_el is not None and text_el.text:
                candidate = text_el.text.strip().upper()
                if candidate in SUPERATOM_EXPANSIONS:
                    superatom_label = candidate
                    stats["superatoms_detected"].append(candidate)

            elem_attr = n.get('Element')
            elem = int(elem_attr) if elem_attr and elem_attr.isdigit() else 6
            charge = int(n.get('Charge', '0'))
            num_radicals = int(n.get('NumRadicals', '0'))
            radical_attr = n.get('Radical', '')
            if not num_radicals and radical_attr in ('Doublet', 'Radical', 'monoradical'):
                num_radicals = 1
            elif not num_radicals and radical_attr in ('Triplet', 'diradical'):
                num_radicals = 2

            p = n.get('p', '0 0').split()
            x = float(p[0]) / 20.0
            y = float(p[1]) / 20.0
            
            nodes[nid] = {
                'elem': elem,
                'charge': charge,
                'num_radicals': num_radicals,
                'superatom': superatom_label,
                'x': x,
                'y': y
            }

        for b in frag.iter('b'):
            bid_start = b.get('B')
            bid_end = b.get('E')
            order = int(b.get('Order', '1'))
            display = b.get('Display', '')
            if bid_start and bid_end and bid_start in nodes and bid_end in nodes:
                bonds.append({
                    'start': bid_start,
                    'end': bid_end,
                    'order': order,
                    'display': display
                })

        if len(nodes) > 0:
            # Estimate rough formula
            elem_counts = {}
            for d in nodes.values():
                if d.get('superatom'):
                    lbl = d['superatom']
                    elem_counts[lbl] = elem_counts.get(lbl, 0) + 1
                else:
                    sym = ob.GetSymbol(d['elem'])
                    elem_counts[sym] = elem_counts.get(sym, 0) + 1
            
            formula_str = "".join([f"{k}{elem_counts[k] if elem_counts[k]>1 else ''}" for k in sorted(elem_counts.keys())])
            
            fragments_data.append({
                'fragment_id': frag.get('id', str(frag_idx + 1)),
                'index': frag_idx,
                'nodes': nodes,
                'bonds': bonds,
                'atom_count': len(nodes),
                'formula': formula_str,
                'is_largest': False
            })

    # Flag largest molecule
    if fragments_data:
        sorted_by_size = sorted(fragments_data, key=lambda x: x['atom_count'], reverse=True)
        sorted_by_size[0]['is_largest'] = True

    stats["total_molecules_found"] = len(fragments_data)
    return fragments_data, stats


def build_3d_mol_from_fragment(frag_data: Dict[str, Any]) -> Optional[pybel.Molecule]:
    nodes = frag_data['nodes']
    bonds = frag_data['bonds']

    # Separate real chemical atoms from superatom shortcut nodes
    real_nodes = {nid: d for nid, d in nodes.items() if not d.get('superatom')}
    super_nodes = {nid: d for nid, d in nodes.items() if d.get('superatom')}

    rw = Chem.RWMol()
    id_to_idx = {}

    for nid, data in real_nodes.items():
        atom = Chem.Atom(data['elem'])
        if data['charge'] != 0:
            atom.SetFormalCharge(data['charge'])
        if data.get('num_radicals', 0) > 0:
            atom.SetNumRadicalElectrons(data['num_radicals'])
        idx = rw.AddAtom(atom)
        id_to_idx[nid] = idx

    bond_order_map = {
        1: Chem.BondType.SINGLE,
        2: Chem.BondType.DOUBLE,
        3: Chem.BondType.TRIPLE,
        4: Chem.BondType.AROMATIC
    }

    # Add bonds between real nodes
    for b in bonds:
        if b['start'] in real_nodes and b['end'] in real_nodes:
            idx1 = id_to_idx[b['start']]
            idx2 = id_to_idx[b['end']]
            btype = bond_order_map.get(b['order'], Chem.BondType.SINGLE)
            b_idx = rw.AddBond(idx1, idx2, btype) - 1

            display = b.get('display', '')
            if display == 'WedgeBegin':
                rw.GetBondWithIdx(b_idx).SetBondDir(Chem.BondDir.BEGINWEDGE)
            elif display == 'WedgedHashBegin':
                rw.GetBondWithIdx(b_idx).SetBondDir(Chem.BondDir.BEGINDASH)

    # Expand any attached superatom groups (Ph, OMe, Boc, etc.)
    for s_nid, s_data in super_nodes.items():
        group_smi = SUPERATOM_EXPANSIONS.get(s_data['superatom'])
        if not group_smi:
            continue
        
        # Locate the real atom to which this superatom node was connected
        connected_real_nid = None
        for b in bonds:
            if b['start'] == s_nid and b['end'] in real_nodes:
                connected_real_nid = b['end']
                break
            elif b['end'] == s_nid and b['start'] in real_nodes:
                connected_real_nid = b['start']
                break

        if connected_real_nid is not None:
            try:
                group_mol = Chem.MolFromSmiles(group_smi)
                if group_mol:
                    base_mol = rw.GetMol()
                    combo = Chem.CombineMols(base_mol, group_mol)
                    rw = Chem.RWMol(combo)

                    # Dummy atom (*) has atomic num 0
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
            except Exception:
                pass

    mol = rw.GetMol()

    try:
        Chem.SanitizeMol(mol, Chem.SANITIZE_ALL ^ Chem.SANITIZE_KEKULIZE)
    except Exception:
        try: Chem.SanitizeMol(mol)
        except Exception: pass

    # If only real nodes were present and coordinates match, assign 2D conformer
    if mol.GetNumAtoms() == len(real_nodes):
        conf = Chem.Conformer(mol.GetNumAtoms())
        for nid, data in real_nodes.items():
            atom_idx = id_to_idx[nid]
            conf.SetAtomPosition(atom_idx, (float(data.get('x', 0.0)), float(data.get('y', 0.0)), 0.0))
        mol.AddConformer(conf, assignId=True)

        try:
            Chem.AssignChiralTypesFromBondDirs(mol, confId=0)
            Chem.SetDoubleBondNeighborDirections(mol, conf)
        except Exception:
            pass

    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except Exception:
        pass

    mol = Chem.AddHs(mol, addCoords=True)
    try:
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    except Exception:
        pass

    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.enforceChirality = True
    params.useExpTorsionAnglePrefs = True
    params.useBasicKnowledge = True
    
    embed_ok = AllChem.EmbedMolecule(mol, params)
    if embed_ok < 0:
        AllChem.EmbedMolecule(mol, randomSeed=42)

    try: AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    except Exception:
        try: AllChem.UFFOptimizeMolecule(mol, maxIters=500)
        except Exception: pass

    mol_block = Chem.MolToMolBlock(mol)
    py_mol = pybel.readstring("mol", mol_block)
    return py_mol


def create_zip_archive_of_all_fragments(fragments_data: List[Dict[str, Any]], target_format: str, export_fn, qm_options: dict) -> Tuple[bytes, str]:
    """
    Converts EVERY molecule in the reaction scheme into the target format
    and bundles them into a clean, ready-to-download .ZIP archive.
    """
    mem_zip = io.BytesIO()
    with zipfile.ZipFile(mem_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for idx, frag in enumerate(fragments_data):
            mol = build_3d_mol_from_fragment(frag)
            if mol:
                res = export_fn(mol, target_format, qm_options)
                out_content, fname = res[0], res[1]
                base_name, ext = os.path.splitext(fname)
                clean_formula = re.sub(r'[^a-zA-Z0-9_]', '', frag.get('formula', 'mol'))
                clean_ext = re.sub(r'[^a-zA-Z0-9\.]', '', ext)
                zip_filename = f"molecule_{idx + 1}_{clean_formula}{clean_ext}"
                zf.writestr(zip_filename, out_content.encode('utf-8'))

    mem_zip.seek(0)
    return mem_zip.getvalue(), f"reaction_scheme_all_molecules.zip"
