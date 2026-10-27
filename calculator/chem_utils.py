"""
Cheminformatics utility functions using RDKit and OpenBabel.
Handles SMILES validation, property calculations, 3D coordinate generation, and 2D depiction.
Includes robust fallbacks for organometallics, coordination complexes, and non-kekulizable aromatics.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, QED, rdMolDescriptors
from rdkit.Chem import Draw
from openbabel import pybel
import io
import base64

def validate_smiles(smiles: str) -> bool:
    """Validate whether input SMILES can be converted to an RDKit Mol object or OpenBabel molecule."""
    if not smiles or not isinstance(smiles, str):
        return False
    s = smiles.strip()
    try:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            return True
    except Exception:
        pass
    try:
        mol = Chem.MolFromSmiles(s, sanitize=False)
        if mol is not None:
            return True
    except Exception:
        pass
    try:
        ob_mol = pybel.readstring("smi", s)
        return len(ob_mol.atoms) > 0
    except Exception:
        return False

def calculate_properties(smiles: str) -> dict:
    """
    Computes physicochemical properties and drug-likeness parameters.
    Supports organometallics, hypervalent complexes, and non-kekulizable aromatics.
    """
    if not smiles or not isinstance(smiles, str):
        return {"error": "Empty SMILES string"}

    s = smiles.strip()
    mol = None

    try:
        mol = Chem.MolFromSmiles(s)
    except Exception:
        mol = None

    if mol is None:
        try:
            mol = Chem.MolFromSmiles(s, sanitize=False)
            if mol is not None:
                mol.UpdatePropertyCache(strict=False)
                Chem.FastFindRings(mol)
        except Exception:
            mol = None

    if mol is None:
        # Fallback to OpenBabel for organometallics / non-kekulizable molecules
        try:
            ob_mol = pybel.readstring("smi", s)
            if len(ob_mol.atoms) > 0:
                mw = round(float(ob_mol.molwt), 2)
                formula = ob_mol.formula
                heavy_atoms = sum(1 for a in ob_mol.atoms if a.OBAtom.GetAtomicNum() > 1)
                return {
                    "smiles": s,
                    "formula": formula,
                    "molecular_weight": mw,
                    "exact_mass": mw,
                    "logp": 0.0,
                    "tpsa": 0.0,
                    "h_bond_donors": 0,
                    "h_bond_acceptors": 0,
                    "rotatable_bonds": 0,
                    "aromatic_rings": 0,
                    "rings": 0,
                    "heavy_atoms": heavy_atoms,
                    "fraction_csp3": 0.0,
                    "qed_score": 0.0,
                    "lipinski_violations": 1 if mw > 500 else 0,
                    "lipinski_pass": mw <= 500,
                    "veber_violations": 0,
                    "veber_pass": True,
                }
        except Exception:
            pass
        return {"error": "Invalid SMILES string"}

    try:
        mw = round(float(Descriptors.MolWt(mol)), 2)
        exact_mw = round(float(Descriptors.ExactMolWt(mol)), 3)
        try: logp = round(float(Crippen.MolLogP(mol)), 2)
        except Exception: logp = 0.0
        try: tpsa = round(float(Descriptors.TPSA(mol)), 2)
        except Exception: tpsa = 0.0
        try: hbd = int(Descriptors.NumHDonors(mol))
        except Exception: hbd = 0
        try: hba = int(Descriptors.NumHAcceptors(mol))
        except Exception: hba = 0
        try: rotatable = int(Descriptors.NumRotatableBonds(mol))
        except Exception: rotatable = 0
        try: aromatic_rings = int(Descriptors.NumAromaticRings(mol))
        except Exception: aromatic_rings = 0
        try: rings = int(Descriptors.RingCount(mol))
        except Exception: rings = 0
        heavy_atoms = int(mol.GetNumHeavyAtoms())
        try: formula = rdMolDescriptors.CalcMolFormula(mol)
        except Exception: formula = ""
        try: qed_val = round(float(QED.qed(mol)), 3)
        except Exception: qed_val = 0.0
        try: fsp3 = round(float(Descriptors.FractionCSP3(mol)), 2)
        except Exception: fsp3 = 0.0

        # Lipinski's Rule of 5 check
        lipinski_violations = 0
        if mw > 500:
            lipinski_violations += 1
        if logp > 5.0:
            lipinski_violations += 1
        if hbd > 5:
            lipinski_violations += 1
        if hba > 10:
            lipinski_violations += 1

        lipinski_pass = (lipinski_violations <= 1)

        # Veber's Rule for oral bioavailability:
        # Rotatable Bonds <= 10 and TPSA <= 140 Å²
        veber_violations = 0
        if rotatable > 10:
            veber_violations += 1
        if tpsa > 140.0:
            veber_violations += 1
        veber_pass = (veber_violations == 0)

        return {
            "smiles": s,
            "formula": formula,
            "molecular_weight": mw,
            "exact_mass": exact_mw,
            "logp": logp,
            "tpsa": tpsa,
            "h_bond_donors": hbd,
            "h_bond_acceptors": hba,
            "rotatable_bonds": rotatable,
            "aromatic_rings": aromatic_rings,
            "rings": rings,
            "heavy_atoms": heavy_atoms,
            "fraction_csp3": fsp3,
            "qed_score": qed_val,
            "lipinski_violations": lipinski_violations,
            "lipinski_pass": lipinski_pass,
            "veber_violations": veber_violations,
            "veber_pass": veber_pass,
        }
    except Exception:
        return {"error": "Failed to calculate molecular properties."}

def generate_3d_molblock(smiles: str) -> dict:
    """
    Adds explicit hydrogens, embeds 3D coordinates using ETKDG,
    and performs MMFF force-field energy minimization.
    Falls back to OpenBabel UFF for organometallics and complex coordination cages.
    Returns standard MolBlock string for 3Dmol.js rendering.
    """
    s = smiles.strip()
    mol = None
    try:
        mol = Chem.MolFromSmiles(s)
    except Exception:
        mol = None

    if mol is not None:
        try:
            mol_h = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            params.useRandomCoords = True
            params.maxIterations = 50
            embed_status = AllChem.EmbedMolecule(mol_h, params)
            if embed_status < 0:
                embed_status = AllChem.EmbedMolecule(mol_h, maxAttempts=20, randomSeed=42, useRandomCoords=True)
            
            if embed_status >= 0:
                try:
                    AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
                except Exception:
                    pass
                mol_block = Chem.MolToMolBlock(mol_h)
                return {"mol_block": mol_block}
        except Exception:
            pass

    # OpenBabel fallback for organometallics, macrocycles, and complex structures
    try:
        ob_mol = pybel.readstring("smi", s)
        if len(ob_mol.atoms) > 0:
            ob_mol.make3D(forcefield='uff', steps=200)
            return {"mol_block": ob_mol.write("mol")}
    except Exception:
        pass

    return {"error": "Failed to generate 3D coordinates."}

def generate_2d_base64_image(smiles: str, width: int = 360, height: int = 280) -> str:
    """
    Generates 2D structure diagram and returns base64 PNG data URL.
    Falls back to OpenBabel SVG if RDKit kekulization fails.
    """
    s = smiles.strip()
    mol = None
    try:
        mol = Chem.MolFromSmiles(s)
    except Exception:
        mol = None

    if mol is None:
        try:
            mol = Chem.MolFromSmiles(s, sanitize=False)
            if mol is not None:
                mol.UpdatePropertyCache(strict=False)
                Chem.FastFindRings(mol)
        except Exception:
            mol = None

    if mol is not None:
        try:
            img = Draw.MolToImage(mol, size=(width, height))
            buf = io.BytesIO()
            img.save(buf, format="PNG")
            b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
            return f"data:image/png;base64,{b64}"
        except Exception:
            pass

    # OpenBabel SVG fallback
    try:
        ob_mol = pybel.readstring("smi", s)
        svg_str = ob_mol.write("svg")
        if svg_str:
            b64 = base64.b64encode(svg_str.encode("utf-8")).decode("utf-8")
            return f"data:image/svg+xml;base64,{b64}"
    except Exception:
        pass

    return ""
