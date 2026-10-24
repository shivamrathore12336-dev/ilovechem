"""
Cheminformatics utility functions using RDKit.
Handles SMILES validation, property calculations, 3D coordinate generation, and 2D depiction.
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, QED, rdMolDescriptors
from rdkit.Chem import Draw
import io
import base64

def validate_smiles(smiles: str) -> bool:
    """Validate whether input SMILES can be converted to an RDKit Mol object."""
    if not smiles or not isinstance(smiles, str):
        return False
    try:
        mol = Chem.MolFromSmiles(smiles.strip())
        return mol is not None
    except Exception:
        return False

def calculate_properties(smiles: str) -> dict:
    """
    Computes physicochemical properties and drug-likeness parameters.
    """
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return {"error": "Invalid SMILES string"}

    try:
        mw = round(float(Descriptors.MolWt(mol)), 2)
        exact_mw = round(float(Descriptors.ExactMolWt(mol)), 3)
        logp = round(float(Crippen.MolLogP(mol)), 2)
        tpsa = round(float(Descriptors.TPSA(mol)), 2)
        hbd = int(Descriptors.NumHDonors(mol))
        hba = int(Descriptors.NumHAcceptors(mol))
        rotatable = int(Descriptors.NumRotatableBonds(mol))
        aromatic_rings = int(Descriptors.NumAromaticRings(mol))
        rings = int(Descriptors.RingCount(mol))
        heavy_atoms = int(mol.GetNumHeavyAtoms())
        formula = rdMolDescriptors.CalcMolFormula(mol)
        qed_val = round(float(QED.qed(mol)), 3)
        fsp3 = round(float(Descriptors.FractionCSP3(mol)), 2)

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

        return {
            "smiles": smiles.strip(),
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
        }
    except Exception as e:
        return {"error": str(e)}

def generate_3d_molblock(smiles: str) -> dict:
    """
    Adds explicit hydrogens, embeds 3D coordinates using ETKDG,
    and performs MMFF force-field energy minimization.
    Returns standard MolBlock string for 3Dmol.js rendering.
    """
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return {"error": "Invalid SMILES string"}

    try:
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        params.useRandomCoords = True
        params.maxIterations = 50
        embed_status = AllChem.EmbedMolecule(mol, params)
        if embed_status < 0:
            # Fallback to standard embed with random coords if ETKDGv3 fails
            embed_status = AllChem.EmbedMolecule(mol, maxAttempts=20, randomSeed=42, useRandomCoords=True)
            if embed_status < 0:
                return {"error": "Failed to generate 3D coordinates."}

        # Energy minimize using MMFF94
        try:
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
        except Exception:
            pass

        mol_block = Chem.MolToMolBlock(mol)
        return {"mol_block": mol_block}
    except Exception as e:
        return {"error": str(e)}

def generate_2d_base64_image(smiles: str, width: int = 360, height: int = 280) -> str:
    """
    Generates 2D structure diagram and returns base64 PNG data URL.
    """
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return ""
    try:
        img = Draw.MolToImage(mol, size=(width, height))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode("utf-8")
        return f"data:image/png;base64,{b64}"
    except Exception:
        return ""
