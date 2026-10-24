"""
Exhaustive Master Registry of Computational Chemistry Formats
Grouping over 40+ production-grade scientific formats into clear, intuitive categories.
"""

COMPLETE_TARGET_FORMATS = {
    "Quantum Chemistry & Electronic Structure": [
        ("orca", "ORCA Input (.inp)", "DFT, TD-DFT, ab initio, coupled cluster"),
        ("gjf", "Gaussian Input (.gjf / .com)", "Gaussian 16/09 electronic structure"),
        ("qchem", "Q-Chem Input (.in)", "Quantum chemical calculations"),
        ("mop", "MOPAC Input (.mop)", "Semiempirical QM (PM6, PM7, RM1)"),
        ("nw", "NWChem Input (.nw)", "High-performance computational chemistry"),
        ("gam", "GAMESS (US) Input (.inp)", "Ab initio quantum chemistry"),
        ("fhiaims", "FHI-aims Input (.in)", "All-electron electronic structure"),
    ],
    "Materials Science, DFT & Crystallography": [
        ("poscar", "VASP POSCAR (Cartesian)", "Vienna Ab initio Simulation Package"),
        ("cif", "Crystallographic Info (.cif)", "X-ray crystal structure exchange"),
        ("pwi", "Quantum ESPRESSO (.pwi)", "Plane-wave pseudopotential DFT"),
        ("lmpdat", "LAMMPS Data (.data)", "Large-scale atomic/molecular massively parallel"),
        ("xtl", "XTL Crystal (.xtl)", "Crystal structure visualization"),
    ],
    "Drug Discovery, Docking & Biochemistry": [
        ("pdbqt", "AutoDock Vina (.pdbqt)", "Flexible docking with partial charges & torsions"),
        ("pdb", "Protein Data Bank (.pdb)", "Standard macromolecular coordinate standard"),
        ("mol2", "Tripos Mol2 (.mol2)", "SYBYL force field & docking format"),
        ("gro", "GROMACS Coordinate (.gro)", "Molecular dynamics simulation box"),
        ("pqr", "PQR (.pqr)", "PDB format with atomic charge & radius for APBS"),
    ],
    "2D Lab Sketching & Cheminformatics": [
        ("cdxml", "ChemDraw XML (.cdxml)", "PerkinElmer ChemDraw exchange format"),
        ("mol", "MDL Molfile (.mol)", "Standard 2D/3D connection table"),
        ("sdf", "Structure-Data File (.sdf)", "High-throughput multi-molecule library"),
        ("smi", "SMILES String (.smi)", "Simplified Molecular Input Line Entry"),
        ("inchi", "InChI (.txt)", "IUPAC International Chemical Identifier"),
        ("inchikey", "InChIKey (.txt)", "Fixed-length 27-character hashed key"),
        ("cml", "Chemical Markup Language (.cml)", "XML-based chemical structure format"),
    ],
    "Standard 3D Coordinates & Graphics": [
        ("xyz", "Standard XYZ (.xyz)", "Universal 3D Cartesian coordinates"),
        ("txyz", "TINKER XYZ (.xyz)", "TINKER molecular mechanics coordinate format"),
        ("svg", "2D Vector Diagram (.svg)", "Publication-quality scalable vector graphic"),
    ]
}
