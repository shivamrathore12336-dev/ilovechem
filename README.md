# iLoveChem 🧪💻
### Universal Computational Chemistry & CADD File Converter Studio

[![Live Web App](https://img.shields.io/badge/Live_App-ilovechem.onrender.com-blue?style=for-the-badge&logo=render)](https://ilovechem.onrender.com/)
[![Python 3.11](https://img.shields.io/badge/Python-3.11-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://www.python.org/)
[![Django](https://img.shields.io/badge/Django-6.1-092E20?style=for-the-badge&logo=django&logoColor=white)](https://www.djangoproject.com/)
[![RDKit](https://img.shields.io/badge/RDKit-Cheminformatics-238636?style=for-the-badge)](https://www.rdkit.org/)
[![OpenBabel](https://img.shields.io/badge/OpenBabel-Chemical_Toolbox-blueviolet?style=for-the-badge)](https://openbabel.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg?style=for-the-badge)](https://opensource.org/licenses/MIT)

**iLoveChem** is a free, fast, and minimalist computational chemistry web platform designed to streamline molecular file prep for drug discovery, quantum mechanics (DFT), molecular docking, and crystallography.

Inspired by the clean simplicity of **iLovePDF**, iLoveChem eliminates the friction of broken spin multiplicities, missing polar hydrogens, and manual command-line format conversions.

🌐 **Live Production App:** [https://ilovechem.onrender.com/](https://ilovechem.onrender.com/)  
📥 **Essential Software Downloads Hub:** [https://ilovechem.onrender.com/downloads/](https://ilovechem.onrender.com/downloads/)

---

## ⚡ Core Features

- **Universal 20+ Format Interconversion:** PDB, PDBQT, SDF, MOL2, CIF, XYZ, CDXML, SMILES, InChI, Gaussian (`.com`/`.gjf`), ORCA (`.inp`), and quantum logs (`.log`/`.out`).
- **Automated 3D Conformers & Minimization:** Embeds 3D coordinates using **ETKDG v3** and relaxes structures with **MMFF94 / UFF** force-field minimization.
- **ORCA Input Generator (`.inp`):** Calculates valence electrons and detects odd-electron radicals to automatically balance singlet/doublet spin states, preventing runtime DFT crashes.
- **Gaussian Route Card Builder (`.com` / `.gjf`):** Generates ready-to-run calculation decks (`# opt freq b3lyp/6-31g(d)`) with automatic net molecular charge resolution.
- **1-Click AutoDock Vina Preparation (`.pdbqt`):** Computes Gasteiger partial charges, merges non-polar hydrogens, and configures rotatable torsions for virtual screening.
- **ChemDraw Unpacker (`.cdxml`):** Filters reaction arrows, solvents, and annotations, cleanly extracting genuine chemical fragments into single files or a consolidated ZIP.
- **SuperAtom & Protecting Group Expander:** Unpacks common organic abbreviations (`Boc`, `Fmoc`, `OTs`, `TFA`, `OMe`) into explicit atomic connections before 3D embedding.
- **Crystallographic CIF to XYZ:** Translates fractional crystal symmetry coordinates into Cartesian coordinates.
- **Counter-Ion & Salt Stripper:** Neutralizes screening libraries by stripping Hydrochloride, Sodium, and TFA salts.
- **Physicochemical Properties:** Real-time calculation of Lipinski's Rule of 5 parameters, Molecular Weight, LogP, TPSA, Rotatable Bonds, and H-Bond Donors/Acceptors.
- **Curated Software Downloads Hub:** Direct, verified links to official downloads for AutoDock Vina, PyMOL, Avogadro 2, UCSF ChimeraX, ORCA 6, GROMACS, VMD, OpenBabel, and NWChem.

---

## 🛡️ Pharma Zero-Disk Invariant (Security & IP Protection)

In chemical and pharmaceutical research, molecular intellectual property is critical. **iLoveChem** is architected with a strict zero-disk policy:
- **100% In-Memory RAM Execution:** All conversions and coordinate embeddings are executed in volatile RAM buffers (`io.BytesIO`, `io.StringIO`).
- **Zero Disk Persistence:** No user molecules, structures, or logs are written to physical disks.
- **Immediate In-Memory Purge:** Buffers are explicitly closed and garbage-collected immediately upon HTTP response delivery.
- **Zero AI Training:** No user data is ever stored, analyzed, or used for model training.

---

## 🚀 Local Development Setup

### 1. Prerequisites
- Python 3.11+
- [OpenBabel](http://openbabel.org/) system libraries
- `pip`

### 2. Clone & Configure
```bash
git clone https://github.com/shivamrathore12336-dev/ilovechem.git
cd ilovechem

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Configure environment variables
cp .env.example .env
```

### 3. Run Migrations & Start Server
```bash
python manage.py migrate
python manage.py runserver 127.0.0.1:8000
```
Open [http://127.0.0.1:8000/](http://127.0.0.1:8000/) in your browser.

---

## 🐳 Docker Deployment

You can build and run iLoveChem in an isolated Docker container with all chemistry C-libraries pre-compiled:

```bash
docker build -t ilovechem .
docker run -p 10000:10000 -e PORT=10000 ilovechem
```

---

## 👨‍🔬 Author & Maintainer

**Shivam Rathore**  
*Medicinal Chemist & CADD Developer*  
- **LinkedIn:** [shivam-rathore-medchem](https://www.linkedin.com/in/shivam-rathore-medchem/)  
- **GitHub:** [@shivamrathore12336-dev](https://github.com/shivamrathore12336-dev)  
- **Email:** [shivamrathore.medchem@gmail.com](mailto:shivamrathore.medchem@gmail.com)

---

## 📄 License

This project is licensed under the [MIT License](LICENSE).
