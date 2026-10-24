# iLoveChem — Requirements & Feature Specifications

## Status Legend
- `[x]` Implemented & Verified in Production
- `[ ]` Backlog / Next Development Wave

---

## 1. File Format Interconversions
- [x] **ORCA Input Generator (`.inp`):** Automated valence counting, radical detection (Singlet vs Doublet), `%maxcore`, `%pal`, method & basis selection.
- [x] **Gaussian Route Card Builder (`.com` / `.gjf`):** Net charge and multiplicity resolution, `%mem`, `%nprocshared`, guaranteed 2 trailing blank lines.
- [x] **AutoDock Vina Preparation (`.pdbqt`):** 1-click Gasteiger partial charges, polar hydrogens, rotatable torsions.
- [x] **ChemDraw Unpacker (`.cdxml`):** Reaction arrow/text filtering, multiple product extraction, multi-molecule ZIP export.
- [x] **SuperAtom Expansion:** Expands Boc, Fmoc, OTs, TFA, OMe, and common protecting groups prior to 3D embedding.
- [x] **3D Conformer & Minimization:** RDKit ETKDG v3 embedding with MMFF94 and UFF force-field relaxation.
- [x] **Crystallography CIF to XYZ:** Orthogonal Cartesian coordinate extraction from periodic fractional unit cells via ASE.
- [x] **Salt & Counter-Ion Stripper:** Automatic removal of HCl, Na+, TFA salts from screening decks.
- [x] **Solid-State POSCAR:** VASP POSCAR generator with automatic 15 Å vacuum padding.
- [x] **Standard Formats:** SDF, MOL2, PDB, XYZ, SMILES, InChI, CML.
- [ ] **Interactive 3D In-Browser Molecular Viewer:** Integration of 3Dmol.js / NGL Viewer for inline rotatable 3D inspection.
- [ ] **Batch Multi-File Converter:** Multi-file drag-and-drop with bulk ZIP download.

---

## 2. Platform Infrastructure & SEO
- [x] **Docker Containerization:** Python 3.11 with pre-compiled OpenBabel & RDKit system libraries.
- [x] **Render Deployment:** Auto-deploying from GitHub `main` branch with dynamic PORT binding.
- [x] **Google Search Console Ownership:** Verified via HTML token `google889b369925c8069a.html`.
- [x] **Dynamic Sitemap & Robots:** Priority-weighted `/sitemap.xml` including `/`, `/downloads/`, `/properties/`.
- [x] **Zero Cold-Start Uptime:** 5-minute UptimeRobot ping monitor keeping Render cloud active 24/7.
- [x] **Standalone Downloads Page:** Dedicated `/downloads/` hub with 10 industry-standard packages and category filters.

---

## 3. Security & IP Protection
- [x] **100% In-Memory RAM Policy:** Zero disk writes for user chemical data.
- [x] **Environment-Variable Secret Management:** `SECRET_KEY`, `DEBUG`, `ALLOWED_HOSTS` loaded via env vars.
- [x] **Deep Content Sniffing:** Prevents corrupted or binary `.chk` uploads from crashing parsers.
