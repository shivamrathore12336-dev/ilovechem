# iLoveChem — Live Project State

**Current Milestone:** Milestone 2.0: Inbuilt Molecular Docking Studio (ChemDock Engine)  
**Active Phase:** Phase 9: Redocking Benchmark & Scientific Verification (Phase 1 & Phase 2 Complete)  
**Last Updated:** 2026-10-25  
**Author:** Shivam Rathore  

---

## 📈 System Health & Metrics
- **Live Production URL:** https://ilovechem.onrender.com/ (HTTP 200 OK)
- **Render Production Service:** `srv-dbe71cavcj2c73buidvg`
- **High-Availability Monitor:** 100% Uptime, zero cold-starts
- **Codebase Mapping:** COMPLETE (7/7 documents in `.planning/codebase/` stamped at `fb7a0df`)
- **Onboarding Status:** Ready (`.planning/onboarding/SUMMARY.md` committed)
- **Test Invariants:**
  - 18 / 18 Django unit & integration tests passing (100%)
  - 22 / 22 GRAVITAS scientific probes passing (100%)
  - 20 / 20 Code review & security hardening findings patched

---

## 🧭 Active Context & Decisions
- **Decision 1 (Pharma Zero-Storage Invariant):** All docking preparations (Receptor PDBQT, Ligand PDBQT, pocket sniffing, docking scoring) run strictly in volatile RAM buffers. No user structures or protein targets touch persistent disk.
- **Decision 2 (Pure In-Process ChemDock Engine):** Zero external binaries (`vina.exe`, etc.). Native vectorized NumPy implementation of AutoDock Vina empirical free energy scoring (Gauss1, Gauss2, Repulsion, Hydrophobic, H-bond, Torsional entropy penalty) combined with quaternion Monte Carlo Simulated Annealing and spatial pocket pruning.
- **Decision 3 (Smart Auto-Pocket Centroid Sniffer):** Co-crystallized ligands in PDB files are auto-detected to calculate the binding site centroid $(X, Y, Z)$ automatically, eliminating manual grid box guesswork for users.
- **Decision 4 (Interactive WebGL 3D Studio):** Dedicated `/docking/` page with 3Dmol.js viewer rendering macromolecular cartoon ribbons, active-site stick residues, real-time 3D Grid Box with sliders, pose table, live dashed H-bonds with distance measurements in Å, and 1-click PDB complex download.

---

## 🎯 Completed Deliverables: Milestone 2.0
1. **`calculator/docking_core.py` (ChemDockEngine):**
   - Automated receptor preparation with crystallographic water and co-crystal stripping, polar hydrogens, and rigid PDBQT formatting.
   - Spatial pocket pruning: filters 2000-atom proteins down to active pocket boundary, delivering a 12x speedup.
   - Flexible ligand preparation: 3D conformer generation (RDKit ETKDGv3 + MMFF94), AutoDock torsion tree, and Gasteiger charges.
   - Vectorized AutoDock Vina scoring and quaternion Monte Carlo simulated annealing.
   - Real-time active-site hydrogen bond contact detection with distance measurements.
2. **`calculator/docking_views.py`:**
   - `docking_studio_view`: Dedicated template view for `/docking/`.
   - `api_fetch_rcsb_pdb`: Official RCSB PDB download with pocket centroid sniffer.
   - `api_run_docking`: Full in-memory docking pipeline with consolidated PDB complex generator.
3. **`calculator/templates/calculator/docking.html`:**
   - Two-column responsive desktop/mobile layout.
   - 1-Click RCSB fetch + benchmark presets (1HSG HIV-1 Protease, 6LU7 SARS-CoV-2 Mpro, 3ERT Estrogen Receptor).
   - Drug presets (Indinavir, Aspirin, Paracetamol, 4-OHT) + custom SMILES and file upload.
   - Interactive 3D Grid Box with center and size sliders synced live with 3Dmol.js.
   - Ranked docking pose table with binding affinities ($\Delta G$), RMSD, and click-to-pose 3D viewer.
   - Detailed H-bond contact inspector with atom pairs and Ångström distances.
   - 1-Click download of Docked PDB Complex and Ligand PDBQT.
4. **Site-Wide Navigation Updated:**
   - Dedicated "Molecular Docking" navigation links in `studio.html` (header, hero banner, mega-menu), `about.html`, and `downloads.html`.
