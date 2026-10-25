# iLoveChem — Product Roadmap

## Milestone 1.0: Foundation & Core Converter Studio (Completed ✅)
- **Phase 1: Chemistry Engine & Guardrails:** RDKit + OpenBabel core conversion pipeline with spin multiplicity, salt stripping, and SuperAtom expander.
- **Phase 2: iLovePDF Minimalist UI:** Header mega-dropdown, 12-tool quick selector, zero-disk security modal, benchmark preset catalog.
- **Phase 3: Production Deployment & SEO:** Dockerized container on Render, Google Search Console verification, dynamic sitemap.
- **Phase 4: Standalone Downloads Hub:** Separated 10-tool academic software directory into `/downloads/`.
- **Phase 5: GRAVITAS Stress Testing & Hardening:** 20 diverse complex specimens tested with 100% pass rate; 20-issue code review patched and live on Render.

---

## Milestone 2.0: Inbuilt Molecular Docking Studio (ChemDock Engine) (In Progress 🚀)
- **Phase 6: Inbuilt ChemDock Engine Kernel & PDBQT Preprocessor (Completed ✅):**
  - Automated Receptor Preparation: PDB clean-up, crystallographic water stripping, polar hydrogen addition, Gasteiger partial charge assignment.
  - Flexible Ligand Preparation: Rotatable bond detection, torsion tree hierarchy, non-polar hydrogen merging, ETKDGv3 3D conformation.
  - Smart Pocket Centroid Sniffer: Auto-detection of co-crystallized ligand coordinates to set default Grid Box center $(X, Y, Z)$ and dimensions.
  - Native ChemDock Engine: 100% in-process vectorized AutoDock Vina empirical free energy scoring function + quaternion Monte Carlo Simulated Annealing with spatial pocket pruning.
  - Pure in-memory execution pipeline adhering strictly to the Pharma Zero-Storage Invariant (zero external binaries, zero disk storage).
- **Phase 7: Interactive 3D Docking Studio UI & 3Dmol.js Viewport (Completed ✅):**
  - Dedicated `/docking/` studio page with intuitive two-panel layout (Inputs & Parameters + 3D Viewport).
  - 1-Click PDB ID fetch (e.g. `1HSG`, `6LU7`, `3ERT`) directly from RCSB Protein Data Bank with auto-pocket sniffing.
  - Interactive 3D Grid Box visualization with real-time XYZ sliders dynamically updating 3Dmol.js viewport.
  - Standard drug presets (Indinavir, Aspirin, Paracetamol, 4-OHT) and custom file upload.
- **Phase 8: Pose Scoring, Contact Profiler & Export (Completed ✅):**
  - Pose ranking table with binding affinities ($\Delta G$ in $\text{kcal/mol}$), RMSD lower/upper bounds.
  - Interactive pose selector with real-time 3D pose switching.
  - Live 3D hydrogen bonding contacts (yellow dashed lines with distance measurements in Ångströms).
  - 1-Click export of Docked Ligand PDBQT and consolidated Receptor-Ligand PDB complex ready for PyMOL and ChimeraX.
- **Phase 9: Redocking Benchmark & Scientific Verification (Next):**
  - Automated redocking test suite verifying RMSD $< 2.0\text{ \AA}$ against crystallographic ground truth (1HSG, 3ERT).

---

## Milestone 3.0: Bulk Automation & Ecosystem Integrations (Planned 🔭)
- **Phase 10: Multi-File Batch Queue:**
  - Multi-file drag-and-drop with async ZIP bundle download.
- **Phase 11: Public REST API & Python SDK:**
  - Clean API documentation (`/about/#rest-api`) with pip package `ilovechem`.
- **Phase 12: ChEMBL & PubChem Direct Fetch:**
  - Search by compound name, CID, or ChEMBL ID to instantly fetch, convert, and dock.
