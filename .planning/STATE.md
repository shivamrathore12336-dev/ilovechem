# iLoveChem — Live Project State

**Current Milestone:** Milestone 2.0: Inbuilt Molecular Docking Studio (AutoDock Vina Engine)  
**Active Phase:** Phase 6: Docking Kernel & PDBQT Preprocessor (Phase 1 of Docking)  
**Last Updated:** 2026-10-25  
**Author:** Shivam Rathore  

---

## 📈 System Health & Metrics
- **Live Production URL:** https://ilovechem.onrender.com/ (HTTP 200 OK)
- **Render Production Service:** `srv-dbe71cavcj2c73buidvg` (Commit `fb7a0df` LIVE)
- **High-Availability Monitor:** 100% Uptime, zero cold-starts
- **Codebase Mapping:** COMPLETE (7/7 documents in `.planning/codebase/` stamped at `fb7a0df`)
- **Onboarding Status:** Ready (`.planning/onboarding/SUMMARY.md` committed)
- **Test Invariants:**
  - 11 / 11 Django unit tests passing
  - 20 / 20 GRAVITAS scientific probes passing (100%)
  - 20 / 20 Code review & security hardening findings patched

---

## 🧭 Active Context & Decisions
- **Decision 1 (Pharma Zero-Storage Invariant):** All docking preparations (Receptor PDBQT, Ligand PDBQT, pocket sniffing, docking scoring) must run in volatile RAM buffers. No user structures or protein targets touch persistent disk.
- **Decision 2 (Rigorous Redocking Benchmark):** Binding affinity scoring must use the AutoDock Vina scoring function with redocking RMSD $< 2.0\text{ \AA}$ against standard crystallographic complexes (e.g. 1HSG, 3ERT).
- **Decision 3 (Smart Auto-Pocket Centroid):** Co-crystallized ligands in PDB files must be auto-detected to calculate the binding site centroid $(X, Y, Z)$ automatically, eliminating manual grid box guesswork for users.

---

## 🎯 Active Execution Plan: Phase 6 (Docking Kernel & PDBQT Preprocessor)
1. **`calculator/docking_core.py`:** Create standalone, high-performance macromolecule preprocessor:
   - `prepare_receptor_pdbqt(pdb_content_str)`: Strip crystallographic waters, detect alternate conformations, add polar hydrogens, assign Gasteiger charges, merge non-polar hydrogens into heavy carbons.
   - `prepare_ligand_pdbqt(mol_input)`: Identify rotatable bonds, generate AutoDock torsion tree (`ROOT`, `BRANCH`, `TORSDOF`), assign partial charges.
   - `detect_binding_pocket_centroid(pdb_content_str)`: Extract co-crystallized ligand HETATMs, compute bounding box & centroid $(X, Y, Z)$ with suggested search volume $(20\times 20\times 20\text{ \AA})$.
2. **AutoDock Vina Engine Binding:** Connect Python Vina engine or lightweight execution bridge with in-memory pipes.
3. **Verification Probes:** Test on HIV-1 Protease (`1HSG`) with Indinavir ligand, validating accurate PDBQT formatting and centroid extraction.
