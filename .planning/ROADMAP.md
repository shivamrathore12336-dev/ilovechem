# iLoveChem — Product Roadmap

## Milestone 1.0: Foundation & Core Converter Studio (Completed ✅)
- **Phase 1: Chemistry Engine & Guardrails:** RDKit + OpenBabel core conversion pipeline with spin multiplicity, salt stripping, and SuperAtom expander.
- **Phase 2: iLovePDF Minimalist UI:** Header mega-dropdown, 12-tool quick selector, zero-disk security modal, benchmark preset catalog.
- **Phase 3: Production Deployment & SEO:** Dockerized container on Render, Google Search Console verification, dynamic sitemap.
- **Phase 4: Standalone Downloads Hub:** Separated 10-tool academic software directory into `/downloads/`.
- **Phase 5: GRAVITAS Stress Testing:** 20 diverse, complex specimens tested with 100% pass rate.

---

## Milestone 2.0: Advanced Scientific Visualization & Bulk Automation (Current Focus 🚀)
- **Phase 6: Interactive 3D Molecular Web Viewer:**
  - Embed lightweight `3Dmol.js` or `NGL Viewer` canvas inside the converter studio.
  - Allow users to rotate, zoom, and inspect generated 3D conformers, docking poses, and electron surfaces directly in the browser before downloading.
- **Phase 7: Multi-File Batch Queue:**
  - Enable multi-file drag-and-drop (e.g. 50 SDF/PDB files at once).
  - Background async conversion with progress bar and 1-click consolidated ZIP bundle download.
- **Phase 8: Additional Quantum Chemistry Engines:**
  - Full support for Q-Chem (`.in`), VASP (`POSCAR`/`INCAR`), CP2K, and Molpro input builders.
  - Explicit basis set compatibility matrix for heavy elements (LANL2DZ ECP).

---

## Milestone 3.0: API & Ecosystem Integrations (Planned 🔭)
- **Phase 9: Public REST API & Python SDK:**
  - Clean API documentation (`/api/docs/`) with rate limiting and API key optionality.
  - PyPI package `pip install ilovechem` for terminal conversions.
- **Phase 10: ChEMBL & PubChem Direct Fetch:**
  - Search by compound name, CID, or ChEMBL ID to instantly fetch, convert, and download without needing local files.
