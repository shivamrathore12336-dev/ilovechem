---
last_mapped_commit: fb7a0df1088938b104f273279d7d7919077adbcf
last_mapped_at: 2026-10-25
---
# Testing Strategy & Test Suites

**Analysis Date:** 2026-10-25  
<!-- refreshed: 2026-10-25 -->

## Test Architecture

The test framework is built on Django's native `django.test.TestCase` supplemented with standalone scientific integration suites.

### 1. Django Test Suite (`calculator/tests.py`)

- **Execution Command:** `python manage.py test`
- **Current Coverage:** 11 active test cases covering:
  - PDBQT Gasteiger partial charge calculation on carboxylic acids (Aspirin)
  - Organometallic metallocene (Ferrocene) descriptor calculation & 3D generation
  - Giant log file reverse streaming parser ($20\text{ MB+}$ logs)
  - PDB residue column alignment and atomic symbol extraction
  - Security defenses: HPC command injection, path traversal, XXE attack prevention, DoS oversized payloads, and Zero-Storage RAM verification
  - Spin multiplicity parity corrections (radicals, anions, cations)
  - Sitemap and Robots.txt canonical origin routing

### 2. GRAVITAS Scientific Integration Probe Suite (`run_20_gravitas_tests.py`)

- **Execution Command:** `python run_20_gravitas_tests.py`
- **Specimens Tested (20 Probes):**
  1. Paclitaxel (Taxol) — Complex Natural Product (11 stereocenters)
  2. Cisplatin — Inorganic / Transition Metal Complex (Pt)
  3. Meso-Tartaric Acid — Stereochemical Symmetry / Achiral
  4. TEMPO Radical — Open-shell odd-electron radical
  5. Ferrocene — Sandwich metallocene complex
  6. Gramicidin S Fragment — Constrained cyclic peptide macrocycle
  7. Boc-L-Alanine TFA Salt — Protecting group + acid salt
  8. Atorvastatin — Polysubstituted pyrrole drug
  9. Multi-Fragment Scheme — Reaction mixture parsing
  10. Aspirin — Classic NSAID docking ligand (PDBQT)
  11. Remdesivir — Phosphoramidate ProTide nucleotide antiviral
  12. C60 Fullerene Segment — Carbon nanomaterial
  13. L-Tryptophan Zwitterion — Amino acid zwitterion
  14. Cetylpyridinium Cation — Permanent cationic ammonium salt
  15. Sulfamethoxazole — Hypervalent sulfur sulfonamide
  16. Sucrose — Complex disaccharide polyol
  17. Coronene — Planar ultra-rigid 7-ring PAH
  18. Levothyroxine (T4) — Tetra-iodinated heavy halogen hormone
  19. D-(+)-Camphor — Bicyclic bridged terpenoid
  20. Rhodamine B Cation — Conjugated xanthene dye

*Automated codebase analysis: 2026-10-25*
