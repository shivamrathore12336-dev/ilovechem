---
last_mapped_commit: fb7a0df1088938b104f273279d7d7919077adbcf
last_mapped_at: 2026-10-25
---
# Technical Debt & Architecture Concerns

**Analysis Date:** 2026-10-25  
<!-- refreshed: 2026-10-25 -->

## 1. Hosting Resource Ceiling (Render Free Tier 512 MB RAM)

- **Status:** Actively Guarded.
- **Risk:** Concurrently loading multiple heavy C++ chemistry libraries (`rdkit`, `openbabel`, `scipy`) can exceed 512 MB memory under high concurrency.
- **Mitigation:**
  - `WEB_CONCURRENCY=1` enforced in `render.yaml`.
  - Streaming reverse-tail seek for large quantum logs ($50\text{ MB+}$) so that the entire log file is never buffered in RAM simultaneously.
  - Zero-storage architecture guarantees immediate garbage collection after response return.

## 2. Windows vs Linux Wheel Divergence (OpenBabel Data Paths)

- **Status:** Mitigated.
- **Risk:** Windows binary wheels for `openbabel-wheel` sometimes lack bundled force-field parameter text files (`UFF.prm`, `ring-fragments-index.txt`), whereas Linux containers have system-installed `/usr/share/openbabel`.
- **Mitigation:**
  - Runtime auto-detection of `BABEL_DATADIR` in `converter_core.py`.
  - Primary 3D conformer embedding delegated to RDKit ETKDGv3 + MMFF94, using OpenBabel UFF strictly as a secondary fallback.

## 3. Quantum Input Format Diversity

- **Status:** Monitored.
- **Risk:** Computational chemistry programs (Gaussian, ORCA, Q-Chem, VASP, CASTEP) frequently update input syntax and basis set naming conventions.
- **Mitigation:**
  - DeepContentSniffer validates coordinate formats with strict Cartesian float checks before export.

*Automated codebase analysis: 2026-10-25*
