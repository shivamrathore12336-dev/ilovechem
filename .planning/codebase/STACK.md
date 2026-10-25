---
last_mapped_commit: fb7a0df1088938b104f273279d7d7919077adbcf
last_mapped_at: 2026-10-25
---
# Technology Stack & Dependencies

**Analysis Date:** 2026-10-25  
<!-- refreshed: 2026-10-25 -->

## Core Runtime & Languages

- **Runtime Environment:** Python 3.11 (`python:3.11-slim` Debian base)
- **Web Server Gateway:** Gunicorn 22.0.0 (`chemsite.wsgi:application`)
- **Web Framework:** Django 5.1 (`django>=5.1`)
- **API Framework:** Django REST Framework 3.15.0 (`djangorestframework>=3.15.0`)
- **Frontend Layer:** HTML5, Tailwind CSS (via CDN), Vanilla ES6 JavaScript, Lucide Icons, 3Dmol.js (WebGL molecular visualizer)

## Scientific & Cheminformatics Kernels

- **RDKit (`rdkit>=2023.9.1`):** Primary organic chemistry engine, SMILES kekulization, ETKDGv3 3D conformer embedding, MMFF94 force-field optimization, Lipinski Rule of 5 and MedChem physicochemical descriptor generation.
- **OpenBabel (`openbabel-wheel>=3.1.1.23` + system `libopenbabel-dev`):** Universal chemical file parser (50+ formats), PDBQT Gasteiger partial charge assignment, UFF force field for organometallics, CDXML reaction parsing.
- **ASE (Atomic Simulation Environment, `ase>=3.29.0`):** Solid-state crystallographic operations, fractional-to-Cartesian matrix conversion, 15 Å vacuum box generation for VASP/Quantum ESPRESSO, Periodic Boundary Unwrapping.
- **cclib (`cclib>=1.8.1`):** Computational quantum chemistry log file parser for Gaussian, ORCA, Q-Chem, and NWChem outputs.
- **NumPy & SciPy (`numpy`, `scipy`):** Vectorized coordinate operations, distance matrices, nuclear clash detection, and Kabsch alignment.

## Hosting & Infrastructure

- **Cloud Platform:** Render (`srv-dbe71cavcj2c73buidvg`)
- **Deployment Mode:** Dockerized Container (`Dockerfile`)
- **Process Model:** Single worker with Gunicorn (`WEB_CONCURRENCY=1` tuned for 512 MB memory constraints on Render Free Tier)
- **Security Protocols:** TLS termination, HSTS enabled (`SECURE_HSTS_SECONDS=31536000`), CSRF cookie secure, Session cookie secure, non-root container user (UID 1000).

*Automated codebase analysis: 2026-10-25*
