---
last_mapped_commit: fb7a0df1088938b104f273279d7d7919077adbcf
last_mapped_at: 2026-10-25
---
# External Integrations & Boundary Interfaces

**Analysis Date:** 2026-10-25  
<!-- refreshed: 2026-10-25 -->

## Network & External Services

- **Database:** SQLite (`db.sqlite3`) for default Django framework state; **Zero Molecular Persistence** for user chemistry (all molecular processing is strictly 100% in-memory / volatile RAM).
- **Authentication Providers:** Stateless public endpoints; no OAuth or third-party auth dependencies.
- **Third-Party CDNs:**
  - Tailwind CSS: `https://cdn.tailwindcss.com`
  - Lucide Icons: `https://unpkg.com/lucide@latest`
  - 3Dmol.js WebGL: `https://3dmol.org/build/3Dmol-min.js`
  - Google Fonts: Plus Jakarta Sans / JetBrains Mono

## Inbound REST API Endpoints

- `POST /api/convert/`: Universal multi-format chemical file conversion, geometry optimization, quantum input generation, and PDBQT preparation.
- `POST /api/calculate/`: Physicochemical property calculation, Lipinski Rule of 5 evaluation, 2D PNG/SVG depiction, and 3D coordinate molblock.
- `GET /`: Studio dashboard.
- `GET /about/`: Platform capabilities, scientific citation (BibTeX), and REST API documentation.
- `GET /downloads/`: Computational chemistry software downloads hub.
- `GET /sitemap.xml` & `GET /robots.txt`: Search engine discovery metadata.

## Zero-Storage Invariant Guarantee

- No user-uploaded chemical files, coordinates, or SMILES strings are ever written to disk or recorded in database tables.
- All conversions execute in volatile memory buffers (`io.BytesIO`, `io.StringIO`) and are immediately garbage-collected upon response completion.

*Automated codebase analysis: 2026-10-25*
