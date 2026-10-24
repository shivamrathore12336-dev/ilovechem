# iLoveChem — Architectural Blueprint & Codebase Map

## System Topology

```
+-----------------------------------------------------------------------+
|                           CLIENT BROWSER                              |
|   Tailwind CSS Minimalist UI · Drag-and-Drop Dropzone · Lucide Icons   |
+-----------------------------------+-----------------------------------+
                                    |
                    HTTP Multipart / JSON REST API
                                    |
                                    v
+-----------------------------------------------------------------------+
|                    DJANGO APPLICATION SERVER (Gunicorn)                |
|                                                                       |
|  chemsite/urls.py ──> calculator/studio_views.py                      |
|                                                                       |
|  [Entry Guards]                                                       |
|  - DeepContentSniffer: Detects binary .chk, empty files, true formats |
|  - SecurityGuard: Strips malicious path traversal characters         |
|  - Stream Parser: Reverse tail scanner for 50MB+ Gaussian/ORCA logs   |
|                                                                       |
|  [Core Chemistry Engines] (In-Memory Volatile RAM Only)                |
|  - calculator/converter_core.py: Universal chemical parsing           |
|  - calculator/scientific_engine.py: Spin multiplicity & parity audit |
|  - RDKit: ETKDG v3 3D coordinate generation, MMFF94/UFF minimization  |
|  - OpenBabel 3.1.1: Multi-format I/O, PDBQT Gasteiger charges        |
|  - ASE: Fractional crystallographic symmetry unrolling for CIF/POSCAR |
|                                                                       |
|  [Output Generation]                                                  |
|  - Clean ASCII/UTF-8 buffer or Base64 ZIP Archive                     |
|  - RAM buffer purged immediately on HTTP response                     |
+-----------------------------------------------------------------------+
```

## Directory Structure
- `calculator/`: Main Django application
  - `studio_views.py`: Primary HTTP views for Studio, Downloads, Sitemap, Robots, Google Verification, and `/api/convert/`
  - `views.py`: Physicochemical properties calculation API (`/api/calculate/`)
  - `converter_core.py`: Universal multi-format exporter and parser
  - `scientific_engine.py`: Scientific guardrails, parity checks, zero-coordinate guard, salt stripper
  - `templates/calculator/`:
    - `studio.html`: Main minimalist Converter Studio (iLovePDF style)
    - `downloads.html`: Standalone Software Downloads Hub (10 curated packages)
- `chemsite/`: Django project configuration (`settings.py`, `urls.py`, `wsgi.py`)
- `Dockerfile`: Debian-slim container with pre-compiled RDKit and OpenBabel libraries
- `.planning/`: GSD autonomous planning directory (`PROJECT.md`, `REQUIREMENTS.md`, `ROADMAP.md`, `STATE.md`)
