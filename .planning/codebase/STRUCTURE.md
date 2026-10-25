---
last_mapped_commit: fb7a0df1088938b104f273279d7d7919077adbcf
last_mapped_at: 2026-10-25
---
# Directory Layout & File Structure

**Analysis Date:** 2026-10-25  
<!-- refreshed: 2026-10-25 -->

## Repository Organization

```
chemsite/
├── calculator/                         # Primary computational chemistry Django app
│   ├── templates/calculator/           # Server-rendered HTML templates
│   │   ├── studio.html                 # Converter Studio UI (drag-and-drop, 3Dmol viewer)
│   │   ├── index.html                  # Molecular properties & descriptors dashboard
│   │   ├── about.html                  # Documentation, capabilities, REST API spec
│   │   └── downloads.html              # Software hub for quantum & modeling tools
│   ├── chem_utils.py                   # RDKit descriptor calculation, 2D/3D depiction
│   ├── converter_core.py               # Multi-format conversion pipeline & format dispatch
│   ├── scientific_engine.py            # Physics guardrails, spin parity, PBC unwrapper
│   ├── phase1_engine.py                # ChemDraw CDXML XML extractor & superatom resolver
│   ├── studio_views.py                 # Views for studio, downloads, about, /api/convert/
│   ├── views.py                        # Views for /properties/ and /api/calculate/
│   ├── urls.py                         # URL routing for calculator app
│   └── tests.py                        # Comprehensive Django unit and integration test suite
├── chemsite/                           # Django project root
│   ├── settings.py                     # Security, ALLOWED_HOSTS, static files, middleware
│   ├── urls.py                         # Global route dispatcher
│   └── wsgi.py                         # WSGI application entry point for Gunicorn
├── .planning/                          # GSD autonomous planning directory
│   ├── codebase/                       # Architectural and technical reference maps
│   ├── PROJECT.md                      # Project definition and core invariants
│   ├── REQUIREMENTS.md                 # Product specifications and milestones
│   ├── ROADMAP.md                      # Phased development roadmap
│   └── STATE.md                        # Current execution state and decisions
├── Dockerfile                          # Production container specification (Python 3.11-slim)
├── render.yaml                         # Render Cloud Infrastructure-as-Code blueprint
├── requirements.txt                    # Pinned Python package dependencies
└── manage.py                           # Django administrative CLI
```

## Key Architectural Locations

- **Entry Request Routing:** `chemsite/urls.py` $\to$ `calculator/urls.py`
- **Conversion Gateway:** `calculator/studio_views.py:api_convert_file` $\to$ `calculator/converter_core.py:convert_content`
- **Scientific Audit & Parity:** `calculator/scientific_engine.py` (Classes: `PhysicalIntegrityChecker`, `ElectronicStateEngine`, `PeriodicBoundaryUnwrapper`, `DeepContentSniffer`)
- **Web Interface:** `calculator/templates/calculator/studio.html` (3Dmol.js WebGL viewport + reactive preset catalog)

*Automated codebase analysis: 2026-10-25*
