---
last_mapped_commit: fb7a0df1088938b104f273279d7d7919077adbcf
last_mapped_at: 2026-10-25
---
# Code Conventions & Engineering Standards

**Analysis Date:** 2026-10-25  
<!-- refreshed: 2026-10-25 -->

## Python & Backend Standards

- **Python Version:** Python 3.11 with modern type hints (`typing.Optional`, `typing.Tuple`, `typing.Dict`, `typing.List`).
- **Defensive Chemical Parsing:** All external file parsing (RDKit, OpenBabel, cclib, ASE) is wrapped in hierarchical try-except fallbacks. Uncaught exceptions must never bubble up as raw HTTP 500 errors to the client; meaningful scientific warnings or structured error JSON are returned instead.
- **Pharma Zero-Storage Invariant:** Strict rule that no molecular data, coordinates, or SMILES may be stored on disk or in the database. All operations must use volatile RAM buffers (`io.BytesIO`, `io.StringIO`, `pybel.readstring`).
- **Variable & Function Naming:**
  - Functions and methods: `snake_case` (e.g. `export_to_target_format`, `validate_smiles`, `_rebuild_conformer_with_rdkit`).
  - Classes: `PascalCase` (e.g. `PhysicalIntegrityChecker`, `ElectronicStateEngine`, `QuantumLogParser`).
  - Constants: `UPPER_SNAKE_CASE` (e.g. `TRANSITION_AND_HEAVY_METALS`, `KNOWN_INPUT_EXTENSIONS`).

## Frontend & JavaScript Standards

- **Frameworkless Architecture:** Vanilla modern ES6+ JavaScript embedded directly in templates to eliminate heavy build-step tooling (npm/webpack) and enable instant deployment.
- **Security Standards (OWASP):**
  - All external anchor links with `target="_blank"` must specify `rel="noopener noreferrer"`.
  - Dynamic user-supplied text injected into DOM must be sanitized via `escapeHtml()` or `textContent` to prevent XSS.
  - AJAX state-mutating requests (`POST`) must attach the `X-CSRFToken` header.
- **Styling Standards:** Tailwind CSS utility classes used consistently across all views with dark mode accents (`slate-900`, `cyan-400`, `blue-600`).

*Automated codebase analysis: 2026-10-25*
