# iLoveChem — Project Definition

## 1. What This Is
**iLoveChem** (`https://ilovechem.onrender.com/`) is a free, high-performance, browser-based computational chemistry & CADD conversion studio inspired by the minimalist design philosophy of iLovePDF. It solves the critical bottleneck of chemical file format incompatibility across drug discovery, molecular docking, quantum mechanics (DFT), and crystallography.

## 2. Core Stakeholders
- **Creator & Lead Architect:** Shivam Rathore (Medicinal Chemist & CADD Developer)
- **Target Audience:** Medicinal chemists, computational biologists, DFT researchers, crystallographers, organic synthesis labs, undergraduate/graduate students, and academic educators worldwide.

## 3. Tech Stack & Infrastructure
- **Backend Framework:** Django 6.1 (Python 3.11)
- **Chemistry Engines:** RDKit (Cheminformatics, ETKDG v3 3D conformers), OpenBabel 3.1.1 (Universal Chemical Toolbox), ASE (Atomic Simulation Environment)
- **Frontend:** Tailwind CSS (Clean iLovePDF-inspired layout), Vanilla JS, Lucide Icons, JetBrains Mono & Inter typography
- **Containerization & Deployment:** Docker (Debian slim, non-root user UID 1000) deployed on Render Cloud (`srv-dbe71cavcj2c73buidvg`) with automated GitHub CI/CD (`shivamrathore12336-dev/ilovechem`).
- **High-Availability:** UptimeRobot 5-minute automated keep-alive ping daemon to eliminate Render free tier cold-starts.
- **Search Engine Optimization (SEO):** Dynamically generated `/robots.txt` and `/sitemap.xml`, verified on Google Search Console via HTML verification token.

## 4. Fundamental Invariants
1. **Pharma Zero-Disk Invariant:** 100% in-memory RAM execution (`io.BytesIO`, `io.StringIO`). No user chemical structures, SMILES, or calculation outputs are ever written to server disk or retained. Zero IP persistence, zero logs, zero AI training.
2. **Spin-Multiplicity Parity Invariant:** Automatic counting of valence electrons and formal charges to prevent Gaussian (`.com`) and ORCA (`.inp`) crashes due to spin multiplicity mismatches.
3. **Deterministic 3D Conformers:** Bound RDKit distance geometry (`maxAttempts=20`, `useRandomCoords=True`) with fallback to OpenBabel and golden-angle dispersion to guarantee non-hanging conformer generation.
4. **Clean Minimalist UX:** Ultra-fast, responsive single-screen converter studio with separate dedicated standalone hub for software downloads (`/downloads/`).

## 5. Repository & Links
- **Production Web Service:** https://ilovechem.onrender.com/
- **Software Downloads Directory:** https://ilovechem.onrender.com/downloads/
- **GitHub Repository:** https://github.com/shivamrathore12336-dev/ilovechem
- **Author LinkedIn:** https://www.linkedin.com/in/shivam-rathore-medchem/
