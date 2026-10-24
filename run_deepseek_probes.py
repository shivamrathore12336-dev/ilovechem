"""API-level probes against the locally running iLoveChem Django app."""
import requests, json, base64, io, zipfile, pathlib, sys
sys.stdout.reconfigure(encoding='utf-8')
BASE = "http://127.0.0.1:8000/api/convert/"
ROOT = pathlib.Path(r"C:\Users\shiva\OneDrive\Documents\deepseek-harness\default-workspace")
OUT = {}

def rec(k, v):
    OUT[k] = v
    print(f"\n=== {k} ===\n{v}")

def post(data=None, files=None):
    r = requests.post(BASE, data=data or {}, files=files, timeout=120)
    try:
        j = r.json()
    except Exception:
        j = {"_text": r.text[:300]}
    return r.status_code, j

# 1. O2 spin default
sc, j = post({"raw_content": "O=O", "target_format": "orca"})
rec("o2_spin", {"status": sc, "mult_meta": j.get("metadata", {}).get("suggested_mult"),
                "state": j.get("metadata", {}).get("electronic_state"),
                "orca_head": j.get("output_content", "")[:220]})

# 2. metal salt stripping
sc, j = post({"raw_content": "[Zn+2].[O-]C(=O)c1ccccc1", "target_format": "orca"})
rec("zn_salt", {"status": sc, "formula": j.get("metadata", {}).get("formula"),
                "charge": j.get("metadata", {}).get("detected_charge"),
                "salts": j.get("metadata", {}).get("detected_salts"),
                "stripped": j.get("metadata", {}).get("salts_stripped"),
                "warns": j.get("metadata", {}).get("warnings")})

# 3. ferrocene XYZ
import math
ferro = "11\nferrocene\nFe 0 0 0\n" + "\n".join(
    f"C {2.0*math.cos(i*2*math.pi/5):.3f} {2.0*math.sin(i*2*math.pi/5):.3f} {1.65 if i<5 else -1.65}" for i in range(10))
sc, j = post({"raw_content": ferro, "target_format": "orca"})
rec("ferrocene", {"status": sc, "formula": j.get("metadata", {}).get("formula"),
                  "metals": j.get("metadata", {}).get("metal_symbols"),
                  "warns": (j.get("metadata", {}).get("warnings") or [])[:3],
                  "orca_head": j.get("output_content", "")[:260]})

# 4. PDBQT upload
pdbqt_path = ROOT/"docking/validation/ligand/1E9H_INR_explicitH.pdbqt"
if pdbqt_path.exists():
    pdbqt = pdbqt_path.read_bytes()
    sc, j = post({"target_format": "orca"}, files={"file": ("lig.pdbqt", pdbqt)})
    rec("pdbqt_upload", {"status": sc, "detected": j.get("detected_input_format"),
                         "formula": j.get("metadata", {}).get("formula"),
                         "atoms": j.get("metadata", {}).get("num_atoms"),
                         "warns": (j.get("metadata", {}).get("warnings") or [])[:3]})
else:
    rec("pdbqt_upload", {"status": "NOT_FOUND"})

# 5. QE log
qe = b"""     Program PWSCF v.7.2 starts
     convergence has been achieved
!    total energy              =     -30.12345678 Ry
     JOB DONE.
"""
sc, j = post({"target_format": "xyz"}, files={"file": ("qe.pwout", qe)})
rec("qe_upload", {"status": sc, "detected": j.get("detected_input_format"),
                  "error": j.get("error"), "text": j.get("_text")})

# 6. CDXML real file single
cdxml_path = pathlib.Path(r"C:\Users\shiva\OneDrive\Desktop\meridianin c like derivative.cdxml")
if cdxml_path.exists():
    cdxml = cdxml_path.read_bytes()
    sc, j = post({"target_format": "xyz", "fragment_index": 0}, files={"file": ("scheme.cdxml", cdxml)})
    meta = j.get("metadata", {})
    rec("cdxml_single", {"status": sc, "detected": j.get("detected_input_format"),
                         "formula": meta.get("formula"), "frag_count": meta.get("fragment_count"),
                         "frag_list": meta.get("fragments_list", [])[:5],
                         "warns": (meta.get("warnings") or [])[:3]})

    # 7. CDXML all_zip
    sc, j = post({"target_format": "xyz", "convert_mode": "all_zip"}, files={"file": ("scheme.cdxml", cdxml)})
    z = {}
    if j.get("is_zip"):
        raw = base64.b64decode(j["zip_base64"])
        z = {"zip_bytes": len(raw), "members": []}
        with zipfile.ZipFile(io.BytesIO(raw)) as zf:
            z["members"] = [(n, zf.getinfo(n).file_size) for n in zf.namelist()][:40]
    rec("cdxml_zip", {"status": sc, "zip": z, "meta": j.get("metadata")})

# 8. method injection
sc, j = post({"raw_content": "CCO", "target_format": "orca", "method": "B3LYP\n%maxcore 999999"})
rec("inject_method", {"status": sc, "head": j.get("output_content", "")[:260]})

# 9. bad nprocs
for val in ("0", "-3", "1000000000"):
    sc, j = post({"raw_content": "CCO", "target_format": "orca", "nprocs": val})
    rec(f"nprocs_{val}", {"status": sc, "error": j.get("error"),
                          "head": j.get("output_content", "")[:240]})

# 10. malformed JSON types
try:
    r = requests.post(BASE, json={"raw_content": "CCO", "target_format": 123}, timeout=30)
    rec("json_type", {"status": r.status_code, "ctype": r.headers.get("content-type"),
                      "text": r.text[:300]})
except Exception as e:
    rec("json_type", f"EXC {e}")

try:
    evidence_dir = ROOT/"audit/evidence"
    evidence_dir.mkdir(parents=True, exist_ok=True)
    json.dump(OUT, open(str(evidence_dir/"probe_api.json"), "w"), indent=2, default=str)
    print("\nwrote probe_api.json")
except Exception as e:
    print(f"\ncould not write to ROOT evidence: {e}")
