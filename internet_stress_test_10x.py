import urllib.request
import json
import time
import sys
import requests

if hasattr(sys.stdout, 'reconfigure'):
    sys.stdout.reconfigure(encoding='utf-8')

BASE_URL = 'http://127.0.0.1:8000/api/convert/'

print("=" * 85)
print("INTERNET-GROUNDED 10-DATASET & 10-ITERATION STRESS TEST HARNESS")
print("=" * 85)

headers = {'User-Agent': 'iLoveChem-Scientific-Validator/2.4'}

def fetch_url(url, timeout=12):
    req = urllib.request.Request(url, headers=headers)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        return resp.read().decode('utf-8', errors='ignore')

# -----------------------------------------------------------------------------
# 1. Fetch Real Data from Internet Repositories
# -----------------------------------------------------------------------------
print("\n[PHASE 1] Fetching live structures from PubChem & RCSB PDB...")
datasets = {}

# 1. PubChem Sildenafil (Viagra) SDF
try:
    sildenafil_sdf = fetch_url("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/sildenafil/SDF")
    datasets["PubChem_Sildenafil_SDF"] = (sildenafil_sdf, "sildenafil.sdf", "orca")
    print("  [OK] Fetched Sildenafil (Sulfonamide drug) from PubChem")
except Exception as e:
    print("  [SKIP] Sildenafil fetch failed:", e)

# 2. PubChem Oxaliplatin (Pt Coordination Complex) SDF
try:
    oxaliplatin_sdf = fetch_url("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/oxaliplatin/SDF")
    datasets["PubChem_Oxaliplatin_Pt_SDF"] = (oxaliplatin_sdf, "oxaliplatin.sdf", "orca")
    print("  [OK] Fetched Oxaliplatin (Pt Coordination complex) from PubChem")
except Exception as e:
    print("  [SKIP] Oxaliplatin fetch failed:", e)

# 3. PubChem Taxol (11 Stereocenters Diterpenoid Natural Product)
try:
    taxol_sdf = fetch_url("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/paclitaxel/SDF")
    datasets["PubChem_Paclitaxel_Taxol"] = (taxol_sdf, "taxol.sdf", "xyz")
    print("  [OK] Fetched Paclitaxel / Taxol (110+ atoms, 11 stereocenters) from PubChem")
except Exception as e:
    print("  [SKIP] Paclitaxel fetch failed:", e)

# 4. PubChem Grubbs 2nd Gen Catalyst (Ruthenium Carbene)
try:
    grubbs_sdf = fetch_url("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/grubbs%20catalyst%202nd%20generation/SDF")
    datasets["PubChem_Grubbs_Ruthenium_SDF"] = (grubbs_sdf, "grubbs.sdf", "orca")
    print("  [OK] Fetched Grubbs 2nd Gen Catalyst (Ru metal carbene) from PubChem")
except Exception as e:
    print("  [SKIP] Grubbs catalyst fetch failed:", e)

# 5. PubChem Dess-Martin Periodinane (Hypervalent Iodine(V))
try:
    dess_martin_sdf = fetch_url("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/dess-martin%20periodinane/SDF")
    datasets["PubChem_DessMartin_Hypervalent_Iodine"] = (dess_martin_sdf, "dess_martin.sdf", "orca")
    print("  [OK] Fetched Dess-Martin Periodinane (Hypervalent Iodine) from PubChem")
except Exception as e:
    print("  [SKIP] Dess-Martin fetch failed:", e)

# 6. PubChem InChI String: Morphine (Pentacyclic Chiral Alkaloid)
morphine_inchi = "InChI=1S/C17H19NO3/c1-18-7-6-17-10-3-5-13(20)16(17)21-15-12(19)4-2-9(14(15)17)8-11(10)18/h2-5,10-11,13,16,19-20H,6-8H2,1H3/t10-,11+,13-,16-,17-/m0/s1"
datasets["PubChem_InChI_Morphine"] = (morphine_inchi, "morphine.inchi", "orca")
print("  [OK] Loaded Morphine InChI string from NIH PubChem")

# 7. RCSB PDB: 1CRN (Crambin Hydrophobic Protein)
try:
    crambin_pdb = fetch_url("https://files.rcsb.org/download/1CRN.pdb")
    datasets["RCSB_PDB_1CRN_Crambin"] = (crambin_pdb, "1CRN.pdb", "xyz")
    print("  [OK] Fetched 1CRN (Crambin 327 atoms) from RCSB PDB")
except Exception as e:
    print("  [SKIP] 1CRN fetch failed:", e)

# 8. RCSB PDB: 1STP (Streptavidin + Bound Biotin Ligand with Crystal Waters)
try:
    stp_pdb = fetch_url("https://files.rcsb.org/download/1STP.pdb")
    # Test ligand extraction mode
    datasets["RCSB_PDB_1STP_Ligand"] = (stp_pdb, "1STP.pdb", "orca")
    print("  [OK] Fetched 1STP (Streptavidin-Biotin Complex) from RCSB PDB")
except Exception as e:
    print("  [SKIP] 1STP fetch failed:", e)

# 9. Real Multi-Model NMR PDB: 1G03 (20 NMR solution models)
try:
    nmr_pdb = fetch_url("https://files.rcsb.org/download/1G03.pdb")
    datasets["RCSB_PDB_1G03_NMR_20Models"] = (nmr_pdb, "1G03.pdb", "gjf")
    print("  [OK] Fetched 1G03 (20-Model NMR Solution Ensemble) from RCSB PDB")
except Exception as e:
    print("  [SKIP] 1G03 NMR fetch failed:", e)

# 10. Real Crystallography CIF: Diamond / Quartz Unit Cell
diamond_cif = """data_diamond
_cell_length_a 3.567
_cell_length_b 3.567
_cell_length_c 3.567
_cell_angle_alpha 90.00
_cell_angle_beta 90.00
_cell_angle_gamma 90.00
_symmetry_space_group_name_H-M 'F d -3 m'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
C1 C 0.0000 0.0000 0.0000
C2 C 0.2500 0.2500 0.2500
"""
datasets["Crystallography_Diamond_CIF"] = (diamond_cif, "diamond.cif", "poscar")
print("  [OK] Loaded Diamond CIF unit cell structure")

print(f"\nTotal real datasets prepared: {len(datasets)}")

# -----------------------------------------------------------------------------
# 2. RUN 10 ITERATIONS ACROSS ALL DATASETS (10 x N conversions)
# -----------------------------------------------------------------------------
print("\n" + "=" * 85)
print("[PHASE 2] Executing 10 Iterations Across All Real-World Datasets...")
print("=" * 85)

NUM_RUNS = 10
dataset_stats = {name: {"passes": 0, "failures": 0, "latencies": []} for name in datasets}

for run_idx in range(1, NUM_RUNS + 1):
    print(f"\n--- [ITERATION {run_idx}/{NUM_RUNS}] ---")
    for name, (content, fname, target_fmt) in datasets.items():
        t0 = time.time()
        try:
            payload = {
                'raw_content': content,
                'target_format': target_fmt,
                'pdb_filter_mode': 'ligand_only' if 'Ligand' in name else 'all',
                'add_hydrogens': 'true',
                'strip_salts': 'true'
            }
            r = requests.post(BASE_URL, data=payload, timeout=20)
            elapsed = time.time() - t0
            res = r.json()
            if r.status_code == 200 and res.get('success') and len(res.get('output_content', '')) > 0:
                dataset_stats[name]["passes"] += 1
                dataset_stats[name]["latencies"].append(elapsed)
                print(f"  ✓ Run {run_idx:>2}: {name:<35} ➔ {target_fmt.upper():<6} ({elapsed*1000:.1f}ms, {res.get('metadata',{}).get('num_atoms')} atoms)")
            else:
                dataset_stats[name]["failures"] += 1
                print(f"  ✗ Run {run_idx:>2}: {name:<35} FAILED: {res.get('error')}")
        except Exception as e:
            dataset_stats[name]["failures"] += 1
            print(f"  ✗ Run {run_idx:>2}: {name:<35} EXCEPTION: {e}")

# -----------------------------------------------------------------------------
# 3. STATISTICAL REPORT & SUMMARY
# -----------------------------------------------------------------------------
print("\n" + "=" * 85)
print("FINAL REPRODUCIBILITY & STRESS TEST SCORECARD (10 RUNS PER DATASET)")
print("=" * 85)
print(f"{'Dataset Name':<38} | {'Success Rate':<14} | {'Avg Latency':<12} | {'Status'}")
print("-" * 85)

total_success = 0
total_attempts = 0

for name, stats in dataset_stats.items():
    passes = stats["passes"]
    total = passes + stats["failures"]
    total_success += passes
    total_attempts += total
    pct = (passes / total * 100) if total > 0 else 0
    avg_lat = (sum(stats["latencies"]) / len(stats["latencies"]) * 1000) if stats["latencies"] else 0
    status_str = "PASS (100%)" if pct == 100 else f"WARN ({pct:.1f}%)"
    print(f"{name:<38} | {passes:>2}/{total:<2} ({pct:>5.1f}%)   | {avg_lat:>8.1f} ms  | {status_str}")

print("=" * 85)
print(f"GRAND TOTAL: {total_success} / {total_attempts} RUNS SUCCESSFUL ({(total_success/total_attempts)*100:.1f}%)")
print("=" * 85)
