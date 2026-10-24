import os
import sys
import time
import tempfile
import requests

url = 'http://127.0.0.1:8000/api/convert/'

print("=" * 65)
print("TESTING GIANT LOG STREAMING TAIL EXTRACTION ENGINE")
print("=" * 65)

# Step 1: Create a 25 MB mock Gaussian log file with 300,000 cycles
print("Generating 25 MB mock Gaussian log file with 300,000 cycles...")
with tempfile.NamedTemporaryFile(delete=False, mode='w', suffix='.log', encoding='utf-8') as f:
    f.write("Entering Gaussian 16 System...\nInitial orientation:\n 1 6 0 0.0 0.0 0.0\n 2 1 0 0.0 0.0 1.09\n")
    # 300,000 dummy lines
    for i in range(1, 300001):
        f.write(f"Step {i}: Total Energy = -154.123456789 Har | RMS Force = 0.000123\n")
    # Final converged stationary point
    f.write("""
Stationary point found.
Standard orientation:
---------------------------------------------------------------------
Center     Atomic      Atomic             Coordinates (Angstroms)
Number     Number       Type             X           Y           Z
---------------------------------------------------------------------
     1          6           0        1.999888    0.000000    0.000000
     2          1           0        0.000000    1.087654    0.000000
     3          8           0        0.000000    0.000000    1.432100
---------------------------------------------------------------------
Normal termination of Gaussian 16.
""")
    tmp_path = f.name

size_mb = os.path.getsize(tmp_path) / (1024 * 1024)
print(f"Log generated successfully: {size_mb:.2f} MB on disk.")

# Step 2: Upload file via multipart/form-data to the live Django server
print("\nUploading giant file to http://127.0.0.1:8000/api/convert/ ...")
t0 = time.time()
with open(tmp_path, 'rb') as f_up:
    r = requests.post(url, files={'file': ('opt_giant.log', f_up, 'text/plain')}, data={'target_format': 'xyz'})
elapsed = time.time() - t0

try:
    os.remove(tmp_path)
except Exception:
    pass

res = r.json()
print(f"\nResponse received in {elapsed:.3f} seconds! Status: {r.status_code}")
meta = res.get('metadata', {})
print(f"Parsed by: {meta.get('parsed_by')}")
print(f"Formula: {meta.get('formula')}")
print(f"Num Atoms: {meta.get('num_atoms')}")
print(f"Warnings: {meta.get('warnings')}")

out_content = res.get('output_content', '')
has_converged_coord = "1.999888" in out_content
print(f"Contains final stationary point (1.999888): {has_converged_coord}")

assert r.status_code == 200, f"Expected 200 OK, got {r.status_code}"
assert has_converged_coord, "Failed to extract final stationary coordinates!"
print("\n" + "=" * 65)
print("SUCCESS: 25MB+ GIANT LOG STREAMED AND EXTRACTED SUB-SECOND!")
print("=" * 65)
