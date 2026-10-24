import sys
sys.stdout.reconfigure(encoding='utf-8')
import requests

url = 'http://127.0.0.1:8000/api/convert/'

tests = [
    ('Methyl Radical', '[CH3]', 'orca', None, None, '* xyz 0 2'),
    ('Acetate Anion', 'CC(=O)[O-]', 'orca', None, None, '* xyz -1 1'),
    ('Ammonium Cation', '[NH4+]', 'orca', None, None, '* xyz 1 1'),
    ('TEMPO Radical', 'CC1(C)CCCC(C)(C)N1[O]', 'orca', None, None, '* xyz 0 2'),
    ('Mistaken User Mult on Radical', '[CH3]', 'orca', 0, 1, '* xyz 0 2')
]

print('=' * 65)
print('TESTING AUTOMATIC RADICAL, CHARGE & SPIN MULTIPLICITY ENGINE')
print('=' * 65)

for name, smi, target, user_charge, user_mult, expected_header in tests:
    payload = {'raw_content': smi, 'input_format': 'smi', 'target_format': target}
    if user_charge is not None: payload['charge'] = user_charge
    if user_mult is not None: payload['mult'] = user_mult

    r = requests.post(url, data=payload)
    res = r.json()
    meta = res.get('metadata', {})
    out = res.get('output_content', '')
    has_expected = expected_header in out

    print(f'{name}:')
    print(f"  Electronic state: {meta.get('electronic_state')}")
    print(f"  Charge: {meta.get('detected_charge')} | Suggested Mult: {meta.get('suggested_mult')} | Total e-: {meta.get('total_electrons')}")
    print(f"  Header matches '{expected_header}': {has_expected}")
    if meta.get('warnings'):
        print(f"  Warnings: {meta.get('warnings')}")
    print()
