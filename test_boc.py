import requests

cdxml_test = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<fragment id="1">
  <n id="10" p="200 200" Element="7"/>
  <n id="11" p="220 200" NodeType="Fragment"><t p="220 200"><s>Boc</s></t></n>
  <b id="12" B="10" E="11" Order="1"/>
</fragment>
</page>
</CDXML>"""

r = requests.post('http://127.0.0.1:8000/api/convert/', data={
    'raw_content': cdxml_test,
    'target_format': 'xyz',
    'fragment_index': 0
})

res = r.json()
print('Boc-Amine Success:', res.get('success'))
print('Boc-Amine Formula:', res.get('metadata', {}).get('formula'))
print('Boc-Amine Atoms:', res.get('metadata', {}).get('num_atoms'))
