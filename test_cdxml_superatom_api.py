import requests

cdxml_with_superatom = """<?xml version="1.0" encoding="UTF-8" ?>
<CDXML>
<page id="1">
<fragment id="1">
  <n id="1" p="100 100" Element="6"/>
  <n id="2" p="120 100" NodeType="Fragment"><t p="120 100"><s>Ph</s></t></n>
  <b id="3" B="1" E="2" Order="1"/>
</fragment>
<fragment id="2">
  <n id="10" p="200 200" Element="7"/>
  <n id="11" p="220 200" NodeType="Fragment"><t p="220 200"><s>Boc</s></t></n>
  <b id="12" B="10" E="11" Order="1"/>
</fragment>
</page>
</CDXML>"""

r = requests.post('http://127.0.0.1:8000/api/convert/', data={
    'raw_content': cdxml_with_superatom,
    'target_format': 'xyz',
    'fragment_index': 0
})

res = r.json()
print("Frag 0 (C-Ph -> Toluene) Success:", res.get('success'))
print("Frag 0 Formula:", res.get('metadata', {}).get('formula'))
print("Frag 0 Atoms:", res.get('metadata', {}).get('num_atoms'))
print("Frag 0 Output:\n", res.get('output_content'))
