"""
Lightweight PubChem PUG-REST Name / CID / InChIKey -> SMILES Resolver.
Uses stdlib urllib with strict timeouts, input sanitization, and bounded LRU cache.
"""
import functools
import json
import re
import urllib.parse
import urllib.request
from typing import Dict, Optional

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


@functools.lru_cache(maxsize=64)
def resolve_compound_identifier(query: str) -> Optional[Dict[str, str]]:
    if not query or not isinstance(query, str):
        return None
    cleaned = query.strip()
    if not cleaned or len(cleaned) > 120:
        return None

    # Detect numeric PubChem CID vs InChIKey vs common drug name
    if re.fullmatch(r"(?:CID[:\s]*)?(\d{1,10})", cleaned, flags=re.IGNORECASE):
        cid = re.sub(r"\D", "", cleaned)
        endpoint = f"{PUBCHEM_BASE}/cid/{cid}/property/CanonicalSMILES,IsomericSMILES,IUPACName,MolecularFormula/JSON"
    elif re.fullmatch(r"[A-Z]{14}-[A-Z]{10}-[A-Z]", cleaned):
        endpoint = f"{PUBCHEM_BASE}/inchikey/{urllib.parse.quote(cleaned)}/property/CanonicalSMILES,IsomericSMILES,IUPACName,MolecularFormula/JSON"
    else:
        safe_name = urllib.parse.quote(cleaned, safe="")
        endpoint = f"{PUBCHEM_BASE}/name/{safe_name}/property/CanonicalSMILES,IsomericSMILES,IUPACName,MolecularFormula/JSON"

    req = urllib.request.Request(endpoint, headers={"User-Agent": "iLoveChem/3.0 (Computational Chemistry Workbench)"})
    try:
        with urllib.request.urlopen(req, timeout=4.5) as resp:
            payload = json.loads(resp.read(262144).decode("utf-8"))
        props = payload["PropertyTable"]["Properties"][0]
        smiles = props.get("IsomericSMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES") or ""
        if not smiles:
            return None
        return {
            "query": cleaned,
            "cid": str(props.get("CID", "")),
            "smiles": smiles,
            "iupac_name": props.get("IUPACName", ""),
            "formula": props.get("MolecularFormula", ""),
            "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{props.get('CID', '')}",
        }
    except Exception:
        return None
