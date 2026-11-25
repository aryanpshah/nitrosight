"""Probe which PubChem property names and lookup routes actually work today."""

import json
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
UA = {"User-Agent": "impurity-scanner/0.1 (research tool; contact via repo)"}


def get(url: str) -> tuple[int, str]:
    req = urllib.request.Request(url, headers=UA)
    try:
        with urllib.request.urlopen(req, timeout=30) as r:
            return r.status, r.read().decode("utf-8", "replace")
    except urllib.error.HTTPError as e:
        return e.code, e.read().decode("utf-8", "replace")[:200]
    except Exception as e:  # noqa: BLE001
        return -1, str(e)


def probe(label: str, url: str) -> None:
    status, body = get(url)
    print(f"  [{status}] {label}")
    if status == 200:
        try:
            props = json.loads(body)["PropertyTable"]["Properties"][0]
            for k, v in props.items():
                print(f"          {k} = {str(v)[:70]}")
        except Exception:
            print(f"          {body[:200]}")
    else:
        print(f"          {body[:160]}")
    print()


if __name__ == "__main__":
    name = urllib.parse.quote("aniline")
    print("Property-name probes (PubChem renamed SMILES properties in 2025):\n")
    probe("new style: SMILES,ConnectivitySMILES",
          f"{BASE}/name/{name}/property/SMILES,ConnectivitySMILES,MolecularFormula/JSON")
    probe("old style: CanonicalSMILES,IsomericSMILES",
          f"{BASE}/name/{name}/property/CanonicalSMILES,IsomericSMILES/JSON")

    print("Lookup-route probes:\n")
    probe("by CAS as name (57365-08-9 / Bromhexine EP Imp C)",
          f"{BASE}/name/57365-08-9/property/SMILES,MolecularFormula,IUPACName/JSON")
    probe("by parent drug name (abacavir)",
          f"{BASE}/name/abacavir/property/SMILES,MolecularFormula/JSON")
    probe("by vendor impurity name (Abacavir Nitroso 1)",
          f"{BASE}/name/{urllib.parse.quote('Abacavir Nitroso 1')}/property/SMILES/JSON")
