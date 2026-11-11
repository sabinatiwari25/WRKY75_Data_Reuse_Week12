"""Download the four GEO inputs and verify the hashes of the analyzed bytes."""
import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

ROOT = Path(__file__).resolve().parents[1]
raw = ROOT / '01_Raw_Data'
raw.mkdir(exist_ok=True)
manifest = json.loads((ROOT/'02_Validation/input_audit.json').read_text(encoding='utf-8'))
for entry in manifest['files']:
    name = entry['name']
    if not name.endswith(('.xlsx','.gz')):
        continue
    accession = name.split('_',1)[0]
    series_group = accession[:-3] + 'nnn'
    url = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{series_group}/{accession}/suppl/{name}'
    target = raw / name
    if target.exists():
        content = target.read_bytes()
    else:
        with urlopen(url,timeout=90) as response:
            content = response.read()
    if hashlib.sha256(content).hexdigest() != entry['sha256']:
        raise SystemExit(f'SHA-256 mismatch: {target}; inspect the file and GEO record before replacing it.')
    if not target.exists():
        target.write_bytes(content)
    print(f'Verified: {name}')
