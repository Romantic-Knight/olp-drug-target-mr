"""Fetch Open Targets tractability for druggable-union + eQTLGen gene universe.
Step 1: mapIds (symbol -> Ensembl).  Step 2: targets(ensemblIds) -> tractability.
Output: data/raw/opentargets/ot_tractability_union_eqtlgen.csv
"""
import urllib.request, json, time, csv, os, sys

URL = 'https://api.platform.opentargets.org/api/v4/graphql'
ROOT = r'D:\AI-Wisp-Science\Drug Target MR'

def gql(query, variables=None, timeout=180):
    body = json.dumps({'query': query, 'variables': variables or {}}).encode()
    req = urllib.request.Request(URL, data=body, headers={'Content-Type': 'application/json'})
    last = None
    for attempt in range(8):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                resp = json.loads(r.read().decode())
            if resp.get('errors'):
                last = resp['errors']
                print(f'  api error (attempt {attempt+1}): {last[:1]}', flush=True)
                time.sleep(6 * (attempt + 1))
                continue
            return resp
        except Exception as e:
            last = e
            print(f'  network error (attempt {attempt+1}): {e}', flush=True)
            time.sleep(5 * (attempt + 1))
    raise RuntimeError(f'persistent failure: {last}')

# ---------- build symbol universe ----------
def load_col(path, col):
    with open(path, encoding='utf-8') as f:
        rdr = csv.DictReader(f)
        return {row[col].strip().upper() for row in rdr if row.get(col) and row[col].strip()}

union = load_col(os.path.join(ROOT, 'data/processed/targets/druggable_genes_final.csv'), 'gene')
eqtl  = load_col(os.path.join(ROOT, 'data/processed/targets/eqtlgen_genes.csv'), 'gene')
universe = sorted(union | eqtl)
print(f'union genes: {len(union)} | eQTLGen genes: {len(eqtl)} | universe: {len(universe)}', flush=True)

# ---------- step 1: mapIds ----------
MAP_Q = """query($terms: [String!]!) {
  mapIds(queryTerms: $terms, entityNames: "Target") {
    mappings { term hits { id name } }
  }
}"""
sym2ensg = {}
for i in range(0, len(universe), 1000):
    chunk = universe[i:i + 1000]
    r = gql(MAP_Q, {'terms': chunk})
    for m in r['data']['mapIds']['mappings']:
        term = m['term'].strip().upper()
        ensg = [h['id'] for h in m.get('hits', []) if h.get('id', '').startswith('ENSG')]
        if ensg:
            sym2ensg[term] = ensg[0]
    print(f'mapIds {min(i + 1000, len(universe))}/{len(universe)} -> mapped {len(sym2ensg)}', flush=True)
    time.sleep(0.5)
print(f'mapped symbols: {len(sym2ensg)}', flush=True)

# ---------- step 2: targets tractability ----------
TGT_Q = """query($ids: [String!]!) {
  targets(ensemblIds: $ids) { approvedSymbol tractability { label modality value } }
}"""
ensg2sym = {v: k for k, v in sym2ensg.items()}
all_ensg = list(ensg2sym.keys())
print(f'unique ensembl ids: {len(all_ensg)}', flush=True)

records = {}
for i in range(0, len(all_ensg), 200):
    chunk = all_ensg[i:i + 200]
    r = gql(TGT_Q, {'ids': chunk})
    for t in r['data']['targets']:
        records[t['approvedSymbol'].strip().upper()] = {
            'ensembl_id': next((e for e in chunk if ensg2sym[e] == t['approvedSymbol'].strip().upper()), ''),
            'tractability': t['tractability'],
        }
    print(f'targets {min(i + 200, len(all_ensg))}/{len(all_ensg)} -> {len(records)} records', flush=True)
    time.sleep(0.4)

# ---------- save ----------
out = os.path.join(ROOT, 'data/raw/opentargets/ot_tractability_union_eqtlgen.csv')
os.makedirs(os.path.dirname(out), exist_ok=True)
with open(out, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['symbol', 'ensembl_id', 'tractability'])
    for sym in sorted(records):
        w.writerow([sym, records[sym]['ensembl_id'],
                    json.dumps(records[sym]['tractability'], ensure_ascii=False)])
print(f'written {out} | records: {len(records)}', flush=True)
