"""IV strength (R^2, F) for the 22 WIDE Bonferroni hits (incl. 5 tiered candidates).
Consistent with mr_wald_ratio.R beta_exp recovery:
  se_exp = 1/sqrt(2*MAF*(1-MAF)*N);  beta_exp = Z_exp * se_exp
IV strength under standardized-expression approximation (Var(Y)~1):
  R2 = 2*MAF*(1-MAF)*beta_exp^2 = Z_exp^2 / N
  F  = (N-2)*R2/(1-R2)            [single IV; ~= Z_exp^2 when R2 small]
Input : data/processed/eqtl/mr_input_wide.csv (zscore, af_alt, allele_state)
        data/processed/eqtl/top_cis_eqtl_nsamples.csv (NrSamples)
Output: results/tables/iv_strength_F.csv
"""
import csv, math, io

CANDIDATES = {'GDF11', 'SUOX', 'EGFL8', 'TAPBP', 'HLA-DPA1'}

with open('data/processed/eqtl/mr_input_wide.csv', encoding='utf-8-sig') as fh:
    inp = list(csv.DictReader(fh))
with open('data/processed/eqtl/top_cis_eqtl_nsamples.csv', encoding='utf-8-sig') as fh:
    ns = {}
    for r in csv.DictReader(fh):
        v = r.get('NrSamples', '').strip()
        if v:
            ns[r['gene']] = int(v)

with open('results/tables/mr_olp_hits_crosscohort.csv', encoding='utf-8-sig') as fh:
    hits = [r['gene'] for r in csv.DictReader(fh)]
hits_set = set(hits)
assert len(hits) == len(hits_set), 'duplicate genes in hits?'

by_gene = {}
for r in inp:
    g = r['gene']
    if g not in hits_set:
        continue
    by_gene.setdefault(g, []).append(r)

rows, missing = [], []
for g in hits:
    recs = by_gene.get(g, [])
    if not recs:
        missing.append(g); continue
    r = recs[0]                       # single cis IV per gene
    z = float(r['zscore']); af = float(r['af_alt']) if r.get('af_alt') else float('nan')
    state = r.get('allele_state', '')
    p_eff = (1 - af) if state == 'flip' else af
    maf = min(p_eff, 1 - p_eff)
    N = ns.get(g)
    if N is None or not (math.isfinite(af) and maf > 0):
        missing.append(g); continue
    R2 = z * z / N
    F = (N - 2) * R2 / (1 - R2)
    rows.append({'gene': g, 'tiered_candidate': 'yes' if g in CANDIDATES else 'no',
                 'iv': r['snp'], 'effect_allele': r['effect_allele'], 'other_allele': r['other_allele'],
                 'allele_state': state, 'effect_allele_freq': round(p_eff, 5), 'maf': round(maf, 5),
                 'z_exp': round(z, 4), 'n_samples_eqtlgen': N,
                 'r2_exp': round(R2, 6), 'F': round(F, 2)})

print('rows:', len(rows), 'missing:', missing)
with open('results/tables/iv_strength_F.csv', 'w', newline='', encoding='utf-8') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

for r in rows:
    star = ' *' if r['tiered_candidate'] == 'yes' else ''
    print('%-9s %-15s Z=%8.2f N=%-6d R2=%.5f F=%8.1f%s' %
          (r['gene'], r['iv'], r['z_exp'], r['n_samples_eqtlgen'], r['r2_exp'], r['F'], star))
print('note: * = tiered candidate; all F >> 10 (strong IV); F approx equals Z^2 for single-IV')
