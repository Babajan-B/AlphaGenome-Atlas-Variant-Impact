import json
from pathlib import Path
import pandas as pd

root = Path(__file__).resolve().parents[1]
out = root / 'results' / 'objective2_benchmark'
ext = out.parent / 'external_predictors'
df = pd.read_parquet(out / 'variant_scores.parquet')
parts = df.canonical_id.str.split('_', expand=True)
df['chrom'] = parts[1]
df['pos'] = pd.to_numeric(parts[2])
df['ref'] = parts[3]
df['alt'] = parts[4]
key = ['chrom', 'pos', 'ref', 'alt']
summary = {}
for name in ('alphamissense', 'cadd', 'revel'):
    path = next(ext.glob(name + '_*_query_matches.csv'))
    x = pd.read_csv(path)
    x['pos'] = pd.to_numeric(x['pos'])
    x = x.drop_duplicates(key)
    value_columns = [c for c in x.columns if c not in key + ['query_key']]
    df = df.merge(x[key + value_columns], on=key, how='left', validate='one_to_one')
    summary[name] = {c: int(df[c].notna().sum()) for c in value_columns}
df.to_parquet(out / 'variant_scores_with_external_predictors.parquet', index=False)
(out / 'external_join_summary.json').write_text(json.dumps({
    'variants': len(df), 'exact_match_nonmissing': summary,
    'join': 'GRCh38 exact chromosome, position, reference and alternate alleles'
}, indent=2) + '\n')
print(json.dumps(summary))
