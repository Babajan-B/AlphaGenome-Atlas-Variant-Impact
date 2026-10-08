"""Bounded Atlas availability probes; score extensions only after successful probes."""
import argparse
import hashlib
import importlib
import json
import fcntl
import importlib.metadata
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from alphagenome.atlas import atlas
from alphagenome.data import genome
import os

pilot = importlib.import_module('07_objective2_atlas_pilot')


def query(rows, output, client):
    records = output / 'records'
    records.mkdir(parents=True, exist_ok=True)
    scores = []
    failures = 0
    for index, row in enumerate(rows, 1):
        ident = row['canonical_id']
        path = records / (hashlib.sha256(ident.encode()).hexdigest()[:24] + '.json')
        if path.exists():
            record = json.loads(path.read_text())
            if record.get('status') == 'annotated' and record['canonical_id'] == ident:
                scores.append(record)
                continue
        record = dict(row, queried_at_utc=pilot.now(), status='api_error',
                      alphagenome_version=importlib.metadata.version('alphagenome'))
        for attempt in range(3):
            try:
                with pilot.deadline(45):
                    response = client.query_variant(genome.Variant(
                        chromosome='chr' + row['chrom'], position=int(row['norm_pos']),
                        reference_bases=row['norm_ref'], alternate_bases=row['norm_alt']),
                        requested_scorers=['AVI_SCORE'])
                data = response['AVI_SCORE']
                values = data.X.toarray() if hasattr(data.X, 'toarray') else np.asarray(data.X)
                if values.size != 1 or not np.isfinite(values).all():
                    record['status'] = 'lookup_not_single_finite_avi'
                else:
                    record.update(status='annotated', avi=float(values.ravel()[0]))
                break
            except Exception as exc:
                record['error_type'] = type(exc).__name__
                if attempt < 2:
                    time.sleep(2 ** attempt)
        pilot.atomic_json(path, record)
        scores.append(record)
        failures = failures + 1 if record['status'] == 'api_error' else 0
        print(f'{output.name} variant={index}/{len(rows)} status={record["status"]}', flush=True)
        if failures >= 3:
            raise RuntimeError('Three consecutive API failures; stopping request budget')
        time.sleep(1)
    pd.DataFrame(scores).to_parquet(output / 'scores.parquet', index=False)
    return all(item['status'] == 'annotated' for item in scores)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    output = args.root / 'manuscript/path_b/atlas_extensions'
    output.mkdir(parents=True, exist_ok=True)
    lock = (args.root / 'logs/path_b_atlas_extensions.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    source = args.root / 'data/processed/longitudinal/objective1_frozen_cohorts_ref_validated.parquet'
    con = duckdb.connect()
    con.execute("SET memory_limit='2GB'")
    con.execute('SET threads=1')
    cursor = con.execute('''SELECT canonical_id, trajectory, chrom, norm_pos, norm_ref, norm_alt,
        molecular_consequence, geneinfo FROM read_parquet(?)
        WHERE trajectory IN ('VUS_TO_B_LB','VUS_TO_P_LP') AND reference_status='valid'
          AND chrom IN ('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y')
          AND length(norm_ref) <> length(norm_alt)
          AND (length(norm_ref)=1 OR length(norm_alt)=1)
          AND regexp_full_match(norm_ref, '[ACGT]+') AND regexp_full_match(norm_alt, '[ACGT]+')
        ORDER BY md5(canonical_id || '20260918')''', [str(source)])
    columns = [col[0] for col in cursor.description]
    indels = [dict(zip(columns, row)) for row in cursor.fetchall()]
    eligibility = args.root / 'annotations/alphagenome/objective2_pilot/eligibility.parquet'
    cursor = con.execute('''SELECT * FROM read_parquet(?)
        WHERE trajectory='VUS_TO_CONFLICTING' AND eligibility_reason='eligible_snv'
        ORDER BY md5(canonical_id || '20260918') LIMIT 20000''', [str(eligibility)])
    columns = [col[0] for col in cursor.description]
    conflicting = [dict(zip(columns, row)) for row in cursor.fetchall()]
    con.close()
    pilot.atomic_json(output / 'design.json', {
        'simple_indel_count': len(indels), 'conflicting_snv_budget': len(conflicting),
        'seed': 20260918, 'sampling': 'Deterministic seeded hash order',
        'boundary': 'Post hoc extension. Published model supports observed indels; query misses do not establish that indels as a class are unsupported. Numerical API quota is unavailable; serial queries stop after three consecutive errors. Persistent VUS excluded.'})
    key = os.environ.get('ALPHA_GENOME_API_KEY')
    if not key:
        raise RuntimeError('API key unavailable')
    client = atlas.create(key)
    results = {}
    for name, rows in [('indels', indels), ('conflicting', conflicting)]:
        destination = output / name
        destination.mkdir(exist_ok=True)
        pilot.atomic_json(destination / 'manifest.json', rows)
        first = rows[:5] if name == 'indels' else rows[:200]
        try:
            probe_complete = query(first, destination, client)
        except RuntimeError:
            results[name] = 'probe_api_errors_full_extension_not_launched'
            pilot.atomic_json(output / 'extension_results.json', results)
            continue
        if not probe_complete:
            results[name] = 'probe_incomplete_coverage_full_extension_not_launched'
            continue
        if name == 'indels' and not query(rows[:100], destination, client):
            results[name] = '100_indel_probe_incomplete_coverage_full_extension_not_launched'
            continue
        complete = query(rows, destination, client)
        results[name] = 'all_annotated' if complete else 'complete_with_lookup_misses'
    pilot.atomic_json(output / 'extension_results.json', results)


if __name__ == '__main__':
    main()
