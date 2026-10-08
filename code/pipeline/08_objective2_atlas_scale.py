"""Bounded, resumable Atlas annotation of reclassified VUS only."""
import argparse
import importlib.util
import json
from pathlib import Path

import duckdb


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, default=Path('/path/to/project_root/AG-VUS'))
    parser.add_argument('--per-cohort', type=int, default=1000)
    parser.add_argument('--output-name', default='objective2_reclassified_batch1')
    args = parser.parse_args()
    if not 1 <= args.per_cohort <= 34623:
        parser.error('Maximum is the complete eligible benign cohort')
    if Path(args.output_name).name != args.output_name:
        parser.error('Output name must be a single directory name')
    root = args.root
    if not (root / 'logs/OBJECTIVE1.COMPLETE').exists():
        raise RuntimeError('Objective 1 completion marker missing')
    spec = importlib.util.spec_from_file_location('pilot', Path(__file__).with_name('07_objective2_atlas_pilot.py'))
    pilot = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(pilot)
    output = root / 'annotations/alphagenome' / args.output_name
    output.mkdir(parents=True, exist_ok=True)
    eligibility = root / 'annotations/alphagenome/objective2_pilot/eligibility.parquet'
    con = duckdb.connect()
    con.execute("SET memory_limit='1GB'")
    con.execute('SET threads=1')
    cursor = con.execute('''
        SELECT * EXCLUDE (rank) FROM (
            SELECT *, row_number() OVER (
                PARTITION BY trajectory ORDER BY md5(canonical_id), canonical_id
            ) AS rank
            FROM read_parquet(?)
            WHERE eligibility_reason = 'eligible_snv'
              AND trajectory IN ('VUS_TO_P_LP', 'VUS_TO_B_LB')
        ) WHERE rank <= ? ORDER BY trajectory, canonical_id
    ''', [str(eligibility), args.per_cohort])
    columns = [item[0] for item in cursor.description]
    rows = [dict(zip(columns, row)) for row in cursor.fetchall()]
    con.close()
    if len(rows) != min(args.per_cohort, 7760) + min(args.per_cohort, 34623):
        raise RuntimeError('Unexpected manifest size')
    manifest = output / 'manifest.json'
    if manifest.exists() and json.loads(manifest.read_text()) != rows:
        raise RuntimeError('Refusing to change existing annotation manifest')
    pilot.atomic_json(manifest, rows)
    # Reject incomplete checkpoints before the shared worker considers reuse.
    for record in (output / 'variant_results').glob('*.json'):
        data = json.loads(record.read_text())
        files = data.get('response_files', {})
        if data.get('requested_scorers') != pilot.SCORERS or any(
            name not in files or not Path(files[name]).is_file() for name in pilot.SCORERS
        ):
            record.rename(record.with_suffix('.invalid.json.bak'))
    qc = json.loads((root / 'reports/objective1_final_qc.json').read_text())
    source = root / 'data/processed/longitudinal/objective1_frozen_cohorts_ref_validated.parquet'
    pilot.run_queries(rows, source, output, qc)
    pilot.atomic_json(output / 'batch_complete.json', {
        'status': 'complete', 'variants': len(rows), 'completed_at_utc': pilot.now(),
        'scope': 'Deterministic balanced discovery batch; not full-cohort validation',
    })


if __name__ == '__main__':
    main()
