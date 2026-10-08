"""Server-side, fail-closed success gates for Atlas batch expansion."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

ROOT = Path('/path/to/project_root/AG-VUS')
BASE = ROOT / 'annotations/alphagenome'
SCRIPT = Path(__file__).with_name('08_objective2_atlas_scale.py')


def verify(directory, expected):
    summary = json.loads((directory / 'pilot_summary.json').read_text())
    if summary['status_counts'] != {'annotated': expected}:
        raise RuntimeError('Annotation status gate failed')
    manifest_path = directory / 'manifest.json'
    rows = json.loads(manifest_path.read_text())
    if len(rows) != expected or len({r['canonical_id'] for r in rows}) != expected:
        raise RuntimeError('Manifest size/uniqueness gate failed')
    for row in rows:
        token = hashlib.sha256(row['canonical_id'].encode()).hexdigest()[:24]
        record = json.loads((directory / 'variant_results' / (token + '.json')).read_text())
        if record['canonical_id'] != row['canonical_id'] or record['status'] != 'annotated':
            raise RuntimeError('Variant identity/status gate failed')
        if record.get('finite_avi_values', 0) < 1:
            raise RuntimeError('Finite AVI gate failed')
        for scorer in ('AVI_SCORE', 'AVI_SCORE_MODEL_FEATURES', 'AVI_SCORE_FEATURE_IMPORTANCE'):
            artifact = Path(record['response_files'][scorer])
            if not artifact.is_file() or artifact.stat().st_size == 0:
                raise RuntimeError('Artifact gate failed')
    return rows


def main():
    import fcntl
    lock = (ROOT / 'logs/objective2_auto_scale.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    first = BASE / 'objective2_reclassified_batch1'
    pid = int((ROOT / 'logs/objective2_reclassified_batch1.pid').read_text())
    while not (first / 'batch_complete.json').exists():
        try:
            os.kill(pid, 0)
        except ProcessLookupError:
            raise RuntimeError('Initial worker ended without successful completion')
        time.sleep(60)
    previous = first
    verify(previous, 2000)
    for limit, name, expected in (
        (5000, 'objective2_reclassified_batch2', 10000),
        (34623, 'objective2_reclassified_full', 42383),
    ):
        target = BASE / name
        records = target / 'variant_results'
        records.mkdir(parents=True, exist_ok=True)
        # Reuse immutable successful artifacts instead of querying them again.
        for record in (previous / 'variant_results').glob('*.json'):
            destination = records / record.name
            if not destination.exists():
                destination.write_bytes(record.read_bytes())
        print(f'starting={name} expected={expected}', flush=True)
        subprocess.run([os.sys.executable, '-u', str(SCRIPT), '--per-cohort', str(limit),
                        '--output-name', name], check=True)
        verify(target, expected)
        print(f'verified={name} variants={expected}', flush=True)
        previous = target
    (ROOT / 'logs/OBJECTIVE2_RECLASSIFIED_ATLAS.COMPLETE').write_text('42383 verified annotations\n')


if __name__ == '__main__':
    main()
