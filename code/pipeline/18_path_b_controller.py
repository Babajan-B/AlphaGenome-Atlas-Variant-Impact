"""HPC-side data acquisition with bounded concurrency and success-gated PBS analyses."""
import concurrent.futures
import fcntl
import json
import subprocess
import sys
import urllib.request
from pathlib import Path

ROOT = Path('/path/to/project_root/AG-VUS')
SCRIPTS = Path(__file__).parent
OUTPUT = ROOT / 'manuscript/path_b'
DOWNLOADS = {
    'am': ('AlphaMissense_hg38.tsv.gz', 'https://zenodo.org/api/records/8208688/files/AlphaMissense_hg38.tsv.gz/content'),
    'revel': ('revel-v1.3_all_chromosomes.zip', 'https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip'),
}


def task(name):
    log = ROOT / 'logs' / f'path_b_{name}_acquisition.log'
    with log.open('a') as handle:
        if name in DOWNLOADS:
            filename, url = DOWNLOADS[name]
            destination = ROOT / 'data/path_b_sources' / filename
            destination.parent.mkdir(parents=True, exist_ok=True)
            if not destination.exists():
                temporary = destination.with_suffix(destination.suffix + '.part')
                subprocess.run(['curl', '--fail', '--location', '--retry', '3', '--connect-timeout', '20',
                                '--max-time', '21600', '--continue-at', '-', '--output', str(temporary), url],
                               stdout=handle, stderr=handle, check=True)
                if temporary.stat().st_size == 0:
                    raise RuntimeError('Empty downloaded source')
                temporary.replace(destination)
            return {'path': str(destination), 'bytes': destination.stat().st_size, 'source_url': url}
        subprocess.run([sys.executable, '-u', str(SCRIPTS / '15_path_b_indexed_extract.py'),
                        '--root', str(ROOT), '--dataset', name], stdout=handle, stderr=handle,
                       check=True, timeout=86400)
        return {'dataset': name, 'extraction': 'exact_allele_extraction_complete'}


def save(report):
    destination = OUTPUT / 'acquisition_status.json'
    tmp = destination.with_suffix('.tmp')
    tmp.write_text(json.dumps(report, indent=2) + '\n')
    tmp.replace(destination)


def submit(report, name, required):
    if name in report['analysis_jobs'] or not all(report['tasks'].get(key, {}).get('status') == 'complete' for key in required):
        return
    result = subprocess.run(['/opt/pbs/bin/qsub', str(SCRIPTS / f'path_b_{name}.pbs')],
                            text=True, capture_output=True, check=True)
    report['analysis_jobs'][name] = result.stdout.strip()
    print(f'submitted={name} job={result.stdout.strip()}', flush=True)


def main():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    lock = (ROOT / 'logs/path_b_controller.lock').open('w')
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    checkpoint = OUTPUT / 'acquisition_status.json'
    report = json.loads(checkpoint.read_text()) if checkpoint.exists() else {'tasks': {}, 'analysis_jobs': {}}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(task, name): name for name in ['am', 'genomes', 'revel', 'joint', 'cadd', 'spliceai']
                   if report['tasks'].get(name, {}).get('status') != 'complete'}
        for future in concurrent.futures.as_completed(futures):
            name = futures[future]
            try:
                report['tasks'][name] = dict(status='complete', **future.result())
            except Exception as exc:
                report['tasks'][name] = {'status': 'failed', 'error_type': type(exc).__name__}
            print(json.dumps({name: report['tasks'][name]}), flush=True)
            save(report)
            submit(report, 'af', ['genomes', 'joint'])
            submit(report, 'comparators', ['am', 'revel', 'cadd', 'spliceai'])
            save(report)
    submit(report, 'af', ['genomes', 'joint'])
    submit(report, 'comparators', ['am', 'revel', 'cadd', 'spliceai'])
    save(report)
    if any(value['status'] != 'complete' for value in report['tasks'].values()):
        raise SystemExit(2)


if __name__ == '__main__':
    main()
