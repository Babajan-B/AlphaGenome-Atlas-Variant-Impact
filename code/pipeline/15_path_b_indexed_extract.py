"""Extract exact GRCh38 alleles from public indexed sources with chromosome checkpoints."""
import argparse
import hashlib
import json
import urllib.request
from pathlib import Path

import pandas as pd
import pysam

FIELDS = {
    'joint': ['AF_joint', 'AF_grpmax_joint', 'fafmax_faf95_max_joint', 'AN_joint'],
    'genomes': ['AF', 'AF_grpmax', 'fafmax_faf95_max', 'AN'],
    'cadd': ['cadd_raw', 'cadd_phred'],
    'spliceai': ['spliceai_max_ds'],
}


def source_url(dataset, chrom):
    if dataset == 'cadd':
        return 'https://kircherlab.bihealth.org/download/CADD/v1.7/GRCh38/whole_genome_SNVs.tsv.gz'
    if dataset == 'spliceai':
        return ('https://ftp.ensembl.org/pub/data_files/homo_sapiens/GRCh38/variation_plugins/'
                'spliceai_scores.raw.snv.ensembl_mane_v1.4.grch38.vcf.gz')
    return (f'https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/{dataset}/'
            f'gnomad.{dataset}.v4.1.sites.chr{chrom}.vcf.bgz')


def windows(positions, gap=1000):
    start = end = positions[0]
    for pos in positions[1:]:
        if pos - end <= gap:
            end = pos
        else:
            yield start, end
            start = end = pos
    yield start, end


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--dataset', choices=FIELDS, required=True)
    parser.add_argument('--chrom')
    parser.add_argument('--limit', type=int)
    args = parser.parse_args()
    source = args.root / 'reports/objective2_benchmark/variant_scores.parquet'
    cohort = pd.read_parquet(source)
    coords = cohort.canonical_id.str.extract(r'^GRCh38_([0-9XY]+)_(\d+)_([ACGT])_([ACGT])$')
    for col, values in zip(['chrom', 'pos', 'ref', 'alt'], [coords[0], coords[1].astype(int), coords[2], coords[3]]):
        cohort[col] = values
    output = args.root / 'annotations/path_b' / args.dataset
    output.mkdir(parents=True, exist_ok=True)
    for chrom, frame in cohort.groupby('chrom'):
        if args.chrom and chrom != args.chrom:
            continue
        frame = frame.sort_values('pos')
        if args.limit:
            frame = frame.head(args.limit)
        destination = output / (f'chr{chrom}' + (f'.smoke{args.limit}' if args.limit else '') + '.parquet')
        if destination.exists():
            old = pd.read_parquet(destination)
            if set(old.canonical_id) != set(frame.canonical_id):
                raise RuntimeError('Checkpoint manifest mismatch')
            print(f'{args.dataset} chr{chrom} reused={len(old)}', flush=True)
            continue
        url = source_url(args.dataset, chrom)
        index = output / (hashlib.sha256(url.encode()).hexdigest()[:16] + '.tbi')
        if not index.exists():
            with urllib.request.urlopen(url + '.tbi', timeout=60) as response:
                content = response.read()
            tmp = index.with_suffix('.tmp')
            tmp.write_bytes(content)
            tmp.replace(index)
        queries = {(r.pos, r.ref, r.alt): r.canonical_id for r in frame.itertuples()}
        hits = {}
        with pysam.TabixFile(url, index=str(index)) as tabix:
            contig = next((value for value in [chrom, 'chr' + chrom] if value in tabix.contigs), None)
            if contig is None:
                raise RuntimeError('Contig unavailable; cannot label variants absent')
            header = '\n'.join(tabix.header)
            if args.dataset not in ['cadd', 'spliceai']:
                for field in FIELDS[args.dataset]:
                    if f'ID={field},' not in header:
                        raise RuntimeError(f'Required field absent from header: {field}')
            positions = sorted(set(frame.pos))
            for number, (start, end) in enumerate(windows(positions), 1):
                for line in tabix.fetch(contig, start - 1, end):
                    parts = line.split('\t')
                    if args.dataset == 'cadd':
                        key = (int(parts[1]), parts[2], parts[3])
                        if key in queries:
                            hits[key] = {'cadd_raw': float(parts[4]), 'cadd_phred': float(parts[5])}
                    else:
                        info = dict(item.split('=', 1) for item in parts[7].split(';') if '=' in item)
                        if args.dataset == 'spliceai':
                            for annotation in info.get('SpliceAI', '').split(','):
                                values = annotation.split('|')
                                if len(values) < 6:
                                    continue
                                key = (int(parts[1]), parts[3], values[0])
                                if key in queries:
                                    numeric = [float(value) for value in values[2:6] if value != '.']
                                    if numeric:
                                        old = hits.get(key, {}).get('spliceai_max_ds', 0)
                                        hits[key] = {'spliceai_max_ds': max(old, max(numeric))}
                            continue
                        alts = parts[4].split(',')
                        for allele, alt in enumerate(alts):
                            key = (int(parts[1]), parts[3], alt)
                            if key not in queries:
                                continue
                            values = {}
                            for field in FIELDS[args.dataset]:
                                raw = info.get(field, '.').split(',')
                                item = raw[allele] if f'ID={field},Number=A,' in header else raw[0]
                                values[field] = None if item == '.' else float(item)
                            values['FILTER'] = parts[6]
                            hits[key] = values
                if number % 100 == 0:
                    print(f'{args.dataset} chr{chrom} windows={number} hits={len(hits)}', flush=True)
        records = []
        for row in frame.itertuples():
            key = row.pos, row.ref, row.alt
            item = dict(canonical_id=row.canonical_id, chrom=chrom, pos=row.pos, ref=row.ref, alt=row.alt,
                        source_url=url, extraction_status='observed' if key in hits else 'absent_after_successful_query')
            item.update({field: None for field in FIELDS[args.dataset]})
            item.update(hits.get(key, {}))
            records.append(item)
        result = pd.DataFrame(records)
        temporary = destination.with_suffix('.tmp.parquet')
        result.to_parquet(temporary, index=False)
        temporary.replace(destination)
        print(f'{args.dataset} chr{chrom} complete={len(result)} observed={len(hits)}', flush=True)
    if not args.chrom and not args.limit:
        all_rows = pd.concat([pd.read_parquet(file) for file in output.glob('chr*.parquet') if '.smoke' not in file.name])
        if len(all_rows) != 42383 or all_rows.canonical_id.nunique() != 42383:
            raise RuntimeError('Extraction count does not reconcile')
        all_rows.to_parquet(output / 'all_variants.parquet', index=False)
        (output / 'COMPLETE.json').write_text(json.dumps({'dataset': args.dataset, 'variants': len(all_rows)}) + '\n')


if __name__ == '__main__':
    main()
