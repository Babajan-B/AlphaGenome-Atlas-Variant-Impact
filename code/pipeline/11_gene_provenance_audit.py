"""Gene-held-out descriptive audit and provenance checklist for Objective 2."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, average_precision_score

ROOT = Path('/path/to/project_root/AG-VUS')
INPUT = ROOT / 'reports/objective2_benchmark/variant_scores.parquet'
OUTPUT = ROOT / 'reports/objective2_benchmark'

def score(frame):
    if len(frame) == 0 or frame.label.nunique() < 2:
        return {'n': int(len(frame)), 'status': 'single_class_or_empty'}
    return {'n': int(len(frame)), 'positive': int(frame.label.sum()),
            'auroc': float(roc_auc_score(frame.label, frame.avi)),
            'auprc': float(average_precision_score(frame.label, frame.avi))}

def main():
    frame = pd.read_parquet(INPUT)
    if len(frame) != 42383 or frame.canonical_id.nunique() != 42383:
        raise RuntimeError('Unexpected benchmark table')
    frame['label'] = frame.trajectory.eq('VUS_TO_P_LP').astype(int)
    # geneinfo is a source field; this does not claim a retrained gene-held-out model.
    frame['gene'] = frame.geneinfo.fillna('UNKNOWN').astype(str).str.split(';').str[0].str.strip()
    rows = []
    for gene, subset in frame.groupby('gene', dropna=False):
        item = score(subset)
        item['gene'] = gene
        rows.append(item)
    by_gene = pd.DataFrame(rows).sort_values('n', ascending=False)
    by_gene.to_csv(OUTPUT / 'gene_group_audit.csv', index=False)
    top = by_gene.head(20).to_dict(orient='records')
    report = {
        'status': 'complete', 'variants': int(len(frame)), 'genes': int(frame.gene.nunique()),
        'overall': score(frame), 'top_gene_groups': top,
        'provenance_audit': {
            'clinvar_label_leakage': 'excluded from AVI evaluation input; trajectory used only as retrospective outcome',
            'predictor_training_overlap': 'not established from downloaded artifacts; requires AlphaGenome and AlphaMissense documentation audit',
            'gene_disjoint_claim': 'not a retrained gene-disjoint model; this is a descriptive held-out-group audit',
            'calibration': 'not assessed because AVI is not assumed to be a probability',
            'clinical_use': 'not supported; research prioritization only',
        },
        'next_required': ['audit predictor training data and terms', 'obtain MANE/UniProt mappings',
                          'run selected AlphaFold residue confidence analysis',
                          'define a true gene-disjoint train/test model if claiming generalization'],
    }
    (OUTPUT / 'gene_provenance_audit.json').write_text(json.dumps(report, indent=2) + '\n')
    (ROOT / 'logs/OBJECTIVE2_GENE_PROVENANCE_AUDIT.COMPLETE').write_text('Gene/provenance audit exported\n')
    print(json.dumps(report['overall']))

if __name__ == '__main__':
    main()
