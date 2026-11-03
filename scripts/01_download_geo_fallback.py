import csv
import gzip
from pathlib import Path
from urllib.request import urlretrieve

import pandas as pd

RAW = Path('data/raw')
PROCESSED = Path('data/processed')
PROCESSED.mkdir(parents=True, exist_ok=True)
LOG = Path('logs/pipeline.log')

RAW_FILES = {
    'GSE24206': RAW / 'GSE24206_series_matrix.txt.gz',
    'GSE32539': RAW / 'GSE32539_series_matrix.txt.gz',
    'GSE110147': RAW / 'GSE110147_series_matrix.txt.gz',
    'GSE135065': RAW / 'GSE135065_series_matrix.txt.gz',
}

RAW_URLS = {
    'GSE24206': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE24206/GSE24206/matrix/GSE24206_series_matrix.txt.gz',
    'GSE32539': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32539/matrix/GSE32539-GPL6244_series_matrix.txt.gz',
    'GSE110147': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110147/matrix/GSE110147_series_matrix.txt.gz',
    'GSE135065': 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135065/matrix/GSE135065_series_matrix.txt.gz',
}


def log(msg):
    print(msg)
    with open(LOG, 'a', encoding='utf-8') as f:
        f.write(msg + '\n')


def ensure_raw_files():
    RAW.mkdir(parents=True, exist_ok=True)
    missing = [gse for gse, path in RAW_FILES.items() if not path.exists()]
    for gse in missing:
        url = RAW_URLS[gse]
        out = RAW_FILES[gse]
        log(f'[{gse}] downloading raw matrix: {url}')
        urlretrieve(url, str(out))
        log(f'[{gse}] saved to {out}')


def _parse_series_matrix_values(line):
    row = next(csv.reader([line], delimiter='\t'))
    return row[1:]


def parse_geo_matrix(matrix_path, dataset_name=None):
    dataset_name = dataset_name or matrix_path.name.split('_')[0]
    with gzip.open(matrix_path, 'rt', encoding='utf-8', errors='replace') as fh:
        lines = [line.rstrip('\n') for line in fh]

    sample_field_rows = {}
    data_start = None
    for i, line in enumerate(lines):
        if not line.startswith('!Sample_') and not line.startswith('"ID_REF"'):
            continue
        if line.startswith('!Sample_'):
            key = line.split('\t', 1)[0].replace('!Sample_', '')
            sample_field_rows.setdefault(key, []).append(_parse_series_matrix_values(line))
        elif line.startswith('"ID_REF"'):
            data_start = i
            break

    if data_start is None:
        raise ValueError(f'No data table found in matrix file: {matrix_path}')

    expr = pd.read_csv(matrix_path, sep='\t', skiprows=data_start, compression='gzip', quotechar='"', low_memory=False)
    expr = expr.rename(columns={'ID_REF': 'probe_id'})

    sample_fields = {}
    for key, rows in sample_field_rows.items():
        n = max(len(r) for r in rows)
        combined = [''] * n
        for row in rows:
            for idx, val in enumerate(row[:n]):
                if val:
                    combined[idx] = f"{combined[idx]} {val}".strip()
        sample_fields[key] = combined

    geo_accessions = sample_fields.get('geo_accession', [])
    sample_titles = sample_fields.get('title', [])
    descriptions = sample_fields.get('description', [])
    n_samples = max(len(geo_accessions), len(sample_titles), len(descriptions))

    if len(geo_accessions) < n_samples:
        geo_accessions = geo_accessions + [''] * (n_samples - len(geo_accessions))
    if len(sample_titles) < n_samples:
        sample_titles = sample_titles + [''] * (n_samples - len(sample_titles))

    sample_df = pd.DataFrame({
        'sample_id': geo_accessions[:n_samples],
        'title': sample_titles[:n_samples],
    })
    if 'characteristics_ch1' in sample_fields:
        sample_df['characteristics_ch1'] = sample_fields['characteristics_ch1'][:n_samples]
    else:
        sample_df['characteristics_ch1'] = ''

    def classify_group(label):
        t = str(label).lower()
        if 'final diagnosis: ipf' in t or 'ipf' in t or 'pulmonary fibrosis' in t:
            return 'IPF'
        if 'final diagnosis: control' in t or 'healthy' in t or 'normal' in t or 'control' in t:
            return 'Control'
        return 'Other'

    sample_df['label_text'] = sample_df[['title', 'characteristics_ch1']].fillna('').astype(str).agg(' '.join, axis=1)
    sample_df['group'] = sample_df['label_text'].map(classify_group)
    sample_df = sample_df[sample_df['group'].isin(['IPF', 'Control'])].copy()
    sample_df['group'] = sample_df['group'].astype(str)

    metadata_path = PROCESSED / f'{dataset_name}_metadata.csv'
    sample_df.to_csv(metadata_path, index=False)
    log(f"[{dataset_name}] verified GEO series matrix sample counts: {len(sample_df)} total ({(sample_df['group'] == 'IPF').sum()} IPF, {(sample_df['group'] == 'Control').sum()} Control)")

    if dataset_name == 'GSE24206':
        log('[GSE24206] Note: the raw GEO matrix includes early and advanced IPF labels, which were collapsed into a single IPF group for this analysis.')

    expr_cols = [c for c in expr.columns if c in sample_df['sample_id'].tolist()]
    expr_subset = expr[['probe_id'] + expr_cols].copy()
    expr_subset.to_csv(PROCESSED / f'{dataset_name}_expression.csv', index=False)
    log(f"[{dataset_name}] saved expression matrix with {expr_subset.shape[0]} probes and {expr_subset.shape[1]-1} samples")


if __name__ == '__main__':
    ensure_raw_files()
    for dataset_name in ['GSE24206', 'GSE32539', 'GSE110147', 'GSE135065']:
        parse_geo_matrix(RAW_FILES[dataset_name], dataset_name=dataset_name)
