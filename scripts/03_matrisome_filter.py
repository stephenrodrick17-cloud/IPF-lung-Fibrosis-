import os
import re
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / 'data' / 'raw'
LOGFILE = ROOT / 'logs' / 'pipeline.log'


def log(msg: str):
    print(msg)
    LOGFILE.parent.mkdir(parents=True, exist_ok=True)
    with open(LOGFILE, 'a', encoding='utf-8') as fh:
        fh.write(msg + '\n')


def download_matrisome_masterlist():
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    candidates = [
        'https://docs.google.com/spreadsheets/d/12OslOlOhQQtHMOUKAFFgAQXL2zYX-g9q/export?format=csv&usp=sharing',
        'https://docs.google.com/spreadsheets/d/12OslOlOhQQtHMOUKAFFgAQXL2zYX-g9q/gviz/tq?tqx=out:csv&sheet=Sheet1',
        'https://sites.google.com/uic.edu/matrisome/matrisome-annotations/homo-sapiens',
        'https://sites.google.com/uic.edu/matrisome/home',
    ]
    for url in candidates:
        try:
            resp = requests.get(url, timeout=60, headers={'User-Agent': 'Mozilla/5.0'})
            resp.raise_for_status()
            text = resp.text
            if 'text/csv' in resp.headers.get('content-type', '').lower() or 'Homo sapiens matrisome masterlist' in text:
                out = RAW_DIR / 'Hs_Matrisome_Masterlist.csv'
                out.write_text(text, encoding='utf-8', errors='ignore')
                log(f'[MATRISOME] downloaded master list from {url}')
                return out
            if 'matrisome' in text.lower() or 'download' in text.lower():
                html_path = RAW_DIR / 'matrisome_masterlist.html'
                html_path.write_text(text, encoding='utf-8', errors='ignore')
                log(f'[MATRISOME] saved HTML snapshot from {url}')
        except Exception as exc:
            log(f'[MATRISOME] candidate failed: {url} ({exc})')

    log('[MATRISOME] Auto-download failed. Manual download required.')
    log('Please download the human matrisome master list from: https://sites.google.com/uic.edu/matrisome/matrisome-annotations/homo-sapiens')
    log('Then save it as data/raw/Hs_Matrisome_Masterlist.csv')
    return None


def read_masterlist():
    candidates = [
        RAW_DIR / 'Hs_Matrisome_Masterlist.csv',
        RAW_DIR / 'Hs_Matrisome_Masterlist.tsv',
        RAW_DIR / 'matrisome_masterlist.csv',
        RAW_DIR / 'matrisome_masterlist.tsv',
        RAW_DIR / 'Hs_Matrisome_Masterlist.xlsx',
        RAW_DIR / 'matrisome_masterlist.xlsx',
        RAW_DIR / 'matrisome_masterlist.html',
    ]
    for path in candidates:
        if not path.exists():
            continue
        try:
            if path.suffix.lower() == '.csv':
                return pd.read_csv(path, skiprows=1)
            if path.suffix.lower() == '.tsv':
                return pd.read_csv(path, sep='\t', skiprows=1)
            if path.suffix.lower() == '.xlsx':
                return pd.read_excel(path)
            if path.suffix.lower() == '.html':
                text = path.read_text(encoding='utf-8', errors='ignore')
                genes = re.findall(r'\b[A-Z][A-Z0-9.\-_]{1,}\b', text)
                return pd.DataFrame({'GeneSymbol': genes[:5000]})
        except Exception as exc:
            log(f'[MATRISOME] failed to read {path}: {exc}')
    return None


def normalize_symbol(s):
    if pd.isna(s):
        return None
    value = str(s).strip()
    if value in {'', '-', 'NA', 'N/A', 'nan'}:
        return None
    return value


def main():
    log('[STEP 3] Loading matrisome master list')
    master = read_masterlist()
    if master is None:
        download_matrisome_masterlist()
        master = read_masterlist()
    if master is None:
        raise FileNotFoundError(
            'Matrisome list not found in data/raw. Please manually download it from '
            'https://sites.google.com/uic.edu/matrisome/downloads and save it as Hs_Matrisome_Masterlist.csv.'
        )

    gene_col = next((c for c in master.columns if 'gene' in c.lower() or 'symbol' in c.lower() or 'protein' in c.lower()), None)
    if gene_col is None:
        raise ValueError(f'Could not identify gene symbol column in matrisome table: {list(master.columns)}')

    mat = master[[gene_col]].dropna().copy()
    mat.columns = ['GeneSymbol']
    mat['GeneSymbol'] = mat['GeneSymbol'].map(normalize_symbol)
    mat = mat[mat['GeneSymbol'].notna()].drop_duplicates().reset_index(drop=True)

    if {'Matrisome Division', 'Matrisome Category'}.issubset(master.columns):
        annotation = master[[gene_col, 'Matrisome Division', 'Matrisome Category']].copy()
        annotation.columns = ['GeneSymbol', 'Matrisome Division', 'Matrisome Category']
        annotation['GeneSymbol'] = annotation['GeneSymbol'].map(normalize_symbol)
        annotation = annotation[annotation['GeneSymbol'].notna()].drop_duplicates()
        mat = mat.merge(annotation, on='GeneSymbol', how='left')
        mat['Category/Division'] = (
            mat['Matrisome Division'].fillna('Unknown').astype(str).str.strip() + ' / ' +
            mat['Matrisome Category'].fillna('Unknown').astype(str).str.strip()
        )
    elif 'Category' in master.columns:
        annotation = master[[gene_col, 'Category']].copy()
        annotation.columns = ['GeneSymbol', 'Category/Division']
        mat = mat.merge(annotation, on='GeneSymbol', how='left')
    elif 'Division' in master.columns:
        annotation = master[[gene_col, 'Division']].copy()
        annotation.columns = ['GeneSymbol', 'Category/Division']
        mat = mat.merge(annotation, on='GeneSymbol', how='left')
    elif 'Class' in master.columns:
        annotation = master[[gene_col, 'Class']].copy()
        annotation.columns = ['GeneSymbol', 'Category/Division']
        mat = mat.merge(annotation, on='GeneSymbol', how='left')
    else:
        mat['Category/Division'] = 'matrisome'

    outpath = RAW_DIR / 'Hs_Matrisome_Masterlist_clean.csv'
    mat.to_csv(outpath, index=False)
    log(f'[MATRISOME] saved cleaned master list to {outpath} ({len(mat)} genes)')

    for gse in ['GSE24206', 'GSE32539', 'GSE110147']:
        deg_path = ROOT / 'results' / 'tables' / f'{gse}_limma_full.csv'
        if not deg_path.exists():
            log(f'[MATRISOME] skipping {gse}: no DEG table found')
            continue
        deg = pd.read_csv(deg_path)
        deg['gene'] = deg['gene'].astype(str).str.strip()
        merged = deg[deg['gene'].isin(mat['GeneSymbol'])].copy()
        merged = merged.merge(mat[['GeneSymbol', 'Category/Division']].rename(columns={'GeneSymbol': 'gene'}), on='gene', how='left')
        out = ROOT / 'results' / 'tables' / f'{gse}_matrisome_DEGs.csv'
        merged.to_csv(out, index=False)
        log(f'[MATRISOME] {gse}: {len(merged)} matrisome genes detected among {len(deg)} genes')

    log('[STEP 3] Matrisome filtering complete.')


if __name__ == '__main__':
    main()
