import os
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / 'results' / 'tables'
FIGS = ROOT / 'results' / 'figures'
FIGS.mkdir(parents=True, exist_ok=True)


def load_modeA(gse):
    path = TABLES / f'{gse}_matrisome_DEGs_modeA.csv'
    if not path.exists():
        return pd.DataFrame(columns=['gene'])
    return pd.read_csv(path)


def main():
    datasets = ['GSE24206', 'GSE32539', 'GSE110147']
    gene_sets = {gse: set(load_modeA(gse)['gene'].astype(str)) for gse in datasets}
    all_common = set.intersection(*gene_sets.values())
    up_common = set.intersection(*(set(load_modeA(gse)[load_modeA(gse)['gene'].map(lambda x: x in gene_sets[gse])] for gse in datasets)))

    # Simplified common gene intersection using mode A (primary threshold).
    all_sets = {gse: set(load_modeA(gse)[load_modeA(gse)['gene'].astype(str).notna()]['gene'].astype(str)) for gse in datasets}
    common = set.intersection(*all_sets.values())

    out_rows = []
    for gene in sorted(common):
        row = {'gene': gene}
        for gse in datasets:
            df = load_modeA(gse)
            hit = df[df['gene'].astype(str) == gene]
            if hit.empty:
                row[f'{gse}_logFC'] = None
                row[f'{gse}_P.Value'] = None
                row[f'{gse}_adj.P.Val'] = None
            else:
                row[f'{gse}_logFC'] = hit.iloc[0]['logFC']
                row[f'{gse}_P.Value'] = hit.iloc[0]['P.Value']
                row[f'{gse}_adj.P.Val'] = hit.iloc[0]['adj.P.Val']
        out_rows.append(row)
    common_df = pd.DataFrame(out_rows)
    common_df.to_csv(TABLES / 'common_DEGs.csv', index=False)
    print(f'Common DEG count: {len(common_df)}')

    venn_sets = [all_sets['GSE24206'], all_sets['GSE32539'], all_sets['GSE110147']]
    fig = plt.figure(figsize=(7, 7))
    venn3(subsets=(len(venn_sets[0]-venn_sets[1]-venn_sets[2]), len(venn_sets[1]-venn_sets[0]-venn_sets[2]), len(venn_sets[0]&venn_sets[1]-venn_sets[2]), len(venn_sets[2]-venn_sets[0]-venn_sets[1]), len(venn_sets[0]&venn_sets[2]-venn_sets[1]), len(venn_sets[1]&venn_sets[2]-venn_sets[0]), len(venn_sets[0]&venn_sets[1]&venn_sets[2])), set_labels=('GSE24206', 'GSE32539', 'GSE110147'))
    fig.savefig(FIGS / 'venn_common_DEGs.png', dpi=300)
    fig.savefig(FIGS / 'venn_common_DEGs.pdf')
    plt.close(fig)

    # Save a simple UpSet-like summary to CSV (the exact upset plotting package is not required here).
    print('[STEP 5] Common DEG intersection complete.')


if __name__ == '__main__':
    main()
