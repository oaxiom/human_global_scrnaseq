"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import matplotlib
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

adata = sc.read('raw_data.h5ad')
print(adata)

sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='replicate', size=0, log=False, rotation=90, cut=0, save='-qc1.pdf')


adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# Filter out the mt genes to stop them getting into the most variable
mask = np.isin(adata.var_names, adata.var['mt'], invert=True, assume_unique=True)
adata = adata[:, mask] # also slices .var[]
# The above causes trouble for some reason in scran

sc.pl.violin(adata, ['n_genes', 'n_counts', 'pct_counts_mt'], groupby='replicate', size=0, log=False, cut=0, rotation=90, save='-qc2.pdf')
sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='-QC-prefilter-MT.pdf')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='-QC-postfilter-countsby_genes.pdf')

# Remove MT
adata = adata[adata.obs.pct_counts_mt < 20, :]

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='-QC-postfilter-MT.pdf')

#sc.pp.filter_genes(adata, min_counts=1)
sc.pp.filter_genes(adata, min_cells=50) # Only filter genes here, don't let MT genes influence.

sc.pl.violin(adata, ['n_genes', 'n_counts', 'pct_counts_mt'], groupby='replicate', size=0, log=False, cut=0, rotation=90, save='-qc3.pdf')

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./filtered.h5ad')

with open('gene_names.filtered.tsv', 'w') as oh:
    for g in adata.var_names:
        oh.write(f'{g}\n')

'''
print('Save CytoTRACE matrix...')
# CytoTRACE version, very slow!
dense = pd.DataFrame(
    adata.X.T.todense(),
    columns=adata.obs_names,
    index=adata.var_names,
    )
dense.to_csv('dense_array.T.csv.gz')
print('Done')
'''