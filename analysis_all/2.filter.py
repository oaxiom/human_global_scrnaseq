"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import matplotlib
import scanpy as sc
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

adata = sc.read('raw_data.h5ad')
print(adata)

#mito_genes = adata.var_names.str.startswith('mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
#adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
#adata.obs['n_counts'] = adata.X.sum(axis=1).A1

sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='replicate', size=0, log=False, rotation=90, cut=0, save='-qc1.pdf')

# Base filtering for QC failures:
sc.pp.filter_cells(adata, min_genes=1000)
sc.pp.filter_cells(adata, max_genes=8000)
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=50000)

#sc.pp.filter_genes(adata, min_counts=1)
sc.pp.filter_genes(adata, min_cells=50) # Only filter genes here;

# Filter out the mt genes to stop them getting into the most variable
adata.var['mt'] = adata.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
#mask = np.isin(adata.var_names, mito_genes, invert=True, assume_unique=True)
#adata = adata[:, mask] # also slices .var[]
# The above causes trouble for some reason in scran

sc.pl.violin(adata, ['n_genes', 'n_counts', 'pct_counts_mt'], groupby='replicate', size=0, log=False, cut=0, rotation=90, save='-qc2.pdf')

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='-QC-prefilter-MT.pdf')
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', save='-QC-postfilter-countsby_genes.pdf')

# Remove MT
adata = adata[adata.obs.pct_counts_mt < 10, :]

sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt', save='-QC-postfilter-MT.pdf')

sc.pl.violin(adata, ['n_genes', 'n_counts', 'pct_counts_mt'], groupby='replicate', size=0, log=False, cut=0, rotation=90, save='-qc3.pdf')

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./filtered.h5ad')

oh = open('gene_names.filtered.tsv', 'w')
for g in adata.var_names:
    oh.write(f'{g}\n')
oh.close()

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