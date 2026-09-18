
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

from sc_utils import compute_sum_factors
import scanpy as sc
#import scranPY as scranpy

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)
sc.settings.autoshow = False

adata = sc.read('./filtered.h5ad')

print(adata)

sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.leiden(adata, key_added='quick_clusters', resolution=0.5, flavor='igraph', n_iterations=2)

compute_sum_factors(adata, clusters='quick_clusters',
                            parallelize=False, # True == broken?
                            algorithm='CVXPY',
                            max_size=8000,
                            min_mean=0.1,
                            plotting=True,
                            lower_bound=0.4,
                            normalize_counts=False,
                            save_plots_dir='./scranpy/')

# Size factors are in adata.obs['size_factors']

sc.pl.scatter(adata, 'size_factors', 'n_counts', show=False, save='size_factors_vs_ncounts.pdf')
sc.pl.scatter(adata, 'size_factors', 'n_genes', show=False, save='size_factors_vs_ngenes.pdf')

# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata

adata.write('./normed.scranpy.h5ad')

"""
library('data.table')
library('scran')
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Sys.setenv('R_MAX_VSIZE'=54000000000) Renviron = 130 Gb

data_mat = data.matrix(fread('dense_array.tsv.gz', sep='\t'))
clusters <- quickCluster(data_mat)
size_factors = calculateSumFactors(data_mat, clusters=clusters, min.mean=0.1)
fwrite(data.frame(size_factors), file='size_factors.csv')

length(size_factors)
"""