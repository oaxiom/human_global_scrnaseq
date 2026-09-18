import logging, matplotlib, os, sys
import scanpy as sc
import scvelo as scv
import scanpy.external as sce
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)
matplotlib.rcParams['font.size'] = 6
sc.settings.autoshow = False
scv.settings.autoshow = False

if __name__ == '__main__': # This solves a crazy weird bug in the joblib!

    adata = sc.read('normed.scranpy.h5ad')
    print(adata)

    sc.pp.log1p(adata)

    adata.uns['log1p']['base'] = None

    print('Number of cells: {:d}'.format(adata.n_obs))

    sc.pp.highly_variable_genes(adata, min_mean=0.1, max_mean=10, min_disp=1)
    sc.pl.highly_variable_genes(adata, save='-highly_variable.pdf')

    # Calculate the visualizations
    sc.pp.scale(adata, max_value=10, zero_center=False)

    sc.pp.pca(adata, n_comps=100, use_highly_variable=True, svd_solver='arpack')
    sc.pl.pca_variance_ratio(adata, log=True, save='-pca_variance.pdf')

    sce.pp.harmony_integrate(adata, 'replicate',
        adjusted_basis='X_pca',
        # decreasing theta and increasing lambda, to achieve a less integrated embedding
        theta=1, # default value is 1.0 ?! Although the docs say it is 2.0 in the R version.
        lamb=1, # default is 1.0
        max_iter_harmony=30,
        )

    sc.pp.neighbors(adata)

    '''
    import bbknn
    bbknn.bbknn(adata)
    #sc.pp.neighbors(adata)
    '''

    sc.tl.tsne(adata)
    sc.tl.umap(adata, min_dist=0.3)

    # Perform clustering - using highly variable genes
    res = [2.0, 0.8, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05]
    for r in res:
        sc.tl.leiden(adata, resolution=r, key_added='leiden_r{0:.2f}'.format(r),
                     flavor='igraph', n_iterations=2)

    #adata.write('./scvelo/learned-scvelo.h5ad')

    #del adata.layers["spliced"] # Remove these to save space
    #del adata.layers["unspliced"]

    adata.write('./learned-noscvelo.h5ad')

    todraw = ['replicate', 'n_counts', 'sample_type', 'organ', 'tissue', 'disease_state', 'disease_label'] + ['leiden_r{0:.2f}'.format(r) for r in res]

    spot_size = 3

    #Visualize the clustering and how this is reflected by different technical covariates
    sc.pl.tsne(adata, color=todraw, size=spot_size, legend_loc='on data', save='-tsne.pdf')
    sc.pl.umap(adata, color=todraw, size=spot_size, legend_loc='on data', save='-umap.pdf')
