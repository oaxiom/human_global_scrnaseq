import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
from glbase3 import glload
import sc_utils

plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=300)
matplotlib.rcParams['font.size'] = 6

sc.settings.figdir = 'organ'

adata = sc.read('learned-noscvelo.h5ad')

spot_size = 3

for cell_type in set(adata.obs['organ']):
    print(cell_type)

    #cells = sc.pl.umap(adata, color='batch', groups=['batch1']

    sc.pl.tsne(adata, color='organ', groups=cell_type, use_raw=False, size=spot_size,
        legend_loc='on data', vmax=3, show=False, save=f'markers-{cell_type}.pdf')

    sc.pl.umap(adata, color='organ', groups=cell_type, use_raw=False, size=spot_size,
        legend_loc='on data', show=False, save=f'markers-{cell_type}.pdf')

