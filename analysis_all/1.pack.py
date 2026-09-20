
"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import logging, matplotlib, os, sys, glob
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import colors
from glbase3 import glload
import anndata

plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

import sc_utils

ensg = glload(os.path.expanduser('~/hg38/gencode/hg38_gencode_ensg_tes_v42.glb'))
ensg_to_symbol = {row['ensg']: row['name'] for row in ensg}

samples = []

dummy_run = False

for f in sorted(list(glob.glob('../tecounts/*/ss.Hs_*.gz'))):
    if 'spliced' in f: continue
    if 'barcode' in f: continue

    bf = os.path.split(f)[1].replace('.tsv.gz', '')
    t = bf.split('.')[1].split('_')
    rep = bf.split('.')[2].replace('rp', '')
    rep_name = '_'.join(t).replace('Hs_', '')

    cell_type = 'XXXXXXXXXXXXXXXXXX'
    tissue =    'XXXXXXXXXXXXXXXXXX'
    organ =     'XXXXXXXXXXXXXXXXXX'

    tissue = f.split('/')[2]

    # TODO: Add age stratification data
    # TODO: Add sc/snRNA-seq method

    if tissue == 'adult_heart':
        organ = 'Heart'
        sample_type = 'Heart post-mortem sample'
        tissue = 'Heart'
        disease_state = 'Normal'
        disease_label = "NA"

    elif tissue == 'adult_eye':
        organ = 'Eye'
        sample_type = 'RPE'
        tissue = 'Eye'
        disease_state = 'Normal'
        disease_label = "NA"

    elif tissue == 'adult_endometrium':
        organ = 'Endometrium'
        sample_type = 'Endometrium biopsy'
        tissue = 'Endometrium'
        disease_state = 'Normal'
        disease_label = "NA"

    elif tissue == 'adult_ileum':
        organ = 'Ileum'
        sample_type = 'Ileum biopsy'
        tissue = 'Ileum'
        if 'inf' in f:
            disease_state = 'Diseased'
            disease_label = "Crohn's disease"
        else:
            disease_state = 'Normal'
            disease_label = "NA"

    elif tissue == 'adult_dorsolateral_prefrontal_cortex':
        organ = 'Brain'
        sample_type = 'Dorsolateral prefrontal cortex'
        tissue = 'Dorsolateral prefrontal cortex'
        if 'Ctrl' in f:
            disease_state = 'Normal'
            disease_label = "NA"
        else:
            disease_state = 'Diseased'
            disease_label = "Parkinson's disease"

    elif tissue == 'adult_ovary':
        organ = 'Ovary'
        sample_type = 'Ovary'
        tissue = 'Ovary'
        disease_state = 'Normal'
        disease_label = "NA"

    # TODO: PBMC
    # TODO: Azheimers
    # TODO: Testes

    obs_add={
        'replicate': f'{rep_name}#{rep}',
        'sample_type': sample_type,
        'organ': organ,
        'tissue': tissue,
        'disease_state': disease_state,
        'disease_label': disease_label,
        }

    print(obs_add)
    if not dummy_run:
        sam = sc_utils.sparsify(f,
            obs_add=obs_add,
            ensg_to_symbol=ensg_to_symbol,
            drop_fusions=True,
            drop_mir=True,
            drop_tes=False,
            drop_ribosomes=True,
            csv=False,
            )

        # Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial
        # uninteresting things and save space.
        sc.pp.filter_cells(sam, min_genes=1000)
        sc.pp.filter_cells(sam, max_genes=8000)
        sc.pp.filter_cells(sam, min_counts=5000)
        sc.pp.filter_cells(sam, max_counts=50000)
        # Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;

        # If the ann data has no cells, don't add it:

        if len(sam.obs) >= 0 and len(sam.var) >= 0:
            samples.append(sam)

if dummy_run:
    1/0

print('Loaded Samples...')

print('Concatenating')
adata = anndata.concat(samples)
adata.obs_names_make_unique()
del samples

adata.X = adata.X.astype('float32', copy=False)

print(adata)

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./raw_data.h5ad')

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()
