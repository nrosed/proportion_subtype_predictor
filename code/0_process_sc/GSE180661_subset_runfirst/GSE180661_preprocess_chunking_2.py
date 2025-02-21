import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from matplotlib.pyplot import rc_context

import os
import argparse

# argparse
parser = argparse.ArgumentParser(description='Step 1 in processing')
parser.add_argument('--indir', dest='data_path',
                    help='Directory where the input H5AD is, and where the written one will be')
parser.add_argument('--index', dest='samp_group_idx', type=int,
                    help='number 0 through 4')

args = parser.parse_args()


# read in params
data_path = args.data_path
samp_group_idx = args.samp_group_idx

in_ad_file = f"{data_path}/GSE180661_subset_{samp_group_idx}.h5ad"
result_ad_file = f"{data_path}/GSE180661_subset2_{samp_group_idx}.h5ad"

# read in the expression
adata = sc.read_h5ad(in_ad_file)

# format metadata
# reformat the cell types to be more granular
cell_type_vec = adata.obs['author_cell_type'].copy().astype(str)
immune_idx = np.where(cell_type_vec.isin(["T.cell", "Monocyte", "Plasma.cell", 
                                          "B.cell", "Dendritic.cell", "Mast.cell"]) )[0]
cell_type_vec[immune_idx] = "immune"

stromal_idx = np.where(cell_type_vec.isin(["Endothelial.cell", "Fibroblast"]) )[0]
cell_type_vec[stromal_idx] = "stromal"

epi_idx = np.where(cell_type_vec.isin(["Ovarian.cancer.cell"]) )[0]
cell_type_vec[epi_idx] = "tumor"
adata.obs['celltype_granular'] = cell_type_vec

adata.obs['scpred_CellType'] = adata.obs['author_cell_type']
adata.obs['cellType'] = adata.obs['author_cell_type']

adata.obs['sample_id'] = adata.obs['author_sample_id'].astype(str)
adata.obs['sample_id'] = adata.obs['sample_id'].tolist()
adata.obs['stim'] = ["CTRL"]*len(adata.obs)

# remove unsorted and ascites
idx_remove = adata.obs.sample_id.str.contains('ASCITES|UNSORTED', regex=True)
adata = adata[~idx_remove, :]

# remove rare cell types
ct_idx = adata.obs.scpred_CellType.isin(["Mast.cell", "Other", "Dendritic.cell"])
inverted_list = [not x for x in ct_idx]
adata = adata[np.where(inverted_list)[0]]

# filter out cells with less than 500 genes and genes expressed in less than 3 cells
sc.pp.filter_cells(adata, min_genes=500)
sc.pp.filter_genes(adata, min_cells=3)



# make sample_id without the CD45 string
adata.obs.sample_id = adata.obs.sample_id.str.replace("_CD45P_", "_")
adata.obs.sample_id = adata.obs.sample_id.str.replace("_CD45N_", "_")

# only keep samples with enough cells in each cell type
keep_samps = []
for curr_samp in adata.obs.sample_id.unique():
    temp = adata.obs.iloc[np.where(adata.obs.sample_id == curr_samp)[0]]
    if np.min(temp.cellType.value_counts()) >= 50:
        keep_samps = keep_samps + [curr_samp]


# filter
idx_keep = adata.obs.sample_id.isin(keep_samps)
adata = adata[np.where(idx_keep)[0], :]


# write the chunk
adata.write(result_ad_file)
