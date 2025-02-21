import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from matplotlib.pyplot import rc_context

import os, gc
import argparse

parser = argparse.ArgumentParser(description='Step 0 in processing')
parser.add_argument('--indir', dest='data_path',
                    help='Directory where the input H5AD is, and where the written one will be')
parser.add_argument('--index', dest='samp_group_idx', type=int,
                    help='number 0 through 4')

args = parser.parse_args()



data_path = args.data_path
samp_group_idx = args.samp_group_idx

run_annot = True
if os.path.exists(f"{data_path}/GSE180661_matrix_wAnnot.h5"):
    run_annot = False


if run_annot:
    # read in the expression
    in_ad_file = f"{data_path}/GSE180661_matrix.h5"
    adata = ad.read_h5ad(in_ad_file, backed='r+')

    # get the obs matrix from the normalized data 
    # then delete the object from memory
    in_ad_file = f"{data_path}/GSE180661.h5ad"
    adata_norm = sc.read_h5ad(in_ad_file, backed='r+')
    meta_cell_df = adata_norm.obs.copy()
    del adata_norm
    gc.collect()

    # now join the metadatas

    # merge the metadata
    merged_df = adata.obs.join(meta_cell_df, how='inner')

    # remove duplicates
    cell_barcode_keep = pd.Series(merged_df.index.to_list()).drop_duplicates(keep=False)
    merged_df = merged_df.iloc[np.where(merged_df.index.isin(cell_barcode_keep))[0]]
    adata = adata[np.where(adata.obs.index.isin(cell_barcode_keep))[0],:]

    #inner: use intersection of keys from both frames, preserve the order of the left keys.
    # so this makes merged df have the correct order of the index
    left_df, right_df = adata.obs.align(merged_df, join="inner", axis=0)

    # double check
    check_val = np.all(right_df.index.to_list() == adata.obs.index.to_list())

    if check_val:
        adata = adata.to_memory()
        adata.obs = right_df
        adata.write(f"{data_path}/GSE180661_matrix_wAnnot.h5")
    else: # don't change
        print("ERROR unable to update OBS")


