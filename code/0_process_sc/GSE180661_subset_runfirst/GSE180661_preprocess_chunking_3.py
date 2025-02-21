import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from matplotlib.pyplot import rc_context

import os
import argparse

parser = argparse.ArgumentParser(description='Step 1 in processing')
parser.add_argument('--indir', dest='data_path',
                    help='Directory where the input H5AD is, and where the written one will be')


args = parser.parse_args()



data_path = args.data_path

data_path = f"{os.getcwd()}/../data/single_cell/GSE180661/"
result_ad_file = f"{data_path}/GSE180661_subset.h5ad"

# read in the expression

adata = None
for curr_samp in range(4):
    in_ad_file = f"{data_path}/GSE180661_subset2_{curr_samp}.h5ad"

    adata_0 = sc.read_h5ad(in_ad_file)
    if adata is None:
        adata = adata_0
    else:
        adata = ad.concat([adata, adata_0])



# write the chunk
adata.write(result_ad_file)
