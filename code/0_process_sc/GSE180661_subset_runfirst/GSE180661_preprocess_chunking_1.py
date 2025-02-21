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
parser.add_argument('--index', dest='samp_group_idx', type=int,
                    help='number 0 through 4')

args = parser.parse_args()



data_path = args.data_path
samp_group_idx = args.samp_group_idx

result_ad_file = f"{data_path}/GSE180661_subset_{samp_group_idx}.h5ad"

all_sample_ids = ['SPECTRUM-OV-107', 'SPECTRUM-OV-031', 'SPECTRUM-OV-003',
                    'SPECTRUM-OV-115', 'SPECTRUM-OV-112', 'SPECTRUM-OV-110',
                    'SPECTRUM-OV-105', 'SPECTRUM-OV-090', 'SPECTRUM-OV-083',
                    'SPECTRUM-OV-081', 'SPECTRUM-OV-082', 'SPECTRUM-OV-080',
                    'SPECTRUM-OV-077', 'SPECTRUM-OV-075', 'SPECTRUM-OV-071',
                    'SPECTRUM-OV-070', 'SPECTRUM-OV-068', 'SPECTRUM-OV-067',
                    'SPECTRUM-OV-065', 'SPECTRUM-OV-024', 'SPECTRUM-OV-025',
                    'SPECTRUM-OV-009', 'SPECTRUM-OV-014', 'SPECTRUM-OV-008',
                    'SPECTRUM-OV-037', 'SPECTRUM-OV-036', 'SPECTRUM-OV-054',
                    'SPECTRUM-OV-053', 'SPECTRUM-OV-002', 'SPECTRUM-OV-049',
                    'SPECTRUM-OV-052', 'SPECTRUM-OV-051', 'SPECTRUM-OV-050',
                    'SPECTRUM-OV-045', 'SPECTRUM-OV-042', 'SPECTRUM-OV-041',
                    'SPECTRUM-OV-026', 'SPECTRUM-OV-022', 'SPECTRUM-OV-007',
                    'SPECTRUM-OV-116', 'SPECTRUM-OV-118']

all_sample_ids = np.array(all_sample_ids)
chunks = np.array_split(range(len(all_sample_ids)),  5)

curr_samples = all_sample_ids[chunks[samp_group_idx]]


# read in the expression
adata = sc.read_h5ad(f"{data_path}/GSE180661_matrix_wAnnot.h5")

# index
idx_keep = adata.obs.donor_id.isin(curr_samples)
adata = adata[np.where(idx_keep)[0], :]


# write the chunk
adata.write(result_ad_file)
