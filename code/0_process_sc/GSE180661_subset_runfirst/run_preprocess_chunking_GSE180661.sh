#!/bin/bash

indir=/YOUR_PATH/data/single_cell/GSE180661/

python preprocess_chunking_0.py --indir ${indir}

idx=0
echo "-------------------"
echo ${idx}
echo "-------------------"
python preprocess_chunking_1.py --indir ${indir} --index ${idx}
python preprocess_chunking_2.py --indir ${indir} --index ${idx}
#rm ${indir}/GSE180661_subset_${idx}.h5ad


idx=1
echo "-------------------"
echo ${idx}
echo "-------------------"
python preprocess_chunking_1.py --indir ${indir} --index ${idx}
python preprocess_chunking_2.py --indir ${indir} --index ${idx}
#rm ${indir}/GSE180661_subset_${idx}.h5ad


idx=2
echo "-------------------"
echo ${idx}
echo "-------------------"
python preprocess_chunking_1.py --indir ${indir} --index ${idx}
python preprocess_chunking_2.py --indir ${indir} --index ${idx}
#rm ${indir}/GSE180661_subset_${idx}.h5ad


idx=3
echo "-------------------"
echo ${idx}
echo "-------------------"
python preprocess_chunking_1.py --indir ${indir} --index ${idx}
python preprocess_chunking_2.py --indir ${indir} --index ${idx}
#rm ${indir}/GSE180661_subset_${idx}.h5ad

idx=4
echo "-------------------"
echo ${idx}
echo "-------------------"
python preprocess_chunking_1.py --indir ${indir} --index ${idx}
python preprocess_chunking_2.py --indir ${indir} --index ${idx}
#rm ${indir}/GSE180661_subset_${idx}.h5ad


python preprocess_chunking_3.py --indir ${indir}

#rm ${indir}/GSE180661_subset2_0.h5ad
#rm ${indir}/GSE180661_subset2_1.h5ad
#rm ${indir}/GSE180661_subset2_2.h5ad
#rm ${indir}/GSE180661_subset2_3.h5ad
#rm ${indir}/GSE180661_subset2_4.h5ad
