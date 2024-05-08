cat ~/PROJECT/NMF/samples|while read i
do
echo "#!/bin/bash

module load anaconda/4.12.0 
conda activate cnmf-env
cd /share/home/NMF/

# step 1. prepare ----
python /share/home/luoylLab/zengyuchen/biosoft/cNMF-master/src/cnmf/cnmf.py prepare --output-dir ./res1/ --name ${i}_cNMF -c ./${i}count.txt -k 3 4 5 6 7 8 9 10 --n-iter 300 --total-workers 3 --numgenes 2000

# step 2. factorize ----
python -W ignore /share/home/luoylLab/zengyuchen/biosoft/cNMF-master/src/cnmf/cnmf.py factorize --output-dir ./res1/ --name ${i}_cNMF --worker-index 0

# step 3. combine ----
python /share/home/luoylLab/zengyuchen/biosoft/cNMF-master/src/cnmf/cnmf.py combine --output-dir ./res1/ --name ${i}_cNMF
rm -f ./res1/001N_cNMF/cnmf_tmp/${i}_cNMF.spectra.k_*.iter_*.df.npz

# step 4. k_selection_plot ----
python /share/home/luoylLab/zengyuchen/biosoft/cNMF-master/src/cnmf/cnmf.py k_selection_plot --output-dir ./res1/ --name ${i}_cNMF 

done



