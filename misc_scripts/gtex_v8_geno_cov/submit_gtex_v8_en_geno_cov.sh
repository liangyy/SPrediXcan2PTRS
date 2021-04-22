# args1: tag=en|ctimp
tag=$1
for tissue in `cat gtex_v8_tissue_list.txt`
do
  echo qsub -v TISSUE=$tissue -N ${tag}_$tissue gtex_v8_${tag}_geno_cov.qsub
done
