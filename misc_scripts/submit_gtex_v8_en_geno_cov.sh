for tissue in `cat gtex_v8_tissue_list.txt`
do
  qsub -v TISSUE=$tissue -N en_$tissue gen_gtex_v8_en_geno_cov.qsub
done