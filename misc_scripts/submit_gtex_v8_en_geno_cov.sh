for tissue in `cat gtex_v8_tissue_list.txt`
do
  echo qsub -v TISSUE=$tissue -N en_$tissue gtex_v8_en_geno_cov.qsub
done
