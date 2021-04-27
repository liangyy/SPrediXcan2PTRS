# args1: tissue
tissue=$1

mkdir -p logs

mkdir -p /scratch/t.cri.yliang/SPrediXcan2PTRS/run_gtex_gwas

while IFS=' ' read -r col1 col2
do 
  echo qsub -v TISSUE=$tissue,GWASTAG=$col1,GWASN=$col2 -N $tissue-$col1 run.sh
done < trait_list.txt
