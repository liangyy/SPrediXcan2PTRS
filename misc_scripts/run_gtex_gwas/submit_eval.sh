# args1: tissue
# args2: nsample
tissue=$1
nsample=$2

mkdir -p logs

mkdir -p /scratch/t.cri.yliang/SPrediXcan2PTRS/run_gtex_gwas

while IFS=$'\t' read -r col1 col2 col3
do 
  echo qsub -v TISSUE=$tissue,GWASTAG=$col1,NSAMPLE=$nsample -N $col1-$nsample eval.sh
done < trait_list.txt
