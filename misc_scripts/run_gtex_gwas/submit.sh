# args1: tissue
tissue=$1

while IFS=' ' read -r col1 col2
do 
  qsub -v TISSUE=$tissue,GWASTAG=$col1,GWASN=$col2 -N $tissue-$col1 run.sh
done < trait_list.txt
