#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb
#PBS -e logs/eval_${TISSUE}_$GWASTAG_$NSAMPLE.err
#PBS -o logs/eval_${TISSUE}_$GWASTAG_$NSAMPLE.out

# ARGS:
# TISSUE
# GWASTAG
# NSAMPLE
# SEED (optional)

if [[ -z $TISSUE ]]
then
  TISSUE=$1
  GWASTAG=$2
  NSAMPLE=$3
  # SEED=$4
  PBS_O_WORKDIR=`pwd`
fi

if [[ $NSAMPLE == "all" ]]
then
  nsample_cmd=""
else
  nsample_cmd="--n_samples $NSAMPLE"
fi

source ~/.bash_profile
source ~/.bashrc

conda activate SPrediXcan2PTRS

outdir=/scratch/t.cri.yliang/SPrediXcan2PTRS/run_gtex_gwas
pxcan_file=$outdir/spredixcan.${GWASTAG}.${TISSUE}.csv
ptrs_prefix=$outdir/spxcan2ptrs.${GWASTAG}.${TISSUE}
ptrs_file=$ptrs_prefix.results.h5
pred_expr=/home/t.cri.yliang/labshare/PTRS/predicted_expression/predicted_expression.ukb_imp_x_ctimp_${TISSUE}.h5
sample_list=${PBS_O_WORKDIR}/sample_list_for_eval.txt

if [[ ! -f $sample_list ]]
then
  echo "Generating sample list"
  echo "Input: ..."
  echo "Output: $sample_list"
  python $PBS_O_WORKDIR/gen_sample_list_for_eval.py
fi

eval_file=$outdir/eval.${GWASTAG}.${TISSUE}.parquet

if [[ ! -f $eval_file ]]
then
  echo "Running evaluation."
  echo "Input: ..."
  echo "Output: $eval_file"
  python $PBS_O_WORKDIR/eval_on_ukb.py \
    --pred_expr $pred_expr \
    --sample_list $sample_list \
    $nsample_cmd \
    --list_of_ptrs ${GWASTAG}:$ptrs_file \
    --list_of_pxcan ${GWASTAG}:$pxcan_file \
    --output_parquet $eval_file
fi

# --seed $SEED \

