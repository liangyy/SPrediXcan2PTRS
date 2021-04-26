#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb
#PBS -e logs/run_${TISSUE}_$GWASTAG.err
#PBS -o logs/run_${TISSUE}_$GWASTAG.out

# ARGS:
# TISSUE
# GWASTAG
# GWASN

if [[ -z $TISSUE ]]
then
  TISSUE=$1
  GWASTAG=$2
  GWASN=$3
  PBS_O_WORKDIR=`pwd`
fi

source ~/.bash_profile
source ~/.bashrc

predict_db=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/ctimp/ctimp_${TISSUE}.db
predict_db_cov=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/ctimp/ctimp_${TISSUE}.txt.gz
gwas=/gpfs/data/im-lab/nas40t2/Data/SummaryResults/imputed_gwas_hg38_1.1/imputed_${GWASTAG}.txt.gz
outdir=/scratch/t.cri.yliang/SPrediXcan2PTRS/run_gtex_gwas


export PYTHONPATH=$PYTHONPATH:/home/t.cri.yliang/labshare/GitHub/SPrediXcan2PTRS
export PYTHONPATH=$PYTHONPATH:/home/t.cri.yliang/labshare/GitHub/transethnic_prs

conda activate imlabtools

# impute beta and se from z
imputeb_gwas=$outdir/imputed_bhat.${GWASTAG}.txt.gz
if [[ ! -f $imputeb_gwas ]]
then
  echo "Imputing effect size of GWAS"
  echo "Input: $gwas"
  echo "Output: $imputeb_gwas"
  python $PBS_O_WORKDIR/impute_b_for_gwas.py \
    --input $gwas \
    --zscore zscore \
    --freq frequency \
    --sample_size sample_size \
    --output $imputeb_gwas
fi

# run s-predixcan
spxcanscript=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/MetaXcan/software/SPrediXcan.py
pxcan_file=$outdir/spredixcan.${GWASTAG}.${TISSUE}.csv
if [[ ! -f pxcan_file ]]
then
  echo "Running S-PrediXcan"
  echo "Input: $imputeb_gwas"
  echo "Output: $pxcan_file"
  python $spxcanscript \
    --gwas_file $imputeb_gwas \
    --snp_column panel_variant_id \
    --effect_allele_column effect_allele \
    --non_effect_allele_column non_effect_allele \
    --beta_column effect_size \
    --se_column standard_error \
    --model_db_path $predict_db \
    --covariance $predict_db_cov \
    --keep_non_rsid \
    --additional_output \
    --model_db_snp_key varID \
    --throw \
    --output_file $pxcan_file
fi
  
  
# run SPrediXcan2PTRS
conda deactivate
conda activate SPrediXcan2PTRS


runscript=/home/t.cri.yliang/labshare/GitHub/SPrediXcan2PTRS/run_pxcan2ptrs.py

geno_cov_file=/home/t.cri.yliang/scratch/SPrediXcan2PTRS/geno_cov/gtex_v8/ctimp_$TISSUE.geno_cov.chr{chr_num}.evd.npz

ptrs_prefix=$outdir/spxcan2ptrs.${GWASTAG}.${TISSUE}
ptrs_file=$ptrs_prefix.results.h5

if [[ -f $ptrs_file ]]
then
  echo "Running SPrediXcan2PTRS"
  echo "Input: $pxcan_file"
  echo "Output: $ptrs_file"
  python $runscript \
    --predixcan $pxcan_file \
    --predictdb $predict_db \
    --geno_cov $geno_cov_file \
    --gwas $gwas \
    --gwas_cols chromosome=chromosome \
      position=position \
      effect_allele=effect_allele \
      non_effect_allele=non_effect_allele \
    --gwas_sample_size $GWASN \
    --output_prefix $ptrs_prefix
fi

