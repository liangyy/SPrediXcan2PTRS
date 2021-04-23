#PBS -S /bin/bash
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=1
#PBS -l mem=32gb
#PBS -e run_ukb_gwas_to_ptrs.err
#PBS -o run_ukb_gwas_to_ptrs.out
#PBS -N run_ukb_gwas_to_ptrs


source ~/.bash_profile
source ~/.bashrc

# run s-predixcan
conda activate imlabtools
SNMK=/gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/mixqtl/bin/snakemake
snmkscript=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/spredixcan_workflow/spredixcan_pipe.snmk


pxcan_file=/scratch/t.cri.yliang/SPrediXcan2PTRS/test_ptrs_ukb_gwas/spredixcan.height_x_Whole_Blood.harmonization_full.csv 
gwas_file=/scratch/t.cri.yliang/SPrediXcan2PTRS/test_ptrs_ukb_gwas/height.harmonization_full.txt.gz
# db_file=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/ctimp/ctimp_Whole_Blood.db
if [[ ! -f $pxcan_file ]]
then
  $SNMK -s $snmkscript --configfile ${PBS_O_WORKDIR}/config.spxcan.yaml -p --cores 1 \
    > ${PBS_O_WORKDIR}/run_ukb_gwas_to_ptrs.spredixcan.log 2>&1
fi

# run SPrediXcan2PTRS

conda deactivate
conda activate SPrediXcan2PTRS

export PYTHONPATH=$PYTHONPATH:/home/t.cri.yliang/labshare/GitHub/SPrediXcan2PTRS
export PYTHONPATH=$PYTHONPATH:/home/t.cri.yliang/labshare/GitHub/transethnic_prs

runscript=/home/t.cri.yliang/labshare/GitHub/SPrediXcan2PTRS/run_pxcan2ptrs.py

geno_cov_file=/home/t.cri.yliang/scratch/SPrediXcan2PTRS/geno_cov/gtex_v8/ctimp_Whole_Blood.geno_cov.chr{chr_num}.evd.npz
db_file=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/ctimp/ctimp_Whole_Blood.db

outdir=/scratch/t.cri.yliang/SPrediXcan2PTRS/test_ptrs_ukb_gwas
outprefix=SPrediXcan2PTRS.height_x_Whole_Blood.harmonization_full 
outprefix2=SPrediXcan2PTRS.height_x_Whole_Blood.harmonization_full.2 
outprefix3=SPrediXcan2PTRS.height_x_Whole_Blood.harmonization_full.3

chainfile=/gpfs/data/im-lab/nas40t2/yanyul/data/MetaXcan/data/liftover/hg19ToHg38.over.chain.gz

# test 1
if [[ $1 == 1 ]]
then
  echo python $runscript \
    --predixcan $pxcan_file \
    --predictdb $db_file \
    --geno_cov $geno_cov_file \
    --gwas $gwas_file \
    --gwas_cols chromosome=chromosome \
      position=position \
      effect_allele=effect_allele \
      non_effect_allele=non_effect_allele \
    --gwas_sample_size 360388 \
    --output_prefix $outdir/$outprefix
fi

# test 2
if [[ $1 == 2 ]]
then
  echo python $runscript \
    --predixcan $pxcan_file \
    --predictdb $db_file \
    --geno_cov $geno_cov_file \
    --gwas $gwas_file \
    --gwas_cols variant={chromosome}:{position}:{non_effect_allele}:{effect_allele} \
    --liftover_chain $chainfile \
    --gwas_sample_size 360388 \
    --output_prefix $outdir/$outprefix2
fi

# test 3
if [[ $1 == 3 ]]
then
  echo python $runscript \
    --predixcan $pxcan_file \
    --predictdb $db_file \
    --geno_cov $geno_cov_file \
    --gwas $gwas_file \
    --gwas_cols chromosome=chromosome \
      position=position \
      effect_allele=effect_allele \
      non_effect_allele=non_effect_allele \
    --gwas_sample_size 360388 \
    --hyperparam_yaml params.yaml \
    --output_prefix $outdir/$outprefix3
fi

