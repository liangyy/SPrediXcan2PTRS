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


outfile=/scratch/t.cri.yliang/SPrediXcan2PTRS/test_ptrs_ukb_gwas/spredixcan.height_x_Whole_Blood.harmonization_full.csv 
# db_file=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/ctimp/ctimp_Whole_Blood.db
if [[ ! -f $outfile ]]
then
  $SNMK -s $snmkscript --configfile ${PBS_O_WORKDIR}/config.spxcan.yaml -p --cores 1 \
    > ${PBS_O_WORKDIR}/run_ukb_gwas_to_ptrs.spredixcan.log 2>&1
fi

# run SPrediXcan2PTRS
# TODO
