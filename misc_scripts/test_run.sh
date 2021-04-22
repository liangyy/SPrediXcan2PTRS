# conda activate SPrediXcan2PTRS

export PYTHONPATH=$PYTHONPATH:/home/t.cri.yliang/labshare/GitHub/SPrediXcan2PTRS
export PYTHONPATH=$PYTHONPATH:/home/t.cri.yliang/labshare/GitHub/transethnic_prs

runscript=/home/t.cri.yliang/labshare/GitHub/SPrediXcan2PTRS/run_pxcan2ptrs.py

pxcan_file=/home/t.cri.yliang/scratch/tcga-survival/spredixcan/spredixcan.BCAC_2020_HER2_x_Breast_Mammary_Tissue.harmonization_full.csv
geno_cov_file=/home/t.cri.yliang/scratch/SPrediXcan2PTRS/geno_cov/gtex_v8/en_Breast_Mammary_Tissue.geno_cov.chr{chr_num}.evd.npz
gwas_file=/home/t.cri.yliang/scratch/tcga-survival/spredixcan/BCAC_2020_HER2.harmonization_full.txt.gz
db_file=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/elastic_net_models/en_Breast_Mammary_Tissue.db

outdir=/gpfs/data/im-lab/nas40t2/yanyul/tmp/SPrediXcan2PTRS/tmp

python $runscript \
  --predixcan $pxcan_file \
  --predictdb $db_file \
  --geno_cov $geno_cov_file \
  --gwas $gwas_file \
         chromosome:chromosome \
         position:position \
         effect_allele:effect_allele \
         non_effect_allele:non_effect_allele \
  --gwas_sample_size 197735 \
  --output $outdir/test_run
