db_file=/gpfs/data/im-lab/nas40t2/Data/PredictDB/GTEx_v8/models_v1/eqtl/ctimp/ctimp_Whole_Blood.db

# run s-predixcan
conda activate imlabtools
SNMK=/gpfs/data/im-lab/nas40t2/yanyul/softwares/miniconda2/envs/snakemake/bin/snakemake
snmkscript=/gpfs/data/im-lab/nas40t2/yanyul/GitHub/misc-tools/spredixcan_workflow/spredixcan_pipe.snmk

$SNMK -s $snmkscript --configfile config.spxcan.yaml -np