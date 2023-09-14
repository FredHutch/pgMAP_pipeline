#!/bin/bash

## un-comment the following line if the Snakemake conda environment 
##   has not been set up yet
# mamba env create -f workflow/envs/snakemake.yaml

## activate conda envt
source activate snakemake 
echo "snakemake env activated"

## set up folders and install idemp if needed
bash workflow/scripts/setup.sh

## run the pipeline
snakemake --snakefile "workflow/Snakefile"  \
  --use-conda --conda-prefix "~/tmp/" --conda-frontend mamba \
  -k -p --reason --jobs 50 --latency-wait 80 --restart-times 3

  ## other useful flags
  ## -n = dry run, -r prints reasons
  ## -R forces snakemake to rerun from a specific rule

echo -e "\npgMAP run completed! See results/pgRNA_counts/ folder for output.\n"

bash workflow/scripts/make_reports.sh