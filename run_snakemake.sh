#!/bin/bash

## configure file paths
CONFIG_FILE="/home/pparrish/pgPEN_pipeline/config/config.yaml"
SNAKE_FILE="/home/pparrish/pgPEN_pipeline/workflow/Snakefile"
# CONDA_ENV="/home/pparrish/pgPEN_pipeline/environment.yaml"

## activate conda envt
source activate snakemake_env

## run the pipeline
##   -p prints out shell commands, -k keeps going w/ independent jobs
##   if one fails
snakemake --snakefile $SNAKE_FILE \
  --configfile $CONFIG_FILE \
  --use-conda -k -p --jobs 50 --reason
  # --cluster "sbatch -o {log}" --jobs 50 \
  # --use-conda \ # --conda-prefix $CONDA_ENV
  # -k -p --restart-times 3 --latency-wait 180

## export PDF and svg visualizations of DAG structure of pipeline steps
## NOTE: this will not work if there are print statements in the pipeline
# echo -e "Exporting pipeline DAG to svg and pdf..."
# snakemake --snakefile $SNAKE_FILE \
#   --configfile $CONFIG_FILE \
#   --dag > images/dag.dot
# dot -Tpdf dag.dot > images/pipeline_dag.pdf
# dot -Tsvg dag.dot > images/pipeline_dag.svg
# rm dag.dot
#
# echo -e "Generating pipeline HTML report..."
# snakemake --configfile $CONFIG_FILE \
#   --snakefile $SNAKE_FILE \
#   --report workflow/logs/pipeline/report.html
#
# echo -e "\nDONE!\n"
