#!/bin/bash

## configure file paths
CONFIG_FILE="/home/pparrish/pgPEN_pipeline/config/AB_config.yaml"
SNAKE_FILE="/home/pparrish/pgPEN_pipeline/workflow/Snakefile"
SBATCH_OUT="/home/pparrish/pgPEN_pipeline/workflow/logs/pipeline/sbatch_out/slurm-"
CONDA_ENV="/home/pparrish/pgPEN_pipeline/workflow/envs/snakemake.yaml"

## activate conda envt
source activate snakemake

## run the pipeline
##   -p prints out shell commands, -k keeps going w/ independent jobs
##   if one fails
snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE \
  -k -p --reason --jobs 50 --use-conda --latency-wait 180 \
  --cluster "sbatch -e ${SBATCH_OUT}%j.err -o ${SBATCH_OUT}%j.out"
  # --use-conda --conda-prefix $CONDA_ENV \
  # --restart-times 3 --conda-prefix $CONDA_ENV
  # force rerun: -R trim_reads

## export PDF and svg visualizations of DAG structure of pipeline steps
## NOTE: this will not work if there are print statements in the pipeline
# echo -e "Exporting pipeline DAG to svg and pdf..."
# snakemake --snakefile $SNAKE_FILE \
#   --configfile $CONFIG_FILE \
#   --dag > images/dag.dot
# dot -Tpdf images/dag.dot > images/pipeline_dag.pdf
# dot -Tsvg images/dag.dot > images/pipeline_dag.svg
# rm images/dag.dot
#
# echo -e "Generating pipeline HTML report..."
# snakemake --configfile $CONFIG_FILE \
#   --snakefile $SNAKE_FILE \
#   --report workflow/logs/pipeline/report.html
#
# echo -e "\nDONE!\n"
