#!/bin/bash

## configure file paths
BASE_PATH="/fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/200122_PC9_screen"
CONFIG_FILE="${BASE_PATH}/config/config.yaml"
SNAKE_FILE="${BASE_PATH}/workflow/Snakefile"
SBATCH_OUT="${BASE_PATH}/workflow/logs/sbatch/slurm-"
CONDA_ENV="${BASE_PATH}/workflow/envs/snakemake.yaml"
REPORT_DIR="${BASE_PATH}/workflow/report"

## activate conda envt
source activate snakemake

## run the pipeline
##   -p prints out shell commands, -k keeps going w/ independent jobs
##   if one fails
snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE \
  --use-conda --conda-prefix "~/tmp/" --conda-frontend mamba \
  -k -p --reason --jobs 50 --latency-wait 180 \
  --cluster "sbatch &> ${SBATCH_OUT}%j.out"
  # --use-conda --conda-prefix $CONDA_ENV \
  # --restart-times 3 --conda-prefix $CONDA_ENV
  # force rerun: -R trim_reads

## export PDF and svg visualizations of DAG structure of pipeline steps
## NOTE: this will not work if there are print statements in the pipeline
echo -e "Exporting pipeline DAG to svg and pdf..."
snakemake --snakefile $SNAKE_FILE \
  --configfile $CONFIG_FILE \
  --dag > "${REPORT_DIR}/dag.dot"
dot -Tpdf "${REPORT_DIR}/dag.dot" > "${REPORT_DIR}/pipeline_dag.pdf"
dot -Tsvg "${REPORT_DIR}/dag.dot" > "${REPORT_DIR}/pipeline_dag.svg"
rm "${REPORT_DIR}/dag.dot"

echo -e "Generating pipeline HTML report..."
snakemake --configfile $CONFIG_FILE \
  --snakefile $SNAKE_FILE \
  --report "${REPORT_DIR}/report.html"

echo -e "\nDONE!\n"
