#!/bin/bash

## configure file paths
## make BASE_PATH a variable entered by the user so run_snakemake can be used for any
## test_snakemake path: /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/test_snakemake
## PC9 screen path: /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/200122_PC9_screen
BASE_PATH="$1" ## gets user input path
# echo "base_path=${BASE_PATH}"
CONFIG_FILE="${BASE_PATH}/config/config.yaml"
SNAKE_FILE="${BASE_PATH}/workflow/Snakefile"
SBATCH_OUT="${BASE_PATH}/workflow/logs/sbatch/slurm-"
CONDA_ENV="${BASE_PATH}/workflow/envs/snakemake.yaml"
REPORT_DIR="${BASE_PATH}/workflow/report"

# mamba create -f $CONDA_ENV

## activate conda envt
source activate snakemake

## make log and report directories
mkdir -p "${BASE_PATH}/workflow/logs/trim_fastqs"
mkdir -p "${BASE_PATH}/workflow/logs/demux_fastqs"
mkdir -p "${BASE_PATH}/workflow/logs/align_reads"
mkdir -p "${BASE_PATH}/workflow/logs/make_sorted_bams"
mkdir -p "${BASE_PATH}/workflow/logs/get_stats"
mkdir -p "${BASE_PATH}/workflow/logs/count_pgRNAs"
mkdir -p "${BASE_PATH}/workflow/logs/combine_counts"
mkdir -p "${BASE_PATH}/workflow/logs/sbatch"
mkdir -p "${REPORT_DIR}"

## run the pipeline
##   -p prints out shell commands, -k keeps going w/ independent jobs
##   if one fails
snakemake --snakefile $SNAKE_FILE --configfile $CONFIG_FILE \
  --use-conda --conda-prefix "~/tmp/" --conda-frontend mamba \
  -k -p --reason --jobs 50 --latency-wait 180 \
  --cluster "sbatch &> ${SBATCH_OUT}%j.out" \
  -R trim_reads ## this forces snakemake to rerun from rule trim_reads
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
