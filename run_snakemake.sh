#!/bin/bash

BASE_PATH=$(pwd)

# bash make_dirs.sh "${BASE_PATH}"

SNAKE_FILE="${BASE_PATH}/workflow/Snakefile"
CONFIG_FILE="${BASE_PATH}/config/config.yaml"
CONDA_ENV="${BASE_PATH}/workflow/envs/snakemake.yaml"

# ## make report dirs
REPORT_DIR="${BASE_PATH}/workflow/report"
mkdir -p "${REPORT_DIR}"

## run this if the Snakemake environment has not been set up yet
# mamba env create -f $CONDA_ENV

## activate conda envt
source activate snakemake && echo "snakemake env activated"

## check if idemp successfully installed
IDEMP=config/idemp/idemp
if test -f "$IDEMP"; then
    echo "idemp has been successfully installed and compiled in config folder."
else 
    echo "$IDEMP does not exist. Installing idemp."
    ## install idemp in config dir
    git -C config clone --quiet https://github.com/yhwu/idemp.git
    echo "Compiling idemp."
    ## compile idemp
    make -C config/idemp  
fi

## make output dirs
mkdir -p "${BASE_PATH}/results/fastq_trimmed"
mkdir -p "${BASE_PATH}/results/fastq_demuxed"
mkdir -p "${BASE_PATH}/results/sam"
mkdir -p "${BASE_PATH}/results/bam"
mkdir -p "${BASE_PATH}/results/bam_sorted"
mkdir -p "${BASE_PATH}/results/pgRNA_counts"

## make log dirs
mkdir -p "workflow/logs/trim_reads"
mkdir -p "workflow/logs/demux_fastqs"
mkdir -p "workflow/logs/align_reads"
mkdir -p "workflow/logs/make_sorted_bams"
mkdir -p "workflow/logs/get_stats"
mkdir -p "workflow/logs/count_pgRNAs"
mkdir -p "workflow/logs/combine_counts"

# ## run the pipeline
snakemake --snakefile "workflow/Snakefile"  \
  --use-conda --conda-prefix "~/tmp/" --conda-frontend mamba \
  -k -p --reason --jobs 50 --latency-wait 80 -n

  ## -n = dry run, -r prints reasons
  ## -R forces snakemake to rerun from a specific rule
  # --cluster "sbatch &> ${SBATCH_OUT}%j.out" \
  # --use-conda --conda-prefix $CONDA_ENV \
  # --use-conda --conda-prefix "~/tmp/"
  # --restart-times 3 --conda-prefix $CONDA_ENV
  # force rerun: -R trim_reads

echo -e "\nDONE!\n"

# echo -e "Exporting pipeline rulegraph to svg and pdf..."
# snakemake --snakefile $SNAKE_FILE \
#   --configfile $CONFIG_FILE \
#   --rulegraph > "${REPORT_DIR}/rulegraph.dot"
# dot -Tpdf "${REPORT_DIR}/rulegraph.dot" > "${REPORT_DIR}/pipeline_rulegraph.pdf"
# dot -Tsvg "${REPORT_DIR}/rulegraph.dot" > "${REPORT_DIR}/pipeline_rulegraph.svg"
# rm "${REPORT_DIR}/rulegraph.dot"

# ## export PDF and svg visualizations of DAG structure of pipeline steps
# ## NOTE: this will not work if there are print statements in the pipeline
# echo -e "Exporting pipeline DAG to svg and pdf..."
# snakemake --snakefile $SNAKE_FILE \
#   --configfile $CONFIG_FILE \
#   --dag > "${REPORT_DIR}/dag.dot"
# dot -Tpdf "${REPORT_DIR}/dag.dot" > "${REPORT_DIR}/pipeline_dag.pdf"
# dot -Tsvg "${REPORT_DIR}/dag.dot" > "${REPORT_DIR}/pipeline_dag.svg"
# rm "${REPORT_DIR}/dag.dot"

# echo -e "\nDONE!\n"
