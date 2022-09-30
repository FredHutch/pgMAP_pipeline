#!/bin/bash
  
## configure file paths
## make cwd a variable entered by the user so run_snakemake can be used for any
# cwd="/fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/test_snakemake"

cwd=$(pwd)

# echo "base_path=${cwd}"
CONFIG_FILE="${cwd}/config/config.yaml"
SNAKE_FILE="${cwd}/workflow/Snakefile"
SBATCH_OUT="${cwd}/workflow/logs/sbatch/slurm-"
CONDA_ENV="${cwd}/workflow/envs/snakemake.yaml"
REPORT_DIR="${cwd}/workflow/report"

# mamba create -f $CONDA_ENV

## activate conda envt
source activate snakemake

## make log directories
mkdir -p "${cwd}/workflow/logs/trim_fastqs"
mkdir -p "${cwd}/workflow/logs/demux_fastqs"
mkdir -p "${cwd}/workflow/logs/align_reads"
mkdir -p "${cwd}/workflow/logs/make_sorted_bams"
mkdir -p "${cwd}/workflow/logs/get_stats"
mkdir -p "${cwd}/workflow/logs/count_pgRNAs"
mkdir -p "${cwd}/workflow/logs/combine_counts"
mkdir -p "${cwd}/workflow/logs/sbatch"
mkdir -p "${REPORT_DIR}"

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
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         
~                                                                                         

