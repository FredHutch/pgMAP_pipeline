#!/bin/bash

## configure file paths
## make BASE_PATH a variable entered by the user so run_snakemake can be used for any
## test_snakemake path: /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/test_snakemake
## PC9 screen path: /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/200122_PC9_screen
# BASE_PATH="$1" ## gets user input path
# BASE_PATH=$(pwd)
# echo "base_path=${BASE_PATH}"
# CONFIG_FILE="${BASE_PATH}/config/config.yaml"
# SNAKE_FILE="${BASE_PATH}/workflow/Snakefile"
# SBATCH_OUT="${BASE_PATH}/workflow/logs/sbatch/slurm-"
# CONDA_ENV="${BASE_PATH}/workflow/envs/snakemake.yaml"
# REPORT_DIR="${BASE_PATH}/workflow/report"

## run this if the Snakemake environment has not been set up yet
# mamba create -f $CONDA_ENV

## activate conda envt
source activate snakemake
echo "snakemake env activated"

# ## make log directories
# mkdir -p "workflow/logs/trim_reads"
# mkdir -p "workflow/logs/demux_fastqs"
# mkdir -p "workflow/logs/build_bowtie_index"
# mkdir -p "workflow/logs/align_reads"
# mkdir -p "workflow/logs/make_sorted_bams"
# mkdir -p "workflow/logs/get_stats"
# mkdir -p "workflow/logs/count_pgRNAs"
# mkdir -p "workflow/logs/combine_counts"


## running snakemake on the cluster:
## - add in --cluster "sbatch {threads} --mem={resources.mem}" etc

# snakemake --snakefile "workflow/Snakefile"  \
#   --use-conda --conda-prefix "~/tmp/" --conda-frontend mamba

echo -e "Exporting pipeline rulegraph to svg and pdf..."
snakemake --snakefile "workflow/Snakefile" \
  --rulegraph > "workflow/report/rulegraph.dot"
dot -Tpdf "workflow/report/rulegraph.dot" > "workflow/report/pipeline_rulegraph.pdf"
dot -Tsvg "workflow/report/rulegraph.dot" > "workflow/report/pipeline_rulegraph.svg"
rm "workflow/report/rulegraph.dot"

## export PDF and svg visualizations of DAG structure of pipeline steps
## NOTE: this will not work if there are print statements in the pipeline
echo -e "Exporting pipeline DAG to svg and pdf..."
snakemake --snakefile "workflow/Snakefile" \
  --dag > "workflow/report/dag.dot"
dot -Tpdf "workflow/report/dag.dot" > "workflow/report/pipeline_dag.pdf"
dot -Tsvg "workflow/report/dag.dot" > "workflow/report/pipeline_dag.svg"
rm "workflow/report/dag.dot"

echo -e "\nDONE!\n"
