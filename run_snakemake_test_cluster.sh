#!/bin/bash

## configure file paths
## make BASE_PATH a variable entered by the user so run_snakemake can be used for any
## test_snakemake path: /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/test_snakemake
## PC9 screen path: /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgPEN_library/200122_PC9_screen
# BASE_PATH="$1" ## gets user input path
BASE_PATH=$(pwd)
echo "base_path=${BASE_PATH}"
CONFIG_FILE="${BASE_PATH}/config/config.yaml"
SNAKE_FILE="${BASE_PATH}/workflow/Snakefile"
SBATCH_OUT="${BASE_PATH}/workflow/logs/sbatch/slurm-"
CONDA_ENV="${BASE_PATH}/workflow/envs/snakemake.yaml"
REPORT_DIR="${BASE_PATH}/workflow/report"

## run this if the Snakemake environment has not been set up yet
# mamba create -f $CONDA_ENV

## activate conda envt
source activate snakemake
echo "snakemake env activated"

# ## make log directories
mkdir -p "${BASE_PATH}/workflow/logs/trim_reads"
mkdir -p "${BASE_PATH}/workflow/logs/demux_fastqs"
mkdir -p "${BASE_PATH}/workflow/logs/build_bowtie_index"
mkdir -p "${BASE_PATH}/workflow/logs/align_reads"
mkdir -p "${BASE_PATH}/workflow/logs/make_sorted_bams"
mkdir -p "${BASE_PATH}/workflow/logs/get_stats"
mkdir -p "${BASE_PATH}/workflow/logs/count_pgRNAs"
mkdir -p "${BASE_PATH}/workflow/logs/combine_counts"

## make output dirs
# mkdir -p "${BASE_PATH}/results/fastq"
# mkdir -p "${BASE_PATH}/results/fastq_trimmed"
# mkdir -p "${BASE_PATH}/results/fastq_demuxed"
# mkdir -p "${BASE_PATH}/results/sam"
# mkdir -p "${BASE_PATH}/results/bam"
# mkdir -p "${BASE_PATH}/results/bam_sorted"
# mkdir -p "${BASE_PATH}/results/pgRNA_counts"

# ## make report dirs
# mkdir -p "${BASE_PATH}/workflow/logs/sbatch"
# mkdir -p "${REPORT_DIR}"
#
# ## run the pipeline
# ##   -p prints out shell commands, -k keeps going w/ independent jobs
# ##   if one fails
snakemake --snakefile $SNAKE_FILE \
  --use-conda \
  --conda-prefix "~/tmp/" \
  --conda-frontend mamba \
  -k -p --reason \
  --jobs 50 --latency-wait 180 --restart-times 3 \
  -R trim_reads \
  --cluster-config config/cluster_slurm.yaml \
  --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {cluster.ncpus} -n {cluster.ntasks} -o {cluster.output}" \
  ## -n = dry run, -r prints reasons
  ## -R forces snakemake to rerun from a specific rule
  # --cluster "sbatch &> ${SBATCH_OUT}%j.out" \
  # --use-conda --conda-prefix $CONDA_ENV \
  # force rerun: -R trim_reads

## running snakemake on the cluster:
## - add in --cluster "sbatch {threads} --mem={resources.mem}" etc

echo -e "\nDONE!\n"
