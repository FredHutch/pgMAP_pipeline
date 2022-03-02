#!/bin/bash

## sbatch parameters
#SBATCH --job-name=HeLa_counter
#SBATCH --output=HeLa_counter_%A_%a.out
#SBATCH --array=1-5
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=98G

## load appropriate R module with necessary packages
module load R/3.6.2-foss-2018b-fh1

## assign variables
script_dir="/home/pparrish/paralog_pgRNA_screen/scripts/"
bam_dir="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200722_HeLa_screen/bam_sorted/"
gRNA1_bam="${bam_dir}HeLa.sample${SLURM_ARRAY_TASK_ID}.gRNA_1.bam"
n_chunks=50

## run script
Rscript "${script_dir}counter_efficient.R" "${gRNA1_bam}" "${n_chunks}"
