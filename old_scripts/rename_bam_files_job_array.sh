#!/bin/bash

## SETUP
# set job parameters
#SBATCH --job-name=rename_bam
#SBATCH --output=rename_bam_%A_%a.out
#SBATCH --array=1-6
#SBATCH --ntasks=1
#SBATCH --mem=8G

# assign directory variable
dir="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200722_HeLa_screen/bam_sorted/"

## gRNA1 commands
# assign variables
gRNA1_bam_in="${dir}PP_pgRNA_HeLa_S1_R1_001_trimmed_sample${SLURM_ARRAY_TASK_ID}_aligned_sorted.bam"
gRNA1_bam_out="${dir}PP_pgRNA_HeLa_S1_R1_001_trimmed_sample${SLURM_ARRAY_TASK_ID}_aligned_sorted_flagstat.txt"

# run samtools
samtools flagstat "$gRNA1_bam" > "$gRNA1_bam_stats"


## gRNA2 commands
# assign variables
gRNA2_bam="${dir}PP_pgRNA_HeLa_S1_I1_001_trimmed_sample${SLURM_ARRAY_TASK_ID}_aligned_sorted.bam"
gRNA2_bam_stats="${dir}PP_pgRNA_HeLa_S1_I1_001_trimmed_sample${SLURM_ARRAY_TASK_ID}_aligned_sorted_flagstat.txt"

# run samtools
samtools flagstat "$gRNA2_bam" > "$gRNA2_bam_stats"
