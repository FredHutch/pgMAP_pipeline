#!/bin/bash

## SETUP
# set job parameters
#SBATCH --job-name=bam_stats
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G

## general stuff
# load samtools
module load SAMtools/1.9-foss-2018b

# assign directory variable
dir="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200122_PC9_screen/"

## gRNA1 commands
# assign variables
gRNA1_bam="${dir}PP_pgRNA_PC9_S1_R1_001_trimmed_sample5_aligned_sorted.bam"
gRNA1_bam_stats="${dir}PP_pgRNA_PC9_S1_R1_001_trimmed_sample5_aligned_sorted_flagstat.txt"

# run samtools
samtools flagstat "$gRNA1_bam" > "$gRNA1_bam_stats"


## gRNA2 commands
# assign variables
gRNA2_bam="${dir}PP_pgRNA_PC9_S1_I1_001_trimmed_sample5_aligned_sorted.bam"
gRNA2_bam_stats="${dir}PP_pgRNA_PC9_S1_I1_001_trimmed_sample5_aligned_sorted_flagstat.txt"

# run samtools
samtools flagstat "$gRNA2_bam" > "$gRNA2_bam_stats"
