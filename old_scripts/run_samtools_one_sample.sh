#!/bin/bash

## SETUP
# set job parameters
#SBATCH --job-name=bam_sort
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=16G

# load samtools
module load SAMtools/1.9-foss-2018b

## general stuff
# assign directory variable
dir="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200122_PC9_screen/"

## gRNA 1 stuff
# assign variables
gRNA1_sam="${dir}PP_pgRNA_PC9_S1_R1_001_trimmed_sample5_aligned.sam"
gRNA1_bam="${dir}PP_pgRNA_PC9_S1_R1_001_trimmed_sample5_aligned.bam"
gRNA1_sorted_bam="${dir}PP_pgRNA_PC9_S1_R1_001_trimmed_sample5_aligned_sorted.bam"
gRNA1_tmp="${dir}tmp1_5"

# run samtools
samtools view -bS -o "$gRNA1_bam" "$gRNA1_sam"

samtools sort -O bam -n "$gRNA1_bam" -T "$gRNA1_tmp" -o "$gRNA1_sorted_bam"


## gRNA2 stuff
# assign variables
gRNA2_sam="${dir}PP_pgRNA_PC9_S1_I1_001_trimmed_sample5_aligned.sam"
gRNA2_bam="${dir}PP_pgRNA_PC9_S1_I1_001_trimmed_sample5_aligned.bam"
gRNA2_sorted_bam="${dir}PP_pgRNA_PC9_S1_I1_001_trimmed_sample5_aligned_sorted.bam"
gRNA2_tmp="${dir}tmp2_5"

# run samtools
samtools view -bS -o "$gRNA2_bam" "$gRNA2_sam"

samtools sort -O bam -n "$gRNA2_bam" -T "$gRNA2_tmp" -o "$gRNA2_sorted_bam"


## Notes on Samtools options:
## view
# -bS: outputs to BAM format; S is no longer used but specifies SAM input format
## sort
# -n: sorts by read names rather than chromosomal coordinates
