#!/bin/bash

#SBATCH --job-name=align
#SBATCH --output=align_%A_%a.out
#SBATCH --array=1-2
#SBATCH --ntasks=1

module load Bowtie/1.2.2-foss-2018b

## GUIDE 1
# set gRNA1 file name variables
gRNA1_idx_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_indices/pgRNA1"
gRNA1_infile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/190522_plasmid/Undetermined.R1_sample${SLURM_ARRAY_TASK_ID}_trimmed.fastq"
gRNA1_outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/test/plasmid_R1_test3_sample${SLURM_ARRAY_TASK_ID}_aligned.sam"

bowtie -q -v 1 --best --strata --all --sam -p 4 \
  "$gRNA1_idx_file" \
  "$gRNA1_infile" \
  "$gRNA1_outfile"


## GUIDE 2
# set gRNA2 file name variables
gRNA2_idx_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_indices/pgRNA2"
gRNA2_infile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/190522_plasmid/Undetermined.R2_sample${SLURM_ARRAY_TASK_ID}_trimmed.fastq"
gRNA2_outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/test/plasmid_R2_test3_sample${SLURM_ARRAY_TASK_ID}_aligned.sam"

bowtie -q -v 1 --best --strata --all --sam -p 4 \
  "$gRNA2_idx_file" \
  "$gRNA2_infile" \
  "$gRNA2_outfile"

