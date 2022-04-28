#!/bin/bash

## SETUP
# set job parameters
#SBATCH --job-name=align
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G

# load bowtie1
module load Bowtie/1.2.2-foss-2018b

## GUIDE 1
# set gRNA1 file name variables
# pretty sure the FASTQ has to be unzipped
gRNA1_idx_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_indices/pgRNA1"
gRNA1_infile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/demultiplex/200122_PC9_screen/PP_pgRNA_PC9_S1_R1_001_trimmed.fastq_sample5.fastq"
gRNA1_outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200122_PC9_screen/PP_pgRNA_PC9_S1_R1_001_trimmed_sample5_aligned.sam"

bowtie -q -v 1 --best --strata --all --sam -p 4 \
  "$gRNA1_idx_file" \
  "$gRNA1_infile" \
  "$gRNA1_outfile"


## GUIDE 2
# set gRNA2 file name variables
# pretty sure the FASTQ has to be unzipped
gRNA2_idx_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_indices/pgRNA2"
gRNA2_infile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/demultiplex/200122_PC9_screen/PP_pgRNA_PC9_S1_I1_001_trimmed.fastq_sample5.fastq"
gRNA2_outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200122_PC9_screen/PP_pgRNA_PC9_S1_I1_001_trimmed_sample5_aligned.sam"

bowtie -q -v 1 --best --strata --all --sam -p 4 \
  "$gRNA2_idx_file" \
  "$gRNA2_infile" \
  "$gRNA2_outfile"
