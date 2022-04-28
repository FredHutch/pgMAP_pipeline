#!/bin/bash

# run idemp - remember it is only installed on my account so this must be run locally and not submitted via sbatch
# can use grabnode with 28 cores to get max power for running this

# assign variables
barcode_ref_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/demultiplex/200122_PC9_screen/screen_barcodes.txt"
index_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/200122_PC9_screen/PP_pgRNA_PC9_S1_I2_001_trimmed.fastq"
gRNA1_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/200122_PC9_screen/PP_pgRNA_PC9_S1_R1_001_trimmed.fastq"
gRNA2_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/200122_PC9_screen/PP_pgRNA_PC9_S1_I1_001_trimmed.fastq"
outfile_path="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/demultiplex/200122_PC9_screen/"

# run command
## -n 1 allows 1 mismatch in the barcodes
idemp -b "$barcode_ref_file" -n 1 -I1 "$index_file" -R1 "$gRNA1_file" -R2 "$gRNA2_file" -o "$outfile_path"
