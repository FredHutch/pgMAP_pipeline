#!/bin/bash

# load module
module load FASTX-Toolkit/0.0.14-foss-2018b

# set variables
## infile should be the standard Illumina I1 (pgRNA1)
infile="/home/pparrish/ngs/illumina/200722_D00300_1006_AHFFK7BCX3/Unaligned/Project_pparrish/PP_pgRNA_HeLa_S1_I1_001.fastq.gz"
outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/200722_HeLa_screen/PP_pgRNA_HeLa_S1_I1_001_trimmed.fastq"

# command to run
zcat "$infile" | fastx_trimmer -f 2 -l 21 -o "$outfile"
