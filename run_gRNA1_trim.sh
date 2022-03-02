#!/bin/bash

# load module
module load FASTX-Toolkit/0.0.14-foss-2018b

# set variables
## infile should be the standard Illumina R1 (pgRNA1)
infile="/home/pparrish/ngs/illumina/200722_D00300_1006_AHFFK7BCX3/Unaligned/Project_pparrish/PP_pgRNA_HeLa_S1_R1_001.fastq.gz"
outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/200722_HeLa_screen/PP_pgRNA_HeLa_S1_R1_001_trimmed.fastq"

# command to run
## recommended to call variables inside double quotes to prevent interpretation of special symbols
zcat "$infile" | fastx_trimmer -f 1 -l 20 -o "$outfile"
