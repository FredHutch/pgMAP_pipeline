#!/bin/bash

# load module
module load FASTX-Toolkit/0.0.14-foss-2018b

# set variables
## infile should be the standard Illumina I2 (index)
infile="/home/pparrish/ngs/illumina/200722_D00300_1006_AHFFK7BCX3/Unaligned/Project_pparrish/PP_pgRNA_HeLa_S1_I2_001.fastq.gz"
outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/200722_HeLa_screen/PP_pgRNA_HeLa_S1_I2_001_trimmed.fastq"

# command to run
## recommended to call variables inside double quotes to prevent interpretation of special symbols
zcat "$infile" | fastx_trimmer -f 1 -l 6 -o "$outfile"








#module load FASTX-Toolkit/0.0.14-foss-2016a

# then run
#zcat /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/fastqs/Undetermined.R3.fastq.gz | \
#  fastx_trimmer -f 1 -l 6 -o /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/Undetermined.R3_trimmed.fastq
