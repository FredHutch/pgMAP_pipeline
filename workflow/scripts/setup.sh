#!/bin/bash

## check if idemp successfully installed
IDEMP=config/idemp/idemp
if test -f "$IDEMP"; then
    echo "idemp has been successfully installed and compiled in config folder."
else 
    echo "$IDEMP does not exist. Installing idemp."
    ## install idemp in config dir
    git -C config clone --quiet https://github.com/yhwu/idemp.git
    echo "Compiling idemp."
    ## compile idemp
    make -C config/idemp  
fi

## make output dirs
mkdir -p "results/fastq_trimmed"
mkdir -p "results/fastq_demuxed"
mkdir -p "results/sam"
mkdir -p "results/bam"
mkdir -p "results/bam_sorted"
mkdir -p "results/pgRNA_counts"

## make log dirs
mkdir -p "workflow/logs/trim_reads"
mkdir -p "workflow/logs/demux_fastqs"
mkdir -p "workflow/logs/align_reads"
mkdir -p "workflow/logs/make_sorted_bams"
mkdir -p "workflow/logs/get_stats"
mkdir -p "workflow/logs/count_pgRNAs"
mkdir -p "workflow/logs/combine_counts"