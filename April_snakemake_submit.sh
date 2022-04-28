#!/bin/bash

mkdir -p logs/picard
mkdir -p logs/indexbam
mkdir -p logs/rMATS
mkdir -p logs/R

snakemake --jobs 50 \
    --cluster " sbatch -n {threads} -o {log} --mem=8096 " \
    --use-conda -k -p \
    --latency-wait 180
