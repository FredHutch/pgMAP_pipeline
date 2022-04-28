#!/bin/bash

## SETUP
# set job parameters
#SBATCH --job-name=trim_reads
#SBATCH --output=trim_reads_%A_%a.out
#SBATCH --array=1-3
#SBATCH --time=3:00:00
#SBATCH --ntasks=1
#SBATCH --mem=8G

# load module
module load FASTX-Toolkit/0.0.14-foss-2018b

# set variables
## R1 = R1, R2 = I1, R3 = I2
infile="/home/pparrish/ngs/illumina/200722_D00300_1006_AHFFK7BCX3/Unaligned/Project_pparrish/test_array/PP_pgRNA_HeLa_S1_R${SLURM_ARRAY_TASK_ID}_001.fastq.gz"
outfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/trimmed_fastqs/200722_HeLa_screen/test_array/PP_pgRNA_HeLa_S1_R${SLURM_ARRAY_TASK_ID}_001_trimmed.fastq"

## make an associative array because bash arrays are zero-indexed so using a normal array
##   was giving me a headache (kept having to subtract 1 b/c the job array is one-indexed)
declare -A start_pos=([1]=1 [2]=2 [3]=1)
declare -A end_pos=([1]=20 [2]=21 [3]=6)

## test code to just print start/end coords
# echo "start pos = " ${start_pos[${SLURM_ARRAY_TASK_ID}]}
# echo "end pos = " ${end_pos[${SLURM_ARRAY_TASK_ID}]}

# command to run
## recommended to call variables inside double quotes to prevent interpretation of special symbols
zcat "$infile" | fastx_trimmer -f ${start_pos[${SLURM_ARRAY_TASK_ID}]} \
  -l ${end_pos[${SLURM_ARRAY_TASK_ID}]} \
  -o "$outfile"
