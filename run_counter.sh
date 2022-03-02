#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=98G

# load appropriate R module with necessary packages
module load R/3.6.2-foss-2018b-fh1

# run script
Rscript /home/pparrish/paralog_pgRNA_screen/scripts/counter.R
