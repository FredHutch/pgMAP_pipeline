#!/bin/bash

# load module
module load samtools

# assign variables
sam_file="/fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/Undetermined.R1_trimmed_aligned.sam"
bam_file="/fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/Undetermined.R1_trimmed_aligned.bam"
sorted_bam="/fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/Undetermined.R1_trimmed_aligned_sorted.bam"

### OPTIONS
## view
# -bS: outputs to BAM format; S is no longer used but specifies SAM input format
## sort
# -n: sorts by read names rather than chromosomal coordinates

samtools view -bS -o "$bam_file" "$sam_file"

samtools sort -O bam -n "$bam_file" -T tmp1 -o "$sorted_bam"
