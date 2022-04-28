#!/bin/bash

module load samtools

# assign variables
sam_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200122_PC9_screen/PP_pgRNA_PC9_S1_I1_001_trimmed_aligned.sam"
bam_file="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200122_PC9_screen/PP_pgRNA_PC9_S1_I1_001_trimmed_aligned.bam"
sorted_bam="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_aligned/200122_PC9_screen/PP_pgRNA_PC9_S1_I1_001_trimmed_aligned_sorted.bam"

### OPTIONS
## view
# -bS: outputs to BAM format; S is no longer used but specifies SAM input format
## sort
# -n: sorts by read names rather than chromosomal coordinates

samtools view -bS -o "$bam_file" "$sam_file"

samtools sort -O bam -n "$bam_file" -T tmp2 -o "$sorted_bam"
