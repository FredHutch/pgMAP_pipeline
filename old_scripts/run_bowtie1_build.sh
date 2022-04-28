#!/bin/bash

module load bowtie

bowtie-build -f /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgRNA_fastas/paralog_pgRNA1.fa \
  /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_indices/pgRNA1

bowtie-build -f /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/pgRNA_fastas/paralog_pgRNA2.fa \
  /fh/fast/berger_a/grp/bergerlab_shared/Projects/paralog_pgRNA/bowtie1_indices/pgRNA2
