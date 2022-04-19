## Phoebe C. R. Parrish
## Berger Lab, FHCRC
## 2022-03-01
## updated 2022-03-07
## snakemake v7.1.0

import os
import sys
import re
import pandas as pd

configfile: "config.yaml"

## set up paths for output files
## do this automatically for any line containing "dir"?
fastq_dir = os.path.abspath(config["fastq_dir"]) + "/"
fastq_trimmed_dir = os.path.abspath(config["fastq_trimmed_dir"]) + "/"
demultiplexed_dir = os.path.abspath(config["demultiplexed_dir"]) + "/"
base_filename = config["base_filename"]
bowtie_idx_dir = os.path.abspath(config["bowtie_idx_dir"]) + "/"
bowtie_idx_ref_dir = os.path.abspath(config["bowtie_idx_ref_dir"]) + "/"
aligned_dir = os.path.abspath(config["aligned_dir"]) + "/"
pgRNA_counts_dir = os.path.abspath(config["pgRNA_counts_dir"]) + "/"
counter_script_dir = os.path.abspath(config["counter_script_dir"]) + "/"


read_to_fastq = {}
fastq_trim_coords = {}

## note: find a way to not hard-code this (in case seq strategy changes)
## convert to a dictionary too (so position isn't key)
trim_coords = [[1, 20], [2, 21], [1, 6]]

fastq_fofn = open(config["fastq_fofn"], "r")
i = 0
for fastq in fastq_fofn:
    fastq = fastq.strip()
    print(fastq)
    print(trim_coords[i])
    key = "R" + str(i+1)
    print("key =", key)
    read_to_fastq[key] = fastq
    fastq_trim_coords[key] = trim_coords[i]

    i += 1
fastq_fofn.close()

print("read_to_fastq=", read_to_fastq)
print("fastq_trim_coords=", fastq_trim_coords)

sample_list = []
barcode_ref_file = open(demultiplexed_dir + "ref/" + config["barcode_ref_file"], "r")
for line in barcode_ref_file:
    ## get sample names from idemp ref file
    sample = line.strip().split("\t")[1]
    # print("sample =", sample)
    sample_list.append(sample)
barcode_ref_file.close()
print("sample_list = ", sample_list)

# print(sample_to_fastq)


# sample_wildcard_to_fastq

# trim_in = config["fastq_dir"] + "PP_pgRNA_HeLa_S1_R3_001.fastq.gz"
# trim_out = config["fastq_trimmed_dir"] + "PP_pgRNA_HeLa_S1_R3_001_trimmed.fastq"

## this is related to the wildcard thing - output file names for each rule
rule all:
    input:
        expand(fastq_trimmed_dir + base_filename + "_{read}_trimmed.fastq",
            read = read_to_fastq.keys()),
        # demultiplexed_dir,
        expand(demultiplexed_dir + base_filename + "_{dmx_read}_trimmed.fastq_{sample}.fastq",
            dmx_read = ["R1", "R2"],
            sample = sample_list),
        expand(demultiplexed_dir + base_filename + "_{dmx_read}_trimmed_{sample}.fastq",
            dmx_read = ["R1", "R2"],
            sample = sample_list),
        expand(bowtie_idx_dir + "pgPEN_{dmx_read}",
            dmx_read = ["R1", "R2"]),
        expand(aligned_dir + "sam/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned.sam",
            dmx_read = ["R1", "R2"],
            sample = sample_list),
        expand(aligned_dir + "bam_sorted/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned_sorted.bam",
            dmx_read = ["R1", "R2"],
            sample = sample_list),
        expand(aligned_dir + "flagstat/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned_sorted_flagstat.txt",
            dmx_read = ["R1", "R2"],
            sample = sample_list)


## add wildcards!
rule trim_reads:
    input:
        # fastq_dir + "PP_pgRNA_HeLa_S1_{read}_001.fastq.gz"
        lambda wildcards: fastq_dir + read_to_fastq[wildcards.read]
    output:
        # config["fastq_trimmed_dir"] + "PP_pgRNA_HeLa_S1_R3_001_trimmed.fastq"
        fastq_trimmed_dir + base_filename + "_{read}_trimmed.fastq"
    params: ## is this the best way to do this?
        start = lambda wildcards: fastq_trim_coords[wildcards.read][0],
        end = lambda wildcards: fastq_trim_coords[wildcards.read][1]
    log:
        "logs/trim_fastqs/{read}.log"
    shell:
        "zcat {input} | fastx_trimmer -f {params.start} -l {params.end} -o {output}"


## note to self: is this the best way to install/run idemp? (reproducible?)
rule demux_fastqs:
    input:
        idx = fastq_trimmed_dir + base_filename + "_R3_trimmed.fastq",
        R1 = fastq_trimmed_dir + base_filename + "_R1_trimmed.fastq",
        R2 = fastq_trimmed_dir + base_filename + "_R2_trimmed.fastq",
        ref = demultiplexed_dir + "ref/" + config["barcode_ref_file"]
    output:
        ## put the expand statement here so this is only run once
        ## (instead of once per output file)
        expand(demultiplexed_dir + base_filename + "_{dmx_read}_trimmed.fastq_{sample}.fastq",
            dmx_read = ["R1", "R2"],
            sample = sample_list)
    params:
        idemp = config["idemp"],
        n_mismatch = 1,
        ## this is a param not output because having a folder as the output messes everything up
        out_dir = demultiplexed_dir
    log:
        "logs/demux_fastqs/idemp.log"
    shell:
        """
        {params.idemp} -b {input.ref} -n {params.n_mismatch} \
        -I1 {input.idx} -R1 {input.R1} -R2 {input.R2} -o {params.out_dir}
        """


## TO DO: resolve conflict b/w R1/R2 and gRNA1/gRNA2 naming schemes...have to fix this
##    somehow so I can rename the idemp output fastqs


## remove the extra ".fastq" that idemp puts in for some reason
rule rename_fastqs:
    input:
        demultiplexed_dir + base_filename + "_{dmx_read}_trimmed.fastq_{sample}.fastq"
    output:
        demultiplexed_dir + base_filename + "_{dmx_read}_trimmed_{sample}.fastq"
    shell:
        "mv {input} {output}"

    # run: # python script
    #     for filename in os.scandir(input):
    #         if filename.is_file():
    #             if "decode" in filename:
    #                 continue
    #             else:
    #                 print("filename =", filename)
    #                 rm_start = filename.index(".fastq")
    #                 rm_end = rm_start + 6
    #                 new_prefix = filename[:rm_start] ## get filename up to the start of .fastq
    #                 new_suffix = filename[rm_end:]
    #                 new_filename = new_prefix + new_suffix
    #                 print("new_filename =", new_filename)
    #                 rename_fastq[new_filename] = filename
    #                 os.rename(new_filename, filename)
    #     print("rename_fastq =", rename_fastq)

rule build_bowtie_index:
    input:
        bowtie_idx_ref_dir + "pgPEN_{dmx_read}.fa"
    output:
        bowtie_idx_dir + "pgPEN_{dmx_read}"
    log:
        "logs/align_reads/{dmx_read}"
    shell:
        "bowtie-build -f {input} {output}"

rule align_reads:
    input:
        fastq = demultiplexed_dir + base_filename + "_{dmx_read}_trimmed_{sample}.fastq",
        idx = bowtie_idx_dir + "pgPEN_{dmx_read}"
    output:
        aligned_dir + "sam/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned.sam"
    resources:
        mem = 8,
        time = 24
    log:
        "logs/align_reads/{dmx_read}_{sample}.log"
    shell:
        "bowtie -q -v 1 --best --strata --all --sam -p 4 {input.idx} {input.fastq} {output}"

rule make_sorted_bam:
    input:
        aligned_dir + "sam/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned.sam"
    output:
        aligned_dir + "bam_sorted/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned_sorted.bam"
    params:
        unsorted_bam = aligned_dir + "bam/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned.bam",
        tmp = aligned_dir + "bam_sorted/tmp_{dmx_read}_{sample}"
    resources:
        mem = 16,
        time = 24
    log:
        "logs/make_sorted_bam/{dmx_read}_{sample}.log"
    shell:
        """
        samtools view -bS -o {params.unsorted_bam} {input}

        samtools sort -O bam -n {params.unsorted_bam} -T {params.tmp} -o {output}
        """

rule get_stats:
    input:
        aligned_dir + "bam_sorted/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned_sorted.bam"
    output:
        aligned_dir + "flagstat/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned_sorted_flagstat.txt"
    shell:
        "samtools flagstat {input} > {output}"

## note: make a new rule c to combine all sample counts => 1 file
## update snakemake envt to include R (& other) dependencies
# run count_pgRNAs:
#     input:
#         ## fix this bc the filename needs to be diff for counter_efficient.R ... or fix that script?
#         expand(aligned_dir + "bam_sorted/" + base_filename + "_{dmx_read}_trimmed_{sample}_aligned_sorted.bam",
#             dmx_read = ["R1"],
#             sample = sample_list[0])
#     output:
#         pgRNA_counts_dir +
#     params:
#         n_chunks =
#         script = counter_script_dir +
#     resources:
#         mem = 98
#         time = 8
#         cpus-per-task = 16
#     log:
#     shell:
