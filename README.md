# pgPEN_pipeline

Install the latest Snakemake version (7.1?) using mambaforge as described here:
https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

Folder setup/running info as described here:
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html

Use ctrl + shift + M to see rendered version in Atom

# To Do

* increase memory/cores for demux_fastqs rule

* get package versions from conda environment => add to envt.yamls

* figure out how to copy stuff in tmux!!!

* reorganize folders on my fast drive - by screen and then sub-directory

* figure out how to make the fastq.fofn within the pipeline

* add step to make directories (for log files and output files) if they don't already exist

* figure out how to access snakemake variables from inside python scripts

* add fastQC step

* add a rule to gzip all files at the end of the pgPEN pipeline

* figure out rule all - what do I actually need to include here? where should I define wildcards if not here?

* remove weird conda files from Snakemake dir

* Decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* replace "." with "_" in variable names for counter_efficient.R

* might have to make log folder first
