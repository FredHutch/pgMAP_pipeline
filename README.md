# pgPEN_pipeline

Install the latest Snakemake version (7.1?) using mambaforge as described here:
https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

Folder setup/running info as described here:
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html

Use ctrl + shift + M to see rendered version in Atom

# To Do

* run test analysis on the cluster instead of on an interactive node

* run on Arnab's data - write to a shared folder that he can also access?

* reorganize folders on my fast drive - by screen and then sub-directory

* figure out rule all - what do I actually need to include here? where should I define wildcards if not here?

* remove weird conda files from Snakemake dir

* Decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* replace "." with "_" in variable names for counter_efficient.R

* add fastQC step?

* read Vince Buffalo's snakemake guide

* might have to make log folder first
