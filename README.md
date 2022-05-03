# pgPEN_pipeline

Install the latest Snakemake version (7.1?) using mambaforge as described here:
https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

Folder setup/running info as described here:
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html

Use ctrl + shift + M to see rendered version in Atom

# To Do

* fix whatever is going on with counter_efficient.R to produce incorrect output
  * confirm I'm not going to have the same issue as above

* read Vince Buffalo's snakemake guide

* figure out rule all - what do I actually need to include here? where should I define wildcards if not here?

* remove weird conda files from Snakemake dir

* Decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* might have to make log folder first

* add fastQC step?

* replace "." with "_" in variable names for counter_efficient.R
