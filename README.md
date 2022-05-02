# pgPEN_pipeline

Install the latest Snakemake version (7.1?) using mambaforge as described here:
https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

Folder setup/running info as described here:
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html

Use ctrl + shift + M to see rendered version in Atom

# To Do

* figure out how to fix moving demuxed fastqs => new name causing it to rerun the whole pipeline because the output from demux_fastqs rule is not there

* fix whatever is going on with counter_efficient.R to produce incorrect output
  * confirm I'm not going to have the same issue as above

* remove weird conda files from Snakemake dir

* Decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* might have to make log folder first

* add fastQC step?
