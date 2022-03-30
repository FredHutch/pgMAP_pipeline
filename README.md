# pgPEN_pipeline

Install the latest Snakemake version (7.1?) using mambaforge as described here:
https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

# To Do
* Decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* get rid of R1 and R2 dicts, convert to list of just sample names
  * expand list and use ad hoc list of R1, R2 to get new file names
* use read_to_fastq to get fastq filenames (lambda function again)

  * function sample_to_wildcard if file name doesn't match sample name
  * might have to make log folder first

  * change pooled => pool in filenames
