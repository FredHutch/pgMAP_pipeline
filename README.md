# pgPEN_pipeline

## Running the pipeline

1. Clone (or fork) the git repo from https://github.com/bergerbio/pgPEN_pipeline. Use the main branch for now.

2. Copy the following config files and edit them to point to your analysis files:
    * `barcode_ref_file.sample.txt` => `barcode_ref_file.txt`
    * `config.sample.yaml` => `config.yaml`
    * `fastqs.sample.fofn` => `fastqs.fofn`

3. **First time only:** make a Snakemake conda environment (defined by `workflow/envs/snakemake.yaml`) using the following steps:
  1. Install the Mamba package manger as described in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) by either:
    * If you do not already have Conda installed, [install Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
    * If you already have Conda installed, run the following command: `$ conda install -n base -c conda-forge mamba`
  2. Then create your Snakemake environment using Mamba by either:
    * Running the command `mamba env create -f workflow/envs/snakemake.yaml`
    * Un-commenting line 16 in the script `run_snakemake.sh` (command: `# mamba create -f $CONDA_ENV`)

4. Run the script `run_snakemake.sh` using the command: `bash run_snakemake.sh`

### Snakemake installation documentation
Install Snakemake v7.1.0 using mambaforge as described here:
https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

Folder setup/running info as described here:
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html


### Phoebe note to self:
Use ctrl + shift + M to see rendered markdown in Atom

# To Do

* increase memory/cores for demux_fastqs rule

* figure out what is going on with conda error about placeholder of length '80' - can I just use mamba to fix the problem? Or do I need to use the --conda-prefix option? (Test without --conda-prefix set as home to see)

* get package versions from conda environment => add to envt.yamls - done for all except Python env

* figure out how to make the fastq.fofn within the pipeline

* add step to make directories (for log files and output files) if they don't already exist

* figure out how to access snakemake variables from inside python scripts

* add fastQC step (once I do this I will have to add the output to rule all)

* add a rule to gzip all files at the end of the pgPEN pipeline

* Decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* replace "." with "_" in variable names for counter_efficient.R

* move idemp install from my home directory => berger_a fast directory
  * paralog_pgRNA/pgPEN_library/installs

* figure out how to copy files even if I am gitignoring them (so they can be shared/modified by others)

* decide about improving config file setup/tracking based on: https://gist.github.com/canton7/1423106

* add step to create log folders at the beginning of the pipeline? Or track them from git?

* save root dir as yaml variable then concat within python

* copy test_snakemake fastqs to folder and change directory for that command

* add idemp to conda? or find some other way to ensure that users can easily install idemp and add it to their bashrc

* make interpretable error messages for n_chunks, etc.

* link git and Atom

* figure out a clearer way to explain conda/mamba installs - both Saksham & Emma had issues with this

* consider directly calling the R2 sorted BAM from the count_pgRNAs rule (so that it shows up on the DAG and stuff)

* check that log files are actually being made correctly
  * see: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files

* figure out what's going on with build_bowtie_index - why is it running so many times? maybe change the Snakemake output so it doesn't exactly match the rule output?

* figure out how to run Snakemake correctly on the cluster
  * https://snakemake.readthedocs.io/en/stable/executing/cluster.html
  * https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads

* figure out how to run atom remotely:
  * https://gist.github.com/NTag/9d9be611e03098c282241652894bda7f
  * https://www.quora.com/What-is-the-best-way-to-use-Atom-io-over-SSH
