# pgMAP_pipeline

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


### Phoebe notes to self:
* Use ctrl + shift + M to see rendered markdown in Atom

* running conditional rules: https://stackoverflow.com/questions/64949149/is-it-possible-to-add-a-conditional-statement-in-snakemakes-rule-all


# To Do

## Necessary before submission

* **Phoebe:** make sure the workflow can be run on the cluster with sbatch
  * https://snakemake.readthedocs.io/en/stable/executing/cluster.html
  * https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads

* **Phoebe:** run the pipeline on Arnab's pgPEN data

* **Phoebe:** get bowtie index rule to work!
  * for now, just copy the bowtie index to the results folder. Pull Daniel's changes. 
  * LATER: figure out what's actually going on and get the code to work.
  * write some kind of lambda statement that will make the bowtie index rule work?
  * figure out what's going on with build_bowtie_index - why is it running so many times? maybe change the Snakemake output so it doesn't exactly match the rule output?

* **Phoebe:** make sure the report.html is actually being ignored by .gitignore

* **Phoebe:** delete current dev branch and make a new one for each update or each group of updates

* **Daniel:** make a new branch and make sure the existing version of the pipeline runs on the test_snakemake downsampled dataset 
  * run on an interactive node first
  * once Phoebe has gotten the scheduler part working, confirm that this works for the downsampled dataset too

* **Daniel/Emma:** confirm that updated conda/mamba/Snakemake install instructions make sense

* update config files, etc. for full pgPEN library

* get Python env package versions from conda environment => add to envt.yaml

* add other stats output apart from counts?

* make a git tracked results/fastq/ dir with a TXT file saying "put your FASTQs here"?

* increase memory/cores for demux_fastqs rule

* add step to make directories (for log files and output files) if they don't already exist - or just add those directories to git tracking but ignore their contents?

* Decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* figure out how to install and successfully idemp (does it have to be added to a user's bashrc?) as part of the pipeline/in a conda environment, or find another demultiplexer that users can download through conda
  * confirm that you get the exact same output from new demultipexer as you do from idemp

* check that log files are actually being made correctly (esp. when being run via sbatch)
  * see: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files


* **Daniel:** run the pipeline on the full PC9 dataset and confirm that we get the same output as what was published in the pgPEN paper

* **Phoebe:** make sure the pipeline can (conceptually) be run on any dual gRNA sequencing dataset

* run the pipeline on data from another dual gRNA sequencing approach 

### Completed - Daniel please confirm
* **Daniel:** figure out how to make the fastq.fofn within the pipeline (look at the CRISPR_pipeline dev Snakefile for reference) (**CONFIRMED**)

* **Daniel:** counter.R: replace "." with "_" in variable names for counter_efficient.R - for dataframes, please make sure the variable names is still "d.[rest_of_var_name]"! (**CONFIRMED**)

* **Daniel:** add fastQC step (once you do this you will have to add the output to rule all) (**CONFIRMED**)


### Completed

* merge Daniel's changes to the Snakefile into Phoebe's dev branch

* figure out why Daniel's changes aren't showing up in the Snakefile

* figure out what is going on with conda error about placeholder of length '80' - can I just use mamba to fix the problem? Or do I need to use the --conda-prefix option? (Test without --conda-prefix set as home to see)

* config.yaml: save root dir as yaml variable then concat within python
  * future note: make the whole pipeline more self-contained - tell users to run the pipeline from within the root dir and only use relative file paths from there?

* figure out how to make folders in a way that makes sense



## Can probably be done after submission

* **Phoebe:** make interpretable error messages for n_chunks, etc.

* add a rule to gzip all large files at the end of the pgPEN pipeline

* add in a "test case" option for people to try out the aligner on the downsampled PC9 dataset

* decide about improving config file setup/tracking based on: https://gist.github.com/canton7/1423106

* make a nicely-formatted report for counts output, alignment stats

### Completed



## Would be nice but probably unnecessary

* figure out how to access snakemake variables from inside python scripts

### Completed


