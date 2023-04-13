# pgMAP_pipeline

## Running the pipeline

1. Clone or fork the git repo from https://github.com/FredHutch/pgMAP_pipeline/ (if you are not sure of the difference between cloning and forking, check out the explainer [here](https://github.com/FredHutch/pgMAP_pipeline/)). Make sure you are on the main branch. 


2. Update the config files as described below: Duplicate the following files and edit them so that the barcodes and samples match those used in your sequencing run and all file paths point to your working directory:
    * Make a copy of `barcode_ref_file.sample.txt` named `barcode_ref_file.txt`. Update the sample and barcode information to match your experimental design and sequencing setup. 
    * Make a copy of `config.sample.yaml` named `config.yaml`. Update the `base_filename` variable, the read coordinates, and the number of chunks to split your BAM files into. 


3. **First time only:** Build a Conda environment for Snakemake (defined by `workflow/envs/snakemake.yaml`, for more detail see the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)). Running your analysis within `snakemake_env` will enable you to use the same version of Snakemake and all other supporting packages as we used when we developed pgMAP. To build the `snakemake_env`, do one of the following: 
   1. Install [Mamba](https://mamba.readthedocs.io/en/latest/), which is essentially a faster version of Conda that is required to run Snakemake, by doing one of the following:
      * If you **do not** already have the Conda package manager installed, [install Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
      * If you **do** already have Conda installed, run the following command: `$ conda install -n base -c conda-forge mamba`
   2. Next, create your Snakemake environment using Mamba by either:
      * Running the command `mamba env create -f workflow/envs/snakemake.yaml`
      * Un-commenting line 16 in the script `run_snakemake.sh` (command: `# mamba create -f $CONDA_ENV`)


4. Run the script `run_snakemake.sh` by entering the command: `bash run_snakemake.sh`


### More info on installing/running/troubleshooting Snakemake
Install Snakemake v7.1.0 using mambaforge as described here:
https://snakemake.readthedocs.io/en/stable/tutorial/setup.html

Folder setup/running info as described here:
https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html


### Phoebe notes to self:
* running conditional rules: https://stackoverflow.com/questions/64949149/is-it-possible-to-add-a-conditional-statement-in-snakemakes-rule-all

* see Mitchell's slideshow on creating smk workflows: https://mrvollger.github.io/SmkTemplate/slides/#2

# To Do

## Necessary before submission

* **Phoebe:** figure out why the pipeline is only running 3 rules??? look at old code to compare what might be missing

* **Daniel:** please confirm the FASTQC rule is running/producing the correct output

* **Daniel:** make sure all of your old/inactive branches have been merged and deleted

* **Daniel:** make a git tracked results/fastq/ dir with a README file saying "put your FASTQs here"

* **Daniel:** make a git 

* **Phoebe:** make sure all of your old/inactive branches have been merged and deleted 

* run full pipeline on test dataset using an interactive node

* run pipeline on test dataset using cluster submit

* **Phoebe:** run the pipeline on Arnab's pgPEN data

* **Phoebe:** make sure the test workflow can be run on the cluster with sbatch
  * https://snakemake.readthedocs.io/en/stable/executing/cluster.html
  * https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#threads

* **Phoebe:** get bowtie index rule to work!
  * LATER: figure out what's actually going on and get the code to work.
  * write some kind of lambda statement that will make the bowtie index rule work?
  * figure out what's going on with build_bowtie_index - why is it running so many times? maybe change the Snakemake output so it doesn't exactly match the rule output?

* update step #2 of **Running the Pipeline** to make it clear what actually needs to be updated (after incorporating Daniel's changes/getting the pipeline to run on PC9 data and Arnab's data)

* make sure the report.html is actually being ignored by .gitignore

* **Phoebe:** get VScode Snakemake extension?

* get rid of "dir" in results_dir_dict key names

* make a new branch and make sure the existing version of the pipeline runs on the test_snakemake downsampled dataset 
  * run on an interactive node first
  * once the scheduler part is working, confirm that this works for the downsampled dataset too

* update config files, etc. for full pgPEN library

* add other stats output apart from counts?

* increase memory/cores for demux_fastqs rule

* add step to make directories (for log files and output files) if they don't already exist - or just add those directories to git tracking but ignore their contents?

* decide on a consistent file naming strategy
  * name based on pgRNA/full sample name or just numbers?

* figure out how to install and successfully idemp (does it have to be added to a user's bashrc?) as part of the pipeline/in a conda environment, or find another demultiplexer that users can download through conda
  * confirm that you get the exact same output from new demultipexer as you do from idemp
  * remove idemp from config.sample.yaml and other config.yaml files

* check that log files are actually being made correctly (esp. when being run via sbatch)
  * see: https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#log-files

* run the pipeline on the full PC9 dataset and confirm that we get the same output as what was published in the pgPEN paper

* **Phoebe:** make sure the pipeline can (conceptually) be run on any dual gRNA sequencing dataset

* run the pipeline on data from another dual gRNA sequencing approach 

* **Daniel/Alice:** confirm that updated conda/mamba/Snakemake install instructions make sense

* run the pipeline on a Bradley lab dataset from the new sequencing strategy
  * Taylor has 2 pools that use the same barcodes - **Phoebe/Emma** figure out what to do about this
  
* **Daniel:** write a vignette for running pgMAP on the downsampled dataset on an interactive node (can just add to GitHub README)
  * example: https://github.com/Atkinson-Lab/Tractor-tutorial


### Completed
* merge Daniel's changes to the Snakefile into Phoebe's dev branch

* **Phoebe:** get bowtie index rule to work!
  * for now, just copy the bowtie index to the results folder. Pull Daniel's changes. 

* figure out why Daniel's changes aren't showing up in the Snakefile

* figure out what is going on with conda error about placeholder of length '80' - can I just use mamba to fix the problem? Or do I need to use the --conda-prefix option? (Test without --conda-prefix set as home to see)

* config.yaml: save root dir as yaml variable then concat within python
  * future note: make the whole pipeline more self-contained - tell users to run the pipeline from within the root dir and only use relative file paths from there?

* figure out how to make folders in a way that makes sense

* **Daniel:** figure out how to make the fastq.fofn within the pipeline (look at the CRISPR_pipeline dev Snakefile for reference) (**CONFIRMED**)

* **Daniel:** counter.R: replace "." with "_" in variable names for counter_efficient.R - for dataframes, please make sure the variable names is still "d.[rest_of_var_name]"! (**CONFIRMED**)

* **Daniel:** add fastQC step (once you do this you will have to add the output to rule all) (**CONFIRMED**)

* **Phoebe:** delete current dev branch and make a new one for each update or each group of updates

* get Python env package versions from conda environment => add to envt.yaml

* ask Daniel about fastqc wrapper - why do all the output files say ".fastq.gz_fastqc.zip"?
  * find out what the wrapper is - I think this is messing w/ the other rules running
  
* copy example FASTQ files from Bradley lab folder => Berger lab folder
  

## Can probably be done after submission

* **Phoebe:** make interpretable error messages for n_chunks, etc.

* add a rule to gzip all large files at the end of the pgPEN pipeline

* add in a "test case" option for people to try out the aligner on the downsampled PC9 dataset

* decide about improving config file setup/tracking based on: https://gist.github.com/canton7/1423106

* make a nicely-formatted report for counts output, alignment stats

* move the results/fastq + README to an input directory


### Completed



## Would be nice but probably unnecessary

* figure out how to access snakemake variables from inside python scripts

### Completed


