#!/bin/bash

## configure file paths
CONFIG_FILE="config.yaml"
SNAKE_FILE="Snakefile"
CONDA_ENV="environment.yaml"

## activate conda envt
source activate snakemake_env

## run the pipeline
snakemake --snakefile $SNAKE_FILE \
  --configfile $CONFIG_FILE \
  --use-conda --conda-prefix $CONDA_ENV \
  --cores --restart-times 3

## export PDF and svg visualizations of DAG structure of pipeline steps
## NOTE: this will not work if there are print statements in the pipeline
echo -e "Exporting pipeline DAG to svg and pdf..."
snakemake --snakefile $SNAKE_FILE \
  --configfile $CONFIG_FILE \
  --dag > snakemake_out/dag.dot
dot -Tpdf dag.dot > snakemake_out/pipeline_dag.pdf
dot -Tsvg dag.dot > snakemake_out/pipeline_dag.svg
rm dag.dot

echo -e "Generating pipeline HTML report..."
snakemake --configfile $CONFIG_FILE \
  --snakefile $SNAKE_FILE \
  --report snakemake_out/report.html

echo -e "\nDONE!\n"
