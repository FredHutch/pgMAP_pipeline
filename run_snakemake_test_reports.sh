#!/bin/bash

## run this if the Snakemake environment has not been set up yet
# mamba create -f $CONDA_ENV

## activate conda envt
source activate snakemake
echo "snakemake env activated"

echo -e "Making html report..."
snakemake --report "workflow/report/report.html"

echo -e "Exporting pipeline rulegraph to svg and pdf..."
snakemake --snakefile "workflow/Snakefile" \
  --rulegraph > "workflow/report/rulegraph.dot"
dot -Tpdf "workflow/report/rulegraph.dot" > "workflow/report/pipeline_rulegraph.pdf"
dot -Tsvg "workflow/report/rulegraph.dot" > "workflow/report/pipeline_rulegraph.svg"
rm "workflow/report/rulegraph.dot"

## export PDF and svg visualizations of DAG structure of pipeline steps
## NOTE: this will not work if there are print statements in the pipeline
echo -e "Exporting pipeline DAG to svg and pdf..."
snakemake --snakefile "workflow/Snakefile" \
  --dag > "workflow/report/dag.dot"
dot -Tpdf "workflow/report/dag.dot" > "workflow/report/pipeline_dag.pdf"
dot -Tsvg "workflow/report/dag.dot" > "workflow/report/pipeline_dag.svg"
rm "workflow/report/dag.dot"

echo -e "\nDONE!\n"
