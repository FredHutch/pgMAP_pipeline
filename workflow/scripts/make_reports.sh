#!/bin/bash

echo -e "Exporting pipeline rulegraph to svg and pdf..."
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --rulegraph > "workflow/report/rulegraph.dot"
dot -Tpdf "workflow/report/rulegraph.dot" > "workflow/report/pipeline_rulegraph.pdf"
dot -Tsvg "workflow/report/rulegraph.dot" > "workflow/report/pipeline_rulegraph.svg"
rm "workflow/report/rulegraph.dot"

## export PDF and svg visualizations of DAG structure of pipeline steps
## NOTE: this will not work if there are print statements in the pipeline
echo -e "Exporting pipeline DAG to svg and pdf..."
snakemake --snakefile workflow/Snakefile \
  --configfile config/config.yaml \
  --dag > "workflow/report/dag.dot"
dot -Tpdf "workflow/report/dag.dot" > "workflow/report/pipeline_dag.pdf"
dot -Tsvg "workflow/report/dag.dot" > "workflow/report/pipeline_dag.svg"
rm "workflow/report/dag.dot"

echo -e "\nReports created! See workflow/report/ folder for output.\n"
