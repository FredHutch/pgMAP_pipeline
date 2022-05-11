import os
import sys
import re
import pandas as pd

def main(annot_file, in_dir, out_file):
    ## read in pgRNA annotations to get IDs
    annot_df = pd.read_csv(annot_file, sep = "\t") ## will this work??

    ## loop through files in pgRNA_counts_dir
    dfs = [annot_df]
    
    ## sorted/lambda stuff is to order files alphabetically
    for file in sorted(os.scandir(in_dir), key=lambda e: e.name):
        df = pd.read_csv(file, sep = "\t")
        counts = df.loc[:, df.columns.str.contains("counts")] ## extract counts col
        dfs.append(counts)

    out_df = pd.concat(dfs, axis = 1)

    ## write to output file
    out_df.to_csv(out_file, sep = "\t", index = False, header = True)

main(sys.argv[1], sys.argv[2], sys.argv[3])
