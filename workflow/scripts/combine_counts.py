import os
import sys
import re
import pandas as pd

def main(annot_file, in_dir, out_file):
    ## read in pgRNA annotations to get IDs
    annot_df = pd.read_csv(annot_file, sep = "\t")

    ## loop through files in pgRNA_counts_dir
    dfs = [annot_df]
    with os.scandir(in_dir) as iterator:
        ## sorted/lambda stuff is to order files alphabetically
        for entry in sorted(iterator, key = lambda e: e.name):
            if not entry.name.startswith('.') and entry.is_file():
                df = pd.read_csv(entry, sep = "\t")
                counts = df.loc[:, df.columns.str.contains("counts")] ## extract counts col
                dfs.append(counts)

    out_df = pd.concat(dfs, axis = 1)

    ## write to output file
    out_df.to_csv(out_file, sep = "\t", index = False, header = True)

main(sys.argv[1], sys.argv[2], sys.argv[3])

sys.stdout.write("Done!")
