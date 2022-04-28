#!/bin/bash

# Calculate MD5 checksums of fastq's for submission to GEO (from April)
# can get list of file names using command $ ls -1a
# checksums must be calculated for unzipped FASTQs

# TODO: Make fofn a positional argument
fofn="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/demultiplex/geo_files/fastq_filenames.txt"
checksumsfile="/home/pparrish/bergerlab_shared/Projects/paralog_pgRNA/demultiplex/geo_files/md5checksums.txt"

readarray -t fastqlist < $fofn
md5sum "${fastqlist[@]}" > $checksumsfile
