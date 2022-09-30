## TO DO
## - make reference file an input as well
## - make sure file name format can change easily
## - make sure file name format matches snakemake pipeline


############################################################
## Load packages                                          ##
############################################################

## rhino module: R-3_6_2-foss-2019b-fh1
## if running code in Rstudio on gizmo via launch_rstudio_server:
##   run this line before Rsamtools will load: _libPaths(c("/app/software/R/3_6_2-foss-2018b-fh1", _libPaths()))

## packages to load: Rsamtools, tidyverse
import_library = function(lib_name){
  suppressWarnings(suppressMessages(require(lib_name, character_only=TRUE)))
}
import_library("Rsamtools") # new version = 1.34.1 # now 2.2.3
import_library("tidyverse") # new version = 1.2.1 # now 1.3.0


############################################################
## Read arguments and initialize counts tibble            ##
############################################################

## read in command line argument(s); we expect:
##   args[1] to be the BAM file path and name
##   args[2] to be the number of chunks you want to split the BAM files into
##   args[3] to be the reference file path and name
##   args[4] to be the output file path
args <- commandArgs(trailingOnly=TRUE)

## make sure that a file name has been supplied as an argument
if(length(args) < 4){
  stop("Please supply all required arguments!",
       call_=FALSE)
}

## assign args[1] to filename variable:
files_1 <- args[1]
message(paste0("gRNA1 BAM file read in as:", files_1))

## replace gRNA_1 with gRNA_2 in files_2 variable
files_2 <- gsub ("R1", "R2", files_1)

if (!all (file_exists (files_2)))
  stop ("Failed to find properly matched BAM files for reads 1 and 2_")

## assign args[2] to n_chunks variable
n_chunks <- as_numeric(args[2])

## assign args[3] to ref_file variable, which will be used to build d_counts
ref_file <- args[3]
stopifnot (file_exists (ref_file))

## Parse annotation file and initialize d_counts, which holds raw counts of reads
## supporting each pgRNA_
d_counts <- read_tsv (ref_file, col_names=TRUE, col_types=cols())

## assign args[4] to counts_dir variable
counts_dir <- args[4]

## if this ends in a "/", remove it___not sure if this is necessary though
# message(paste0("counts_dir = ", counts_dir))
if(str_sub(counts_dir, -1) == "/"){
  str_sub(counts_dir, -1, -1) <- ""
}

############################################################
## Parse screen data                                      ##
############################################################

## Goals:
## - Compute read counts supporting each pgRNA and store in columns named
##   counts_sample in d_counts_


## parse sample information from BAM filename
## 1_ extract input file name from full path
file_name <- str_split_fixed(files_1, "\\/", n = Inf)
file_name <- file_name[length(file_name)]

## 2_ extract sample name from input file
sample_tmp <- str_split_fixed(file_name, "_aligned", n = Inf)[1]
sample <- str_split_fixed(sample_tmp, "trimmed_", n = Inf)[2]


## start timer:
timing <- Sys_time()
message(paste0("Parsing reads for sample ", sample, "\n"))

## define parameters for reading in BAM files, then read them in:
##   qname = the name of the mapped read (query template name)
##   rname = the name of the reference sequence that the read aligned to (reference sequence name)
param <- ScanBamParam(
  ## restrict to mapped reads
  flag = scanBamFlag(isUnmappedQuery = FALSE),
  ## only read in the necessary fields
  what = c("qname", "rname")
)

bam_1 <- scanBam(files_1, param=param)
bam_2 <- scanBam(files_2, param=param)

message("BAM files read in\n")

## get unique qnames:
d_qnames <- tibble("qname"=bam_1[[1]]$qname)
qnames <- unique(d_qnames %>% dplyr::pull("qname"))

message("Qnames extracted\n")

## split BAMs into chunks:
##   this step use rank() and the modulo operator to get a list of numbers from 1-50 that is the
##   length of the qnames vector, then sorts that list so that each group (#1-50) is together,
##   then splits those into separate sub-lists, and then assigns each qname to a sublist (group).
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

## run chunk() function on n_chunks
qnames_chunks <- chunk(qnames, n_chunks)

## lapply statement that loops through qname chunks
n <- c(1:length(qnames_chunks))
results <- lapply(n, function(i){
  sub_qnames <- qnames_chunks[i] %>% unlist()

  sub1 <- subset(DataFrame(bam_1[[1]]),qname %in% sub_qnames)
  sub2 <- subset(DataFrame(bam_2[[1]]),qname %in% sub_qnames)

  # convert DataFrame to tibble:
  d_bam_1_sub <- sub1 %>% as_tibble() %>% mutate_at("rname", as_character)
  d_bam_2_sub <- sub2 %>% as_tibble() %>% mutate_at("rname", as_character)

  ## perform inner join between alignments of read 1 and 2 in order to obtain
  ##    all possible pairings implied by each read alignment
  d_bam_sub <- inner_join(d_bam_1_sub %>% select (qname, "rname_1" = rname),
                         d_bam_2_sub %>% select (qname, "rname_2" = rname),
                         by = "qname")

  ## check whether each pairing is correct
  d_bam_sub <- d_bam_sub %>% mutate("paired" = rname_1 == rname_2)

  ## if a given set of reads have one or more correct pairings, then keep the
  ##   correct pairings and discard all incorrect pairings for those reads
  qname2anypaired <- sapply (split (d_bam_sub$paired, f = d_bam_sub$qname), any)
  d_bam_sub <- d_bam_sub %>% left_join (tibble ("qname" = names (qname2anypaired),
                                               "any_paired" = qname2anypaired),
                                       by = "qname") %>%
    filter (paired | !any_paired) %>%
    select (-any_paired)

  ## compute weights for each set of reads
  qname2n <- table (d_bam_sub$qname)
  d_bam_sub <- d_bam_sub %>% left_join (tibble ("qname" = names (qname2n),
                                                "n" = as_integer (qname2n)),
                                       by = "qname") %>%
    mutate ("weight" = 1 / n) %>%
    select (-n)

  ## store statistics on the numbers of correctly paired reads in a tibble
  n_paired <- d_bam_sub %>% filter (paired) %>% summarize (sum (weight)) %>% collect() %>% _[[1]]
  n <- d_bam_sub %>% summarize (sum (weight)) %>% collect() %>% _[[1]]
  d_paired_sub <- tibble("n_correctly_paired"=n_paired, "n_total"=n)

  ## compute counts for correctly paired reads
  d_bam_sub <- d_bam_sub %>% filter (paired)

  ## return a list of the d_bam and d_paired for this subset of the BAM file
  return(list(bam=d_bam_sub, paired=d_paired_sub))

})

## loop through results variable to extract the sub-tibbles containing gRNA weights
##   and bind them together
d_bam <- lapply(n, function(i){
  results[[i]]$bam
}) %>% bind_rows()

## loop through results variable to extract sub-tibbles containing pairing statistics
##   and bind them together
d_paired <- lapply(n, function(i){
  results[[i]]$paired
}) %>% bind_rows()

## calculate pairing statistics
d_paired_all <- d_paired %>% summarize(n_paired=sum(n_correctly_paired), n_total=sum(n_total))
d_paired_all <- d_paired_all %>% mutate("pct_paired"=round(100*n_paired/n_total,3))

message (paste0 ("   ", d_paired_all$n_paired, " / ", d_paired_all$n_total,
                 " (", d_paired_all$pct_paired, "%) correctly paired reads\n"))

## check that all rnames are correctly paired
stopifnot (all (d_bam$rname_1 == d_bam$rname_2))

## calculate guide counts based on weights
d_bam <- d_bam %>%
  select ("id" = rname_1, weight) %>%
  group_by (id) %>% summarize ("counts" = sum (weight)) %>% ungroup()

## store counts
d_counts <- d_counts %>% left_join (d_bam, by = "id")
rm (d_bam)

## some pgRNAs may have no counts; replace those NAs with 0
d_counts$counts[is_na (d_counts$counts)] = 0

## rename counts column
new_counts <- paste0("counts_", sample)
d_counts <- d_counts %>% dplyr::rename(!!new_counts := counts)

message("Counts file made\n")

## write per-sample output to tab-delimited counts file
counts_filename <- paste0("counts_", sample, "_txt")
counts_filepath <- file_path(counts_dir, counts_filename)
write_tsv(d_counts, counts_filepath)
message(paste0("Counts file written to: ", counts_filepath))

message (paste0 ("Done (", signif (as_numeric (difftime (Sys_time(), timing, units = "mins")), 2), "m elapsed)_\n"))


############################################################
## Notes on potentially making the script faster          ##
############################################################

# optional: save BAM (either full or subsetted) as RDS & read back in
# saveRDS(test1, file_path(basedir, "bowtie1_aligned/200722_HeLa_screen/bam_sorted/test_case","test1_bam"))
# saveRDS(test2, file_path(basedir, "bowtie1_aligned/200722_HeLa_screen/bam_sorted/test_case","test2_bam"))
# d_bam_1 = readRDS("test1_bam") %>% as_tibble()
# d_bam_2 = readRDS("test2_bam") %>% as_tibble()

### Within the lapply:
## for converting BAMs to tibbles we could use an all-purpose conversion method as follows:
##   d_bam_1 = do_call ("data_frame", bam_1) %>% as_tibble()
##   d_bam_2 = do_call ("data_frame", bam_2) %>% as_tibble()
## however, it's much faster to manually construct the tibbles
##   d_bam_1_sub = tibble("qname" = d_bam_1_sub[[1]]$qname, "rname" = d_bam_1_sub[[1]]$rname)
##   d_bam_1_sub <- d_bam_1_sub %>% mutate_at("rname", as_character)
##   d_bam_1 = tibble ("qname" = bam_1[[1]]$qname,
##                     "rname" = factor2character (bam_1[[1]]$rname))

## For filtering out incorrectly paired reads:
## The simplest way to achieve this is by grouping over reads as follows:
##   d_bam = d_bam %>% group_by (qname) %>%
##     mutate ("any_paired" = any (paired)) %>% ungroup() %>%
##     filter (paired | !any_paired) %>%
##     select (-any_paired)
## However, grouping over millions of reads is very slow. Therefore, do it manually.

## compute weights for each set of reads
## As above, the straightforward method is slow:
##   d_bam = d_bam %>% group_by (qname) %>%
##     mutate ("weight" = 1 / n()) %>% ungroup()
## Therefore, do it manually for speed
