############################################################
## Load packages                                          ##
############################################################

## rhino module: R-3.6.2-foss-2019b-fh1
## if running code in Rstudio on gizmo via launch_rstudio_server:
##   run this line before Rsamtools will load: .libPaths(c("/app/software/R/3.6.2-foss-2018b-fh1", .libPaths()))

## packages to load: Rsamtools, tidyverse
import_library = function(lib_name){
  suppressWarnings(suppressMessages(require(lib_name, character.only=TRUE)))
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
       call.=FALSE)
}

## assign args[1] to filename variable:
files.1 <- args[1]
message(paste0("gRNA1 BAM file read in as:", files.1))

## replace gRNA_1 with gRNA_2 in files.2 variable
files.2 <- gsub ("R1", "R2", files.1)

if (!all (file.exists (files.2)))
  stop ("Failed to find properly matched BAM files for reads 1 and 2.")

## assign args[2] to n.chunks variable
n.chunks <- as.numeric(args[2])

## assign args[3] to ref.file variable, which will be used to build d.counts
ref.file <- args[3]
stopifnot (file.exists (ref.file))

## Parse annotation file and initialize d.counts, which holds raw counts of reads
## supporting each pgRNA.
d.counts <- read_tsv (ref.file, col_names=TRUE, col_types=cols())

## assign args[4] to counts.dir variable
counts.dir <- args[4]

## if this ends in a "/", remove it...not sure if this is necessary though
# message(paste0("counts.dir = ", counts.dir))
if(str_sub(counts.dir, -1) == "/"){
  str_sub(counts.dir, -1, -1) <- ""
}

############################################################
## Parse screen data                                      ##
############################################################

## Goals:
## - Compute read counts supporting each pgRNA and store in columns named
##   counts_sample in d.counts.


## parse sample information from BAM filename
## 1. extract input file name from full path
file_name <- str_split_fixed(files.1, "\\/", n = Inf)
file_name <- file_name[length(file_name)]

## 2. extract sample name from input file
sample_tmp <- str_split_fixed(file_name, "_aligned", n = Inf)[1]
sample <- str_split_fixed(sample_tmp, "trimmed_", n = Inf)[2]


## start timer:
timing <- Sys.time()
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

bam.1 <- scanBam(files.1, param=param)
bam.2 <- scanBam(files.2, param=param)

message("BAM files read in\n")

## get unique qnames:
d.qnames <- tibble("qname"=bam.1[[1]]$qname)
qnames <- unique(d.qnames %>% dplyr::pull("qname"))

message("Qnames extracted\n")

## split BAMs into chunks:
##   this step use rank() and the modulo operator to get a list of numbers from 1-50 that is the
##   length of the qnames vector, then sorts that list so that each group (#1-50) is together,
##   then splits those into separate sub-lists, and then assigns each qname to a sublist (group).
chunk <- function(x,n) split(x, factor(sort(rank(x)%%n)))

## run chunk() function on n.chunks
qnames_chunks <- chunk(qnames, n.chunks)

## lapply statement that loops through qname chunks
n <- c(1:length(qnames_chunks))
results <- lapply(n, function(i){
  sub_qnames <- qnames_chunks[i] %>% unlist()

  sub1 <- subset(DataFrame(bam.1[[1]]),qname %in% sub_qnames)
  sub2 <- subset(DataFrame(bam.2[[1]]),qname %in% sub_qnames)

  # convert DataFrame to tibble:
  d.bam.1.sub <- sub1 %>% as_tibble() %>% mutate_at("rname", as.character)
  d.bam.2.sub <- sub2 %>% as_tibble() %>% mutate_at("rname", as.character)

  ## perform inner join between alignments of read 1 and 2 in order to obtain
  ##    all possible pairings implied by each read alignment
  d.bam.sub <- inner_join(d.bam.1.sub %>% select (qname, "rname.1" = rname),
                         d.bam.2.sub %>% select (qname, "rname.2" = rname),
                         by = "qname")

  ## check whether each pairing is correct
  d.bam.sub <- d.bam.sub %>% mutate("paired" = rname.1 == rname.2)

  ## if a given set of reads have one or more correct pairings, then keep the
  ##   correct pairings and discard all incorrect pairings for those reads.
  qname2anypaired <- sapply (split (d.bam.sub$paired, f = d.bam.sub$qname), any)
  d.bam.sub <- d.bam.sub %>% left_join (tibble ("qname" = names (qname2anypaired),
                                               "any_paired" = qname2anypaired),
                                       by = "qname") %>%
    filter (paired | !any_paired) %>%
    select (-any_paired)

  ## compute weights for each set of reads
  qname2n <- table (d.bam.sub$qname)
  d.bam.sub <- d.bam.sub %>% left_join (tibble ("qname" = names (qname2n),
                                                "n" = as.integer (qname2n)),
                                       by = "qname") %>%
    mutate ("weight" = 1 / n) %>%
    select (-n)

  ## store statistics on the numbers of correctly paired reads in a tibble
  n.paired <- d.bam.sub %>% filter (paired) %>% summarize (sum (weight)) %>% collect() %>% .[[1]]
  n <- d.bam.sub %>% summarize (sum (weight)) %>% collect() %>% .[[1]]
  d.paired.sub <- tibble("n.correctly.paired"=n.paired, "n.total"=n)

  ## compute counts for correctly paired reads
  d.bam.sub <- d.bam.sub %>% filter (paired)

  ## return a list of the d.bam and d.paired for this subset of the BAM file
  return(list(bam=d.bam.sub, paired=d.paired.sub))

})

## loop through results variable to extract the sub-tibbles containing gRNA weights
##   and bind them together
d.bam <- lapply(n, function(i){
  results[[i]]$bam
}) %>% bind_rows()

## loop through results variable to extract sub-tibbles containing pairing statistics
##   and bind them together
d.paired <- lapply(n, function(i){
  results[[i]]$paired
}) %>% bind_rows()

## calculate pairing statistics
d.paired.all <- d.paired %>% summarize(n.paired=sum(n.correctly.paired), n.total=sum(n.total))
d.paired.all <- d.paired.all %>% mutate("pct.paired"=round(100*n.paired/n.total,3))

message (paste0 ("   ", d.paired.all$n.paired, " / ", d.paired.all$n.total,
                 " (", d.paired.all$pct.paired, "%) correctly paired reads\n"))

## check that all rnames are correctly paired
stopifnot (all (d.bam$rname.1 == d.bam$rname.2))

## calculate guide counts based on weights
d.bam <- d.bam %>%
  select ("id" = rname.1, weight) %>%
  group_by (id) %>% summarize ("counts" = sum (weight)) %>% ungroup()

## store counts
d.counts <- d.counts %>% left_join (d.bam, by = "id")
rm (d.bam)

## some pgRNAs may have no counts; replace those NAs with 0
d.counts$counts[is.na (d.counts$counts)] = 0

## rename counts column
new.counts <- paste0("counts_", sample)
d.counts <- d.counts %>% dplyr::rename(!!new.counts := counts)

message("Counts file made\n")

## write per-sample output to tab-delimited counts file
counts.filename <- paste0("counts_", sample, ".txt")
counts.filepath <- file.path(counts.dir, counts.filename)
write_tsv(d.counts, counts.filepath)
message(paste0("Counts file written to: ", counts.filepath))

message (paste0 ("Done (", signif (as.numeric (difftime (Sys.time(), timing, units = "mins")), 2), "m elapsed).\n"))
