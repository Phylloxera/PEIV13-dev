#!/usr/bin/env Rscript

#The PEIV13_ctrl2counts_rscript.R Rscript file. Tested via system() in RStudio, 
#inside of the peiv13cont apptainer container, and via Rscript

beg_time <- Sys.time();cat("Running c2c!\n") #Rscript time accounting

#commandArgs() test
if (length(commandArgs(trailingOnly = T)) != 1) {
  stop("This script requires a controlfile argument; exiting", call. = F)}
#regular commandArgs() to get self and find source dir
script.name <- sub("--file=", "", commandArgs()[grep("--file=", commandArgs())])
#this should just give the same value on win. used only in mac.
script.name <- gsub("\\~\\+\\~", " ", script.name)
Rsrcdir <- file.path(dirname(dirname(script.name)), "src")
source(paste0(Rsrcdir, "/PEIV13funcs.R")) #get funcs

#trailingOnly TRUE to get controlfile path from commandArgs()
uservars <- getc(commandArgs(trailingOnly = T)[1])

checkfs(); checkenv() #check files/folders & packages

#get & order metadata md. check assignTaxonomy files/parameters. set up results
check_taxdb()
cat("Checking metadata after comparing file names against library!\n")
md <- read.csv(pattern2string("metadatafile"), tryLogical = F); check_mdat()
md <- md[with(md, order(as.character(run), as.character(library) )), ]
setup_peiv13_results(); saveRDS(md, "results/metadata1.rds")

#define R1, R2, refine R1, R2 & md, samples and test
fns <- sort(list.files(pattern2string("inputfolder"), full.names = T, 
                       recursive = T)); fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]; fnFs <- unique(batch_ilike_subset_vec(fnFs, md[
  , "library"])); fnRs <- unique(batch_ilike_subset_vec(fnRs, md[, "library"]))

#subset and save final metadata file
md <- subset_md(); cat("Saving final metadata file!\n"); saveRDS(
  md, "results/metadata2.rds"); f <- list.files(
    "results", "metadata", full.names = T); if (nrow(readRDS(f[1])) == nrow(
      readRDS(f[2]))) {unlink(f[1]); file.rename(f[2], "results/metadata.rds")
      if (file.exists(f[1])) stop("failed to remove intermediate files")}

#extract trunc,maxerror,trimLeft,run variables from uservars
varMessages(); forwardmaxerror <- dada2var_from_ctrl("forwardmaxerror")
reversemaxerror <- dada2var_from_ctrl("reversemaxerror")
forwardtrunc <- dada2var_from_ctrl("forwardtrunc")
reversetrunc <- dada2var_from_ctrl("reversetrunc")
forwardtrimLeft <- dada2var_from_ctrl("forwardtrimLeft")
reversetrimLeft <- dada2var_from_ctrl("reversetrimLeft")
runs <- unique(md[, "run"]); word <- "run"; vct <- length(runs)
if (vct == 1L) {word <- word} else {word <- paste0(word, "s")}
message("--", vct, " ", word, " found in metadata.--")

#check environment variables (allow for optional raw2seqtab variables)
pl <- as.character(pattern2string("pool")); if (length(pl) > 0) {
  #For a single optional variable, pl is character vector length 1.
  #For multiple optional variables, this might need to be a list.
  pl <- PEIV13varcheck()} else {remove(pl); PEIV13varcheck()}

#split libraries by runs. get ready to run dada2.
FileList <- split(md["library"], md[, "run"]); suppressPackageStartupMessages(
  library(dada2))

#run dada2 on each run (optional variable aware)
if (exists('pl')) {for (r in runs) {raw2seqtab(plv = pl)} 
} else {for (r in runs) {raw2seqtab()}}

st.all <- multiseqtab() #merge sequence tables

st.mr <- mat_sum_repeat() #merge repeat rows

st.chim <- Bn(); taxa <- aT() #remove chimeras, assign taxonomy

#subset taxa and remove nonbacterial sequences
#(match to track reads summaries if necessary)
st.bac <- rmBac(); taxa.bac <- rmBac2()

#save files (written as rds not csv for container version 1.1)
cat("Writing sequence table and taxonomy table!\n")
saveRDS(st.bac, "results/SequenceTable.rds")
saveRDS(taxa.bac, "results/TaxaTable.rds")

make_ps() #make a phyloseq (tree-less) and save

pkg <- "package:dada2"; detach(pkg, character.only = T) #unload dada2

rm_peiv13_int() #cleanup intermediate files/folders

#Rscript time accounting
end_time_h <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours")
cat("c2c completed in", paste0(elapsed_time, "!\n"))
