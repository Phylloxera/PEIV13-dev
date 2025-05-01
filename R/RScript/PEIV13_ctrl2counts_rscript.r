#!/usr/bin/env Rscript

#The PEIV13_ctrl2counts_rscript.R Rscript file. intended to be run via 
#system() in RStudio or inside of the PEIV13 apptainer container
#tested mac/pc but too much data for use with mac (onedrive problem or system?)
#not yet tested linux.

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

#check files/folders & packages
checkfs(); checkenv()

#get & order metadata md. check assignTaxonomy files/parameters. set up results
check_taxdb()
cat("first check of metadata after comparing file names against library\n")
md <- read.csv(pattern2string("metadatafile"), tryLogical = F); check_mdat()
md <- md[with(md, order(as.character(run), as.character(library) )), ]
setup_peiv13_results()

#define R1, R2, refine R1, R2 & md, samples and test
fns <- sort(list.files(pattern2string("inputfolder"), full.names = T, 
                       recursive = T)); fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]; fnFs <- unique(batch_ilike_subset_vec(fnFs, md[
  , "library"])); fnRs <- unique(batch_ilike_subset_vec(fnRs, md[, "library"]))

#subset md
md <- subset_md()

#extract trunc,maxerror,trimLeft,run variables from uservars
varMessages(); forwardmaxerror <- dada2var_from_ctrl("forwardmaxerror")
reversemaxerror <- dada2var_from_ctrl("reversemaxerror")
forwardtrunc <- dada2var_from_ctrl("forwardtrunc")
reversetrunc <- dada2var_from_ctrl("reversetrunc")
forwardtrimLeft <- dada2var_from_ctrl("forwardtrimLeft")
reversetrimLeft <- dada2var_from_ctrl("reversetrimLeft")
runs <- unique(md[, "run"]); cat(length(runs), "runs found in metadata\n")

#check environment variables
PEIV13varcheck()

#split libraries by runs. get ready to run dada2.
FileList <- split(md["library"], md[, "run"]); suppressPackageStartupMessages(
  library(dada2))

#run dada2 on each run in loop.
for (r in runs) {raw2seqtab()} #did not complete on mac. try connected to
#screen downstairs!

#merge sequence tables. test updated with subsetted md/runs
st.all <- multiseqtab()

#merge repeat rows. ""
st.mr <- mat_sum_repeat()

#remove chimeras. ""
st.chim <- Bn()

#assign taxonomy. ""
taxa <- aT()

#subset taxa and st.bac (tm if necessary)
st.bac <- rmBac(); taxa.bac <- rmBac2()

#save st.bac, taxa.bac
message("Writing sequence table and taxonomy table")
write.csv(st.bac, "results/SequenceTable.csv")
write.csv(taxa.bac, "results/TaxaTable.csv")

#make a phyloseq (tree-less) and save
make_ps()

#unload dada2
pkg <- "package:dada2"; detach(pkg, character.only = T)

#cleanup intermediate
rm_peiv13_int()

#i can add in fasttree2 and ribovore (others?) BUT only inside container and
#with separate rscripts