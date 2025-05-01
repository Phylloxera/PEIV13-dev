#!/usr/bin/env Rscript

#The PEIV13_ctrl2rqf_rscript.R Rscript file. can be run via 
#system() in RStudio or via Rscript. Possibly inside container.

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
mf <- pattern2string("metadatafile")
if (!file.exists(mf)) stop(mf, " not found"); rawfol <- pattern2string(
  "inputfolder"); if (!dir.exists(rawfol)) stop(rawfol, " not found")
if (!("dada2" %in% installed.packages())) stop(
  "dada2 is reqired for PEIV13 but is not an installed package in your R insta",
  "llation")

#get & order metadata md
cat("first check of metadata after comparing file names against library\n")
md <- read.csv(mf, tryLogical = F); check_mdat()
md <- md[with(md, order(as.character(run), as.character(library) )), ]

#define R1, R2, refine R1, R2 & md, samples and test
fns <- sort(list.files(pattern2string("inputfolder"), full.names = T, 
                       recursive = T)); fnFs <- fns[grepl("R1", fns)]
fnRs <- fns[grepl("R2", fns)]; fnFs <- unique(batch_ilike_subset_vec(fnFs, md[
  , "library"])); fnRs <- unique(batch_ilike_subset_vec(fnRs, md[, "library"]))
md <- subset_md()
cat("second check of metadata after comparing file names against library\n")
check_mdat()
runs <- unique(md[, "run"]); cat(length(runs), "runs found in metadata\n")

#split libraries by runs. get ready to run dada2.
FileList <- split(md["library"], md[, "run"]); suppressPackageStartupMessages(
  library(dada2))

#create and save read quality profiles for each run in loop.
for (r in runs) {raw2pQP()}

#unload dada2
pkg <- "package:dada2"; detach(pkg, character.only = T)