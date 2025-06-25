#PEIV13funcs.R for PEIV13-dev

#getc takes controlfile path as character string and returns a dataframe
#from the file at that path
getc <- function(cp) {if (!file.exists(cp)) stop(cp, ",\nnot found")
  uv <- read.csv(cp, header = F); if (!is.list(uv)) stop(
    cp, ",\nnot properly coerced to a list"); return(uv)}
#e.g., uservars <- getcf("ctrl_file_path")

#pattern2string(pattern, uv = uservars) takes a patern string and a peiv13
#controlfile uservars and returns the user supplied value from column 2
pattern2string <- function(pattern, uv = uservars) {
  out <- uv[which(uv[, 1] == pattern), 2]}
#e.g., pattern2string("metadatafile")

#checkfs takes uv (default uservars) and performs confirms existence of inputs
checkfs <- function(uv = uservars) {mf <- pattern2string("metadatafile")
if (!file.exists(mf)) stop(mf, ",\nnot found"); rawfol <- pattern2string(
  "inputfolder"); if (!dir.exists(rawfol)) stop(rawfol, ",\nnot found")
dbf <- pattern2string("databasefile"); if (!file.exists(dbf)) stop(
  dbf, ",\nnot found")}
  #never reproduced or found source of env: Rscript\r: No such file or directory
#e.g., checkfs()

#checkenv stops if a required R package is not found in the R installation
checkenv <- function() {if (!("dada2" %in% installed.packages())) stop(
  "dada2 is reqired for PEIV13 but is not an installed package in your R\ninst",
  "allation"); if (!("data.table" %in% installed.packages())) stop(
  "data.table is reqired for PEIV13 but is not an installed package in your R",
  "\ninstallation"); if (!("dplyr" %in% installed.packages())) stop(
  "dplyr is reqired for PEIV13 but is not an installed package in your R\ninst",
  "allation"); if (!("phyloseq" %in% installed.packages())) stop(
  "phyloseq is reqired for PEIV13 but is not an installed package in your R\ni",
  "nstallation")}
#e.g., checkenv()

#check_mdat function takes PEIV13 metadata data.frame default = md and prints 
#test results
check_mdat <- function (m = md) {if ("library" %in% colnames(
  m) && "run" %in% colnames(m)) {cat(
    "Required columns were found in metadata!\n")
  } else {stop("required columns not found in metadata")}
  #ENFORCE unique(<library column>)
  if (nrow(m) > length(unique(m[, "library"]))) {stop(
    "library column MUST be unique for PEIV13")
    #set conditions for (optional) repeat column
    } else if ("repeat." %in% colnames(m) && length(unique(m[
      , "repeat."])) < nrow(m)) {message("--", length(unique(m[
        , "repeat."])), " samples from ", nrow(
          m), " rows accounting for repeats were found in metadata.\nPer-varia",
        "nt sequence counts will be summed for repeats.--")
      } else {message("--No technical replicates in metadata.--")}}
#e.g., check_mdat()

#dada2var_from_ctrl function for PEIV13 extracts a named integer based on 
#variable labels in a PEIV13 controlfile dataframe. ordering by names in 
#as.character() order more robust
dada2var_from_ctrl <- function(pattern, uv = uservars) {namedint <- as.integer(
  uv[which(grepl(pattern, uv[, 1])), 2]); names(namedint) <- sub(
    pattern, "", uv[which(grepl(pattern, uv[, 1])), 1])
  namedint <- namedint[order(names(namedint))]; return(namedint)}
#e.g., forwarderror <- dada2var_from_ctrl("forwarderror")

#subset_md function takes md and fnFs and returns a subsetted md
subset_md <- function(m = md, fF = fnFs, fR = fnRs) {unused <- unused_patterns(
  fF, m[, "library"]); if (length(unused) > 0) m <- m[-unused, ]; libids <- m[
    , "library"]
    #stop the run in length(fnFs) != length(fnRs) or if length(libids) is zero
  if (length(fF) != length(fR)) {stop(
    "equal numbers of forward and reverse raw read files needed for paired-end",
    " protocol")}; if (length(libids) == 0L) {stop(
      "no samples to process. check compatibility between library and\ninputfo",
      "lder")}; return(m)}
#e.g., subset_md()

#varMessages function lacks arguments. just posts messages
varMessages <- function() {message(
  "--Reading maxerror parameters from controlfile.  If you want to speed up",
  "\ndownstream computation, consider tightening (decreasing) these.  If too f",
  "ew\nreads are passing the filter, consider relaxing these, perhaps especial",
  "ly\nreversemaxerror and reducing the trunc parameters to remove low quality",
  " tails.\nRemember though, when choosing trunc you must maintain overlap aft",
  "er truncation\nin order to merge R1 and R2 pairs.--")
  message(
  "--Reading forwardtrunc and reversetrunc parameters from controlfile.  You",
  "r reads\nmust still overlap after truncation in order to merge them later! ",
  " When using a\nless-overlapping primer set, like V1-V3, these parameters mu",
  "st be large enough\nto maintain 20 + biological-length-variation nucleotide",
  "s of overlap between\nthem.--")
  message(
  "--Reading trimLeft parameters from controlfile.  Constant length primers ",
  "at the\nstart of your reads is a common scenario.  The forward and reverse ",
  "parameters\nwill be applied to your R1 and R2 files respectively.  If your ",
  "forward primers\nare in R2, set tryRC to TRUE for assignTaxonomy.  Please d",
  "ouble-check that your\nprimers have been removed by dada2!--")}
#e.g., varMessages()

#multieq function for PEIV13 takes 3 or more R objects of the same class and 
#tests identity returning TRUE or FALSE
multieq <- function(A, B, C, ...) {classes <- c(class(A), class(B), class(C))
input_list <- list(...); elipsecount <- length(input_list); if (
  elipsecount > 0) {for (i in 1:elipsecount) {classes <- c(classes, class(
    input_list[[i]]))}}; if (elipsecount > 0) {for (i in 1:elipsecount) {assign(
      LETTERS[3 + i], input_list[[i]])}}; tcount <- elipsecount + 3
if (all(sapply(list(classes[1:(tcount - 1)]), function(x) x == classes[
  tcount])) == F) {stop("arguments do not share the same class")}; all(sapply(
    mget(LETTERS[1:(tcount - 1)]), function(x) x == get(LETTERS[tcount])))}
#e.g., multieq(names(forwarderror), names(reverseerror), as.character(runs))

#check_taxdb function takes uservars. includes tryRC check
check_taxdb <- function(uv = uservars) {cat(
  "Checking Taxonomy database files/parameters!\n")
  trc <- as.logical(pattern2string("tryRC")); if (is.na(trc)) stop(
    "tryRC needs to be TRUE or FALSE"); d <- pattern2string("databasefile")
  #what happens with below when multiple trailing ';'
  L1 <- readLines(d, n = 1); if (length(unlist(strsplit(L1, ';'))) != 6) stop(
    d, "\ndoes not appear to have 6 taxonomic ranks as expected for Genus-leve",
    "l identification"); if (substring(L1,1,1) != ">") stop(
      d, "\nis not a FASTA file")}
#e.g., check_taxdb()

#setup_peiv13_results function takes uservars.
#message() is more about what the results might mean
setup_peiv13_results <- function(uv = uservars) {ow <- as.logical(
  pattern2string("overwrite")); if (is.na(ow)) stop(
    "overwrite needs to be TRUE or FALSE")
  #if results exists and overwrite is FALSE exit
if (dir.exists("results") && ow == F) stop(
  "the PEIV13 results directory,\n", getwd(),
  "/results,\nwill not be overwritten since overwrite is FALSE")
#if results exists and overwrite is TRUE delete and recreate else create
#force = T needed on windows one drive
if (dir.exists("results") && ow == T) {message(
  "--Overwriting,\n", getwd(), "/results.--"); unlink(
  "results", recursive = T, force = T); if (dir.exists("results")) {stop(
    "the PEIV13 results directory,\n", getwd(), "/results,\ncould not be overw",
    "ritten")}; dir.create("results"); dir.create("results/filt")
  dir.create("results/seqtab"); dir.create("results/plotErrors")} else {
    message("--Creating\n", getwd(), "/results.--"); dir.create("results")
    dir.create("results/filt"); dir.create("results/seqtab"); dir.create(
      "results/plotErrors")}}
#e.g., setup_peiv13_results()

#check_peiv13_results function
check_peiv13_results <- function(uv = uservars) {ow <- as.logical(
  pattern2string("overwrite")); if (is.na(ow)) stop(
    "overwrite needs to be TRUE or FALSE")
#if results does not exist exit
if (!dir.exists("results")) stop(
  "the PEIV13 results directory,\n", getwd(), "/results,\ndoes not exist")
#if seqtab does not exist exit
if (!dir.exists("results/seqtab")) {stop(
  "the PEIV13 results/seqtab directory,\n", getwd(), 
  "/results/seqtab,\ndoes not exist")}}
#e.g., check_peiv13_results()

#the unused_patterns function takes a long vector (subject vector) and a
#pattern vector (query vector) (no defaults) and returns the indices of the 
#long/subject, which do not match any of the patterns
unused_patterns <- function(longvec, pvec) {up <- integer(0)
suppressPackageStartupMessages(library(data.table)); for (i in 1:length(pvec)) {
  if (length(which(longvec %ilike% pvec[i])) == 0L) {up <- c(up, i)}}
pkg <- "package:data.table"; detach(pkg, character.only = T); return(up)}
#e.g., unused <- unused_patterns(fnFs, md[, "library"])

#the batch_ilike_subset_vec function takes a subset (subject) vector and a 
#pattern vector (query) (no defaults) and returns the subject with entries
#removed that are not %ilike% (package data.table) any of the patterns.
batch_ilike_subset_vec <- function(ssvector, pvector) {index <- integer()
suppressPackageStartupMessages(library(data.table)); for (i in 1:length(
  pvector)) {index <- c(index, which(ssvector %ilike% pvector[i]))}
pkg <- "package:data.table"; detach(pkg, character.only = T)
out <- ssvector[index]}
#e.g., fnRs <- unique(batch_ilike_subset_vec(fnRs, md[, "library"]))

#PEIV13varcheck function. no arguments. performs final checks of existing
#R objects before running dada2
PEIV13varcheck <- function() {if (is.integer(reversetrimLeft) == F) stop(
  "reversetrimLeft needs to be an integer"); if (is.integer(
    forwardtrimLeft) == F) stop("forwardtrimLeft needs to be an integer")
  if (is.integer(reversetrunc) == F) stop("reversetrunc needs to be an integer")
  if (is.integer(reversemaxerror) == F) stop(
    "reversemaxerror needs to be an integer"); if (is.integer(
      forwardtrunc) == F) stop("forwardtrunc needs to be an integer")
  if (is.integer(forwardmaxerror) == F) stop(
    "forwardmaxerror needs to be an integer")
  if (length(forwardtrimLeft) > 1) stop("too many forwardtrimLeft values found")
  if (length(reversetrimLeft) > 1) stop("too many reversetrimLeft values found")
  if (multieq(names(forwardmaxerror), names(reversemaxerror), names(
    forwardtrunc), names(reversetrunc), as.character(runs)) == T) {cat(
      "User supplied dada2 trunc and error variables match!\n")
    } else {stop("check that error and trunc controlfile labels match your met",
                 "adata run column")}
  #new container V 1.1 text for new ... variables that can be passed to
  #raw2seqtab(). if add a second one, pl will need to be returned as a list
  pl <- as.character(pattern2string("pool")); if (length(pl) > 0) {
    if (pl != "pseudo") pl <- as.logical(casefold(pl, upper = T)); if (is.na(
      pl)) stop('pool needs to be TRUE, FALSE or "pseudo"')} 
  cat("Performing second check of metadata after comparing file names against",
      "library!\n"); check_mdat()
  message(
    "--PEIV13 will now run the dada2 pipeline for paired-end bacterial 16s v", 
    "ariable\nregion 1-3 data.  Error plots will be in\nresults/plotErrors/<di",
    "rection>_<batch>.pdf showing the error rates for each\npossible transitio",
    "n (A→C, A→G, ...).  Points are the observed error rates for\neach consens",
    "us quality score.  The black line shows the estimated error rates\nafter ",
    "convergence of the machine-learning algorithm and the red line shows the",
    "\nerror rates expected under the nominal definition of the Q-score.  If t",
    "he\nestimated error rates (black line) are a good fit to the observed rat",
    "es(points)\nand the error rates drop with increased quality as expected y",
    "ou can have\nconfidence in the error model.  Example plots that might jus",
    "tify enforcing\nmonotonicity in the fitted error model are here:\nhttps:/",
    "/github.com/benjjneb/dada2/issues/1156.  During merging, sequences are\no",
    "nly output if the forward and reverse reads overlap by at least 12 bases,",
    " and\nare identical to each other in the overlap region.--")
  #new text for container version 1.1. return optional variable. needs to have
  #assignment inside PEIV_ctrl2counts_rscript.r
  if (length(pl) > 0) return(pl)}
#e.g., PEIV13varcheck()

#raw2pQP requires metadata and controlfile. makes a folder called qualProf if
#not present and writes plots to pdf, labeling according to library and run.
raw2pQP <- function(fnF = fnFs, fnR = fnRs, FL = FileList) {r <- as.character(r)
res <- "qualProf"; if (!dir.exists(res)) dir.create(res); if (!dir.exists(
  res)) stop("PEIV13 was unable to create a quality plots folder"); s <- unlist(
    FL[r]); Fs <- unique(batch_ilike_subset_vec(fnF, s)); Rs <- unique(
      batch_ilike_subset_vec(fnR, s))
#plot quality profiles
cat("Printing quality profile plots for run", r, "!\n")
beg_time <- Sys.time(); if (length(Fs) == 1) {pqF <- suppressWarnings(
  plotQualityProfile(Fs)); warning(
    "a single R1 file is available to plot a quality profile for run ", r)
} else if (length(Fs) > 1) {pqF <- suppressWarnings(plotQualityProfile(Fs[1:2]))
} else {stop("no R1 libraries were found from run ", r)}; pdf(file.path(
  res, paste0(r, "_R1.pdf"))); suppressWarnings(print(pqF))
garbage <- dev.off()
if (length(Rs) == 1) {pqR <- suppressWarnings(plotQualityProfile(Rs))
warning("a single R2 file is available to plot a quality profile for run ", r)
} else if (length(Rs) > 1) {pqR <- suppressWarnings(plotQualityProfile(Rs[1:2]))
} else {stop("no R2 libraries were found from run ", r)}; pdf(file.path(
  res, paste0(r, "_R2.pdf"))); suppressWarnings(print(pqR))
garbage <- dev.off(); end_time_m <- Sys.time(); elapsed_time <- paste0(round(
  as.numeric(difftime(time1 = end_time_m, time2 = beg_time, units = "mins")), 
    2), " Minutes"); cat("The first 2 read quality profiles for run", r, 
                         "completed in", paste0(elapsed_time, "!\n"))}
#e.g., raw2pQP()

#raw2seqtab requires 10 default defined variables. As of container version 1.1,
#it allows for optional variables. the first optional variable is for pooling
#default FALSE as in dada2::dada(). This function writes a dada2 sequence
#table as a rds file in file.path(getwd(), "results/seqtab"). It runs 
#filterAndTrim, learnErrors, dada, mergePairs, and makeSequenceTable functions
#from the dada2 r package for a single run of paired end raw data (only miseq 
#tested)
raw2seqtab <- function(
    fe = forwardmaxerror, re = reversemaxerror, ft = forwardtrunc, 
    rt = reversetrunc, ftl = forwardtrimLeft, rtl = reversetrimLeft, fnF = fnFs,
    fnR = fnRs, FL = FileList, ...) {args_list <- list(...)
plv <- args_list$plv; seed <- 100; multithread_val = T; maxN_val <- 0
truncQ_val <- 2; r <- as.character(r); maxEE_val <- c(fe[r], re[r])
truncspec <- c(ft[r], rt[r]); res <- "results"; filt <- file.path(res, "filt")
s <- unlist(FL[r]); filtFs <- file.path(filt, r, paste0(s, "_F_filt.fastq.gz"))
filtRs <- file.path(filt, r, paste0(s, "_R_filt.fastq.gz"))
#batch_ilike_subset not necessary but some subset is.
Fs <- unique(batch_ilike_subset_vec(fnF, s))
Rs <- unique(batch_ilike_subset_vec(fnR, s))
if (.Platform$OS.type == "windows") multithread_val = F; set.seed(seed)
#filterAndTrim()
beg_time <- Sys.time(); cat("Starting filterAndTrim for run", paste0(r, "!\n"))
out <- filterAndTrim(Fs, filtFs, Rs, filtRs, truncLen = truncspec, trimLeft = c(
  ftl, rtl), maxN = maxN_val, maxEE = maxEE_val, truncQ = truncQ_val, 
  rm.phix = T, compress = T, multithread = multithread_val)
end_time_h <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours")
cat("filterAndTrim for run", r, "completed in", paste0(elapsed_time, "\n"))
cat("Learning errors for run", paste0(r, "!\n")); beg_time <- Sys.time()
#learnErrors()
errF <- learnErrors(
  filtFs, nbases = 1e8, multithread = multithread_val, randomize = T)
pEF <- suppressWarnings(plotErrors(errF, nominalQ = T)); pdf(paste0(
  "results/plotErrors/F_", r, ".pdf")); suppressWarnings(print(pEF))
garbage <- dev.off(); errR <- learnErrors(
  filtRs, nbases = 1e8, multithread = multithread_val, randomize = T)
pER <- suppressWarnings(plotErrors(errR, nominalQ = T)); pdf(paste0(
  "results/plotErrors/R_", r, ".pdf")); suppressWarnings(print(pER))
garbage <- dev.off(); end_time_h <- Sys.time(); elapsed_time <- paste0(round(
  as.numeric(difftime(time1 = end_time_h, time2 = beg_time, units = "hours")), 
  2), " Hours"); cat("learnErrors for run", r, "completed in", paste0(
    elapsed_time, "!\n")); if (length(plv) > 0) {} else {plv = F}; cat(
      "Running dada for run", r, "with", plv, "pooling!\n") #dada()
beg_time <- Sys.time(); dadaFs <- dada(
  filtFs, err = errF, multithread = multithread_val, verbose = F, pool = plv)
dadaRs <- dada(
  filtRs, err = errR, multithread = multithread_val, verbose = F, pool = plv)
end_time_h <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours"); cat(
    "dada for run", r, "completed in", paste0(elapsed_time, "!\n")); cat(
      "Merging pairs for run", paste0(r, "!\n")); beg_time <- Sys.time()
#mergePairs()
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
#additional || nrow(mergers[[1]]) == 0L category needed windows.
if(isTRUE(nrow(mergers) == 0L || nrow(mergers[[1]]) == 0L)) { 
  #report the column 1 reduction in out
  qualfiltper <- round((sum(out[, 1]) - sum(out[, 2])) / sum(out[, 1]) * 100, 2)
  stop("There are no remaining sequences in run ", r, " after mergePairs.\n",
       qualfiltper, "% of reads from run ", r, " were removed by quality filte",
       "ring\nIf this is >>50, your PEIV13 controlfile forward and reversetrun",
       "c parameters\nmight be too small.\nLikewise, your maxerror parameters.",
       "\nRemember that when using a a less-overlapping primer set, like V1-V3",
       ",\ntrunc parameters must be large enough to maintain 20 + biological-l",
       "ength-\nvariation nucleotides of overlap between them.")}
end_time_s <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_s, time2 = beg_time, units = "secs")), 2), " Seconds")
cat("mergePairs for run", r, "completed in", paste0(elapsed_time, "!\n"))
cat("Making sequence table and initiating read-loss tracking for run", paste0(
  r, "!\n")) #makeSequenceTable()
beg_time <- Sys.time(); seqtab <- makeSequenceTable(mergers)
ct_filt400 <- sum(nchar(colnames(seqtab)) < 400); if (ct_filt400 > 0) {
  message("--Removing ", ct_filt400, " variants shorter than 400.  See doi:1",
          "10.1371/journal.pone.0129174--")}; seqtab <- seqtab[, nchar(
            colnames(seqtab)) > 399, drop = F]
hist <- hist(nchar(getSequences(seqtab)), plot = F); hist[["xname"]] <- paste0(
  "Batch ", r, " Variants"); png(paste0("results/seqtab/hist_", r, ".png"))
plot(hist, xlab = 'Lengths'); garbage <- dev.off(); getN <- function(x) sum(
  getUniques(x)); if (nrow(out) > 1L) {track <- cbind(out, sapply(
    dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(
      seqtab))} else {track <- cbind(out, sum(getUniques(dadaFs)), sum(
        getUniques(dadaRs)), sum(getUniques(mergers)), sum(seqtab))}
colnames(track) <- c(
  "input", "filtered", "denoisedF", "denoisedR", "merged", "filter_400")
minDenois <- min(sum(track[, 4]), sum(track[, 3]))
fremain <- round(sum(track[, 2])/sum(track[, 1]) * 100, 2)
dremain <- round(minDenois/sum(track[, 1]) * 100, 2)
mremain <- round(sum(track[, 5])/sum(track[, 1]) * 100, 2)
remain400 <- round(sum(track[, 6])/sum(track[, 1]) * 100, 2); message(
  "--For run ", r, ", ", fremain, "%, ", dremain, "%, ", mremain, "%, and ",
  remain400, "% of input reads remain after\nquality filtering, denoising, mer",
  "ging, and 400 length filtering respectively.\nOutside of filtering, there s",
  "hould be no step in which a majority of reads are\nlost.--")
qualfiltcomp <- (sum(track[, 1]) - sum(track[, 2])) / sum(track[, 1])
if (qualfiltcomp > 0.5) warning(
  "A majority of reads from run ", r, " were removed by quality filtering!\nYo",
  "u might consider relaxing (increasing) the PEIV13 controlfile forward and\n",
  "reversemaxerror parameters used in filterAndTrim, however, your raw data mi",
  "ght\nhave too many errors if a maxerror of 5 is still insufficient.")  
denoisecomp <- (sum(track[, 2]) - minDenois) / sum(track[, 2])
if (denoisecomp > 0.5) warning(
  "A majority of quality filtered reads from run ", r, " were removed by\ndeno",
  "ising!\nYou might consider relaxing (increasing) the PEIV13 controlfile for",
  "ward and\nreversemaxerror parameters used in filterAndTrim, however, your r",
  "aw data might\nhave too many errors if a maxerror of 5 is still insufficien",
  "t."); mergecomp <- (minDenois - sum(track[, 5])) / minDenois
if (mergecomp > 0.5) warning(
  "A majority of quality filtered and denoised reads from run ", r, " were rem",
  "oved by\nmerging!  You should revisit the PEIV13 controlfile forwardtrunc a",
  "nd\nreversetrunc parameters used in filterAndTrim and make sure that the tr",
  "uncated\nreads span your amplicon.")
filter400comp <- (sum(track[, 5]) - sum(track[, 6])) / sum(track[, 5])
if (filter400comp > 0.5) warning(
  "User specified raw reads are probably not V1V3 data!\nA majority of reads f",
  "from run ", r, " were removed by 400 length filtering!")
rownames(track) <- s; rownames(seqtab) <- s; saveRDS(seqtab, paste0(
  "results/seqtab/seqtab", r, ".rds")); saveRDS(track, paste0(
    "results/seqtab/track", r, ".rds")); end_time_s <- Sys.time()
elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_s, time2 = beg_time, units = "secs")), 2), " Seconds"); cat(
    "sequence table and initial read-loss tracking for run", r, "completed in",
    paste0(elapsed_time, "!\n"), fill = T)}
#e.g., raw2seqtab(plv = "pseudo")

#the multiseqtab() function. binds multiple seqtabs.
multiseqtab <- function(r = runs) {seqtabs <- character()
trackfiles <- character(); for (run in r) {run <- as.character(run); seqtabs[
  run] <- c(paste0("results/seqtab/seqtab", run, ".rds")); trackfiles[run] <- c(
    paste0("results/seqtab/track", run, ".rds"))}; cat(
      "Combining track reads!\n"); tracks <- do.call(rbind, lapply(
        trackfiles, readRDS)); saveRDS(tracks, "results/seqtab/trackmult.rds")
if (length(r) == 1L) {s <- readRDS(seqtabs)} else {seed <- 100; cat(
  "Merging Sequence Tables!\n"); set.seed(seed); s <- mergeSequenceTables(
    tables = seqtabs)}; return(s)}
#e.g., multiseqtab()

#rs10k function takes a sequence count matrix s (no default) and returns an 
#index of samples to remove from the sequence count matrix AND track reads 
#matrix. Only be run if any(rowSums(s) < 10000)
rs10k <- function(s) {rs <- which(rowSums(s) < 10000); warning(
  "Samples ", toString(names(rs)), "\nwill be removed due to sequence depth <",
  " 10,000. See doi:10.1093/jas/skab346."); return(rs)}
#e.g., rs10k(st.bac)

#rsSpecify function for debugging as a drop-in replacement for rs10k.
#Added paste() to warning 23 June, 2025
rsSpecify <- function(s, i) {i <- as.integer(i); if (!(is.integer(i))) stop(
  "i is not an integer"); rs <- which(rowSums(s) < i); warning(
    "Samples\n", toString(names(rs)), "\nwill be removed due to arbitrary sequ",
    "ence depth < ", i); return(rs)}
#e.g., rsSpecify(st.bac, 100000)

#rt0 function takes a sequence count matrix s (no default) and returns an index 
#of taxa (sequences) to remove from the sequence count matrix AND taxa 
#identification/classification matrix. Run after rs10k has been run if 
#any(colSums(s) == 0)
rt0 <- function(s) {if (nrow(s) == 0){stop(
  "No samples remain in the sequence table after filtering rows with a\nsequen",
  "ce depth < 10,000.")} else if (nrow(s) > 1L) {rt <- which(colSums(s) == 0)
  } else {rt <- which(s[1, ] == 0)}; warning(length(
    rt), " variants will also be removed as they were only found in the\nsampl",
    "es with sequence depth < 10,000."); return(rt)}
#e.g., rt0(st.bac)

#name_compare2 function takes a sequence count matrix s (no default) and a
#track changes object t (default tm) and performs 2 dim names identity checks
name_compare2 <- function(s, t = tm) {
  mes <- "a problem occurred when removing non Bacterial sequences"
  if (nrow(s) > 1L) {if (any(colSums(s) == 0)) stop(mes)
    } else {if (any(s[1, ] == 0)) stop(mes)}
  if (all(rownames(t) != rownames(s))) stop(mes); return(mes)}
#e.g., name_compare2(st.chim)

#name_compare3 function takes a sequence count matrix s2 (no default) and a
#track changes object t2 (default tm) and a dada2 taxa object tax
#(default taxa.bac) and performs names_compare2 + 1 more dim name identity check
name_compare3 <- function(s2, t2 = trm, tax = taxbac) {mes <- name_compare2(
  s = s2, t = t2); if (all(rownames(tax) != colnames(s2))) stop(mes)}
#e.g., name_compare3(st.bac)

#mat_sum_repeat function takes PEIV13 metadata data.frame default = md and 
#mergeSequenceTable() matrix result default = st.all and merges repeats via
#colSums(). had been getting colSums() errors when nrow was 1L, but then,
#stopped. not replicated.
mat_sum_repeat <- function (m = md, r = runs, FL = FileList, s = st.all) {
  libvec <- unlist(m["library"]); if (all(rownames(s) != libvec)) stop(
    "unable to match seqtab to metadata")
  tracks <- readRDS("results/seqtab/trackmult.rds"); if (all(rownames(
    tracks) != libvec)) stop("unable to match read-loss tracking to metadata")
  saveRDS(tracks, "results/TrackReads_lib.rds"); nfl <- names(FL)
  lr <- length(r); if (nrow(s) > 1L) {cat("Calculating shared variants!\n")
  message("--Few shared variants can indicate technical problems such as inexac", 
          "t region\namplification, incomplete primer removal, heterogeneity sp", 
          "acers and batch\neffects.--")
  for (i in 1:lr) {ss <- s[unlist(FL[which(nfl == as.character(r[
    i]))]), , drop = F]; if (nrow(ss) > 1L) {ss_shared <- ss; ss_shared[
      ss_shared > 0.1] <- 1; shared_vec <- which(colSums(ss_shared) >= 2)
      word <- "variant"; vct <- length(shared_vec); if (vct == 1L) {word <- word
      } else {word <- paste0(word, "s")}; message("--", round(sum(colSums(ss)[
        shared_vec]) / sum(ss), 4) * 100, " % of reads were in ", vct, 
        " shared ", word, "\namong ", nrow(ss), " samples within batch ", r[i],
        ".--")}}}
  if (lr > 1) {combns <- combn(r, 2); for (i in 1:ncol(combns)) {
    j <- as.character(combn(r, 2)[1, i]); k <- as.character(combns[2, i])
    ssj <- s[unlist(FL[which(nfl == j)]), , drop = F]; ssk <- s[unlist(FL[which(
      nfl == k)]), , drop = F]; ssj_shared <- ssj
    ssj_shared[ssj_shared > 0.1] <- 1; ssk_shared <- ssk; ssk_shared[
      ssk_shared > 0.1] <- 1; shared_vec <- which(colSums(
        ssk_shared) >= 1 & colSums(ssj_shared) >= 1); word <- "variant"
    vct <- length(shared_vec); if (vct == 1L) {word <- word} else {
      word <- paste0(word, "s")}; message("--", round(sum(colSums(ssj)[
      shared_vec]) / sum(ssj), 4) * 100, "% of run ", j, " reads and ", round(
        sum(colSums(ssk)[shared_vec]) / sum(ssk), 4) * 100, "% of run ", k, "r",
        "eads are in ", vct, " shared ", word, "\namong batches.--")}
  } #combining repeats.
  if ("repeat." %in% colnames(m)) {suppressPackageStartupMessages(library(
    dplyr)); repeats_by_lib <- m %>% group_by(repeat.) %>%
      select(repeat., library) %>% group_nest()
    cat("Summing counts for repeats!\n"); sm <- matrix(nrow = nrow(
      repeats_by_lib), ncol = ncol(s), dimnames = list(
        repeats_by_lib$repeat., colnames(s))); tm <- matrix(nrow = nrow(
          repeats_by_lib), ncol = ncol(tracks), dimnames = list(
            repeats_by_lib$repeat., colnames(tracks)))
    for (i in 1:nrow(repeats_by_lib)) {
      if (length(repeats_by_lib$data[[i]][[1]]) == 1) {
        sm[rownames(sm)[i], ] <- s[repeats_by_lib$data[[i]][[1]][1], ]
        tm[rownames(tm)[i], ] <- tracks[repeats_by_lib$data[[i]][[1]][1], ]
        } else {
          sm[rownames(sm)[i], ] <- colSums(s[repeats_by_lib$data[[i]][[1]], ])
          tm[rownames(tm)[i], ] <- colSums(tracks[
            repeats_by_lib$data[[i]][[1]], ])}} #error check below
    if (sum(s) != sum(sm)) {stop(
      "an unexpected error occurred while merging repeats for the sequence tab",
      "le")}; if (sum(tracks) != sum(tm)) {stop(
        "an unexpected error occurred while merging repeats for tracking reads")
      }; pkg <- "package:dplyr"; detach(pkg, character.only = T)
    } else {sm <- s; tm <- tracks}
  #10k filtering. run only AFTER combining repeats. 
  #calculate shared variants among repeats. 
  #but do so internally to repeat if(). need a test case.
  if (any(rowSums(sm) < 10000)) {rr <- rs10k(
      sm); sm <- sm[-rr, ]; tm <- tm[-rr, ]; if (any(colSums(sm) == 0)) {
        rc <- rt0(sm); sm <- sm[, -rc]}; name_compare2(s = sm, t = tm)}
  #save sm and write tm
  mode(sm) <- "integer"; saveRDS(tm, "results/seqtab/tm.rds"); return(sm)}
#e.g., mat_sum_repeat()

#Bn function takes PEIV13 sequence table matrix default = st.mr and 
#identifies and removes chimeras via removeBimeraDenovo(method = "consensus")
#bn reports decimals of total counts for schim, relative to st.mr
Bn <- function(s = st.mr) {seed <- 100; multithread_val = T
if (.Platform$OS.type == "windows") multithread_val = F; set.seed(seed)
cat("Removing chimeric sequences with the dada2 consensus method!\n")
beg_time <- Sys.time(); schim <- removeBimeraDenovo(
  s, method = "consensus", multithread = multithread_val); tm <- readRDS(
    "results/seqtab/tm.rds"); tm <- cbind(tm, rowSums(schim)); colnames(
      tm)[7] <- "nonchim"; if (any(rowSums(schim) < 10000)) {
        rr <- rs10k(schim); schim <- schim[-rr, ]; tm <- tm[-rr, ]
  if (any(colSums(schim) == 0)) {rc <- rt0(schim); schim <- schim[, -rc]}
  #made this explicit like in mat_sum_repeat on 11Apr25
  name_compare2(s = schim, t = tm)}; saveRDS(tm, "results/seqtab/tm.rds")
end_time_h <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours")
cat("Chimera removal completed in", paste0(elapsed_time, "!\n")); message(
  "--", round(sum(schim)/sum(st.mr) * 100, 2), "% of merged sequence reads a",
  "nd ", round(dim(schim)[2] / dim(st.mr)[2] * 100, 2), "% of merged sequence ",
  "variants remain\nin the sequence table with chimeras removed.  Most of your",
  " reads should remain\nafter chimera removal (it is not uncommon for a major",
  "ity of sequence variants to\nbe removed though).  If most of your reads wer",
  "e removed as chimeric, upstream\nprocessing may need to be revisited.  In a",
  "lmost all cases this is caused by\nunremoved primers interfering with chime",
  "ra identification.  PEIV13 controlfile\nprimer removing parameters for R1 a",
  "nd R2 raw reads are forwardtrimLeft and\nreversetrimLeft respectively.--")
return(schim)}
#e.g., Bn()

#aT function takes PEIV13 sequence table matrix default = st.chim and a PEIV13
#controlfile default uservars. returns a taxonomy sequence table.
aT <- function(s = st.chim, uv = uservars) {d <- pattern2string("databasefile")
trc <- as.logical(pattern2string("tryRC")); seed <- 100; multithread_val = T
if (.Platform$OS.type == "windows") multithread_val = F; set.seed(seed)
beg_time <- Sys.time(); cat("Taxonomy assignment started!\n")
taxa <- assignTaxonomy(s, d, multithread = multithread_val, tryRC = trc)
end_time_h <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours"); cat(
    "Taxonomy assignment completed in", paste0(elapsed_time, "!\nKingdom"),
    "summary:"); print(table(taxa[, 1], useNA = "ifany")); return(taxa)}
#e.g., aT()

#rmBac function takes PEIV13 sequence table (default st.chim) and otu table
#default(taxa) and returns a new bacteria only sequence table
rmBac <- function(s = st.chim, t = taxa) {suppressPackageStartupMessages(
  library(data.table)) 
  #https://stackoverflow.com/questions/58327399/negate-like-in-data-table-package-r/58327444#58327444
  #for negating. t[!, 1] %ilike% "Bacteria" failed to load via source() with 
  #error. a notilike function worked
  `%notilike%` <- Negate(`%ilike%`); message("--Removing ", length(which(t[
    , 1] %notilike% "Bacteria")), " non Bacterial taxa.--"); sb <- s[, which(
      t[, 1] %ilike% "Bacteria"), drop = F]; pkg <- "package:data.table"
  detach(pkg, character.only = T); return(sb)}
#e.g., st.bac <- rmBac()

#rmBac2 function takes a PEIV13 sequence table (default st.bac) and otu table
#default(taxa) and returns a new bacteria only otu table. final tm update
rmBac2 <- function(sbac = st.bac, schim = st.chim, tax = taxa) {
  suppressPackageStartupMessages(library(data.table)); taxbac <- tax[which(tax[
    , 1] %ilike% "Bacteria"), ]; trm <- readRDS("results/seqtab/tm.rds")
    name_compare3(sbac, trm, taxbac); if (any(rowSums(sbac) < 10000)) {
      rr <- rs10k(sbac); sbac <- sbac[-rr, ]; trm <- trm[-rr, ]
      if (any(colSums(sbac) == 0)) {rc <- rt0(sbac); sbac <- sbac[, -rc]
      taxbac <- taxbac[-rc, ]}; name_compare3(sbac, trm, taxbac)}
    trm <- cbind(trm, rowSums(sbac)); colnames(trm)[8] <- "bac"; saveRDS(
      trm, "results/seqtab/tm.rds"); taxareadloss <- (sum(schim) - sum(
        sbac)) / sum(schim)
    if (taxareadloss > 0.5) warning(
      "Your reads do not seem to be appropriately assigned!\nCheck for lots of",
      " 16S sequences assigned as Kingdom Eukaryota NA NA NA rather\nthan Bact",
      "eria.  Your reads may be in the opposite orientation as the reference\n",
      "database.  Use a different PEIV13 parameter tryRC and see if this fixes",
      "the\nassignments."); message("--", round(
        taxareadloss, 4) * 100, "% of reads removed as non Bacteria.--")
    cat("Saving final TrackReads file"); saveRDS(
      trm, "results/TrackReads_rep.rds"); pkg <- "package:data.table"
    detach(pkg, character.only = T); return(taxbac)}
#e.g., rmBac2()
							    
#make_ps function takes PEIV13 metadata, sequence table, and taxa and saves a
#phyloseq rds
make_ps <- function(m = md, s = st.bac, t = taxa.bac) {cat(
  "Making phyloseq object from metadata, sequence table and taxonomy table!\n")
  if ("repeat." %in% colnames(m)) {samdf <- m[, -which(colnames(
    m) == "library")]; samdf <- samdf[!duplicated(samdf$repeat.), ]; rownames(
      samdf) <- samdf$repeat.} else {samdf <- m
  rownames(samdf) <- samdf$library}; library(phyloseq); ps <- phyloseq(
    otu_table(s, taxa_are_rows = F), sample_data(samdf), tax_table(t))
  saveRDS(ps, "results/ps.rds"); pkg <- "package:phyloseq"
  detach(pkg, character.only = T)}
#e.g., make_ps()

#rm_peiv13_int function cleans up intermediate files.
rm_peiv13_int <- function() {f <- list.files("results/seqtab", "rds")
cat("Removing intermediate files!\n"); for (i in 1:length(f)) unlink(file.path(
  "results/seqtab", f[i]), force = T); f <- list.files("results/seqtab", "rds")
if (length(f) > 0) stop("failed to remove intermediate files"); unlink(
  "results/filt", recursive = T, force = T); if (dir.exists(
  "results/filt")) stop("failed to remove intermediate files"); file.rename(
    "results/seqtab", "results/varLenHistograms"); if (dir.exists(
      "results/seqtab")) stop("failed too create histograms directory")
f <- list.files("results", "TrackReads", full.names = T); if (nrow(readRDS(
  f[1])) == nrow(readRDS(f[2]))) {unlink(f[1]); if (file.exists(f[1])) stop(
  "failed to remove intermediate files")}}
#e.g., rm_peiv13_int()
