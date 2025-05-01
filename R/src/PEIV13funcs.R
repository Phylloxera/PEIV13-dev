#PEIV13funcs.R being developed to run a single Rscript,
#PEIV13_ctrl2counts_rscript.r

#getc takes controlfile path as character string and returns a dataframe
#from the file at that path
getc <- function(cp) {if (!file.exists(cp)) stop(cp, " not found") #this is
  #an informative AND exiting error
  uv <- read.csv(cp, header = F); if (!is.list(uv)) stop(
    cp, " not properly coerced to a list"); return(uv)}
#e.g. uservars <- getcf("ctrl_file_path")

#pattern2string(pattern, uv = uservars) takes a patern string and a peiv13
#controlfile uservars and returns the user supplied value from column 2
pattern2string <- function(pattern, uv = uservars) {
  out <- uv[which(uv[, 1] == pattern), 2]}
#e.g. pattern2string("metadatafile")

#checkfs takes uv (default uservars) and performs confirms existence of inputs
checkfs <- function(uv = uservars) {mf <- pattern2string("metadatafile")
if (!file.exists(mf)) stop(mf, " not found"); rawfol <- pattern2string(
  "inputfolder"); if (!dir.exists(rawfol)) stop(rawfol, " not found")
dbf <- pattern2string("databasefile"); if (!file.exists(dbf)) stop(
  dbf, " not found")}
  #never reproduced or found source of env: Rscript\r: No such file or directory
#e.g. checkfs()

#checkenv stops if a required R package is not found in the R installation
checkenv <- function() {if (!("dada2" %in% installed.packages())) stop(
  "dada2 is reqired for PEIV13 but is not an installed package in your R insta",
  "llation"); if (!("data.table" %in% installed.packages())) stop(
  "data.table is reqired for PEIV13 but is not an installed package in your R ",
  "installation"); if (!("dplyr" %in% installed.packages())) stop(
  "dplyr is reqired for PEIV13 but is not an installed package in your R insta",
  "llation"); if (!("phyloseq" %in% installed.packages())) stop(
  "phyloseq is reqired for PEIV13 but is not an installed package in your R in",
  "stallation")}
#e.g. checkenv()

#check_mdat function takes PEIV13 metadata data.frame default = md and prints 
#test results
check_mdat <- function (m = md) {if ("library" %in% colnames(
  m) && "run" %in% colnames(m)) {cat("required columns found in metadata\n")
} else {stop("required columns not found in metadata")}
  #MUST ENFORCE UNIQUENESS OF library COLUMN FOR FileList TO WORK.
  if (nrow(m) > length(unique(m[, "library"]))) {stop(
    "library column MUST be unique for PEIV13")
	#set conditions for (optional) repeat column
	} else if ("repeat." %in% colnames(m) && length(unique(m[
      , "repeat."])) < nrow(m)) {cat(length(unique(m[
        , "repeat."])), "samples from", nrow(
          m), "rows accounting for repeats in metadata.\n"); message(
            "per-variant sequence counts will be summed for repeats")
    } else {cat("no technical replicates in metadata\n")}}
#e.g. check_mdat()

#dada2var_from_ctrl function for PEIV13 extracts a named integer based on 
#variable labels in a PEIV13 controlfile dataframe. ordering by names in 
#as.character() order more robust
dada2var_from_ctrl <- function(pattern, uv = uservars) {namedint <- as.integer(
  uv[which(grepl(pattern, uv[, 1])), 2]); names(namedint) <- sub(
    pattern, "", uv[which(grepl(pattern, uv[, 1])), 1])
  namedint <- namedint[order(names(namedint))]; return(namedint)}
#e.g. forwarderror <- dada2var_from_ctrl("forwarderror")

#subset_md function takes md and fnFs and returns a subsetted md
subset_md <- function(m = md, fF = fnFs, fR = fnRs) {unused <- unused_patterns(
  fF, m[, "library"]); if (length(unused) > 0) m <- m[-unused, ]; libids <- m[
    , "library"]
    #stop the run in length(fnFs) != length(fnRs) or if length(libids) is zero
  if (length(fF) != length(fR)) {stop(
    "equal numbers of forward and reverse raw read files needed for paired-end",
    " protocol")}; if (length(libids) == 0L) {stop(
      "no samples to process. check compatibility between library and inputfol",
      "der")}; return(m)}

#varMessages function lacks arguments. just posts messages
varMessages <- function() {
  message(
  "Reading maxerror parameters from controlfile.\nIf you want to speed up down",
  "stream computation, consider tightening (decreasing)\nthese.\nIf too few re",
  "ads are passing the filter, consider relaxing these, perhaps\nespecially re",
  "versemaxerror and reducing the trunc parameters to remove low\nquality tail",
  "s.\nRemember though, when choosing trunc you must maintain overlap after tr",
  "uncation\nin order to merge R1 and R2 pairs.")
  message(
  "Reading forwardtrunc and reversetrunc parameters from controlfile.\nYour re",
  "ads must still overlap after truncation in order to merge them later!\nWhen",
  " using a a less-overlapping primer set, like V1-V3, these parameters must b",
  "e\nlarge enough to maintain 20 + biological-length-variation nucleotides of",
  " overlap\nbetween them.")
  message(
  "Reading trimLeft parameters from controlfile.\nConstant length primers at t",
  "he start of your reads is a common scenario.\nThe forward and reverse param",
  "eters will be applied to your R1 and R2 files\nrespectively.\nIf your forwa",
  "rd primers are in R2, set tryRC to TRUE for assignTaxonomy.\nPlease double-",
  "check that your primers have been removed by dada2!")}

#multieq function for PEIV13 takes 3 or more R objects of the same class and 
#tests identity returning TRUE or FALSE
multieq <- function(A, B, C, ...) {classes <- c(class(A), class(B), class(C))
input_list <- list(...); elipsecount <- length(input_list); if (
  elipsecount > 0) {for (i in 1:elipsecount) {classes <- c(classes, class(
    input_list[[i]]))}}; if (elipsecount > 0) {for (i in 1:elipsecount) {assign(
      LETTERS[3 + i], input_list[[i]])}}; tcount <- elipsecount + 3; if (all(
        sapply(list(classes[1:(tcount - 1)]), function(
    x) x == classes[tcount])) == F) {stop(
      "arguments do not share the same class")}; all(sapply(mget(LETTERS[1:(
        tcount - 1)]), function(x) x == get(LETTERS[tcount])))}
#e.g. multieq(names(forwarderror), names(reverseerror), as.character(runs))

#check_taxdb function takes uservars. includes tryRC check
check_taxdb <- function(uv = uservars) {message(
  "checking Taxonomy database files/parameters"); trc <- as.logical(
    pattern2string("tryRC")); if (is.na(trc)) stop(
      "tryRC needs to be TRUE or FALSE"); d <- pattern2string("databasefile")
      #what happens with below when multiple trailing ';'
      L1 <- readLines(d, n = 1); if (length(unlist(strsplit(
        L1, ';'))) != 6) stop(
          d, " does not appear to have 6 taxonomic ranks as expected for Genus",
  "-level identification"); if (substring(L1,1,1) != ">") stop(
    d, " is not a FASTA file")}
#e.g. check_taxdb()

#setup_peiv13_results function takes uservars.
#cat() and message() syntax seem right. generally "ing" indicates message()
setup_peiv13_results <- function(uv = uservars) {ow <- as.logical(
  pattern2string("overwrite")); if (is.na(ow)) stop(
    "overwrite needs to be TRUE or FALSE")
  #if results exists and overwrite is FALSE exit
if (dir.exists("results") && ow == F) stop(
  "the PEIV13 results directory ", getwd(),
  "/results will not be overwritten since overwrite is FALSE")
#if results exists and overwrite is TRUE delete and recreate else create
#force = T needed on windows one drive
if (dir.exists("results") && ow == T) {
  message("overwriting ", getwd(), "/results"); unlink(
    "results", recursive = T, force = T); if (dir.exists("results")) {
      stop("the PEIV13 results directory ", getwd(), "/results could not be ov",
           "erwritten")}; dir.create("results"); dir.create("results/filt")
  dir.create("results/seqtab"); dir.create("results/plotErrors")} else {
    message("creating ", getwd(), "/results"); dir.create("results")
    dir.create("results/filt"); dir.create("results/seqtab"); dir.create(
      "results/plotErrors")}}
#e.g. setup_peiv13_results()

#check_peiv13_results function
check_peiv13_results <- function(uv = uservars) {ow <- as.logical(
  pattern2string("overwrite")); if (is.na(ow)) stop(
    "overwrite needs to be TRUE or FALSE")
#if results does not exist exit
if (!dir.exists("results")) stop(
  "the PEIV13 results directory ", getwd(), "/results does not exist")
#if seqtab does not exist exit
if (!dir.exists("results/seqtab")) {stop(
  "the PEIV13 results/seqtab directory ", getwd(), 
  "/results/seqtab does not exist")}}
#e.g. check_peiv13_results()

#the unused_patterns function takes a long vector (subject vector) and a
#pattern vector (query vector) (no defaults) and returns the indices of the 
#long/subject, which do not match any of the patterns
unused_patterns <- function(longvec, pvec) {
  up <- integer(0); suppressPackageStartupMessages(library(data.table))
  for (i in 1:length(pvec)) {if (length(which(
    longvec %ilike% pvec[i])) == 0L) {up <- c(up, i)}}
  pkg <- "package:data.table"; detach(pkg, character.only = T); return(up)}
#e.g. unused <- unused_patterns(fnFs, md[, "library"])

#the batch_ilike_subset_vec function takes a subset (subject) vector and a 
#pattern vector (query) (no defaults) and returns the subject with entries
#removed that are not %ilike% (package data.table) any of the patterns.
batch_ilike_subset_vec <- function(ssvector, pvector) {
  index <- integer(); suppressPackageStartupMessages(library(data.table))
  for (i in 1:length(pvector)) {index <- c(index, which(
    ssvector %ilike% pvector[i]))}; pkg <- "package:data.table"; detach(
      pkg, character.only = T); out <- ssvector[index]}
#e.g. fnRs <- unique(batch_ilike_subset_vec(fnRs, md[, "library"]))

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
      "user supplied dada2 trunc and error variables match\n")
    } else {stop("check that error and trunc controlfile labels match your met",
             "adata run column")}
  cat("second check of metadata after comparing file names against library\n")
  #PEIV13 message
  check_mdat(); message(
    "PEIV13 will now run the dada2 pipeline for paired-end bacterial 16s varia", 
    "ble\nregion 1-3 data.\nError plots will be in results/plotErrors/<directi",
    "on>_<batch>.pdf showing the\nerror rates for each possible transition (A→",
    "C, A→G, ...).\nPoints are the observed error rates for each consensus qua",
    "lity score.\nThe black line shows the estimated error rates after converg",
    "ence of the machine-\nlearning algorithm.\nThe red line shows the error r",
    "ates expected under the nominal definition of the\nQ-score.\nIf the estim",
    "ated error rates (black line) are a good fit to the observed rates\n(poin",
    "ts) and the error rates drop with increased quality as expected you can\n",
    "have confidence in the error model.\nExample plots that might justify enf",
    "orcing monotonicity in the fitted error\nmodel are here: https://github.c",
    "om/benjjneb/dada2/issues/1156.\nDuring merging, sequences are only output",
    " if the forward and reverse reads\noverlap by at least 12 bases, and are ",
    "identical to each other in the overlap\nregion.")}
#e.g. PEIV13varcheck()

#raw2pQP requires metadata and controlfile. makes a folder called qualProf if
#not present and writes plots to pdf, labeling according to library and run.
raw2pQP <- function(fnF = fnFs, fnR = fnRs, FL = FileList) {r <- as.character(r)
res <- "qualProf"; if (!dir.exists(res)) dir.create(res); if (!dir.exists(
  res)) stop("PEIV13 was unable to create a quality plots folder"); s <- unlist(
    FL[r]); Fs <- unique(batch_ilike_subset_vec(fnF, s)); Rs <- unique(
      batch_ilike_subset_vec(fnR, s))
#plot quality profiles
cat("printing quality profile plots for run", r, "\n"); beg_time <- Sys.time()
if (length(Fs) == 1) {pqF <- suppressWarnings(plotQualityProfile(Fs))
warning("a single R1 file is available to plot a quality profile")
} else if (length(Fs) > 1) {pqF <- suppressWarnings(plotQualityProfile(Fs[1:2]))
} else {stop("no R1 libraries were found from run ", r)}; pdf(file.path(
  res, paste0(r, "_R1.pdf"))); suppressWarnings(print(pqF))
garbage <- dev.off()
if (length(Rs) == 1) {pqR <- suppressWarnings(plotQualityProfile(Rs))
warning("a single R2 file is available to plot a quality profile")
} else if (length(Rs) > 1) {pqR <- suppressWarnings(plotQualityProfile(Rs[1:2]))
} else {stop("no R2 libraries were found from run ", r)}; pdf(file.path(
  res, paste0(r, "_R2.pdf"))); suppressWarnings(print(pqR))
garbage <- dev.off()
end_time_m <- Sys.time(); elapsed_time <- paste0(round(
    as.numeric(difftime(time1 = end_time_m, time2 = beg_time, units = "mins")), 
    2), " Minutes"); cat("the first 2 read quality profiles for run", r, 
      "completed in", elapsed_time, "\n")}

#raw2seqtab requires 10 default defined variables. and writes a dada2 sequence
#table as a rds file in file.path(getwd(), "results/seqtab"). runs 
#filterAndTrim, learnErrors, dada, mergePairs, and makeSequenceTable functions
#from the dada2 r package for a single run of raw data (only miseq paired end 
#tested)
raw2seqtab <- function(
    fe = forwardmaxerror, re = reversemaxerror, ft = forwardtrunc, 
    rt = reversetrunc, ftl = forwardtrimLeft, rtl = reversetrimLeft, fnF = fnFs,
    fnR = fnRs, FL = FileList) {
  #set up
  seed <- 100; multithread_val = T; maxN_val <- 0; truncQ_val <- 2
  r <- as.character(r); maxEE_val <- c(fe[r], re[r]); truncspec <- c(
    ft[r], rt[r]); res <- "results"; filt <- file.path(res, "filt")
  #for batches/runs put in different folders
  s <- unlist(FL[r]); filtFs <- file.path(filt, r, paste0(
    s, "_F_filt.fastq.gz")); filtRs <- file.path(filt, r, paste0(
      s, "_R_filt.fastq.gz"))
  #batch_ilike_subset not necessary but some subset is.
  Fs <- unique(batch_ilike_subset_vec(fnF, s))
  Rs <- unique(batch_ilike_subset_vec(fnR, s))
  if (.Platform$OS.type == "windows") multithread_val = F; set.seed(seed)
  #filterAndTrim()
  beg_time <- Sys.time(); cat("filterAndTrim for run", r, "\n")
  out <- filterAndTrim(
    Fs, filtFs, Rs, filtRs, truncLen = truncspec, trimLeft = c(ftl, rtl), 
    maxN = maxN_val, maxEE = maxEE_val, truncQ = truncQ_val, rm.phix = T, 
    compress = T, multithread = multithread_val); end_time_h <- Sys.time()
  elapsed_time <- paste0(round(as.numeric(difftime(
    time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours")
  cat("filterAndTrim for run", r, "completed in", elapsed_time, "\n"); message(
    "learning errors for forward libraries"); beg_time <- Sys.time()
  #learnErrors()
  errF <- learnErrors(
    filtFs, nbases = 1e8, multithread = multithread_val, randomize = T)
  pEF <- suppressWarnings(plotErrors(errF, nominalQ = T)); pdf(paste0(
    "results/plotErrors/F_", r, ".pdf")); suppressWarnings(print(pEF))
  garbage <- dev.off(); message("learning errors for reverse libraries")
  errR <- learnErrors(
    filtRs, nbases = 1e8, multithread = multithread_val, randomize = T)
  pER <- suppressWarnings(plotErrors(errR, nominalQ = T))
  pdf(paste0("results/plotErrors/R_", r, ".pdf")); suppressWarnings(print(pER))
  garbage <- dev.off(); end_time_h <- Sys.time(); elapsed_time <- paste0(round(
    as.numeric(difftime(time1 = end_time_h, time2 = beg_time, units = "hours")), 
    2), " Hours"); cat("learnErrors for run", r, "completed in", elapsed_time, 
                       "\n"); message("running dada for forward libraries")
  beg_time <- Sys.time()
  #dada()
  dadaFs <- dada(filtFs, err = errF, multithread = multithread_val, verbose = F)
  message("running dada for reverse libraries")
  dadaRs <- dada(filtRs, err = errR, multithread = multithread_val, verbose = F)
  end_time_h <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
    time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours"); cat(
      "dada for run", r, "completed in", elapsed_time, "\n"); message(
        "Merging pairs"); beg_time <- Sys.time() 
  #mergePairs()
  mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
  #additional || nrow(mergers[[1]]) == 0L category needed windows.
  if(isTRUE(nrow(mergers) == 0L || nrow(mergers[[1]]) == 0L)) { 
    #report the column 1 reduction in out
    qualfiltper <- round((sum(out[, 1]) - sum(out[, 2])) / sum(
      out[, 1]) * 100, 2)
    stop(
      "There are no remaining sequences in run\n", r, " after mergePairs.\n",
      qualfiltper, "% of reads from run ", r, " were removed by quality filter",
      "ing\nIf this is >>50, your PEIV13 controlfile forward and reversetrunc",
      "parameters\nmight be too small.\nLikewise, your maxerror parameters.\n",
      "Remember that when using a a less-overlapping primer set, like V1-V3,\n",
      "trunc parameters must be large enough to maintain 20 + biological-lengt",
      "h-\nvariation nucleotides of overlap between them.")}
  end_time_s <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
    time1 = end_time_s, time2 = beg_time, units = "secs")), 2), " Seconds")
  cat("mergePairs for run", r, "completed in", elapsed_time, "\n"); message(
      "making sequence table and initiating read-loss tracking")
  #makeSequenceTable()
  beg_time <- Sys.time(); seqtab <- makeSequenceTable(mergers)
  ct_filt400 <- sum(nchar(colnames(seqtab)) < 400); if (ct_filt400 > 0) {
    message("Removing ", ct_filt400, " variants shorter than 400.\nSee doi:10.",
            "1371/journal.pone.0129174")}
  seqtab <- seqtab[, nchar(colnames(seqtab)) > 399, drop = F]; hist <- hist(
    nchar(getSequences(seqtab)), plot = F); hist[["xname"]] <- paste0(
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
  remain400 <- round(sum(track[, 6])/sum(track[, 1]) * 100, 2); cat(
    "For run", r, paste0(
      fremain, "%, ", dremain, "%, ", mremain, "%, ", remain400, "%"), "of", 
    "input reads remain after\nquality filtering, denoising, merging, and 400",
    "length filtering respectively.\nOutside of filtering, there should be no",
    "step in which a majority of reads are lost.\n")
  qualfiltcomp <- (sum(track[, 1]) - sum(track[, 2])) / sum(track[, 1])
  if (qualfiltcomp > 0.5) warning(
    "A majority of reads from run ", r, " were removed by quality filtering!\n",
    "You might consider relaxing (increasing) the PEIV13 controlfile forward a",
    "nd reversemaxerror parameters used in filterAndTrim, however, your raw da",
    "ta might have too many errors if a maxerror of 5 is still insufficient.")  
  denoisecomp <- (sum(track[, 2]) - minDenois) / sum(track[, 2])
  mergecomp <- (minDenois - sum(track[, 5])) / minDenois
  if (mergecomp > 0.5) warning(
    "A majority of quality filtered reads from run ", r, " were removed by den",
    "oising!\nYou might consider relaxing (increasing) the PEIV13 controlfile ",
    "forward and reversemaxerror parameters used in filterAndTrim, however, yo",
    "ur raw data might have too many errors if a maxerror of 5 is still insuff",
    "icient.")
  if (mergecomp > 0.5) warning(
    "A majority of quality filtered and denoised reads from run ", r, " were r",
    "emoved by merging!\nYou should revisit the PEIV13 controlfile forwardtrun",
    "c and reversetrunc parameters used in filterAndTrim and make sure that th",
    "e truncated reads span your amplicon.")
  filter400comp <- (sum(track[, 5]) - sum(track[, 6])) / sum(track[, 5])
  if (filter400comp > 0.5) warning(
    "User specified raw reads are robably not V1V3 data!\nA majority of reads ",
    "from run ", r, " were removed by 400 length filtering!")
  rownames(track) <- s; rownames(seqtab) <- s; saveRDS(seqtab, paste0(
    "results/seqtab/seqtab", r, ".rds")); saveRDS(track, paste0(
      "results/seqtab/track", r, ".rds")); end_time_s <- Sys.time()
  elapsed_time <- paste0(round(as.numeric(difftime(
    time1 = end_time_s, time2 = beg_time, units = "secs")), 2), " Seconds")
  cat("sequence table and initial read-loss tracking for run", r, "completed",
      "in", elapsed_time, "\n")}
#e.g. raw2seqtab()

#i think works but should be tested when length(runs) == 1.
multiseqtab <- function(r = runs) {seqtabs <- character()
trackfiles <- character(); for (run in r) {run <- as.character(run); seqtabs[
  run] <- c(paste0("results/seqtab/seqtab", run, ".rds")); trackfiles[run] <- c(
    paste0("results/seqtab/track", run, ".rds"))}; message(
      "reading track reads"); tracks <- do.call(rbind, lapply(
        trackfiles, readRDS)); saveRDS(tracks, "results/seqtab/trackmult.rds")
if (length(r) == 1L) {s <- readRDS(seqtabs)} else {seed <- 100; message(
  "merging Sequence Tables"); set.seed(seed); s <- mergeSequenceTables(
    tables = seqtabs)}; return(s)}
#e.g. multiseqtab()

#rs10k function takes a sequence count matrix s (no default) and returns an 
#index of samples to remove from the sequence count matrix AND track reads 
#matrix. was not being applied correctly. 
#MUST only be run if any(rowSums(s) < 10000)
rs10k <- function(s) {rs <- which(rowSums(s) < 10000); warning(
  "Samples ", names(rs), " will be removed due to sequence depth < 10,000.\nSe",
  "e doi:10.1093/jas/skab346."); return(rs)}
#e.g. rs10k(st.bac)

#rsSpecify function for debugging as a drop-in replacement for rs10k.
#Takes a sequence count matrix s (no default) and an integer (no default) and 
#returns an index of samples to remove
rsSpecify <- function(s, i) {i <- as.integer(i); if (!(is.integer(i))) stop(
  "i is not an integer"); rs <- which(rowSums(s) < i); warning(
    "Samples ", names(rs), " will be removed due to arbitrary sequence depth <",
    " ", i); return(rs)}

#rt0 function takes a sequence count matrix s (no default) and returns an index 
#of taxa (sequences) to remove from the sequence count matrix AND taxa 
#identification/classification matrix. 
#MUST be run after rs10k has been run AND if any(colSums(s) == 0)
rt0 <- function(s) {if (nrow(s) == 0){stop(
  "No samples remain in the sequence table after filtering rows with a sequenc",
  "e depth < 10,000.")} else if (nrow(s) > 1L) {rt <- which(colSums(s) == 0)
  } else {rt <- which(s[1, ] == 0)}; warning(length(
    rt), " variants will also be removed as they were only found in the sample",
    "s with sequence depth < 10,000."); return(rt)}
#e.g. rt0(st.bac)

#name_compare2 function takes a sequence count matrix s (no default) and a
#track changes object t (default tm) and performs 2 dim names identity checks
#1lib case needs testing
name_compare2 <- function(s, t = tm) {
  mes <- "a problem occurred when removing non Bacterial sequences"
  if (nrow(s) > 1L) {if (any(colSums(s) == 0)) stop(mes)
    } else {if (any(s[1, ] == 0)) stop(mes)}
  if (all(rownames(t) != rownames(s))) stop(mes); return(mes)}
#e.g. name_compare2(st.chim)

#name_compare3 function takes a sequence count matrix s2 (no default) and a
#track changes object t2 (default tm) and a dada2 taxa object tax
#(default taxa.bac) and performs names_compare2 + 1 more dim name identity check
name_compare3 <- function(s2, t2 = trm, tax = taxbac) {mes <- name_compare2(
  s = s2, t = t2); if (all(rownames(tax) != colnames(s2))) stop(mes)}
#e.g. name_compare3(st.bac)

#mat_sum_repeat function takes PEIV13 metadata data.frame default = md and 
#mergeSequenceTable() matrix result default = st.all and merges repeats via
#colSums(). had been getting colSums() errors when nrow was 1L, but then,
#stopped. confused.
mat_sum_repeat <- function (m = md, r = runs, FL = FileList, s = st.all) {
  libvec <- unlist(m["library"]); if (all(rownames(s) != libvec)) stop(
    "unable to match seqtab to metadata")
  tracks <- readRDS("results/seqtab/trackmult.rds"); if (all(rownames(
    tracks) != libvec)) stop("unable to match read-loss tracking to metadata")
  write.csv(tracks, "results/TrackReads_lib.csv", quote = F)
  #above: keep tracks-by-library for reference
  #below: shared variants within run.
  nfl <- names(FL); lr <- length(r); if (nrow(s) > 1L) {cat(
    "calculating shared variants. few shared variants can indicate technical", 
    "problems such as inexact region amplification, incomplete primer removal,", 
    "heterogeneity spacers and batch effects", fill = T); for (i in 1:lr) {
      ss <- s[unlist(FL[which(nfl == as.character(r[i]))]), , drop = F]
      if (nrow(ss) > 1L) {ss_shared <- ss; ss_shared[ss_shared > 0.1] <- 1
      shared_vec <- which(colSums(ss_shared) >= 2); word <- "variant"
      vct <- length(shared_vec); if (vct == 1L) {word <- word
      } else {word <- paste0(word, "s")}; cat(paste0(round(sum(colSums(ss)[
        shared_vec]) / sum(ss), 4) * 100, "%"), "of reads were in", vct, 
        "shared", word, "among", nrow(ss), "samples within batch", r[i],
        fill = T)}}} #shared variants across runs.
  if (lr > 1) {combns <- combn(r, 2); for (i in 1:ncol(combns)) {
    j <- as.character(combn(r, 2)[1, i]); k <- as.character(combns[2, i])
    ssj <- s[unlist(FL[which(nfl == j)]), , drop = F]; ssk <- s[unlist(FL[which(
      nfl == k)]), , drop = F]; ssj_shared <- ssj
    ssj_shared[ssj_shared > 0.1] <- 1; ssk_shared <- ssk; ssk_shared[
      ssk_shared > 0.1] <- 1; shared_vec <- which(colSums(
        ssk_shared) >= 1 & colSums(ssj_shared) >= 1); word <- "variant"
    vct <- length(shared_vec); if (vct == 1L) {word <- word} else {
      word <- paste0(word, "s")}; cat(paste0(round(sum(colSums(ssj)[
      shared_vec]) / sum(ssj), 4) * 100, "%"), "of run", j, "reads and", paste0(
        round(sum(colSums(ssk)[shared_vec]) / sum(ssk), 4) * 100, "%"), "of", 
      "run", k, "reads are in", vct, "shared", word, "among batches", fill = T)}
  } #combining repeats.
  if ("repeat." %in% colnames(m)) {suppressPackageStartupMessages(library(
    dplyr)); repeats_by_lib <- m %>% group_by(repeat.) %>%
      select(repeat., library) %>% group_nest()
    message("summing counts for repeats"); sm <- matrix(nrow = nrow(
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
#e.g. mat_sum_repeat()

#Bn function takes PEIV13 sequence table matrix default = st.mr and 
#identifies and removes chimeras via removeBimeraDenovo(method = "consensus")
#bn reports decimals of total counts for schim, relative to st.mr
Bn <- function(s = st.mr) {seed <- 100; multithread_val = T
if (.Platform$OS.type == "windows") multithread_val = F; set.seed(seed)
message("Removing chimeric sequences with the dada2 consensus method")
beg_time <- Sys.time(); cat("removeBimeraDenovo started\n")
schim <- removeBimeraDenovo(
  s, method = "consensus", multithread = multithread_val); tm <- readRDS(
    "results/seqtab/tm.rds"); tm <- cbind(tm, rowSums(schim)); colnames(
      tm)[7] <- "nonchim"; if (any(rowSums(schim) < 10000)) {
        rr <- rs10k(schim); schim <- schim[-rr, ]; tm <- tm[-rr, ]
  if (any(colSums(schim) == 0)) {rc <- rt0(schim); schim <- schim[, -rc]}
  #made this explicit like in mat_sum_repeat on 11Apr25
        name_compare2(s = schim, t = tm)}; saveRDS(tm, "results/seqtab/tm.rds")
end_time_h <- Sys.time(); elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours")
cat("removeBimeraDenovo completed in", paste0(elapsed_time, ".\n", round(sum(
  schim)/sum(st.mr) * 100, 2), "% of merged sequence reads and ", round(dim(
    schim)[2] / dim(st.mr)[2] * 100, 2), "% of merged sequence variants remain",
    " in the sequence table with chimeras removed.\nMost of your reads should ",
    "remain after chimera removal (it is not uncommon for a majority of sequen",
    "ce variants to be removed though).\nIf most of your reads were removed as",
    " chimeric, upstream processing may need to be revisited.\nIn almost all c",
    "ases this is caused by unremoved primers interfering with chimera identif",
    "cation.\nforwardtrimLeft and reversetrimLeft are the PEIV13 controlfile p",
    "rimer removing parameters for R1 and R2 raw reads respectively."),
  fill = T); return(schim)}
#e.g. Bn()

#aT function takes PEIV13 sequence table matrix default = st.chim and a PEIV13
#controlfile default uservars. returns a taxonomy sequence table.
aT <- function(s = st.chim, uv = uservars) {d <- pattern2string("databasefile")
trc <- as.logical(pattern2string("tryRC")); seed <- 100; multithread_val = T
if (.Platform$OS.type == "windows") multithread_val = F; set.seed(seed)
beg_time <- Sys.time(); cat("assignTaxonomy started\n"); taxa <- assignTaxonomy(
  s, d, multithread = multithread_val, tryRC = trc); end_time_h <- Sys.time()
elapsed_time <- paste0(round(as.numeric(difftime(
  time1 = end_time_h, time2 = beg_time, units = "hours")), 2), " Hours"); cat(
    "assignTaxonomy completed in", paste0(elapsed_time, ".\nKingdom"),
    "summary:\n"); print(table(taxa[, 1], useNA = "ifany")); return(taxa)}
#e.g. aT()

#rmBac function takes PEIV13 sequence table (default st.chim) and otu table
#default(taxa) and returns a new bacteria only sequence table
rmBac <- function(s = st.chim, t = taxa) {message("Removing ", sum(length(which(
  t[, 1] != "Bacteria")), length(which(is.na(t[, 1])))), " non Bacterial taxa")
sb <- s[, which(t[, 1] == "Bacteria"), drop = F]; return(sb)}
#e.g. st.bac <- rmBac()

#rmBac2 function takes a PEIV13 sequence table (default st.bac) and otu table
#default(taxa) and returns a new bacteria only otu table. final tm update
rmBac2 <- function(sbac = st.bac, schim = st.chim, tax = taxa) {
  taxbac <- tax[which(tax[, 1] == "Bacteria"), ]; trm <- readRDS(
    "results/seqtab/tm.rds"); name_compare3(sbac, trm, taxbac); if (any(rowSums(
      sbac) < 10000)) {rr <- rs10k(sbac); sbac <- sbac[-rr, ]; trm <- trm[-rr, ]
      if (any(colSums(sbac) == 0)) {rc <- rt0(sbac); sbac <- sbac[, -rc]
      taxbac <- taxbac[-rc, ]}; name_compare3(sbac, trm, taxbac)}; trm <- cbind(
        trm, rowSums(sbac)); colnames(trm)[8] <- "bac"; saveRDS(
          trm, "results/seqtab/tm.rds")
  taxareadloss <- (sum(schim) - sum(sbac)) / sum(schim)
  if (taxareadloss > 0.5) warning(
    "Your reads do not seem to be appropriately assigned!\n Check for lots of ",
    "16S sequences assigned as Kingdom Eukaryota NA NA NA rather than Bacteria",
    ".\nYour reads may be in the opposite orientation as the reference databas",
    "e.\nUse a different PEIV13 parameter tryRC and see if this fixes the assi",
    "gnments."); cat(paste0(round(
      taxareadloss, 4) * 100, "% of reads removed as non Bacteria\n"))
  message("Saving final TrackReads file"); write.csv(
    trm, "results/TrackReads_rep.csv"); return(taxbac)}

#make_ps function takes PEIV13 metadata, sequence table, and taxa and saves a
#phyloseq rds
make_ps <- function(m = md, s = st.bac, t = taxa.bac) {message(
  "Making phyloseq object from metadata, sequence table and taxonomy table")
  if ("repeat." %in% colnames(m)) {samdf <- m[, -which(colnames(
    m) == "library")]; samdf <- samdf[!duplicated(samdf$repeat.), ]; rownames(
      samdf) <- samdf$repeat.} else {samdf <- m; rownames(
        samdf) <- samdf$library}; library(phyloseq); ps <- phyloseq(otu_table(
          s, taxa_are_rows = F), sample_data(samdf), tax_table(t)); saveRDS(
            ps, "results/ps.rds"); pkg <- "package:phyloseq"; detach(
              pkg, character.only = T)}
#e.g. make_ps()

#rm_peiv13_int function cleans up intermediate files.
rm_peiv13_int <- function() {f <- list.files("results/seqtab", "rds")
message("Removing intermediate files"); for (i in 1:length(f)) unlink(file.path(
  "results/seqtab", f[i]), force = T); f <- list.files("results/seqtab", "rds")
if (length(f) > 0) stop("failed to remove intermediate files"); f <- list.files(
  "results/seqtab", "csv"); for (i in 1:length(f)) unlink(file.path(
    "results/seqtab", f[i]), force = T); f <- list.files(
      "results/seqtab", "csv"); if (length(f) > 0) stop(
        "failed to remove intermediate files")
unlink("results/filt", recursive = T, force = T); if (dir.exists(
  "results/filt")) stop("failed to remove intermediate files"); file.rename(
    "results/seqtab", "results/varLenHistograms"); if (dir.exists(
      "results/seqtab")) stop("failed too create histograms directory")
f <- list.files("results", "TrackReads", full.names = T); if (nrow(read.csv(
  f[1], row.names = 1)) == nrow(read.csv(f[2], row.names = 1))) {unlink(f[1])
  if (file.exists(f[1])) stop("failed to remove intermediate files")}}
#e.g. rm_peiv13_int()