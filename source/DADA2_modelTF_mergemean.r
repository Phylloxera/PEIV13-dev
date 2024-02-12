set.seed(seed)
filtFs <- file.path(filt_path, basename(fnFs[substring(fnFs, start, end) %in%
                                               LibList[[run]]]))
filtRs <- file.path(filt_path, basename(fnRs[substring(fnRs, start, end) %in%
                                               LibList[[run]]]))
Fs <- fnFs[substring(fnFs, start, end) %in% LibList[[run]]]
Rs <- fnRs[substring(fnRs, start, end) %in% LibList[[run]]]
out <- filterAndTrim(Fs, filtFs, Rs, filtRs, truncLen = trimspec,
                     maxN = maxN_val, maxEE = maxEE_val, truncQ = truncQ_val,
                     rm.phix = T, compress = T, multithread = multithread_val)
errF <- learnErrors(filtFs, nbases = 1e8, multithread = multithread_val,
                    randomize = T)
errR <- learnErrors(filtRs, nbases = 1e8, multithread = multithread_val,
                    randomize = T)
dadaFs <- dada(filtFs, err = errF, multithread = multithread_val)
dadaRs <- dada(filtRs, err = errR, multithread = multithread_val)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = T)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = multithread_val, verbose = T)
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN),
               sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nochim")
print(paste0("the average number of merged reads per library in run ", run,
             " for the user provided parameter values is ", 
             round(mean(track[, 5]), 0)))
