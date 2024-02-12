filtFs <- file.path(filt_path, basename(fnFs[substring(fnFs, start, end) %in% LibList[[run]]]))
filtRs <- file.path(filt_path, basename(fnRs[substring(fnRs, start, end) %in% LibList[[run]]]))
Fs <- fnFs[substring(fnFs, start, end) %in% LibList[[run]]]; Rs <- fnRs[substring(fnRs, start, end) %in% LibList[[run]]]
out <- filterAndTrim(Fs, filtFs, Rs, filtRs, truncLen = trimspec, maxN = maxN_val, maxEE = maxEE_val, truncQ = truncQ_val, rm.phix = T, compress = T, multithread = T)
if (exists(quote(outdf))) {outdf <- rbind(outdf, data.frame(out))} else {outdf <- data.frame(out)}
errF[[run]] <- learnErrors(filtFs, nbases = 1e8, multithread = T, randomize = T); errR[[run]] <- learnErrors(filtRs, nbases = 1e8, multithread = T, randomize = T)
dadaFs <- dada(filtFs, err = errF[[run]], multithread = T); dadaRs <- dada(filtRs, err = errR[[run]], multithread = T)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = T); seqtab <- makeSequenceTable(mergers)
stchim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = T, verbose = T); print(dim(stchim)); print(sum(stchim)/sum(seqtab))
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(stchim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "seqtab", "nochim"); sample.names <- LibList[[run]]
rownames(track) <- sample.names; tracks[[run]] <- track; saveRDS(track, paste0("track", run, ".rds")); saveRDS(stchim, paste0("stchim", run, ".rds"))
