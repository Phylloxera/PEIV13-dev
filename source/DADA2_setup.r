.cran_packages <- c("ggplot2", "gridExtra", "knitr",
                    "readxl")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn",
                    "BiocStyle")
setwd(wd); sapply(c(.cran_packages, .bioc_packages), require,
                  character.only = T)
filtered_data_folder <- "filt"
miseq_path <- file.path(raw_data_folder)
start <- nchar(miseq_path) + 2; end <- nchar(miseq_path) + 10
fns <- sort(list.files(miseq_path, full.names = T))
fnFs <- fns[grepl("R1", fns)]; fnRs <- fns[grepl("R2", fns)]
metafile <- read_excel(metadata)
metafile <- metafile[order(metafile[[libcol]]), ]
if (length(fnFs) > nlibs) {
  fnFs <- fnFs[substring(fnFs, start, end) %in% metafile[[libcol]]]
  fnRs <- fnRs[substring(fnRs, start, end) %in% metafile[[libcol]]]}
filt_path <- file.path(filtered_data_folder)
if(!file_test("-d", filt_path)) dir.create(filt_path)
LibList <- split(metafile[[libcol]], apply(metafile[, splitcol], 1,
                                           paste, collapse = ""))
nbatch <- length(LibList); print(paste0("there are ", nbatch,
                                        " batches in the data"))