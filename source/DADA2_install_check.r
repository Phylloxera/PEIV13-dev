{
.cran_packages <- c("ggplot2", "gridExtra", "knitr",
                    "readxl")
.bioc_packages <- c("dada2", "phyloseq", "DECIPHER", "phangorn",
                    "BiocStyle")
.inst <- .cran_packages %in% installed.packages()
if(any(!.inst)) {install.packages(.cran_packages[!.inst])}
.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {BiocManager::install(.bioc_packages[!.inst],
                                      ask = F)}
}
