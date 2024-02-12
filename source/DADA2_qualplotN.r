set.seed(seed)
flibs <- sample(fnFs[substring(fnFs, start,
                               end) %in% LibList[[run]]], plot_ct)
rlibs <- sample(fnRs[substring(fnRs, start,
                               end) %in% LibList[[run]]], plot_ct)
for(lib in flibs) { print(plotQualityProfile(lib) + ggtitle("Fwd")) }
for(lib in rlibs) { print(plotQualityProfile(lib) + ggtitle("Rev")) }