##Each model run may take several hours on a personal computer. Much faster with
##high performance computing (hpc) resources

#user must set variables according to THEIR data and system
source <- "/Volumes/New Volume/Tara/R/source/" #set R source script directory.
#trailing slash required.
wd <- "/Volumes/New Volume/Tara/16s_purdue/" #set working directory. working 
#directory must contain a raw data folder and an excel formatted metadata table.
raw_data_folder <- "fqgz_purge" #set raw data folder
metadata <- "Purdue_Library_Info_8_31_23.xlsx" #set metadata file name
libcol <- 6 #specify which metadata column has the illumina library identifiers
splitcol <- c(8, 9) #specify which metadata column(s) separate batches/runs
nlibs <- 317 #specify illumina library count
#user supplied variables for inspecting read quality profiles
seed <- 99 #specify the seed. Guarantees identical, reproducible results
plot_ct <- 3 #specify how many plots for each forward/reverse you want to produce
run <- 2 #specify which batch/run you want to plot. follows from splitcol
#Additional user supplied variables for modeling a parameter set
maxN_val <- 0 #specify how many Ns DADA2 will allow
maxEE_val <- c(2, 2) #specify the maximum number of “expected errors” allowed in a 
#read, which is a better filter than simply averaging quality scores.
truncQ_val <- 2 #specify the truncation quality parameter
multithread_val <- T #specify the multithread parameter
trimspec <- c(295, 250) #specify the data trimming (forward, reverse) specifications.
#this can be one of the most sensitive parameters.

#The data modeling workflow: https://benjjneb.github.io/dada2/tutorial.html
source(paste0(source,"DADA2_install_check.r")) #check DADA2 prerequisite packages
source(paste0(source,"DADA2_setup.r")) #setup the dada2 environment
source(paste0(source,"DADA2_modelTF_mergemean.r"))
