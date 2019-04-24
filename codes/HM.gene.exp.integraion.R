# Qiwen HU 2019
# This script is used to integrate HM signals with gene expression data and generate new feature files for modelling
setwd("/home/qiwen/Desktop/data/git/HM_splicing")
source("./util/core.R")

input.dir <- file.path("data", "different_timepoints")
gene.exp.dir <- file.path("data", "tissue.gene.count.data")
output.dir <- file.path("data", "HM_features", "new")

file.list <- list.files(input.dir)
#remove allsample.reads.txt from file.list
file.list <- file.list[-1]

time_point <- c("12.5", "13.5", "14.5", "15.5", "16.5")
tissue <- c("forebrain.mixed", "heart.mixed", "hindbrain.mixed", 
            "limb.mixed", "neuraltube.mixed", "liver.mixed", "midbrain.mixed")

unique_file_header <- as.vector(outer(tissue, time_point, paste, sep="."))

unique_file_header <- unique_file_header[-which(unique_file_header == "limb.mixed.16.5" 
                                                | unique_file_header == "neuraltube.mixed.16.5")]
# get gene length information
gene.length <- read.table(file.path(gene.exp.dir, "mm10.encode.gene.length.txt"), sep ="\t", header = TRUE)

for(file.iter in unique_file_header){
  HM_files <- file.list[grep(file.iter, file.list)]
  # read gene expression data
  gene.exp.file <- paste(gsub(".mixed", "", file.iter), "com.txt", sep = ".")
  gene.exp <- read.table(file.path(gene.exp.dir, gene.exp.file), header = TRUE)
  
  HM_flanking_features <- list()
  for(i in 1:length(HM_files)){
    total_reads <- get_total_reads(HM_files[i])
    sample_info <- get_sample_info(HM_files[i])
    HM_file <- read.table(file.path(input.dir, HM_files[i]), sep = "\t", header = TRUE)
    HM_flanking_features[[i]] <- cal_HM_signal_flanking(HM_file, sample_info[2], total_reads, exp.flag = 1, 
                                                        exp_file = gene.exp, gene.length=gene.length)
    if(i > 1){
      HM_flanking_features[[i]] <- HM_flanking_features[[i]][, -c(1, 6)]
    }
  }
  HM_flanking_features <-  dplyr::bind_cols(HM_flanking_features)
  write.table(HM_flanking_features, quote = FALSE, col.names = TRUE, 
              row.names = FALSE, sep="\t", 
              file = file.path(output.dir, paste(sample_info[1], sample_info[3], 
                                                 "HMfeatures.txt", sep = ".")))
}

