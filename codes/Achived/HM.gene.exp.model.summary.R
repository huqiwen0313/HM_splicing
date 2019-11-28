# Qiwen Hu - 2019
# Perform random forest and logistic regression model by integrating gene expression and HM features

setwd("/home/qiwen/Desktop/data/git/HM_splicing")
source("./util/core.R")
source("./util/model.R")

#input.dir <- file.path("data", "HM_features", "new")
input.dir <- file.path("data", "model_results", "different_timepoints")
output.dir <- input.dir

rf.perf.all.gain.loss <- list()
rf.perf.all.high.low <- list()

logit.perf.all.gain.loss <- list()
logit.perf.all.high.low <- list()

time_point <- c("12.5day", "13.5day", "14.5day", "15.5day", "16.5day")
tissue <- c("forebrain", "heart", "hindbrain",
            "limb", "neuraltube", "liver", "midbrain")
unique_file_header <- as.vector(outer(tissue, time_point, paste, sep="."))
unique_file_header <- unique_file_header[-which(unique_file_header == "limb.16.5day"
                                                | unique_file_header == "neuraltube.16.5day")]
for(i in 1:length(unique_file_header)){
  
  sample.info <- unlist(strsplit(unique_file_header[i], "[.]"))
  timepoints <- paste(sample.info[2], sample.info[3], sep=".")
  tissue <- sample.info[1]
  
  rf.perf.gain.loss <- read.table(file.path(input.dir, paste(unique_file_header[i], "rf.gain.loss.perf.txt", sep=".")),
                        sep="\t", header = TRUE)
  rf.perf.all.gain.loss[[i]] <- round(rf.perf.gain.loss[6, ], 3)
  rf.perf.all.gain.loss[[i]]$timepoints <- timepoints
  rf.perf.all.gain.loss[[i]]$tissue <- tissue
  
  rf.perf.high.low <- read.table(file.path(input.dir, paste(unique_file_header[i], "rf.high.low.perf.txt", sep=".")),
                                  sep="\t", header = TRUE)
  rf.perf.all.high.low[[i]] <- round(rf.perf.high.low[6, ], 3)
  rf.perf.all.high.low[[i]]$timepoints <- timepoints
  rf.perf.all.high.low[[i]]$tissue <- tissue
  
}

all.rf.perf.gain.vs.loss <- dplyr::bind_rows(rf.perf.all.gain.loss)
all.rf.perf.high.vs.low <- dplyr::bind_rows(rf.perf.all.high.low)

write.table(all.rf.perf.gain.vs.loss, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.rf.perf.gain.vs.loss.txt"))
write.table(all.rf.perf.high.vs.low, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.rf.perf.high.vs.low.txt"))

all.logit.perf.gain.vs.loss <- dplyr::bind_rows(logit.perf.all.gain.loss)
all.logit.perf.high.vs.low <- dplyr::bind_rows(logit.perf.all.high.low)

write.table(all.logit.perf.gain.vs.loss, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.logit.perf.gain.vs.loss.txt"))
write.table(all.logit.perf.high.vs.low, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.logit.perf.high.vs.low.txt"))
