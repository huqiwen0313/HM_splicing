# Qiwen Hu - 2019
# Perform random forest and logistic regression model by integrating gene expression and HM features

setwd("/home/qiwen/Desktop/data/git/HM_splicing")
source("./util/core.R")
source("./util/model.R")

input.dir <- file.path("data", "HM_features")
output.dir <- file.path("data", "model_results", "different_timepoints")

files <- list.files(input.dir)
rf.perf.all <- list()
logit.perf.all <- list()

for(i in 1:length(files)){
  out.file <- gsub("HMfeatures.txt", "", files[i])
  data <- read.table(file.path(input.dir, files[i]),
                     sep = "\t", header = TRUE)
  gain.vs.loss <- data[data$class == 0 | data$class == 1, ]
  high.vs.low <- data[data$class == 2 | data$class == 3, ]
  gain.vs.loss.perf <- model(gain.vs.loss)
  write.table(gain.vs.loss.perf[[1]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "rf.gain.loss.perf.txt", sep = "")))
  write.table(gain.vs.loss.perf[[2]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "logit.gain.loss.perf.txt", sep = "")))
  write.table(gain.vs.loss.perf[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "impscore.gain.loss.txt", sep = "")))
  
  high.vs.low.perf <- model(high.vs.low)
  write.table(high.vs.low.perf[[1]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "rf.high.low.perf.txt", sep = "")))
  write.table(high.vs.low.perf[[2]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "logit.high.low.perf.txt", sep = "")))
  write.table(high.vs.low.perf[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "impscore.high.low.txt", sep = "")))
  ##combine performance of all models
  rf.perf.all.gain.vs.loss[[i]] <- gain.vs.loss.perf[[1]][6, ]
  #rownames(rf.perf.all[[i]]) <- out.file
  rf.perf.all.high.vs.low[[i]] <- high.vs.low.perf[[1]][6, ]
  #rownames(logit.perf.all[[i]]) <- out.file
}

names <- gsub(".HMfeatures.txt", "", files)
all.rf.perf.gain.vs.loss <- dplyr::bind_rows(rf.perf.all.gain.vs.loss)
all.rf.perf.high.vs.low <- dplyr::bind_rows(rf.perf.all.high.vs.low)
rownames(all.rf.perf.gain.vs.loss) <- names
rownames(all.rf.perf.high.vs.low) <- names

write.table(all.rf.perf.gain.vs.loss, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.rf.perf.gain.vs.loss.txt"))
write.table(all.rf.perf.high.vs.low, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.rf.perf.high.vs.low.txt"))
