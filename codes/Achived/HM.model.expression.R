# this is a tempory script
setwd("/home/qiwen/Desktop/data/git/HM_splicing")
source("./util/core.R")
source("./util/model.R")

input.dir <- file.path("data", "HM_features", "new")
output.dir <- file.path("data", "model_results", "different_timepoints_rpkm")

files <- list.files(input.dir)
output.dir <- file.path("data", "model_results", "gene_expression")
rf.perf.all <- list()
for(i in 1:length(files)){
  out.file <- gsub(".HMfeatures.txt", "", files[i])
  data <- read.table(file.path(input.dir, files[i]),
                     sep = "\t", header = TRUE)
  model.perf <- model_exp(data)

  write.table(model.perf[[1]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "rf.perf.rsme.txt", sep = "")))
  write.table(model.perf[[2]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "impscore.exp.txt", sep = "")))
  ##combine performance of all models

  rf.perf.all[[i]] <- data.frame(mean.rmse = round(mean(model.perf[[1]]), 2),
                                 var.mse = round(var(model.perf[[1]]), 2))
  rownames(rf.perf.all[[i]]) <- out.file
}

rf.perf.all <- dplyr::bind_rows(rf.perf.all)
write.table(rf.perf.all, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "rf.regression.perf.txt"))
