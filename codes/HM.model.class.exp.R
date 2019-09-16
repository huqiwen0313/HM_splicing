setwd("/home/qiwen/Desktop/data/git/HM_splicing")
source("./util/core.R")
source("./util/model.R")

input.dir <- file.path("data", "HM_features", "new")
output.dir <- file.path("data", "model_results", "different_timepoints_rpkm")

files <- list.files(input.dir)
output.dir <- file.path("data", "model_results", "categories")

rf.perf.gain.vs.loss.low <- list()
rf.perf.gain.vs.loss.media <- list()
rf.perf.gain.vs.loss.high <- list()
rf.perf.high.vs.low.low <- list()
rf.perf.high.vs.low.media <- list()
rf.perf.high.vs.low.high <- list()

for(i in 1:length(files)){
  out.file <- gsub(".HMfeatures.txt", "", files[i])
  data <- read.table(file.path(input.dir, files[i]),
                     sep = "\t", header = TRUE)
  
  # classify data into different gene expression category
  exp.quantile <- quantile(data$exp)
  exp.class <- c(exp.quantile[2], exp.quantile[4])
  exp.data.low <- data[data$exp < exp.class[1], ]
  exp.data.low <- exp.data.low[, -1*which(colnames(exp.data.low) == "exp")]
  exp.data.media <- data[data$exp >= exp.class[1] & data$exp < exp.class[2], ]
  exp.data.media <- exp.data.media[, -1*which(colnames(exp.data.media) == "exp")]
  exp.data.high <- data[data$exp >= exp.class[2], ]
  exp.data.high <- exp.data.high[, -1*which(colnames(exp.data.high) == "exp")]

  elow.gain.vs.loss <- exp.data.low[exp.data.low$class == 0 | exp.data.low$class == 1, ]
  elow.high.vs.low <- exp.data.low[exp.data.low$class == 2 | exp.data.low$class == 3, ]
  emed.gain.vs.loss <- exp.data.media[exp.data.media$class == 0 | exp.data.media$class == 1, ]
  emed.high.vs.low <- exp.data.media[exp.data.media$class == 2 | exp.data.media$class == 3, ]
  ehigh.gain.vs.loss <- exp.data.high[exp.data.high$class == 0 | exp.data.high$class == 1, ]
  ehigh.high.vs.low <- exp.data.high[exp.data.high$class == 2 | exp.data.high$class == 3, ]
  
  gain.vs.loss.model.perf.low <- model(elow.gain.vs.loss)
  gain.vs.loss.model.perf.media <- model(emed.gain.vs.loss)
  gain.vs.loss.model.perf.high <- model(ehigh.gain.vs.loss)
  
  write.table(gain.vs.loss.model.perf.low[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "gain.vs.loss.low.impscore.exp.txt", sep = "")))
  write.table(gain.vs.loss.model.perf.media[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "gain.vs.loss.media.impscore.exp.txt", sep = "")))
  write.table(gain.vs.loss.model.perf.high[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "gain.vs.loss.high.impscore.exp.txt", sep = "")))
  
  high.vs.low.model.perf.low <- model(elow.high.vs.low)
  high.vs.low.model.perf.media <- model(emed.high.vs.low)
  high.vs.low.model.perf.high <- model(ehigh.high.vs.low)
  
  #write.table(model.perf.low[[1]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
  #            file.path(output.dir, paste(out.file, "rf.perf.rsme.txt", sep = "")))
  write.table(high.vs.low.model.perf.low[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "high.vs.low.low.impscore.exp.txt", sep = "")))
  write.table(high.vs.low.model.perf.media[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "high.vs.low.media.impscore.exp.txt", sep = "")))
  write.table(high.vs.low.model.perf.high[[3]], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file, "high.vs.low.high.impscore.exp.txt", sep = "")))
  
  print("done .....")
  #combine performance of all models
  rf.perf.gain.vs.loss.low[[i]] <- data.frame(mean.accuracy = round(mean(gain.vs.loss.model.perf.low[[1]]$accuracy), 2),
                                 var.mse = round(var(gain.vs.loss.model.perf.low[[1]]$accuracy), 2))
  rownames(rf.perf.gain.vs.loss.low[[i]]) <- out.file
  rf.perf.gain.vs.loss.media[[i]] <- data.frame(mean.accuracy = round(mean(gain.vs.loss.model.perf.media[[1]]$accuracy), 2),
                                              var.mse = round(var(gain.vs.loss.model.perf.media[[1]]$accuracy), 2))
  rownames(rf.perf.gain.vs.loss.media[[i]]) <- out.file
  rf.perf.gain.vs.loss.high[[i]] <- data.frame(mean.accuracy = round(mean(gain.vs.loss.model.perf.high[[1]]$accuracy), 2),
                                                var.mse = round(var(gain.vs.loss.model.perf.high[[1]]$accuracy), 2))
  rownames(rf.perf.gain.vs.loss.high[[i]]) <- out.file
  
  rf.perf.high.vs.low.low[[i]] <- data.frame(mean.accuracy = round(mean(high.vs.low.model.perf.low[[1]]$accuracy), 2),
                                             var.mse = round(var(high.vs.low.model.perf.low[[1]]$accuracy), 2))
  rownames(rf.perf.high.vs.low.low[[i]]) <- out.file
  rf.perf.high.vs.low.media[[i]] <- data.frame(mean.accuracy = round(mean(high.vs.low.model.perf.media[[1]]$accuracy), 2),
                                             var.mse = round(var(high.vs.low.model.perf.media[[1]]$accuracy), 2))
  rownames(rf.perf.high.vs.low.media[[i]]) <- out.file
  rf.perf.high.vs.low.high[[i]] <- data.frame(mean.accuracy = round(mean(high.vs.low.model.perf.high[[1]]$accuracy), 2),
                                               var.mse = round(var(high.vs.low.model.perf.high[[1]]$accuracy), 2))
  rownames(rf.perf.high.vs.low.high[[i]]) <- out.file
}

rf.perf.gain.vs.loss.low <- dplyr::bind_rows(rf.perf.gain.vs.loss.low)
rf.perf.gain.vs.loss.media <- dplyr::bind_rows(rf.perf.gain.vs.loss.media)
rf.perf.gain.vs.loss.high <- dplyr::bind_rows(rf.perf.gain.vs.loss.high)
rf.perf.high.vs.low.low <- dplyr::bind_rows(high.vs.low.model.perf.low)
rf.perf.high.vs.low.media <- dplyr::bind_rows(high.vs.low.model.perf.media)
rf.perf.high.vs.low.high <- dplyr::bind_rows(high.vs.low.model.perf.high)

write.table(rf.perf.gain.vs.loss.low, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "rf.regression.perf.gain.vs.loss.low.txt"))
write.table(rf.perf.gain.vs.loss.media, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "rf.regression.perf.gain.vs.loss.media.txt"))
write.table(rf.perf.gain.vs.loss.high, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "rf.regression.perf.gain.vs.loss.high.txt"))
write.table(rf.perf.high.vs.low.low, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "rf.regression.perf.gain.vs.loss.low.txt"))
write.table(rf.perf.high.vs.low.media, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "rf.regression.perf.gain.vs.loss.media.txt"))
write.table(rf.perf.high.vs.low.high, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "rf.regression.perf.gain.vs.loss.high.txt"))

