library(cvTools)
library(nnet)
library(caret)
library(pROC)
library(glmulti)
library(iRF)
library(AUC)

# function to set up seeds for caret package 
setSeeds <- function(method = "cv", numbers = 1, repeats = 1, tunes = NULL, seed = 1237) {
  # This function is used to setup caret training random seeds for model reproducibility
  #B is the number of resamples and integer vector of M (numbers + tune length if any)
  B <- if (method == "cv") numbers
  else if(method == "repeatedcv") numbers * repeats
  else NULL
  
  if(is.null(length)) {
    seeds <- NULL
  } else {
    set.seed(seed = seed)
    seeds <- vector(mode = "list", length = B)
    seeds <- lapply(seeds, function(x) sample.int(n = 1000000, 
                                                  size = numbers + ifelse(is.null(tunes), 0, tunes)))
    seeds[[length(seeds) + 1]] <- sample.int(n = 1000000, size = 1)
  }
  return(seeds)
}

cal_performance <- function(pred, pred_prob, class, cat){
  
  confusion <- as.matrix(table(class, pred, deparse.level = 0))
  
  if(cat > 2){
    n <- sum(confusion) # number of instances
    nc <- nrow(confusion) # number of classes
    diag <- diag(confusion) # number of correctly classified instances per class 
    accuracy <- sum(diag)/sum(confusion)
    rowsums <- apply(confusion, 1, sum) # number of instances per class
    colsums <- apply(confusion, 2, sum) # number of predictions per class
    p <- rowsums / n # distribution of instances over the actual classes
    q <- colsums / n # distribution of instances over the predicted classes
    precision <- diag / colsums 
    recall <- diag / rowsums 
    f1 <- 2 * precision * recall / (precision + recall)
    macroPrecision <- mean(precision, na.rm = T)
    macroRecall <- mean(recall, na.rm = T)
    macroF1 <- mean(f1, na.rm = T)
    return(data.frame(accuracy = accuracy, macroPrecision = macroPrecision, 
                      macrorecall = macroRecall, macrof1 = macroF1))
  } else {
    accuracy <- (confusion[1,1] + confusion[2,2]) / sum(confusion)
    precision <- confusion[1,1] / (confusion[1,1] + confusion[1,2])
    recall <- confusion[1,1] / (confusion[1,1] + confusion[2,1])
    f1 <- 2 * precision * recall / (precision + recall)
    auc <- pROC::auc(pROC::roc(class, pred_prob))[1]
    return(data.frame(accuracy = accuracy, Precision = precision, 
                      recall = recall, f1 = f1, auc = auc))
  }
}

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

rf_perf_vis <- function(files){
  input.dir <- file.path("data", "model_results", "different_timepoints")
  perf.all <- list()
  for(j in 1:length(files)){
    perf <- read.table(file.path(input.dir, files[j]), sep = "\t", header = TRUE)
    perf <- perf[, c(1, 5)]
    time_point <- gsub(paste(tissue[i], ".", sep = ""), "", files[j])
    time_point <- gsub(".rf.gain.loss.perf.txt", "", time_point)
    time_point <- gsub(".rf.high.low.perf.txt", "", time_point)
    perf$tissue <- tissue[i]
    perf$timepoint <- time_point
    perf.all[[j]] <- perf
  }
  perf.all <- dplyr::bind_rows(perf.all)
  perf.all <- reshape2::melt(perf.all)
  perf.all <- data_summary(perf.all, varname="value", 
                                 groupnames=c("variable", "timepoint"))
  plot <- ggplot2::ggplot(perf.all, aes(x=timepoint, y=value, group=variable, color = variable)) +
    geom_line() +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                  position=position_dodge(0.05)) + 
    theme_classic() +
    scale_color_manual(values=c('#999999','#E69F00')) + 
    theme(plot.title = element_text(hjust = 0.5), legend.position="top") + 
    ggplot2::ggtitle(tissues[i])
  return(plot)
}

model <- function(data) {
  
  #set seed for external cross-validation
  set.seed(11223)
  # set seed for caret training
  logit_caretSeeds <- setSeeds(method = "cv", numbers = 3, repeats = 1)
  rf_caretSeeds <- setSeeds(method = "cv", numbers = 3, repeats = 1, tunes = 10)
  
  data$class <- factor(data$class)
  
  k <- 5
  cv <- cvTools::cvFolds(nrow(data), K=k, R=1)
  
  rf_performance <- list()
  logit_performance <- list()

  for(i in 1:k){
    test <- data[cv$subset[which(cv$which == i)],]
    train <- data[cv$subset[-1*which(cv$which == i)],]
    
    #m <- multinom(class ~ ., data = train)
    logit.lass <- caret::train(train[, -1*which(colnames(train) == "class")],
                        as.factor(train$class), 
                        method = "glmnet", 
                        trControl = trainControl(method = "cv",number = 3,
                                                 seeds = logit_car
                                                 setSeeds), 
                        family="binomial",
                        tuneGrid = expand.grid(.alpha=1, .lambda=seq(0, 10, by = 0.01)))
    
    rf <- caret::train(train[, -1*which(colnames(train) == "class")],
                as.factor(train$class), method="rf", 
                trControl = trainControl(method = "cv",number = 3, seeds = rf_caretSeeds), 
                tuneLength=10)
    #mod_fit <- train(class ~ .,  data=train, method="glm", family="binomial")
    rf_predict <- predict(rf, test[, -1*which(colnames(test) == "class")])
    logit_predict <- predict(logit.lass, test[, -1*which(colnames(test) == "class")])
    logit_predict_prob <- predict(logit.lass, newdata = test[, -1*which(colnames(test) == "class")],
                                  "prob")
    rf_predict_prob <- predict(rf, test[, -1*which(colnames(test) == "class")], "prob")
    
    rf_performance[[i]] <- cal_performance(rf_predict, rf_predict_prob[, 1], test$class, 2)
    logit_performance[[i]] <- cal_performance(logit_predict, logit_predict_prob[, 1], test$class, 2)
    
    impscore = varImp(rf$finalModel)
    impscore$variable = rownames(impscore)
    colnames(impscore)[1] = i
    if(i == 1){
      impscore_all = impscore
    } else{
      impscore_all = merge(impscore_all, impscore, by=c("variable"))
    }
  }
  
  rf_performance <- dplyr::bind_rows(rf_performance)
  rf_performance <- rbind(rf_performance, colMeans(rf_performance))
  logit_performance <- dplyr::bind_rows(logit_performance)
  logit_performance <- rbind(logit_performance, colMeans(logit_performance))
  impscore_all$ave = rowMeans(impscore_all[,2:6])
  
  return(list(rf_performance, logit_performance, impscore_all))
}

iRF_model <- function(data, tissue, timepoint, n.iter = 20, n.interaction = 10, bootstrap = 30){
  set.seed(11223)
  data$class <- factor(data$class)
  k <- 5
  cv <- cvTools::cvFolds(nrow(data), K=5, R=1)
  interaction.rank <- list()
  for(i in 1:k){
    test <- data[cv$subset[which(cv$which == i)],]
    train <- data[cv$subset[-1*which(cv$which == i)],]
    train.X <- train[, -which(colnames(train) == "class")]
    train.Y <- train$class
    test.X <- test[, -(colnames(test) == "class")]
    test.Y <- test$class
    ff <- iRF::iRF(x = as.matrix(train.X), 
              y=train.Y, 
              xtest = as.matrix(test.X), 
              ytest = test.Y, 
              n.iter = n.iter, 
              interactions.return = n.interaction,
              n.bootstrap = bootstrap
    )
    iter.interaction <- ff$interaction[[n.interaction]]
    interaction.rank[[i]] <- data.frame(interaction = names(iter.interaction), 
                                       score = as.vector(iter.interaction))
    names(interaction.rank[[i]])[2] <- paste("score", i, sep = ".")
    #roc.info <- AUC::roc(ff$rf.list[[20]]$test$votes[,2], test.Y)
    #rauc <- round(100*auc(roc.info), 2)
  }
  interaction.rank <- Reduce(function(x, y) merge(x, y, by = "interaction"), interaction.rank)
  interaction.rank$ave <- rowMeans(interaction.rank[, 2:(k+1)])
  interaction.rank <- interaction.rank[order(interaction.rank$ave, decreasing = TRUE), ]
  
  # Generate interaction plot
  interaction.list <- interaction.rank[, c(1, ncol(interaction.rank))]
  p <- HM_interaction_plot(interaction.list, tissue = tissue, timepoint = timepoint)
  
  return(list(interaction.rank, p))
}

# Main ----
# Model performance for different tissues in same timepoints
input.dir <- file.path("data", "processed", "HM_features")
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

###combine all performance results for logistic regression (removed in the future)##
files <- list.files(output.dir)
files <- files[grep("logit", files)]
gain.vs.loss <- list()
high.vs.low <- list()
for(i in 1:length(files)){
  file <- read.table(file.path(output.dir, files[i]), sep = "\t", header = TRUE)
  if(grepl("gain", files[i])){
    gain.vs.loss[[i]] <- file[6, ]
  } else{
    high.vs.low[[i]] <- file[6, ]
  }
}
all.logit.perf.gain.vs.loss <- dplyr::bind_rows(gain.vs.loss)
all.logit.perf.high.vs.low <- dplyr::bind_rows(high.vs.low)
write.table(all.logit.perf.gain.vs.loss, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.logit.perf.gain.vs.loss.txt"))
write.table(all.logit.perf.high.vs.low, quote = FALSE, col.names = TRUE, row.names = TRUE, sep="\t",
            file.path(output.dir, "alldata.logit.perf.high.vs.low.txt"))
  
#visulize model performance for different tissues ###
input.dir <- file.path("data", "model_results", "different_timepoints")
output.dir <- file.path("data", "figures", "model")
rf.high.vs.low <- read.table(file.path(input.dir, "alldata.rf.perf.high.vs.low.txt"), 
                             sep = "\t", header = TRUE)
logit.high.vs.low <- read.table(file.path(input.dir, "alldata.logit.perf.high.vs.low.txt"), 
                             sep = "\t", header = TRUE)
rf.gain.vs.loss <- read.table(file.path(input.dir, "alldata.rf.perf.gain.vs.loss.txt"), 
                             sep = "\t", header = TRUE)
logit.gain.vs.loss <- read.table(file.path(input.dir, "alldata.logit.perf.gain.vs.loss.txt"), 
                                sep = "\t", header = TRUE)


high.vs.low.pl <- plot_model_performance(rf.high.vs.low, logit.high.vs.low) 
gain.vs.loss.pl <- plot_model_performance(rf.gain.vs.loss, logit.gain.vs.loss)

pdf(file.path(output.dir, "high.vs.low.model.perf.plot.pdf"), width = 5, height = 15)
cowplot::plot_grid(plotlist = high.vs.low.pl, ncol = 1)
dev.off()

pdf(file.path(output.dir, "gain.vs.loss.model.perf.plot.pdf"), width = 5, height = 15)
cowplot::plot_grid(plotlist = gain.vs.loss.pl, ncol = 1)
dev.off()

# Visulize random forest performance
input.dir <- file.path("data", "model_results", "different_timepoints")
output.dir <- file.path("data", "figures", "model")
rf.high.vs.low <- read.table(file.path(input.dir, "alldata.rf.perf.high.vs.low.txt"), 
                             sep = "\t", header = TRUE)
rf.gain.vs.loss <- read.table(file.path(input.dir, "alldata.rf.perf.gain.vs.loss.txt"), 
                              sep = "\t", header = TRUE)
rf.high.vs.low.pl <- plot_rf_performance(rf.high.vs.low) 
rf.gain.vs.loss.pl <- plot_rf_performance(rf.gain.vs.loss)


pdf(file.path(output.dir, "rf.high.vs.low.model.perf.plot.pdf"), width = 10, height = 20)
cowplot::plot_grid(plotlist = rf.high.vs.low.pl, ncol = 2)
dev.off()

pdf(file.path(output.dir, "rf.gain.vs.loss.model.perf.plot.pdf"), width = 10, height = 20)
cowplot::plot_grid(plotlist = rf.gain.vs.loss.pl, ncol = 2)
dev.off()

# visulize random forest performance, add variation
input.dir <- file.path("data", "model_results", "different_timepoints")
output.dir <- file.path("data", "figures", "model")
tissue <- c("forebrain", "midbrain", "hindbrain", "heart", "limb", "liver", "neuraltube")
files <- list.files(input.dir)
gain.loss.plot <- list()
high.low.plot <- list()
for(i in 1:length(tissue)){
  rf.file.gain.loss <- files[grep(paste(tissue[i], ".*[.]rf.gain.loss.*", sep = ""), files)]
  rf.file.high.low <- files[grep(paste(tissue[i], ".*[.]rf.high.low.*", sep = ""), files)]
  gain.loss.plot[[i]] <- rf_perf_vis(rf.file.gain.loss)
  high.low.plot[[i]] <- rf_perf_vis(rf.file.high.low)
}
  
pdf(file.path(output.dir, "rf.high.vs.low.model.perf.plot.pdf"), width = 10, height = 20)
cowplot::plot_grid(plotlist = high.low.plot, ncol = 2)
dev.off()

pdf(file.path(output.dir, "rf.gain.vs.loss.model.perf.plot.pdf"), width = 10, height = 20)
cowplot::plot_grid(plotlist = gain.loss.plot, ncol = 2)
dev.off()


#HM interations identification and visulization

input.dir <- file.path("data", "processed", "HM_features")
output.dir <- file.path("data", "model_results", "different_timepoints", "feature_interaction")

files <- list.files(input.dir)

gain.vs.loss.interaction.plot <- list()
high.vs.low.interaction.plot <- list()

for(i in 1:length(files)){
  out.file <- gsub("HMfeatures.txt", "", files[i])
  sample_info <- unlist(strsplit(out.file, split = "[.]"))
  tissue <- sample_info[1]
  timepoint <- paste(sample_info[2], sample_info[3], sep = ".")
  
  data <- read.table(file.path(input.dir, files[i]),
                     sep = "\t", header = TRUE)
  gain.vs.loss <- data[data$class == 0 | data$class == 1, ]
  high.vs.low <- data[data$class == 2 | data$class == 3, ]
  
  gain.vs.loss.interaction <- iRF_model(gain.vs.loss, tissue = tissue, timepoint = timepoint)
  gain.vs.loss.interaction.plot[[i]] <-  gain.vs.loss.interaction[2]
  write.table(gain.vs.loss.interaction[1], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file,"gain.loss.interaction.txt", sep = "")))
 
  
  high.vs.low.interaction <- iRF_model(high.vs.low, tissue = tissue, timepoint = timepoint)
  high.vs.low.interaction.plot[[i]] <-  high.vs.low.interaction[2]
  write.table(high.vs.low.interaction[1], quote = FALSE, col.names = TRUE, row.names = FALSE, sep="\t",
              file.path(output.dir, paste(out.file,"high.low.interaction.txt", sep = "")))
}

pdf(file.path(output.dir, "gain.vs.loss.interaction.pdf"), width = 20, height = 15)
cowplot::plot_grid(plotlist = gain.vs.loss.interaction.plot, ncol = 5)
dev.off()

pdf(file.path(output.dir, "high.vs.low.interaction.pdf"), width = 20, height = 15)
cowplot::plot_grid(plotlist = high.vs.low.interaction.plot, ncol = 5)
dev.off()
