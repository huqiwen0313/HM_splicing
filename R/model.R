# functions for modeling

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

modelPerformance <- function(pred, pred_prob, class, cat){

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

data_summary <- function(data, varname, groupnames){
  # Function to calculate the mean and the standard deviation
  # for each group
  # data : a data frame
  # varname : the name of a column containing the variable
  #to be summariezed
  # groupnames : vector of column names to be used as
  # grouping variables

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

    logit.lass <- caret::train(train[, -1*which(colnames(train) == "class")],
                               as.factor(train$class),
                               method = "glmnet",
                               trControl = caret::trainControl(method = "cv",number = 3,
                                                        seeds = logit_caretSeeds,
                                                        setSeeds),
                               family="binomial",
                               tuneGrid = expand.grid(.alpha=1, .lambda=seq(0, 10, by = 0.01)))

    rf <- caret::train(train[, -1*which(colnames(train) == "class")],
                       as.factor(train$class), method="rf",
                       trControl = caret::trainControl(method = "cv",number = 3, seeds = rf_caretSeeds),
                       tuneLength=10)
    rf_predict <- predict(rf, test[, -1*which(colnames(test) == "class")])
    logit_predict <- predict(logit.lass, test[, -1*which(colnames(test) == "class")])
    logit_predict_prob <- predict(logit.lass, newdata = test[, -1*which(colnames(test) == "class")],
                                  "prob")
    rf_predict_prob <- predict(rf, test[, -1*which(colnames(test) == "class")], "prob")

    rf_performance[[i]] <- modelPerformance(rf_predict, rf_predict_prob[, 1], test$class, 2)
    logit_performance[[i]] <- modelPerformance(logit_predict, logit_predict_prob[, 1], test$class, 2)

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

  return(list(rf.performance=rf_performance, logit.performance=logit_performance, impscore=impscore_all))
}

model_exp <- function(data) {
  # random forest regression model to identify HMs asscociated with expression
  # input: data - same as the input in model function

  #set seed for external cross-validation
  set.seed(11223)

  # now we are look at expression, remove exon type classification label
  data <- data[, -1]

  k <- 5
  cv <- cvTools::cvFolds(nrow(data), K=k, R=1)

  rf_performance <- c()

  for(i in 1:k){
    test <- data[cv$subset[which(cv$which == i)],]
    train <- data[cv$subset[-1*which(cv$which == i)],]

    rf.fit <- caret::train(exp ~ ., data = train,  method = 'ranger',
                           tuneLength = 10,
                           trControl = trainControl(method = "cv",number = 3),
                           num.trees = 700,
                           importance = "permutation")

    testPred <- predict(rf.fit , test[, -1*which(colnames(test) == "exp")])
    rmse.score <- round(RMSE(testPred,test$exp), 2)

    rf_performance <- c(rf_performance, rmse.score)
    impscore = varImp(rf.fit)$importance
    impscore$variable = rownames(impscore)
    colnames(impscore)[1] = i
    if(i == 1){
      impscore_all = impscore
    } else{
      impscore_all = merge(impscore_all, impscore, by=c("variable"))
    }
  }

  impscore_all$ave = rowMeans(impscore_all[,2:6])
  return(list(rf_performance, impscore_all))
}

iRFModel <- function(data, tissue, timepoint, n.iter = 20, n.interaction = 10, bootstrap = 30){
  # run interative random forest model
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
  }
  interaction.rank <- Reduce(function(x, y) merge(x, y, by = "interaction"), interaction.rank)
  interaction.rank$ave <- rowMeans(interaction.rank[, 2:(k+1)])
  interaction.rank <- interaction.rank[order(interaction.rank$ave, decreasing = TRUE), ]

  # Generate interaction plot
  interaction.list <- interaction.rank[, c(1, ncol(interaction.rank))]
  p <- HM_interaction_plot(interaction.list, tissue = tissue, timepoint = timepoint)

  return(list(interaction.rank, p))
}
