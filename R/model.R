# functions for modeling

performance.eval <- function(pred, pred_prob, class, cat=2){
  # evaluate model performance
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

rf.model <- function(HM_count_file, AS_pattern, cv=5) {
  # random forest model to predict HMs that associated with different splicing pattern
  # Input:
  #  HM_count_file: file contains normalized HM count matrix in each exon flanking region
  #  AS_patter: a correspondent vector contain AS pattern information
  #  cv: k-fold cross validation to train the model, default=5
  # return a list contains performance evaluation and important score of each HM marker
  data <- cbind(AS_pattern, HM_count_file)
  names(data)[1] <- "class"

  k <- cv
  cv <- cvTools::cvFolds(nrow(data), K=k, R=1)

  rf_performance <- list()
  for(i in 1:k){
    test <- data[cv$subset[which(cv$which == i)],]
    train <- data[cv$subset[-1*which(cv$which == i)],]
    rf <- caret::train(train[, -1*which(colnames(train) == "class")],
                       as.factor(train$class), method="rf",
                       trControl = trainControl(method = "cv",number = 3),
                       tuneLength=10)
    rf_predict <- predict(rf, test[, -1*which(colnames(test) == "class")])
    rf_predict_prob <- predict(rf, test[, -1*which(colnames(test) == "class")], "prob")

    rf_performance[[i]] <- performance.eval(rf_predict, rf_predict_prob[, 1], test$class, 2)

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
  rf_performance_mean <- colMeans(rf_performance)
  impscore_all$ave = rowMeans(impscore_all[,2:6])
  names(impscore_all$ave) <- impscore_all$variable

  perf <- list(rf_performance_mean, data.frame(marker=impscore_all$variable, impscore = impscore_all$ave))
  names(perf) <- c("performance", "marker_importance")
  return(perf)
}

logit.model <- function(HM_count_file, AS_pattern, cv=5){
  # logistic regression to predict predict HMs that associated with different splicing pattern

}
