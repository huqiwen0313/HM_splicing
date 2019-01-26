# Qiwen Hu - 2017 - 2019

library(cvTools)
library(nnet)
library(caret)
library(pROC)
library(glmulti)
library(iRF)
library(AUC)

setSeeds <- function(method = "cv", numbers = 1, repeats = 1, tunes = NULL, seed = 1237) {
  # This function is used to setup caret training random seeds for model reproducibility
  # B is the number of resamples and integer vector of M (numbers + tune length if any)
  B <- if (method == "cv") numbers
  else if(method == "repeatedcv") numbers * repeats
  else NULL
  
  if(is.null(length)) {
    seeds <- NULL
  } else {
    set.seed(seed = seed)
    seeds <- vector(mode = "list", length = B)
    seeds <- lapply(seeds, function(x) sample.int(n = 1000000, 
                                                  size = numbers + ifelse(is.null(tunes), 
                                                                          0, tunes)))
    seeds[[length(seeds) + 1]] <- sample.int(n = 1000000, size = 1)
  }
  return(seeds)
}

cal_performance <- function(pred, pred_prob, class, cat){
  # Calculate model performance 
  # Input:
  #  pred: vector contains predicted value
  #  pred_prob: vector contains probability of predicted value
  #  class: class label
  #  cat: category
  # Return:
  #  returns performance table
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
  } else{
    accuracy <- (confusion[1,1] + confusion[2,2]) / sum(confusion)
    precision <- confusion[1,1] / (confusion[1,1] + confusion[1,2])
    recall <- confusion[1,1] / (confusion[1,1] + confusion[2,1])
    f1 <- 2 * precision * recall / (precision + recall)
    auc <- pROC::auc(pROC::roc(class, pred_prob))[1]
    return(data.frame(accuracy = accuracy, Precision = precision, 
                      recall = recall, f1 = f1, auc = auc))
  }
}

