# Qiwen Hu - 2017 - 2019

library(cvTools)
library(nnet)
library(caret)
library(pROC)
library(glmulti)
library(iRF)
library(AUC)

get_sample_info <- function(sample_file_name){
  sample_info <- unlist(strsplit(sample_file_name, split = "[.]"))
  tissue <- sample_info[1]
  histone_marker <- sample_info[5]
  time_point <- paste(sample_info[3], sample_info[4], sep = ".")
  return(c(tissue, histone_marker, time_point))
}

get_total_reads <- function(sample_file_name){
  input.dir <- file.path("data", "different_timepoints")
  sample_reads <- read.table(file.path(input.dir, "allsample.reads.txt"),
                             sep = "\t", header = FALSE)
  sample <- gsub(".1.bam.sam.hm.signal", "", sample_file_name)
  total_reads <- sample_reads[sample_reads[ ,1] == sample, 2]
  return(total_reads)
}

cal_HM_signal_flanking <- function(HM_file, marker, total_reads, exp.flag=1, exp_file = NULL){
  # caclulate chip-seq signals in the exon flanking regions
  # Input:  exp.flag - 1, incoprate gene expression in the feature
  
  if(is.null(exp_file) & exp.flag == 1){
    stop("please provide corespondent gene expression file")
  }
  
  if(exp.flag == 1){
    # normalization
    gene.exp$cpm <- (gene.exp$count/sum(gene.exp$count)) * 1e6
    gene.exp <- gene.exp[, c(1, 4)]
    names(gene.exp) <- c("GeneID", "cpm")
    HM_file <- merge(HM_file, gene.exp, by=c("GeneID"))
  }
  
  HM_signal <- list()
  for(i in 1:nrow(HM_file)){
    if(HM_file[i,]$strand == "+"){
      chip_left <- as.numeric(unlist(strsplit(as.character(HM_file[i, ]$chip_left), ",")))
      chip_right <- as.numeric(unlist(strsplit(as.character(HM_file[i, ]$chip_right), ",")))
    } else{
      chip_left <- as.numeric(unlist(strsplit(as.character(HM_file[i, ]$chip_right), ",")))
      chip_left <- rev(chip_left)
      chip_right <- as.numeric(unlist(strsplit(as.character(HM_file[i, ]$chip_left), ",")))
      chip_right <- rev(chip_right)
    }
    chip_left_intron <- (sum(chip_left[1:10]) / total_reads) * 1e6
    chip_left_exon <- (sum(chip_left[11:20]) / total_reads) * 1e6
    chip_right_intron <- (sum(chip_right[11:20]) / total_reads) * 1e6
    chip_right_exon <- (sum(chip_left[1:10]) / total_reads) * 1e6
    
    if(exp.flag == 1){
      HM_signal[[i]] <- data.frame(class = HM_file[i, ]$class, chip_left_intron = chip_left_intron,
                                   chip_left_exon = chip_left_exon, chip_right_intron = chip_right_intron,
                                   chip_right_exon = chip_right_exon, exp=HM_file[i, ]$cpm)
    } else{
      HM_signal[[i]] <- data.frame(class = HM_file[i, ]$class, chip_left_intron = chip_left_intron,
                                   chip_left_exon = chip_left_exon, chip_right_intron = chip_right_intron,
                                   chip_right_exon = chip_right_exon)
    }
  }
  HM_signal <- dplyr::bind_rows(HM_signal)
  names(HM_signal)[-1] = paste(marker, names(HM_signal)[-1], sep="_")
  return(HM_signal)
}

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

