install.packages("caret")
install.packages("cvTools")
install.packages("ggplot2")
install.packages("e1071")
install.packages("reshape2")
install.packages("pROC")
install.packages("randomForest")
install.packages("iRF")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")
