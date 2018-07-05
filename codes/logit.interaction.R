

logit.model <- function(x, y) {
  
  #set seed for external cross-validation
  set.seed(11223)
  # set seed for caret training
  logit_caretSeeds <- setSeeds(method = "cv", numbers = 3, repeats = 1)
  
  k <- 5
  cv <- cvTools::cvFolds(nrow(x), K=k, R=1)
  
  logit_performance <- list()
  
  for(i in 1:k){
    train.x <- x[cv$subset[-1*which(cv$which == i)], ]
    train.y <- y[cv$subset[-1*which(cv$which == i)]]
    test.x <- x[cv$subset[which(cv$which == i)], ]
    test.y <- y[cv$subset[which(cv$which == i)]]
    
    logit.lass <- caret::train(train.x,
                               train.y, 
                               method = "glmnet", 
                               trControl = trainControl(method = "cv",number = 3,
                                                        seeds = logit_caretSeeds), 
                               family="binomial",
                               tuneGrid = expand.grid(.alpha=1, .lambda=seq(0, 10, by = 0.01)))
    
    logit_predict <- predict(logit.lass, test.x)
    logit_predict_prob <- predict(logit.lass, newdata = test.x,
                                  "prob")
    logit_performance[[i]] <- cal_performance(logit_predict, logit_predict_prob[, 1], test.y, 2)
    
  }
  
  logit_performance <- dplyr::bind_rows(logit_performance)
  logit_performance <- rbind(logit_performance, colMeans(logit_performance))
  return(logit_performance)
}


setwd("E:/work/penn/greene/HM_splicing_motifs/data/processed/HM_features")
file <- read.table("heart.16.5day.HMfeatures.txt", sep = "\t", header = T)
file <- file[file$class == 0 | file$class == 1, ]

x <- as.matrix(file[, -which(colnames(file) == "class")])
class <- as.factor(file$class)
f <- as.formula(paste("class ~ ", paste(c(colnames(x), "H3K27ac_chip_right_intron:H3K36me3_chip_right_intron")
                      , collapse= "+")))
f <- as.formula(class ~ .*.) 
x.interaction <- model.matrix(f, file)[, -1]
logit.performance <- logit.model(x.interaction, class)

logit.performance.orignal <- logit.model(x, class)
