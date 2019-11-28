
setwd("E:/work/penn/greene/HM_splicing_motifs/data/processed/HM_features")
data = read.table("forebrain.12.5day.HMfeatures.txt", sep = "\t", header = TRUE)
data = data[data$class == 0 | data$class == 1, ]

n <- 200
p <- 250
X <- matrix(rnorm(n * p), nrow=n)

Y <- (X[,1] > 0.35 & X[,2] > 0.35) | (X[,5] > 0.35 & X[,7] > 0.35)
Y <- as.factor(as.numeric(Y > 0))

train.id <- 1:(n / 2)
test.id <- setdiff(1:n, train.id)

ff <- iRF(x=X[train.id,], 
          y=Y[train.id], 
          xtest=X[test.id,], 
          ytest=Y[test.id], 
          n.iter=5, 
          interactions.return=5,
          n.bootstrap=10
)
ff$interaction

toplot <- rev(ff$interaction[[3]])
dotchart(as.matrix(toplot[189:199]), xlab='Stability Score', 
         main='Prevalent Features/Interactions \n on Decision paths')


library(AUC)
plot(0:1, 0:1, type='l', lty = 2, xlab = 'FPR', ylab = 'TPR', main='ROC Curve')
for (iter in 1:4){
  # performance on test set
  cat(paste('iter = ', iter, ':: '))
  roc.info <- roc(ff$rf.list[[iter]]$test$votes[,2], Y[test.id])
  lines(roc.info$fpr, roc.info$tpr, type='l', col=iter, lwd=2)
  cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep=''))
} 
legend('bottomright', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')