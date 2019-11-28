# Qiwen Hu - 2019
# This script is used to calculate the correlation between gene expression and psi
#  for each splicing event

input.dir <- "/data/tissue.gene.count.data/"
setwd(input.dir)

# list file names
files <- dir(input.dir)

cor.table <- data.frame()
for(i in 1:length(files)){
  f <- read.table(files[i], sep = "\t", header = TRUE)
  
  # get sample name
  name <- unlist(strsplit(files[i], split = "[.]"))
  tissue <- name[1]
  timepoint <- paste(name[2], name[3], sep = ".")
  
  # Normalized gene expression CPM
  f$count <- (f$count/sum(f$count)) * 1e6
  correlation <- round(cor(f$count, f$psi), 2)
  p.value <- cor.test(f$count, f$psi)$p.value
  
  if(i == 1){
    cor.table <- data.frame(tissue = tissue, timepoint = timepoint, cor = correlation, p.value = p.value)
  } else{
    table <- data.frame(tissue = tissue, timepoint = timepoint, cor = correlation, p.value = p.value)
    cor.table <- rbind(cor.table, table)
  }
}

# plot data
pvalues <- cor.table$p.value
names(pvalues) <- paste(cor.table$tissue, cor.table$timepoint, sep = " ")

pdf("correlation.plot.pdf")
barplot(pvalues, las = 2, ylab = "p-value", cex.names = 0.7)
abline(h = 0.01, col = "red")
dev.off()

write.table(cor.table, quote = FALSE, col.names = TRUE,
            row.names = FALSE, sep = "\t", file = "exp.psi.cor.table")