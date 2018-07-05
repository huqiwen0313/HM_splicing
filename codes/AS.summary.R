library(ggplot2)

##summarize AS statistics

input.dir <- file.path("data", "processed", "different_timepoints")
output.dir <- file.path("data", "processed", "AS_statistics_summary")

files <- list.files(input.dir)

summary.stat <- data.frame()
lineage.specific.per <- c()
SE.alldata.gain.vs.loss <- list()
SE.alldata.high.vs.low <- list()
timepoints <- c("12.5day", "13.5day", "14.5day", "15.5day", "16.5day")
for(i in 1:length(timepoints)){
  #SE.high.vs.low <- list()
  SE.gain.vs.loss <- list()
  SE.high.vs.low <- list()
  file <- files[grep(paste(timepoints[i], "H3K27ac", sep = "."), files)]
  for(j in 1:length(file)){
    data <- read.table(file.path(input.dir, file[j]), sep = "\t", header = T)
    tissue <- unlist(strsplit(file[j], "[.]"))[1]
    gain.vs.loss <- data[data$class == 0 | data$class == 1, ]
    high.vs.low <- data[data$class == 2 | data$class == 3, ]
    gain.num <- nrow(data[data$class == 0, ])
    loss.num <- nrow(data[data$class == 1, ])
    high.num <- nrow(data[data$class == 2, ])
    low.num <- nrow(data[data$class == 3, ])
    
    gain.vs.loss.gene.num <- length(unique(gain.vs.loss$geneSymbol))
    high.vs.low.gene.num <- length(unique(high.vs.low$geneSymbol))
    
    if( i == 1 & j == 1){
      summary.stat <- data.frame(timepoint = timepoints[i], tissue = tissue, 
               total_AS_number = nrow(data), gain_num = gain.num, loss_num = loss.num,
               gain.loss.gene.num = gain.vs.loss.gene.num,
               high.low.gene.num = high.vs.low.gene.num,
               high_num = high.num, low_num = low.num)
    } else {
      tmp <- data.frame(timepoint = timepoints[i], tissue = tissue, 
                        total_AS_number = nrow(data), gain_num = gain.num, loss_num = loss.num,
                        gain.loss.gene.num = gain.vs.loss.gene.num,
                        high.low.gene.num = high.vs.low.gene.num,
                        high_num = high.num, low_num = low.num)
      summary.stat <- rbind(summary.stat, tmp)
    }
    SE.gain.vs.loss[[j]] <- unique(gain.vs.loss[, 2:3])
    SE.gain.vs.loss[[j]]$tissue <- tissue
    SE.high.vs.low[[j]] <- unique(high.vs.low[, 2:3])
    SE.high.vs.low[[j]]$tissue <- tissue
    #write.table(SE.gain.vs.loss[[j]], quote = FALSE, col.names = TRUE, 
    #            row.names = FALSE, sep="\t", 
    #            file = file.path(output.dir, paste(tissue, timepoints[i], "sig.SE.txt",
    #                                               sep = ".")))
  }
  SE.all.gain.vs.loss <- dplyr::bind_rows(SE.gain.vs.loss)
  SE.alldata.gain.vs.loss[[i]] <- SE.all.gain.vs.loss
  SE.all.high.vs.low <- dplyr::bind_rows(SE.high.vs.low)
  SE.alldata.high.vs.low[[i]] <- SE.all.high.vs.low
  for( k in 1:length(SE.gain.vs.loss)){
    all.AS <- unique(SE.all.gain.vs.loss[-which(SE.all.gain.vs.loss$tissue == unique(SE.gain.vs.loss[[k]]$tissue)), 
                            ]$geneSymbol)
    comm <- intersect(unique(SE.gain.vs.loss[[k]]$geneSymbol), all.AS)
    lineage.specific <- 1 - (length(comm)/length(unique(SE.gain.vs.loss[[k]]$geneSymbol)))
    lineage.specific.per <- c(lineage.specific.per, lineage.specific)
  }
}
summary.stat$lineage.specific.percentage = lineage.specific.per
write.table(summary.stat, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep="\t", 
            file = file.path(output.dir, "AS.summary.stat.txt"))

###
SE.alldata <- dplyr::bind_rows(SE.alldata.gain.vs.loss)
SE.alldata <- unique(SE.alldata)
write.table(SE.alldata, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep="\t", 
            file = file.path(output.dir, "all.gain.vs.loss.SE.gene.txt"))

SE.alldata <-  dplyr::bind_rows(SE.alldata.high.vs.low)
SE.alldata <- unique(SE.alldata)
write.table(SE.alldata, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep="\t", 
            file = file.path(output.dir, "all.high.vs.low.SE.gene.txt"))

tissues <- unique(SE.alldata$tissue)
lineage.specific.per <- c()
for(i in 1:length(tissues)){
  all.AS <- unique(SE.alldata[-which(SE.alldata$tissue == tissues[i]), 
                          ]$geneSymbol)
  tissue.AS <- SE.alldata[-which(SE.alldata$tissue == tissues[i]), ]$geneSymbol
  comm <- intersect(tissue.AS, all.AS)
  lineage.specific <- 1 - (length(comm)/length(tissue.AS))
  lineage.specific.per <- c(lineage.specific.per, lineage.specific)
}
lineage.percentage <- data.frame(tissue = tissues, percentage = round(lineage.specific.per*100, 2))
#plot
ggplot2::ggplot(lineage.percentage, aes(x = tissue, y = percentage, fill = "blue")) + 
  ggplot2::geom_bar(stat="identity", color="black", position=position_dodge()) +
  ggplot2::scale_fill_brewer(palette="Paired") + ggplot2::theme_minimal() + 
  geom_text(aes(label = percentage), vjust= -0.2, color="black",
                                       position = position_dodge(0.9), size = 3)

###top function##
setwd("E:/work/penn/greene/HM_splicing_motifs/data/processed/AS_statistics_summary/function")
data <- read.table("forebrain.top20.high.vs.low.txt", sep = "\t", header = TRUE)
data$pvalue <- -1*log(data$pvalue)/log(10)

ggplot2::ggplot(data, aes(x = reorder(GO_function, pvalue), y = pvalue, fill = "blue")) + 
  ggplot2::geom_bar(stat="identity", color="black", position=position_dodge()) +
  ggplot2::scale_fill_brewer(palette="Paired") + ggplot2::theme_minimal() + 
  coord_flip() + ggplot2::theme(legend.position="none") + 
  ggplot2::xlab("GO function") + ggplot2::ylab("-log10 p-value") 

# pairwsie lineage overlap heatmap
input.dir <- file.path("data", "processed", "AS_statistics_summary", "gain.vs.loss.summary")
tissues <- c("forebrain", "hindbrain", "midbrain", "neuraltube", "heart", "limb", "liver")

overlap.matrix <- matrix(0, nrow = length(tissues), ncol = length(tissues))
for(i in 1:length(tissues)){
  tissue_1 <- read.table(file.path(input.dir, paste(tissues[i], "all.as.txt", sep = ".")), 
                                   sep = "\t", header = TRUE)
  for(j in 1:length(tissues)){
    tissue_2 <- read.table(file.path(input.dir, paste(tissues[j], "all.as.txt", sep = ".")), 
                           sep = "\t", header = TRUE)
    all.events <- unique(c(as.vector(tissue_1$V2), as.vector(tissue_2$V2)))
    intersect.events <- intersect(unique(tissue_1$V2), unique(tissue_2$V2))
    overlap.rate <- round(length(intersect.events) / length(all.events), 2)
    overlap.matrix[i, j] <- overlap.rate
  }
}

rownames(overlap.matrix) <- tissues
colnames(overlap.matrix) <- tissues

my_palette <- colorRampPalette(c("yellow", "red"))(n = 50)
gplots::heatmap.2(overlap.matrix,
          cellnote = overlap.matrix,
          col=my_palette, 
          main = "",
          notecol="black",     
          density.info="none",  
          trace="none",        
          margins =c(10,9),    
          #col=my_palette,      
          #breaks=col_breaks,   
          dendrogram="none",     
          Colv="NA", Rowv = "NA") 

# Plot lineage-specific events at different time points for the same tissue
input.dir <- file.path("data", "processed", "AS_statistics_summary")
output.dir <- file.path("data", "figures", "AS_summary")
tissues <- c("forebrain", "midbrain", "hindbrain", "heart", "limb", "liver", "neuraltube")

AS_summary_file <- read.table(file.path(input.dir, "AS.summary.stat.txt"), 
                              sep = "\t", header = TRUE)
tissue_specific_plot <- list()
for(i in 1:length(tissues)){
  timepoint <- AS_summary_file[AS_summary_file$tissue == tissues[i], ]
  timepoint$lineage.specific.percentage <- round(timepoint$lineage.specific.percentage, 2)
  
  tissue_specific_plot[[i]] <- ggplot2::ggplot(timepoint, aes(x = timepoint, 
                                                              y = lineage.specific.percentage, 
                                                              fill = "blue")) + 
    ggplot2::geom_bar(stat="identity", color="black", position=position_dodge()) + 
    ggplot2::scale_fill_brewer(palette="Paired") + ggplot2::theme_minimal() + 
    geom_text(aes(label = lineage.specific.percentage), vjust= -0.2, color="black",
              position = position_dodge(0.9), size = 3) + ggtitle(tissues[i]) + 
    theme(plot.title = element_text(hjust = 0.5), legend.position="none")
}

pdf(file.path(output.dir, "tissue.specific.plot.pdf"), width = 5, height = 25)
cowplot::plot_grid(plotlist = tissue_specific_plot, ncol = 1)
dev.off()




