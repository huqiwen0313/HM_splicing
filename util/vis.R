# Qiwen Hu - 2019
# Core functions for data visulization

library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)
library(core.R)

plot_flanking_HM_distribution <- function(HM_file, total_reads, sample_info){
  # This function is used to visualize the chip-seq distribution of flanking regions
  # Input paramters:
  #  HM_file: files contain the chip-seq signals in the flanking region of AS exons and AS exon class
  #  total_reads: number of total reads in the sample
  #  sample_info: vector contains sample information (tissue, markers and time points)
  # Output: 
  #  plots
  
  # get sample information
  tissue <- sample_info[1]
  histone_marker <- sample_info[2]
  time_point <- sample_info[3]
  
  # Group AS exons into 4 categories
  gain <- HM_file[HM_file$class == 0, ]
  loss <- HM_file[HM_file$class == 1, ]
  High <- HM_file[HM_file$class == 2, ]
  Low <- HM_file[HM_file$class == 3, ]
  
  # Caculate average Chip-seq signals in the flanking regions of AS exons
  gain_ave_chip <- as_ave_chip_signal(gain, total_reads)
  loss_ave_chip <- as_ave_chip_signal(loss, total_reads)
  high_ave_chip <- as_ave_chip_signal(High, total_reads)
  low_ave_chip <- as_ave_chip_signal(Low, total_reads)
  
  chip_all <- data.frame(signal = gain_ave_chip, group = "Gain")
  chip_all <- rbind(chip_all, data.frame(signal = loss_ave_chip, group = "Loss"))
  chip_all <- rbind(chip_all, data.frame(signal = high_ave_chip, group = "high"))
  chip_all <- rbind(chip_all, data.frame(signal = low_ave_chip, group = "low"))
  
  # caculate annova p-value
  mod <- lm(signal ~ group, data = chip_all)
  mod.pval <- signif(anova(mod)$Pr[1], 3)
  
  # subsample canonical exons
  set.seed(1236)
  canonical.dir <- file.path("data", "processed", "canonical")
  canonical <- read.table(file.path(canonical.dir, 
                                    paste(tissue,histone_marker, "rep1.12.5day.bam.sam.canonical.signal",
                                          sep=".")), sep = "\t")
  canonical <- canonical[sample(nrow(canonical), 
                                size = max(nrow(gain), nrow(loss), nrow(High), nrow(Low)), replace = F), ]
  canonical_ave <- canonical_ave_chip_signal(canonical, total_reads)
  
  # plot average chip-seq distribution for alternative spliced exons
  plot(loss_ave_chip[1:20], type="l", xlim = c(0,55), 
       frame.plot = FALSE, ylim=c(0,max(max(gain_ave_chip), max(loss_ave_chip), 
                                        max(high_ave_chip), max(low_ave_chip)) + 0.005), 
       xlab=paste("p-value:", mod.pval, sep=" "), ylab="chip signal", xaxt='n', 
       lwd=2, main = as.character(paste(histone_marker, tissue, time_point, sep = " ")), 
       cex.main=1)
  Axis(side=1, at=c(1, 10, 20, 30, 40, 50), labels=c("-150", "0", "+150", "-150", "0", "+150"))
  lines(30:49, loss_ave_chip[22:41], lwd=2)
  lines(1:20, gain_ave_chip[1:20], col="red", lwd=2)
  lines(30:49, gain_ave_chip[22:41], col="red", lwd=2)
  abline(v = c(10, 40), lty=2, col="grey")
  
  lines(1:20, high_ave_chip[1:20], col="blue", lwd=2)
  lines(30:49, high_ave_chip[21:40], col="blue", lwd=2)
  lines(1:20, low_ave_chip[1:20], col="green", lwd=2)
  lines(30:49, low_ave_chip[21:40], col="green", lwd=2)
  
  # plot average chip-seq distribution for canonical exons
  lines(1:20, canonical_ave[1:20], col="purple", lwd=1, lty = 2)
  lines(30:49, canonical_ave[21:40], col="purple", lwd=1, lty = 2)
  
  
  legend("top", ncol=3, 
         legend = c("inclusion gain", "inclusion loss", "high", "low", "canonical"), 
         col = c("red", "black", "blue", "green", "purple"), 
         bty = "n",lty=c(1, 1, 1, 1, 2), cex = 0.75, xjust = 0
  )
  
}

plot_model_performance <- function(rf_perf, logit_perf){
  # used to visualize the performance of random forest and logistic regression model
  # Input:
  #  rf_perf: performance table of random forest model
  #  logit_perf: performance table of logistic regression
  # Output:
  #  list of bargraphs for model performance
  
  rf_perf <- rf_perf[, c(1:2, 7)]
  rf_perf$method <- "rf"
  logit_perf <- logit_perf[, c(1:2, 7)]
  logit_perf$method <- "logit"
  
  all.perf <- rbind(rf_perf, logit_perf)
  tissues <- unique(all.perf$tissue)
  
  pl.list <- list()
  for(i in 1:length(tissues)){
    perf <- all.perf[all.perf$tissue == tissues[i], ]
    perf$auc = round(perf$auc, 2)
    pl.list[[i]] <- ggplot2::ggplot(perf, aes(x = timepoint, y = auc, fill = method)) + 
      ggplot2::geom_bar(stat="identity", color="black", position=position_dodge()) +
      ggplot2::scale_fill_brewer(palette="Paired") + ggplot2::theme_minimal() + 
      ggplot2::ggtitle(tissues[i]) + ggplot2::theme(plot.title = element_text(hjust = 0.5)) + 
      ggplot2::ylim(c(0, 0.9)) + geom_text(aes(label = auc), vjust= -0.2, color="black",
                                           position = position_dodge(0.9), size = 3)
  }
  return(pl.list)
}

plot_rf_performance <- function(rf_perf){
  # used to plot the performance of random forest model
  # Input:
  #  rf_perf: performance table of random forest model
  # Output:
  #  lists of plots 
  
  rf_perf <- rf_perf[, c(1:3, 7)]
  tissues <- unique(rf_perf$tissue)
  
  pl.list <- list()
  
  for(i in 1:length(tissues)){
    perf <- rf_perf[rf_perf$tissue == tissues[i], ]
    perf$accuracy <- round(perf$accuracy, 2)
    perf$auc <- round(perf$auc, 2)
    perf <- reshape2::melt(perf)
    colnames(perf) <- c("tissue", "timepoint", "measurement", "value")
    pl.list[[i]] <- ggplot2::ggplot(perf, aes(x=timepoint, y=value, group=measurement)) +
      geom_line(aes(color=measurement)) +
      geom_point(aes(color=measurement)) +
      theme(plot.title = element_text(hjust = 0.5), legend.position="top") + 
      ylim(0.5, 0.8) +
      ggplot2::ggtitle(tissues[i])
  } 
  return(pl.list)
}

HM_interaction_plot <- function(interaction.list, cutoff = 0.5, tissue, timepoint){
  # plot interaction of histone markers
  # Input:
  #  interaction.list: table contain interaction information
  #  cutoff: cutoff value for visualization of interactions
  #  tissue: tissue type
  #  timepoint: data timepoint
  # Output:
  #  interaction plot
  
  HM.interaction <- interaction.list[interaction.list$ave >= cutoff, ]
  
  HM.interaction$interaction <- gsub("_chip_left_exon", " 5' D", HM.interaction$interaction)
  HM.interaction$interaction <- gsub("_chip_left_intron", " 5' U", HM.interaction$interaction)
  HM.interaction$interaction <- gsub("_chip_right_intron", " 3' D", HM.interaction$interaction)
  HM.interaction$interaction <- gsub("_chip_right_exon", " 3' U", HM.interaction$interaction)
  
  ggplot2::theme_set(theme_bw())
  plot <- ggplot2::ggplot(HM.interaction, aes(reorder(interaction, ave), ave)) + 
    ggplot2::geom_point(stat='identity', aes(col = "tomato2"), size=6)  +
    ggplot2::geom_text(color="white", size=2, label = round(HM.interaction$ave, 2)) +
    ggplot2::ylim(0, 1) + ggplot2::theme(legend.position="none") + 
    ggplot2::xlab("interactions") + ggplot2::ylab("stability score") + 
    ggplot2::ggtitle(paste(tissue, timepoint)) + ggplot2::coord_flip() +
    ggplot2::theme(plot.title = element_text(hjust = 0.5)) 
  
  return(plot)
}