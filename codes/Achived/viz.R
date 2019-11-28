# Qiwen Hu - 2017
# This script visulization HM distribution in the flanking regions of exons

library(dplyr)
library(reshape2)
library(ggplot2)
library(gridExtra)

as_ave_chip_signal <- function(HM_file, total_reads){
  
  # This function returns mean HM signal in the flanking regions of alternative spliced exons
  #
  # Args: 
  #  HM_file: processed HM rMAST file with one splice code category
  #  total_reads: total number of aligned reads
  #
  # Returns: 
  #  Mean HM signal in the flanking regions
  
  left_signal <- list()
  right_signal <- list()
  for(i in 1:nrow(HM_file)){
    if(HM_file[i,]$strand == "+"){
      # get chip-seq signal in the flanking region
      left_region <- as.numeric(unlist(strsplit(as.character(HM_file[i,]$chip_left), ",")))
      right_region <- as.numeric(unlist(strsplit(as.character(HM_file[i,]$chip_right), ",")))
    } else {
      # if in minus strand, reverse the flanking region
      left_region <- as.numeric(unlist(strsplit(as.character(HM_file[i,]$chip_left), ",")))
      left_region <- rev(left_region)
      right_region <- as.numeric(unlist(strsplit(as.character(HM_file[i,]$chip_right), ",")))
      right_region <- rev(right_region)
    }
    left_signal[[i]] <- as.data.frame(left_region)
    right_signal[[i]] <- as.data.frame(right_region)
  }
  
  left_region_signal <- dplyr::bind_cols(left_signal)
  right_region_signal <- dplyr::bind_cols(right_signal)
  
  left_region_mean = rowMeans(left_region_signal / total_reads) * 1e6
  right_region_mean = rowMeans(right_region_signal / total_reads) * 1e6
  return(c(left_region_mean, right_region_mean))
}

canonical_ave_chip_signal <- function(HM_file, total_reads){
  left_signal <- list()
  right_signal <- list()
  for(i in 1:nrow(HM_file)){
    if(HM_file[i, 4] == "+"){
      left_reads <- as.numeric(unlist(strsplit(as.character(HM_file[i, 6]), ",")))
      right_reads <- as.numeric(unlist(strsplit(as.character(HM_file[i, 7]), ",")))
      
    } else {
      left_reads <- as.numeric(unlist(strsplit(as.character(HM_file[i, 6]), ",")))
      left_reads <- rev(left_reads)
      right_reads <- as.numeric(unlist(strsplit(as.character(HM_file[i, 7]), ",")))
      right_reads <- rev(right_reads)
    }
    left_signal[[i]] <- as.data.frame(left_reads)
    right_signal[[i]] <- as.data.frame(right_reads)
  }
  left_signal <- dplyr::bind_cols(left_signal)
  right_signal <- dplyr::bind_cols(right_signal)
  
  left_region_mean = rowMeans(left_signal / total_reads) * 1e6
  right_region_mean = rowMeans(right_signal / total_reads) * 1e6
  return(c(left_region_mean, right_region_mean))
}

get_sample_info <- function(sample_file_name){
  sample_info <- unlist(strsplit(sample_file_name, split = "[.]"))
  tissue <- sample_info[1]
  histone_marker <- sample_info[5]
  time_point <- paste(sample_info[3], sample_info[4], sep = ".")
  return(c(tissue, histone_marker, time_point))
}

plot_flanking_HM_distribution <- function(HM_file, total_reads, sample_info){
  
  tissue <- sample_info[1]
  histone_marker <- sample_info[2]
  time_point <- sample_info[3]
  
  gain <- HM_file[HM_file$class == 0, ]
  loss <- HM_file[HM_file$class == 1, ]
  High <- HM_file[HM_file$class == 2, ]
  Low <- HM_file[HM_file$class == 3, ]
  
  gain_ave_chip <- as_ave_chip_signal(gain, total_reads)
  loss_ave_chip <- as_ave_chip_signal(loss, total_reads)
  high_ave_chip <- as_ave_chip_signal(High, total_reads)
  low_ave_chip <- as_ave_chip_signal(Low, total_reads)
  
  chip_all <- data.frame(signal = gain_ave_chip, group = "Gain")
  chip_all <- rbind(chip_all, data.frame(signal = loss_ave_chip, group = "Loss"))
  chip_all <- rbind(chip_all, data.frame(signal = high_ave_chip, group = "high"))
  chip_all <- rbind(chip_all, data.frame(signal = low_ave_chip, group = "low"))
  
  # annova p-value
  mod <- lm(signal ~ group, data = chip_all)
  mod.pval <- signif(anova(mod)$Pr[1], 3)
  
  # subsample canonical exons
  #set.seed(1236)
  #canonical.dir <- file.path("data", "processed", "canonical")
  #canonical <- read.table(file.path(canonical.dir, 
  #                         paste(tissue,histone_marker, "rep1.12.5day.bam.sam.canonical.signal",
  #                         sep=".")), sep = "\t")
  #canonical <- canonical[sample(nrow(canonical), 
  #                        size = max(nrow(gain), nrow(loss), nrow(High), nrow(Low)), replace = F), ]
  #canonical_ave <- canonical_ave_chip_signal(canonical, total_reads)
  
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
  #lines(1:20, canonical_ave[1:20], col="purple", lwd=1, lty = 2)
  #lines(22:41, canonical_ave[21:40], col="purple", lwd=1, lty = 2)
  
  
  legend("bottom", ncol=3, 
         legend = c("inclusion gain", "inclusion loss", "high", "low"), 
         col = c("red", "black", "blue", "green"), 
         bty = "n",lty=1, cex = 0.75, xjust = 0
   )
  
 
}

# model performance bar graph
plot_model_performance <- function(rf_perf, logit_perf){
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
  
  
# main ---------
input.dir <- file.path("data", "processed", "different_timepoints")
output.dir <- file.path("data", "figures", "HM_signal_plot", "new")

# plot average HM signal in flanking regions for different timepoints 
file.list <- list.files(input.dir)
# remove allsample.reads.txt from file.list
file.list <- file.list[-1]
sample_reads <- read.table(file.path(input.dir, "allsample.reads.txt"),
                           sep = "\t", header = FALSE)

for(i in 1:length(file.list)) {
  sample <- gsub(".1.bam.sam.hm.signal", "", file.list[i])
  sample_info <- get_sample_info(sample)
  
  HM_file <- read.table(file.path(input.dir, file.list[i]), sep="\t", header=T)
  
  total_reads <- sample_reads[sample_reads[ ,1] == sample, 2]
  
  if(i == 1) {
    pdf(file.path(output.dir, paste(sample_info[1], sample_info[3], "pdf", sep = ".")), 
        width = 10, height = 13)
    par(mfrow = c(4, 2))
  } 
  
  if(i %% 8 == 1){
    dev.off()
    pdf(file.path(output.dir, paste(sample_info[1], sample_info[3], "pdf", sep = ".")),
        width = 10, height = 13)
    par(mfrow = c(4, 2))
  }
  
  plot_flanking_HM_distribution(HM_file, total_reads, sample_info)
  
} 
dev.off()

###interaction plot##
input.dir <- file.path("data", "model_results", "different_timepoints", "feature_interaction")
output.dir <- file.path("data", "figures", "interaction_plot")

file.list <- list.files(input.dir)
tissues <- c("forebrain", "heart", "hindbrain", "limb", "neuraltube", "liver", "midbrain")
for(i in 1:length(tissues)){
  gain.vs.loss.files <- file.list[grep(paste(tissues[i],".*gain.loss.*", sep = ""), file.list)]
  high.vs.low.files <- file.list[grep(paste(tissues[i],".*high.low.*", sep = ""), file.list)]
  gain.vs.loss.interaction.plot <- list()
  high.vs.low.interaction.plot <- list()

  for(j in 1:length(gain.vs.loss.files)){
    timepoint <- gsub(paste(tissues[i], ".", sep = ""), "", gain.vs.loss.files[j])
    timepoint <- gsub("day.gain.loss.interaction.txt", "", timepoint)
    
    gain.vs.loss <- read.table(file.path(input.dir, gain.vs.loss.files[j]),
                               sep = "\t", header = TRUE)
    #gain.vs.loss$interaction <- gsub("chip_", "", gain.vs.loss$interaction)
    gain.vs.loss.interaction.plot[[j]] <- HM_interaction_plot(gain.vs.loss, 
                                                              tissue = tissues[i], 
                                                              timepoint = timepoint)
  }
  
  for(j in 1:length(high.vs.low.files)){
    timepoint <- gsub(paste(tissues[i], ".", sep = ""), "", high.vs.low.files[j])
    timepoint <- gsub("day.gain.loss.interaction.txt", "", timepoint)
    
    high.vs.low <- read.table(file.path(input.dir, high.vs.low.files[j]),
                               sep = "\t", header = TRUE)
    #high.vs.low$interaction <- gsub("chip_", "", high.vs.low$interaction)
    high.vs.low.interaction.plot[[j]] <- HM_interaction_plot(high.vs.low, 
                                                              tissue = tissues[i], 
                                                              timepoint = timepoint)
  }
  
  # gain v.s. loss plot
  cowplot::plot_grid(plotlist = gain.vs.loss.interaction.plot, ncol = 2)
  ggsave(file.path(output.dir, paste(tissues[i], 
                                     "gain.loss.difftimepoints.interaction.pdf", sep = ".")), 
         width = 20, height = 20)
  # high v.s. low plot
  cowplot::plot_grid(plotlist = high.vs.low.interaction.plot, ncol = 2)
  ggsave(file.path(output.dir, paste(tissues[i], 
                                     "high.low.difftimepoints.interaction.pdf", sep = ".")), 
         width = 20, height = 20)
}