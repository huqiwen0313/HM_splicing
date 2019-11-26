# functions for data visulization

plotCoverage <- function(region.cov, groups=NULL){
  # This function is used to visualize the chip-seq distribution of flanking regions
  # Input paramters:
  #  region.cov: vector of read coverage in each individual position
  #  groups: group information
  # Output:
  #  plots



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
