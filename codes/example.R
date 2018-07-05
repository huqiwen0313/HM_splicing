
### look at correlation between HM and psi for some example genes

total_chip_signal <- function(HM_file, total_reads = 1){
  HM_file$total_signal <- 0
  for(i in 1:nrow(HM_file)){
    chip_left <- as.numeric(unlist(strsplit(as.character(HM_file[i, ]$chip_left), ",")))
    chip_right <- as.numeric(unlist(strsplit(as.character(HM_file[i, ]$chip_right), ",")))
    total_signal <- sum(chip_left) + sum(chip_right)
    HM_file[i, ]$total_signal <- total_signal
    
  }
  HM_file <- HM_file[,c(3:7,23,26)]
  return(HM_file)
}

sample_reads <- function(HM, timepoints, tissue = "forebrain"){
  input.dir <- file.path("data", "examples")
  allsample.reads <- read.table(file.path(input.dir,"allsample.reads.txt"), sep = "\t", header = FALSE)
  reads <- c()
  for(i in 1:length(timepoints)){
    sample <- paste(tissue, "mixed", timepoints[i], HM, sep = ".")
    read <- allsample.reads[allsample.reads[, 1] == sample, 2]
    reads <- c(reads, read)
  }
  return(reads)
}


# select AS examples
input.dir <- file.path("data", "processed", "AS_statistics_summary", "gain.vs.loss.summary")
AS.example <- read.table(file.path(input.dir, "gene.examples.txt"), header = FALSE, sep = "\t")
AS.example <- AS.example[, c(1, 4)]
colnames(AS.example) <- c("genename", "length")
tissue <- c("forebrain", "heart", "hindbrain", "limb", "liver", "midbrain", "neuraltube")
for(i in 1:length(tissue)){
  AS.file <- read.table(file.path(input.dir, paste(tissue[i], "all.as.txt", sep = ".")), 
                                  sep = "\t", header = TRUE)
  AS.file$length <- AS.file[, 6] - AS.file[, 5]
  colnames(AS.file)[2] <- "genename"
  AS.file <- merge(AS.file, AS.example, by=c("genename", "length"))
  write.table(AS.file, quote = FALSE, col.names = TRUE, 
              row.names = FALSE, sep="\t", 
              file = file.path(input.dir, paste(tissue[i], "as.example.txt"), sep = "."))
}


# get chip-seq count in flanking regions
input.dir <- file.path("data", "raw", "different_timepoints")
file.list <- list.files(input.dir)

time_point <- c("12.5", "13.5", "14.5", "15.5", "16.5")
tissue <- c("forebrain.mixed", "heart.mixed", "hindbrain.mixed", 
            "limb.mixed", "neuraltube.mixed", "liver.mixed", "midbrain.mixed")
all.signal <- list()
for(i in 1:length(time_point)){
  file <- file.list[grep(paste(tissue, time_point[i], sep = "."), file.list)]
  HM.signal <- list()
  for( j in 1:length(file)){
    sample_info <- get_sample_info(file[j])
    HM_file <- read.table(file.path(input.dir, file[j]), sep = "\t", header = T)
    HM_file <- total_chip_signal(HM_file)
    colnames(HM_file)[6] = sample_info[2]
    if(j == 1){
      HM.signal[[j]] <- HM_file
      HM.signal[[j]]$timepoint <- sample_info[3]
    }
    else{
      HM.signal[[j]] <- data.frame(HM_file[6])
    }
  }
  all.signal[[i]] <- dplyr::bind_cols(HM.signal)
}
all.signal <- dplyr::bind_rows(all.signal)

output.dir <- file.path("data", "examples")
write.table(all.signal, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep="\t", 
            file = file.path(output.dir, "forebrain.all.chip.signal.txt"))

###
input.dir <- file.path("data", "examples")
example.gene <- read.table(file.path(input.dir, "forebrain.as.example.txt"), sep = "\t", header = TRUE)
chip <- read.table(file.path(input.dir, "forebrain.all.chip.signal.txt"), sep = "\t", header = TRUE)
chip$geneSymbol <- as.character(chip$geneSymbol)
example.gene$genename <- as.character(example.gene$genename)
hm.cor.all <- c()

for(i in 1:nrow(example.gene)){
  gene.chip <- chip[chip$geneSymbol == example.gene[i,]$genename & 
                      chip$exonStart_0base == example.gene[i,]$V5 & chip$exonEnd == example.gene[i,]$V6, ]
  gene.hm <- unique(gene.chip[, c(6,8:14)])
  hm.cor <- c()
  for(j in 1:ncol(gene.hm)){
    reads <- sample_reads(colnames(gene.hm)[j], unique(gene.chip$timepoint))
    HM <- gene.hm[, j]/reads * 1e6
    psi <- -1 * as.numeric(example.gene[i, 12:ncol(example.gene)])
    cor <- cor(HM, psi)
    hm.cor <- c(hm.cor, cor)
    
    if(abs(cor) > 0.5){
      gene <- data.frame(HM = HM, psi = psi, timepoint = unique(gene.chip$timepoint))
      
      ggplot2::ggplot(gene, aes(x=HM, y=psi, color=timepoint)) +
        geom_point(size=3) + geom_smooth(method=lm, se=FALSE, color="black") + 
        ggplot2::theme_classic() + 
        ggplot2::xlab(paste(example.gene[i, ]$genename, colnames(gene.hm[j]))) +
        annotate(geom="text", x=0.15, y=0.7, 
                 label=paste("cor =", round(cor, 2)),
                 color="red")
      ggsave(file.path(input.dir, paste(example.gene[i, ]$genename,
                                        colnames(gene.hm[j]), "eps", sep = ".")))
    }
  }
  names(hm.cor) <- colnames(gene.hm)
  hm.cor.all <- rbind(hm.cor.all, hm.cor)
}
hm.gene <- cbind(example.gene, hm.cor.all)
write.table(hm.gene, quote = FALSE, col.names = TRUE, 
            row.names = FALSE, sep="\t", 
            file = file.path(input.dir, "forebrain.example.gene.cor.txt"))

####
flna <- all.signal[all.signal$geneSymbol == "Flna" & all.signal$exonStart_0base == 74240815, ]
flna$psi <- abs(c(0,-0.126,-0.146,-0.232,-0.313))
sample_reads_k36 <- c(45224195, 50466837, 25907176, 69516474, 20400228)
sample_reads_k4 <- c(54178271, 52930402, 31267709, 53292662, 23026590)
sample_reads_k9 <- c(21064618, 38094512, 21068187, 28206586, 19902788)
flna$H3K36me3 <- flna$H3K36me3/sample_reads_k36 * 1e6
flna.k36 <- flna[,c(7,9,15)]

flna$H3K9ac <- flna$H3K9ac/sample_reads_k9 * 1e6
flna.k9 <- flna[,c(7,13,15)]

flna$H3K4me1 <- flna$H3K4me1/sample_reads_k4 * 1e6
flna.k4 <- flna[,c(7,10,15)]

ggplot2::ggplot(flna.k9, aes(x=psi, y=H3K9ac, color=timepoint)) +
  geom_point(size=3) + geom_smooth(method=lm, se=FALSE, color="black") + 
  ggplot2::theme_classic() + 
  annotate(geom="text", x=0.05, y=0.5, 
           label=paste("cor =", round(cor(flna$psi, flna$H3K9ac), 2)),
           color="red")
