
cal_HM_signal_flanking <- function(HM_file, marker, total_reads){
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
    
    HM_signal[[i]] <- data.frame(class = HM_file[i, ]$class, chip_left_intron = chip_left_intron,
                            chip_left_exon = chip_left_exon, chip_right_intron = chip_right_intron,
                            chip_right_exon = chip_right_exon)
  }
  HM_signal <- dplyr::bind_rows(HM_signal)
  names(HM_signal)[-1] = paste(marker, names(HM_signal)[-1], sep="_")
  return(HM_signal)
}

get_sample_info <- function(sample_file_name){
  sample_info <- unlist(strsplit(sample_file_name, split = "[.]"))
  tissue <- sample_info[1]
  histone_marker <- sample_info[5]
  time_point <- paste(sample_info[3], sample_info[4], sep = ".")
  return(c(tissue, histone_marker, time_point))
}

get_total_reads <- function(sample_file_name){
  input.dir <- file.path("data", "processed", "different_timepoints")
  sample_reads <- read.table(file.path(input.dir, "allsample.reads.txt"),
                             sep = "\t", header = FALSE)
  sample <- gsub(".1.bam.sam.hm.signal", "", sample_file_name)
  total_reads <- sample_reads[sample_reads[ ,1] == sample, 2]
  return(total_reads)
}

# main ---------
input.dir <- file.path("data", "processed", "different_timepoints")
output.dir <- file.path("data", "processed", "HM_features")

file.list <- list.files(input.dir)
#remove allsample.reads.txt from file.list
file.list <- file.list[-1]

time_point <- c("12.5", "13.5", "14.5", "15.5", "16.5")
tissue <- c("forebrain.mixed", "heart.mixed", "hindbrain.mixed", 
            "limb.mixed", "neuraltube.mixed", "liver.mixed", "midbrain.mixed")

unique_file_header <- as.vector(outer(tissue, time_point, paste, sep="."))

unique_file_header <- unique_file_header[-which(unique_file_header == "limb.mixed.16.5" 
                                                | unique_file_header == "neuraltube.mixed.16.5")]

for(file.iter in unique_file_header){
  HM_files <- file.list[grep(file.iter, file.list)]
  HM_flanking_features <- list()
  for(i in 1:length(HM_files)){
    total_reads <- get_total_reads(HM_files[i])
    sample_info <- get_sample_info(HM_files[i])
    HM_file <- read.table(file.path(input.dir, HM_files[i]), sep = "\t", header = TRUE)
    HM_flanking_features[[i]] <- cal_HM_signal_flanking(HM_file, sample_info[2], total_reads)
    if(i > 1){
      HM_flanking_features[[i]] <- HM_flanking_features[[i]][, -1]
    }
  }
  HM_flanking_features <-  dplyr::bind_cols(HM_flanking_features)
  write.table(HM_flanking_features, quote = FALSE, col.names = TRUE, 
              row.names = FALSE, sep="\t", 
              file = file.path(output.dir, paste(sample_info[1], sample_info[3], 
                                                 "HMfeatures.txt", sep = ".")))
}



