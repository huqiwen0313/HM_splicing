# Qiwen Hu - 2017
# This script is used to process the rMAST/rMAST_HM files and classify each alternative spliced exon
# into different categories
# Usage:
#  Rscript SE.files.R

classify_splicing_code <- function(SE_file){
  
  # This function take the input of rMATS (http://rnaseq-mats.sourceforge.net/user_guide.htm#output) 
  # or rMATS HM files and calssify alternative spliced exons into different 
  # categories (splicing codes): gain (0), loss (1), High (2) and Low (3)
  #
  # Args: 
  #  SE_file: rMAT result or rMAST result with HM signal at the end
  #
  # Returns: 
  #  File with a class label indicates splicing codes
  
  inclevel <- colsplit(as.character(SE_file$IncLevel2), split = ",", names = c("s1", "s2"))
  ave_inclevel <- rowMeans(inclevel, na.rm = T)
  SE_file$ave_inclevel <- ave_inclevel
  gain <- SE_file[SE_file$PValue <= 0.05 & SE_file$FDR <= 0.1 
                  & SE_file$IncLevelDifference >= 0.1, ]
  gain <- as.data.frame(gain)
  gain$class <- 0
  loss <- SE_file[SE_file$PValue <= 0.05 & SE_file$FDR <= 0.1 
                  & SE_file$IncLevelDifference <= -0.1, ]
  loss <- as.data.frame(loss)
  loss$class = 1
  High <- SE_file[SE_file$PValue > 0.5 & abs(SE_file$IncLevelDifference) < 0.1
                 & SE_file$ave_inclevel >= quantile(ave_inclevel, 0.75, na.rm = T), ]
  High <- as.data.frame(High)
  High$class <- 2
  Low <- SE_file[SE_file$PValue > 0.5 & abs(SE_file$IncLevelDifference) < 0.1
                & SE_file$ave_inclevel <= quantile(ave_inclevel, 0.25, na.rm = T), ]
  Low <- as.data.frame(Low)
  Low$class <- 3
  out_file <- rbind(gain, loss, High, Low)
  return(out_file)
}

# Processing all rMAST HM results from different timepoints and 
# save them to data/processed/different_timepoints

input.dir <- file.path("data", "raw", "different_timepoints")
output.dir <- file.path("data", "processed", "different_timepoints")

file.list <- list.files(input.dir)
for(f in file.list) {
  HM_SE_file <- read.table(file.path(input.dir, f), sep = "\t", header = TRUE)
  processed_HM_SE_file <- classify_splicing_code(HM_SE_file)
  write.table(processed_HM_SE_file, quote=FALSE, col.names=TRUE, 
              row.names=FALSE,sep="\t", file.path(output.dir, f))
}






