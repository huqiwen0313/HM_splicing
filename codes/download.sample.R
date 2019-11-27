# This scriput is used to download the bam files from encode database (https://www.encodeproject.org/)

download_bam <- function(file){
  for(i in 1:nrow(file)){
    if(file[i,]$Output.type == "alignments" & file[i,]$File.Status == "released"){
     if(file[i,]$Assay == "RNA-seq" | file[i,]$Assay == "DNase-seq"){
        destiny = paste(file[i,]$Biosample.term.name, file[i,]$Biosample.sex, file[i,]$Biosample.Age,
                        file[i,]$Assay, file[i,]$Assembly, file[i,]$Biological.replicate.s., "bam", sep=".")
       download.file(file[i,]$File.download.URL, destiny, mode="wb", method = "wget", quiet = T)
     } else if(file[i,]$Assay == "ChIP-seq"){
       destiny = paste(file[i,]$Biosample.term.name, file[i,]$Biosample.sex, file[i,]$Biosample.Age,
                       file[i,]$Experiment.target, file[i,]$Assembly, file[i,]$Biological.replicate.s., "bam", sep=".")
       download.file(file[i,]$File.download.URL, destiny, mode="wb", method = "wget", quiet = T)
      }
   }
  }
}

input.dir <- "./data/encode_data_table"
filelist = read.table(file.path(input.dir, "list.txt"), header = F, stringsAsFactors = FALSE)
for(j in 1:nrow(filelist)){
  file = read.table(filelist[j,], sep="\t", header=T, fill=T, stringsAsFactors = FALSE)
  download_bam(file)
}

