# functions for data processing

covsum <- function(rl.table, region, ws=1){
  # caculate cumulative coverage based on a sliding window in a region

  cov.table <-rl.table[[which(names(rl.table) %in% region[1])]]
  region.pos <- seq(region[2], region[3])
  cov <- rep(0, length(region.pos))
  names(cov) <- region.pos

  cov.pos <- intersect(dimnames(cov.table)[[1]], region.pos)
  cov.pos.table <- cov.table[which(as.numeric(dimnames(cov.table)[[1]]) %in% cov.pos)]
  cov[names(cov) %in% cov.pos] <- cov.pos.table
  return(cov)
}

getCoverage <- function(bamfile, region, flanking.ws = 150, read.tag.names=F){
  # caculate read coverage in the flanking region of target exons
  # input: bamfile - bam file of chip-seq data (bowtie alignment)
  #        region - a dataframe indicate chromosome, target exon start and target exon end position
  #        flanking.ws - window size to look at exon flanking region
  # returns: list contain: 1. vector of read coverage of each individual position in the left flanking region
  #                        2. vector of read coverage of each individual position in the right flanking region
  #                        3. total library size of bam file

  if(!is.element("Rsamtools", installed.packages()[, 1])) {
    stop("Rsamtools Bioconductor package is now required for BAM file support. Please install")
  }
  ww <- c("flag","rname","pos","isize","strand","mapq","qwidth"); if(read.tag.names) { ww <- c(ww,"qname") }
  bam <- Rsamtools::scanBam(bamfile, param=Rsamtools::ScanBamParam(what=ww,tag="XC",
                                                                   flag=Rsamtools::scanBamFlag(isUnmappedQuery=FALSE)))[[1]]
  # get read position for each chromosome
  strm <- as.integer(bam$strand=="+")
  rl <- list(chr=tapply(1:length(bam$pos),bam$rname,
                        function(ii) bam$pos[ii]*strm[ii] + (1-strm[ii])*(bam$pos[ii]+bam$qwidth[ii])))
  library.size <- sum(unlist(lapply(rl$chr, length)))

  rl.table <- lapply(rl$chr, table)
  # flanking region of left exons
  exon.left.flanking.start <- region[, 2] - flanking.ws
  exon.left.flanking.end <- region[, 2] + flanking.ws
  left.flanking.region <-data.frame(chr=region$chr, start=exon.left.flanking.start, end=exon.left.flanking.end)
  # flanking region of right exons
  exon.right.flanking.start <- region[, 3] - flanking.ws
  exon.right.flanking.end <- region[, 3] + flanking.ws
  right.flanking.region <-data.frame(chr=region$chr, start=exon.right.flanking.start, end=exon.right.flanking.end)

  cov.left<- apply(left.flanking.region, 1, function(x) covsum(rl.table, region=x))
  cov.right <- apply(right.flanking.region, 1, function(x) covsum(rl.table, region=x))


  return(list(left.flanking=t(cov.left), right.flanking=t(cov.right), size=library.size))
}

