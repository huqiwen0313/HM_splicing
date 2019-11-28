## Data exploration and modelling

### Plot distributions of hPTMs in the exon flanking region

Below is an example that we plot the distribution of H3K4me1 in the exon flanking regions
for different exon splicing patterns. In this example, signals from constitutive exons were not sampled.

```R
hPTM <- read.table("./example_data/forebrain.mixed.12.5day.H3K4me1.1.bam.sam.hm.signal", sep="\t", header=TRUE)
plotCoverage(hPTM, total_reads=22497119, c("forebrain", "H3K4me1", "12.5 day"), subsampleCanonical=F, CanonicalFile=NULL)
```
<img src="../figs/H3K4me1.hm.distr.png" width="509" height="362" /> 

Constitutive exons can be used as a control to see if the hPTM pattern is enriched in alternative exons, 
so we also provide an option to sample signals from constitutive exons. 
In this case, a file contains all constitutive exons and their correspondent hPTM signals 
in the flanking regions is needed (example file at ./example_data/H3K4me1.canonical.exon.signal).

```R
plotCoverage(hPTM, total_reads=22497119, c("forebrain", "H3K4me1", "12.5 day"), subsampleCanonical=T, CanonicalFile="./example_data/H3K4me1.canonical.exon.signal")
```
<img src="../figs/H3K4me1.sampled.distr.png" width="509" height="362" /> 

### Run plain logistic and random forest model

```R
allFeatures <- read.table("forebrain.12.5day.HMfeatures.txt", sep="\t", header=TRUE)
gain.loss.features <- allFeatures[allFeatures$class==0 | allFeatures$class==1, ]
gain.loss.modelResult <- model(gain.loss.features)
str(gain.loss.modelResult)
```
It will generate the performance measurements for logistic regression and random forest model and
provide information for important markers.

```
List of 3
 $ rf.performance   :'data.frame':	6 obs. of  5 variables:
  ..$ accuracy : num [1:6] 0.576 0.608 0.544 0.56 0.584 ...
  ..$ Precision: num [1:6] 0.723 0.606 0.552 0.716 0.75 ...
  ..$ recall   : num [1:6] 0.573 0.672 0.578 0.571 0.593 ...
  ..$ f1       : num [1:6] 0.639 0.637 0.565 0.636 0.662 ...
  ..$ auc      : num [1:6] 0.692 0.643 0.585 0.649 0.631 ...
 $ logit.performance:'data.frame':	6 obs. of  5 variables:
  ..$ accuracy : num [1:6] 0.48 0.544 0.576 0.536 0.544 0.536
  ..$ Precision: num [1:6] 0.723 0.775 0.657 1 1 ...
  ..$ recall   : num [1:6] 0.5 0.573 0.595 0.536 0.544 ...
  ..$ f1       : num [1:6] 0.591 0.659 0.624 0.698 0.705 ...
  ..$ auc      : num [1:6] 0.536 0.584 0.578 0.5 0.5 ...
 $ impscore         :'data.frame':	32 obs. of  7 variables:
  ..$ variable: chr [1:32] "H3K27ac_chip_left_exon" "H3K27ac_chip_left_intron" "H3K27ac_chip_right_exon" "H3K27ac_chip_right_intron" ...
  ..$ 1       : num [1:32] 7.68 5.92 6 8.8 9.97 ...
  ..$ 2       : num [1:32] 7.53 6.49 6.75 9.27 8.21 ...
  ..$ 3       : num [1:32] 8.55 5.02 5.29 10.55 9.33 ...
  ..$ 4       : num [1:32] 6.62 7.43 7.35 8.51 8.32 ...
  ..$ 5       : num [1:32] 7.25 6.27 6.84 10.17 8.37 ...
  ..$ ave     : num [1:32] 7.53 6.23 6.44 9.46 8.84 ...
```

 