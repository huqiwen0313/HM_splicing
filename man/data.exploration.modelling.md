## Data exploration and modelling

### Plot distributions of hPTMs in the exon flanking region

Below is an example that we plot the distribution of H3K4me1 in the exon flanking regions
for different exon splicing patterns. In this example, signals from constitutive exons were not sampled.

```R
hPTM <- read.table("./example_data/forebrain.mixed.12.5day.H3K4me1.1.bam.sam.hm.signal", sep="\t", header=TRUE)
plotCoverage(hPTM, 22497119, c("forebrain", "H3K4me1", "12.5 day"), subsampleCanonical=F, CanonicalFile=NULL)
```
<img src="../figs/H3K4me1.hm.distr.png" width="509" height="362" /> 


Constitutive exons can be used as a control to see if the hPTM pattern is enriched in alternative exons, 
so we also provide an option to sample signals from constitutive exons. 
In this case, a file contains all constitutive exons and their correspondent hPTM signals 
in the flanking regions is needed (example file at ./example_data/H3K4me1.canonical.exon.signal).

```R
plotCoverage(hPTM, 22497119, c("forebrain", "H3K4me1", "12.5 day"), subsampleCanonical=T, CanonicalFile="./example_data/H3K4me1.canonical.exon.signal")
```
<img src="../figs/H3K4me1.sampled.distr.png" width="509" height="362" /> 

