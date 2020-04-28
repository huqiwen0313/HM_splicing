## Specific histone modifications associate with alternative exon selection during mammalian development

**Qiwen Hu, Casey Greene, Elizabeth Heller.** [Specific histone modifications associate with alternative exon selection during mammalian development](https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gkaa248/5823706)

## Summary

Alternative splicing (AS) is frequent during early mouse embryonic development. 
Specific histone post-translational modifications (hPTMs) have been shown to regulate exon splicing by either directly recruiting splice machinery or indirectly modulating transcriptional elongation. 
In this study, we hypothesized that hPTMs regulate expression of alternatively spliced genes for specific processes during differentiation.

To address this notion, we trained supervised machine learning models to relate global hPTM enrichment to AS regulation during mammalian tissue development.
We found that specific histone modifications play a dominant role in skipped exon selection among all the tissues and developmental time points examined. 

In addition, we trained iterative random forest [Github](https://github.com/sumbose/iRF) model to identify interactions of hPTMs that associated with skipped exon selection.
Collectively, our data demonstrated a link between hPTMs and alternative splicing which will drive further experimental studies on the functional relevance of these modifications to alternative splicing.

## Data
Data used in this analysis were downloaded from ENCODE database (https://www.encodeproject.org/). We selected mouse embryonic tissue developmental data from forebrain, hindbrain, midbrain, neural tube, heart, liver and limb from 6 time points (E11.5 - E16.5 day)

## Approach
![overview](https://github.com/huqiwen0313/HM_splicing/blob/master/figs/cover.figure.png)

## Installation

### Requirements
* R >= 3.6.0
* Dependencies: Rscript ./codes/dependencies.R

```r
install.packages("devtools")
devtools::install_github("huqiwen0313/HM_splicing")
library(HMSplicing)
```
## Analysis examples
* [Basic analysis](https://github.com/huqiwen0313/HM_splicing/blob/master/man/basic.analysis.md)

## License 

This repository is dual licensed as **[BSD 3-Clause](https://github.com/huqiwen0313/HM_splicing/blob/master/LICENSE_BSD-3.md)** (source code) and **[CC0 1.0](https://github.com/huqiwen0313/HM_splicing/blob/master/license_CC0.md)** (figures, documentation, and our arrangement of the facts contained in the underlying data)

## Citation
Qiwen Hu, Casey S Greene, Elizabeth A Heller, Specific histone modifications associate with alternative exon selection during mammalian development, Nucleic Acids Research, gkaa248, https://doi.org/10.1093/nar/gkaa248