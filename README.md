# HM_splicing

### Specific histone modifications associate with alternative exon selection during mammalian development

**Qiwen Hu, Casey Greene, Elizabeth Heller.** [Specific histone modifications associate with alternative exon selection during mammalian development](https://doi.org/10.1101/361816). _bioRxiv._ 2018.

## Summary

Alternative splicing (AS) is frequent during early mouse embryonic development. 
Specific histone post-translational modifications (hPTMs) have been shown to regulate exon splicing by either directly recruiting splice machinery or indirectly modulating transcriptional elongation. 
In this study, we hypothesized that hPTMs regulate expression of alternatively spliced genes for specific processes during differentiation.

To address this notion, we trained supervised machine learning models to relate global hPTM enrichment to AS regulation during mammalian tissue development.
We found that specific histone modifications play a dominant role in skipped exon selection among all the tissues and developmental time points examined. 

In addition, we trained iterative random forest model to identify interactions of hPTMs that associated with skipped exon selection.
Collectively, our data demonstrated a link between hPTMs and alternative splicing which will drive further experimental studies on the functional relevance of these modifications to alternative splicing.

## Data
Data used in this analysis were downloaded from ENCODE database (https://www.encodeproject.org/). We selected mouse embryonic tissue developmental data from forebrain, hindbrain, midbrain, neural tube, heart, liver and limb from 6 time points (E11.5 - E16.5 day)

## Installation

```r
devtools::install_github("huqiwen0313/HM_splicing")
```

## License 

This repository is dual licensed as **[BSD 3-Clause](https://github.com/huqiwen0313/HM_splicing/blob/master/LICENSE_BSD-3.md)** (source code) and **[CC0 1.0](https://github.com/huqiwen0313/HM_splicing/blob/master/license_CC0.md)** (figures, documentation, and our arrangement of the facts contained in the underlying data)

