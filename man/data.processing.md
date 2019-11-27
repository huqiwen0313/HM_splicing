# Analysis pipelines

## Data preprocessing
### Step 1: Data download
Data is downloaded from encode database (https://www.encodeproject.org/). 
We selected mouse embryonic tissue developmental data from forebrain, hindbrain, midbrain, neural tube, heart, liver and limb from 6 time points 
(E11.5 - E16.5 day). 
The detail datasets we used are at data/encode_data_table.

The datasets used in this manuscript can be downloaded using the follow R script:

```R
Rscript ./data/download.sample.R
```

### Step 2: Run rMATs (http://rnaseq-mats.sourceforge.net/) to detect alternative splicing events
For each sample, we used rMATs to detect potential alternative spliced exons. 
All the jobs were ran on PMACS sever at University of Pennsylvania using customized scripts.
Example of rMATs command: 

```R
RNASeq-MATS.py -b1 input_bam_timepoint_1_rep1,input_bam_timepoint_1_rep2 -b2 input_bam_timepoint_2_rep1,input_bam_timepoint_2_rep2 
-gtf gencode.vM4.annotation.gtf -t single -len 100 -c 0.0001 -analysis U -novelSS 1 -keepTemp -o rMATs.out
```
### Step 3: Identify spliced code from rMATs output
Example usage for individual sample:

```R
library(HMSplicing)
rmats.result <- read.table("./example_data/forebrain.12.5.vs.11.5.rmast.out", sep="\t", header=T)
rmatsClassified <- classify_splicing_code(rmats.result)
head(rmatsClassified)
```
```
      ID                GeneID geneSymbol  chr strand exonStart_0base   exonEnd upstreamES upstreamEE downstreamES
6  16886 ENSMUSG00000040797.12     Iqsec3 chr6      -       121378959 121379909  121372932  121376531    121379994
10 11233 ENSMUSG00000038671.11     Arfrp1 chr2      -       181361000 181361066  181357832  181359561    181361390
12 11235 ENSMUSG00000038671.11     Arfrp1 chr2      -       181361000 181361069  181359314  181359561    181361390
13 12590  ENSMUSG00000059974.6        Ntm chr9      -        29006253  29006289   29000704   29000737     29009225
15 16697  ENSMUSG00000037211.8      Spry1 chr3      +        37640519  37640808   37639957   37640105     37642554
16 30964  ENSMUSG00000059742.6      Kcnh7 chr2      -        62845314  62845338   62837065   62837280     62850349
   downstreamEE  ID.1 IC_SAMPLE_1 SC_SAMPLE_1 IC_SAMPLE_2 SC_SAMPLE_2 IncFormLen SkipFormLen       PValue          FDR
6     121380044 16886         4,1         0,0        5,12         8,2       1049          99 5.108369e-11 3.184898e-07
10    181361472 11233         3,1         0,0         0,0         2,3        165          99 4.071635e-10 1.523117e-06
12    181361472 11235         2,1         0,0         0,0         2,3        168          99 1.787005e-09 4.774877e-06
13     29009377 12590         1,1         0,0         3,0         6,3        135          99 1.786725e-09 4.774877e-06
15     37644598 16697        8,17       24,16         0,0       35,35        388          99 2.081150e-09 5.190110e-06
16     62850778 30964         1,1         0,0         0,0        2,12        123          99 2.467113e-09 5.768110e-06
     IncLevel1   IncLevel2 IncLevelDifference ave_inclevel class
6      1.0,1.0 0.056,0.362              0.791        0.209     0
10     1.0,1.0     0.0,0.0              1.000        0.000     0
12     1.0,1.0     0.0,0.0              1.000        0.000     0
13     1.0,1.0   0.268,0.0              0.866        0.134     0
15 0.078,0.213     0.0,0.0              0.145        0.000     0
16     1.0,1.0     0.0,0.0              1.000        0.000     0
```



