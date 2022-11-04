

## Overview

IFAA is a novel approach to make inference on the association of covariates with the absolute abundance (AA) of microbiome in an ecosystem. 

## Installation
```r
# install from GitHub:
devtools::install_github("gitlzg/IFAA")
```

```r
# install from Bioconductor:
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("IFAA")
```r
## Usage

Use sample datasets to run `IFAA()` function.
```r
# Detailed instructions on the package are provided in the manual and vignette
library(IFAA)
library(SummarizedExperiment)

data(dataM)
dim(dataM)
dataM[1:5, 1:8]

data(dataC)
dim(dataC)
dataC[1:3, ]

## Merge microbiome data and covariate data by id, to avoid unmatching observations.
data_merged<-merge(dataM,dataC,by="id",all=FALSE)

## Seperate microbiome data and covariate data, drop id variable from microbiome data
dataM_sub<-data_merged[,colnames(dataM)[!colnames(dataM)%in%c("id")]]
dataC_sub<-data_merged[,colnames(dataC)]

## Create SummarizedExperiment object
test_dat<-SummarizedExperiment(assays=list(MicrobData=t(dataM_sub)), colData=dataC_sub)

## If you already have a SummarizedExperiment format data, you can
## ignore the above steps.


results <- IFAA(experiment_dat = test_dat,
                testCov = c("v1"),
                ctrlCov = c("v2","v3"),
                fdrRate = 0.05)

```


Once the analysis is done, you can extract the regression coefficients along with 95% confidence intervals using this command:
```r
summary_res<-results$full_results
```



Use sample datasets to run `MZILN()` function.
```r
results <- MZILN(experiment_dat=test_dat,
                 targetTaxa = "rawCount18",
                 refTaxa=c("rawCount11"),
                 allCov=c("v1","v2","v3"),
                 fdrRate=0.15)
                 ```
Regression results including confidence intervals can be extracted in the following way:
```r
results$full_results
```

## References 
- Zhigang Li, Lu Tian, A. James O'Malley, Margaret R. Karagas, Anne G. Hoen, Brock C. Christensen, Juliette C. Madan, Quran Wu, Raad Z. Gharaibeh, Christian Jobin, Hongzhe Li (2020) IFAA: Robust association identification and Inference For Absolute Abundance in microbiome analyses. arXiv:1909.10101v3

- Zhigang Li, Katherine Lee, Margaret Karagas, Juliette Madan, Anne Hoen, James O’Malley and Hongzhe Li (2018 ) Conditional regression based on a multivariate zero-inflated logistic normal model for modeling microbiome data. Statistics in Biosciences  10(3):587-608
