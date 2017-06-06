# KAMBLUP [![](https://img.shields.io/badge/Issues-1%2B-brightgreen.svg)](https://github.com/YinLiLin/R-KAMBLUP/issues) [![](https://img.shields.io/badge/Release-v1.0.1-ff69b4.svg)](https://github.com/YinLiLin/R-KAMBLUP/commits/master)

## Kinship Adjusted Multi-locus Best Linear Unbiased Prediction

## Contents
* [GETTING STARTED](#getting-started)
  - [Installation](#installation)
* [INPUT](#input)
  - [Phenotype file](#phenotype-file)
  - [Covariate file](#covariate-file)
  - [Kinship file](#kinship-file) 
  - [Genotype file](#genotype-file)
    * [Hapmap](#hapmap)
    * [PLINK Binary](#plink-binary)
    * [Numeric](#numeric)
* [USAGE](#usage)
  - [Basic](#basic)
  - [Advanced](#advanced)
* [OUTPUT](#output)
* [FAQ AND HINTS](#faq-and-hints)

---
## GETTING STARTED
`KAMBLUP` is compatible with both [R](https://www.r-project.org/) and [Microsoft R Open](https://mran.microsoft.com/open/), We highly recommend **MRO** instead of **R** for running `KAMBLUP`. **MRO** is the enhanced distribution of **R**, it includes multi-threaded math libraries. These libraries make it possible for so many common R operations, ***such as matrix multiply/inverse, matrix decomposition, and some higher-level matrix operations***, to compute in parallel and use all of the processing power available to [reduce computation times](https://mran.microsoft.com/documents/rro/multithread/#mt-bench).

### Installation
`KAMBLUP` is not available on CRAN, but can be installed using the R package **"devtools"**. There are two packages should be installed beforehand, **"snpStats"** and **"rfunctions"**. `KAMBLUP` can be installed with the following R code:
```r
#if "devtools" isn't installed, please "install.packages(devtools)" first.
devtools::install_version('RcppEigen', version = "0.3.2.9.0")
devtools::install_github("Bioconductor-mirror/snpStats")
devtools::install_github("jaredhuling/rfunctions")
devtools::install_github("YinLiLin/R-KAMBLUP")
```
---
## INPUT
### Phenotype file
> `testPheno.txt`

| Trait1 | Trait2 | Trait3 | ... | TraitN |
| :---: | :---: |  :---: |  :---: |  :---: |
| 3.20 | NA | 1 | ... | -0.31 |
| 0.58 | 58.2 | NA| ... | NA |
| NA | 35.7 | 0 | ... | NA |
| -0.25 | 42.4 | 1| ... | 0.47 |
| NA | NA | 0| ... | -1.01 |
| ... | ... | ... | ... | ... |
| 1.24 | 25.6 | 1 | ... | 0.25 |

### Covariate file
> `testCV.txt`

| female | 168 | 1 | ... | 55 |
| :---: | :---: |  :---: |  :---: |  :---: |
| female | 178 | 2| ... | 57 |
| male | 187 | 2 | ... | 62 |
| male | 156 | 1| ... | 75 |
| male | 148 | 1| ... | 45 |
| ... | ... | ... | ... | ... |
| female | 150 | 1 | ... | 80 |

### Kinship file
> `testKin.txt`

### Genotype file
> `test.geno.desc, test.geno.bin`

#### Hapmap
| rs | alleles | chrom | pos | strand | assembly | center | protLSID | assayLSID | panel | QCcode | Ind1 | Ind2 | Ind3 | Ind4 | Ind5 |
| :---: | :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |  :---: |
| PZB00859.1 | A/C | 1 | 157104 | + | AGPv1 | Panzea | NA | NA | maize282 | NA | CC | CC | CC | CC | AA | 
| PZA01271.1 | C/G | 1 | 1947984 | + | AGPv1 | Panzea | NA | NA | maize282 | NA | CC | GG | CC | GG | CC | 
| PZA03613.2 | G/T | 1 | 2914066 | + | AGPv1 | Panzea | NA | NA | maize282 | NA | GG | GG | GG | GG | GG | 
| PZA03613.1 | A/T | 1 | 2914171 | + | AGPv1 | Panzea | NA | NA | maize282 | NA | TT | TT | TT | TT | TT | 
| PZA03614.2 | A/G | 1 | 2915078 | + | AGPv1 | Panzea | NA | NA | maize282 | NA | GG | GG | GG | GG | GG | 

```r
KAMBLUP.Data(hfile="", out="testGeno")
```

#### PLINK Binary
```r
KAMBLUP.Data(bfile="", out="testGeno")
```

#### Numeric
```r
KAMBLUP.Data(numfile="", mapfile="", out="testGeno")
```

---
## USAGE
### Basic
```r
KAMBLUP(pfile="./testPheno.txt", pheno=1, gfile="./testGeno")
```
```r
KAMBLUP(pfile="./testPheno.txt", pheno=1, gfile="./testGeno", cfile="./testCV.txt", kfile="./testKin.txt")
```
### Advanced
<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-KAMBLUP/master/Manhattan_Fpro.jpg">
<img src="/Manhattan_Fpro.jpg" height="250px" width="800px">
</a>
</p>

---
## OUTPUT

---
## FAQ and Hints

:sos: **Question1:** Failing to install "devtools":

***ERROR: configuration failed for package ‘git2r’***

***removing ‘/Users/acer/R/3.4/library/git2r’***

***ERROR: dependency ‘git2r’ is not available for package ‘devtools’***

***removing ‘/Users/acer/R/3.4/library/devtools’***

:yum: **Answer:** Please type the following codes in terminal.
```ssh
apt-get install libssl-dev/unstable
```
---
:sos: **Question2:** When installing packages from Github with "devtools", there is a error:
 
 ***Error in curl::curl_fetch_disk(url, x$path, handle = handle): Problem with the SSL CA cert (path? access rights?)***
 
:yum: **Answer:** Please type the following codes and than try agian.
```r
library(httr)
set_config(config(ssl_verifypeer = 0L))
```

**Questions, suggestions, and bug reports are welcome and appreciated.** [:arrow_right:](https://github.com/YinLiLin/R-KAMBLUP/issues)

