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
`KAMBLUP` is compatible with both [R](https://www.r-project.org/) and [Microsoft R Open](https://mran.microsoft.com/open/), We strongly recommend **MRO** instead of **R** for running `KAMBLUP`. **MRO** is the enhanced distribution of **R**, it includes multi-threaded math libraries. These libraries make it possible for so many common R operations, ***such as matrix multiply/inverse, matrix decomposition, and some higher-level matrix operations***, to compute in parallel and use all of the processing power available to [reduce computation times](https://mran.microsoft.com/documents/rro/multithread/#mt-bench).

### Installation
`KAMBLUP` is not available on CRAN, but can be installed using the R package **"devtools"**. There are two packages should be installed beforehand, **"snpStats"** and **"rfunctions"**. `KAMBLUP` can be installed with the following R code:
```r
#if "devtools" isn't installed, please "install.packages(devtools)" first.
devtools::install_version('RcppEigen', version = "0.3.2.9.0")
devtools::install_github("Bioconductor-mirror/snpStats")
devtools::install_github("jaredhuling/rfunctions")
devtools::install_github("YinLiLin/R-KAMBLUP/KAMBLUP")
```
After installed successfully, `KAMBLUP` can be loaded by typing ```library(KAMBLUP)```. Typing `?KAMBLUP` could get the details of all parameters.

---
## INPUT
### Phenotype file

> `mouse.Pheno.txt`

| Trait1 | Trait2 | Trait3 | Trait4 | Trait5 | Trait6 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: |
| 0.224992 | 0.224991 | NA | 1 | NA | -0.285427 |
| -0.974543 | -0.974542 | NA | 0 | NA | -2.333531 |
| 0.195909 | 0.195909 | NA | 1 | NA | 0.046818 |
| NA | NA | NA | NA | NA | NA |
| NA | NA | NA | NA | NA | NA |
| ... | ... | ... | ... | ... | ... |
| NA | NA | NA | NA | NA | 0.720009 |


### Covariate file
> `CV.txt` ***(optional)***

| female | group1 | 1 | ... | 55 |
| :---: | :---: |  :---: |  :---: |  :---: |
| female | group2 | 1| ... | 57 |
| male | group2 | 2 | ... | 62 |
| male | group3 | 2| ... | 75 |
| male | group2 | 2| ... | 45 |
| ... | ... | ... | ... | ... |
| female | group3 | 1 | ... | 80 |

### Kinship file
`KAMBLUP` requires a n×n relatedness matrix. By default, it is automatically calculated using choosed one type of three algorithms(
***"scale","center","vanraden"***
), but can be supplied in file by users. If in this case, the order of individuals for each row and each column in the file must correspond to phenotype file, no column and row names.

> `mouse.Kin.txt` ***(optional)***

| 0.3032 | -0.0193 | 0.0094 | 0.0024 | 0.0381 | ... | -0.0072 |
| :---: | :---: |  :---: |  :---: |  :---: |  :---: |  :---: |
| -0.0193 | 0.274 | -0.0243 | 0.0032 | -0.0081 | ... | 0.0056 |
| 0.0094 | -0.0243 | 0.3207 | -0.0071 | -0.0045 | ... | -0.0407 |
| 0.0024 | 0.0032 | -0.0071 | 0.321 | -0.008 | ... | -0.0093 |
| 0.0381 | -0.0081 | -0.0045 | -0.008 | 0.3498 | ... | -0.0238 |
| ... | ... | ... | ... | ... | ... | ... | 
| -0.0072 | 0.0056 | -0.0407 | -0.0093 | -0.0238 | ... | 0.3436 |

### Genotype file
is a n×n matrix, where each row and each column corresponds to individuals in the same order as in the

> `mouse.map, mouse.geno.desc, mouse.geno.bin`

#### Hapmap

> `mouse.hmp.txt`

| rs# | alleles | chrom | pos | strand | assembly# | center | protLSID | assayLSID | panelLSID | QCcode | A048005080 | A048006063 | A048006555 | A048007096 | A048010273 | ... | A084292044 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| rs3683945 | G/A | 1 | 3197400 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs3707673 | A/G | 1 | 3407393 | + | NA | NA | NA | NA | NA | NA | GA | GA | AA | GA | AA | ... | GG |
| rs6269442 | G/A | 1 | 3492195 | + | NA | NA | NA | NA | NA | NA | AG | GG | GG | AG | GG | ... | AA |
| rs6336442 | G/A | 1 | 3580634 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs13475699 | G | 1 | 3860406 | + | NA | NA | NA | NA | NA | NA | GG | GG | GG | GG | GG | ... | GG |

```r
KAMBLUP.Data(hfile="", out="testGeno")
```

#### PLINK Binary

`mouse.fam`, `mouse.bim`, `mouse.bed`

```r
KAMBLUP.Data(bfile="", out="testGeno")
```

#### Numeric

> `mouse.Numeric.txt`

| 1 | 1 | 2 | 1 | 2 | … | 0 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: | :---: |
| 1 | 1 | 0 | 1 | 0 | … | 2 |
| 1 | 2 | 2 | 1 | 2 | … | 0 |
| 1 | 1 | 2 | 1 | 2 | … | 0 |
| 0 | 0 | 0 | 0 | 0 | … | 0 |

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
<a href="https://raw.githubusercontent.com/YinLiLin/R-KAMBLUP/master/figures/Manhattan_Fpro.jpg">
<img src="/figures/Manhattan_Fpro.jpg" height="250px" width="800px">
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

