# KAML [![](https://img.shields.io/badge/Issues-1%2B-brightgreen.svg)](https://github.com/YinLiLin/R-KAML/issues) [![](https://img.shields.io/badge/Release-v1.0.1-important.svg)](https://github.com/YinLiLin/R-KAML/commits/master) [![](https://img.shields.io/badge/license-GPL-blue.svg)](https://github.com/YinLiLin/KAML/blob/master/LICENSE)

## *[K](https://github.com/YinLiLin/R-KAML)inship [A](https://github.com/YinLiLin/R-KAML)djusted [M](https://github.com/YinLiLin/R-KAML)ultiple [L](https://github.com/YinLiLin/R-KAML)oci Best Linear Unbiased Prediction*

## Contents

* [OVERVIEW](#overview)
* [CITATION](#citation)
* [GETTING STARTED](#getting-started)<img src="https://raw.githubusercontent.com/YinLiLin/R-KAML/master/figures/KAML_log.png" height="250" align="right" />
  - [Installation](#installation)
  - [Test Datasets](#test-datasets)
* [INPUT](#input)
  - [Phenotype file](#phenotype-file)/[Covariate file](#covariate-file)/[Kinship file](#kinship-file) 
  - [Genotype file](#genotype-file)
    * [Hapmap](#hapmap)/[VCF](#vcf)/[PLINK Binary](#plink-binary)/[Numeric](#numeric)
* [USAGE](#usage)
  - [Basic](#basic)
  - [Advanced](#advanced)
* [OUTPUT](#output)
* [FAQ AND HINTS](#faq-and-hints)

---
## OVERVIEW
***`KAML`*** is designed to predict genetic values using genome-wide or chromosome-wide SNPs for either simple traits that controlled by a limited number of major genes or complex traits that influenced by many polygenes with minor effects. In brief, ***`KAML`*** provides a flexible assumption to accommodate traits of various genetic architectures and incorporates pseudo-QTNs as fixed effect terms and a trait-specific random effect term under the LMM framework. The model parameters are optimized using the information of bootstrap strategy based GWAS results in a parallel accelerated machine learning procedure combing cross-validation, grid search and bisection algorithms.

<!-- <p align="center">
<a target="_blank" href="https://camo.githubusercontent.com/6af206ddeaf9b5fdd9b944e82e7b5d698f29ccbf/687474703a2f2f6c617465782e636f6465636f67732e636f6d2f6769662e6c617465783f253742793d58622b253742253543636f6c6f7225374252656425374451712537442b5a253742253543636f6c6f722537425265642537442535436d752537442b653b2535436d7525354373696d2537424e25374428302c253742253543636f6c6f722537425265642537444b2537442535437369676d612535452537423225374429253744"><img src="https://camo.githubusercontent.com/6af206ddeaf9b5fdd9b944e82e7b5d698f29ccbf/687474703a2f2f6c617465782e636f6465636f67732e636f6d2f6769662e6c617465783f253742793d58622b253742253543636f6c6f7225374252656425374451712537442b5a253742253543636f6c6f722537425265642537442535436d752537442b653b2535436d7525354373696d2537424e25374428302c253742253543636f6c6f722537425265642537444b2537442535437369676d612535452537423225374429253744" alt="equation" data-canonical-src="http://latex.codecogs.com/gif.latex?%7By=Xb+%7B%5Ccolor%7BRed%7DQq%7D+Z%7B%5Ccolor%7BRed%7D%5Cmu%7D+e;%5Cmu%5Csim%7BN%7D(0,%7B%5Ccolor%7BRed%7DK%7D%5Csigma%5E%7B2%7D)%7D" style="max-width:100%;">
</a>
</p>

<p align="center">
<a target="_blank" href="https://camo.githubusercontent.com/9f1a2f99488778d874d22d6720f1a5bd946cccf1/687474703a2f2f6c617465782e636f6465636f67732e636f6d2f6769662e6c617465783f2537424b5f253742696a2537443d25354366726163253742312537442537426d25374425354373756d5f2537426b3d312537442535452537426d25374425354366726163253742284d5f253742696b2537442d32705f2537426b25374429253742253543636f6c6f7225374252656425374425354378695f2537426b253744253744284d5f2537426a6b2537442d32705f2537426b2537442925374425374232705f2537426b25374428312d705f2537426b25374429253744253744"><img src="https://camo.githubusercontent.com/9f1a2f99488778d874d22d6720f1a5bd946cccf1/687474703a2f2f6c617465782e636f6465636f67732e636f6d2f6769662e6c617465783f2537424b5f253742696a2537443d25354366726163253742312537442537426d25374425354373756d5f2537426b3d312537442535452537426d25374425354366726163253742284d5f253742696b2537442d32705f2537426b25374429253742253543636f6c6f7225374252656425374425354378695f2537426b253744253744284d5f2537426a6b2537442d32705f2537426b2537442925374425374232705f2537426b25374428312d705f2537426b25374429253744253744" alt="equation" data-canonical-src="http://latex.codecogs.com/gif.latex?%7BK_%7Bij%7D=%5Cfrac%7B1%7D%7Bm%7D%5Csum_%7Bk=1%7D%5E%7Bm%7D%5Cfrac%7B(M_%7Bik%7D-2p_%7Bk%7D)%7B%5Ccolor%7BRed%7D%5Cxi_%7Bk%7D%7D(M_%7Bjk%7D-2p_%7Bk%7D)%7D%7B2p_%7Bk%7D(1-p_%7Bk%7D)%7D%7D" style="max-width:100%;">
</a>
</p>

<p align="center">
<a target="_blank" href="https://camo.githubusercontent.com/844a26b674cc2acbc1f0c3ef08dd55c074ffb475/687474703a2f2f6c617465782e636f6465636f67732e636f6d2f6769662e6c617465783f25374225354378695f2537426b2537442535436d696428253543616c7068612c253543626574612925354373696d253543626567696e25374263617365732537443126253543746578742537423b253744312d25354362657461264c6f675f253742253543616c70686125374428505f2537426d253543626574612b31253744292d4c6f675f253742253543616c70686125374428505f2537426d253543626574612537442926253543746578742537423b25374425354362657461253543656e642537426361736573253744253744"><img src="https://camo.githubusercontent.com/844a26b674cc2acbc1f0c3ef08dd55c074ffb475/687474703a2f2f6c617465782e636f6465636f67732e636f6d2f6769662e6c617465783f25374225354378695f2537426b2537442535436d696428253543616c7068612c253543626574612925354373696d253543626567696e25374263617365732537443126253543746578742537423b253744312d25354362657461264c6f675f253742253543616c70686125374428505f2537426d253543626574612b31253744292d4c6f675f253742253543616c70686125374428505f2537426d253543626574612537442926253543746578742537423b25374425354362657461253543656e642537426361736573253744253744" alt="equation" data-canonical-src="http://latex.codecogs.com/gif.latex?%7B%5Cxi_%7Bk%7D%5Cmid(%5Calpha,%5Cbeta)%5Csim%5Cbegin%7Bcases%7D1&amp;%5Ctext%7B;%7D1-%5Cbeta&amp;Log_%7B%5Calpha%7D(P_%7Bm%5Cbeta+1%7D)-Log_%7B%5Calpha%7D(P_%7Bm%5Cbeta%7D)&amp;%5Ctext%7B;%7D%5Cbeta%5Cend%7Bcases%7D%7D" style="max-width:100%;">
</a>
</p>

 -->

***`KAML`*** is developed by [***Lilin Yin***](https://github.com/YinLiLin), [***Haohao Zhang***](https://github.com/hyacz), and [***Xiaolei Liu***](https://github.com/XiaoleiLiuBio)**\*** at the [***Huazhong(Central China) Agricultural University***](http://www.hzau.edu.cn/en/HOME.htm).

If you have any bug reports or questions please send an email to **Xiaolei Liu** at **xiaoleiliu@mail.hzau.edu.cn**

---
## CITATION
(waiting for update)

---
## GETTING STARTED
***`KAML`*** is compatible with both [R](https://www.r-project.org/) and [Microsoft R Open](https://mran.microsoft.com/open/), **WE STRONGLY RECOMMEND TO INSTALL ***`KAML`*** ON Microsoft R Open(https://mran.microsoft.com/download/)**. **MRO** is the enhanced distribution of **R** from Microsoft Corporation, and it includes the state-of-the-art parallel accelerated math libraries. Those libraries would [reduce the time consumption significantly](https://mran.microsoft.com/documents/rro/multithread/#mt-bench) as many matrix operations could be computed in parallel by using all the available processing power.

### Installation
***`KAML`*** can be installed with **"devtools"** by using the following R codes:
```r
#if "devtools" isn't installed, please "install.packages(devtools)" first.
devtools::install_github("YinLiLin/R-KAML")
```
If you get trouble in installing the **"devtools"**, please download [KAML_1.0.1.tar.gz](https://github.com/YinLiLin/KAML/releases/download/1.0.1/KAML_1.0.1.tar.gz), and then try the following steps:
```r
pkg <- setdiff(c("RcppEigen", "bigmemory"), installed.packages()[,c("Package")])
install.packages(pkg)
install.packages("KAML_1.0.1.tar.gz", repos=NULL)
```
After installed successfully, the ***`KAML`*** package can be loaded by typing
```r
library(KAML)
```
Typing `?KAML` could get the details of all parameters. Some related functions for ***`KAML`*** could be refered to our developed GWAS package rMVP (https://github.com/XiaoleiLiuBio/rMVP), and integrated analysis extending ***`KAML`*** to BLUP could be found at our flexible GP/GS package HIBLUP (https://hiblup.github.io).

### Test Datasets

The example data can be downloaded by typing:
```bash
wget https://github.com/YinLiLin/KAML/releases/download/1.0.1/Data_example.zip
unzip example.zip
```

**Or** by clicking the [example.zip](https://github.com/YinLiLin/KAML/releases/download/1.0.1/example.zip) in your browser. After downloaded, unzip the file and change the workplace by ``` setwd("")```.

---

## INPUT
### Phenotype file
The file must contain a header row, which may represents the trait names. The missing values should be denoted by NA, which will be treated as candidates. Notice that only the numeric values are allowed and the characters will not be recognized. However, if a phenotype takes only values of 0, 1 (or only two levels), ***`KAML`*** would consider it to be a case-control trait, and the predicted value could be directly interpreted as the probability of being a case. <br>

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
The **Covariate file** is **optional**, in order to fit the model for raw phenotype prediction, users can provide the covariates. Please attention that NAs are not allowed in the **Covariate file**,  and all individuals should be in the same order with phenotype file. In ***`KAML`***, there are two types of covariates: **`dcovfile`** and **`qcovfile`**.

**`dcovfile`*****(optional)***: discrete covariates, e.g. *dcov.txt*. Each discrete covariate is recognized as a categorical factor with several levels. Each level can be denoted by a single character, a word, or a numerical number. 

**`qcovfile`*****(optional)***: quantitative covariates, e.g. *qcov.txt*. Each quantitative covariate is recognized as a continuous variable.


***NOTE:*** the design matrix of the mean, which is a vector of all ones, in the model is always a linear combination of the design matrix of a discrete covariate so that not all the effects of the levels (or classes, e.g. male and female) of a discrete covariate are estimable. ***`KAML`*** will always constrain the effect of the first level to be zero and the effect of any other level represents its difference in effect compared to the first level.

<table>
<tbody>
<tr>
<td align="center"><em><strong><code>dcov.txt</code></strong></em></td>
<td align="center"><em><strong><code>qcov.txt</code></strong></em></td>
</tr>
<tr>
<td align="center">

<table>
<tbody>
<tr>
<td align="center">F</td>
<td align="center">H1</td>
<td align="center">1</td>
<td align="center">pop1</td>
</tr>
<tr>
<td align="center">F</td>
<td align="center">H3</td>
<td align="center">0</td>
<td align="center">pop2</td>
</tr>
<tr>
<td align="center">M</td>
<td align="center">H2</td>
<td align="center">0</td>
<td align="center">pop2</td>
</tr>
<tr>
<td align="center">F</td>
<td align="center">H3</td>
<td align="center">1</td>
<td align="center">pop4</td>
</tr>
<tr>
<td align="center">M</td>
<td align="center">H3</td>
<td align="center">1</td>
<td align="center">pop4</td>
</tr>
<tr>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
</tr>
<tr>
<td align="center">M</td>
<td align="center">H5</td>
<td align="center">0</td>
<td align="center">pop5</td>
</tr></tbody></table>

</td>

<td align="center">
<table>
<tbody>
<tr>
<td align="center">12</td>
<td align="center">0.01</td>
<td align="center">0.13</td>
</tr>
<tr>
<td align="center">5</td>
<td align="center">-0.05</td>
<td align="center">0.25</td>
</tr>
<tr>
<td align="center">7</td>
<td align="center">0.05</td>
<td align="center">-0.36</td>
</tr>
<tr>
<td align="center">13</td>
<td align="center">0.16</td>
<td align="center">0.28</td>
</tr>
<tr>
<td align="center">2</td>
<td align="center">0.07</td>
<td align="center">0.95</td>
</tr>
<tr>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
</tr>
<tr>
<td align="center">10</td>
<td align="center">-0.12</td>
<td align="center">0.35</td>
</tr></tbody></table>

</td>
</tr></tbody></table>


### Kinship file 
***`KAML`*** requires a n×n matrix that represents the relationship among individuals. By default, it could be automatically calculated by using one of the three algorithms (***"scale","center","vanraden"***) that implemented in ***`KAML`*** package. It could be also supplied by the users, however, in this case, the order of individuals in either row or column should be the same as phenotype file, the column and row names are not needed.

> `mouse.Kin.txt` ***(optional)***

<table>
<tbody>
<tr>
<td align="center">0.3032</td>
<td align="center">-0.0193</td>
<td align="center">0.0094</td>
<td align="center">0.0024</td>
<td align="center">0.0381</td>
<td align="center">...</td>
<td align="center">-0.0072</td>
</tr>
<tr>
<td align="center">-0.0193</td>
<td align="center">0.274</td>
<td align="center">-0.0243</td>
<td align="center">0.0032</td>
<td align="center">-0.0081</td>
<td align="center">...</td>
<td align="center">0.0056</td>
</tr>
<tr>
<td align="center">0.0094</td>
<td align="center">-0.0243</td>
<td align="center">0.3207</td>
<td align="center">-0.0071</td>
<td align="center">-0.0045</td>
<td align="center">...</td>
<td align="center">-0.0407</td>
</tr>
<tr>
<td align="center">0.0024</td>
<td align="center">0.0032</td>
<td align="center">-0.0071</td>
<td align="center">0.321</td>
<td align="center">-0.008</td>
<td align="center">...</td>
<td align="center">-0.0093</td>
</tr>
<tr>
<td align="center">0.0381</td>
<td align="center">-0.0081</td>
<td align="center">-0.0045</td>
<td align="center">-0.008</td>
<td align="center">0.3498</td>
<td align="center">...</td>
<td align="center">-0.0238</td>
</tr>
<tr>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
<td align="center">...</td>
</tr>
<tr>
<td align="center">-0.0072</td>
<td align="center">0.0056</td>
<td align="center">-0.0407</td>
<td align="center">-0.0093</td>
<td align="center">-0.0238</td>
<td align="center">...</td>
<td align="center">0.3436</td>
</tr></tbody></table>

### Genotype file
With the increasing number of SNPs, the genotype file becomes very big and it is not possible to read it into the memory directly with a memory-limited machine. Hence, ***`KAML`*** is integrated with ***`bigmemory`*** package, which designed a specific data format named ***big.matrix*** for saving the memory usage.<br>
A total of two files should be provided: `*.geno.bin` and `*.geno.desc`, both files should use the same prefix. `*.geno.bin` is the numeric m(number of markers) by n(number of individuals) genotype file in ***big.matrix*** format, and `*.geno.desc` is the description file of `*.geno.bin`. Actually, users could manually make those files, but time-consuming and error prone, so ***`KAML`*** provides a function ***`KAML.Data()`*** for genotype format transformation. In the released version, ***`KAML`*** could accept the genotype file in four formats, including the ***Hapmap*** format, the ***VCF*** format, the ***PLINK Binary*** format, and the ***Numeric*** format. When transforming, missing genotypes are replaced by the selected methods (Left, Middle, "Right") of a given marker. After transformed, ***`KAML`*** can read it on-the-fly without a memory attenuation.<br>
***NOTE***: No matter what type of genotype format, the order of individuals in columns should be the same as phenotype file.

#### Hapmap
Hapmap is one of the commonly used data format for storing genotype. As the example shown below, the SNP information is stored in rows while the individual information is stored in columns. The first 11 columns showed attributes of the SNPs and the remaining columns show the nucleotides information that genotyped at each SNP for all individuals.

> `mouse.hmp.txt`

| rs# | alleles | chrom | pos | strand | assembly# | center | protLSID | assayLSID | panelLSID | QCcode | A048005080 | A048006063 | A048006555 | A048007096 | A048010273 | ... | A084292044 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| rs3683945 | G/A | 1 | 3197400 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs3707673 | A/G | 1 | 3407393 | + | NA | NA | NA | NA | NA | NA | GA | GA | AA | GA | AA | ... | GG |
| rs6269442 | G/A | 1 | 3492195 | + | NA | NA | NA | NA | NA | NA | AG | GG | GG | AG | GG | ... | AA |
| rs6336442 | G/A | 1 | 3580634 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs13475699 | G | 1 | 3860406 | + | NA | NA | NA | NA | NA | NA | GG | GG | GG | GG | GG | ... | GG |

The genotype in ***Hapmap*** format can be transformed to ***big.matrix*** format by the following codes:

```r
KAML.Data(hfile="mouse.hmp.txt", out="mouse")
```

If the genotype in ***Hapmap*** format is stored in multiple files for many chromosomes, it can be transformed to ***big.matrix*** format by the following codes:

```r
KAML.Data(hfile=c("mouse.chr1.hmp.txt", "mouse.chr2.hmp.txt",...), out="mouse")
```

#### VCF
***VCF*** (Variant Call Format) file has been developed with the advent of large-scale genotyping and DNA sequencing projects, such as the 1000 Genomes Project, it is one of the most widely used genotype format. An VCF file example is shown below: 

```
##fileformat=VCFv4.2
##fileDate=20171105
##source=PLINKv1.90
##contig=<ID=1,length=2>
##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional reference allele, may not be based on real reference genome">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	-9_CZTB0004	-9_CZTB0006	-9_CZTB0008	-9_CZTB0010	-9_CZTB0011	-9_CZTB0012
1	1	10000235	A	C	.	.	PR	GT	0/1	0/0	0/0	0/0	0/0	0/1
1	1	10000345	A	G	.	.	PR	GT	0/0	0/0	0/0	0/0	1/1	1/1
1	1	10004575	G	.	.	.	PR	GT	0/0	0/0	0/0	0/0	0/0	0/0
1	1	10006974	C	T	.	.	PR	GT	0/0	0/0	0/1	1/1	0/1	1/1
1	1	10006986	A	G	.	.	PR	GT	0/0	0/0	0/1	./.	1/1	1/1
```

The genotype in ***VCF*** format can be transformed to ***big.matrix*** format by the following codes:

```r
KAML.Data(vfile="mouse.vcf", out="mouse")
```

If the genotype in ***VCF*** format is stored in multiple files for many chromosomes, it can be transformed to ***big.matrix*** format by the following codes:

```r
KAML.Data(vfile=c("mouse1.vcf", "mouse2.vcf",...), out="mouse")
```

#### PLINK Binary
The ***PLINK Banary*** format is derived from Plink software. This format requires three files: \*.bed, \*.bim and
\*.fam, all with the same prefix. ***`KAML`*** only use the \*.bed and the \*.bim file.<br>
***NOTE:*** the SNP ID must be unique.

>`mouse.fam`, `mouse.bim`, `mouse.bed`

The genotype in ***PLINK Banary*** format can be transformed to ***big.matrix*** format by the following codes:

```r
KAML.Data(bfile="mouse", out="mouse")
```

#### Numeric
***`KAML`*** also accepts the ***Numeric*** format. The homozygote should be coded as 0, 2 while the heterozygote is coded as 1. The SNP information is stored in the rows and individual information is stored in the columns, it means that the dimension of the numeric matrix is m by n, **the order of individuals in columns must correspond to the phenotype file in rows**. Additionally, this format does not contain the chromosome and position of the SNPs. Therefore, two separate files should be provided, including one file contains the numeric genotype data, and the other file contains the position of each SNP.<br>
***NOTE:*** Row names and column names are not allowed, the number of row and the order of SNPs must be same in the two files.

<table>
<tbody>
<tr>
<td align="center"><em><strong><code>mouse.Numeric.txt</code></strong></em></td>
<td align="center"><em><strong><code>mouse.map</code></strong></em></td>
</tr>
<tr>
<td align="center">

<table>
<tbody>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">1</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">2</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">1</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">1</td>
<td align="center">2</td>
<td align="center">…</td>
<td align="center">0</td>
</tr>
<tr>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">0</td>
<td align="center">…</td>
<td align="center">0</td>
</tr></tbody></table>

</td>

<td align="center">
<table>
<tbody>
<tr>
<td align="center">rs3683945</td>
<td align="center">1</td>
<td align="center">3197400</td>
</tr>
<tr>
<td align="center">rs3707673</td>
<td align="center">1</td>
<td align="center">3407393</td>
</tr>
<tr>
<td align="center">rs6269442</td>
<td align="center">1</td>
<td align="center">3492195</td>
</tr>
<tr>
<td align="center">rs6336442</td>
<td align="center">1</td>
<td align="center">3580634</td>
</tr>
<tr>
<td align="center">rs13475699</td>
<td align="center">1</td>
<td align="center">3860406</td>
</tr></tbody></table>

</td>
</tr></tbody></table>

This type of file can be transformed by the following codes:
```r
KAML.Data(numfile="mouse.Numeric.txt", mapfile="mouse.map", out="mouse")
```

After transformed from one of four types of format above, two needed files will be generated, all with the same prefix which is assigned by users in the *`KAML.Data`* function.

> *The example mouse datasets:* `mouse.geno.desc, mouse.geno.bin`

---
## USAGE
### Basic

To run ***`KAML`***, you should provide two basic files: the phenotype file (values for training, NAs for predictors) and the genotype file. By default, the first column of phenotype will be analyzed, if there are more than one trait, please specify which column of trait is to be analyzed with the parameter "pheno=". For example: *`KAML(..., pheno=3)`* means the trait in third column would be analyzed. For the genotype, only the prefix need to be assigned, ***`KAML`*** could automatically attach the `*.geno.bin` and `*.geno.desc` files. <br>
***Note again:*** ***`KAML`*** has no function for adjusting the order of individuals. So please make sure the same order of individuals between phenotype and genotype.

```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse")
# pfile: phenotype file
# pheno: the column number of predicted phenotype
# gfile: transformed genotype file

# multiple threads computation
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", cpu=30)
# cpu: number pf threads
```

**Run ***`KAML`*** with the provided covariate file ***cfile*** and kinship file ***kfile***:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", cfile="CV.txt", kfile="mouse.Kin.txt")
# cfile: covariates file
# kfile: kinship file
```
**Set the sample number ***sample.num*** and validation number ***crv.num*** for cross_validation:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", sample.num=2, crv.num=5)
# sample.num: the number of replicates on cross-validation 
# crv.num: fold of cross-validation
```
**Change the top selected number of SNPs ***Top.num*** and GWAS model ***(the options are "MLM", "GLM", "RR")***:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", Top.num=15, GWAS.model="MLM")
# Top.num: max number of top LD-filtered SNPs before pseudo QTNs optimization
# GWAS.model: select the model of Genome-Wide Association Study
```
**Change the methods of variance components estimation ***vc.method*** ***(the options are "brent", "emma", "he", "ai")***:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", GWAS.model="MLM", vc.method="brent")
# vc.method: select the method of variance components estimation
```
**Change the start value of grid search procedure of Kinship optimization ***Top.perc*** & ***Logx***:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", GWAS.model="MLM", vc.method="brent",
            Top.perc=c(0.0001,0.001,0.01), Logx=c(0.01,0.05,0.1,0.5,1,5,10,15))
# Top.perc: prior value of the percentage of SNPs which would be given more weights
# Logx: prior value of the base of log function
# Note: More levels of start values will lead to much more calculation burden.
```
**Change the maximum iteration number of bisection algorithm ***bisection.loop***:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", GWAS.model="MLM", vc.method="brent",
            Top.perc=c(0.0001,0.001,0.01), Logx=c(0.01,0.05,0.1,0.5,1,5,10,15), bisection.loop=8)
# bisection.loop: the max number of iteration for bisection algorithm
# Note: if bisection.loop=0, the bisection procedure will not run.
```

### Advanced

**Only to optimize a trait-specific (weighted) kinship matrix:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", Top.num=NULL)
```
**Only to add the selected pseudo QTNs with big effects as covariates:**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", Top.perc=NULL)
```
**Switch KAML to LMM (GBLUP)**
```r
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", Top.num=NULL, Top.perc=NULL)
```
**Integrate some previously validated causal SNPs of the trait as covariates directly:**
```r
# directly predict by LM (QTN)
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", prior.QTN=c(9358, 9375), prior.model="QTN")

# predict by MLM with QTNs and a Kinship optimization procedure (QTN + weighted Kinship)
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", prior.QTN=c(9358, 9375), prior.model="QTN+K")

# predict by MLM with QTNs and a Kinship (QTN + standard Kinship)
> mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", prior.QTN=c(9358, 9375), 
            prior.model="QTN+K", Top.perc=NULL)
```
In realistic analysis, we don't know its actual genetic architecture of unknow traits, which could be obtained from a machine learning strategy of  ***`KAML`.*** Although those procedures could be speeded up by parallel computation, it's still time-consuming with limited computation resources. So it would be a better choice to run ***`KAML`*** within a smaller population to obtain the parameters, and then apply the optimized parameters to greater populations, which has been proved to be more efficient but generate similar prediction performance in our numbers of tests.

```r
mykaml <- KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", prior.QTN=c(9358, 9375), prior.model="QTN+K",
          Top.perc=0.0276, Logx=3.1094)
```
**Integrate KAML to SSBLUP (SSKAML)**

The weighted Kinship matrix can be directly applied into SSBLUP model to improve the prediction accuracy for both genotyped and non-genotyped individuals. To run SSBLUP model, please install our developed tool [HIBLUP](https://hiblup.github.io).
```r
library(hiblup)
library(KAML)
# extract the weighted kinship from KAML
G <- mykaml$K
G.ind <- read.table("mouse.Pheno.txt", head=F)[,1]
# fit SSBLUP
fit <- hiblup(pheno=pheno, pedigree=pedigree, geno=NULL, G=G, geno.ind=G.ind)
```

---
## OUTPUT
***`KAML`*** returns total 8 lists of results including: predicted phenotype, coefficients of fixed effects, predicted GEBVs, pseudo QTNs, the used model, the optimized top percentage of GWAS results, the optimized base value of Log and the optimized Kinship matrix.
```r
> str(mykaml)
List of 8
 $ y       : num [1:1940, 1] -0.109 -0.362 -0.149 -1.101 0.114 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr "Trait1"
 $ beta    : num -0.0801
 $ gebv    : num [1:1940, 1] -0.0286 -0.2823 -0.0691 -1.0205 0.1943 ...
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : NULL
  .. ..$ : chr "Trait1"
 $ qtn     : num [1:2] 9358 9375
 $ model   : chr "QTN+K"
 $ top.perc: num 0.0277
 $ logx    : num 3.11
 $ K       : num [1:1940, 1:1940] 0.8989 -0.0996 0.0538 0.0372 0.1102 ...
```

## FAQ and Hints

:sos: **Question1:** Failing to install "devtools":

***ERROR: configuration failed for package ‘git2r’***

***removing ‘/Users/acer/R/3.4/library/git2r’***

***ERROR: dependency ‘git2r’ is not available for package ‘devtools’***

***removing ‘/Users/acer/R/3.4/library/devtools’***

:yum: **Answer:** Please type the following codes in terminal.
```bash
apt-get install libssl-dev/unstable
```
---
:sos: **Question2:** When installing packages from Github with "devtools", there is a error:
 
 ***Error in curl::curl_fetch_disk(url, x$path, handle = handle): Problem with the SSL CA cert (path? access rights?)***
 
:yum: **Answer:** Please type the following codes and then try again.
```r
library(httr)
set_config(config(ssl_verifypeer = 0L))
```

**Questions, suggestions, and bug reports are welcome and appreciated.** [:arrow_right:](https://github.com/YinLiLin/R-KAML/issues)
