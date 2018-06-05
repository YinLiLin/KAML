# KAML [![](https://img.shields.io/badge/Issues-1%2B-brightgreen.svg)](https://github.com/YinLiLin/R-KAML/issues) [![](https://img.shields.io/badge/Release-v1.0.1-ff69b4.svg)](https://github.com/YinLiLin/R-KAML/commits/master)

## *[K](https://github.com/YinLiLin/R-KAML)inship [A](https://github.com/YinLiLin/R-KAML)djusted [M](https://github.com/YinLiLin/R-KAML)ultiple [L](https://github.com/YinLiLin/R-KAML)ocus Best Linear Unbiased Prediction*

## Contents
* [OVERVIEW](#overview)
* [CITATION](#citation)
* [GETTING STARTED](#getting-started)
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
***`KAML`*** is originally designed to predict phenotypic value using genome- or chromosome-wide SNPs for sample traits which are controled by limited major markers or complex traits that are influenced by many minor-polygene. In brief, ***`KAML`*** incorporates pseudo QTNs as fixed effects and a trait-specific K matrix as random effect in a mixed linear model. Both pseudo QTNs and trait-specific K matrix are optimized using a parallel-accelerated machine learning strategy.


<p align="center">
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

***`KAML`*** is developed by [***Lilin Yin***](https://github.com/YinLiLin), [***Xiaolei Liu***](https://github.com/XiaoleiLiuBio)**\*** at the [***Huazhong(Centra of China) agriculture University***](http://www.hzau.edu.cn/en/HOME.htm).

If you have any bug reports or questions please send an email to [**Lilin Yin**](https://github.com/YinLiLin) at **ylilin@163.com**

---
## CITATION
(waiting for updating)

---
## GETTING STARTED
***`KAML`*** is compatible with both [R](https://www.r-project.org/) and [Microsoft R Open](https://mran.microsoft.com/open/), We strongly recommend [Microsoft R Open](https://mran.microsoft.com/open/) instead of [R](https://www.r-project.org/) for running ***`KAML`*** . **MRO** is the enhanced distribution of **R**, it includes multi-threaded math libraries. These libraries make it possible for so many common R operations, ***such as matrix multiply/inverse, matrix decomposition, and some higher-level matrix operations***, to compute in parallel and use all of the processing power available to [reduce computation times](https://mran.microsoft.com/documents/rro/multithread/#mt-bench).

### Installation
***`KAML`*** is not available on CRAN, but can be installed using the R package **"devtools"**. There are two packages should be installed beforehand, **"snpStats"** and **"rfunctions"**. ***`KAML`*** can be installed with the following R code:
```r
#if "devtools" isn't installed, please "install.packages(devtools)" first.
devtools::install_version('RcppEigen', version = "0.3.2.9.0")
devtools::install_github("hclimente/snpStats")
devtools::install_github("jaredhuling/rfunctions")
devtools::install_github("YinLiLin/R-KAML/KAML")
```
After installed successfully, ***`KAML`*** can be loaded by typing
```r
library(KAML)
```
Typing `?KAML` could get the details of all parameters.

### Test Datasets

The example data is available for Linux by:
```bash
wget https://raw.githubusercontent.com/YinLiLin/R-KAML/master/example/example.zip
unzip example.zip
cd example
```

**Or** click [here](https://raw.githubusercontent.com/YinLiLin/R-KAML/master/example/example.zip) in your browser to download for windows, after downloaded, unzip the file and change the workplace by ``` setwd("")``` in R.

---

## INPUT
### Phenotype file
The file must contain a header row. Missing values should be denoted by NA, which will be treated as candidates. Notice that only numeric values are allowed and characters will not be recognized. However, if a phenotype takes only values 0, 1(or only two levels), then ***`KAML`*** would consider it to be a case-control trait, and the predicted value could be directly interpreted as the probability of being a case. <br>
When a phenotype file contains more than one trait, users should specify which to analyse using the option "pheno=N", for example: *`KAML(..., pheno=1)`* means the trait in first column would be predicted. 

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
Generally, there are no covariates when predicting candidates in most cases, especially genomic selection of animal breeding, because the predicted values are not original phenotypes but the (genomic) estimated breeding value(GEBV/EBV), which have been corrected by covariates. So **Covariate file** is **optional**, in order to fit the model for original phenotype prediction, users can provide the covariates in file. If provided, NAs are not allowed in the file, **the order of all individuals must be corresponding to phenotype file**. Actually, there are two parameters for covariates: **`dcovfile`** and **`qcovfile`**;

**`dcovfile`*****(optional)***: Input discrete covariates from a plain text file, e.g. *dcov.txt*. Each discrete covariate is recognized as a categorical factor with several levels. The levels of each factor can be represented by a single character, word or numerical number. NOTE: the design matrix of the mean in the model (which is a vector of all ones) is always a linear combination of the design matrix of a discrete covariate so that not all the effects of the levels (or classes, e.g. male and female) of a discrete covariate are estimable. ***`KAML`*** will always constrain the effect of the first level to be zero and the effect of any other level represents its difference in effect compared to the first level.

**`qcovfile`*****(optional)***: Input quantitative covariates from a plain text file, e.g. *qcov.txt*. Each quantitative covariate is recognized as a continuous variable.

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
***`KAML`*** requires a n×n relatedness matrix. By default, it is automatically calculated using choosed one type of three algorithms(
***"scale","center","vanraden"***
), but can be supplied in file by users. If in this case, **the order of individuals for each row and each column in the file must correspond to phenotype file, no column and row names**.

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

With the increasing number of SNPs through whole genome, the storage of genotype is a big problem. Obviously, it's not a good choice to read it into memory with a memory-limited PC directly. Here, ***`KAML`*** is integrated with a memory-efficient tool named ***bigmemory*** and could obtain the genotype information from disk instead, which can save much of memory to do more analysis.<br>
By default, total two files should be provided: `*.geno.bin`, `*.geno.desc` and **all with the same prefix**. `*.geno.bin` is the numeric m(number of markers) by n(number of individuals) genotype file in *big.matrix* format, and `*.geno.desc` is the description file of `*.geno.bin`. Actually, users could manually make those files, but time-consuming and error prone, so ***`KAML`*** provides a function ***`KAML.Data()`*** for genotype format transformation. Currently, genotype file could be in four formats, either in the ***Hapmap*** format, ***VCF*** format, ***PLINK binary ped*** format, or the m by n ***Numeric*** format. When transforming, missing genotypes are replaced by the mean genotype value of a given marker. After transformed, ***`KAML`*** can read it on-the-fly without a memory attenuation.<br>
***NOTE***: **No matter what type of format of genotype, the order of individuals in columns of the file must correspond to phenotype file in rows.**

#### Hapmap
Hapmap is the most popular used format for storing genotype data. As the example below, the SNP information is stored in the rows and individuals information is stored in the columns. The first 11 columns display attributes of the SNPs and the remaining columns show the nucleotides genotyped at each SNP for all individuals.

> `mouse.hmp.txt`

| rs# | alleles | chrom | pos | strand | assembly# | center | protLSID | assayLSID | panelLSID | QCcode | A048005080 | A048006063 | A048006555 | A048007096 | A048010273 | ... | A084292044 |
| :---: | :---: |  :---: |  :---: |  :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| rs3683945 | G/A | 1 | 3197400 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs3707673 | A/G | 1 | 3407393 | + | NA | NA | NA | NA | NA | NA | GA | GA | AA | GA | AA | ... | GG |
| rs6269442 | G/A | 1 | 3492195 | + | NA | NA | NA | NA | NA | NA | AG | GG | GG | AG | GG | ... | AA |
| rs6336442 | G/A | 1 | 3580634 | + | NA | NA | NA | NA | NA | NA | AG | AG | GG | AG | GG | ... | AA |
| rs13475699 | G | 1 | 3860406 | + | NA | NA | NA | NA | NA | NA | GG | GG | GG | GG | GG | ... | GG |

This type of file can be transformed by the following codes:

```r
KAML.Data(hfile="mouse.hmp.txt", out="mouse")
```

Normally, all chromosomes are stored in one file, but can be stored in separated files for different chromosomes. If in this case, it can be transformed by following codes:

```r
KAML.Data(hfile=c("mouse.chr1.hmp.txt", "mouse.chr2.hmp.txt",...), out="mouse")
```

#### VCF

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


#### PLINK Binary
The **PLINK Banary** format is derived from Plink software. This format requires three files: \*.bed, \*.bim and
\*.fam, all with the same prefix. ***`KAML`*** only use \*.bed and \*.bim file. ***NOTE*** that the id of SNPs must be unique.

>`mouse.fam`, `mouse.bim`, `mouse.bed`

This type of file can be transformed by the following codes:

```r
KAML.Data(bfile="mouse", out="mouse")
```

#### Numeric
***`KAML`*** also accepts the numeric format. All nucleotides have been coded as 0, 1, 2. The SNP information is stored in the rows and individuals information is stored in the columns, it means that the dimension of the numeric matrix is m by n, **the order of individuals in columns must correspond to phenotype file in rows**. Additionally, this format does not contain the chromosome and position of the SNPs. Therefore, two separate files must be provided. One file contains the numeric genotypic data, and the other contains the position of each SNP. ***NOTE:*** **Row names and column names are not allowed, the number of row and the order of SNPs must be same in those two files.**

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

After transformed from one of three type of format above, three file will be generated, all with the same prefix which is assigned by users in the *`KAML.Data`* function.

> *The example mouse datasets:*
> `mouse.map, mouse.geno.desc, mouse.geno.bin`

---
## USAGE
### Basic

```r
KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse")
```
```r
KAML(pfile="mouse.Pheno.txt", pheno=1, gfile="mouse", cfile="CV.txt", kfile="mouse.Kin.txt")
```
### Advanced

<p align="center">
<a href="https://raw.githubusercontent.com/YinLiLin/R-KAML/master/figures/Trait1.jpg">
<img src="/figures/Trait1.jpg" height="300px" width="840px"/>
</a>
</p>

---
## OUTPUT

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
 
:yum: **Answer:** Please type the following codes and then try agian.
```r
library(httr)
set_config(config(ssl_verifypeer = 0L))
```

**Questions, suggestions, and bug reports are welcome and appreciated.** [:arrow_right:](https://github.com/YinLiLin/R-KAML/issues)
