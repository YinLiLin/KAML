\name{KAML-package}
\alias{KAML-package}
\alias{KAML}
\docType{package}
\title{
KAML: Kinship Adjusted Multiple Loci Best Linear Unbiased Prediction
}
\description{
KAML is originally designed to predict phenotypic value using genome- or chromosome-wide SNPs for sample traits which are controled by limited major markers or complex traits that are influenced by many minor-polygene. In brief, KAML incorporates pseudo QTNs as fixed effects and a trait-specific K matrix as random effect in a mixed linear model. Both pseudo QTNs and trait-specific K matrix are optimized using a parallel-accelerated machine learning strategy.}
\usage{
KAML(pfile="", gfile="", kfile=NULL, dcovfile=NULL, qcovfile=NULL,
  pheno=1, GWAS.model=c("MLM","GLM", "RR"), GWAS.npc=NULL,
  prior.QTN=NULL, prior.model=c("QTN+K", "QTN", "K"),
  vc.method=c("brent", "he", "emma", "ai"),
  K.method=c("center", "scale", "vanraden"),
  Top.perc=c(1e-4, 1e-3, 1e-2, 1e-1), Top.num=10,
  Logx=c(1.01, 1.11, exp(1), 10), qtn.model=c("MR", "SR", "BF"),
  BF.threshold=NULL, binary=FALSE, bin.size=1000000, max.nQTN=TRUE,
  sample.num=2, SNP.filter=NULL, crv.num=5, cor.threshold=0.3,
  count.threshold=0.9, bisection.loop=10, ref.gwas=TRUE,
  theSeed=666666, file.output=TRUE, cpu=10
)
}
\arguments{
\item{pfile}{phenotype file, one column for a trait, the name of each column must be provided(NA is allowed)}
\item{gfile}{genotype files, including "gfile.geno.desc", "gfile.geno.bin" and "gfile.map"}
\item{kfile}{n*n, optional, provided KINSHIP file for all individuals}
\item{dcovfile}{n*x, optional, the provided discrete covariates file}
\item{qcovfile}{n*x, optional, the provided quantitative covariates file}
\item{pheno}{specify phenotype column in the phenotype file(default 1)}
\item{GWAS.model}{which model will be used for GWAS(only "GLM" and "MLM" can be selected presently)}
\item{GWAS.npc}{the number of PC that will be added as covariance to control population structure}
\item{prior.QTN}{the prior QTNs which will be added as covariants, if provided prior QTNs,  KAML will not optimize QTNs and model during cross-validation}
\item{prior.model}{the prior Model for the prior.QTN that added as covariants}
\item{vc.method}{method for variance components estimation("emma" or "gemma")}
\item{K.method}{which algorithm will be applied("center", "scale", "vanraden")}
\item{Top.perc}{a vector, a subset of top SNPs for each iteration are amplified when calculating KINSHIP}
\item{Top.num}{a number, a subset of top SNPs for each iteration are used as covariants}
\item{Logx}{a vector, the base for LOG}
\item{qtn.model}{the strategy of selecting pseudo QTNs. c("MR", "SR", "BF")}
\item{BF.threshold}{the threshold of BF method}
\item{binary}{whether the phenotype is case-control}
\item{bin.size}{the size of each bin}
\item{max.nQTN}{whether limits the max number of Top.num}
\item{sample.num}{the sample number of cross-validation}
\item{SNP.filter}{the SNPs whose P-value below this threshold will be deleted}
\item{crv.num}{the cross number of cross-validation}
\item{cor.threshold}{if the top SNP which is used as covariant is in high correlation with others, it will be deleted}
\item{count.threshold}{if the count of selected SNP for all iteration >= sample.num*crv.num*count.threshold, than it will be treated as covariance in final predicting model}
\item{bisection.loop}{the max loop(iteration) number of bisection algorithm}
\item{ref.gwas}{whether to do GWAS for reference population(if not, KAML will merge all GWAS results of cross-validation by mean)}
\item{theSeed}{the random seed}
\item{file.output}{whether to write the predicted values in file}
\item{cpu}{the number of CPU for calculation}
}
\details{
\tabular{ll}{
Package: \tab KAML\cr
Type: \tab Package\cr
Version: \tab 0.99.1\cr
Date: \tab 2019-7-11\cr
License: \tab GPL(>=2)\cr
}
}
\author{
Lilin Yin, Haohao Zhang and Xiaolei Liu\cr
Maintainer: \cr
Lilin Yin <\bold{ylilin@163.com}> \cr
Xiaolei Liu <\bold{xiaoleiliu@mail.hzau.edu.cn}> \cr
}
\examples{
Please see at: https://github.com/YinLiLin/R-KAML
}
