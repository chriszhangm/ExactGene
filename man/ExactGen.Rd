% package: ExactGene
\name{ExactGen}
\docType{package}
\title{
  Exact Inference
}
\description{
An exact meta-analysis method constructs confidence intervals, estimates, and p-values for risk difference based on exact method with identity, arcsin square root, and inverse normal CDF transformations, and Mantel-Haenszel method with and without correction.
}
\usage{
ExactGen(Xt,Xc,Nt,Nc,pcupper,pclower,Thetagrid,Pgrid,cov_prob,B,dev_region_num)
}
\arguments{
  \item{Xt}{number of events in treatment/case group(s).}
  \item{Xc}{number of events in control group(s).}
  \item{Nt}{sample size in treatment/case group(s).}
  \item{Nc}{sample size in control group(s).}
  \item{pcupper}{initial upper bounds of pc (event rate of the control group).}
  \item{pclower}{initial lower bounds of pc (event rate of the control group).}
  \item{Thetagrid}{number of grid values for risk difference. Set 200 as default.}
  \item{Pgrid}{number of grid values set for pc. Set 20 as default.}
  \item{cov_prob}{a value set for coverage probability of confidence interval.}
  \item{B}{number of iteration used to generate random uniform variables.}
  \item{dev_region_num}{a number (2 to 10) used to accelerate algorithm by dividing elements in the step of calculating individual p-value functions. Set 6 as default.}

}
\examples{
  xt = c(1:5)
  xc = c(6:10)
  nt = c(5:9)*10
  nc = c(5:9)*10
  #get pcupper and pclower by the method proposed by Tian and et.al.
  pcupper=pclower=numeric(length(xt))
  alpha=0.99 #set as the largest value
  for(i in 1:length(xt)){
  pcupper[i] = binom.confint(xc[i],nc[i],alpha,method="exact")$upper
  pclower[i] = binom.confint(xc[i],nc[i],alpha,method="exact")$lower
  }
  Thetagrid=200
  Pgrid=20
  cov_prob=0.95
  B=5e+4
  dev_region_num=6
  eee = ExactGen(Xt,Xc,Nt,Nc,pcupper,pclower,Thetagrid,Pgrid,cov_prob,B,dev_region_num)
  #row1-row5: identity,square root, and inverse normal CDF transformations, and Mantel-Haenszel method with and without correction.
  ##Need estimated Risk Difference:
  eee$Estimate
  ###p-value
  eee$pvalue
  ###Confidence Upper and Lower bound
  eee$CILower;eee$CIUpper
  ###Odds ratio
  eee#OR
}
\author{
\packageAuthor{ExactGene}
Maintainer: \packageMaintainer{ExactGene}
}
\references{
Tian, L., Cai, T., Pfeffer, M. A., Piankov, N., Cremieux, P. Y., & Wei, L. J. (2009). Exact and efficient inference procedure for meta-analysis and its application to the analysis of independent 2*2 tables with all available data but without artificial continuity correction. Biostatistics, 10(2), 275-281.

Liu, D., Liu, R. Y., & Xie, M. G. (2014). Exact meta-analysis approach for discrete data and its application to 2*2 tables with rare events. Journal of the American Statistical Association, 109(508), 1450-1465.
}


