% package: ExactGene
\name{individual_p_cpp}
\docType{package}
\title{
  Generate individual p value
}
\description{
This function helps generate right-tailed and left-tailed individual p-value functions for each study.
}
\usage{
individual_p_cpp(xt, xc, nt, nc, theta_grd, Pgrid, dev_region_num)
}
\arguments{
  \item{xt}{number of events in treatment/case group.}
  \item{xc}{number of events in control group.}
  \item{nt}{sample size in treatment/case group.}
  \item{nc}{sample size in control group.}
  \item{theta_grd}{a grid values set for theta (risk difference between the treatment and control group.)}
  \item{Pgrid}{number of grid values set for pc (event rate of the control group). Set 20 as default.}
  \item{dev_region_num}{a number (2 to 10) used to accelerate algorithm by dividing elements in the step of calculating individual p-value functions. Set 6 as default.}

}
\examples{
  xt = 10
  xc = 5
  nt = 500
  nc = 600
  theta_grd = seq(-0.1,0.1,length.out=1000)
  Pgrid = 20
  dev_region_num=6
  ind = individual_p_cpp(xt, xc, nt, nc, theta_grd, Pgrid, dev_region_num)
  #get the right tailed p value functions:
  upper = ind$p_up
  #get the left tailed p value functions:
  upper = ind$p_low
}
\author{
\packageAuthor{ExactGene}
Maintainer: \packageMaintainer{ExactGene}
}


