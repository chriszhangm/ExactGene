# ExactGene
Exact inference for meta-analysis of rare events based on the fixed-effect model
## Installation of the package

To install our package, you may execute the following codes:

```{r, eval = FALSE}
install.packages("devtools")
devtools::install_github("chriszhangm/ExactGene")
library(ExactGene)
```
For Mac Users who cannot compile the code, please refer [this answer](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).
## A example of using _ExactGene_
We show a real application to apply the function `ExactGen`, which constructs confidence intervals, estimates, and p-values for risk difference based on exact method with identity, arcsin square root, and inverse normal CDF transformations, and Mantel-Haenszel method with and without correction. The dataset we use is a meta-analysis conducted by Paolo Zanoni and et.al (2016) about the gene scavenger receptor class B member 1 (SCARB1) encodes scavenger receptor BI (SR-BI).
```{r, eval = FALSE}
------------------------------------------------------------------------------------------------------------------------------------------------------
#input data
xt = c(6,1,1,0,0,0,1,3,4,8,1,4,4,0,1,0)
xc = c(10,0,0,4,1,1,1,0,0,29,1,3,2,0,0,0)
nt = c(4587,2833,1568,2483,4464,1024,728,683,5255,2860,2020,8079,9810,2153,640,659)
nc = c(16546,5912,2772,8085,2886,2267,808,156,2921,14929,6087,10367,10970,2118,638,687)
------------------------------------------------------------------------------------------------------------------------------------------------------
#find pcupper and pclower
pcupper=pclower=numeric(length(xt))
for (i in 1:length(xt)) {
  pcupper[i]=binom.confint(xc[i],nc[i],.99,"exact")$upper
  pclower[i]=binom.confint(xc[i],nc[i],.99,"exact")$lower
}
------------------------------------------------------------------------------------------------------------------------------------------------------
#define configurations
Thetagrid=500
Pgrid=20
cov_prob=0.95
B=5e+4
dev_region_num=6
set.seed(1234)
#run exact method
eee = ExactGen(Xt=xt,Xc=xc,Nt=nt,Nc=nc,pcupper,pclower,Thetagrid,Pgrid,cov_prob,B,dev_region_num)
output = rbind(t(eee$Estimate), t(eee$OR), t(eee$CILower),t(eee$CIUpper),t(eee$pvalue))
rownames(output)=c("RD Estimate", "Converted OR", "CI Lower Bound", "CI Upper Bound", "P-value")
colnames(output)=c(" Exact-Identity","  Exact-Asin_sqrt", "  Exact-IN-CDF", "  MH w/o correction", "  MH w/ correction")
options(digits = 4)
------------------------------------------------------------------------------------------------------------------------------------------------------
#show the output
output
------------------------------------------------------------------------------------------------------------------------------------------------------
                Exact-Identity   Exact-Asin_sqrt   Exact-IN-CDF   MH w/o correction   MH w/ correction
RD Estimate          2.868e-04         2.741e-04      2.358e-04           3.207e-04          3.281e-04
Converted OR         1.704e+00         1.673e+00      1.579e+00           1.787e+00          1.805e+00
CI Lower Bound      -6.374e-06        -4.462e-05     -8.286e-05           3.458e-05          2.002e-05
CI Upper Bound       6.183e-04         6.183e-04      6.438e-04           6.068e-04          6.361e-04
P-value              7.156e-02         1.653e-01      4.103e-01           2.803e-02          3.686e-02
```

## reference
Tian, L., Cai, T., Pfeffer, M. A., Piankov, N., Cremieux, P. Y., & Wei, L. J. (2009). Exact and efficient inference procedure for meta-analysis and its application to the analysis of independent 2*2 tables with all available data but without artificial continuity correction. Biostatistics, 10(2), 275-281.
Liu, D., Liu, R. Y., & Xie, M. G. (2014). Exact meta-analysis approach for discrete data and its application to 2*2 tables with rare events. Journal of the American Statistical Association, 109(508), 1450-1465.
Paolo Zanoni, Sumeet A Khetarpal, Daniel B Larach, William F Hancock-Cerutti, John S Millar, Marina Cuchel, Stephanie DerOhannessian, Anatol Kontush, Praveen Surendran, Danish Saleheen, et al. Rare variant in scavenger receptor bi raises hdl cholesterol and increases risk of coronary heart disease. Science, 351(6278):1166–1171, 2016.
