\name{binom.confint}
\alias{binom.confint}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Tian's method to calculate restricted search
}

\usage{
binom.confint(x, n, alpha, method = "exact")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
number of events in control group.
}
  \item{n}{
sample size in control group.
}
  \item{alpha}{
Set 0.99 (conservative enough) as default.
}
  \item{method}{
Set "Exact" as default.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}

\references{
Tian, L., Cai, T., Pfeffer, M. A., Piankov, N., Cremieux, P. Y., & Wei, L. J. (2009). Exact and efficient inference procedure for meta-analysis and its application to the analysis of independent 2*2 tables with all available data but without artificial continuity correction. Biostatistics, 10(2), 275-281.
}
\author{
%%  ~~who you are~~
}

\examples{
binom.confint(10, 600, .99, method = "exact")
}
