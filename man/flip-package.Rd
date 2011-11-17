\name{flip-package}
\alias{flip-package}
\docType{package}
\title{
The library is devoted to permutation-based inferential methods.
}
\description{
The library is devoted to Univariate and Multivariate permutation-based inferential methods.

It allows for multivariate one sample, C>=2 samples and any regression tests. 
Also the use of covariates (to be fitted in the model but) not under test is allowed.

Permutation and Rotation method to get the null distribution of the test statistic(s) are available.
}
\details{
\tabular{ll}{
Package: \tab flip\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2012-02-05\cr
License: \tab GPL <=2\cr
LazyLoad: \tab yes\cr
Depends: \tab methods, foreach, e1071, someMTP, cherry\cr
}
}
\author{
livio finos

Maintainer: livio finos <livio@stat.unipd.it>
}
\references{
Pesarin (2001) Multivariate Permutation Tests with Applications in Biostatistics. Wiley, Mew York.

I widely borrowed the structure of \code{flip} from library \code{globaltest} of (http://www.bioconductor.org/packages/release/bioc/html/globaltest.html, Authors: J. Goeman and J. Oosting).
}
\keyword{ package }

\examples{

Y=data.frame(matrix(rnorm(50),10,5))
names(Y)=LETTERS[1:5]
Y[,1:2]=Y[,1:2]
x=rep(0:1,5)
data=data.frame(x=x, Z=rnorm(10))
res = flip(Y+matrix(x*2,10,5),~x,~Z,data=data)
res

plot(res)

p2=npc(res,"fisher",subsets=list(c1=c("A","B"),c2=names(Y)))
p2
}
