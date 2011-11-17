\name{flip-package}
\alias{flip-package}
\docType{package}
\title{
The library is devoted to permutation-based inferential methods.
}
\description{
The library is devoted to permutation-based inferential methods. Univariate and multivariate methods are developed.
}
\details{
\tabular{ll}{
Package: \tab flip\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-06-20\cr
License: \tab GPL <=2\cr
LazyLoad: \tab yes\cr
Depends: \tab methods\cr
}
}
\author{
livio finos

Maintainer: livio finos <livio@stat.unipd.it>
}
\references{
Pesarin (2001) Multivariate Permutation Tests with Applications in Biostatistics. Wiley, Mew York.
}
\keyword{ package }
\seealso{
}
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
