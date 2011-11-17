\name{flip}
\alias{flip}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The main function for univariate tests under permutation framework
}
\description{
The general function for permutation (or rotation) test.

It allows for multivariate one sample, C>=2 samples and any regression tests. Also the use of covariates (to be fitted in the model but) not under test is allowed.
}
\usage{
flip(Y, X=NULL, Z=NULL, data=NULL, tail = 0, perms = 1000, flipReturn, rotationTest = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{Y}{The response vector of the regression model. May be
      supplied as a vector or as a \code{\link[stats:formula]{formula}}
      object. In the latter case, the right hand side of \code{Y}
      is passed on to \code{alternative} if that argument is missing, or
      otherwise to \code{null}.}

    \item{X}{The part of the design matrix corresponding to
      the alternative hypothesis. The covariates of the null model do
      not have to be supplied again here. May be given as a half
      \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}). In
      that case the intercept is always suppressed.}
    \item{Z}{The part of the design matrix corresponding to the null hypothesis. May be given as a design matrix or as a half \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}). The default for \code{Z} is \code{~1}, i.e. only an intercept. This intercept may be suppressed, if desired, with \code{Z = ~0}.}    
	\item{data}{Only used when \code{Y}, \code{X}, or \code{Z} is given in formula form. An optional data frame, list or environment containing the variables used in the formulae. If the variables in a formula are not found in \code{data}, the variables are taken from environment(formula), typically the environment from which \code{gt} is called.}
  \item{tail}{Vector of values -1, 0 or 1 indicating the tail to be used in the test for each column of \code{Y}. I the length of \code{tail} is smaller than number of coulmns of \code{Y}, the values are recycled.  }
\item{perms}{The number of permutations to use. The default is \code{perms = 1000}. Alternatively it can be a matrix (i.e. the permutation space) or a list with elements \code{number} and \code{seed} }
\item{flipReturn}{list of objects indicating what will be included in the output. e.g. \code{list(permP=TRUE,permT=TRUE,data=TRUE)}.}
\item{rotationTest}{Logical value. If \code{TRUE} resampling is performed through random linear combinations (i.e. a rotation test is performed). This option is usefull when only few permutations are availables, that is, minimun reacheable significance is hight.}
\item{\dots}{ Further parameter such as (default is reported here) \code{permT.return = TRUE, permP.return = FALSE, permSpace.return = FALSE, permY.return = FALSE} }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\note{If \code{null} is supplied as a \code{\link[stats:formula]{formula}} object, an intercept is automatically included. As a consequence \code{gt(Y, X, Z)} will usually give a different result from \code{gt(Y, X, ~Z)}. The first call is equivalent to \code{gt(Y, X, ~0+Z)}, whereas the second call is equivalent to \code{gt(Y, X, cbind(1,Z))}.

P-values from the asymptotic distribution are accurate to at least two decimal places up to a value of around \code{1e-12}. Lower p-values are numerically less reliable.

Missing values in the \code{Y} matrix are are dealed with approximate solution proposed by Pesarin (2001) or by weighted regression assuming with (almost) infinite variance for missing values.
Missing values in the \code{X} matrix are imputed conservatively (i.e. under the null hypothesis).}
\value{The function returns an object of class \code{flip.object}. 
Several operations and plots can be made from this object. See also \code{\link{flip.object-class}}.}
\author{
livio finos livioATstatDOTunipdDOTit
}

\references{
For the general framework of univariate and multivariate permutation tests see: 
Pesarin, F. (2001) Multivariate Permutation Tests with Applications in Biostatistics. Wiley, Mew York.

for Rotation tests see: 
Langsrud, O. (2005) Rotation tests, Statistics and Computing, 15, 1, 53-60
}

\seealso{
The permutation spaces on which the test is based: \code{\link{permutationSpace}} function and useful functions associated with that object.

Multiplicity correction and Overall test: \code{\link{npc}}.
}

\examples{
Y=data.frame(matrix(rnorm(50),10,5))
names(Y)=LETTERS[1:5]
Y[,1:2]=Y[,1:2]+1.5
X=rep(0:1,5)
data=data.frame(X=X, Z=rnorm(10))
#testing dependence among Y and x
res = flip(Y+matrix(X*2,10,5),~X,~Z,data=data)
res

#testing symmetry of Y around 0
Y[,1:2]=Y[,1:2] +2
res = flip(Y)
res
plot(res)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
