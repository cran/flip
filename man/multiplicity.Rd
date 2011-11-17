\name{npc}
\alias{npc}
\alias{flip.adjust}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for multiplicity corrections
}
\description{
\code{npc} provides overall tests (i.e. weak FWER control), while \code{flip.adjust} provides adjusted p-values.
}
\usage{
npc(permTP, comb.funct = c("Fisher", "Liptak", "minP", "maxT", "sumT", "sumT2"),subsets=NULL,weights=NULL,...)
flip.adjust(permTP, method = "maxT", maxalpha = 1, weights = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{permTP}{
  A permutation space (B times m matrix) or an \code{flip.object} as produced by \code{\link{flip}}.
}
  \item{comb.funct}{ A cobining function among: "Fisher", "Liptak", "minP" (i.e. Tippet), "maxT", "sumT" (i.e. direct) , "sumT2" (sum of T^2)
  see also Pesarin (2001) for further details.
}
  \item{method}{ A method among:
"maxT", "maxTstd" (i.e. max T on \code{scale(permTP)}) ,"minP" or any \code{p.adjust.methods} of function \code{p.adjust}
}
  \item{maxalpha}{ Adjusted p-vlaues greather than \code{maxalpha} are forced to 1. It saves computational time when there are many hypotheses under test.
}
  \item{weights}{
Optional argument that can be used to give certain variables greater weight in the combined test. 
Can be a vector or a list of vectors. In the latter case, a separate test will be performed for each weight vector. 
If both \code{subsets} and \code{weights} are specified as a list, they must have the same length. In that case, \code{weights} vectors may have either the same length as the number of covariates in \code{alternative}, or the same length as the corresponding subset vector. Weights can be negative; the sign has no effect unless \code{directional} is \code{TRUE}.
It works for \code{npc} and \code{flip.adjust} with \code{method}= "maxT", "maxTstd" or "minP"
}
  \item{subsets}{Optional argument that can be used to test one or more subsets of variables. 
  Can be a vector of column names or indices of a \code{\link{flip.object-class}}, or a list of such vectors. 
  In the latter case, a separate test will be performed for each subset.}
  \item{...}{further arguments }
}


\value{
 a \code{\link{flip.object-class}}
}
\author{
livio finos
}

\examples{

Y=data.frame(matrix(rnorm(50),10,5))
names(Y)=LETTERS[1:5]
Y[,1:2]=Y[,1:2]+1.5
res = flip(Y,perms=10000)

########npc
p2=npc(res,"tippet")
summary(p2)
p2=npc(res,"fisher",subsets=list(c1=c("A","B"),c2=names(Y)))
summary(p2)
p2=npc(res,"fisher",subsets=list(c1=c("A","B"),c2=names(Y)),weights=1:5)
summary(p2)



flip.adjust(res,"maxT")

flip.adjust(res,"BH")
##same as 
#p.adjust(res,"BH")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
