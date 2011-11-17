\name{npc}
\alias{npc}
\alias{flip.adjust}
\alias{flip.npc.methods}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions for multiplicity corrections
}
\description{
\code{npc} provides overall tests (i.e. weak FWER control), while \code{flip.adjust} provides adjusted p-values (i.e. strong FWER control).
}
\usage{
npc(permTP, comb.funct = c(flip.npc.methods, p.adjust.methods), subsets=NULL,weights=NULL, stdSpace=FALSE, ...)
flip.adjust(permTP, method = flip.npc.methods, maxalpha = 1, weights = NULL, stdSpace=FALSE, ...)

flip.npc.methods 
#    c("Fisher", "Liptak", "minP", "maxT", "maxTstd", "sumT", "sumT2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{permTP}{
  A permutation space (B times m matrix) or an \code{flip.object} as produced by \code{\link{flip}}. Alternativelly it can be a \code{\link{flip.object-class}} resulting, for example from a call of function \code{\link{flip}}.
}
  \item{comb.funct}{ A cobining function \code{flip.npc.methods}: "Fisher", "Liptak", "minP" (i.e. Tippet), "maxT", "sumT" (i.e. direct) , "sumT2" (sum of T^2).
  \code{"Fisher"} combining function is the default.
  See also the section \code{Details}.
}
  \item{method}{ A method among \code{\link{flip.npc.methods}} or \code{\link{p.adjust.methods}}. By default \code{"maxT"} is used. 
  See also the section \code{Details}.
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
 
 \item{stdSpace}{Ask if the permutation distribution of the test statistic should be standardized or not. The default is \code{FALSE}. The option is applied only 
 if \code{comb.funct} or \code{method} is equal to \code{"maxT"} or \code{"sumT"}, it becomes usefull when test statistics are of different nature (e.g. chisquare and t-test).}
  
  \item{...}{further arguments }
  
}

\details{
\code{npc} combine the p-values using the combinig functions (and the method) described in Pesarin (2001). It make use of the join space of the permutations. This is usually derived from a call of \code{flip} function.

\code{flip.adjust} adjusts the p-value for multiplicity (FamilyWise Error Rate).
When \code{method} is equal to \code{"maxT"}, \code{"maxTstd"} (i.e. max T on \code{scale(permTP)}) or "minP" (i.e. Tippett) it performs the step-down method of Westfall and Young (1993).
 For any other element of \code{flip.npc.methods} (i.e. "Fisher", "Liptak", "sumT" (i.e. direct) or "sumT2" (sum of T^2)) a call to \code{noc} together with a closed testing procedure is used (it make use of  \code{\link[cherry:closed]{cherry:closed}}).
 When \code{method} is any among \code{p.adjust.methods} the function \code{stats:p.adjust} or -if weights are provided- \code{someMTP:p.adjust.w} is used.
}


\value{
 The function returns an object of class \code{\link{flip.object-class}}.
}

\author{
livio finos
}

\references{
Pesarin (2001) Multivariate Permutation Tests with Applications in Biostatistics. Wiley, Mew York.

P. H. Westfall and S. S. Young (1993). Resampling-based multiple testing: Examples and methods for p-value adjustment. John Wiley & Sons.
}

\examples{

Y=data.frame(matrix(rnorm(50),10,5))
names(Y)=LETTERS[1:5]
Y[,1:2]=Y[,1:2]+1.5
res = flip(Y,perms=10000)

########npc
p2=npc(res,"minP")
summary(p2)
p2=npc(res,"Fisher",subsets=list(c1=c("A","B"),c2=names(Y)))
summary(p2)
p2=npc(res,"Fisher",subsets=list(c1=c("A","B"),c2=names(Y)),weights=1:5)
summary(p2)



flip.adjust(res,"maxT")

flip.adjust(res,"BH")
##same as 
#p.adjust(res,"BH")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
