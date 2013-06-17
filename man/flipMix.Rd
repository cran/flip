\name{flipMix}
\alias{flipMix}
\alias{obs2coeffWithin}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The main function for testing mixed models under a permutation (and rotation) framework
}
\description{
It allows to test fixed effect in mixed models. You can test within-unit effects, between-unit and interactions of the two.
The response can be uni- or multi-variate.
See also examples below.
}
\usage{
flipMix(modelWithin, X = NULL, Z = NULL, units, perms = 1000, data = NULL, tail = NULL,
                    statTest = NULL, flipReturn, testType = "permutation", 
                    Su = NULL, equal.se = FALSE, se = NA, replaceNA.coeffWithin = NA, 
                    replaceNA.coeffWithin.se = 0, ...)

obs2coeffWithin(modelWithin, units, X = NULL, Z = NULL, data = NULL, equal.se = FALSE, 
          se = NA, replaceNA.coeffWithin = NA, replaceNA.coeffWithin.se = Inf, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{modelWithin}{ When it is a  \code{\link[stats:formula]{formula}} object, a (possibly multivariate) multiple linear model is fitted. 
     Responses are on the left, while the right part contains ONLY within-unit variables. 
     In this case \code{data} must be supplied.
     Alternatively, it can be a \code{glm}, a \code{lm} or \code{vgam} (\code{library(VGAM)}) (i.e. \code{vglm}) object.
     The \code{modelWithin} have to be performed using only variables within-unit, without using \code{units} indicator (in this case the argument \code{data} is not used).
     It can be  also  a list of models.
It can be null if data is provided in the right format (see below).
}
\item{X}{The part of the design matrix corresponding to the between-unit effect that are not null under 
      the alternative hypothesis. If it is a matrix or a data.frame it must have a number of rows equal to the number of units or equal to the total number of observations (in the latest case all elements of the same units must have the same values since they are between-unit effects).
The non-null between-unit covariates of null model are defined in \code{Z} (see argument below) and do
      not have to be supplied again here.
See also the function \code{\link{flip}}.
}
    \item{Z}{The part of the design matrix corresponding to the non-null between-unit covariates of the model under the null hypothesis. May be given as a design matrix or as a half \code{\link[stats:formula]{formula}} object (e.g. \code{~a+b}). See also the function \code{\link{flip}}.
 If it is a matrix or a data.frame it must have a number of rows equal to the number of units or equal to the total number of observations (in the latest case all elements of the same units must have the same values since they are between-unit effects).
}
  \item{units}{
Vector of units IDs. May be given as a vector or as a half \code{\link[stats:formula]{formula}} object (e.g. \code{~subj}).
}
\item{perms}{The number of permutations to use. The default is \code{perms = 1000}. Alternatively it can be a matrix (i.e. the permutation space) or a list with elements \code{number} and \code{seed}. See also the function \code{\link{flip}}.}
\item{data}{Same as in the function \code{\link{flip}}. If can also be the results of \code{obs2coeffWithin}.
}
  \item{tail}{
Same as in the function \code{\link{flip}}.
}
  \item{statTest}{ Choose among \code{"t"}, \code{"F"} and \code{"Trace"}.
}
  \item{flipReturn}{
Same as in the function \code{\link{flip}}.
}
  \item{testType}{
See also the function \code{\link{flip}}.
}
  \item{Su}{Usually \code{NULL}. It is the covariance matrix of the random effects. If not supplied, it is estimated by iterative least square algorithm.
}
  \item{equal.se}{
Logical. If \code{TRUE} it force the unit to have the same variance of errors (like it is usually assumed in the lmer methods).
}
  \item{se}{Usually \code{NULL}. It is a matrix of unit-specific standard errors. If not supplied it is estimated by the algorithm.
}
  \item{replaceNA.coeffWithin}{ deafult is \code{NA}
}
  \item{replaceNA.coeffWithin.se}{deafult is \code{Inf}
}
\item{\dots}{Further parameters. 
\code{test.coeffWithin} Vector of names or IDs of within-unit variables that have to be tested (and reported). Note that variables not in the list are used in the model (i.e. the they play the role of nuisance parameters).
See also the function \code{\link{flip}} for other parameters.
}
}
\value{An object of class \code{flip.object}. 
Several operations and plots can be made from this object. See also \code{\link{flip.object-class}}.
}

\references{
L. Finos and D. Basso (2013) Permutation Tests for Between-Unit Fixed Effects in Multivariate Generalized Linear Mixed Models. Statistics and Computing. Accepted.
}

\author{
Livio Finos and Dario Basso
}

\seealso{
 \code{\link{flip}}
}

\examples{
n=10
toyData=  data.frame(subj=rep(1:n,rep(4,n)), Within=rep(1:2,n*2),
          XBetween= rep(1:2,rep(n/2*4,2)),ZBetween= rep(rnorm(n/2),rep(8,n/2)))
toyData=cbind(Y1=rnorm(toyData$subj+toyData$ZBetween+toyData$XBetween+toyData$Within),
              Y2=rnorm(toyData$subj+toyData$ZBetween+toyData$Within),toyData)
(toyData)

#####################
###Testing Between-unit effects
(flipMix(modelWithin=as.matrix(toyData[,c("Y1","Y2")])~Within,data=toyData, 
      X=~XBetween,Z=~ZBetween,units=~subj,perms=1000,testType="permutation",statTest="t"))
#same as:
modelWithin <- lm(as.matrix(toyData[,c("Y1","Y2")])~Within,data=toyData)
(flipMix(modelWithin=modelWithin,data=toyData, X=~XBetween,Z=~ZBetween,units= ~subj, 
        perms=1000,testType="permutation",statTest="t"))


### Note that this is different from:
modelWithin <- list(Y1=lm(Y1~Within,data=toyData),Y2=lm(Y2~Within,data=toyData))
(flipMix(modelWithin=modelWithin,data=toyData, X=~XBetween,Z=~ZBetween,units= ~subj,
        perms=1000,testType="permutation",statTest="t"))

################################
###Testing Within-unit effects
(flipMix(modelWithin=as.matrix(toyData[,c("Y1","Y2")])~Within,data=toyData, units= ~subj, 
        perms=1000,testType="permutation",statTest="t"))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}