flip <- function(Y, X=NULL, Z=NULL, data=NULL, tail = 0, perms = 1000, rotationTest = FALSE, ...) {
  if(!exists("flipReturn") || is.null(flipReturn)) 
  flipReturn=list(permT=TRUE,permP=FALSE,permSpace=FALSE,data=FALSE)
  
  # store the call
  call <- match.call()
  data=.getXY(Y,X,Z,data)

  
  ### residualize by Z and orthogonalize X and Y
 #forse non serve:  data$Y = as.matrix(data$Y)
  if ((dim(data$Z)[2]>=1) && ((any(var(data$Z)>0)) || (rotationTest)) ) {
    ##if possible fill missing data (conservatively)
	data$X <- .fillX(data$X,data$Z)
    ##orthogocanlize the residuals on Z
	data <- .orthoZ(data$Y,data$X,data$Z)
	}

	
  if(!missing(rotationTest)&&(rotationTest==TRUE)){
	res = .rotation.nptest(data$Y, data$X ,B=1000,...)
  } else {
	if(length(unique(data$X))>1){ # if X is not a constant go on
		data$X=data$X[,apply(data$X,2,var)>0]
		res= .dependence.nptest(Y=data$Y, X=data$X, perms=perms,  tail = tail,...)
	} else { #otherwise perform a symmetry test
		res= .symmetry.nptest(Y=data$Y, perms=perms,  tail = tail,...)
	}
 }
  
	p=t2p(res$permT,obs.only=TRUE,tail=tail)
	
	
	#build the flip-object
	out=.getOut(permSpace=res$permSpace,permP=p,permT=res$permT, data=data,tail=tail, call=call, flipReturn=flipReturn,tstat=res$tStats)
  return(out)
}

