flip <- function(Y, X=NULL, Z=NULL, data=NULL, tail = NULL, perms = 1000,...) {
  
  # store the call
  call <- match.call()
  data=.getXY(Y,X,Z,data)

  
  ### residualize by Z and orthogonalize X and Y
 #forse non serve:  data$Y = as.matrix(data$Y)
  if (dim(data$Z)[2]>1) if(  any(var(data$Z)>0) ){
    ##if possible fill missing data (conservatively)
	data$X <- .fillX(data$X,data$Z)
    ##orthogocanlize the residuals on Z
	data <- .orthoZ(data$Y,data$X,data$Z)
	}
	
  if(length(unique(data$X))>1){ # if X is not a constant go on
	out= .dependence.nptest(Y=data$Y, X=data$X, perms=perms,  tail = tail,...)
  } else { #otherwise perform a symmetry test
	out= .symmetry.nptest(Y=data$Y, perms=perms,  tail = tail,...)
  }
  out@call=call
  return(out)
}

