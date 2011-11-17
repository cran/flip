flip <- function(Y, X=NULL, Z=NULL, data=NULL, tail = 0, perms = 1000, flipReturn, rotationTest=FALSE, ...) {

  if(missing(flipReturn)||is.null(flipReturn)) 
  flipReturn=list(permT=TRUE,permP=FALSE,permSpace=FALSE,data=FALSE)
  
  # store the call
  call <- match.call()
  # get matrices from inputs
  data=.getXY(Y,X,Z,data)

  
  ##########GIA' inserito in getXY
  # # if Z have at least 1 column fill the missing values of X
  # if (dim(data$Z)[2]>=1  )  {
	# if possible fill missing data (conservatively)
	# data$X <- .fillX(data$X,data$Z)
	# }
  
  
  ################ solution for missing values if any NA in Y AND not a symmetry test (for symmetry test uses standard solution)
  if ( any(is.na(data$Y)) && (length(unique(data$X))>1)   ){
	#applica soluzione missing values pesarin
	#################
	##################assicurarsi che X sia di variabili dummy!
	###################
	if( (ncol(data$Z)==0) && (!rotationTest) )  { 
		##only 1 grouping variable is allowed, compact if X has more than 1 column
		res <- .dependenceNA.nptest(Y=data$Y, X=data$X, perms=perms,  tail = tail, ...)
	} else {
	#applica regr pesata
		data$W=10^(-5*(apply(data$Y,2,is.na)))
		data$X=data$X[,apply(data$X,2,var)>0,drop=FALSE]
		data$Y[is.na(data$Y)]=0
		res=.dependenceW.nptest(Y=data$Y, X=data$X, Z=data$Z, W=data$W,  perms=perms,  tail = tail, ...)
		#OR ROTATION TEST TO BE ADDED. e cambia anche .othoZ (vedi righe 20-23), che accetti anche NA.
	######MANCA	
	}  
  ##############################################
  } else ##############################################
  ############## Standard -  (NOT missing values in Y)
  {
   #if (Z have at least 1 column and is not constant ) OR (rotationTest)	
   if ((dim(data$Z)[2]>=1) && any(apply(data$Z,2,var)>0))  {
		##orthogocanlize the residuals on Z
		data <- .orthoZ(data$Y,data$X,data$Z)
		}	  
   
   #rotation test
   if(!missing(rotationTest)&&(rotationTest==TRUE)){
		res = .rotation.nptest(data$Y, data$X ,perms=perms,...)
	} else {
   #permutation test
		if(length(unique(data$X))>1){ # if X is not a constant perform dependence.nptest
			data$X=data$X[,apply(data$X,2,var)>0,drop=FALSE]
			res= .dependence.nptest(Y=data$Y, X=data$X, perms=perms,  tail = tail,...)
		} else { #otherwise perform a symmetry test
			res= .symmetry.nptest(Y=data$Y, perms=perms,  tail = tail,...)
		}
	}
 }
  
	#build the flip-object
	out=.getOut(permSpace=res$permSpace,permT=res$permT, data=data,tail=tail, call=call, flipReturn=flipReturn)
  return(out)
}

