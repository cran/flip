############################

############################

npc <- function(permSpace, comb.funct = c("Fisher", "Liptak", "minP", "maxT", "sumT", "sumT2"),permP.return=FALSE,permT.return=TRUE,subsets=NULL,weights=NULL){
	
	### just in analogy with gt(). to be implemented as flip-options
	trace=TRUE
	
	#### arrange comb.funct
	comb.funct <- match.arg(tolower(comb.funct[1]), c("fisher", "liptak", "minp", "tippett", "maxt", "sumt", "direct", "sumt2"))
	if(comb.funct == "tippett") comb.funct ="minp"
	if(comb.funct == "direct") comb.funct ="sumt"
	
	if(is(permSpace,"flip.object")){
		nperms = permSpace@call$perms
		if(comb.funct %in% c("fisher", "liptak", "minp")) {
			if(!is.null(permSpace@permP)){ 
				permSpace=permSpace@permP			
			} else { 
				if(!is.null(permSpace@permT)) permSpace=t2p(permSpace@permT,obs.only=FALSE)
			} 
		} else if(comb.funct %in% c("maxt", "sumt", "sumt2")) permSpace=permSpace@permT
	}
	
	if(!is.matrix(permSpace)) permSpace=as.matrix(permSpace)
	if(!exists("nperm"))  nperms = list(number=dim(permSpace)[1],seed=NA)
	
	if(comb.funct=="fisher"){permSpace = -log(permSpace)}
	if(comb.funct=="liptak"){permSpace = -qnorm(permSpace)}
	if(comb.funct=="minp"){permSpace = 1/permSpace}
	if(comb.funct=="sumt2"){permSpace = permSpace^2}
	
	
	
	############
	temp=.getSubsetWeights(weights,subsets,colnames(permSpace))
	weights=temp$weights
	subsets=temp$subsets
	many.weights=temp$many.weights
	many.subsets=temp$many.subsets
	one.weight=temp$one.weight
	
	rm(temp)
	
	# prepare progress info
   ###  if (missing(trace)) trace <- gt.options()$trace && (many.weights || many.subsets)
  if (trace && (many.subsets || many.weights)) {
    if (many.subsets) 
      K <- length(subsets) 
    else 
      K <- length(weights) 
    digitsK <- trunc(log10(K))+1
  }

  
  # weight
    if (one.weight) {
      if (length(weights) != ncol(permSpace)) 
        stop("length of \"weights\" does not match column count of \"permSpace\"")
      all.weights <- weights
    } else
      all.weights <- rep(1, ncol(permSpace))
	  names(all.weights) = colnames(permSpace)
  
  
  
	if(comb.funct %in% c("fisher", "liptak", "sumt", "sumt2")){
		  test= function(subset=NULL,weights=NULL){ 
		  permT = matrix(if(is.null(subset)) permSpace%*%all.weights else permSpace[,subset,drop=FALSE]%*%all.weights[subset]) ;permT}
	} else 	if(comb.funct %in% c("minp", "maxt"))
		   test= function(subset=NULL,weights=NULL){ #browser()
		   permT = matrix(apply(if(is.null(subset)) t(all.weights*t(permSpace)) else 
		   t(all.weights[subset]*t(permSpace[,subset,drop=FALSE])) , 1, max))  ;permT}

	
	
	# Do the test
  if ((!many.subsets) && (!many.weights)) {           # single weighting; single subset
    permT <- test()
	nVar=dim(permSpace)[2]
  } else {     
    L <- if (many.subsets) length(subsets) else length(weights)                 
    permT <- sapply(1:L, function (i) { 
      if (trace && L>1) {
        cat(rep("\b", 2*digitsK+3), i, " / ", K, sep="")
        flush.console()
      }
      if (!many.weights) {                                           # single weighting; many subsets
        uit <- test(subset=subsets[[i]]) 
      } else if (!many.subsets) {                                    # many weightings; single subset
        uit <- test(weights=weights[[i]]) 
      } else {                                                      # many weightings; many subsets
        uit <- test(subset=subsets[[i]], weights=weights[[i]])
      } 
	  uit
    })
    if (many.subsets && !is.null(names(subsets))){
      colnames(permT) <- names(subsets)
	  nVar=sapply(subsets,length)
	}
    else if (many.weights && !is.null(names(weights))){
      colnames(permT) <- names(weights)
	  nVar=sapply(weights,length)
	}
  }
  if (trace && (many.subsets || many.weights) && L>1) cat("\n")
	
	
	
#	out=eval.parent(makeFlipObject(),n=2)
#makeFlipObject <- function(){
	# get p-values
	p=t2p(permT,obs.only=TRUE)
	# test statistic and std dev
	stat=permT[1,]
	stDev=apply(permT,2,sd)
	
	#build the results table
	res=data.frame(Stat=as.vector(stat),Std.dev=as.vector(stDev), Z=as.vector(stat/stDev), p=as.vector(p),combFunct=comb.funct, nvar=nVar)
	colnames(res)[colnames(res)=="nvar"]="#Vars"
	colnames(res)[colnames(res)=="p"]="p-value"
	rownames(res)=colnames(permT)

	out <- new("flip.object")  
    out @res = res
	out @nperms = nperms
	out @call = match.call()
	out @permP=if(permP.return) t2p(permT, obs.only=FALSE)
	out @permT=if(permT.return) permT
#	out
#}
	return(out)
}