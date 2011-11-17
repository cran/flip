############################

############################

npc <- function(permTP, comb.funct = c("Fisher", "Liptak", "minP", "maxT", "sumT", "sumT2"),subsets=NULL,weights=NULL,...){
	
	### just in analogy with gt(). to be implemented as flip-options
	trace=TRUE
	
	#### arrange comb.funct
	comb.funct <- match.arg(tolower(comb.funct[1]), c("fisher", "liptak", "minp", "tippett", "maxt", "sumt", "direct", "sumt2"))
	if(comb.funct == "tippett") comb.funct ="minp"
	if(comb.funct == "direct") comb.funct ="sumt"
	
	if(is(permTP,"flip.object")){
		nperms = permTP@call$perms
		if(comb.funct %in% c("fisher", "liptak", "minp")) {
			if(!is.null(permTP@permP)){ 
				permTP=permTP@permP			
			} else { 
				if(!is.null(permTP@permT)) { permTP=t2p(permTP@permT,obs.only=FALSE,tail=permTP@tail)}
			} 
		} else if(comb.funct %in% c("maxt", "sumt", "sumt2")) {permTP=.setTail(permTP@permT,tail=.fitTail(permTP@permT,permTP@tail) ) }
	}
	
	if(!is.matrix(permTP)) permTP=as.matrix(permTP)
	if(!exists("nperm"))  nperms = list(number=dim(permTP)[1],seed=NA)
	
	if(comb.funct=="fisher"){permTP = -log(permTP)}
	if(comb.funct=="liptak"){permTP = -qnorm(permTP)}
	if(comb.funct=="minp"){permTP = 1/permTP}
	if(comb.funct=="sumt2"){permTP = permTP^2}
	
	
	
	############
	temp=.getSubsetWeights(weights,subsets,colnames(permTP))
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
      if (length(weights) != ncol(permTP)) 
        stop("length of \"weights\" does not match column count of \"permTP\"")
      all.weights <- weights
    } else
      all.weights <- rep(1, ncol(permTP))
	  names(all.weights) = colnames(permTP)
  
  
  
	if(comb.funct %in% c("fisher", "liptak", "sumt", "sumt2")){
		  test= function(subset=NULL,weights=NULL){ 
		  permT = matrix(if(is.null(subset)) permTP%*%all.weights else permTP[,subset,drop=FALSE]%*%all.weights[subset]) ;permT}
	} else 	if(comb.funct %in% c("minp", "maxt"))
		   test= function(subset=NULL,weights=NULL){ #browser()
		   permT = matrix(apply(if(is.null(subset)) t(all.weights*t(permTP)) else 
		   t(all.weights[subset]*t(permTP[,subset,drop=FALSE])) , 1, max))  ;permT}

	
	
	# Do the test
  if ((!many.subsets) && (!many.weights)) {           # single weighting; single subset
    permT <- test()
	nVar=dim(permTP)[2]
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

	#build the flip-object
	  if(!exists("flipReturn") || is.null(flipReturn)) 
			flipReturn=list(permT=TRUE,permP=FALSE)
  
	out=.getOut(type="npc",permSpace=NULL,permP=p,permT=permT, data=list(),tail=NULL, call=match.call(), flipReturn=flipReturn,comb.funct=comb.funct,nVar=nVar)
	
	return(out)
}