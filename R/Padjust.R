
flip.adjust <- function (permSpace, method = "maxT", maxalpha=1, weights=NULL,...) {
    
	##TODO: here is case sensitive, while in npc it is not. here is so for compatibility with p.adjust. shall we change?
	method <- match.arg(method,c(p.adjust.methods,"maxT","minP","maxTstd"))  
	
	##TODO: add some check for names(weights) compatibility. 
	##TODO: exclude missing names (eg permSpace=permSpace[names(weights)]
	if(!is.null(weights)) if(is.null(names(weights))) names(weights)<- colnames(permSpace)
	
	if(is(permSpace,"flip.object")){
		if((method%in%p.adjust.methods) & (is.null(weights))){ 
			#run the standard function
			adjs=p.adjust(p.value(permSpace), method = method)
		} else if((method%in%p.adjust.methods) & (!is.null(weights))){
			#standard function but weigthed
			adjs=p.adjust.w(p.value(permSpace), method = method, w=weights)
		} else{		#perform permutation specific procedures
			if(method=="minP") {
				adjs=.maxt.adjust(if(is.null(permSpace@permP)) -t2p(permSpace@permT,obs.only=FALSE) else -permSpace@permP, method = method, maxalpha,weights=weights)
			} else if(method=="maxTstd") { 
				adjs=.maxt.adjust(scale(permSpace@permT), method = method, maxalpha,weights=weights)
			} else {
				adjs=.maxt.adjust(permSpace@permT, method = method, maxalpha,weights=weights)
			}
		}
		#fit it to the flip-object
		permSpace@res=cbind(permSpace@res,adjs)
		colnames(permSpace@res)[length(colnames(permSpace@res))]=paste("Adjust:",method,sep="")	
	
		return(permSpace)
			
	} else { # not a flip.object, return a vector
		if((method%in%p.adjust.methods) & (is.null(weights))){ 
			return(p.adjust(permSpace, method = method))
		} else if((method%in%p.adjust.methods) & (!is.null(weights))){
			#standard function but weigthed
			return(p.adjust.w(permSpace, method = method, w=weights))
		} else{		#perform permutation specific procedures
			if(method=="minP") {
				return(.maxt.adjust(-permSpace, method = method, maxalpha,weights=weights))
			} else if(method=="maxTstd") { 
				return(.maxt.adjust(scale(permSpace), method = method, maxalpha,weights=weights))
			} else {
				return(.maxt.adjust(permSpace, method = method, maxalpha,weights=weights))
			}
			
		}
	}
	
}


################
.maxt.adjust <- function(permT,method,maxalpha=1,weights=NULL) {
	
	#numb of hypos
	m=dim(permT)[2]
	
	#get colnames to 
	if(is.null(colnames(permT))) colnames(permT)=1:m
	
	#define the order of testing hypotheses
	steps=names(sort(-permT[1,]))
	
	if(!is.null(weights)) permT=t(weights*t(permT))
	
	#set of NOT yet rejected p-values
	notrejs=rep(TRUE,m)
	names(notrejs)=colnames(permT)
	
	#set of adjusted p-values
	Padjs=rep(1,m)
	names(Padjs)=colnames(permT)
	
	i <- 1
	while((i<=m) & ifelse(i>1,Padjs[steps[i-1]] <= maxalpha,TRUE)){
		Padjs[steps[i]]=max( t2p(c(permT[1,steps[i]], apply(permT[-1,notrejs,drop=FALSE],1,max) )) ,Padjs[steps[i-1]] ) #first max ensure monocinity
		notrejs[steps[i]]=FALSE
		
		#avoid to compute max if test statistic are equal (specially useful in minp)
		while(permT[1,steps[i]]==ifelse(i==m,Inf,permT[1,steps[i+1]])){  
			Padjs[steps[i+1]]=Padjs[steps[i]]
			notrejs[steps[i+1]]=FALSE
			i=i+1
		}
		
		i=i+1
	}

	return(Padjs)
}
