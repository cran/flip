
flip.adjust <- function (permTP, method = flip.npc.methods, maxalpha=1, weights=NULL, stdSpace=FALSE, ...) {
    
	##TODO: here is case sensitive, while in npc it is not. here is so for compatibility with p.adjust. shall we change?
	method <- match.arg(method,c(p.adjust.methods,flip.npc.methods))
	
	
	##TODO: add some check for names(weights) compatibility. 
	##TODO: exclude missing names (eg permTP=permTP[names(weights)]
	if(!is.null(weights)) if(is.null(names(weights))) names(weights)<- colnames(permTP)
	
	if(is(permTP,"flip.object")){
		if((method%in%p.adjust.methods) & (is.null(weights))){ 
			#run the standard function
			adjs=p.adjust(p.value(permTP), method = method)
		} else if((method%in%p.adjust.methods) & (!is.null(weights))){
			#standard function but weigthed
			if(method=="hommel") {print("weighted \"hommel\" method not allowed."); return()}
			adjs=p.adjust.w(p.value(permTP), method = method, w=weights)
		} else {		#perform permutation specific procedures
			if(method=="minP") {
				adjs=.maxt.adjust(if(is.null(permTP@permP)) -t2p(permTP@permT,obs.only=FALSE,tail=permTP@tail) else -permTP@permP,  maxalpha,weights=weights)
			} else if(method=="maxT") {
				if(stdSpace) {permTP@permT = .t2stdt(permTP@permT,FALSE)}
				adjs=.maxt.adjust(.setTail(permTP@permT, tail=permTP@tail), maxalpha,weights=weights)
			} else {
			#otherwise perform closed testing
			if(stdSpace & (method %in% c("sumT", "sumT2"))) {permTP@permT = .t2stdt(permTP@permT,FALSE)}
			# Define the local test to be used in the closed testing procedure
			mytest <- function(hyps) {p.value(npc(permTP[hyps],method,weights=weights[hyps]))}
			cl <- closed(mytest, names(permTP),alpha=NA)
			adjs=sapply(names(permTP),function(id) adjusted(cl,id))
			}
		}
		#fit it to the flip-object
		permTP@res=cbind(permTP@res,adjs)
		colnames(permTP@res)[length(colnames(permTP@res))]=paste("Adjust:",method,sep="")	
	
		return(permTP)
			
	} else { # not a flip.object, return a vector
		if((method%in%p.adjust.methods) & (is.null(weights))){ 
			return(p.adjust(permTP, method = method))
		} else if((method%in%p.adjust.methods) & (!is.null(weights))){
			#standard function but weigthed
			if(method=="hommel") {print("weighted \"hommel\" method not allowed."); return()}
			return(p.adjust.w(permTP, method = method, w=weights))
		} else{		#perform permutation specific procedures
			if(method=="minp") {
				return(.maxt.adjust(-permTP, maxalpha,weights=weights))
			} else if(method=="maxtstd") { 
				return(.maxt.adjust(scale(permTP),  maxalpha,weights=weights))
			} else if(method=="maxT") {
				return(.maxt.adjust(permTP, maxalpha,weights=weights))
			} else {
			#then perform closed testing
			# Define the local test to be used in the closed testing procedure
			mytest <- function(hyps) {p.value(npc(permTP[,hyps],method,weights=weights[hyps]))}
			cl <- closed(mytest, colnames(permTP),adjust=TRUE)
			adjs=cl@adjusted[1:ncol(permTP)]
			}
		}
	}
	
}


################
.maxt.adjust <- function(permT,maxalpha=1,weights=NULL,m=ncol(permT)) {
	
	#m=numb of hypos

	#get colnames to 
	if(is.null(colnames(permT))) colnames(permT)=1:m
	
	#define the order of testing hypotheses
	steps=names(sort(permT[1,],decreasing = TRUE))
	
	if(!is.null(weights)) permT=t(weights*t(permT))
	
	#set of NOT yet rejected p-values
	notrejs=rep(TRUE,m)
	names(notrejs)=colnames(permT)
	
	#set of adjusted p-values
	Padjs=rep(1,m)
	names(Padjs)=colnames(permT)
	
	i <- 1
	while((i<=m) & ifelse(i>1,Padjs[steps[i-1]] <= maxalpha,TRUE)){
		Padjs[steps[i]]=max( t2p(c(permT[1,steps[i]], apply(permT[-1,notrejs,drop=FALSE],1,max) )) ,Padjs[steps[i-1]] ) #first max ensures monocinity
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
