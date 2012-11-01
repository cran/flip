############################
# permutation Test of dependence 
# Y is the Nxp matrix of responces
# X is the Nxq matrix of predictors
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of permutations (IDs) be returned?
############################

.dependence.nptest <- function(data, perms=5000,  tail = NULL,statTest="t",separatedX=TRUE,testType="permutation",...){
	
	if(is.null(statTest)) 
		if(is.null(separatedX) || separatedX) 
			statTest="t" else
			statTest="F"
	
	if(is.function(statTest)) {
		test <- statTest
	} else {
		### mettere qui la pulizia degli statTest per dependence test
		if(statTest%in%c("sum","t")){
			if(testType=="rotation"){
				test <- .t.rotation.nptest
			} else # not a rotation test
				test <- .t.dependence.nptest		
		} else  
		if(statTest=="F"){ #ANOVAtype test, 1 column for each column of Y summarizing the dependence with all Xs			
			if(testType=="rotation"){
				test <- .F.rotation.nptest
			} else #not a rotation test
				test <- .F.dependence.nptest
		} else
		if(statTest=="Fisher"){ 
			if(testType=="rotation"){
				warning("Rotations are not allowed for Fisher exact test, permutations will be used instead.")
				testType="permutation"
			} else #not a rotation test
				test <- .fisher.dependence.nptest
		} else if(statTest%in%c("Wilcoxon","Kruskal-Wallis","ranks")){
			if(testType=="rotation"){
				warning("Rotations are not allowed for Wilcoxon (i.e. ranks) test, permutations will be used instead.")
				testType="permutation"
			} else ## permutation test
				test <- .rank.dependence.nptest
		} else if(statTest%in%c("chisq","chisq.separated")){
			test <- .chisq.dependence.nptest
		} else if(statTest%in%"Kolmogorov-Smirnov"){
			test <- .kolmogorov.dependence.nptest
		} else	{stop("This test statistic is not valid, nothing done."); return()}
		if(statTest%in%c("Fisher","Wilcoxon","Kruskal-Wallis","ranks","chisq","Kolmogorov-Smirnov")) 
			if(length(unique(data$Z))>1 ) warning("Covariates Z can not be used in this test. Use strata instread.")
	}
  environment(test) <- sys.frame(sys.nframe())
  out <- sys.frame(sys.nframe())
  return(out)
}

###########################################

.t.dependence.nptest <- function(){
	#data <- .orthoZ(data) #if Z is.null, it doesn't make anything
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE)
	#search for intercept in the model
	intercept=.getIntercept(data$X)
	if(any(intercept)) {
		data$intercept=TRUE
		data$X=data$X[,!intercept,drop=FALSE]
	}
		#data$X=scale(data$X,scale=FALSE)
	permT=.prod.perms(data,perms,testType)
	
	if(statTest=="sum") 
		permT= .prod2sum(permT,data)
	else if(statTest=="t") {
		permT= .prod2t(permT,data)
		if(any(intercept) ) permT =permT * sqrt((N-1-sum(intercept))/(N-1))
		}

	colnames(permT) = .getTNames(data$Y,data$X,permT=permT)
	rownames(permT)=.getTRowNames(permT)		
	
	if(statTest=="sum") {
		center=as.vector(t(colSums(data$X))%*%colSums(data$Y))/N
		names(center)=colnames(permT)
		attributes(tail)$center=center
	}
	res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test=statTest))
}
			
.F.dependence.nptest <- function(){
	#data <- .orthoZ(data)
	
	N=nrow(data$Y)	
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE)
	#search for intercept in the model
	intercept=.getIntercept(data$X) 
	if(any(intercept)) {
		data$X=data$X[,!intercept,drop=FALSE]
		data$X=scale(data$X,scale=FALSE)
		data$Y=scale(data$Y,scale=FALSE)
	}
	
	q=ncol(data$X)#-sum(intercept) #just a trick, it should be subtracetd to N, not to q
	
#	environment(.prod2F) <- sys.frame(sys.nframe())
	P=.get.eigenv.proj.mat(data)
	permT=.prod.perms.P(data,perms,testType,P)
	permT=.prod2F(permT,data)
	if(any(intercept)) permT =permT * (N-q-sum(intercept))/(N-q)
	
	colnames(permT) = .getTNames(data$Y,permT=permT)
	rownames(permT)=.getTRowNames(permT)		
	res=list(permT=permT,perms=perms,tail=1,extraInfoPre=list(Test="F"))
}

#################################
.rank.dependence.nptest <- function(){
	data$Y=apply(data$Y,2,rank)
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE)
	
	#search for intercept in the model
	intercept=.getIntercept(data$X)
	if(any(intercept)) data$X=data$X[,!intercept,drop=FALSE]

	if(statTest=="ranks")
			statTest=ifelse(ncol(data$X)>1,"Kruskal-Wallis","Wilcoxon")

	if(statTest=="Wilcoxon"){
#		data$X=scale(data$X,scale=FALSE)
		permT=.prod.perms(data,perms)
		permT=scale(permT,center=rep(colSums(data$X)*(nrow(data$X)+1)/2,ncol(data$Y)))
		colnames(permT)=.getTNames(data$Y,data$X)
		rownames(permT)=.getTRowNames(permT)
		res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="Wilcoxon"))
	} else { #kruskall-wallis (anova)
		data$X=scale(data$X,scale=FALSE)
		data$Y=scale(data$Y,scale=FALSE)
		P=.get.eigenv.proj.mat(data)
		permT=.prod.perms.P(data,perms,testType="permutation",P)
		permT=.prod2F(permT,data)
		permT=permT *ncol(data$X)/(nrow(data$Y)-ncol(data$X)) 
		permT = (1+ permT^-1)^-1 *(nrow(data$Y)-1)
		res=list(permT=permT,perms=perms,tail=1,extraInfoPre=list(Test="Kruskal-Wallis")) 
		res
	}
	
}
#################################
.fisher.dependence.nptest <- function(){
	if(any(apply(data$Y,2,function(y) length(unique(y))!=2 ))) stop("Only factor or dichotomous variables are allowed with Fisher exact test. Nothing done.")
	
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE)
	
	#search for intercept in the model
	intercept=.getIntercept(data$X)
	if(any(intercept)) data$X=data$X[,!intercept,drop=FALSE]

#	data$X=scale(data$X,scale=FALSE)
	permT=.prod.perms(data,perms)
	
	colnames(permT)=.getTNames(data$Y)
	rownames(permT)=.getTRowNames(permT)		
	center=as.vector(colSums(data$X)%*%t(colSums(data$Y)))/N
	names(center)=colnames(permT)
	attributes(tail)<-list(center=center)
	res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="Fisher"))
}
################################
# .kolmogorov.nptest<-function(){
#   
#   }
####################################
.chisq.dependence.nptest<-function(){ 
	if(any(apply(data$Y,2,function(y) length(unique(y))!=2 ))) stop("Only factor or dichotomous variables are allowed with Chi Squared test. Nothing done.")
	
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,return.permIDs=TRUE)
	
	#search for intercept in the model
	intercept=.getIntercept(data$X)
	if(any(intercept)) { 
		attrsX=attributes(data$X)[c("assign","factors")]
		data$X=data$X[,!intercept,drop=FALSE]
		attrsX$assign=attrsX$assign[attrsX$assign>0]
		attributes(data$X)[c("assign","factors")] <- attrsX
	}
#	data$X=scale(data$X,scale=FALSE)
	permT=.prod.perms(data,perms)
	expected=as.vector(colSums(data$X)%*%t(colSums(data$Y)))/N
	permT=scale(permT,center=expected,scale=sqrt(expected))
	if(statTest=="chisq.separated"){
		colnames(permT)=.getTNames(data$Y,data$X) 
	} else {
		newNames=apply(expand.grid(attributes(data$X)$factors,attributes(data$Y)$factors),1,paste,collapse="_|_")
		permT=sapply(unique(attributes(data$Y)$assign),
				function(idy) {
					whichCol= as.vector(outer(1:ncol(data$X),(which(attributes(data$Y)$assign==idy)-1)*ncol(data$X),"+"))
					rowSums(permT[,whichCol,drop=FALSE]^2)
					})
		colnames(permT)=newNames
	}
	rownames(permT)=.getTRowNames(permT)		
	res=list(permT=permT,perms=perms,tail=ifelse(statTest=="chisq",1,tail),extraInfoPre=list(Test="Chi Squared"))
}