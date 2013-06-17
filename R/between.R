.between.nptest <- function(data, perms=5000,  tail = NULL,statTest="t",testType="permutation",...){
	
	if(is.null(statTest)) 
		if(is.null(otherParams$separatedX) || otherParams$separatedX)
			statTest="t" else
			statTest="F"
	
	if(!is.null(otherParams$rotationTest)) {
    if(otherParams$rotationTest) { testType="rotation"; rotationTest=TRUE} else rotationTest=FALSE
	} else 
	
	if(is.function(statTest)) {
		test <- statTest
	} else
	if(statTest=="t"){
		test <- .t.between.nptest		
	} else  
	if(statTest=="F"){ #ANOVAtype test, 1 column for each column of Y summarizing the dependence with all Xs			
	    test <- .F.between.nptest
	} else  
	if(statTest=="Trace"){ #ANOVAtype test, 1 column for each column of Y summarizing the dependence with all Xs			
		test <- .trace.between.nptest
	} else {warning("This test is not valid, nothing done."); test <- function() return()}

  environment(test) <- sys.frame(sys.nframe())
  out <- sys.frame(sys.nframe())
  return(out)
}

###########################################
.t.between.nptest <- function(){
	if(ncol(data$Y)==1) {
		data$W=as.matrix(1/sqrt(data$Su + (data$covs)))
	} else {
		data=.getLinCombAndWeights(data,otherParams)
	}
	
	#colnames(data$W)=colnames(data$Y)
    perms <- make.permSpace(1:(nrow(data$Y) - ncol(data$Z)),perms,return.permIDs=TRUE,testType=testType)
	uni.test <- function(i,data){	
		data <- .orthoZ(list(X=data$W[,i]*data$X,Y=data$W[,i]*data$Y[,i,drop=FALSE],Z=data$W[,i]*data$Z,intercept = FALSE))
		#environment(.prod2t) <- sys.frame(sys.parent())
   	permT=.prod.perms(data,perms,testType=testType)
		permT=.prod2t(permT,data)
		permT
	}
	permT = matrix(, perms$B+1,ncol(data$Y))
  for(i in 1:ncol(data$Y)) permT[,i]=uni.test(i,data)
  colnames(permT)=.getTNames(data$Y)
  rownames(permT)=.getTRowNames(permT)
	return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="w-t")))
}
	
	
# #########################################
.trace.between.nptest <- function(){

	if(ncol(data$Y)==1) { warning("Only one response variable is inputed, use F statistic. Trace test is not performed.")
		return(NULL)
	} else {
	
  covs=data$covs
  for(i in 1:nrow(data$covs)) covs=data$covs[i,,] +data$Su
  
	perms <- make.permSpace(covs,perms,testType="Simulation")	
	
	permT=.trace.sim(data,perms)
	rownames(permT)=.getTRowNames(permT)
	return(list(permT=permT,perms=perms,tail=1,extraInfoPre=list(Test="F-trace")))
	}
}

.trace.sim <- function(data,perms){
	
	PZX=cbind(data$Z,data$X)%*%solve(t(cbind(data$Z,data$X))%*%cbind(data$Z,data$X))%*%t(cbind(data$Z,data$X))
	PZ=data$Z%*%solve(t(data$Z)%*%data$Z)%*%t(data$Z)
	PZXPZ= PZX-PZ
	HZX= diag(nrow(PZX)) - PZX
	rm(PZ,PZX)
	dfratio=(nrow(data$Y)- ncol(data$X) -ncol(data$Z) )/ncol(data$X)
	.stat <- function(y) .tr(t(y)%*%PZXPZ %*% y)/.tr(t(y)%*%HZX%*%y) * dfratio
  
	permT=matrix(,perms$B+1,1)
  permT[1,]=.stat(data$Y)
	for(i in 1:perms$B) permT[i+1,]=.stat(perms$rotFunct())
	colnames(permT)="F-trace"
	permT
}

#################################################
.F.between.nptest <- function(){
  
  if(ncol(data$Y)==1) 
    data$W=as.matrix(1/sqrt(data$Su + (data$covs))) else 
      data=.getLinCombAndWeights(data, otherParams)
  
  perms <- make.permSpace(1:(nrow(data$Y) - ncol(data$Z)),perms,return.permIDs=TRUE,testType=testType)
  uni.test <- function(i,data){	
    data <- .orthoZ(list(X=data$W[,i]*data$X,Y=data$W[,i]*data$Y[,i,drop=FALSE],Z=data$W[,i]*data$Z,intercept = FALSE))
    P=.get.eigenv.proj.mat(data)
    permT=.prod.perms.P(data,perms,testType=testType,P)
    permT=.prod2F(permT,data)
    permT
  }
  permT = matrix(,perms$B+1,ncol(data$Y))
  for(i in 1:ncol(data$Y)) permT[,i]= uni.test(i,data)
  colnames(permT)=.getTNames(data$Y)
  rownames(permT)=.getTRowNames(permT)
  return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="w-F")))
}


#####################################
.getLinCombAndWeights<- function(data, otherParams){
	if(is.null(otherParams$fastSumCombination)){
		if(is.null(otherParams$linComb)) {
			if(is.null(otherParams$whichpcs)) whichpcs=1:(1+(ncol(data$Y)==2)) else whichpcs=otherParams$whichpcs
			# standard deviations
			sd2s=diag(apply(data$covs,c(2,3),mean))+diag(data$Su)
			linComb=1/sqrt(sd2s)
			CORR=array(t(apply(data$covs,1,function(cov){ diag(linComb)%*%(cov+data$Su)%*%diag(linComb)})),dim(data$covs))
			CORR=apply(CORR,c(2,3),mean)
			linComb = (linComb)*prcomp(CORR)$rotation[,whichpcs,drop=FALSE]
			linComb = cbind(linComb,sum=1)
			
		} else 
			linComb=otherParams$linComb
			
		data$Y =cbind(data$Y,comb=data$Y%*%linComb)
		data$W=array(,dim(data$Y))
		rownames(data$W)=rownames(data$Y)  
		for(j in 1:nrow(data$covs)) {
		  cov= data$Su + data$covs[j,,]
      data$W[j,]=1/sqrt(c(diag(cov), diag(t(linComb)%*%cov%*%linComb)))
		}	
	} else{
		if(is.null(otherParams$linComb)) {
			if(is.null(otherParams$whichpcs)) whichpcs=1:(1+(ncol(data$Y)==2)) else whichpcs=otherParams$whichpcs
			meanCov=apply(data$covs,c(2,3),mean)
			linComb= diag(1/sqrt(diag(meanCov)))%*%prcomp(meanCov,scale. = TRUE)$rotation[,whichpcs,drop=FALSE]
			linComb = cbind(linComb,sum=1)
		} else 
			linComb=cbind(otherParams$linComb,1)
			
		comb=data$Y%*%linComb
		covMulti=array(t(apply(data$covs,1,function(cov){ diag(t(linComb)%*%cov%*%linComb)})),c(nrow(data$covs),ncol(linComb)))
		SuMulti=sapply(1:ncol(comb),function(i) .estimateSuMultiILS(Y=comb[,i,drop=FALSE],Z=as.matrix(cbind(data$X,data$Z)), S=array(covMulti[,i,drop=FALSE],c(nrow(comb),1,1))))
		data$Y =cbind(data$Y,comb=comb)
		data$W=array(,dim(data$Y))
		rownames(data$W)=rownames(data$Y)  
		for(j in 1:nrow(data$covs)) {
		  cov= data$Su + data$covs[j,,]
		  data$W[j,]=1/sqrt(c(diag(cov), covMulti[j,]+SuMulti))
		}	
		colnames(data$W)=colnames(data$Y)
	}
data
}
