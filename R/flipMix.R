flipMix <- function(modelWithin,X=NULL,Z=NULL,units, perms=1000, data=NULL, tail=NULL,
                    statTest=NULL,flipReturn, testType="permutation", 
                    Su=NULL, equal.se=FALSE,se=NA,replaceNA.coeffWithin=NA,
                    replaceNA.coeffWithin.se=0, ...) {

  otherParams= list(...)
  if(missing(flipReturn)||is.null(flipReturn)) 
  flipReturn=list(permT=TRUE,permP=FALSE,permSpace=FALSE,data=FALSE)
  
  if(is.null(statTest) ) if(is.null(otherParams$separatedX)   || otherParams$separatedX)   { statTest="t" } else statTest="F"
  
  if(testType=="permutation") 
  rotationTest=FALSE  else if(testType=="rotation") rotationTest=TRUE else {
	if(is.null(otherParams$rotationTest) || (!otherParams$rotationTest) ) {testType="permutation"; rotationTest=FALSE } else { testType="rotation"; rotationTest=TRUE}
  }
  
  if(is.null(statTest)) statTest="t"
  
  # store the call
  call <- match.call()
  
  if(!(is.list(data) && (!is.data.frame(data)))) {
     data<-obs2coeffWithin(modelWithin,X=X,Z=Z,units=units, data=data,equal.se=equal.se,se=se,
                        replaceNA.coeffWithin=replaceNA.coeffWithin,replaceNA.coeffWithin.se=replaceNA.coeffWithin.se,...)
 
  }
  rm(Z,X,modelWithin)
  #########
	N = nrow(data$coeffWithin)
	p = ncol(data$coeffWithin)
	############################## Estimate of random effects
	if(is.null(data$covs)) {data$covs=array(,c(N,p,p)); for( id in 1:N) data$covs[id,,]=diag(data$se[id,]^2)}
	if(is.null(Su)){
   		data$Su=.estimateSuMultiILS(Y=data$coeffWithin,Z=as.matrix(cbind(data$X,data$Z)), S=data$covs)
	 } else {
		data$Su=Su; rm(Su) 
		}
	names(data)[names(data)=="coeffWithin"]="Y"
		if(length(unique(unlist(data$X)))>1){ # if X is not a constant perform dependence.nptest
		  res=.between.nptest(data, perms=perms, statTest=statTest[1], tail = tail, testType=testType, ...)
		  out=res$test()
		  #una pezza estetica:
		  if(ncol(data$X)==1)
		    colnames(out$permT)=paste(colnames(out$permT), "_|_",colnames(data$X),sep="")
		  res=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,test=ressub)
      if(length(statTest)>1){
        for(i in 2:length(statTest)) {
          ressub=.between.nptest(data, perms=perms, statTest=statTest[1], tail = tail, testType=testType, ...)
          out=ressub$test()
          #una pezza estetica:
          if(ncol(data$X)==1)
            colnames(out$permT)=paste(colnames(out$permT), "_|_",colnames(data$X),sep="")
          out=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,test=ressub)
          res=cFlip(res,out)          
        }
      }
		} else { #otherwise perform a symmetry test
		  data$W=data$Y
		  data$W[,]=NA
		  for(j in 1:nrow(data$covs)){
		    data$W[j,]=1/sqrt(diag(data$Su) + diag(data$covs[j,,]))
		  }
      for(j in 1:nrow(data$covs)) {  data$W[j,] = diag(data$Su) + diag(data$covs[j,,])} 
		  res=.symmetry.nptest(data, perms=perms, statTest=statTest,  tail = tail,testType="t",...)
		  out=res$test()
		  res=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=res)
      
      if(length(statTest)>1){
			for (i in 2:length(statTest)){
				data$W=matrix(,dim(data$Y))
        dimnames(data$W)=dimnames(data$Y)
		    for(j in 1:nrow(data$covs)){  
          data$W[j,]=1/sqrt(diag(data$Su) + diag(data$covs[j,,]))
		    }
			  ressub=.symmetry.nptest(data, perms=perms, statTest=statTest,  tail = tail,testType="t",...)
				out=ressub$test()
				out=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=res)
        res=cFlip(res,out)
			}
		}
	}
return(res)
}