flipMix <- function(modelWithin,X=NULL,Z=NULL,units, perms=1000, data=NULL, tail=NULL,
                    statTest=NULL,flipReturn, testType="permutation", 
                    Su=NULL, equal.se=FALSE,se=NA,replaceNA.coeffWithin=NA,
                    replaceNA.coeffWithin.se=0, ...) {

  otherParams= list(...)
  require(foreach)
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
			res=foreach (i = 1:length(statTest),.combine=cFlip) %do% {
				ressub=.between.nptest(data, perms=perms, statTest=statTest[i], tail = tail, testType=testType, ...)
				out=ressub$test()
        #una pezza estetica:
        if(ncol(data$X)==1)
          colnames(out$permT)=paste(colnames(out$permT), "_|_",colnames(data$X),sep="")
				out=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,test=ressub)
			}
		} else { #otherwise perform a symmetry test
			res=foreach (i = 1:length(statTest),.combine=cFlip) %do% {
			  sd=sqrt(   foreach(j= 1:nrow(data$covs),.combine=rbind) %do% {  cov= diag(data$Su) + diag(data$covs[j,,])} )
				data$Y <- data$Y/sd
			  ressub=.symmetry.nptest(data, perms=perms, statTest=statTest,  tail = tail,testType="t",...)
				out=ressub$test()
				out=.getOut(res=out,data=data, call=call, flipReturn=flipReturn,call.env=ressub)
			}
		}
return(out)
}