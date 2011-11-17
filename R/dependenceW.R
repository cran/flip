############################
# permutation Test of dependence 
# Y is the Nxp matrix of responces
# X is the Nxq matrix of predictors
# Z is the Nx? matrix of covariates
# W is the Nxp matrix of weights to be apply to Y X, and Z (separated columns Y)
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of permutations (IDs) be returned?
############################

.dependenceW.nptest <- function(Y, X, Z, W=NULL, perms=5000,  tail = 0 ,subsets, rotationTest=FALSE, flipReturn=list(permT=TRUE)){
	require(foreach)
	if(!is.matrix(Y)) Y=as.matrix(Y)
	if(!is.matrix(X)) X=as.matrix(X)
	
	colnames(Y)=.getYNames(Y)	
	colnames(X)=.getXNames(X)	
	if(is.null(colnames(W))) colnames(W) = colnames(Y)
	
	if(!missing(Z)){
		if(ncol(Z)==0) Z=matrix(rep(1,nrow(Y)))
		if(!is.matrix(Z)) Z=as.matrix(Z)
		if(is.null(colnames(Z))) colnames(Z) = paste("Z",sep="",if(ncol(Z)>1)1:ncol(Z))
	}
	
	if(missing(subsets)|| is.null(subsets)) {subsets=as.list(colnames(Y)); names(subsets)=colnames(Y)}
	
	N=dim(Y)[1] - ifelse(missing(Z),0,dim(Z)[2]) ###occhio, questo è un po' pericoloso
	
	if(!rotationTest)	{
		permSpace <- make.permSpace(1:N,perms); rm(perms)
		permT <- matrix(NA,nrow(permSpace$permID),length(subsets))
	} else {
		perms
		permT <- matrix(NA,perms,length(subsets))
	}
	colnames(permT) = names(subsets)
	rownames(permT)=.getTRowNames(permT)

	for(i in 1:length(subsets)){
		diagW=if(length(subsets[[i]])>1) diag( apply(W[,subsets[[i]]]^-2,1,sum)^-.5 ) else diag( W[,subsets[[i]]])
		if(!missing(Z)){
		  D=.orthoZ(X=diagW%*%X,Y=diagW%*%Y[,subsets[[i]] ,drop=FALSE],Z=diagW%*%Z)
		} else{
		  D=list(X=diagW%*%X,Y=diagW%*%Y[,subsets[[i]],drop=FALSE])
		}
		if(!rotationTest) {
			res=.dependence.nptest(D$Y, D$X, perms=permSpace,  tail = tail, separatedX=FALSE)
		} else res=.rotation.nptest(D$Y, D$X, perms=perms,  tail = NULL, separatedX=FALSE)
		
		permT[,i]=res$permT
		
	}
	permT=round(permT,10)  #########occhio qui, con numeri molto piccoli possono verificarsi problemi	
	if(!rotationTest) StatType=rep("pseudoF",ncol(Y))  else StatType=rep("pseudoF-rotation",ncol(Y))
		
	return(list(permT=permT,permSpace=res$permSpace,extraInfoPre=list(StatType=StatType)))
}