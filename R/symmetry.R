############################
# permutation Test of symmetry (i.e. one sample test) 
# Y is the Nxp matrix of respoonces
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of signs be returned?
############################
.symmetry.nptest <- function(Y, perms=5000,  tail = NULL, W=NULL){
	if(!is.matrix(Y)) Y=as.matrix(Y)
	
	Ns=apply(!is.na(Y),2,sum)
	Ms=apply(Y,2,mean,na.rm=TRUE)
	M2s=apply(Y^2,2,mean,na.rm=TRUE)
	Y[is.na(Y)]=0
	permSpace <- make.signSpace(nrow(Y),perms)
	
	if(is.null(W)){
		permT <- permSpace$permID %*% Y
	} else { ### W can be a vector of length nrow(Y) or a matrix of same dim of Y
		permT <- permSpace$permID %*% (W * Y)
	
	}
	
	colnames(permT) = .getTNames(Y)
	permT = rbind(permT,-permT[nrow(permT):1,])
	rownames(permT)=.getTRowNames(permT)
	permT=permT/t(sqrt((M2s-t(permT^2))/(Ns-1)))
	StatType=rep("t",ncol(Y))

#	out=eval.parent(makeFlipObject(),n=2)
#makeFlipObject <- function(){
	# get p-values
	return(list(permT=permT,permSpace=permSpace,extraInfoPre=list(StatType=StatType)))
}
