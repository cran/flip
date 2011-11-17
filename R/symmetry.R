############################
# permutation Test of symmetry (i.e. one sample test) 
# Y is the Nxp matrix of respoonces
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of signs be returned?
############################
.symmetry.nptest <- function(Y, perms=5000,  tail = NULL){
	if(!is.matrix(Y)) Y=as.matrix(Y)
	if(is.null(colnames(Y))) colnames(Y) = paste("Y",sep="",if(ncol(Y)>1)1:ncol(Y))
	
	permSpace <- make.signSpace(dim(Y)[1],perms)
	
    permT <- permSpace$permID %*% Y
	colnames(permT)=colnames(Y)
    permT = rbind(permT,-permT)
	tStats=list(t=permT[1,]/sqrt(apply(Y,2,var)/(dim(Y)[1]-1)))

#	out=eval.parent(makeFlipObject(),n=2)
#makeFlipObject <- function(){
	# get p-values
	return(list(permT=permT,permSpace=permSpace,tStats=tStats))
}
