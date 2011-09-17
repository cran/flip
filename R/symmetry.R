############################
# permutation Test of symmetry (i.e. one sample test) 
# Y is the Nxp matrix of respoonces
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of signs be returned?
############################
.symmetry.nptest <- function(Y, perms=5000,  tail = NULL,permT.return=TRUE,permP.return=FALSE,permSpace.return=FALSE){
	if(!is.matrix(Y)) Y=as.matrix(Y)
	if(is.null(colnames(Y))) colnames(Y) = paste("Y",sep="",if(ncol(Y)>1)1:ncol(Y))
	
	permSpace <- make.signSpace(dim(Y)[1],perms)
	
    permT <- permSpace$permID %*% Y
	colnames(permT)=colnames(Y)
    
    permT = .setTail(rbind(permT,-permT),tail)

#	out=eval.parent(makeFlipObject(),n=2)
#makeFlipObject <- function(){
	# get p-values
	p=t2p(permT,obs.only=TRUE)
	
	# test statistic and std dev
	stat=permT[1,]
	stDev=apply(permT,2,sd)
	
	# tails of the test
	dirNum=sign(c(1,-1)%*%.setTail(matrix(c(1,-1),2,dim(permT)[2]),tail))
	dir=rep("",length(dirNum))
	stat[dirNum<0]=-stat[dirNum<0]
	dir[dirNum==1]=">"
	dir[dirNum==-1]="<"
	dir[dirNum==0]="><"
	#build the results table
	res=data.frame(Stat=as.vector(stat),Std.dev=as.vector(stDev), Z=as.vector(stat/stDev), p=as.vector(p),tail=as.vector(dir))
	rownames(res)=colnames(permT)
	colnames(res)[colnames(res)=="p"]="p-value"

	out <- new("flip.object")  
    out @res = res
	out @nperms = permSpace[c("number","seed")]
	out @call = match.call()
	out @call$perms = permSpace[c("seed","number")]
	out @permP=if(permP.return) t2p(permT, obs.only=FALSE)
	out @permT=if(permT.return) permT
	out @permSpace=if(permSpace.return) permSpace$permID
	out @tail = tail
#	out
#}
	
	return(out)
}
