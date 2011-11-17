############################
# permutation Test of dependence 
# Y is the Nxp matrix of responces
# X is the Nxq matrix of predictors
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of permutations (IDs) be returned?
############################

.dependence.nptest <- function(Y, X, perms=5000,  tail = NULL,separatedX=TRUE){
	require(foreach)
	if(!is.matrix(Y)) Y=as.matrix(Y)
	if(!is.matrix(X)) X=as.matrix(X)
	if(is.null(colnames(Y))) colnames(Y) = paste("Y",sep="",if(ncol(Y)>1)1:ncol(Y))
	if(is.null(colnames(X))) colnames(X) = paste("X",sep="",if(ncol(X)>1)1:ncol(X))
	
	N=dim(Y)[1]
	permSpace <- make.permSpace(1:N,perms); rm(perms)
    
	#permY=array(NA,c(dim(permSpace$permID)[1],dim(Y)[1],dim(Y)[2]))
	
	if(separatedX){
		permT=matrix(NA,dim(permSpace$permID)[1],dim(Y)[2]*dim(X)[2])
		colnames(permT) = paste(rep(colnames(Y),rep(ncol(X),ncol(Y))),colnames(X),sep="_|_")
		X=scale(X)
		Y=scale(Y)
		for(i in colnames(Y)){
			#permY[,,i] <- matrix(Y[permSpace$permID,i],nrow=dim(permSpace$permID)[1])
			permT[,paste(i,colnames(X),sep="_|_")] <-(matrix(Y[permSpace$permID,i],nrow=dim(permSpace$permID)[1]) %*% X)/(N-1) #permY[,,i] %*% X
			## transformed in a t-statistic
#			permT[,paste(i,colnames(X),sep="_|_")] <- permT[,paste(i,colnames(X),sep="_|_")]/sqrt(1-permT[,paste(i,colnames(X),sep="_|_")]^2)*sqrt(N-2)
		}
		permT=round(permT,10)  #########occhio qui, con numeri molto piccoli possono verificarsi problemi
		#paste(sort(rep(colnames(data$Y),2)),colnames(data$X),sep="_|_")
		#colnames(permT)=colnames(Y)
		#permT = .setTail(permT,tail)
		#get output
		#	out=eval.parent(makeFlipObject(),n=2)
		#makeFlipObject <- function(){
		# get p-values
		tStats=list(t=permT[1,]/sqrt(rep(apply(Y,2,var),rep(dim(X)[2],dim(Y)[2]))-permT[1,])*sqrt(dim(Y)[1]-2))
		
	} else {
		P = X%*%solve(t(X)%*%X)%*%t(X)
		#IminusP= diag(N)-P
		permT=sapply(1:dim(Y)[2],function(col) {
				STATS=foreach( i =1:dim(permSpace$permID)[1] ,.combine=c)	%do% {
				t(Y[permSpace$permID[i,],col])%*%P%*%Y[permSpace$permID[i,],col]
				}
				STATS=STATS/(sum(Y[,col]^2)-STATS)
			}
		)
		permT=permT*(N-dim(X)[2])/dim(X)[1]
		permT=round(permT,10)  #########occhio qui, con numeri molto piccoli possono verificarsi problemi
		colnames(permT) = colnames(Y)
		tStats=list(F=permT[1,])
	}
	
	return(list(permT=permT,permSpace=permSpace,tStats=tStats))	
}