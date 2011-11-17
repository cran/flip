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
	colnames(Y)=.getYNames(Y)	
	
	N=nrow(Y)
	permSpace <- make.permSpace(1:N,perms); rm(perms)
    
	#permY=array(NA,c(dim(permSpace$permID)[1],dim(Y)[1],dim(Y)[2]))
	
	if(separatedX){
		permT=matrix(NA,dim(permSpace$permID)[1],dim(Y)[2]*dim(X)[2])
		colnames(X)=.getYNames(X)
		colnames(permT) = .getTNames(Y,X)

		X=scale(X)
		Y=scale(Y)
		for(i in 1:ncol(Y)){
			#permY[,,i] <- matrix(Y[permSpace$permID,i],nrow=dim(permSpace$permID)[1])
			permT[,.getTNames(Y[,i,drop=FALSE],X) ] <-(matrix(Y[permSpace$permID,i],nrow=dim(permSpace$permID)[1]) %*% X) #permY[,,i] %*% X
		}
		## transformed in a t-statistic
		permT=permT/N
		permT=permT/sqrt(1-(permT^2))*sqrt(N)
		permT=round(permT,10)  #########occhio qui, con numeri molto piccoli possono verificarsi problemi
		#paste(sort(rep(colnames(data$Y),2)),colnames(data$X),sep="_|_")
		#colnames(permT)=colnames(Y)
		#permT = .setTail(permT,tail)
		#get output
		#	out=eval.parent(makeFlipObject(),n=2)
		#makeFlipObject <- function(){
		# get p-values
		StatType=rep("t",ncol(Y)*ncol(X))
		
	} else { #ANOVAtype test, 1 colums for each columns of Y summarizing the dependence with all Xs
		P = X%*%solve(t(X)%*%X)%*%t(X)
		#IminusP= diag(N)-P
		permT=sapply(1:ncol(Y),function(col) {
				STATS=foreach( i =1:dim(permSpace$permID)[1] ,.combine=c)	%do% {
				.tr(t(Y[permSpace$permID[i,],col])%*%P%*%Y[permSpace$permID[i,],col])
				}
				STATS=STATS/(sum(sum(Y[,col]^2))-STATS)
			}
		)
		permT=permT*(N-dim(X)[2])/dim(X)[1]
		permT=round(permT,10)  #########occhio qui, con numeri molto piccoli possono verificarsi problemi
		colnames(permT)=.getTNames(Y,X)
		rownames(permT)=.getTRowNames(permT)
		StatType=rep("F",ncol(Y))
	}
	
	return(list(permT=permT,permSpace=permSpace,extraInfoPre=list(StatType=StatType)))
}