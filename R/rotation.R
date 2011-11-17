.rotation.nptest <- function (Y, X ,B=1000,separatedX=TRUE,...){
  
	if(!is.matrix(Y)) Y=as.matrix(Y)
	if(!is.matrix(X)) X=as.matrix(X)
	if(is.null(colnames(Y))) colnames(Y) = paste("Y",sep="",if(ncol(Y)>1)1:ncol(Y))
	if(is.null(colnames(X))) colnames(X) = paste("X",sep="",if(ncol(X)>1)1:ncol(X))
	
	
	if(separatedX){
		permT=.randRotate(Y,X,B)
		tStats=list(t=permT[1,]/sqrt( as.vector(outer(apply(X^2,2,sum),apply(Y^2,2,sum)))-(permT[1,]^2))*sqrt(dim(Y)[1]-2))
		return(list(permT=permT,tStats=tStats))	
	} else {
	    print("Option \"separatedX=FALSE\" is not implemented (yet) for rotation tests!")
		
		# P = X%*%solve(t(X)%*%X)%*%t(X)
		# #IminusP= diag(N)-P
		# permT=sapply(1:dim(Y)[2],function(col) {
				# STATS=foreach( i =1:dim(permSpace$permID)[1] ,.combine=c)	%do% {
				# t(Y[permSpace$permID[i,],col])%*%P%*%Y[permSpace$permID[i,],col]
				# }
				# STATS=STATS/(sum(Y[,col]^2)-STATS)
			# }
		# )
		# permT=permT*(N-dim(X)[2])/dim(X)[1]
		# permT=round(permT,10)  #########occhio qui, con numeri molto piccoli possono verificarsi problemi
		# colnames(permT) = colnames(Y)
		# tail="F"
	}
	
}


################ the core of the procedure
.randRotate <- function(Y,X,B){
require(foreach)
	T=matrix(NA,B+1,dim(Y)[2]*dim(X)[2])
	T[1,]=as.vector(t(X)%*%Y)
	for (i in 1:B) { 
		# R is random matrix of independent standard-normal entries 
		R <- matrix(rnorm(dim(Y)[1]^2),ncol=dim(Y)[1]) 
		R <- qr.Q(qr(R, LAPACK = TRUE)) 
		# Z shall be a random matrix with the same mean and covariance structure as Y 
		T[i+1,] <- as.vector(t(X)%*%R%*%Y)
    }
	T
}