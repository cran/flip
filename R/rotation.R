.rotation.nptest <- function (Y, X ,perms=1000,tail, separatedX=TRUE,...){
  
	if(!is.matrix(Y)) Y=as.matrix(Y)
	if(!is.matrix(X)) X=as.matrix(X)
	
	
	if(separatedX){
		permT=.randRotate(Y,X,perms)
		StatType=rep("t-rotation",ncol(Y)*ncol(X))
		return(list(permT=permT,extraInfoPre=list(StatType=StatType)))
	} else {
	    print("Option \"separatedX=FALSE\" is not implemented (yet) for rotation tests! (hint: use  \"separatedX=TRUE\" and cobine the results with npc() )")
		
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
.randRotate <- function(Y,X,perms){
	require(foreach)
	perms=.RotSpaceMatchInput(perms)
	if(!is.null(perms$seed)) set.seed(perms$seed)

	n=nrow(Y)
	permT=matrix(NA,perms$number+1,ncol(Y)*ncol(X))
	colnames(permT) = .getTNames(Y,X)
	rownames(permT)=.getTRowNames(permT)

	X=scale(X,center=FALSE)*sqrt(n/(n-1))
	Y=scale(Y,center=FALSE)*sqrt(n/(n-1))
	permT[1,]=as.vector(t(X)%*%Y)
	for (i in 1:perms$number) { 
		# R is random matrix of independent standard-normal entries 
		R <- perms$rotFunct(n)
		# Z shall be a random matrix with the same mean and covariance structure as Y 
		permT[i+1,] <- as.vector(t(X)%*%R%*%Y)
    }
	permT=permT/nrow(Y)
	permT=permT/sqrt(1-(permT^2))*sqrt(nrow(Y)-1)
}