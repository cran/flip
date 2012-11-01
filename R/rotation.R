.t.rotation.nptest <- function(){
	#data <- .orthoZ(data)
	N=nrow(data$Y)
	perms <- make.permSpace(1:N,perms,testType=testType)
	
	permT=.prod.perms.rotation(data,perms)
	if(statTest=="sum") permT= .prod2sum(permT,data) else
	if(statTest=="t") permT= .prod2t(permT,data)
	
	res=list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="t"))
}

		
.F.rotation.nptest <- function(){
	#data <- .orthoZ(data)
	N=nrow(data$Y)
	m=ncol(data$Y)
	q=ncol(data$X)
	perms <- make.permSpace(1:N,perms,testType=testType)
	
	P=.get.eigenv.proj.mat(data)
	permT=.prod.perms.P.rotation(data,perms,P)
	permT=.prod2F(permT,data)
	
	res=list(permT=permT,perms=perms,tail=1,extraInfoPre=list(Test="F"))
}


.t.rotation.nptest.1sample <- function(){
				data <- .orthoZ(data)
				
				naRows=apply(is.na(data$Y),1,sum)
				if(any(naRows>0)) {
					warning("Some NA on Y. In rotationTest observations are excluded row-wise")
					Y=Y[naRows==0,,drop=FALSE]
				}
				Ns=nrow(data$Y)
				M2s=apply(data$Y^2,2,sum)
				
				perms <- make.permSpace(1:Ns,perms,testType=testType)
				
				permT=rbind(apply(data$Y,2,sum),
					foreach(i = 1:perms$B,.combine=rbind) %do% { 
						# R is random matrix of independent standard-normal entries 
						# Z shall be a random matrix with the same mean and covariance structure as Y 
						apply(perms$rotFunct(),2,sum)
					}
				)
				colnames(permT) = .getTNames(data$Y)
				rownames(permT)=.getTRowNames(permT)
				permT=permT/t(sqrt((M2s-t((permT)^2)/Ns)*((Ns)/(Ns-1))))
				return(list(permT=permT,perms=perms,tail=tail,extraInfoPre=list(Test="t")))
}	


.prod.perms.rotation <-function(data,perms){
	permT=rbind(as.vector(t(data$X)%*%data$Y),
			foreach(i = 1:perms$B,.combine=rbind) %do% { 
				# R is random matrix of independent standard-normal entries 
				# Z shall be a random matrix with the same mean and covariance structure as Y 
				as.vector(t(data$X)%*%perms$rotFunct())
			}
		)
	colnames(permT)=.getTNames(data$Y,data$X)
	rownames(permT)=.getTRowNames(permT)
	permT
}

.prod.perms.P.rotation <-function(data,perms,P=P){
	require(foreach)				
	.Fvector <- function(Y,P) {
		apply((t(Y)%*%P)^2,1,sum)
	}
	permT=rbind(.Fvector(data$Y,P),
					foreach(i = 1:perms$B,.combine=rbind) %do% { 
						# R is random matrix of independent standard-normal entries 
						# perms$rotFunct() #the function already multiplies the rotation matrix by  Y 
						# Z shall be a random matrix with the same mean and covariance structure as Y 
						.Fvector(perms$rotFunct(),P)
					}
				)				
	colnames(permT)=.getTNames(data$Y,data$X)
	rownames(permT)=.getTRowNames(permT)
	permT

}