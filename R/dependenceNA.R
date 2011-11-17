############################
# permutation Test of dependence 
# Y is the Nxp matrix of responces
# X is the Nxq matrix of predictors
# perms: number of permutations
# tail : vector of tails 1, -1 or 0
# permP.return, permT.return, permSpace.return : logical: shoul space of p-values, of statistic and of permutations (IDs) be returned?
############################

.dependenceNA.nptest <- function(Y, X, perms=5000,  tail = NULL){
	require(foreach)
	if(!is.matrix(Y)) Y=as.matrix(Y)
	colnames(Y)=.getYNames(Y)
	#contr=contrasts(levels(X), contrasts =FALSE)
	#contr[contr==0]=-1
	X=data.frame(X=factor(apply(X,1,paste,collapse="_")))
	
	if (nlevels(X[,1])==2) nlevs2=TRUE else nlevs2=FALSE
	X=.makeContrasts(~.-1,data=X,excludeIntercept=TRUE,excludeRefCat=FALSE)
	
	### when only 2 categories are prestent, the first is omittted
	if(nlevs2) X=X[,-1,drop=FALSE]
	
	colnames(X)=.getXNames(X)
	
	N=nrow(Y)
	permSpace <- make.permSpace(1:N,perms); rm(perms)
    
	notnaY=!is.na(Y)
	Y[is.na(Y)]=0	
	sumNotNaYs= apply(notnaY,2,sum)
	sumYs= apply(Y,2,sum)
	
	Tcontrast=function(x,i){
		ni1s=t(x)%*%notnaY[permSpace$permID[i,],]
		coefs1= as.vector(sqrt((sumNotNaYs-ni1s)/ni1s))
		sumYs1 =t(x)%*%Y[permSpace$permID[i,],]
		sumYs1*coefs1 - (sumYs-sumYs1)*(coefs1^-1) 
		######alternativa, fuori: sumT= apply(Y,2,sum)
		# ni1s=t(x)%*%notnaY[permSpace$permID[i,],]
		# coefs1= as.vector(sqrt((notNAys-ni1s)/ni1s))
		# sumY1=apply(Y[permSpace$permID[i,][x==1],],2,sum)
		# return(sumY1%*%(coefs1) - (sumY-sumY1)*(coefs1^-1) )
		}
		i=1
	permT=foreach( i =1:nrow(permSpace$permID) ,.combine=rbind)	%do% {
			as.vector(t(apply(X,2,Tcontrast,i)))
			}
			
	permT=scale(permT,scale=FALSE)
	colnames(permT)=.getTNames(Y,X)
	rownames(permT)=.getTRowNames(permT)
	tStats=list(F=permT[1,])

	
	StatType=rep("meanDiff-NA",ncol(Y)*ncol(X))

	return(list(permT=permT,permSpace=permSpace,extraInfoPre=list(StatType=StatType)))
}