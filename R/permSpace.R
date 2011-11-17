
######## match the setting for permutation 
.PermSpaceMatchInput <- function(perms) {
if (!is.list(perms)) {
		#the whole matrix of random permutaitons is provided
		if (is.matrix(perms)) return(list(permID = perms, number=dim(perms)[1] ,seed=NA))		
		#only the number of random permutaitons is provided
		if (is.numeric(perms)) return(list(permID = NULL, number=perms ,seed=NA))
		
      } else return(perms)  #hopefully there are 3 elements, migliorare la funzione qui
}	  



############################
# Calculates signs flips of a vector of N elements.
# perms is the number of flips; if perms > number of all possible flips, then compute the complete space
############################
make.signSpace <- function(N,perms) {
    perms=.PermSpaceMatchInput(perms)
	if(is.null(perms$permID)){
		if (2^(N-1) <= perms$number) {
		    # all permutations if possible and if no stratas
			#random <- FALSE
			if(N>1){
				require(e1071)
				perms$permID <-cbind(0,bincombinations(N-1))
				perms$permID [which(perms$permID ==0)] <- 1/N
				perms$permID [which(perms$permID ==1)] <- -1/N
			} else if(N==1) {
				perms$permID <-rbind(1,-1)
			}	else if(N==0) {
				perms$permID <-matrix(0,2,0)
			}
				perms$seed=NA
				perms$number=2^(N)
		} else {
			#otherwise random permutations
			if (is.na(perms$seed)) perms$seed <- round(runif(1)*1000)
			set.seed(perms$seed)
			perms$permID <- rbind(rep(1, N), matrix(1 - 2 * rbinom(N * (perms$number%/%2),1, 0.5), (perms$number%/%2), N))/N
		}	
	}
	perms
}

############################
# Calculates permutations space of a vector Y. 
# perms is the number of permutations; if perms > number of all permutations, then compute the complete space
############################
make.permSpace <- function(Y,perms) {
    perms=.PermSpaceMatchInput(perms)
	if(is.null(perms$permID)){
		allperms=npermutations(Y)
		# all permutations if possible
		if ( allperms <= perms$number) {
			#random <- FALSE
			perms$permID <- t(allpermutations(Y))
			perms$seed=NaN
			perms$number=allperms
		} else {
			# otherwise random permutations
			if (is.na(perms$seed)) perms$seed <- round(runif(1)*1000)            
			set.seed(perms$seed)
			N <- length(Y)
			perms$permID <- t(cbind(unlist(Y), replicate(perms$number, Y[sample(N)])))
		}
  }
  perms
}

############################
# rotation space of a vector Y. 
# perms is the number of permutations; seed the seed for random number generation, rotFunct the function to generate the random rotations
############################
.RotSpaceMatchInput <- function(perms) {

	rotFunct <- function(n){
		R <- matrix(rnorm(n^2),ncol=n) 
		R <- qr.Q(qr(R, LAPACK = TRUE)) 
		return(R)
	}

	if (is.list(perms)) {
		perms <- perms[c("number","seed","rotFunct")]
		if(is.null(perms$number))  perms$number <- 1000
		if(is.null(perms$seed) || is.na(perms$seed) )  perms$seed <- round(runif(1)*1000)
		if(is.null(perms$rotFunct))  perms$rotFunct <- rotFunct
		 return(perms)
	}
	
	#only the number of random permutaitons is provided
	if (is.numeric(perms)) return(list(rotFunct = rotFunct, number=perms, seed= round(runif(1)*1000)))	
}	  


		
############################
# Iterative function calculates all permutations of a vector
# values: vector of all unique values
# multiplicity: multiplicity of each value
############################
.allpermutations <- function(values, multiplicity) {

  if (length(values)==1) {
    out <- values
  } else {
    n <- sum(multiplicity)
    out <- matrix(0 , n, .npermutations(multiplicity))
    where <- 0
    for (i in 1:length(values)) {
      if (multiplicity[i] == 1) {
        newmult <- multiplicity[-i]
        newvals <- values[-i]
      } else {
        newmult <- multiplicity
        newmult[i] <- newmult[i] - 1
        newvals <- values
      }
      range <- where + seq_len(.npermutations(newmult))
      where <- range[length(range)]
      out[1,range] <- values[i]
      out[2:n, range] <- .allpermutations(newvals, newmult)
    }
  }
  out
}

############################
# Iterative function counts all permutations of a vector
# values: vector of all unique values
# multiplicity: multiplicity of each value
############################
.npermutations <- function(multiplicity) {
  round(exp(lfactorial(sum(multiplicity)) - sum(lfactorial(multiplicity))))
}

############################
# Counts all permutations of a vector y
# user-friendly version of .npermutations()
############################
npermutations <- function(Y) {
  .npermutations(table(Y))
}

############################
# Calculates all permutations of a vector y
# user-friendly version of .allpermutations()
############################
allpermutations <- function(Y) {
  values <- unique(Y)
  multiplicity <- colSums(outer(Y, values, "=="))
  .allpermutations(values, multiplicity)
}


##########
#compute p-value space P from statistic space T (the percentile of the statistic T column-wise)
t2p<-function(T,obs.only=TRUE,tail = 1){  
    
	if(!missing(tail))	T = .setTail(T,tail)
	
	if(!is.matrix(T)) {T<-as.matrix(T)}
	if(obs.only) { 
		P=matrix(apply(T,2,function(permy)mean(permy>=permy[1])),1,dim(T)[2])
		rownames(P)="p-value"
	}
	else{
		#oth<-seq(1:length(dim(T)))[-1]
		B<-dim(T)[1]
		P=apply(-T,2,rank,ties.method ="max")/B
		P=as.matrix(P)
		rownames(P)=c("p-obs",paste("p-*",1:(B-1),sep=""))
	}
	colnames(P)=colnames(T)
	return(P)
}


############### standardize permT space. used in maxTstd and
.t2stdt <- function(permT,obs.only=TRUE){ return(t(t((permT[1:(nrow(permT)^(!obs.only)),]))/apply(permT,2,sd,na.rm=TRUE)))}
