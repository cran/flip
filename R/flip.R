#######################
flip.statTest <-
    c("t", "F", "ANOVA",
	"Wilcoxon","Kruskal-Wallis", "kruskal", "rank", "Mann-Whitney",
	"chisq","chisq.separated", "Fisher",
	#"KS", "kolmogorow", "Kolmogorow-Smirnov", "ad",
  "McNemar", "Sign","sum")

.get.statTest <- function(statTest){ 
	if(is(statTest,"function")) return(statTest) else
	
	statTest <- match.arg(tolower(statTest[1]),tolower(flip.statTest))
	statTest= flip.statTest[which(statTest==tolower(flip.statTest))]
	#synonyms
	if(statTest=="ANOVA") 
		statTest="F" else
	if(statTest=="kruskal") 
		statTest="Kruskal-Wallis" else
	if(statTest=="Mann-Whitney")
		statTest="Wilcoxon" else
# 	if(statTest%in%c("KS", "kolmogorow"))
# 		statTest="Kolmogorow-Smirnov"
# 		
	statTest
}

#########################
flip <- function(Y, X=NULL, Z=NULL, data=NULL, tail = 0, perms = 1000, statTest=NULL, Strata=NULL, flipReturn, testType=c("permutation","rotation"), ...) {

  if(missing(flipReturn)||is.null(flipReturn)) 
  flipReturn=list(permT=TRUE,permP=FALSE,permSpace=FALSE,data=FALSE)

  if(is.null(statTest) ) if(is.null(list(...)$separatedX)   || list(...)$separatedX)   { statTest="t" } else statTest="F"
    statTest <- .get.statTest(statTest)
	
  if(is.null(testType)){
	if(is.null(list(...)$rotationTest) || (!list(...)$rotationTest) ) {testType="permutation"; rotationTest=FALSE } else { testType="rotation"; rotationTest=TRUE} 
  } 
  testType=match.arg(testType,c("permutation","rotation"))
  rotationTest= (testType=="rotation")



  # store the call
  call <- match.call()
  # get matrices from inputs
  data <- .getXY(Y,X,Z,data,rotationTest=rotationTest,dummyfy=list(...)$dummyfy,statTest=statTest,Strata=Strata)
  rm(X,Y,Z,Strata)
  
  symmetryTest= is.null(data$X) || (length(unique(data$X))==1)

  #check if the problem can be set as one sample problem
  if(!symmetryTest)
	if(statTest%in% c("t","sum","rank","Wilcoxon","McNemar","Sign"))
	  if(  !is.null(data$Strata) ){#is.null(data$Z)|| ncol(data$Z)==0)  &
			keep=setdiff(1:ncol(data$X),.getIntercept(data$X))
			if( (length(unique(data$X[,keep]))==2) && 
				(ncol(data$X[,keep,drop=FALSE])==1) )
					if(all(table(data$X[,keep],unlist(data$Strata))==1)){
						attrsYassign=attributes(data$Y)$assign
						attrsYfactors=attributes(data$Y)$factors
            
						data$X=data$X[,keep,drop=FALSE]
						levs=unique(data$X)
						data$Y=t(sapply(unique(unlist(data$Strata)), function(ids){
              data$Y[(data$Strata==ids)&(data$X==levs[2]),]-
              data$Y[(data$Strata==ids)&(data$X==levs[1]),]}))
            
						attributes(data$Y)$assign=attrsYassign
						attributes(data$Y)$factors=attrsYfactors
						data$X=NULL
						data$Strata=NULL
						data$Z=NULL
						symmetryTest=TRUE
					}	
		}
  
  # if symmetry.nptest
  if(symmetryTest){
  		test= .symmetry.nptest(data, perms=perms, statTest=statTest,  tail = tail,testType=testType,...)
  ##dependence.nptest
  } else 
	if ( !(any(is.na(data$Y))|| ifelse(is.null(data$X),TRUE,any(is.na(data$X))))){ # standard solutions, not missing data
		test= .dependence.nptest(data, perms=perms,statTest=statTest,  tail = tail,testType=testType,...)
	} else {	stop("Warning: NA values are not allowed, nothing done.")	}
	res <- test$test()
	#build the flip-object
	res=.getOut(res=res,data=data, call=call, flipReturn=flipReturn)
  return(res)
}