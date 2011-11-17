#==========================================================
# CLASS DEFINITION *** CLASS DEFINITION *** CLASS DEFINITION
#==========================================================

#setClassUnion("matrixOrNULL", c("matrix", "NULL"))
setClassUnion("numericOrmatrixOrNULL", c("numeric","matrix", "NULL"))
setClassUnion("arrayOrNULL", c("array", "NULL"))
setClassUnion("data.frameOrNULL", c("data.frame", "NULL"))
setClassUnion("numericOrmatrixOrcharacterOrNULL", c("numeric","matrix", "NULL","character"))

#############da togliere per compilazione (esistno gia in someMTP)
#setClassUnion("numericOrNULL", c("numeric", "NULL"))
#setClassUnion("listOrNULL", c("list", "NULL"))


setClass("flip.object", 
  representation(
    res = "data.frameOrNULL",
    call = "call", 
	permP="arrayOrNULL",
	permT="arrayOrNULL",
	permSpace="arrayOrNULL",
	permY="arrayOrNULL",
	nperms="listOrNULL",
    #functions = "environment",#"list",
    #subsets = "listOrNULL",
    #structure = "listOrNULL",
    #weights = "listOrNULL",
    tail = "numericOrmatrixOrcharacterOrNULL",
    #Z = "matrixOrNULL",
    #directional = "logical",
    data = "listOrNULL"
    #model = "character"
  ),
  prototype = list(
    res = NULL,
	nperms=NULL,
	permP=NULL,
	permT=NULL,
	permSpace=NULL,
	permY=NULL,
	data=NULL
  )
)

#==========================================================
# Function "show" prints a "flip.object" object
#==========================================================
setMethod("show", "flip.object", function(object)
{
  result(object)
})

setGeneric("summary")
setMethod("summary", "flip.object", function(object, ...)
{
  nperms= as.list(object@call$perms)
  cat(" \"flip.object\" object of package flip\n")
  cat(" Call:\n ")
  cat(deparse(object@call), "\n")
  cat(ifelse(is.null(nperms$seed),"all",""), nperms$number, ifelse(is.finite(nperms$seed),"random",""), "permutations.\n",   
	ifelse(is.finite(nperms$seed),paste("(seed: ",is.finite(nperms$seed)," )",sep=""),"")) 
  cat("\n")
  show(object)
})


#==========================================================
# Functions to extract relevant information from 
# a flip.object object
#==========================================================
setGeneric("result", function(object, ...) standardGeneric("result"))
setMethod("result", "flip.object",
  function(object) { 
  print(object@res, digits = 3)  
})

# setGeneric(".result", function(object, ...) standardGeneric(".result"))
# setMethod(".result", "flip.object",
  # function(object) 
    # object@result
# )

# setGeneric("weights") 
# setMethod("weights", "flip.object", function(object) {
  
  # # base weights based on covariate variance
  # X <- object@functions$getX()
  # weights <- colSums(X*X)
  # if (object@model == "multinomial")
    # weights <- rowSums(matrix(weights, object@functions$df()[3]))
  # names(weights) <- object@functions$cov.names()
                                  
  # # find weights for specific weights and subsets chosen
  # if (length(object@subsets) > 0) {
    # weights <- lapply(object@subsets, function(set) weights[set])
  # }
  # if (length(object@weights) > 0) 
    # if (is.list(weights)) {
      # weights <- lapply(as.list(1:length(weights)), function(i)
        # weights[[i]] * object@weights[[i]])
      # names(weights) <- names(object@weights)  
    # } else
      # weights <- lapply(object@weights, function(wts) weights * wts)
                  
  # # set the max weight to 1  
  # if (is.list(weights))
    # weights <- lapply(weights, function(wts) wts / max(wts))
  # else
    # weights <- weights / max(weights)
      
  # # reduce a list of length 1 to a vector    
  # if (is.list(weights) && length(weights) == 1)
    # weights <- weights[[1]]

  # weights      
# })

# #==========================================================
# setGeneric("subsets", function(object, ...) standardGeneric("subsets"))
# setMethod("subsets", "flip.object", function(object, ...) {
  # object@subsets
# })


#==========================================================
setGeneric("p.value", function(object, ...) standardGeneric("p.value"))
setMethod("p.value", "flip.object",
  function(object) {
    x=object@res[,"p-value"] 
	names(x)=names(object)
	x
	
  }
)


#==========================================================
setGeneric("size", function(object, ...) standardGeneric("size"))
setMethod("size", "flip.object",
  function(object) {
    dim(object@res,1)
  }
)
# ==========================================================
# setGeneric("dim", function(object, ...) standardGeneric("dim"))
# setMethod("dim", "flip.object",
  # function(object) {
    # c(object@nperms$number , dim(object@res)[1])
  # }
# )


#==========================================================
# The subsetting methods for "flip.object"
#==========================================================
setMethod("[", "flip.object", 
            function(x, i, j,...,drop) 
{
	if(is.character(i) && !all(i %in% names(x))){ 
		search=which(!(i %in% names(x)))
		extended= lapply(i, function(ii) names(x)[if(ii %in% names(x)) ii else grep(ii, names(x))] )
		i=unlist(extended)
	}
  if (all(i %in% names(x)) || 
          all(i %in% 1:length(x)) ||
          all(i %in% -1:-length(x)) ||
          (is.logical(i) && (length(i)== length(x)))) {
    x@res <- x@res[i, ,drop=FALSE]
    if (!is.null(x@permP)) x@permP <- x@permP[,i,drop=FALSE]
    if (!is.null(x@permT)) x@permT <- x@permT[,i,drop=FALSE]
    if (!is.null(x@permY)) x@permY <- x@permY[,,i,drop=FALSE]
	if (!is.null(x@tail)) x@tail <- x@tail[ifelse(length(as.vector(x@tail))==1,1,i)]
    #if (!is.null(x@weights)) x@weights <- x@weights[i]
    x
  } else {
    stop("invalid index set", call. = FALSE)
  }
})            

setMethod("[[", "flip.object", 
            function(x, i, j,...,exact) 
{
   x[i]
})

#==========================================================
# The length method for "flip.object"
#==========================================================
setMethod("length", "flip.object", 
            function(x) 
{
  dim(x@res)[1]
})



#==========================================================
# The names and alias methods for "flip.object" 
# (applies to pathwaynames)
#==========================================================
setMethod("names", "flip.object", 
            function(x) 
{
  rownames(x@res)
})      


setMethod("names<-", "flip.object", 
            function(x, value) 
{
  rownames(x@res) <- value
  if (!is.null(x@permP)) colnames(x@permP) <- value
  if (!is.null(x@permT)) colnames(x@permT) <- value
  x
})            



#==========================================================
# A sort method for "flip.object"
#==========================================================
setGeneric("sort") 
setMethod("sort", "flip.object",
  function(x, decreasing = FALSE ) {
      ix <- order(p.value(x), decreasing=decreasing)
    x[ix]
  }
)

# #==========================================================
# # Model.matrix extracts the model matrix (only if x=TRUE)
# #==========================================================
# setMethod("model.matrix", matchSignature(signature(object = "flip.object"), model.matrix),
  # function(object, ... ) {
    # list(tail = object@tail, Z = object@Z)
  # }
# )


#==========================================================
# Multiple testing correction for "flip.object" object
#==========================================================

setGeneric("p.adjust", function(p, method = p.adjust.methods, n = length(p)) standardGeneric("p.adjust"))
setMethod("p.adjust", matchSignature(signature(p = "flip.object"), p.adjust),
  function(p, method = p.adjust.methods, n = length(p)) {
    method <- method[1]
    method <- p.adjust.methods[grep(method, p.adjust.methods, ign=T)]
    if(length(method)==(0))   # this is just to get a good error message
      method <- match.arg(method)
	if (missing(n))
      p@res <- cbind(p@res, p.adjust(p.value(p), method=method))
    else
      p@res <- cbind(p@res, p.adjust(p.value(p), method=method, n=n))
	
	colnames(p@res)[length(colnames(p@res))]=paste("Adjust:",method,sep="")	
	
    p
  }
)



#==========================================================
# Histogram method to visualize permutations
#==========================================================
setGeneric("hist", function(x,...) standardGeneric("hist"))
#setMethod("hist", matchSignature(signature(x = "flip.object"), hist), 
setMethod("hist", "flip.object", function(x, ...)  {

  flip.hist <- function(x, breaks=20, main="", xlab = "Permutation test statistics", ...) {

     if (length(x) > 1)
     stop("length(object) > 1. Please reduce to a single test result")
  
    # if (is.null(x@weights)) 
      # weights <- rep(1, size(x))
    # else
      # weights <- x@weights[[1]]
    # if (is.null(x@subsets))
      # subset <- seq_len(size(x))
    # else
      # subset <- x@subsets[[1]]
  
    #recalculate <- x@functions$permutations(subset, weights)
    Qs <- x@permT[-1,] #recalculate$permS
    Q <- x@permT[1,]
    nperm <- length(Qs)
    hst <- hist(Qs, xlim = c(1.1 * min(0, Qs, Q), 1.1 * max(Qs, Q)), breaks = breaks, 
      main = main, xlab = xlab, ...)
    h <- max(hst$counts)
    arrows( Q, h/2, Q, 0 , lwd=2)
    text( Q, h/2, 'Observed\ntest\nstatistic' , pos=3)
  
    # No output
    invisible(list(statistic = Q, histogram = hst))
  }
  
  flip.hist(x,...)
})        


# #==========================================================
# # Graph plot for flip-oject
# #==========================================================


setGeneric("plot", function(x, y, ...) standardGeneric("plot"))
setMethod("plot", "flip.object", 
 function(x, y, ...) {
#setMethod("plot", "flip.object", function(x, y, ...) {
  if(!exists("main")) main=NULL 
  if(!exists("xlab")) xlab = NULL
  if(!exists("ylab")) ylab=NULL 
  
   plot.flip <- function(x, y=NULL, main, xlab, ylab,...){
   #draw <- function(x, main, xlab, ylab,...){
 if (length(x)==1 ){
	hist(x, breaks=20, ...)
 } else if (length(x)==2 ){
	plot(x@permT[,1],x@permT[,2],lwd=1,pty="o",xlab=colnames(x@permT)[1],ylab=colnames(x@permT)[2])
	points(x@permT[1,1],x@permT[1,2],col="red",lwd=3)
	title("Bivariate Permutation Space") 
 } else { 
	pc=prcomp(x@permT[,apply(x@permT,2,var)>0],scale. =TRUE,center=FALSE)
	pc$rotation[,1]=pc$rotation[,1]*sign(pc$x[1,1]) 
	pc$rotation[,2]=pc$rotation[,2]*sign(pc$x[1,2]) 
	pc$x[,1]=pc$x[,1]*sign(pc$x[1,1])
	pc$x[,2]=pc$x[,2]*sign(pc$x[1,2]) 
	biplot(pc,xlabs=c("obs",rep("*",dim(pc$x)[1]-1)),
	xlab=paste("PC1 (",round(pc$ sdev [1]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),
	ylab=paste("PC2 (",round(pc$ sdev [2]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),
	main= "PCA of Permutation Space" )
    # lam <- pc$sdev[1:2] #* sqrt(dim(pc$x)[1])
    # #plot(pc$x[, 1:2]/lam)
	# pc$x[,1:2]=pc$x[,1:2] / lam
	# pc$rotation[,1:2]=pc$rotation[,1:2]*lam
	# plot(pc$x[,1],pc$x[,2],lwd=1,pty="o",xlim=range(pc$x[,1])*1.2,ylim=range(pc$x[,2])*1.2,
	# xlab=paste("PC1 (",round(pc$ sdev [1]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),
	# ylab=paste("PC2 (",round(pc$ sdev [2]^2 /sum(pc$ sdev ^2) *100,2)," %)",sep=""),col="gray",pch=21,bg="gray")
	# points(pc$x[1,1],pc$x[1,2],col="red",lwd=3,pch=21,bg="red")
	# text(pc$x[1,1]*1.1,pc$x[1,2]*1.1,col="red","Obs")
	# arrows( 0, 0, 2*pc$rotation[,1], 2*pc$rotation[,2], lwd=1,col="gray")
	# text(2.1*pc$rotation[,1], 2.1*pc$rotation[,2], rownames(pc$rotation), cex=1.5,col="black")
	# title("PCA of Permutation Space") 
	}
  }
   plot.flip(x,y=NULL, main=main, xlab=xlab, ylab=ylab,...)
})

# #==========================================================
# # Graph plot for focus level and inheritance procedures
# #==========================================================
# draw <- function(object, alpha=0.05, type = c("focuslevel","inheritance"), names=FALSE, sign.only = FALSE, interactive = FALSE) {

  # # check availablity of packages
  # require("Rgraphviz") || stop("package \"Rgraphviz\" is not available.")

  # # find ancestors and offspring if missing
  # if (is.null(object@structure$ancestors)) {     # Infer from sets
    # sets <- object@subsets
    # ancestors <- new.env(hash=TRUE)
    # offspring <- new.env(hash=TRUE)
    # for (i in 1:length(sets)) {
      # namei <- names(sets)[i]
      # for (j in 1:length(sets)) {
        # namej <- names(sets)[j]
        # if (i != j && length(sets[[i]]) <= length(sets[[j]]) && all(sets[[i]] %in% sets[[j]])) {
          # ancestors[[namei]] <- c(ancestors[[namei]], namej)
          # offspring[[namej]] <- c(offspring[[namej]], namei)
        # }
      # }
    # }
    # object@structure$ancestors <- as.list(ancestors)
    # object@structure$offspring <- as.list(offspring)
  # }

  # # find type if missing
  # if (missing(type)) 
    # type <- c("focuslevel","inheritance")
  # else
    # type <- match.arg(type)  
  # type <- type[type %in% names(object@extra)]
  # if (length(type) < 1)
    # stop("no focus level or inheritance p-values in object.")
  # if (length(type) > 1)
    # stop("both focus level and inheritance p-values in object. Please specify type.")

  # ps <- object@extra[[type]]
  # significant <- names(object)[ps <= alpha]

  # parents <- ancestors2parents(object@structure$ancestors)

  # # make the graph object
  # graph <- as.matrix(sapply(names(object), function(node) names(object) %in% parents[[node]]))
  # rownames(graph) <- colnames(graph) <- names(object)
  # if (sign.only) {
    # graph <- graph[significant, significant,drop=FALSE]
    # if (length(significant)==0)
      # stop("no significant nodes to plot.")
  # }
  # graph <- as(graph, "graphNEL")
  
  # nodes <- buildNodeList(graph)
  # edges <- buildEdgeList(graph)
  # nAttrs <- list()
  # eAttrs <- list()

  # # color significant nodes
  # if (!sign.only) {
    # signode <- sapply(nodes, name) %in% significant
    # names(signode) <- names(nodes)
    # sigedge <- (sapply(edges, from) %in% significant) & (sapply(edges, to) %in% significant)
    # names(sigedge) <- names(edges)
    # nodecolor <- ifelse(signode, "black", "#BBBBBB")
    # nAttrs$color <- nodecolor
    # nAttrs$fontcolor <- nodecolor
    
    # edgecolor <- ifelse(sigedge, "black", "#BBBBBB")
    # eAttrs$color <- edgecolor
  # }
  
  # # if no names, give the plot numbers in order
  # # this requires plotting the graph twice
  # if (!names) {
    # nAttrs$label <- 1:length(names(nodes))
    # names(nAttrs$label) <- names(nodes)
    # pg <- agopen(graph, name="pg", nodeAttrs = nAttrs, edgeAttrs = eAttrs)
    # x <- getNodeXY(pg)$x
    # y <- getNodeXY(pg)$y
    # ordering <- sort.list(order(-y, x))
    # nAttrs$label <- ordering
    # names(nAttrs$label) <- names(nodes)
    # plot(graph, attrs = list(node=list(shape="rectangle")), nodeAttrs = nAttrs, edgeAttrs = eAttrs)
  # } else
    # plot(graph, attrs = list(node=list(shape="rectangle")), nodeAttrs = nAttrs, edgeAttrs = eAttrs)
    
  # # Make the plot interactive if asked
  # if (interactive) {
    # cat("Click in the plot to see name and alias. Press escape to return.\n")
    # flush.console()
    # repeat {
      # p <- locator(n = 1) 
      # if (is.null(p))
        # break()
      # pg <- plot(graph, attrs = list(node=list(shape="rectangle")), nodeAttrs = nAttrs, edgeAttrs = eAttrs)
      # x <- getNodeXY(pg)$x
      # y <- getNodeXY(pg)$y
      # distance <- abs(p$x - x) + abs(p$y - y)
      # idx <- which.min(distance)
      # legend("topleft", legend = paste(nAttrs$label[idx], names(object)[idx], alias(object)[idx]), bg = "white", box.lty=0)
    # }
  # }
  
  # # return a legend if the plot has numbers
  # if (!names) {
    # legend <- cbind(name=names(object), alias=alias(object))[sort.list(ordering),]
    # legend <- data.frame(legend, stringsAsFactors=FALSE)
  # } else
    # legend <- NULL
    
  # invisible(legend)
# }
