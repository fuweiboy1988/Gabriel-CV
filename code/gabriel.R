# code/gabriel.R
#
# Depends:
#   library("bcv")
#   library("MASS")
#   library("nnet")
#   source("classify.R")
#   source("cluster.R")
#


cv.kmeans.gabriel <- function(x, krow = 2, kcol = 2, maxcenters = 10,
                              classify.method = c("nearest", "lda-equal", "lda-proportions", "svm"))
{
    classify.method <- match.arg(classify.method)
    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    if (n < 2) 
        stop("x should have at least two rows")
    if (p < 2) 
        stop("x should have at least two columns")
    if ((krow > n) || (krow <= 1)) 
        stop("krow outside allowable range")
    if ((kcol > p) || (kcol <= 1)) 
        stop("kcol outside allowable range")
    if (maxcenters <= 0) 
        stop("maxcenters should be positive")

    krow.o <- krow
    krow <- bcv:::round.fold(n, krow)
    kcol.o <- kcol
    kcol <- bcv:::round.fold(p, kcol)
    if (krow != krow.o) 
        warning("krow has been set to ", krow)
    if (kcol != kcol.o) 
        warning("kcol has been set to ", kcol)

    s.r <- bcv:::choose.sets(n, krow)
    s.c <- bcv:::choose.sets(p, kcol)
    n0 <- n - max(table(s.r))
    p0 <- p - max(table(s.c))
    maxcenters.o <- maxcenters
    maxcenters <- min(n0, round(maxcenters.o))
    if (!missing(maxcenters) && maxcenters != maxcenters.o) 
        warning("maxcenters has been set to ", maxcenters)

    if (classify.method == "nearest") {
        classify <- classify_nearest
    } else if (classify.method == "lda-equal") {
        classify <- function(x, grouping)
            classify_lda(x, grouping, prior="equal")
    } else if (classify.method == "lda-proportions") {
        classify <- function(x, grouping)
            classify_lda(x, grouping, prior="proportions")
    } else if (classify.method == "svm") {
        stop("not yet implemented")
    }

    msep <- matrix(NA, krow * kcol, maxcenters)

    for (k in seq_len(krow)) {
        for (l in seq_len(kcol)) {
            test <- s.r == k
            train <- !test
            response <- s.c == l
            predictor <- !response

            for (centers in seq_len(maxcenters)) {
                if (centers == 1) {
                    fit <- colMeans(x[train,response,drop=FALSE])
                    err <- scale(x[test,response,drop=FALSE], center=fit, scale=FALSE)
                } else {
                    cl <- cluster_kmeans(x[train,response,drop=FALSE], centers=centers)
                    cluster <- factor(cl$cluster, levels=seq_len(centers))
                    fit <- classify(x[train,predictor,drop=FALSE], cluster)
                    pred <- predict(fit, x[test,predictor,drop=FALSE])$class
                    err <- x[test,response,drop=FALSE] - cl$centers[pred,,drop=FALSE]
                }

                msep[k + (l - 1) * krow, centers] <- mean(err^2)
            }
        }
    }

    msep.mean <- colMeans(msep)
    centers <- which.min(msep.mean)
    list(msep, centers=centers)
}

###
Cov <- function(dim = 2, rho = 0){
	COV <- diag(dim)
	for(i in 1:dim){
		for(j in 1:dim){
			if(i!= j){COV[i,j] = rho}
		}
	}
	return (COV)
}

#########---------------------------------

Rotate <- function(Data){
  	data = Data 
	#generate a random rotation
        dim <- ncol(data)
        z <- matrix(rnorm(dim * dim), dim, dim)
        qr <- qr(z)
        q <- qr.Q(qr)
        sign <- sample(c(-1, 1), dim, replace=TRUE)
        rot <- q %*% diag(sign, dim)
	#rotate the columns of the data matrix
        x_rot <- data %*% rot 
 	 return(x_rot)
}

#########---------------------------------

Rotate_center <- function(Data){
	Center <- colMeans(Data)
	Centers <- matrix(Center, nrow(Data), ncol(Data), byrow= T)

	  data = Data - Centers
        x_rot <- Rotate(data) + Centers
	  return(x_rot)
}
####----decorrelate a matrix, i.e.diagonize its covariance matrix-------
Uncorrelate <- function(Data, K){

	Fit <- kmeans(Data, K, nstart = 100)
	ClusterID <- unique(Fit$cluster)
	
	for( j in ClusterID){
		ID <- which(Fit$cluster == j)
		DATA <- Data[ID,]
		Center <- colMeans(DATA)
		DATA <- DATA - matrix(Center, nrow(DATA), ncol(DATA), byrow= T)
		if( !exists("Pool") ){
			Pool <- DATA
		}else{
			Pool <- rbind(Pool,DATA)
		}
	}

	Decomp <- eigen(cov(Pool))
	Diag_matrix <- Decomp$vectors
	eigen_value <- sqrt(1/Decomp$values)

	result <- Data%*%Diag_matrix%*%diag(eigen_value)

	##-----------------------------------------
	# result <- Rotate(result)

	return(result)
} 

Uncorrelate2 <- function(Data,K){
	Fit <- kmeans(Data, K, nstart = 100)
	Pool <- Data - Fit$center[Fit$cluster,]

	SVD <- svd(Pool)
	E = 10^(-5)
	result <- Data%*%SVD$v%*%(diag(1/(SVD$d+E)))
	result <- Rotate(result)

	return(result)
} 
###-------------------Gabriel correlation correction-------------------------###########

gabriel_cor_correct <- function(data, maxcenters, type = 1){

	K.original <- cv.kmeans.gabriel(data, 5, 2, maxcenters, classify.method="nearest")$centers
	
	if(type == 1){
		DATA <- Uncorrelate(data,K.original)
	}else if(type == 2){
		DATA <- Uncorrelate2(data,K.original)
	}
	
	K2 <- cv.kmeans.gabriel(DATA, 5, 2, maxcenters, classify.method="nearest")$centers
	
	return(K2)
}

###
