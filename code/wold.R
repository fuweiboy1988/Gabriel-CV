# code/wold.R
#


Impute <- function(Data){
	ID <- which(is.na(Data)==1)
	for(j in 1:length(ID)){
		if(is.na(median(Data[,ceiling(ID[j]/nrow(Data))],na.rm=TRUE))){
			Data[ID[j]] <- 0
		} else {
			Data[ID[j]] <- median(Data[,ceiling(ID[j]/nrow(Data))],na.rm=TRUE)
		}
	}
	return(Data)
}
###===============================================
WoldIter <- function( DATA = Data.impute, Ids = ID.na, k = K, MaxIteration = 1000, ErrorTol = 0.001){
	Oldsegment <- DATA[Ids]
	Continue <- TRUE
      iteration <- 0

      while(Continue){
		iteration <- iteration + 1
		##-------------------------------
		Cond <- TRUE
		while(Cond){
			Kmean <- try(cluster_kmeans(DATA,k),silent=T)
			if(class(Kmean)!="try-error" ){Cond <- FALSE}
		}
		##------------------------------

		DATA_new <- Kmean$center[Kmean$cluster,]
		Newsegment <- DATA_new[Ids]
		Error <- sum((Oldsegment - Newsegment)^2)/length(Oldsegment)
	Continue <- (Error > ErrorTol & iteration <= MaxIteration)
		DATA[Ids] <- Newsegment
		Oldsegment <- Newsegment
	}
	if(iteration >= MaxIteration){print("Convergence doesn't reach")}
	return(Oldsegment)
}
###===============================================
Wold_holdout <- function(data = Data, CV = 5, Errortol = 0.001, max.k=10){
	IDs <- seq(1,nrow(data)*ncol(data),1)

	if(round(length(IDs)/CV) != length(IDs)/CV){
		print("Total observations number is not multiple of CV")
	}

	Id <- sample(IDs,length(IDs))
	Step = length(IDs)/CV
	Error.K <- matrix(NA,CV,max.k)


      for (K in 1:max.k){
		for(T in 1:CV){
			data2 <- as.matrix(data)
			Left.id <- Id[(1+Step*(T-1)):(T*Step)]
			Left.true <- data2[Left.id]
			data2[Left.id] <- NA
			Data.impute <- Impute(data2)
			Left.converge <- WoldIter(Data.impute, Left.id, k=K, ErrorTol = Errortol)
			Error.K[T,K] <- sum((Left.true-Left.converge)^2)/length(Left.id)
		}
	}
	means <- colMeans(Error.K)
	which.min(means)
}

