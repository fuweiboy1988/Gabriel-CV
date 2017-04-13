##Type=1 means null case (single-cluster in 10 dim)
##Type=4 means 4 cluster in 4 dimension cases
##Type=10 means 4 cluster in 10 dimension cases

library(mclust)
library(cluster)
library(MASS)
library(e1071)
library(LICORS)
library(fpc)
library(NbClust)
##------------------------------Jump-----------------------------------------------------
"jump" <-
function(data=NULL,K=10,y=NULL,plotjumps=T,rand=10,fits=NULL,B=0,dist=NULL,trace=F){
  if (!is.null(data)){
    # Compute the kmeans fit to the data
  if (is.null(fits))
  fits <- kmeans.rndstart(data,K,rand)
  if (is.null(y))
    y <- dim(data)[2]/2
  n <- nrow(data)
  p <- ncol(data)
  # Compute the distortion associated with the kmeans fit
  dist<- fits/(n*p)
}
  # Call the compute.jump function to produce plots and calculate
  # maximum jump for each value of Y
  jump.results <- compute.jump(dist,y,plotjumps)
  jump.results$fits <- fits
  # Implement bootstrap routine
  if (B>0 & !is.null(data)){
  n <- nrow(data)
  boot.results <- matrix(0,length(y),K)
  bootdist <- rep(0,K)
  for (b in 1:B){
    if (trace)
    print(paste("Bootstrap Iteration ",b))
    # Make bootstrap data
    bootdata <- data[sample(1:n,replace=T),]
    # Get kmeans fit to the bootstrap data
    bootfits <- kmeans.rndstart(bootdata,K,rand)
    # Compute bootstrap distortion and maximum jumps
    for (k in 1:K)
      bootdist[k] <- sum(bootfits[[k]]$within)/(n*p)
    bootmaxjumps <- compute.jump(bootdist,y,plotjumps=F,printresults=F)$maxjump
    for (j in 1:length(y))
      boot.results[j,bootmaxjumps[j]] <-  boot.results[j,bootmaxjumps[j]]+1
  }
  # Calculate proportions of each number of clusters chosen
  jump.results$boot.result <- round(boot.results/B,3)
  for (j in 1:length(y))
     print(paste(jump.results$boot.result[j,jump.results$maxjump[j]]*100,"% of bootstrap iterations corresponding to ",jump.results$maxjump[j], "clusters with Y=",y[j]))
     }
  jump.results}

"compute.jump" <-
function(dist,y,plotjumps=T,printresults=T){
    K <- length(dist)
    numb.y <- length(y)
    numbclust <- rep(0,numb.y)
    transdist <- matrix(0,numb.y,K+1)
    jumps <- matrix(0,numb.y,K)
    if (plotjumps)
      par(mfrow=c(numb.y,3))
    for (i in 1:numb.y){
      # Compute the transformed distortion
      transdist[i,] <- c(0,dist^(-y[i]))
      # Compute the jumps in transformed distortion
      jumps[i,] <- diff(transdist[i,])
      # Compute the maximum jump
      numbclust[i] <- order(-jumps[i,])[1]
      # Plot distortion, transformed distortion and jumps
      if (plotjumps){
        plot(1:K,dist,type='l',
             xlab="Number of Clusters",ylab="Distortion",main=paste("Y = ",y[i]))
        plot(0:K,transdist[i,],type='l',
             xlab="Number of Clusters",ylab="Transformed Distortion",main=paste("Y = ",y[i]))
        plot(1:K,jumps[i,],type='l',
             xlab="Number of Clusters",ylab="Jumps",main=paste("Y = ",y[i]))
        # Plot line and point to indicate maximum jump
        lines(rep(numbclust[i],2),c(0,jumps[i,numbclust[i]]),lty=3,lwd=3)
        points(numbclust[i],jumps[i,numbclust[i]],col=2,pch=19,cex=1.5)}
      # Report maximum jump
      if (printresults)
    print(paste("The maximum jump occurred at ",numbclust[i], "clusters with Y=",y[i]))
    }
    list(maxjump=numbclust,dist=dist,transdist=transdist[,-1],jumps=jumps)}

"kmeans.rndstart" <-
function(x, K, rand = 10)
{
  fits <- sum((t(x) - apply(x, 2, mean))^2)
  iter.max <- 10
  # Run kmeans for 2 to K clusters
  for (k in 2:K){
  Z=kmeans(x,k,nstart=rand)
  fits=c(fits,sum(Z$withinss))
}
fits
}

##-----------------------------Gabriel----------------------------
Gabriel_holdout <- function(data = Data, CV_r = 5, CV_l = 2, method = "SVM",max.k=10){
	L <- ncol(data)
	R <- nrow(data)
	if( round(L/CV_l) != L/CV_l){print(" Column is not multiply of CV number")}
	if( round(R/CV_r) != R/CV_r){print(" Row is not multiply of CV number")}
	
	Step = R/CV_r
	
	if(CV_l == 2){
		Error.K <- matrix(NA,CV_r,max.k)
		Error.K2 <- matrix(NA,CV_r,max.k)

		for(T in 1:CV_r){ #k=1
			Xtest <- data[((T-1)*Step+1):(T*Step),1:(L/2)]
			Ytest <- data[((T-1)*Step+1):(T*Step),-(1:(L/2))]
			Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),1:(L/2)]
			Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(1:(L/2))]    
			##-------------------------------
			Cond <- TRUE
			while(Cond){
				Kmean <- try(kmeans(Xtrain, 1, nstart = 100),silent=T)
				if(class(Kmean)!="try-error" ){Cond <- FALSE}
			}
		      ##------------------------------
			Predict.test <- Kmean$center[rep(1,nrow(Xtest)),]
			Error.K[T,1] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
                ####=================================================================
			##-------------------------------
			Cond <- TRUE
			while(Cond){
				Kmean2 <- try(kmeans(Ytrain, 1,nstart = 100),silent=T)
				if( class(Kmean2)!="try-error" ){Cond <- FALSE}
			}
		      ##------------------------------
			Predict.test2 <- Kmean2$center[rep(1,nrow(Ytest)),]
			Error.K2[T,1] <- sum((Predict.test2-Ytest)^2)/(nrow(Ytest)*ncol(Ytest))
                ####=================================================================
		}
		
		for(k in 2:max.k){
	    		for(T in 1:CV_r){
				Xtest <- data[((T-1)*Step+1):(T*Step),1:(L/2)]
				Ytest <- data[((T-1)*Step+1):(T*Step),-(1:(L/2))]
				Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),1:(L/2)]
				Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(1:(L/2))]    
				##-------------------------------
				Cond <- TRUE
				while(Cond){
					Kmean <- try(kmeans(Xtrain, k, nstart = 100),silent=T)
					if( class(Kmean)!="try-error" ){Cond <- FALSE}
				}
		      	##------------------------------
				##Kmean <- kmeanspp(Xtrain, k)
				Center = Kmean$cluster
				
				if(method == "LDA"){
					LDA <- lda(Ytrain,factor(Center))
					Predict.test <- Kmean$center[predict(LDA,Ytest)$class,]
				} else{
					SVM.fit <- svm(Ytrain,factor(Center))
					Predict.test <- Kmean$center[predict(SVM.fit,Ytest),]
				}

				Error.K[T,k] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
			   ####=================================================================
				##-------------------------------
				Cond <- TRUE
				while(Cond){
					Kmean2 <- try(kmeans(Ytrain, k, nstart = 100),silent=T)
					if( class(Kmean2)!="try-error" ){Cond <- FALSE}
				}
		      	##------------------------------
				##Kmean2 <- kmeanspp(Ytrain, k) # the number gonna be variable k later
				Center2 = Kmean2$cluster
				
				if(method == "LDA"){
					LDA2 <- lda(Xtrain,factor(Center2))
					Predict.test2 <- Kmean2$center[predict(LDA2,Xtest)$class,]
				} else{
					SVM.fit2 <- svm(Xtrain,factor(Center2))
					Predict.test2 <- Kmean2$center[predict(SVM.fit2,Xtest),]
				}

				Error.K2[T,k] <- sum((Predict.test2-Ytest)^2)/(nrow(Ytest)*ncol(Ytest))
			  ####=================================================================
			}
		}
		
		means <- colMeans((Error.K+Error.K2)/2)
		result <- which(means == min(means))
		
	} else {
			MATRIX <- list(matrix(NA,CV_r,max.k),matrix(NA,CV_r,max.k))
			
			for(t in 3:CV_l){
				MATRIX[[t]] <- matrix(NA,CV_r,max.k)
			}
			
			Step_l = L/CV_l
			
			#=================================================
			for(LL in 1:CV_l){
				for(T in 1:CV_r){ #k=1
					Xtest <- data[((T-1)*Step+1):(T*Step),((LL-1)*Step_l+1):(LL*Step_l)]
					Ytest <- data[((T-1)*Step+1):(T*Step),-(((LL-1)*Step_l+1):(LL*Step_l))]
					Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),((LL-1)*Step_l+1):(LL*Step_l)]
					Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(((LL-1)*Step_l+1):(LL*Step_l))]    
					##-------------------------------
					Cond <- TRUE
					while(Cond){
						Kmean <- try(kmeans(Xtrain, 1,nstart = 100),silent=T)
						if( class(Kmean)!="try-error" ){Cond <- FALSE}
					}
		      		##------------------------------
					##Kmean <- kmeans(Xtrain, 1) # the number gonna be variable k later
					Predict.test <- Kmean$center[rep(1,nrow(Xtest)),]
					MATRIX[[LL]][T,1] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
				}
				#------------------------------------------------------
				for(k in 2:max.k){
					for(T in 1:CV_r){
						Xtest <- data[((T-1)*Step+1):(T*Step),((LL-1)*Step_l+1):(LL*Step_l)]
						Ytest <- data[((T-1)*Step+1):(T*Step),-(((LL-1)*Step_l+1):(LL*Step_l))]
						Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),((LL-1)*Step_l+1):(LL*Step_l)]
						Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(((LL-1)*Step_l+1):(LL*Step_l))]    
						##-------------------------------
						Cond <- TRUE
						while(Cond){
							Kmean <- try( kmeans(Xtrain, k,nstart = 100) ,silent=T)
							if( class(Kmean)!="try-error"){Cond <- FALSE}
						}
		      			##------------------------------
						##Kmean <- kmeanspp(Xtrain, k) # the number gonna be variable k later
						Center = Kmean$cluster
				
						if(method == "LDA"){
							LDA <- lda(Ytrain,factor(Center))
							Predict.test <- Kmean$center[predict(LDA,Ytest)$class,]
						} else{
							SVM.fit <- svm(Ytrain,factor(Center))
							Predict.test <- Kmean$center[predict(SVM.fit,Ytest),]
						}

						MATRIX[[LL]][T,k] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
					}
				}
				#--------------------------------------------------------
				
			}
		means <- colMeans(MATRIX[[1]])
		for (M in 2:CV_l){
		means <- rbind(means, colMeans(MATRIX[[M]]))
		}
		means <- colMeans(as.matrix(means))
		result <- which(means == min(means))	
	}
	
	return(result)
}

##------------------------  Wold ---------------------------------
###=============================================
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
			Kmean <- try(kmeans(DATA,k, nstart = 100),silent=T)
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
	result <- which(means == min(means))
      return(result)
}
#########====================log normal ================================

Data4lognormal <- function(di = 16, dd = 1.2){ 
	condition <- TRUE
  	while(condition){
 		del = 1
 		nn=sample(c(30,60),size=4,replace=T)
 		cl=c(rep(1,nn[1]),rep(2,nn[2]),rep(3,nn[3]),rep(4,nn[4]))
 		
 		c1=dd*rnorm(di)
		u1=matrix(c1,nrow=nn[1],ncol=di,byrow=T)##
 		x1=matrix(c1,nrow=nn[1],ncol=di,byrow=T) + matrix(rlnorm(nn[1]*di, 0, 0.5),ncol=di) - exp(0.25/2)
 		c2=dd*rnorm(di)
		u2=matrix(c2,nrow=nn[2],ncol=di,byrow=T)##
 		x2=matrix(c2,nrow=nn[2],ncol=di,byrow=T) + matrix(rlnorm(nn[2]*di, 0, 0.5),ncol=di) - exp(0.25/2)
 		c3=dd*rnorm(di)
		u3=matrix(c3,nrow=nn[3],ncol=di,byrow=T)##
 		x3=matrix(c3,nrow=nn[3],ncol=di,byrow=T) + matrix(rlnorm(nn[3]*di, 0, 0.5),ncol=di)- exp(0.25/2)
 		c4=dd*rnorm(di)
		u4=matrix(c4,nrow=nn[4],ncol=di,byrow=T)##
 		x4=matrix(c4,nrow=nn[4],ncol=di,byrow=T) + matrix(rlnorm(nn[4]*di, 0, 0.5),ncol=di)- exp(0.25/2)
 
 		x=rbind(x1,x2,x3,x4)
 		Mu=rbind(u1,u2,u3,u4) 

 		ss=dist(rbind(x,c1,c2,c3,c4))
 		d=matrix(0,nrow=nrow(x)+4,ncol=nrow(x)+4)
 		d[row(d)>col(d)]=ss
 		DD=d[(nrow(x)+1):(nrow(x)+4),1:nrow(x)]
 
 		for(i in 1:ncol(DD)){ DD[cl[i],i]=DD[cl[i],i]+del}
 
 		ff=apply(DD,2,which.min)
       	if(sum(ff==cl)==nrow(x)){
 			condition <- FALSE
 		}
    	}

    ID <- sample(1:nrow(x),nrow(x),replace = FALSE)
    Data <- x[ID,]
    Mean <- Mu[ID,]
    result <- list(data = Data, mean = Mean)
    return(result)
}
#------------------------------------------------------------------
location <- function(M.data,M.center){
   if(ncol(M.data)!=ncol(M.center)){return("dimention don't match")}
   center<-rep(NA,nrow(M.data))
   K <- nrow(M.center)
   for (i in 1:length(center)){
     dist <- rep(NA,K)
     for(j in 1:K){
         dist[j]<- sqrt(sum((M.center[j,]-M.data[i,])^2))
     }
     center[i]<- which(dist == min(dist))
   }
   return(center)
}
#-----------------------------------------------------
Data2center <- function(number=50){
	u1 = c(1,0,0,1)
	u2 = c(1,3.5,3.5,1)
      Centers <- matrix(c(1,1,0,3.5,0,3.5,1,1),2,4)

	condition <- TRUE
  	while(condition){
	      obs1 <- mvrnorm(n = number, u1, diag(4))
		obs2 <- mvrnorm(n = number, u2, diag(4))	

		if(sum(location(obs1,Centers)==1)==number & sum(location(obs2,Centers)==2)==number){
			condition <- FALSE
		}
	}

	Data <- rbind(obs1,obs2)
	Mu <- rbind(matrix(u1,nrow=number,ncol=4,byrow=T),matrix(u2,nrow=number,ncol=4,byrow=T))
	ID <- sample(1:nrow(Data),nrow(Data),replace = FALSE)
	Data <- Data[ID,]

	Mean <- Mu[ID,]
    	result <- list(data = Data, mean = Mean)
    	return(result)
}	
##------------------------------------------------------------
##------------------------------------------------------------
Data3center <- function(number=40){
     Centers <- mvrnorm(3, rep(0,20), 19*diag(20))
     u1 <- Centers[1,]
     u2 <- Centers[2,]
     u3 <- Centers[3,]

     condition <- TRUE	
     while(condition){
	      obs1 <- matrix(u1,nrow=number,ncol=20,byrow=T)+ matrix(rexp(number*20, rate = 1),ncol=20)-1
		U1 <- matrix(u1,nrow=number,ncol=20,byrow=T)
		obs2 <- matrix(u2,nrow=number,ncol=20,byrow=T)+ matrix(rexp(number*20, rate = 1/2),ncol=20)-2
		U2 <- matrix(u2,nrow=number,ncol=20,byrow=T)
		obs3 <- matrix(u3,nrow=number,ncol=20,byrow=T)+ matrix(rexp(number*20, rate = 1/5),ncol=20)-5
		U3 <-  matrix(u3,nrow=number,ncol=20,byrow=T)

		if(sum(location(obs1,Centers)==1)==number & sum(location(obs2,Centers)==2)==number & sum(location(obs3,Centers)==3)==number){
			condition <- FALSE
		}
	}
	Data <- rbind(obs1,obs2,obs3)
	Mu <- rbind(U1,U2,U3)
	ID <- sample(1:nrow(Data),nrow(Data),replace = FALSE)
	Data <- Data[ID,]

	Mean <- Mu[ID,]
    	result <- list(data = Data, mean = Mean)
    	return(result)
}
##------------------------------------------------------------
Data10center <- function(di = 100, dd = 0.72){ 
	
	condition <- TRUE
  	while(condition){
 		del = 1
		nn=sample(c(50,100),size=10,replace=T)
 		cl=c(rep(1,nn[1]),rep(2,nn[2]),rep(3,nn[3]),rep(4,nn[4]),rep(5,nn[5]),rep(6,nn[6]),rep(7,nn[7]),rep(8,nn[8]),rep(9,nn[9]),rep(10,nn[10]))
 
 		c1=dd*rnorm(di)
		U1=matrix(c1,nrow=nn[1],ncol=di,byrow=T) 
 		x1=matrix(c1,nrow=nn[1],ncol=di,byrow=T) + matrix(rnorm(nn[1]*di),ncol=di)
 		c2=dd*rnorm(di)
		U2=matrix(c2,nrow=nn[2],ncol=di,byrow=T)
 		x2=matrix(c2,nrow=nn[2],ncol=di,byrow=T) + matrix(rnorm(nn[2]*di),ncol=di)
 		c3=dd*rnorm(di)
		U3=matrix(c3,nrow=nn[3],ncol=di,byrow=T)
 		x3=matrix(c3,nrow=nn[3],ncol=di,byrow=T) + matrix(rnorm(nn[3]*di),ncol=di)
 		c4=dd*rnorm(di)
		U4=matrix(c4,nrow=nn[4],ncol=di,byrow=T)
 		x4=matrix(c4,nrow=nn[4],ncol=di,byrow=T) + matrix(rnorm(nn[4]*di),ncol=di)
 		c5=dd*rnorm(di)
		U5=matrix(c5,nrow=nn[5],ncol=di,byrow=T)
		x5=matrix(c5,nrow=nn[5],ncol=di,byrow=T) + matrix(rnorm(nn[5]*di),ncol=di)
		c6=dd*rnorm(di)
		U6=matrix(c6,nrow=nn[6],ncol=di,byrow=T)
 		x6=matrix(c6,nrow=nn[6],ncol=di,byrow=T) + matrix(rnorm(nn[6]*di),ncol=di)
 		c7=dd*rnorm(di)
		U7=matrix(c7,nrow=nn[7],ncol=di,byrow=T)
 		x7=matrix(c7,nrow=nn[7],ncol=di,byrow=T) + matrix(rnorm(nn[7]*di),ncol=di)
 		c8=dd*rnorm(di)
		U8=matrix(c8,nrow=nn[8],ncol=di,byrow=T)
 		x8=matrix(c8,nrow=nn[8],ncol=di,byrow=T) + matrix(rnorm(nn[8]*di),ncol=di)
 		c9=dd*rnorm(di)
		U9=matrix(c9,nrow=nn[9],ncol=di,byrow=T)
 		x9=matrix(c9,nrow=nn[9],ncol=di,byrow=T) + matrix(rnorm(nn[9]*di),ncol=di)
 		c10=dd*rnorm(di)
		U10=matrix(c10,nrow=nn[10],ncol=di,byrow=T)
		x10=matrix(c10,nrow=nn[10],ncol=di,byrow=T) + matrix(rnorm(nn[10]*di),ncol=di)


 		x=rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)
 		Mu=rbind(U1,U2,U3,U4,U5,U6,U7,U8,U9,U10)

 		ss=dist(rbind(x,c1,c2,c3,c4,c5,c6,c7,c8,c9,c10))
 		d=matrix(0,nrow=nrow(x)+10,ncol=nrow(x)+10)
 		d[row(d)>col(d)]=ss
 		DD=d[(nrow(x)+1):(nrow(x)+10),1:nrow(x)]
 
 		for(i in 1:ncol(DD)){ DD[cl[i],i]=DD[cl[i],i]+del}
 
 		ff=apply(DD,2,which.min)
       	if(sum(ff==cl)==nrow(x)){
 			condition <- FALSE
 		}
    	}
    
    ID <- sample(1:nrow(x),nrow(x),replace = FALSE)
    Data <- x[ID,]
    Mean <- Mu[ID,]

	result <- list(data = Data, mean = Mean)
    	return(result)
}

##----------------------------------------------------------------
Data4center <- function(di = 10, dd = 1.9){ #Another pair is di=4 and dd=5
	condition <- TRUE
  	while(condition){
 		del = 1
		if(di == 100){
			nn=sample(c(100,150),size=4,replace=T)
		}else{
			nn=sample(c(25,50),size=4,replace=T)
		}
 		cl=c(rep(1,nn[1]),rep(2,nn[2]),rep(3,nn[3]),rep(4,nn[4]))
 
 		c1=dd*rnorm(di)
		U1=matrix(c1,nrow=nn[1],ncol=di,byrow=T) 
 		x1=matrix(c1,nrow=nn[1],ncol=di,byrow=T) + matrix(rnorm(nn[1]*di),ncol=di)
 		c2=dd*rnorm(di)
		U2=matrix(c2,nrow=nn[2],ncol=di,byrow=T)
 		x2=matrix(c2,nrow=nn[2],ncol=di,byrow=T) + matrix(rnorm(nn[2]*di),ncol=di)
 		c3=dd*rnorm(di)
		U3=matrix(c3,nrow=nn[3],ncol=di,byrow=T)
 		x3=matrix(c3,nrow=nn[3],ncol=di,byrow=T) + matrix(rnorm(nn[3]*di),ncol=di)
 		c4=dd*rnorm(di)
		U4=matrix(c4,nrow=nn[4],ncol=di,byrow=T)
 		x4=matrix(c4,nrow=nn[4],ncol=di,byrow=T) + matrix(rnorm(nn[4]*di),ncol=di)
 
 		x=rbind(x1,x2,x3,x4)
		Mu = rbind(U1,U2,U3,U4) 	

 		ss=dist(rbind(x,c1,c2,c3,c4))
 		d=matrix(0,nrow=nrow(x)+4,ncol=nrow(x)+4)
 		d[row(d)>col(d)]=ss
 		DD=d[(nrow(x)+1):(nrow(x)+4),1:nrow(x)]
 
 		for(i in 1:ncol(DD)){ DD[cl[i],i]=DD[cl[i],i]+del}
 
 		ff=apply(DD,2,which.min)
       	if(sum(ff==cl)==nrow(x)){
 			condition <- FALSE
 		}
    	}
    
    ID <- sample(1:nrow(x),nrow(x),replace = FALSE)
    Data <- x[ID,]
    Mean <- Mu[ID,]

	result <- list(data = Data, mean = Mean)
    	return(result)
}

#---------------------------------------------------------------------------

Picks <- function(vector){
	KM <- rep(NA,15)
	for(h in 1:15){
 		KM[h] <- sum(vector == h,na.rm=T)
	}
	return(KM) 
}
#--------------------------------------------------------------------------
##----------------------- Comprehensive function ----------------

Cluster.master <- function(Type = 1, Cv = 5, Cv_l = 2, Method = "LDA", Kmax = 10){ ##Method can be "LDA" for Gabriel
	Ora.K <- rep(NA,100)
	Ora.error <- rep(NA,100)
	Gap.K <- rep(NA,100)
	Gap.error <- rep(NA,100)
	BIC.K <- rep(NA,100)
	BIC.error <- rep(NA,100)
	CH.K <- rep(NA,100)
	CH.error <- rep(NA,100)
	Hart.K <- rep(NA,100)
	Hart.error <- rep(NA,100)
	Jump.K <- rep(NA,100)
	Jump.error <- rep(NA,100)
	Ps.K <- rep(NA,100)
	Ps.error <- rep(NA,100)
	Stab.K <- rep(NA,100)
	Stab.error <- rep(NA,100)
	GG.K <- rep(NA,100)
	GG.error <- rep(NA,100)
	Wold.K <- rep(NA,100)
	Wold.error <- rep(NA,100)


	SIMULATE <- list()
	PE.rtv <-list()

	set.seed(0)
	
	for (W in 1:100){
		if( Type==1){
			Simulate <- matrix(NA,200,10)
			for(s in 1:10){
				Simulate[,s]<-runif(200,0,1)
			}
			Truth <- matrix(rep(0.5,10),nrow=200,ncol=10,byrow=T)
		}else if(Type==2){
			SIM <- Data2center(50)
			Simulate <- SIM$data
			Truth <- SIM$mean
		}else if(Type==3){
			SIM <- Data3center(40)
			Simulate <- SIM$data
			Truth <- SIM$mean
		}else if(Type==100){
			SIM <- Data4center(di = 100, dd = 0.65)
			Simulate <- SIM$data
			Truth <- SIM$mean
		}else if(Type=="10center"){
			SIM <- Data10center(di = 100, dd = 0.72)
			Simulate <- SIM$data
			Truth <- SIM$mean
		}else if (Type == "Lognormal"){
			SIM <- Data4lognormal(16,1.2)
			Simulate <- SIM$data
			Truth <- SIM$mean
		}else{
			print("Type is not correct!")
			break
		}

		SIMULATE[[W]] <- Simulate

		##----------------------------------------------
			PE <- rep(NA,Kmax)
			for(i in 1:Kmax){
				KMean <- kmeans(Simulate, i, nstart = 100)
				PE.pred <- KMean$centers[KMean$cluster,]
				PE[i] <- sum((PE.pred-Truth)^2)/(nrow(Truth)*ncol(Truth))	
			}
			Ora.K[W] <- which(PE == min(PE))
			PE <- PE/min(PE)
			Ora.error[W]<-PE[Ora.K[W]]
			PE.rtv[[W]] <- PE
            ##----------------------------------------------
		
	}
	
	for(Z in 1:100){
			set.seed(1)
			GG.K[Z] <- Gabriel_holdout(data = SIMULATE[[Z]], CV_r = Cv, CV_l = Cv_l, method = Method, max.k = Kmax )
			GG.error[Z] <- PE.rtv[[Z]][GG.K[Z]]
			
			set.seed(1)
			Wold.K[Z] <- Wold_holdout(data = SIMULATE[[Z]], CV = Cv, Errortol = 0.01,  max.k = Kmax)
			Wold.error[Z] <- PE.rtv[[Z]][Wold.K[Z]]
		
			set.seed(1)
			Gap <- clusGap(SIMULATE[[Z]], FUN = kmeans, K.max = Kmax)
			Gap.K[Z] <- which(Gap[[1]][,3] == max(Gap[[1]][,3]))
			Gap.error[Z] <- PE.rtv[[Z]][Gap.K[Z]]
	
			set.seed(1)
			mcluster <- Mclust(SIMULATE[[Z]], G = 1:Kmax)
			BIC.K[Z] <- mcluster$G
			BIC.error[Z] <- PE.rtv[[Z]][BIC.K[Z]]				

			set.seed(1)
			Ch <- NbClust(SIMULATE[[Z]], min.nc = 2, max.nc = Kmax,method = "kmeans", index = "ch")
			CH.K[Z] <- Ch$Best.nc[1]
			CH.error[Z]<-PE.rtv[[Z]][CH.K[Z]]	
	
			set.seed(1)
			Hartigan <- NbClust(SIMULATE[[Z]], min.nc = 2, max.nc = Kmax,method = "kmeans", index = "hartigan")
			Hart.K[Z] <- Hartigan$Best.nc[1]
			Hart.error[Z]<-PE.rtv[[Z]][Hart.K[Z]]

			set.seed(1)
			Jump <- jump(SIMULATE[[Z]],plotjumps=F,rand=10,trace=F)
			Jump.K[Z] <- Jump$maxjump
			Jump.error[Z]<-PE.rtv[[Z]][Jump.K[Z]]

			set.seed(1)
			PS <- prediction.strength(SIMULATE[[Z]], Gmin=2, Gmax=Kmax)
			Ps.K[Z] <- PS$optimalk
			Ps.error[Z] <- PE.rtv[[Z]][Ps.K[Z]]

			set.seed(1)
			SB <- nselectboot(SIMULATE[[Z]],clustermethod=kmeansCBI,classification="centroid",krange=2:Kmax)
			Stab.K[Z] <- SB$kopt
			Stab.error[Z]<-PE.rtv[[Z]][Stab.K[Z]]
	}
	result <- list(N1=Ora.K, E1=Ora.error, N2=Gap.K, E2=Gap.error, N3=BIC.K, E3=BIC.error, N4=CH.K, E4=CH.error, N5=Hart.K, E5=Hart.error, N6=Jump.K, E6=Jump.error,N7=Ps.K, E7=Ps.error, N8=Stab.K, E8=Stab.error,N9=GG.K, E9=GG.error ,N10=Wold.K, E10=Wold.error)
	return(result)
}
##-----------------------End of Comprehensive function ----------------
##---------------------------------------------------------------------
##---------------------------------------------------------------------
##---------------------------------------------------------------------
##---------------------------------------------------------------------
##---------------------------------------------------------------------
##----------Example to get result of 4th simulation situation (10 clusters in 100 dimensions)

RT.10 <- Cluster.master(Type = "10center",Kmax = 15)

RT <- RT.10
#------------------------------------------------------

	print("##----Oracle-----##")
	print(Picks(RT$N1))
	print(mean(RT$E1))
	print(sd(RT$E1))
      print("##----Gap-----##")
	print(Picks(RT$N2))
	print(mean(RT$E2))
	print(sd(RT$E2))
	print("##----BIC-----##")
	print(Picks(RT$N3))
	print(mean(RT$E3))
	print(sd(RT$E3))
	print("##----CH-----##")
	print(Picks(RT$N4))
	print(mean(RT$E4))
	print(sd(RT$E4))
	print("##----Hartingan-----##")
	print(Picks(RT$N5))
	print(mean(RT$E5))
	print(sd(RT$E5))
	print("##----Jump-----##")
	print(Picks(RT$N6))
	print(mean(RT$E6))
	print(sd(RT$E6))
	print("##----Ps-----##")
	print(Picks(RT$N7))
	print(mean(RT$E7))
	print(sd(RT$E7))
	print("##----Stability-----##")
	print(Picks(RT$N8))
	print(mean(RT$E8))
	print(sd(RT$E8))
	print("##----Gabriel-----##")
	print(Picks(RT$N9))
	print(mean(RT$E9))
	print(sd(RT$E9))
	print("##----Wold-----##")
	print(Picks(RT$N10))
	print(mean(RT$E10))
	print(sd(RT$E10))
	print(''#------------------#'')

#------------------------------------------------------
#--------------------------------------------------------------------------