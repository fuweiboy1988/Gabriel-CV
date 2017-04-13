Gabriel_holdout <- function(data = Data, CV_r = 5, CV_l = 2, method = "SVM",max.k=10){
	L <- ncol(data)
	R <- nrow(data)
	if( round(L/CV_l) != L/CV_l){print(" Column is not multiply of CV number")}
	if( round(R/CV_r) != R/CV_r){print(" Row is not multiply of CV number")}
	
	Step = floor(R/CV_r)
	
	if(CV_l == 2){
		Error.K <- matrix(NA,CV_r,max.k)
		Error.K2 <- matrix(NA,CV_r,max.k)
		ML <- floor(L/2)
		
		for(T in 1:(CV_r-1)){ #k=1
			Xtest <- data[((T-1)*Step+1):(T*Step),1:ML]
			Ytest <- data[((T-1)*Step+1):(T*Step),-(1:ML)]
			Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),1:ML]
			Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(1:ML)]    
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
		##========================== T= CV_r==============================
                  Xtest <- data[((CV_r-1)*Step+1):R,1:ML]
			Ytest <- data[((CV_r-1)*Step+1):R,-(1:ML)]
			Xtrain  <- data[-(((CV_r-1)*Step+1):R),1:ML]
			Ytrain  <- data[-(((CV_r-1)*Step+1):R),-(1:ML)]    
			##-------------------------------
			Cond <- TRUE
			while(Cond){
				Kmean <- try(kmeans(Xtrain, 1, nstart = 100),silent=T)
				if(class(Kmean)!="try-error" ){Cond <- FALSE}
			}
		      ##------------------------------
			Predict.test <- Kmean$center[rep(1,nrow(Xtest)),]
			Error.K[CV_r,1] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
                ####=================================================================
			##-------------------------------
			Cond <- TRUE
			while(Cond){
				Kmean2 <- try(kmeans(Ytrain, 1,nstart = 100),silent=T)
				if( class(Kmean2)!="try-error" ){Cond <- FALSE}
			}
		      ##------------------------------
			Predict.test2 <- Kmean2$center[rep(1,nrow(Ytest)),]
			Error.K2[CV_r,1] <- sum((Predict.test2-Ytest)^2)/(nrow(Ytest)*ncol(Ytest))
                ####=================================================================

		##=======================================================================
		for(k in 2:max.k){
	    		for(T in 1:(CV_r-1)){
				Xtest <- data[((T-1)*Step+1):(T*Step),1:ML]
				Ytest <- data[((T-1)*Step+1):(T*Step),-(1:ML)]
				Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),1:ML]
				Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(1:ML)]    
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
	            ##========================== T= CV_r==============================
                        Xtest <- data[((CV_r-1)*Step+1):R,1:ML]
				Ytest <- data[((CV_r-1)*Step+1):R,-(1:ML)]
				Xtrain  <- data[-(((CV_r-1)*Step+1):R),1:ML]
				Ytrain  <- data[-(((CV_r-1)*Step+1):R),-(1:ML)]    
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

				Error.K[CV_r,k] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
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

				Error.K2[CV_r,k] <- sum((Predict.test2-Ytest)^2)/(nrow(Ytest)*ncol(Ytest))

                ####=================================================================
		}
		
		means <- colMeans((Error.K+Error.K2)/2)
		result <- which(means == min(means))
		
	} else {
			MATRIX <- list(matrix(NA,CV_r,max.k),matrix(NA,CV_r,max.k))
			
			for(t in 3:CV_l){
				MATRIX[[t]] <- matrix(NA,CV_r,max.k)
			}
			
			Step_l = floor(L/CV_l)
			
			#=================================================
			for(LL in 1:(CV_l-1)){
				for(T in 1:(CV_r-1)){ #k=1
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
				      Xtest <- data[((CV_r-1)*Step+1):R,((LL-1)*Step_l+1):(LL*Step_l)]
					Ytest <- data[((CV_r-1)*Step+1):R,-(((LL-1)*Step_l+1):(LL*Step_l))]
					Xtrain  <- data[-(((CV_r-1)*Step+1):R),((LL-1)*Step_l+1):(LL*Step_l)]
					Ytrain  <- data[-(((CV_r-1)*Step+1):R),-(((LL-1)*Step_l+1):(LL*Step_l))]    
					##-------------------------------
					Cond <- TRUE
					while(Cond){
						Kmean <- try(kmeans(Xtrain, 1,nstart = 100),silent=T)
						if( class(Kmean)!="try-error" ){Cond <- FALSE}
					}
		      		##------------------------------
					##Kmean <- kmeans(Xtrain, 1) # the number gonna be variable k later
					Predict.test <- Kmean$center[rep(1,nrow(Xtest)),]
					MATRIX[[LL]][CV_r,1] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
				#------------------------------------------------------
				for(k in 2:max.k){
					for(T in 1:(CV_r-1)){
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
					#--------------------------------------------------------
					      Xtest <- data[((CV_r-1)*Step+1):R,((LL-1)*Step_l+1):(LL*Step_l)]
						Ytest <- data[((CV_r-1)*Step+1):R,-(((LL-1)*Step_l+1):(LL*Step_l))]
						Xtrain  <- data[-(((CV_r-1)*Step+1):R),((LL-1)*Step_l+1):(LL*Step_l)]
						Ytrain  <- data[-(((CV_r-1)*Step+1):R),-(((LL-1)*Step_l+1):(LL*Step_l))]    
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
						MATRIX[[LL]][CV_r,k] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
					#--------------------------------------------------------
				}
				#--------------------------------------------------------
			}#end of LL loop
				for(T in 1:(CV_r-1)){ #k=1
					Xtest <- data[((T-1)*Step+1):(T*Step),((CV_l-1)*Step_l+1):L]
					Ytest <- data[((T-1)*Step+1):(T*Step),-(((CV_l-1)*Step_l+1):L)]
					Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),((CV_l-1)*Step_l+1):L]
					Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(((CV_l-1)*Step_l+1):L)]    
					##-------------------------------
					Cond <- TRUE
					while(Cond){
						Kmean <- try(kmeans(Xtrain, 1,nstart = 100),silent=T)
						if( class(Kmean)!="try-error" ){Cond <- FALSE}
					}
		      		##------------------------------
					##Kmean <- kmeans(Xtrain, 1) # the number gonna be variable k later
					Predict.test <- Kmean$center[rep(1,nrow(Xtest)),]
					MATRIX[[CV_l]][T,1] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
				}
				#------------------------------------------------------
				      Xtest <- data[((CV_r-1)*Step+1):R,((CV_l-1)*Step_l+1):L]
					Ytest <- data[((CV_r-1)*Step+1):R,-(((CV_l-1)*Step_l+1):L)]
					Xtrain  <- data[-(((CV_r-1)*Step+1):R),((CV_l-1)*Step_l+1):L]
					Ytrain  <- data[-(((CV_r-1)*Step+1):R),-(((CV_l-1)*Step_l+1):L)]    
					##-------------------------------
					Cond <- TRUE
					while(Cond){
						Kmean <- try(kmeans(Xtrain, 1,nstart = 100),silent=T)
						if( class(Kmean)!="try-error" ){Cond <- FALSE}
					}
		      		##------------------------------
					##Kmean <- kmeans(Xtrain, 1) # the number gonna be variable k later
					Predict.test <- Kmean$center[rep(1,nrow(Xtest)),]
					MATRIX[[CV_l]][CV_r,1] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
				#------------------------------------------------------
				for(k in 2:max.k){
					for(T in 1:(CV_r-1)){
						Xtest <- data[((T-1)*Step+1):(T*Step),((CV_l-1)*Step_l+1):L]
						Ytest <- data[((T-1)*Step+1):(T*Step),-(((CV_l-1)*Step_l+1):L)]
						Xtrain  <- data[-(((T-1)*Step+1):(T*Step)),((CV_l-1)*Step_l+1):L]
						Ytrain  <- data[-(((T-1)*Step+1):(T*Step)),-(((CV_l-1)*Step_l+1):L)]    
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

						MATRIX[[CV_l]][T,k] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
					}
					#--------------------------------------------------------
					      Xtest <- data[((CV_r-1)*Step+1):R,((CV_l-1)*Step_l+1):L]
						Ytest <- data[((CV_r-1)*Step+1):R,-(((CV_l-1)*Step_l+1):L)]
						Xtrain  <- data[-(((CV_r-1)*Step+1):R),((CV_l-1)*Step_l+1):L]
						Ytrain  <- data[-(((CV_r-1)*Step+1):R),-(((CV_l-1)*Step_l+1):L)]    
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
						MATRIX[[CV_l]][CV_r,k] <- sum((Predict.test-Xtest)^2)/(nrow(Xtest)*ncol(Xtest))
					#--------------------------------------------------------
				}
			#-------------------End of iteration-------------------------------
		means <- colMeans(MATRIX[[1]])
		for (M in 2:CV_l){
		means <- rbind(means, colMeans(MATRIX[[M]]))
		}
		means <- colMeans(as.matrix(means))
		result <- which(means == min(means))	
	}
	
	return(result)
}
